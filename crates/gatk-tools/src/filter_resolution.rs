//! `GATKReadFilterPluginDescriptor`: which read filters a command line ends up with.
//!
//! Four arguments decide it, and a port reading two of them silently ignores the other two:
//! `--read-filter` and `--disable-tool-default-read-filters` were read, `--disable-read-filter`
//! and `--inverted-read-filter` were not.
//!
//! # The order is defaults, then enabled, then inverted
//!
//! Disabling removes from the DEFAULTS only, and before anything is added. The user's enabled
//! filters follow in the order given, and the inverted ones after them. Nothing is sorted.
//!
//! # An enabled filter already among the defaults is not added twice
//!
//! `getResolvedInstances` filters with `contains()` on the INSTANCES, and every library filter is
//! a singleton, so the test is identity and it agrees with the name here. A filter enabled twice
//! never reaches that test: it is refused first.
//!
//! # An inverted filter is a different filter
//!
//! `negate()` answers a `ReadFilterNegate` wrapping the original, and it is APPENDED rather than
//! replacing anything, so a list can hold a filter and its negation at once. Inverting a filter
//! the tool takes by default is refused outright unless the defaults were disabled, "so we do not
//! inadvertently filter all reads from the input".
//!
//! # One of the refusals names the wrong set
//!
//! The check for a filter that is both enabled and inverted builds `enabledAndInverted` and then
//! formats `enabledAndDisabled`, so its message lists the empty set: `The read filter(s):  are
//! both enabled and inverted`, with two spaces. That is the reference's behaviour and this
//! reproduces it rather than repairing it.
//!
//! Ported from
//! `org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor`.

/// One filter of the resolved list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedFilter {
    pub name: String,
    /// Whether it arrived through `--inverted-read-filter`, which wraps it in a `ReadFilterNegate`.
    pub negated: bool,
}

/// What the descriptor refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FilterResolutionError {
    EnabledTwice {
        names: Vec<String>,
    },
    DisabledTwice {
        names: Vec<String>,
    },
    EnabledAndDisabled {
        names: Vec<String>,
    },
    /// The one whose message lists the empty set, because the reference formats the wrong variable.
    EnabledAndInverted,
    DisabledDoesNotExist {
        name: String,
    },
    InvertedIsADefault {
        names: Vec<String>,
    },
    Unrecognized {
        name: String,
    },
}

impl FilterResolutionError {
    pub fn java_class(&self) -> &'static str {
        match self {
            FilterResolutionError::EnabledTwice { .. }
            | FilterResolutionError::DisabledTwice { .. }
            | FilterResolutionError::DisabledDoesNotExist { .. } => {
                "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
            }
            _ => "org.broadinstitute.barclay.argparser.CommandLineException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            FilterResolutionError::EnabledTwice { names } => format!(
                "Illegal argument value: The read filter(s) are enabled more than once: {}",
                names.join(", ")
            ),
            FilterResolutionError::DisabledTwice { names } => format!(
                "Illegal argument value: The read filter(s) are disabled more than once: {}",
                names.join(", ")
            ),
            FilterResolutionError::EnabledAndDisabled { names } => format!(
                "The read filter(s): {} are both enabled and disabled",
                names.join(", ")
            ),
            // The empty list is the reference's own: it formats `enabledAndDisabled`, which is
            // empty by the time this branch is reached.
            FilterResolutionError::EnabledAndInverted => {
                "The read filter(s):  are both enabled and inverted".to_string()
            }
            FilterResolutionError::DisabledDoesNotExist { name } => {
                format!("Illegal argument value: Disabled filter ({name}) does not exist")
            }
            // `%s` on a `Set`, which is Java's `[a, b]`.
            FilterResolutionError::InvertedIsADefault { names } => format!(
                "The read filter(s): [{}] exist as tool default filters and were inverted, \
                 disable tool default read filters in order to use this filter",
                names.join(", ")
            ),
            FilterResolutionError::Unrecognized { name } => {
                format!("Unrecognized read filter name: {name}")
            }
        }
    }
}

/// The names that appear more than once, in first-seen order.
fn duplicated(names: &[String]) -> Vec<String> {
    let mut seen = Vec::new();
    let mut twice = Vec::new();
    for name in names {
        if seen.contains(name) {
            if !twice.contains(name) {
                twice.push(name.clone());
            }
        } else {
            seen.push(name.clone());
        }
    }
    twice
}

/// `validateAndResolvePlugins` followed by `getResolvedInstances`.
///
/// `known` is every filter name the library carries, which is what an unknown name is measured
/// against; `defaults` is the tool's own list, in the order it declares them.
pub fn resolve(
    defaults: &[&str],
    known: &[&str],
    enabled: &[String],
    disabled: &[String],
    inverted: &[String],
    disable_defaults: bool,
) -> Result<Vec<ResolvedFilter>, FilterResolutionError> {
    let twice = duplicated(enabled);
    if !twice.is_empty() {
        return Err(FilterResolutionError::EnabledTwice { names: twice });
    }
    let twice = duplicated(disabled);
    if !twice.is_empty() {
        return Err(FilterResolutionError::DisabledTwice { names: twice });
    }

    // `getAllUserEnabledReadFilterNames` is the enabled ones AND the inverted ones.
    let all_enabled: Vec<String> = enabled.iter().chain(inverted.iter()).cloned().collect();
    let both: Vec<String> = all_enabled
        .iter()
        .filter(|name| disabled.contains(name))
        .cloned()
        .collect();
    if !both.is_empty() {
        return Err(FilterResolutionError::EnabledAndDisabled { names: both });
    }
    if enabled.iter().any(|name| inverted.contains(name)) {
        return Err(FilterResolutionError::EnabledAndInverted);
    }

    for name in disabled {
        if !known.contains(&name.as_str()) {
            return Err(FilterResolutionError::DisabledDoesNotExist { name: name.clone() });
        }
        // A filter the tool never enabled is only a warning, which is not an observable here.
    }

    if !disable_defaults {
        let redundant: Vec<String> = defaults
            .iter()
            .filter(|name| inverted.contains(&(*name).to_string()))
            .map(|name| (*name).to_string())
            .collect();
        if !redundant.is_empty() {
            return Err(FilterResolutionError::InvertedIsADefault { names: redundant });
        }
    }

    for name in &all_enabled {
        if !known.contains(&name.as_str()) && !defaults.contains(&name.as_str()) {
            return Err(FilterResolutionError::Unrecognized { name: name.clone() });
        }
    }

    let mut resolved: Vec<ResolvedFilter> = if disable_defaults {
        Vec::new()
    } else {
        defaults
            .iter()
            .filter(|name| !disabled.contains(&(*name).to_string()))
            .map(|name| ResolvedFilter {
                name: (*name).to_string(),
                negated: false,
            })
            .collect()
    };

    for name in enabled {
        let candidate = ResolvedFilter {
            name: name.clone(),
            negated: false,
        };
        if !resolved.contains(&candidate) {
            resolved.push(candidate);
        }
    }
    for name in inverted {
        // A negated filter is never equal to the filter it negates, so it is always appended.
        let candidate = ResolvedFilter {
            name: name.clone(),
            negated: true,
        };
        if !resolved.contains(&candidate) {
            resolved.push(candidate);
        }
    }
    Ok(resolved)
}
