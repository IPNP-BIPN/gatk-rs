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
//!
//! # Reading this file without reading Rust
//!
//! - `&[String]` and `&[&str]` are borrowed, read-only views of a list of strings; the caller keeps
//!   the list. `&str` is a borrowed string, `String` an owned one, and `.to_string()` or
//!   `.clone()` makes an owned copy.
//! - `Result<T, E>` is "a value or an error": `Ok(list)` is the resolved list, `Err(..)` is the
//!   `CommandLineException` the descriptor would throw.
//! - an `enum` variant with `{ names }` carries data, like a small final class per exception case;
//!   `match` must cover every variant, and `_` is "any other".
//! - `.iter().filter(..).map(..).collect()` is Java's `stream().filter().map().collect()`.

/// One filter of the resolved list: a `ReadFilter` instance of `getResolvedInstances`, identified
/// by its name and whether it is the `ReadFilterNegate` wrapper.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedFilter {
    /// The filter's simple class name, as given on the command line.
    pub name: String,
    /// Whether it arrived through `--inverted-read-filter`, which wraps it in a `ReadFilterNegate`.
    pub negated: bool,
}

/// What `validateAndResolvePlugins` refuses, one variant per `throw`, in the order the Java
/// checks them.
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
    /// The class of the exception the Java throws: `CommandLineException.BadArgumentValue` for the
    /// two duplicate checks and the unknown disabled filter, plain `CommandLineException` for the
    /// rest.
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

    /// The exception's message. `BadArgumentValue` prefixes `Illegal argument value: ` to the
    /// formatted text, which is why only its three variants carry the prefix here.
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

/// `Utils.getDuplicatedItems`: the names that appear more than once, in first-seen order.
///
/// The Java collects them into a `HashSet`, so it joins them in hash order, not first-seen
/// order. The two differ only when two or more names are duplicated.
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
///
/// The checks run in the Java's order, so a command line with two faults reports the first one
/// the Java reaches. The resolved list is then the defaults minus the disabled ones, followed by
/// the enabled and the inverted filters not already present:
///
/// ```java
/// userEnabledFilters.stream()
///         .filter(f -> !finalFilters.contains(f)) // remove redundant filters
///         .forEach(f -> finalFilters.add(f));
/// ```
///
/// The Java's warnings (a disabled filter the tool never enabled, an enabled filter that is
/// already a default) change nothing in the output, so they are not reproduced.
pub fn resolve(
    defaults: &[&str],
    known: &[&str],
    enabled: &[String],
    disabled: &[String],
    inverted: &[String],
    disable_defaults: bool,
) -> Result<Vec<ResolvedFilter>, FilterResolutionError> {
    // `return Err(..)` is the `throw`: the function stops and the caller gets the error.
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
    // `|name| ...` is a closure, Java's lambda `name -> ...`.
    if enabled.iter().any(|name| inverted.contains(name)) {
        return Err(FilterResolutionError::EnabledAndInverted);
    }

    // `allDiscoveredReadFilters.containsKey(s)`. `as_str()` borrows the `String` as a `&str` so
    // it can be compared with the `&str` entries of `known`.
    for name in disabled {
        if !known.contains(&name.as_str()) {
            return Err(FilterResolutionError::DisabledDoesNotExist { name: name.clone() });
        }
        // A filter the tool never enabled is only a warning, which is not an observable here.
    }

    // `redundantInv`, checked only when the defaults are kept. The Java's set is a `HashSet`, so
    // with two or more names its order is hash order; this lists them in the tool's order.
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

    // `Unrecognized read filter name`: a name neither in the library nor among the tool's
    // defaults.
    for name in &all_enabled {
        if !known.contains(&name.as_str()) && !defaults.contains(&name.as_str()) {
            return Err(FilterResolutionError::Unrecognized { name: name.clone() });
        }
    }

    // `if` is an expression in Rust: the list is whichever branch runs, like Java's `?:` here.
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

    // `getUserEnabledInstances` adds the enabled filters first, then the inverted ones.
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
        // A negated filter is never equal to the filter it negates, so it is appended even when
        // the plain filter is already in the list. The same name inverted twice is appended once
        // here, but twice by the Java: `negate()` answers a new `ReadFilterNegate` each time and
        // nothing checks the inverted list for duplicates.
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
