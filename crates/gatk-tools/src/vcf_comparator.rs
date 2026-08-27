//! `VCFComparator`: what counts as two VCFs disagreeing.
//!
//! The tool's output is its exceptions. It walks two files in step and throws on the first
//! difference it is not told to tolerate, so what is ported is which differences it notices, what
//! it says about each, and which argument silences it.
//!
//! Reading the VCFs and the overlap grouping are not ported. The comparisons are.

/// One variant, reduced to what the comparisons read.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    /// `actual` or `expected`, which is the file's TAG rather than its order.
    pub source: String,
    pub contig: String,
    pub start: i32,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    pub qual: Option<f64>,
    /// Empty means no filter was applied at all, which is not the same as `PASS`.
    pub filters: Vec<String>,
    pub attributes: Vec<(String, String)>,
    pub genotype_qualities: Vec<i32>,
}

impl Variant {
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// `PASS` is a filter that was applied and passed; an empty list is one that never ran.
    pub fn filters_were_applied(&self) -> bool {
        !self.filters.is_empty()
    }

    pub fn has_gq_zero(&self) -> bool {
        self.genotype_qualities.contains(&0)
    }
}

/// The tolerances, with the defaults the tool declares.
#[derive(Debug, Clone, PartialEq)]
pub struct Tolerances {
    pub qual_change_allowed: f64,
    pub ignore_quals: bool,
    pub ignore_filters: bool,
    pub ignore_annotations: bool,
    pub ignored_attributes: Vec<String>,
    pub ignore_dbsnp: bool,
    pub ignore_gq0: bool,
    pub allow_extra_alleles: bool,
    pub allow_new_stars: bool,
    pub allow_missing_stars: bool,
    pub positions_only: bool,
}

impl Default for Tolerances {
    fn default() -> Self {
        Tolerances {
            qual_change_allowed: DEFAULT_QUAL_CHANGE_ALLOWED,
            ignore_quals: false,
            ignore_filters: false,
            ignore_annotations: false,
            ignored_attributes: Vec::new(),
            ignore_dbsnp: false,
            ignore_gq0: false,
            allow_extra_alleles: false,
            allow_new_stars: false,
            allow_missing_stars: false,
            positions_only: false,
        }
    }
}

/// The QUAL tolerance's default, which is not zero.
pub const DEFAULT_QUAL_CHANGE_ALLOWED: f64 = 0.001;

/// One complaint.
#[derive(Debug, Clone, PartialEq)]
pub enum Complaint {
    /// Not wrapped with a position: this one names the file and the position itself.
    UnmatchedVariant {
        source: String,
        contig: String,
        start: i32,
    },
    QualDiffers {
        difference: f64,
        tolerance: f64,
    },
    FiltersNotApplied,
    FiltersDiffer {
        expected: String,
        actual: String,
    },
    AllelesMismatched {
        actual: Vec<String>,
        expected: Vec<String>,
    },
    AttributesDiffer {
        key: String,
        actual: String,
        expected: String,
    },
    DbsnpIdsDiffer,
}

impl Complaint {
    /// The bare message, before the position is wrapped around it.
    pub fn message(&self) -> String {
        match self {
            Complaint::UnmatchedVariant {
                source,
                contig,
                start,
            } => {
                format!("Unmatched variant in {source} at position {contig}:{start}")
            }
            Complaint::QualDiffers {
                difference,
                tolerance,
            } => format!(
                "qual scores differ by {}, which is more than {tolerance}",
                gatk_engine::tsv_table::java_double_to_string(*difference)
            ),
            Complaint::FiltersNotApplied => {
                " filters were not applied to both variants".to_string()
            }
            Complaint::FiltersDiffer { expected, actual } => format!(
                "variants have different filters: expected has {expected} and actual has {actual}"
            ),
            Complaint::AllelesMismatched { actual, expected } => format!(
                "Alleles are mismatched at {{position}}: actual has [{}] and expected has [{}]",
                actual.join(", "),
                expected.join(", ")
            ),
            Complaint::AttributesDiffer {
                key,
                actual,
                expected,
            } => format!(
                "Variant contexts have different attribute values for {key}: actual has {actual} \
                 and expected has {expected}"
            ),
            Complaint::DbsnpIdsDiffer => "dbsnp IDs differ for VCs".to_string(),
        }
    }

    /// `wrapWithPosition`, which puts the position in FRONT of the message rather than inside it.
    /// The unmatched-variant complaint is not wrapped: it carries its own position already.
    pub fn wrapped(&self, contig: &str, start: i32) -> String {
        match self {
            Complaint::UnmatchedVariant { .. } => self.message(),
            Complaint::AllelesMismatched { actual, expected } => format!(
                "At position {contig}:{start} Alleles are mismatched at {contig}:{start}: actual \
                 has [{}] and expected has [{}]",
                actual.join(", "),
                expected.join(", ")
            ),
            _ => format!("At position {contig}:{start} {}", self.message()),
        }
    }
}

/// `actualHasNewAlleles(a, b)`: does `a` carry an alternate that `b` does not?
pub fn has_new_alleles(a: &Variant, b: &Variant) -> bool {
    a.alternates
        .iter()
        .any(|allele| !b.alternates.contains(allele))
}

/// The unmatched-variant test, which is GUARDED ON A GENOTYPE QUALITY OF ZERO: a variant present on
/// one side alone and confidently called passes in silence.
pub fn unmatched(variant: &Variant, tolerances: &Tolerances) -> Option<Complaint> {
    if tolerances.ignore_gq0 || !variant.has_gq_zero() {
        return None;
    }
    Some(Complaint::UnmatchedVariant {
        source: variant.source.clone(),
        contig: variant.contig.clone(),
        start: variant.start,
    })
}

/// The whole comparison of one matched pair, in the order the tool makes it.
///
/// The ATTRIBUTE comparison comes before the allele one, which is why an extra allele's effect on
/// `AC` hides the allele question entirely until that key is ignored.
///
/// The allele check is guarded by `actualHasNewAlleles(expected, actual)`, with its arguments the
/// OTHER WAY ROUND from the check it guards: it asks whether EXPECTED has an allele actual lacks.
/// An allele ADDED to actual is therefore never checked, and `--allow-extra-alleles` exists for a
/// direction this guard never reaches.
pub fn compare(actual: &Variant, expected: &Variant, tolerances: &Tolerances) -> Option<Complaint> {
    if tolerances.positions_only {
        return None;
    }

    if !tolerances.ignore_quals {
        if let (Some(a), Some(e)) = (actual.qual, expected.qual) {
            let difference = (a - e).abs();
            if difference > tolerances.qual_change_allowed {
                return Some(Complaint::QualDiffers {
                    difference,
                    tolerance: tolerances.qual_change_allowed,
                });
            }
        }
    }

    if !tolerances.ignore_filters {
        if actual.filters_were_applied() != expected.filters_were_applied() {
            return Some(Complaint::FiltersNotApplied);
        }
        if actual.filters != expected.filters {
            return Some(Complaint::FiltersDiffer {
                expected: expected.filters.join(", "),
                actual: format!("[{}]", actual.filters.join(", ")),
            });
        }
    }

    if !tolerances.ignore_annotations {
        for (key, expected_value) in &expected.attributes {
            if tolerances.ignored_attributes.contains(key) {
                continue;
            }
            let actual_value = actual.attribute(key).unwrap_or("");
            if actual_value != expected_value {
                return Some(Complaint::AttributesDiffer {
                    key: key.clone(),
                    actual: actual_value.to_string(),
                    expected: expected_value.clone(),
                });
            }
        }
    }

    // The reversed guard, and behind it the check itself.
    if has_new_alleles(expected, actual) {
        if let Some(complaint) = check_alleles(actual, expected, tolerances) {
            return Some(complaint);
        }
    }

    if !tolerances.ignore_dbsnp
        && actual.id != expected.id
        && actual.alternates.len() == expected.alternates.len()
    {
        return Some(Complaint::DbsnpIdsDiffer);
    }

    None
}

/// The spanning-deletion allele.
pub const SPAN_DEL: &str = "*";

fn has_new_star(actual: &Variant, expected: &Variant) -> bool {
    actual.alternates.iter().any(|a| a == SPAN_DEL)
        && !expected.alternates.iter().any(|a| a == SPAN_DEL)
}

fn has_missing_star(actual: &Variant, expected: &Variant) -> bool {
    expected.alternates.iter().any(|a| a == SPAN_DEL)
        && !actual.alternates.iter().any(|a| a == SPAN_DEL)
}

/// `checkAlleles`.
///
/// A chain of three branches with a CATCH-ALL `else` that throws unconditionally, which is what
/// makes the flags weaker than they look: `--allow-extra-alleles` suppresses only the FIRST
/// branch, and its condition is the one the guard in front has already ruled out for a missing
/// allele. Such a pair therefore falls straight through to the catch-all and is refused with the
/// flag set.
///
/// The two star branches are ported but not measured: the fixture carries no spanning deletion.
pub fn check_alleles(
    actual: &Variant,
    expected: &Variant,
    tolerances: &Tolerances,
) -> Option<Complaint> {
    let mismatch = || {
        Some(Complaint::AllelesMismatched {
            actual: actual.alternates.clone(),
            expected: expected.alternates.clone(),
        })
    };
    // The first two branches produce the same complaint from different causes, which is why they
    // are written as one condition here and as two in the reference.
    if (!tolerances.allow_extra_alleles && has_new_alleles(actual, expected))
        || (!tolerances.allow_new_stars && has_new_star(actual, expected))
    {
        mismatch()
    } else if !tolerances.allow_missing_stars && has_missing_star(actual, expected) {
        // The reference lets this one through when the remainder is exactly the star.
        let remainder: Vec<&String> = expected
            .alternates
            .iter()
            .filter(|allele| !actual.alternates.contains(allele))
            .collect();
        if remainder.len() > 1 || !remainder.iter().any(|allele| *allele == SPAN_DEL) {
            mismatch()
        } else {
            None
        }
    } else {
        // The catch-all, which throws whatever was or was not allowed.
        mismatch()
    }
}

/// What the tool refuses about its inputs, before a record is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputError {
    WrongNumberOfInputs,
    NoExpectedInput,
    MutuallyExclusive { argument: String, other: String },
}

impl InputError {
    pub fn message(&self) -> String {
        match self {
            InputError::WrongNumberOfInputs => {
                "Bad input: VCFComparator expects exactly two inputs -- one actual and one \
                 expected."
                    .to_string()
            }
            InputError::NoExpectedInput => {
                "Bad input: Tool requires exactly one expected input file".to_string()
            }
            InputError::MutuallyExclusive { argument, other } => format!(
                "Argument '{argument}' cannot be used in conjunction with argument(s) {other}"
            ),
        }
    }
}

/// `onTraversalStart`'s two checks, in its own order: the COUNT first, then the tag.
pub fn check_inputs(tags: &[String]) -> Result<(), InputError> {
    if tags.len() != 2 {
        return Err(InputError::WrongNumberOfInputs);
    }
    if tags.iter().filter(|tag| *tag == "expected").count() != 1 {
        return Err(InputError::NoExpectedInput);
    }
    Ok(())
}
