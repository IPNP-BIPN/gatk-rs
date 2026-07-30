//! `RawGtCount`, ported from `org.broadinstitute.hellbender.tools.walkers.annotator.RawGtCount`.
//!
//! This is the odd member of the counting family: its `annotate` returns **null**, not an empty
//! map. It is a `ReducibleAnnotation`, so the only path that produces anything is `combineRawData`,
//! which merges the `RAW_GT_COUNT` strings of several GVCFs.
//!
//! Three things in the merge are not what the key suggests.
//!
//! The parser is a regex split on `", *"` after stripping square brackets:
//!
//! ```java
//! final String[] parsed = rawDataString.trim().replaceAll(BRACKET_REGEX, "").split(", *");
//! ```
//!
//! A comma followed by spaces splits; a **space followed by a comma** does not, so `"1 , 2, 3"`
//! parses as two fields and is refused for having the wrong arity rather than for the space. The
//! brackets are stripped anywhere in the string, not just at the ends, because
//! `BRACKET_REGEX = "\\[|\\]"` and `replaceAll` is global.
//!
//! The written value drops one of the three counts it just summed:
//!
//! ```java
//! return "." + SEPARATOR + perAlleleData.get(Allele.NO_CALL).get(1) + SEPARATOR + perAlleleData.get(Allele.NO_CALL).get(2);
//! ```
//!
//! The hom-ref total is replaced by a literal `.`, so a round trip through this annotation is
//! lossy: combining `[1, 2, 3]` with itself gives `.,4,6` and not `2,4,6`. The combined value is
//! still declared `Integer` with `Number=3` in the header, so the `.` sits in an integer field.
//!
//! And the summing is under `Allele.NO_CALL` rather than under any real allele, because the counts
//! are per-site and the reducible machinery is per-allele. That key never reaches a file; it is the
//! map slot the accumulator happens to use.

use crate::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_engine::context::ReferenceContext;
use htsjdk_vcf::variant::VariantContext;

/// `GATKVCFConstants.RAW_GENOTYPE_COUNT_KEY`.
pub const RAW_GENOTYPE_COUNT_KEY: &str = "RAW_GT_COUNT";

const SEPARATOR: &str = ",";

/// What `parseRawDataString` refuses, and how it says so. Both messages are `UserException.BadInput`
/// and differ only in their text, which is what a user sees.
#[derive(Debug, Clone, PartialEq)]
pub enum RawGtCountError {
    /// The field did not have exactly three values.
    WrongArity { found: usize, raw: String },
    /// One of the three was not an integer.
    NotAnInteger { raw: String },
}

impl RawGtCountError {
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }

    pub fn message(&self) -> String {
        match self {
            RawGtCountError::WrongArity { found, raw } => format!(
                "A USER ERROR has occurred: Bad input: Raw value for {RAW_GENOTYPE_COUNT_KEY} has {found} values, expected 3. Annotation value is {raw}"
            ),
            RawGtCountError::NotAnInteger { raw } => format!(
                "A USER ERROR has occurred: Bad input: malformed {RAW_GENOTYPE_COUNT_KEY} annotation: {raw}"
            ),
        }
    }
}

/// `String.trim()`: strips every character at or below `U+0020` from both ends, which is not the
/// same set as Rust's `trim`, since it drops no other Unicode whitespace and does drop control
/// characters.
fn java_trim(text: &str) -> &str {
    let bytes = text.as_bytes();
    let mut start = 0;
    let mut end = bytes.len();
    while start < end && bytes[start] <= b' ' {
        start += 1;
    }
    while end > start && bytes[end - 1] <= b' ' {
        end -= 1;
    }
    // Both ends moved over bytes that are ASCII, so the slice is still on character boundaries.
    &text[start..end]
}

/// `split(", *")`: a comma followed by any number of spaces. Java drops trailing empty fields.
fn split_on_comma_spaces(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let mut parts = Vec::new();
    let mut field_start = 0;
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] == b',' {
            parts.push(&text[field_start..index]);
            index += 1;
            while index < bytes.len() && bytes[index] == b' ' {
                index += 1;
            }
            field_start = index;
        } else {
            index += 1;
        }
    }
    parts.push(&text[field_start..]);
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `parseRawDataString`: strip brackets anywhere, split, and parse three integers.
pub fn parse_raw_data_string(raw: &str) -> Result<[i32; 3], RawGtCountError> {
    let trimmed = java_trim(raw);
    let stripped: String = trimmed.chars().filter(|c| *c != '[' && *c != ']').collect();
    let parsed = split_on_comma_spaces(&stripped);
    if parsed.len() != 3 {
        return Err(RawGtCountError::WrongArity {
            found: parsed.len(),
            raw: raw.to_string(),
        });
    }
    let mut counts = [0i32; 3];
    for (slot, field) in counts.iter_mut().zip(&parsed) {
        // `Integer.parseInt` takes a leading `+` and refuses everything Rust's parser refuses plus
        // the underscore separators Rust also refuses, so the two agree on this alphabet.
        *slot = field
            .parse::<i32>()
            .map_err(|_| RawGtCountError::NotAnInteger {
                raw: raw.to_string(),
            })?;
    }
    Ok(counts)
}

/// `combineRawData`: sum the three counts across the inputs and render the result.
///
/// The rendering discards the hom-ref total, writing `.` in its place, so this is not associative
/// with itself in the way the field name implies.
pub fn combine_raw_data(raw_values: &[String]) -> Result<String, RawGtCountError> {
    let mut combined: Option<[i32; 3]> = None;
    for raw in raw_values {
        let counts = parse_raw_data_string(raw)?;
        combined = Some(match combined {
            None => counts,
            Some(sum) => [sum[0] + counts[0], sum[1] + counts[1], sum[2] + counts[2]],
        });
    }
    // With no inputs at all the reference dereferences a null and throws; a caller reaching that
    // has nothing to combine, so it is reported as the arity failure rather than reproduced.
    let sum = combined.ok_or(RawGtCountError::WrongArity {
        found: 0,
        raw: String::new(),
    })?;
    Ok(format!(".{SEPARATOR}{}{SEPARATOR}{}", sum[1], sum[2]))
}

pub struct RawGtCount;

impl InfoFieldAnnotation for RawGtCount {
    fn key_names(&self) -> Vec<&'static str> {
        vec![RAW_GENOTYPE_COUNT_KEY]
    }

    /// The reference returns `null` here, not an empty map. Every consumer of this annotation goes
    /// through [`combine_raw_data`] instead, and a caller that treated the null as a map would have
    /// thrown rather than written nothing.
    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        _vc: &VariantContext,
    ) -> Vec<(String, AnnotationValue)> {
        Vec::new()
    }
}
