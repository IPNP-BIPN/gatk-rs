//! Ported from `org.broadinstitute.hellbender.utils.variant.VariantContextGetters`
//! (GATK 4.6.2.0): how an INFO attribute becomes a `double[]`.
//!
//! Everything that reads a numeric annotation off a `VariantContext` goes through here, including
//! `Mutect2FilteringEngine.getTumorLogOdds`, and the function is held together by an accident.
//!
//! # The `Iterable` branch is unreachable
//!
//! ```java
//! } else if (value.getClass().isAssignableFrom(Iterable.class)) {
//! ```
//!
//! The test is **backwards**. `isAssignableFrom` asks whether the argument's type can be assigned
//! to the receiver's, so this asks whether an `Iterable` can be stored in an `ArrayList` variable,
//! which is false. Every list therefore falls through to the branch labelled "as a last resort":
//!
//! ```java
//! return Stream.of(String.valueOf(value).trim().replaceAll("\\[|\\]", "").split(","))
//! ```
//!
//! `String.valueOf(List.of(1.0, 2.0))` is `"[1.0, 2.0]"`, so after the brackets are stripped and
//! the string is split on the comma the second field is `" 2.0"`, **with a leading space**. That
//! parses only because `Double.parseDouble` trims. A port that split on `", "`, or that used a
//! parser refusing leading whitespace, would agree on a one-element list and diverge on every
//! longer one.
//!
//! # `.` is not a parse failure
//!
//! A field equal to `VCFConstants.MISSING_VALUE_v4` becomes the caller's `missingValue`, which for
//! the `getAttributeAsDoubleArray(vc, key)` overload is `-1`. So a `TLOD` of `.` reaches
//! `getTumorLogOdds` as `-1` and then, after the log conversion, as `-2.302...`; nothing reports a
//! missing value as missing.
//!
//! # Anything else is a `GATKException`, not a `NumberFormatException`
//!
//! The parse failure is caught and rethrown with the key in the message, so a caller sees a GATK
//! exception naming the annotation rather than Java's.

use htsjdk_vcf::genotype_likelihoods::parse_java_double;
use htsjdk_vcf::variant::{Value, VariantContext};

/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";

/// The `GATKException` the converter throws for a field that is neither `.` nor a double.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NonDoubleValue {
    pub key: String,
    pub value: String,
}

/// A value this port will not render, as distinct from one the reference refuses.
///
/// The last-resort branch renders its input with `String.valueOf`, so an attribute already held as
/// a `Double` goes through `Double.toString`, which is its own algorithm and is not ported. Rather
/// than emit a plausible rendering and compare it to the oracle's, the port says so. Nothing that
/// comes off a parsed VCF reaches this: htsjdk's decoder keeps INFO values as strings, so the
/// `Double` case only arises for a context built in memory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnrenderableDouble;

impl NonDoubleValue {
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.GATKException"
    }

    pub fn message(&self) -> String {
        format!(
            "INFO annotation '{}' contains a non-double value '{}'",
            self.key, self.value
        )
    }
}

/// The `ToDoubleFunction` inside `attributeValueToDoubleArray`.
fn to_double(text: &str, key: &str, missing_value: f64) -> Result<f64, NonDoubleValue> {
    if text == MISSING_VALUE {
        return Ok(missing_value);
    }
    // `Double.parseDouble`, which trims, takes type suffixes and hexadecimal floats, and refuses
    // the `inf`/`nan` spellings Rust's own parser takes.
    parse_java_double(text).ok_or_else(|| NonDoubleValue {
        key: key.to_string(),
        value: text.to_string(),
    })
}

/// `String.valueOf` of an attribute, as the last-resort branch sees it.
///
/// A list renders as `AbstractCollection.toString`, which is `[a, b]`: square brackets and `", "`
/// between the elements. That rendering is what the branch then takes apart, which is why the
/// separator it splits on leaves a space at the front of every field but the first.
fn java_string_of(value: &Value) -> Result<String, UnrenderableDouble> {
    Ok(match value {
        // `String.valueOf(null)` is the string "null", which then fails to parse; the reference
        // reaches that through an attribute explicitly set to a null element.
        Value::Missing => MISSING_VALUE.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Int(number) => number.to_string(),
        // `Double.toString` is not ported; see `UnrenderableDouble`.
        Value::Double(_) => return Err(UnrenderableDouble),
        Value::Str(text) => text.clone(),
        Value::List(values) => {
            let parts: Result<Vec<String>, UnrenderableDouble> =
                values.iter().map(java_string_of).collect();
            format!("[{}]", parts?.join(", "))
        }
    })
}

/// `attributeValueToDoubleArray`.
///
/// `None` for an attribute that is not there, which is the `defaultResult.get()` of the
/// `getAttributeAsDoubleArray(vc, key)` overload, whose supplier returns null.
pub fn attribute_value_to_double_array(
    value: Option<&Value>,
    key: &str,
    missing_value: f64,
) -> Result<Option<Vec<f64>>, NonDoubleValue> {
    let Some(value) = value else {
        return Ok(None);
    };

    // The Iterable branch is unreachable, so a list takes the same path as a scalar: rendered with
    // String.valueOf, stripped of brackets, split on the comma.
    let rendered =
        java_string_of(value).expect("no Double-valued attribute; see UnrenderableDouble");
    let trimmed = rendered.trim_matches(|c: char| c <= ' ');
    let stripped: String = trimmed.chars().filter(|c| *c != '[' && *c != ']').collect();

    let mut out = Vec::new();
    for field in java_split_on_comma(&stripped) {
        out.push(to_double(field, key, missing_value)?);
    }
    Ok(Some(out))
}

/// `String.split(",")`: trailing empty fields dropped, and a string the delimiter never matches
/// returned whole as a single field.
fn java_split_on_comma(text: &str) -> Vec<&str> {
    if !text.contains(',') {
        return vec![text];
    }
    let mut parts: Vec<&str> = text.split(',').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `getAttributeAsDoubleArray(vc, key)`: the overload with a null default and a missing value of
/// `-1`, which is the one every caller in the annotator package uses.
pub fn get_attribute_as_double_array(
    vc: &VariantContext,
    key: &str,
) -> Result<Option<Vec<f64>>, NonDoubleValue> {
    let value = vc
        .attributes
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value);
    attribute_value_to_double_array(value, key, -1.0)
}

/// `MathUtils.maxElementIndex(array)`: the index of the largest element, ties going to the first.
///
/// A `NaN` loses every comparison, so an array whose first element is `NaN` reports index 0 unless
/// something later beats it, and an array that is all `NaN` reports 0.
pub fn max_element_index(values: &[f64]) -> Option<usize> {
    if values.is_empty() {
        return None;
    }
    let mut best = 0usize;
    for (index, value) in values.iter().enumerate().skip(1) {
        if *value > values[best] {
            best = index;
        }
    }
    Some(best)
}

/// `Mutect2FilteringEngine.getTumorLogOdds(vc)`: the `TLOD` attribute, converted from log10 to
/// natural log **in place**.
///
/// The conversion is applied to whatever the getter produced, including the `-1` a `.` field
/// becomes, so a missing element arrives as `-2.302...` and not as anything a caller can
/// recognise as missing.
pub fn get_tumor_log_odds(vc: &VariantContext) -> Result<Option<Vec<f64>>, NonDoubleValue> {
    Ok(
        get_attribute_as_double_array(vc, TUMOR_LOG_10_ODDS_KEY)?.map(|values| {
            values
                .into_iter()
                .map(crate::allele_likelihoods::log10_to_log)
                .collect()
        }),
    )
}

/// `GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY`.
pub const TUMOR_LOG_10_ODDS_KEY: &str = "TLOD";
