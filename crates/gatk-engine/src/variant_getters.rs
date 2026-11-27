//! Ported from `org.broadinstitute.hellbender.utils.variant.VariantContextGetters`
//! (GATK 4.7.0.0; `VariantContextGetters.java` is unchanged since 4.6.2.0): how an INFO attribute
//! becomes a `double[]`.
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
//! return Stream.of(String.valueOf(value).trim().replaceAll("\\[|\\]", "")
//!         .split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR))
//!         .mapToDouble(doubleConverter).toArray();
//! ```
//!
//! `INFO_FIELD_ARRAY_SEPARATOR` is `","` in htsjdk 5.0.0's `VCFConstants`.
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
//! missing value as missing. That holds for a `.` that is the whole field; a `.` after the first
//! element of a list IS a parse failure, as the next section's `list-with-missing` row shows.
//!
//! # Anything else is a `GATKException`, not a `NumberFormatException`
//!
//! The parse failure is caught and rethrown with the key in the message, so a caller sees a GATK
//! exception naming the annotation rather than Java's.
//!
//! # What comes back, by input
//!
//! Rows of `tests/data/variant_getters.txt.gz`:
//!
//! - attribute absent (`absent`): `None`, the overload's null default;
//! - a field exactly equal to `.` (`one-missing`, the string `1.5,.,2.5`): the caller's
//!   `missingValue`, `-1`, silently;
//! - a `.` inside a LIST attribute (`list-with-missing`, `["1.5", ".", "2.5"]`): `Err`. The list
//!   renders as `[1.5, ., 2.5]`, the field is `" ."` with the space `toString` put there, and the
//!   missing test is `s.equals(...)` on the untrimmed field, so it fails and the parse throws.
//!   A multi-valued INFO field read from a VCF is such a list (htsjdk's `parseInfo` splits it),
//!   so `TLOD=3.2,.` throws where `TLOD=.` gives `-1`;
//! - a field `Double.parseDouble` refuses (`not-a-double`, `empty-string`): `Err`, the
//!   `GATKException`;
//! - a trailing empty field (`trailing-comma`, `1.5,2.5,`): dropped by `String.split` before it is
//!   parsed, so two values.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a `T` or nothing", Java's nullable reference; `Result<T, E>` is "a `T` or an
//!   error `E`", how the `GATKException` is modelled; `?` returns the error to the caller at once;
//! - `let Some(x) = y else { ... }` unwraps `y` or takes the early exit in the braces;
//! - `.iter().map(..).collect()` is Java's `stream().map(..).collect(..)`, and `|c| ...` is a
//!   lambda;
//! - `matches!(item, Value::Double(_))` is an `instanceof` test on one case of an enum.

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
///
/// ```java
/// final String s = String.valueOf(o);
/// if (s.equals(VCFConstants.MISSING_VALUE_v4)) {
///     return missingValue;
/// ```
///
/// The comparison is exact and on the field as it stands, leading space included; only the
/// parse after it trims. The reference's first branch, `o == null` giving `missingValue`, has no
/// counterpart because every field reaching here is a substring of a rendered string.
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
fn java_string_of(value: &Value) -> String {
    match value {
        // `Value::Missing` is htsjdk's `VCFConstants.MISSING_VALUE_v4` object, the String ".",
        // not a Java null: htsjdk-rs's `record_parse.rs` assigns it to a bare non-Flag key and to
        // `key=`. So it renders as `.`, and the converter then maps it to the missing value.
        Value::Missing => MISSING_VALUE.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Int(number) => number.to_string(),
        // `Double.toString`, as the pinned Java 17 prints it (htsjdk-rs decision 0046). Nothing
        // off a parsed VCF reaches this, since htsjdk's decoder keeps INFO values as strings; a
        // context built in memory does.
        Value::Double(number) => jmath::to_string::double_to_string(*number),
        Value::Str(text) => text.clone(),
        Value::List(values) => format!(
            "[{}]",
            values
                .iter()
                .map(java_string_of)
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// `attributeValueToDoubleArray`.
///
/// `None` for an attribute that is not there, which is the `defaultResult.get()` of the
/// `getAttributeAsDoubleArray(vc, key)` overload, whose supplier returns null.
///
/// ```java
/// } else if (value.getClass().isArray()) {
///     final double[] result = new double[Array.getLength(value)];
///     for (int i = 0; i < result.length; i++) {
///         result[i] = doubleConverter.applyAsDouble(String.valueOf(Array.get(value, i)));
///     }
/// ```
///
/// `Value` has no separate array case, so a non-empty list whose items are all doubles stands in
/// for the `double[]`. An empty list is not taken as an array: it renders as `[]`, which becomes
/// one empty field and an `Err`, where an empty Java `double[]` would give an empty array.
pub fn attribute_value_to_double_array(
    value: Option<&Value>,
    key: &str,
    missing_value: f64,
) -> Result<Option<Vec<f64>>, NonDoubleValue> {
    let Some(value) = value else {
        return Ok(None);
    };

    // An in-memory `double[]` (Mutect2's `TLOD`, built with `.toArray()`) takes the array branch,
    // `Double.parseDouble(String.valueOf(d))` per element, which gives back `d` exactly: Java's
    // `Double.toString` always round-trips. So the doubles are the answer and no text is needed.
    if let Value::List(items) = value {
        if !items.is_empty() && items.iter().all(|item| matches!(item, Value::Double(_))) {
            return Ok(Some(
                items
                    .iter()
                    // `match` must cover every case of `Value`; the `_` arm is the one the `all`
                    // test above has ruled out, and `unreachable!` panics if that ever breaks.
                    .map(|item| match item {
                        Value::Double(d) => *d,
                        _ => unreachable!("every item is a double"),
                    })
                    .collect(),
            ));
        }
    }

    // The Iterable branch is unreachable, so a list takes the same path as a scalar: rendered with
    // String.valueOf, stripped of brackets, split on the comma.
    let rendered = java_string_of(value);
    // Java's `String.trim()`, which cuts every char up to and including the space (control
    // characters too); Rust's own `trim` cuts Unicode whitespace, a different set.
    let trimmed = rendered.trim_matches(|c: char| c <= ' ');
    // `replaceAll("\\[|\\]", "")`: every bracket anywhere, not only the outer pair.
    let stripped: String = trimmed.chars().filter(|c| *c != '[' && *c != ']').collect();

    let mut out = Vec::new();
    for field in java_split_on_comma(&stripped) {
        out.push(to_double(field, key, missing_value)?);
    }
    Ok(Some(out))
}

/// `String.split(",")`: trailing empty fields dropped, and a string the delimiter never matches
/// returned whole as a single field.
///
/// Rust's `split` keeps trailing empty fields, hence the loop that pops them. The two edge cases
/// differ: `""` is one empty field (no match, returned whole), which then fails to parse, while
/// `","` is no fields at all, so a string attribute `,` gives an empty array.
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
///
/// ```java
/// return getAttributeAsDoubleArray(variantContext, attribute, () -> null, -1);
/// ```
///
/// The attributes are a list of pairs here, and `.find` takes the first whose name matches,
/// where the reference looks the key up in a map.
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
/// A `NaN` loses every comparison, in both directions: `x > NaN` is false for every `x`, so a
/// `NaN` at index 0 is never replaced and reports 0 whatever follows (golden row `with-nan`,
/// `[NaN, 1.0]`, is `0`), and a `NaN` later is never chosen (`nan-last`, `[1.0, NaN]`, is `0`).
///
/// ```java
/// Utils.validateArg(array.length > 0, "array may not be empty");
/// ```
///
/// An empty array is that `IllegalArgumentException` in the reference and `None` here.
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
///
/// The reference reads the attribute with
/// `VariantContextGetters.getAttributeAsDoubleArray(vc, GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY)`
/// and returns `tumorLog10Odds == null ? null :` the array passed through
/// `MathUtils.applyToArrayInPlace(tumorLog10Odds, MathUtils::log10ToLog)`, where `log10ToLog`
/// is `log10 * LOG_10`. `Option::map` is the null test: it converts only a `Some`.
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
