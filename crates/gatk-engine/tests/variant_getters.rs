//! Conformance for how an INFO attribute becomes a `double[]`, against the oracle.
//!
//! Golden from `tools/genotyper-conformance/VariantGettersDump.java`.
//!
//! The row that pays for the whole suite:
//!
//! ```text
//! array  list-with-missing  E:...GATKException:INFO annotation 'X' contains a non-double value ' .'
//! ```
//!
//! The list `["1.5", ".", "2.5"]` renders as `[1.5, ., 2.5]`, so after the split the middle field
//! is `" ."`, with the leading space `List.toString` inserted. The missing-value test is
//! `s.equals(".")` on the **untrimmed** string, while the parse that follows trims. The space
//! therefore rescues every number and destroys the missing value: the same `.` that is accepted in
//! a string attribute throws in a list attribute. Nothing in either signature says so.
//!
//! The doubles are compared as raw bits, because `getTumorLogOdds` multiplies them by
//! `Math.log(10)`.

use std::io::Read;

use gatk_engine::variant_getters::{
    get_attribute_as_double_array, get_tumor_log_odds, max_element_index,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/variant_getters.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn context_with(key: &str, value: Option<Value>) -> VariantContext {
    let mut vc = VariantContext::new(
        "chr1",
        100,
        vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", false).expect("an allele"),
            Allele::from_str("G", false).expect("an allele"),
        ],
    );
    vc.stop = 100;
    if let Some(value) = value {
        vc.attributes.push((key.to_string(), value));
    }
    vc
}

fn text(value: &str) -> Value {
    Value::Str(value.to_string())
}

fn list_of_strings(values: &[&str]) -> Value {
    Value::List(values.iter().map(|v| text(v)).collect())
}

/// The dump's `array` fixtures, in its order.
fn array_cases() -> Vec<(&'static str, Option<Value>)> {
    vec![
        ("absent", None),
        ("one-string", Some(text("1.5"))),
        ("two-strings", Some(text("1.5,2.5"))),
        ("string-with-spaces", Some(text("1.5, 2.5"))),
        ("string-with-brackets", Some(text("[1.5, 2.5]"))),
        ("negative", Some(text("-1.5,-2.5"))),
        ("scientific", Some(text("1.5e3,1.5E-3"))),
        ("all-missing", Some(text("."))),
        ("one-missing", Some(text("1.5,.,2.5"))),
        ("not-a-double", Some(text("abc"))),
        ("empty-string", Some(text(""))),
        ("trailing-comma", Some(text("1.5,2.5,"))),
        ("type-suffix", Some(text("1.5f,2.5d"))),
        ("hexadecimal", Some(text("0x1p3"))),
        ("spelled-infinity", Some(text("Infinity,-Infinity"))),
        ("lower-case-inf", Some(text("inf"))),
        // The list cases the unreachable branch sends through the text path.
        ("list-of-strings", Some(list_of_strings(&["1.5", "2.5"]))),
        (
            "list-with-missing",
            Some(list_of_strings(&["1.5", ".", "2.5"])),
        ),
        ("boxed-integer", Some(Value::Int(3))),
    ]
}

fn rendered(
    result: Result<Option<Vec<f64>>, gatk_engine::variant_getters::NonDoubleValue>,
) -> String {
    match result {
        Ok(None) => "null".to_string(),
        Ok(Some(values)) => values
            .iter()
            .map(|value| (value.to_bits() as i64).to_string())
            .collect::<Vec<_>>()
            .join(","),
        Err(error) => format!("E:{}:{}", error.class(), error.message()),
    }
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let needle = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&needle))
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))[needle.len()..]
        .to_string()
}

#[test]
fn every_attribute_becomes_the_reference_double_array() {
    let golden = golden();
    for (label, attribute) in array_cases() {
        let vc = context_with("X", attribute);
        let ours = rendered(get_attribute_as_double_array(&vc, "X"));
        assert_eq!(
            ours,
            value(&golden, "array", label),
            "the array for {label}"
        );
    }
}

#[test]
fn every_tlod_matches_the_reference() {
    let golden = golden();
    let cases: Vec<(&str, Option<Value>)> = vec![
        ("absent", None),
        ("one", Some(text("1.5"))),
        ("two", Some(text("1.5,2.5"))),
        ("missing", Some(text("."))),
    ];
    for (label, attribute) in cases {
        let vc = context_with("TLOD", attribute);
        let ours = rendered(get_tumor_log_odds(&vc));
        assert_eq!(ours, value(&golden, "tlod", label), "the TLOD for {label}");
    }
}

#[test]
fn every_max_element_index_matches_the_reference() {
    let golden = golden();
    let cases: Vec<(&str, Vec<f64>)> = vec![
        ("one", vec![1.0]),
        ("ascending", vec![1.0, 2.0, 3.0]),
        ("descending", vec![3.0, 2.0, 1.0]),
        ("tie", vec![2.0, 2.0, 1.0]),
        ("tie-at-the-end", vec![1.0, 2.0, 2.0]),
        ("negatives", vec![-3.0, -1.0, -2.0]),
        ("with-nan", vec![f64::NAN, 1.0]),
        ("nan-last", vec![1.0, f64::NAN]),
        ("all-nan", vec![f64::NAN, f64::NAN]),
        ("with-infinity", vec![1.0, f64::INFINITY, 2.0]),
    ];
    for (label, values) in cases {
        let ours = max_element_index(&values).expect("a non-empty array");
        assert_eq!(
            ours.to_string(),
            value(&golden, "maxindex", label),
            "maxElementIndex for {label}"
        );
    }
}

/// The rows a port gets wrong by reading either signature.
#[test]
fn the_rows_that_the_signatures_hide() {
    let golden = golden();

    // The list rendering's space rescues every number and destroys the missing value: the same
    // "." that is accepted inside a string attribute throws inside a list attribute.
    assert_eq!(
        value(&golden, "array", "one-missing"),
        format!(
            "{},{},{}",
            (1.5f64).to_bits() as i64,
            (-1.0f64).to_bits() as i64,
            (2.5f64).to_bits() as i64
        )
    );
    assert!(value(&golden, "array", "list-with-missing").contains("non-double value ' .'"));

    // A missing element is not reported as missing anywhere downstream: it is -1, and after the
    // log conversion it is ln(10) times that.
    let missing_tlod = value(&golden, "tlod", "missing")
        .parse::<i64>()
        .expect("raw bits");
    assert_eq!(
        f64::from_bits(missing_tlod as u64),
        -std::f64::consts::LN_10
    );

    // Double.parseDouble's alphabet, which is not Rust's: type suffixes and hexadecimal floats
    // parse, and the lower-case "inf" spelling does not.
    assert!(!value(&golden, "array", "type-suffix").starts_with("E:"));
    assert!(!value(&golden, "array", "hexadecimal").starts_with("E:"));
    assert!(!value(&golden, "array", "spelled-infinity").starts_with("E:"));
    assert!(value(&golden, "array", "lower-case-inf").starts_with("E:"));

    // A tie in maxElementIndex goes to the first, and a NaN loses every comparison.
    assert_eq!(value(&golden, "maxindex", "tie"), "0");
    assert_eq!(value(&golden, "maxindex", "with-nan"), "0");
    assert_eq!(value(&golden, "maxindex", "all-nan"), "0");
}
