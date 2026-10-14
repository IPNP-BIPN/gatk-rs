//! Conformance for `String.format("%.Nf")` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/StringFormatDump.java`.
//!
//! # What this suite is for
//!
//!  * **the conversion rounds the shortest digits, not the exact expansion**, so `2.675` at two
//!    places is `2.68` where C's printf answers `2.67`;
//!  * **a large value is padded with zeros rather than expanded**;
//!  * **the rounding is HALF_UP on those digits**, not half-even;
//!  * **and the three non-finite spellings carry no sign on `NaN`.**
//!
//! Every row is compared. `Double.toString(4.9E-324)` is `4.9E-324`, longer than the shortest
//! `5.0E-324`, and at 330 places the formatting ends `...4900000`: the conversion passes
//! `Double.toString`'s digits through rather than looking at the value (#399).

use gatk_corpus as corpus;
use gatk_engine::java_format::format_decimals;
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/string_format.txt.gz"),
    )
}

fn rows() -> Vec<(String, String, String)> {
    golden()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut fields = line.splitn(3, '\t');
            (
                fields.next().expect("a kind").to_string(),
                fields.next().expect("a label").to_string(),
                fields.next().expect("a payload").to_string(),
            )
        })
        .collect()
}

fn value(label: &str) -> f64 {
    match label {
        "two-point-six-seven-five" => 2.675,
        "one-sixteenth" => 0.0625,
        "eight-thousandths" => 0.008,
        "one-point-zero-zero-five" => 1.005,
        "one-e-twenty-one" => 1e21,
        "one-e-three-hundred" => 1e300,
        "max-value" => f64::MAX,
        "min-normal" => f64::MIN_POSITIVE,
        "min-value" => 5e-324,
        "two-sevenths-percent" => 100.0 * 2.0 / 7.0,
        "ten-sixty-fifths-percent" => 100.0 * 10.0 / 65.0,
        "all-of-it" => 100.0 * 7.0 / 7.0,
        "nine-nine-nine" => 9.999,
        "zero-nine-nine-nine" => 0.999,
        "zero" => 0.0,
        "negative-zero" => -0.0,
        "negative" => -2.675,
        "nan" => f64::NAN,
        "infinity" => f64::INFINITY,
        "negative-infinity" => f64::NEG_INFINITY,
        other => panic!("no value named {other}"),
    }
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 54, "the golden's row count");
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            // `Double.toString`, which is what the conversion starts from.
            "tostring" => assert_eq!(
                java_double_to_string(value(label)),
                *payload,
                "tostring {label}"
            ),
            "format" => {
                let (places, expected) = payload.split_once('=').expect("a formatting");
                let places: usize = places.parse().expect("a scale");
                assert_eq!(
                    format_decimals(value(label), places),
                    expected,
                    "format {label} at {places}"
                );
            }
            other => panic!("no row kind {other}"),
        }
    }
}
