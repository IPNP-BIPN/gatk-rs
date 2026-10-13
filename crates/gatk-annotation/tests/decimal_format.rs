//! Conformance for `DecimalFormat` against the oracle.
//!
//! Golden from `tools/annotation-conformance/DecimalFormatDump.java`, as bit patterns with the two
//! formatted strings. A decimal literal in the golden would be re-parsed by this file and could
//! name a different double than the one measured, which for a suite about the last printed digit
//! is the whole question.
//!
//! # Every row matches, past 2^53 too
//!
//! Above 2^53 the reference stops printing the shortest decimal form: sometimes it prints the
//! double's exact value, sometimes that value rounded to eighteen significant digits. Those are
//! Java 17's pre-Schubfach digits, and this suite used to accept them as divergences of that shape.
//! The digits now come from `jmath::to_string`, which reproduces Java 17's from measured behaviour
//! (htsjdk-rs decision 0046), so every row must match, and the corpus must still reach past 2^53,
//! where the difference from the shortest form lives.

use std::io::Read;

use gatk_annotation::decimal_format::{DEPTH_FORMAT, FRACTION_FORMAT};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/decimal_format.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

#[test]
fn every_row_matches_including_past_two_to_the_fifty_three() {
    let text = golden();
    let mut compared = 0usize;
    let mut past = 0usize;
    let mut diverging = Vec::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields[0] != "f" {
            continue;
        }
        let value = f64::from_bits(u64::from_str_radix(fields[1], 16).expect("hex bits"));
        let (two, four) = (DEPTH_FORMAT.format(value), FRACTION_FORMAT.format(value));
        compared += 2;
        if value.is_finite() && value.abs() >= 9_007_199_254_740_992.0 {
            past += 1;
        }
        if two != fields[2] || four != fields[3] {
            diverging.push(format!(
                "{}: ours {two:?}/{four:?}, reference {:?}/{:?}",
                fields[1], fields[2], fields[3]
            ));
        }
    }

    assert!(
        diverging.is_empty(),
        "{} value(s) diverge:\n{}",
        diverging.len(),
        diverging.join("\n")
    );
    assert!(
        past > 0,
        "the corpus no longer reaches past 2^53, where Java 17 leaves the shortest form"
    );
    println!("DecimalFormat: {compared} strings compared, {past} values past 2^53, all match");
}
