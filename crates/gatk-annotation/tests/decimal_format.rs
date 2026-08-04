//! Conformance for `DecimalFormat` against the oracle.
//!
//! Golden from `tools/annotation-conformance/DecimalFormatDump.java`, as bit patterns with the two
//! formatted strings. A decimal literal in the golden would be re-parsed by this file and could
//! name a different double than the one measured, which for a suite about the last printed digit
//! is the whole question.
//!
//! # The suite declares its own domain, and the corpus goes past it
//!
//! The port reproduces the reference exactly for values below 2^53 whose shortest decimal form
//! needs at most fifteen significant digits, which is the whole of what `AllelePseudoDepth` can
//! produce. Past that the divergence is in Java's *digit generation* rather than in its rounding,
//! and it has two shapes, both deliberately present in the corpus:
//!
//!  * **sixteen significant digits.** `6.985838094673373e14` is exactly `698583809467337.25`, so
//!    the two sixteen-digit forms are equidistant; Java picks the even one, Rust the larger;
//!  * **above 2^53.** Java stops printing the shortest form and starts printing digits from the
//!    value itself.
//!
//! Every row is compared, including those. What the test asserts is not that the divergences are
//! few but that they are all of that **shape**: a divergence below 2^53 with fewer than sixteen
//! significant digits would mean a rounding rule is wrong, which is a different failure and the one
//! this suite exists to catch. A corpus trimmed to what already passes would report a hundred per
//! cent and mean nothing, and a quarantine listed value by value would keep passing after the rule
//! behind it stopped being true.

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

/// Whether the value is outside what this port claims, and therefore expected to diverge.
///
/// Computed from the value, not looked up: a list of bit patterns would keep passing after the
/// rule behind it stopped being true.
fn past_the_declared_domain(value: f64) -> bool {
    if !value.is_finite() {
        return false;
    }
    if value.abs() >= 9_007_199_254_740_992.0 {
        return true;
    }
    let significant = format!("{value:e}")
        .split('e')
        .next()
        .expect("mantissa")
        .chars()
        .filter(char::is_ascii_digit)
        .count();
    significant >= 16
}

#[test]
fn every_divergence_is_digit_generation_and_not_rounding() {
    let text = golden();
    let mut compared = 0usize;
    let mut diverging = Vec::new();
    let mut rounding = Vec::new();

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
        if two == fields[2] && four == fields[3] {
            continue;
        }
        let row = format!(
            "{}: ours {two:?}/{four:?}, reference {:?}/{:?}",
            fields[1], fields[2], fields[3]
        );
        // Every row is compared. What is allowed to differ is a *shape*, not a list of bit
        // patterns: a list would keep passing after the rule behind it stopped being true.
        if past_the_declared_domain(value) {
            diverging.push(row);
        } else {
            rounding.push(row);
        }
    }

    assert!(
        rounding.is_empty(),
        "{} value(s) diverge inside the declared domain, which means a rounding rule is wrong \
         rather than a digit string being unavailable:\n{}",
        rounding.len(),
        rounding.join("\n")
    );
    // The corpus must still reach past the domain, or the limit is asserted and not measured.
    assert!(
        !diverging.is_empty(),
        "no value diverged at all, so the corpus no longer contains the cases that locate the \
         edge of what this port reproduces"
    );
    println!(
        "DecimalFormat: {compared} strings compared, {} value(s) diverge and all of them are \
         Java emitting digits the shortest form does not have:\n  {}",
        diverging.len(),
        diverging.join("\n  ")
    );
}
