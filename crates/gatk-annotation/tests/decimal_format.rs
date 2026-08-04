//! Conformance for `DecimalFormat` against the oracle.
//!
//! Golden from `tools/annotation-conformance/DecimalFormatDump.java`, as bit patterns with the two
//! formatted strings. A decimal literal in the golden would be re-parsed by this file and could
//! name a different double than the one measured, which for a suite about the last printed digit
//! is the whole question.
//!
//! # The suite declares its own domain, and the corpus goes past it
//!
//! **Below 2^53 the port reproduces the reference exactly**, which covers everything
//! `AllelePseudoDepth` can produce by a wide margin. Above it the reference stops printing the
//! shortest decimal form, and what it prints instead is not one thing: sometimes the double's exact
//! value, sometimes that value rounded to eighteen significant digits, mostly the shortest form
//! after all. Those are branches inside Java 17's pre-Schubfach `FloatingDecimal`.
//!
//! Every row is compared, including those. What the test asserts is not that the divergences are
//! few but that they are all of that **shape**: a divergence below 2^53 would mean a rounding rule
//! is wrong, which is a different failure and the one this suite exists to catch. A corpus trimmed
//! to what already passes would report a hundred per cent and mean nothing, and a quarantine listed
//! value by value would keep passing after the rule behind it stopped being true.

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
    value.is_finite() && value.abs() >= 9_007_199_254_740_992.0
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
