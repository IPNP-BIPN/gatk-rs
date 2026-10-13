//! Conformance for `Double.toString` against GATK 4.6.2.0, over 1059 values.
//!
//! Golden from `tools/readfilter-conformance/DoubleToStringDump.java`.
//!
//! # What this suite is for
//!
//! The pre-JDK19 `FloatingDecimal` does not always emit the shortest digits, and while
//! `java_double_to_string` deferred to Rust's shortest-round-trip formatter, fifteen of these 1059
//! values differed (#399): the smallest subnormal, `1e23`, and values where Java emits a sixteenth
//! or seventeenth digit. It now takes Java 17's digits from `jmath::to_string`, reconstructed from
//! the oracle's behaviour (htsjdk-rs decision 0046), and every row matches.

use gatk_corpus as corpus;
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/double_to_string.txt.gz"),
    )
}

#[test]
fn every_row_matches() {
    let text = golden();
    let rows: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(rows.len(), 1059, "the golden's row count");

    let mut divergences: Vec<String> = Vec::new();
    for row in &rows {
        let mut fields = row.split('\t');
        assert_eq!(fields.next(), Some("tostring"));
        let bits = fields.next().expect("the bits");
        let expected = fields.next().expect("a rendering");
        let ours = java_double_to_string(f64::from_bits(
            u64::from_str_radix(bits, 16).expect("sixteen hex digits"),
        ));
        if ours != expected {
            divergences.push(format!("{bits}: java={expected}, ours={ours}"));
        }
    }
    assert!(
        divergences.is_empty(),
        "renderings differ:\n{}",
        divergences.join("\n")
    );
}
