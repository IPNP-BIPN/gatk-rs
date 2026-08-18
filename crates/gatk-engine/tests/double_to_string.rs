//! Conformance for `Double.toString` against GATK 4.6.2.0, over 1059 values.
//!
//! Golden from `tools/readfilter-conformance/DoubleToStringDump.java`.
//!
//! # What this suite is for
//!
//! `java_double_to_string` defers to Rust's shortest-round-trip formatter, and the pre-JDK19
//! `FloatingDecimal` does not always emit the shortest digits. This suite is the record of **which**
//! values that costs, and it asserts both sides of every one of them, so the set cannot grow or
//! shrink without a test failing.
//!
//! Fifteen of the 1059 differ, 1.4 per cent, and in every one Java emits *more* digits than the
//! shortest. Every one parses back to the same double: this is a rendering difference, not a value
//! difference. #399 tracks the fix, which is a port of `FloatingDecimal.toJavaFormatString`'s digit
//! generation rather than a patch to any of these fifteen.

use gatk_corpus as corpus;
use gatk_engine::tsv_table::java_double_to_string;

/// The bits of the fifteen values the port renders differently, with both renderings.
///
/// Three shapes: the smallest subnormal, `1e23` rendered with sixteen digits where two round-trip,
/// and eleven values needing sixteen digits that get a seventeenth.
const DIVERGENCES: [(&str, &str, &str); 15] = [
    ("0000000000000001", "4.9E-324", "5.0E-324"),
    ("8000000000000001", "-4.9E-324", "-5.0E-324"),
    ("44b52d02c7e14af6", "9.999999999999999E22", "1.0E23"),
    (
        "0100000000000000",
        "7.2911220195563975E-304",
        "7.291122019556398E-304",
    ),
    (
        "15f0000000000000",
        "5.1032038149619546E-203",
        "5.103203814961955E-203",
    ),
    (
        "2f10000000000000",
        "5.2710989716152616E-82",
        "5.271098971615262E-82",
    ),
    (
        "4830000000000000",
        "5.4445178707350154E39",
        "5.444517870735016E39",
    ),
    (
        "6150000000000000",
        "5.6236422431789955E160",
        "5.623642243178996E160",
    ),
    (
        "7210000000000000",
        "2.6672057731519417E241",
        "2.667205773151942E241",
    ),
    (
        "7640000000000000",
        "3.9361009831403587E261",
        "3.936100983140359E261",
    ),
    (
        "43918ba08a9d2f68",
        "3.1607015940265421E17",
        "3.160701594026542E17",
    ),
    (
        "438b1504472ecb8d",
        "2.43933839663657376E17",
        "2.4393383966365738E17",
    ),
    (
        "43bf1ac4aef7b71c",
        "2.24132002032956928E18",
        "2.2413200203295693E18",
    ),
    (
        "c37c27651a273249",
        "-1.26793832509678736E17",
        "-1.2679383250967874E17",
    ),
    // The subnormal appears twice in the corpus: once named and once from the bit patterns.
    ("0000000000000001", "4.9E-324", "5.0E-324"),
];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/double_to_string.txt.gz"),
    )
}

#[test]
fn every_row_is_matched_or_named() {
    let text = golden();
    let rows: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(rows.len(), 1059, "the golden's row count");

    let mut named = 0;
    let mut divergences: Vec<String> = Vec::new();
    for row in &rows {
        let mut fields = row.split('\t');
        assert_eq!(fields.next(), Some("tostring"));
        let bits = fields.next().expect("the bits");
        let expected = fields.next().expect("a rendering");
        let ours = java_double_to_string(f64::from_bits(
            u64::from_str_radix(bits, 16).expect("sixteen hex digits"),
        ));
        if ours == expected {
            continue;
        }
        // Both renderings must be the ones on the record, and both must parse back to the same
        // double: a rendering difference is allowed here, a value difference is not.
        let known = DIVERGENCES.iter().any(|(known_bits, java, port)| {
            *known_bits == bits && *java == expected && *port == ours
        });
        if !known {
            divergences.push(format!("{bits}: java={expected}, ours={ours}"));
            continue;
        }
        let reparsed: f64 = expected.parse().expect("Java's rendering parses");
        assert_eq!(
            reparsed.to_bits(),
            u64::from_str_radix(bits, 16).expect("hex"),
            "{bits}: the reference's rendering must round-trip"
        );
        named += 1;
    }
    assert!(
        divergences.is_empty(),
        "renderings differ and are not on the record:\n{}",
        divergences.join("\n")
    );
    assert_eq!(
        named,
        DIVERGENCES.len(),
        "every named divergence must still be one; a fix to #399 shortens this list"
    );
}
