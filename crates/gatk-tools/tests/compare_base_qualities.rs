//! Conformance for `CompareBaseQualities` against GATK 4.6.2.0, compared as its report text.
//!
//! Golden from `tools/readfilter-conformance/CompareBaseQualitiesDump.java`. The nine input BAMs
//! travel in full, base64, so the port reads the same records.
//!
//! # What this suite is for
//!
//!  * **there are no read filters**, so a duplicate and a vendor failure are counted;
//!  * **a secondary or supplementary read is skipped in each file independently**;
//!  * **the two refusals about order and count are different messages**;
//!  * **ragged qualities name the two lengths and not the read**;
//!  * **the summary collapses onto diagonals** and the percentage is `%.4f`;
//!  * **the diff is QRead1 - QRead2**, so swapping the inputs flips every sign;
//!  * **quantization can make the two halves of one report disagree**;
//!  * **and `--round-down-quantized` alone is refused before anything is read**.
//!
//! The `strict-validation` run is not reproduced here: it fails inside htsjdk's reader under STRICT
//! stringency, which belongs to the reader and not to this tool. The golden holds it so the reason
//! is on record.

use gatk_corpus as corpus;
use gatk_tools::compare_base_qualities::{compare_base_qualities, CompareArguments, CompareError};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/compare_base_qualities.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The records of one fixture, read back off the golden's own bytes.
fn fixture(text: &str, label: &str) -> Vec<BamRecord> {
    let encoded = rows(text, "fixture")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no fixture {label}"))[1]
        .to_string();
    let bytes = corpus::decode_base64(&encoded);
    let decompressed = htsjdk_bgzf::read::decompress_all(&bytes).expect("the fixture is BGZF");
    let reader = htsjdk_bam::reader::BamReader::new(&decompressed).expect("the fixture opens");
    reader.map(|record| record.expect("a record")).collect()
}

/// The two fixtures and the arguments each labelled run used.
fn configuration(label: &str) -> (&str, &str, CompareArguments) {
    let default = CompareArguments::default();
    match label {
        "identical" => ("same-a", "same-b", default),
        "shifted" => ("same-a", "shifted", default),
        "shifted-reversed" => ("shifted", "same-a", default),
        "renamed" => ("same-a", "renamed", default),
        "shorter" => ("same-a", "shorter", default),
        "ragged" => ("same-a", "ragged", default),
        "with-secondary" => ("same-a", "with-secondary", default),
        "flagged" => ("same-a", "flagged", default),
        "quantized" => (
            "same-a",
            "shifted",
            CompareArguments {
                static_quantization_quals: vec![10, 20, 30],
                ..CompareArguments::default()
            },
        ),
        "quantized-round-down" => (
            "same-a",
            "shifted",
            CompareArguments {
                static_quantization_quals: vec![10, 20, 30],
                round_down: true,
                ..CompareArguments::default()
            },
        ),
        "round-down-alone" => (
            "same-a",
            "shifted",
            CompareArguments {
                round_down: true,
                ..CompareArguments::default()
            },
        ),
        "throw-on-diff" => (
            "same-a",
            "shifted",
            CompareArguments {
                throw_on_diff: true,
                ..CompareArguments::default()
            },
        ),
        "throw-on-same" => (
            "same-a",
            "same-b",
            CompareArguments {
                throw_on_diff: true,
                ..CompareArguments::default()
            },
        ),
        other => panic!("no run {other}"),
    }
}

#[test]
fn every_report_is_the_reference() {
    let text = golden();
    let reports = rows(&text, "report");
    let results = rows(&text, "result");
    assert_eq!(reports.len(), 8, "eight runs finish");

    for row in &reports {
        let label = row[0];
        let (first, second, arguments) = configuration(label);
        let ours =
            compare_base_qualities(&fixture(&text, first), &fixture(&text, second), &arguments)
                .unwrap_or_else(|error| panic!("{label} was refused: {}", error.message()));

        assert_eq!(ours.report, unescape(row[1]), "report/{label}");

        let expected: i32 = results
            .iter()
            .find(|other| other[0] == label)
            .expect("every finished run dumps its return value")[1]
            .parse()
            .expect("a number");
        assert_eq!(ours.exit_code, expected, "result/{label}");
    }
}

#[test]
fn every_refusal_this_port_owns_is_the_reference() {
    let text = golden();

    // `strict-validation` belongs to htsjdk's reader, not to this tool, so it is skipped here.
    let ours = |label: &str| -> CompareError {
        let (first, second, arguments) = configuration(label);
        compare_base_qualities(&fixture(&text, first), &fixture(&text, second), &arguments)
            .expect_err("this run is refused")
    };

    for row in rows(&text, "error") {
        let label = row[0];
        if label == "strict-validation" {
            continue;
        }
        let error = ours(label);
        assert!(
            row[1].ends_with(&error.message()),
            "error/{label}: {} against {}",
            error.message(),
            row[1]
        );
    }
}

/// Quantization can make the two halves of one report contradict each other.
#[test]
fn the_binned_half_can_disagree_with_the_raw_one() {
    let text = golden();
    let report = rows(&text, "report")
        .into_iter()
        .find(|row| row[0] == "quantized")
        .expect("the quantized run")[1]
        .to_string();
    let report = unescape(&report);

    let (raw, binned) = report
        .split_once("-----------CompareMatrix-binned summary------------")
        .expect("two halves");
    assert!(
        raw.contains("diff\tcount\t%total"),
        "the raw half has differences"
    );
    assert!(
        binned.contains("all 8 quality scores are the same"),
        "the binned half has none"
    );
}
