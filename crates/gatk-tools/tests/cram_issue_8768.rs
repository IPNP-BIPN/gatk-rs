//! Conformance for `CRAMIssue8768Detector` against GATK 4.6.2.0, compared as the whole report of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/CRAMIssue8768DetectorDump.java`.
//!
//! # What this suite is for
//!
//!  * **the suspect being the container AFTER the one that opens at position 1**;
//!  * **the count beside the rate being the total base count**, not the mismatch count;
//!  * **the fifth good container of a contig being the first one dropped**, and `--verbose`
//!    being the only thing that shows it;
//!  * **a foreign CRAM leaving a report with no body and a return code of 0**;
//!  * **an average over no containers printing as NaN**;
//!  * **and the TSV being the only place a contig name is ever resolved**.

use gatk_corpus as corpus;
use gatk_tools::cram_issue_8768::{
    analyse, report, tsv, ContainerMeta, CramHeaderInfo, RefContext,
};

/// Every CRAM here was written by the same htsjdk, which stamps this version.
const CRAM_VERSION: &str = "3.0";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_issue_8768.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field<'a>(line: &'a str, name: &str) -> &'a str {
    line.split('\t')
        .find_map(|part| part.strip_prefix(&format!("{name}=")))
        .unwrap_or_else(|| panic!("the row carries {name}"))
}

/// Every container of one fixture, in file order. The dump walks until the EOF container, so the
/// last row of a fixture is that container.
fn containers(text: &str, fixture: &str) -> Vec<ContainerMeta> {
    let prefix = format!("container\t{fixture}\t");
    let rows: Vec<&str> = text
        .lines()
        .filter(|line| line.starts_with(prefix.as_str()))
        .collect();
    rows.iter()
        .enumerate()
        .map(|(index, line)| {
            let context = match field(line, "context") {
                "UNMAPPED_UNPLACED" => RefContext::UnmappedUnplaced,
                "MULTIPLE_REFERENCE" => RefContext::Multiple,
                other => RefContext::Single(
                    other
                        .strip_prefix("SINGLE_REFERENCE: ")
                        .expect("a single-reference context")
                        .parse()
                        .expect("a reference id"),
                ),
            };
            ContainerMeta {
                context,
                start: field(line, "start").parse().expect("a start"),
                span: field(line, "span").parse().expect("a span"),
                slices: field(line, "slices").parse().expect("a slice count"),
                reference_required: match field(line, "reference-required") {
                    "none" => None,
                    other => Some(other.parse().expect("a boolean")),
                },
                embedded_reference: field(line, "embedded").parse().expect("a content id"),
                bases: field(line, "bases").parse().expect("a base count"),
                mismatches: field(line, "mismatches").parse().expect("a mismatch count"),
                is_eof: index + 1 == rows.len(),
            }
        })
        .collect()
}

fn escaped(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn code(text: &str, label: &str) -> i32 {
    let prefix = format!("code\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries code/{label}"))
        .parse()
        .expect("a return code")
}

/// The sequence dictionary, in the order the ids index it.
fn dictionary(text: &str) -> Vec<String> {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix("dict\tcram8768="))
            .expect("the golden carries the dictionary"),
    )
    .lines()
    .filter(|line| line.starts_with("@SQ"))
    .map(|line| {
        line.split('\t')
            .find_map(|part| part.strip_prefix("SN:"))
            .expect("a sequence name")
            .to_string()
    })
    .collect()
}

/// `java.util.Base64.getEncoder()`, for the twenty bytes the CRAM header carries.
fn encode_base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::new();
    for chunk in bytes.chunks(3) {
        let mut block = [0u8; 3];
        block[..chunk.len()].copy_from_slice(chunk);
        let value = u32::from(block[0]) << 16 | u32::from(block[1]) << 8 | u32::from(block[2]);
        for index in 0..4 {
            if index <= chunk.len() {
                out.push(ALPHABET[(value >> (18 - 6 * index)) as usize & 0x3f] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// The CRAM id htsjdk writes: the file name, zero padded to twenty bytes.
fn header_for(file_name: &str) -> CramHeaderInfo {
    let mut id = [0u8; 20];
    id[..file_name.len()].copy_from_slice(file_name.as_bytes());
    CramHeaderInfo {
        file_name: file_name.to_string(),
        version: CRAM_VERSION.to_string(),
        id_base64: encode_base64(&id),
    }
}

/// label, fixture, file name, verbose, threshold, echo.
const RUNS: &[(&str, &str, &str, bool, f64, bool)] = &[
    ("default", "input", "reads.cram", false, 0.05, false),
    ("tsv", "input", "reads.cram", false, 0.05, false),
    ("verbose", "input", "reads.cram", true, 0.05, false),
    ("threshold-low", "input", "reads.cram", false, 0.01, false),
    ("threshold-high", "input", "reads.cram", false, 0.9, false),
    ("echo", "input", "reads.cram", false, 0.05, true),
    ("clean", "clean", "clean.cram", false, 0.05, false),
    ("foreign", "foreign", "foreign.cram", false, 0.05, false),
    ("empty", "empty", "empty.cram", false, 0.05, false),
];

#[test]
fn every_report_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, fixture, file_name, verbose, threshold, echo) in RUNS {
        let analysis = analyse(&containers(&text, fixture), *verbose).expect("an analysis");
        let produced = report(&header_for(file_name), &analysis, *threshold, *echo);
        assert_eq!(produced.text, escaped(&text, "report", label), "{label}");
        assert_eq!(produced.stdout, escaped(&text, "stdout", label), "{label}");
        assert_eq!(produced.code, code(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's runs");
}

#[test]
fn the_tsv_matches_the_golden() {
    let text = golden();
    let analysis = analyse(&containers(&text, "input"), false).expect("an analysis");
    assert_eq!(
        tsv(&analysis, "reads.cram", &dictionary(&text)),
        escaped(&text, "tsv", "tsv")
    );
}

/// The container that opens at position 1 is reported GOOD, and the one after it is the suspect.
#[test]
fn the_suspect_is_the_container_after_the_one_at_position_one() {
    let text = golden();
    let containers = containers(&text, "input");
    let analysis = analyse(&containers, false).expect("an analysis");

    // chr1 and chr3 both open at 1, so both are bad contigs, keyed by reference id and not by name.
    assert_eq!(
        analysis
            .bad_contigs
            .iter()
            .map(|(name, _)| name.as_str())
            .collect::<Vec<&str>>(),
        vec!["SINGLE_REFERENCE: 0", "SINGLE_REFERENCE: 2"]
    );
    for (_, suspects) in &analysis.bad_contigs {
        assert_eq!(suspects.len(), 1, "one suspect per contig that opens at 1");
        assert_eq!(suspects[0].container_ordinal, 2);
        assert_ne!(suspects[0].alignment_start, 1, "the suspect's own start");
    }
    // And every container that does open at 1 is in the good list.
    for stats in &analysis.good {
        assert!(!stats.is_bad);
    }
    assert_eq!(
        analysis
            .good
            .iter()
            .filter(|stats| stats.alignment_start == 1)
            .count(),
        2
    );
    // chr2 opens at 5, so nothing on it is ever suspected, whatever its rates.
    assert!(!analysis
        .bad_contigs
        .iter()
        .any(|(name, _)| name == "SINGLE_REFERENCE: 1"));
}

/// The number printed after the slash is the container's TOTAL BASE COUNT, and the golden reports
/// the real mismatch counts beside it.
#[test]
fn the_count_beside_the_rate_is_the_base_count() {
    let text = golden();
    let containers = containers(&text, "input");
    let analysis = analyse(&containers, true).expect("an analysis");

    let measured: Vec<(i64, i64)> = containers
        .iter()
        .filter(|meta| meta.context.is_mapped_single_ref())
        .map(|meta| (meta.bases, meta.mismatches))
        .collect();
    // The mismatches are not all equal, so a count that never varies cannot be the mismatch count.
    assert!(measured.iter().any(|(_, mismatches)| *mismatches != 10));
    for stats in analysis.good.iter().chain(
        analysis
            .bad_contigs
            .iter()
            .flat_map(|(_, containers)| containers),
    ) {
        assert_eq!(stats.mismatch_count, 30, "every container holds 30 bases");
    }
    assert!(escaped(&text, "report", "verbose").contains("Mismatch Rate/Count: 0.733333/30"));
    // 22 of 30, which is the rate. The 22 is nowhere in the report.
    assert!(measured.contains(&(30, 22)));
    assert!(!escaped(&text, "report", "verbose").contains("/22"));
}

/// chr2 carries five containers and only four are reported, because the counter that caps them is
/// set to 1 by the new-context branch rather than to 0.
#[test]
fn the_fifth_good_container_of_a_contig_needs_verbose() {
    let text = golden();
    let containers = containers(&text, "input");
    let on_chr2 = |verbose: bool| {
        analyse(&containers, verbose)
            .expect("an analysis")
            .good
            .iter()
            .filter(|stats| stats.context == RefContext::Single(1))
            .map(|stats| stats.container_ordinal)
            .collect::<Vec<i32>>()
    };
    assert_eq!(on_chr2(false), vec![1, 2, 3, 4]);
    assert_eq!(on_chr2(true), vec![1, 2, 3, 4, 5]);
    // The fixture really does hold five containers on chr2, so the fifth was dropped rather than
    // absent.
    assert_eq!(
        containers
            .iter()
            .filter(|meta| meta.context == RefContext::Single(1))
            .count(),
        5
    );
}

/// A two-slice container stops the walk from inside the loop, so the report keeps its header lines
/// and gains nothing else, and the run still returns 0.
#[test]
fn a_foreign_cram_leaves_a_report_with_no_body() {
    let text = golden();
    let analysis = analyse(&containers(&text, "foreign"), false).expect("an analysis");
    assert_eq!(
        analysis.foreign.as_deref(),
        Some("Multi-slice container detected. This file was not written by GATK or Picard.")
    );
    assert!(analysis.good.is_empty());
    let produced = report(&header_for("foreign.cram"), &analysis, 0.05, false);
    assert_eq!(produced.code, 0, "nothing was judged, and nothing failed");
    assert!(!produced.text.contains("Average mismatch rate"));
    assert!(!produced.text.contains("Presumed GOOD Containers"));
    assert_eq!(
        produced.stdout, "",
        "not even the summary reaches the console"
    );
}

/// A CRAM holding nothing but its EOF container has no good containers, so the average is 0.0/0.
#[test]
fn an_empty_cram_averages_over_nothing_and_prints_nan() {
    let text = golden();
    let containers = containers(&text, "empty");
    assert_eq!(containers.len(), 1);
    assert!(containers[0].is_eof);
    // The EOF container carries no compression header, which is why that field is optional.
    assert_eq!(containers[0].reference_required, None);

    let analysis = analyse(&containers, false).expect("an analysis");
    assert!(analysis.good.is_empty() && analysis.bad_contigs.is_empty());
    let produced = report(&header_for("empty.cram"), &analysis, 0.05, false);
    assert!(produced
        .text
        .contains("Average mismatch rate for presumed good containers: NaN"));
    assert_eq!(produced.code, 0);
}

/// The unmapped container advances the ordinal and is never recorded, so it is the container that
/// closes chr3 without appearing anywhere in the report.
#[test]
fn an_unmapped_container_is_counted_but_never_recorded() {
    let text = golden();
    let containers = containers(&text, "input");
    let unmapped: Vec<&ContainerMeta> = containers
        .iter()
        .filter(|meta| meta.context == RefContext::UnmappedUnplaced)
        .collect();
    assert_eq!(
        unmapped.len(),
        2,
        "the unmapped reads, and the EOF container"
    );
    assert!(!unmapped[0].is_eof && unmapped[1].is_eof);

    let analysis = analyse(&containers, true).expect("an analysis");
    for stats in &analysis.good {
        assert!(stats.context.is_mapped_single_ref());
    }
    // chr3's suspect was closed by that container rather than by the end of the loop.
    assert_eq!(
        analysis.bad_contigs.last().expect("a contig").0,
        "SINGLE_REFERENCE: 2"
    );
}

/// The id in the header is the file name zero padded to twenty bytes, which is what makes the
/// three header lines derivable rather than copied out of the golden.
#[test]
fn the_cram_id_is_the_file_name_padded_to_twenty_bytes() {
    let text = golden();
    for (_, _, file_name, _, _, _) in RUNS {
        let expected = format!("CRAM ID Contents: {}", header_for(file_name).id_base64);
        assert!(
            escaped(&text, "report", "default").contains("CRAM ID Contents: ")
                && text.contains(&expected.replace('\n', "\\n")),
            "{file_name}"
        );
    }
    assert_eq!(
        corpus::decode_base64(&header_for("reads.cram").id_base64),
        b"reads.cram\0\0\0\0\0\0\0\0\0\0".to_vec()
    );
}
