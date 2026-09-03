//! Conformance for the header `SelectVariants` writes, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SelectVariantsHeaderDump.java`.
//!
//! # What this suite is for
//!
//!  * **AN, AC, AF and DP are replaced, not merged**, so an input's own descriptions for those
//!    four are gone and an AF appears in a file that never had one;
//!  * **the keep-original arguments add lines** whether or not a record carries them;
//!  * **the drops run last**, so `--drop-info-annotation AC` removes the standard line just added;
//!  * **the sample columns are the selected set**, sorted, and sites-only leaves none;
//!  * **and `--add-output-vcf-command-line false` removes the `##source` line as well.**
//!
//! # What is compared, and what is not
//!
//! The `##GATKCommandLine` line carries the run's own wall-clock time and the dump elides it, so
//! this compares every OTHER header line and asserts that the elided line is where the reference
//! put it. The line's own construction is the runner's, and it is measured by the
//! `expanded-command-line` suite.

use gatk_corpus as corpus;
use gatk_tools::select_variants_header::{output_header, HeaderArguments};
use htsjdk_vcf::reader::read_vcf;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_header.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The rows of one kind, as `(label, value)`.
fn rows(text: &str, kind: &str) -> Vec<(String, String)> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            let found = parts.next()?;
            if found != kind {
                return None;
            }
            Some((
                parts.next()?.to_string(),
                parts.next().unwrap_or("").to_string(),
            ))
        })
        .collect()
}

/// The input the harness wrote, which every case reads.
fn input(text: &str) -> String {
    let (_, escaped) = rows(text, "input")
        .into_iter()
        .next()
        .expect("the input row");
    unescape(&escaped)
}

/// The header lines the reference wrote for one case, in its own order.
fn reference_lines(text: &str, case: &str) -> Vec<String> {
    rows(text, "header")
        .into_iter()
        .filter(|(label, _)| label == case)
        .map(|(_, line)| unescape(&line))
        .collect()
}

fn reference_samples(text: &str, case: &str) -> Vec<String> {
    let (_, value) = rows(text, "samples")
        .into_iter()
        .find(|(label, _)| label == case)
        .unwrap_or_else(|| panic!("no samples row for {case}"));
    if value.is_empty() {
        Vec::new()
    } else {
        value.split(',').map(str::to_string).collect()
    }
}

/// The arguments each case ran with, in the harness's own order.
fn arguments(case: &str) -> HeaderArguments {
    let mut arguments = HeaderArguments {
        keep_original_chr_counts: false,
        keep_original_depth: false,
        info_annotations_to_drop: Vec::new(),
        genotype_annotations_to_drop: Vec::new(),
        add_output_vcf_command_line: false,
        tool_command_line: None,
        samples: vec!["s0".to_string(), "s1".to_string(), "s2".to_string()],
    };
    match case {
        "plain" => {}
        "keep-original-ac" => arguments.keep_original_chr_counts = true,
        "keep-original-dp" => arguments.keep_original_depth = true,
        "keep-original-both" => {
            arguments.keep_original_chr_counts = true;
            arguments.keep_original_depth = true;
        }
        "drop-info-qd" => arguments.info_annotations_to_drop = vec!["QD".to_string()],
        "drop-info-ac" => arguments.info_annotations_to_drop = vec!["AC".to_string()],
        "drop-genotype" => arguments.genotype_annotations_to_drop = vec!["XX".to_string()],
        "drop-absent" => arguments.info_annotations_to_drop = vec!["NOPE".to_string()],
        "subset" => arguments.samples = vec!["s0".to_string(), "s2".to_string()],
        "sites-only" => arguments.samples = Vec::new(),
        "command-line" => {
            // The line itself is elided in the golden, so the port's own value would not be
            // comparable; what this case measures is that the `##source` line comes WITH it.
            arguments.add_output_vcf_command_line = true;
        }
        other => panic!("no arguments for {other}"),
    }
    arguments
}

const CASES: &[&str] = &[
    "plain",
    "keep-original-ac",
    "keep-original-dp",
    "keep-original-both",
    "drop-info-qd",
    "drop-info-ac",
    "drop-genotype",
    "drop-absent",
    "subset",
    "sites-only",
    "command-line",
];

/// The port's header for one case, as the lines a file would carry.
fn written(text: &str, case: &str) -> Vec<String> {
    let file = read_vcf(&input(text)).expect("the input parses");
    let header = output_header(&file.header, &arguments(case));
    header
        .write()
        .lines()
        .filter(|line| line.starts_with("##"))
        .map(str::to_string)
        .collect()
}

#[test]
fn every_case_writes_the_reference_s_header_lines() {
    let text = golden();
    for case in CASES {
        let expected: Vec<String> = reference_lines(&text, case)
            .into_iter()
            // The command line holds the run's own time and the dump elided it; the port cannot
            // reproduce an instant, and the `##source` line beside it is what this case measures.
            .filter(|line| !line.starts_with("##GATKCommandLine="))
            .collect();
        let produced: Vec<String> = written(&text, case)
            .into_iter()
            .filter(|line| !line.starts_with("##GATKCommandLine="))
            .collect();
        assert_eq!(produced, expected, "case {case}");
    }
}

#[test]
fn the_four_standard_lines_replace_whatever_the_input_declared() {
    let text = golden();
    let plain = reference_lines(&text, "plain");
    for (id, description) in [
        (
            "AC",
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        ),
        (
            "AF",
            "Allele Frequency, for each ALT allele, in the same order as listed",
        ),
        ("AN", "Total number of alleles in called genotypes"),
        (
            "DP",
            "Approximate read depth; some reads may have been filtered",
        ),
    ] {
        let wanted = format!("##INFO=<ID={id},");
        let line = plain
            .iter()
            .find(|line| line.starts_with(&wanted))
            .unwrap_or_else(|| panic!("the reference wrote no INFO line for {id}"));
        assert!(
            line.contains(description),
            "the reference's {id} is not htsjdk's: {line}"
        );
        assert!(
            written(&text, "plain").iter().any(|own| own == line),
            "the port did not replace {id}"
        );
    }
}

#[test]
fn dropping_ac_removes_the_line_the_replacement_had_just_added() {
    let text = golden();
    let reference = reference_lines(&text, "drop-info-ac");
    assert!(
        !reference
            .iter()
            .any(|line| line.starts_with("##INFO=<ID=AC,")),
        "the reference kept an AC line"
    );
    assert_eq!(written(&text, "drop-info-ac"), reference);
}

#[test]
fn the_source_line_goes_with_the_command_line() {
    let text = golden();
    assert!(
        !reference_lines(&text, "plain")
            .iter()
            .any(|line| line.starts_with("##source=")),
        "--add-output-vcf-command-line false left a source line"
    );
    assert!(
        reference_lines(&text, "command-line")
            .iter()
            .any(|line| line == "##source=SelectVariants"),
        "the command-line case has no source line"
    );
    assert!(
        written(&text, "command-line")
            .iter()
            .any(|line| line == "##source=SelectVariants"),
        "the port writes no source line"
    );
}

#[test]
fn the_sample_columns_are_the_selected_set() {
    let text = golden();
    assert_eq!(reference_samples(&text, "subset"), vec!["s0", "s2"]);
    assert!(reference_samples(&text, "sites-only").is_empty());
    for case in CASES {
        let file = read_vcf(&input(&text)).expect("the input parses");
        let header = output_header(&file.header, &arguments(case));
        assert_eq!(
            header.samples,
            reference_samples(&text, case),
            "case {case}"
        );
    }
}
