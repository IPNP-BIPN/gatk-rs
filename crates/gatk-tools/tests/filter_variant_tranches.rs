//! Conformance for `FilterVariantTranches` against GATK 4.6.2.0, compared as the FILTER of every
//! written record, the FILTER header lines, and all five refusals.
//!
//! Golden from `tools/readfilter-conformance/FilterVariantTranchesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the cutoff is a truncated index into the resource scores sorted descending**;
//!  * **a score exactly on the cutoff is filtered**;
//!  * **the band is found afterwards**, and a score below every cutoff takes the last one;
//!  * **the tranche list is deduplicated and sorted**, so the same two given out of order and
//!    repeated produce the same output;
//!  * **a record with no score is written and passes**;
//!  * **and `--invalidate-previous-filters` clears the header's FILTER lines too**.

use gatk_corpus as corpus;
use gatk_tools::filter_variant_tranches::{
    cutoffs, filter_string_from_score, filters_of, is_tranche_filtered, score,
    tranche_header_lines, validate_tranches, FilterVariantTranchesError, INDEL_STRING, SNP_STRING,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_variant_tranches.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

/// One record, decoded as far as this tool looks at it.
struct Record {
    start: i32,
    reference: String,
    alternates: Vec<String>,
    filters: Vec<String>,
    score: Option<f64>,
    line: String,
}

impl Record {
    /// `isSNP()` for the biallelic records this fixture carries.
    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|alt| alt.len() == 1)
    }
}

fn records(whole: &str) -> Vec<Record> {
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            Record {
                start: field[1].parse().expect("a position"),
                reference: field[3].to_string(),
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                filters: match field[6] {
                    "." | "PASS" => Vec::new(),
                    list => list.split(';').map(|name| name.to_string()).collect(),
                },
                score: field[7]
                    .split(';')
                    .find_map(|entry| entry.strip_prefix("SCORE="))
                    .and_then(score),
                line: line.to_string(),
            }
        })
        .collect()
}

fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

fn filter_lines(text: &str, run: &str) -> Vec<String> {
    rows(text, "filter")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// Every record line of one run, rebuilt from the input's own columns.
fn ours(
    text: &str,
    resource: &str,
    snp_tranches: &[f64],
    indel_tranches: &[f64],
    invalidate_previous: bool,
) -> Vec<String> {
    let variants = records(&input(text, "variants"));
    let resources = records(&input(text, resource));
    let snp_tranches = validate_tranches(snp_tranches).expect("valid");
    let indel_tranches = validate_tranches(indel_tranches).expect("valid");

    // The first pass: a record's own score, once, if it overlaps a resource.
    let (mut snp_scores, mut indel_scores) = (Vec::new(), Vec::new());
    let (mut scored_snps, mut scored_indels) = (0usize, 0usize);
    for variant in &variants {
        let Some(value) = variant.score else { continue };
        if variant.is_snp() {
            scored_snps += 1;
        } else {
            scored_indels += 1;
        }
        let overlaps = resources.iter().any(|resource| {
            resource.start == variant.start
                && resource.reference == variant.reference
                && variant
                    .alternates
                    .iter()
                    .any(|alt| resource.alternates.contains(alt))
        });
        if overlaps {
            if variant.is_snp() {
                snp_scores.push(value);
            } else {
                indel_scores.push(value);
            }
        }
    }

    let (snp_cutoffs, indel_cutoffs) = cutoffs(
        &snp_scores,
        &indel_scores,
        scored_snps,
        scored_indels,
        &snp_tranches,
        &indel_tranches,
        "SCORE",
    )
    .expect("this run reaches the second pass");

    variants
        .iter()
        .map(|variant| {
            let added = variant.score.and_then(|value| {
                let (tranches, cutoffs, class) = if variant.is_snp() {
                    (&snp_tranches, &snp_cutoffs, SNP_STRING)
                } else {
                    (&indel_tranches, &indel_cutoffs, INDEL_STRING)
                };
                is_tranche_filtered(value, cutoffs)
                    .then(|| filter_string_from_score("SCORE", class, value, tranches, cutoffs))
            });
            let mut filters = filters_of(&variant.filters, invalidate_previous, added);
            // htsjdk writes a record's filters sorted.
            filters.sort();
            let mut field: Vec<String> = variant.line.split('\t').map(|f| f.to_string()).collect();
            field[6] = filters.join(";");
            field.join("\t")
        })
        .collect()
}

#[test]
fn every_written_record_matches_the_golden_byte_for_byte() {
    let text = golden();
    for (run, snp, indel, invalidate) in [
        ("two-tranches", &[50.0, 99.0][..], &[50.0][..], false),
        (
            "unsorted-and-repeated",
            &[99.0, 50.0, 99.0][..],
            &[50.0][..],
            false,
        ),
        (
            "invalidate-previous-filters",
            &[50.0, 99.0][..],
            &[50.0][..],
            true,
        ),
        ("one-tranche", &[50.0][..], &[50.0][..], false),
    ] {
        assert_eq!(
            ours(&text, "resource", snp, indel, invalidate),
            written(&text, run),
            "{run}"
        );
    }
}

#[test]
fn a_score_on_the_cutoff_is_filtered_and_the_band_comes_after() {
    let text = golden();
    let lines = written(&text, "two-tranches");
    let filter_of = |start: &str| -> String {
        lines
            .iter()
            .find(|line| line.starts_with(&format!("chr1\t{start}\t")))
            .expect("a record")
            .split('\t')
            .nth(6)
            .expect("a FILTER")
            .to_string()
    };
    // Five SNP scores 5..1: the cutoffs are 3.0 and 2.0.
    assert_eq!(filter_of("200"), "PASS");
    assert_eq!(filter_of("300"), "SCORE_SNP_Tranche_50.00_99.00");
    assert_eq!(filter_of("400"), "SCORE_SNP_Tranche_99.00_100.00");
    // Three indel scores 9, 8, 7: the cutoff is 8.0, and the record scoring exactly 8.0 is filtered.
    assert_eq!(filter_of("600"), "PASS");
    assert_eq!(filter_of("700"), "SCORE_INDEL_Tranche_50.00_100.00");
}

#[test]
fn the_tranche_list_is_deduplicated_and_sorted() {
    let text = golden();
    // The same two tranches out of order and one of them twice: the same output, record for record.
    assert_eq!(
        written(&text, "two-tranches"),
        written(&text, "unsorted-and-repeated")
    );
    assert_eq!(
        validate_tranches(&[99.0, 50.0, 99.0]).expect("valid"),
        vec![50.0, 99.0]
    );
}

#[test]
fn a_record_with_no_score_is_written_and_passes() {
    let text = golden();
    let line = written(&text, "two-tranches")
        .into_iter()
        .find(|line| line.starts_with("chr1\t900\t"))
        .expect("the unscored record");
    assert_eq!(line.split('\t').nth(6), Some("PASS"));
    assert_eq!(line.split('\t').nth(7), Some("."));
}

#[test]
fn invalidating_previous_filters_clears_the_header_too() {
    let text = golden();
    assert!(filter_lines(&text, "two-tranches")
        .iter()
        .any(|line| line.contains("ID=weak,")));
    assert!(!filter_lines(&text, "invalidate-previous-filters")
        .iter()
        .any(|line| line.contains("ID=weak,")));
    // And the record that arrived filtered keeps only what the tool gave it.
    let line = written(&text, "invalidate-previous-filters")
        .into_iter()
        .find(|line| line.starts_with("chr1\t500\t"))
        .expect("the filtered record");
    assert_eq!(
        line.split('\t').nth(6),
        Some("SCORE_SNP_Tranche_99.00_100.00")
    );
}

#[test]
fn the_header_lines_are_the_ports_own() {
    let text = golden();
    let mut ours: Vec<String> = tranche_header_lines("SCORE", SNP_STRING, &[50.0, 99.0])
        .into_iter()
        .chain(tranche_header_lines("SCORE", INDEL_STRING, &[50.0]))
        .map(|(id, description)| format!("##FILTER=<ID={id},Description=\"{description}\">"))
        .collect();
    ours.sort();
    let mut theirs: Vec<String> = filter_lines(&text, "invalidate-previous-filters");
    theirs.sort();
    assert_eq!(ours, theirs);
}

#[test]
fn every_refusal_carries_the_references_class_and_words() {
    let text = golden();
    let refusal = |label: &str| -> (String, String) {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let (class, message) = row[1].split_once(':').expect("class and message");
        (class.to_string(), unescape(message))
    };

    let none: [f64; 0] = [];
    for (label, error) in [
        (
            "tranche-of-a-hundred",
            validate_tranches(&[100.0]).expect_err("a hundred"),
        ),
        (
            "nothing-scored",
            cutoffs(&none, &none, 0, 0, &[50.0], &[50.0], "SCORE").expect_err("nothing scored"),
        ),
        (
            "no-overlap",
            cutoffs(&none, &none, 5, 3, &[50.0], &[50.0], "SCORE").expect_err("no overlap"),
        ),
        (
            "snps-without-snp-resources",
            cutoffs(&none, &[9.0], 5, 3, &[50.0], &[50.0], "SCORE").expect_err("indels only"),
        ),
    ] {
        let (class, message) = refusal(label);
        assert_eq!(error.class(), class, "{label}");
        assert_eq!(error.message(), message, "{label}");
    }

    // And the one that is checked before the traversal starts.
    let (class, message) = refusal("info-key-not-in-header");
    let error = FilterVariantTranchesError::InfoKeyNotInHeader("MISSING".to_string());
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
}
