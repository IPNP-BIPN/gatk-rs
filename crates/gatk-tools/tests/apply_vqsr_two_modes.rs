//! Conformance for `ApplyVQSR` run twice against GATK 4.6.2.0, compared as the FILTER the second
//! run gives a mixed site and as what the first run's header tells it.
//!
//! Golden from `tools/readfilter-conformance/ApplyVqsrTwoModesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the leniency regex is greedy and loses a digit**, so a tranche starting at 90 beats one
//!    starting at 5 and the site is filtered with the less lenient name;
//!  * **the prefix test is case-insensitive and the mode letter is not**, so a lowercase tranche
//!    name leaves the run believing no other mode has been applied;
//!  * **an interval that is two non-numbers is a refusal**;
//!  * **and a mixed site is filtered only once both modes have run**.

use gatk_corpus as corpus;
use gatk_engine::tranches::read_tranches;
use gatk_tools::apply_vqsr::{
    keep, parse_filter_lower_limit, previous_runs, site_filter_from_alleles,
    tranche_interval_is_valid, ApplyVqsrError, PreviousRuns,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/apply_vqsr_two_modes.txt.gz"),
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

/// The whole text of one input.
fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

/// The record lines of one run's output.
fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// The `##FILTER` IDs of one run's output, which is what the next run reads.
fn filter_ids(text: &str, run: &str) -> Vec<String> {
    rows(text, "filter")
        .into_iter()
        .filter(|row| row[0] == run)
        .filter_map(|row| {
            unescape(row[1])
                .strip_prefix("##FILTER=<ID=")
                .and_then(|rest| rest.split(',').next())
                .map(|id| id.to_string())
        })
        .collect()
}

/// One record line's FILTER and its `AS_FilterStatus`.
fn filter_and_status(line: &str) -> (String, String) {
    let field: Vec<&str> = line.split('\t').collect();
    let status = field[7]
        .split(';')
        .find_map(|entry| entry.strip_prefix("AS_FilterStatus="))
        .expect("an AS_FilterStatus")
        .to_string();
    (field[6].to_string(), status)
}

#[test]
fn a_mixed_site_is_filtered_only_once_both_modes_have_run() {
    let text = golden();
    let (first, status) = filter_and_status(&written(&text, "first-snp")[0]);
    assert_eq!(first, ".");
    assert_eq!(status, "VQSRTrancheSNP90.00to99.00,NA");
    let (second, status) = filter_and_status(&written(&text, "second-indel")[0]);
    assert_ne!(second, ".");
    // The second run's list is the first run's entry beside its own.
    assert_eq!(
        status,
        "VQSRTrancheSNP90.00to99.00,VQSRTrancheINDEL99.00to100.00"
    );
}

#[test]
fn the_second_run_takes_the_most_lenient_name_of_both_modes() {
    let text = golden();
    let (filter, status) = filter_and_status(&written(&text, "second-indel")[0]);
    // The INDEL run's own answer for a LOD of 0.0 is the widest tranche; the SNP run left a
    // narrower one, and the narrower one is what wins on both readings of the names.
    let tranches =
        read_tranches("tranches-indel", &input(&text, "tranches-indel")).expect("a good file");
    let ours = site_filter_from_alleles(
        true,
        false,
        previous_runs(&filter_ids(&text, "first-snp"))
            .expect("well formed")
            .both_modes_were_run(false),
        Some(&status),
        0.0,
        &keep(&tranches, 0.0),
    );
    assert_eq!(ours, filter);
}

#[test]
fn the_leniency_regex_is_greedy_and_loses_a_digit() {
    let text = golden();
    let (filter, status) = filter_and_status(&written(&text, "second-indel-inverted")[0]);
    // The two names in play, and what the regex makes of them.
    assert_eq!(
        status,
        "VQSRTrancheSNP5.00to90.00,VQSRTrancheINDEL90.00to99.00"
    );
    assert_eq!(parse_filter_lower_limit("VQSRTrancheSNP5.00to90.00"), 5.0);
    assert_eq!(
        parse_filter_lower_limit("VQSRTrancheINDEL90.00to99.00"),
        0.0
    );
    // So the tranche whose interval truly starts at 90 is taken as the more lenient of the two.
    assert_eq!(filter, "VQSRTrancheINDEL90.00to99.00");
}

#[test]
fn the_prefix_test_is_case_insensitive_and_the_mode_letter_is_not() {
    let text = golden();
    // The header of the lowercase run carries a name that looks like a tranche.
    let ids = filter_ids(&text, "lowercase-tranche");
    assert!(ids.iter().any(|id| id == "vqsrtranchesnp90.00to99.00"));
    // And the run behaved as though no SNP mode had ever been applied: the site is unfiltered.
    let (filter, _) = filter_and_status(&written(&text, "lowercase-tranche")[0]);
    assert_eq!(filter, ".");
    assert_eq!(
        previous_runs(&ids).expect("no refusal"),
        PreviousRuns {
            snp: false,
            indel: true
        }
    );
}

#[test]
fn an_interval_of_two_non_numbers_is_the_references_refusal() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "malformed-tranche")
        .expect("a refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");
    let error = previous_runs(&["VQSRTrancheSNPaatobb".to_string()]).expect_err("two non-numbers");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
    assert_eq!(error, ApplyVqsrError::PoorlyFormattedTrancheName);
    // The input header of that run is where the name came from.
    assert!(rows(&text, "input")
        .into_iter()
        .any(|row| row[0] == "malformed-tranche" && row[1].contains("VQSRTrancheSNPaatobb")));
    assert!(tranche_interval_is_valid("90.00to99.00").expect("a pair of numbers"));
}
