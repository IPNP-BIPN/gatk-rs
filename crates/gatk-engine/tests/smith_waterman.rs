//! Conformance for `SmithWatermanJavaAligner` and `CigarUtils.calculateCigar` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SmithWatermanDump.java`, which ran every Smith-Waterman
//! implementation the pinned container could build over 46 reference and alternate pairs, under the
//! four `SWParameters` sets and the four `SWOverhangStrategy` values, and printed the alignment
//! offset and CIGAR of each; and then `CigarUtils.calculateCigar` over the same pairs with the two
//! strategies GATK passes it. The port targets the pure-Java aligner, which is the one the oracle
//! contract pins.
//!
//! # What this suite is for
//!
//!  * **every offset and CIGAR of the matrix, including the tie-breaks of the traceback**: the
//!    homopolymer and tandem-repeat indels whose placement is ambiguous;
//!  * **the exact-match shortcut, which takes the last occurrence and only under two strategies**;
//!  * **the overhang strategies, on an alternate longer than the reference**;
//!  * **the two empty inputs, which the Java refuses**;
//!  * **and `calculateCigar`'s padding, trimming and left-alignment, with the leading and trailing
//!    deletions it puts back.**

use gatk_corpus as corpus;
use gatk_engine::alignment_utils::AlignmentError;
use gatk_engine::cigar_utils::calculate_cigar;
use gatk_engine::smith_waterman::{
    SmithWatermanAligner, SmithWatermanJavaAligner, SwOverhangStrategy, SwParameters,
    ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS, NEW_SW_PARAMETERS, ORIGINAL_DEFAULT, STANDARD_NGS,
};

const JAVA: &str = "SmithWatermanJavaAligner";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/smith_waterman.txt.gz"),
    )
}

/// The reference's exception class, without its package.
fn exception(error: &AlignmentError) -> String {
    format!(
        "exception:{}",
        error.class().rsplit('.').next().expect("a class name")
    )
}

const PAIRS: &[(&str, &str, &str)] = &[
    (
        "identical",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
    ),
    (
        "sub-middle",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAAAGTCCAGTTGACCATGCAAGT",
    ),
    (
        "sub-first-base",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "CATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
    ),
    (
        "sub-last-base",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGA",
    ),
    (
        "two-adjacent-subs",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGGTAGT",
    ),
    (
        "two-distant-subs",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACTGCCTAGGCTTAACGTCCAGTTGACCATGCAAGC",
    ),
    (
        "insertion-middle",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCGGGAGTTGACCATGCAAGT",
    ),
    (
        "deletion-middle",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGCAGTTGACCATGCAAGT",
    ),
    (
        "insertion-near-start",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTTTTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
    ),
    (
        "insertion-near-end",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGTTTTGT",
    ),
    (
        "deletion-near-start",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
    ),
    (
        "deletion-near-end",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCA",
    ),
    (
        "insertion-and-deletion",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGGGCTTAACGTCCAGTGACCATGCAAGT",
    ),
    (
        "sub-then-deletion",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGTCGAGTGACCATGCAAGT",
    ),
    (
        "homopolymer-insertion",
        "CGTCAAAAAAAGTCGATCGTAGCTA",
        "CGTCAAAAAAAAGTCGATCGTAGCTA",
    ),
    (
        "homopolymer-deletion",
        "CGTCAAAAAAAGTCGATCGTAGCTA",
        "CGTCAAAAAAGTCGATCGTAGCTA",
    ),
    (
        "homopolymer-deletion-2",
        "CGTCAAAAAAAGTCGATCGTAGCTA",
        "CGTCAAAAAGTCGATCGTAGCTA",
    ),
    (
        "tandem-insertion",
        "TTGCACACACACACGGATCCAGTT",
        "TTGCACACACACACACGGATCCAGTT",
    ),
    (
        "tandem-deletion",
        "TTGCACACACACACGGATCCAGTT",
        "TTGCACACACGGATCCAGTT",
    ),
    (
        "tandem-trinucleotide-deletion",
        "ATGCAGCAGCAGCAGCAGTTGCA",
        "ATGCAGCAGCAGTTGCA",
    ),
    (
        "tandem-trinucleotide-insertion",
        "ATGCAGCAGCAGTTGCA",
        "ATGCAGCAGCAGCAGCAGTTGCA",
    ),
    (
        "homopolymer-at-start-deletion",
        "AAAAAGTCCGATTGCA",
        "AAAAGTCCGATTGCA",
    ),
    (
        "homopolymer-at-end-insertion",
        "GTCCGATTGCAAAAA",
        "GTCCGATTGCAAAAAA",
    ),
    (
        "alt-longer-both-ends",
        "ACGTCCAGTTGACCAT",
        "TTTTTACGTCCAGTTGACCATGGGGG",
    ),
    (
        "alt-overhangs-left",
        "ACGTCCAGTTGACCAT",
        "GGGGGACGTCCAGTTGACCAT",
    ),
    (
        "alt-overhangs-right",
        "ACGTCCAGTTGACCAT",
        "ACGTCCAGTTGACCATGGGGG",
    ),
    (
        "alt-overhangs-left-mismatch",
        "ACGTCCAGTTGACCAT",
        "GGGGGACGTCCAGATGACCAT",
    ),
    (
        "alt-overhangs-right-mismatch",
        "ACGTCCAGTTGACCAT",
        "ACGTCCAGTTGACGATGGGGG",
    ),
    (
        "alt-is-prefix",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "GATTACAGCCTAGGCTTAACGT",
    ),
    (
        "alt-is-suffix",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "TTAACGTCCAGTTGACCATGCAAGT",
    ),
    (
        "alt-is-interior",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "CCTAGGCTTAACGTCCAGTTGA",
    ),
    (
        "alt-occurs-twice",
        "TTGACCAGTCATTGACCAGTCATTTGAC",
        "TTGACCAGTCA",
    ),
    (
        "alt-occurs-twice-with-mismatch",
        "TTGACCAGTCATTGACCAGTCATTTGAC",
        "TTGACCAGACA",
    ),
    ("one-base-same", "A", "A"),
    ("one-base-different", "A", "C"),
    ("one-base-ref-long-alt", "A", "ACGTACGT"),
    (
        "long-ref-one-base-alt",
        "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT",
        "T",
    ),
    ("two-bases-swapped", "ACGT", "AGCT"),
    ("nothing-in-common", "AAAAAAAAAAAA", "CCCCCCCCCCCC"),
    ("poly-n", "NNNNNNNNNN", "NNNNNNNN"),
    ("empty-reference", "", "ACGT"),
    ("empty-alternate", "ACGT", ""),
    (
        "haplotype-deletion-and-snv",
        "CTGAACGTTAGCCATGCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGATCGTAC",
        "CTGAACGTTAGCCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGTTCGTAC",
    ),
    (
        "haplotype-insertion-and-snv",
        "CTGAACGTTAGCCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGATCGTAC",
        "CTGAACGTTAGCCATGCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGTTCGTAC",
    ),
    (
        "haplotype-leading-deletion",
        "AAAAAAGTCCGATTGCAGGCTAAC",
        "AAAAGTCCGATTGCAGGCTAAC",
    ),
    (
        "haplotype-trailing-deletion",
        "GTCCGATTGCAGGCTAACTTTTTT",
        "GTCCGATTGCAGGCTAACTTTT",
    ),
];

const PARAMETERS: [(&str, SwParameters); 4] = [
    ("ORIGINAL_DEFAULT", ORIGINAL_DEFAULT),
    ("STANDARD_NGS", STANDARD_NGS),
    ("NEW_SW_PARAMETERS", NEW_SW_PARAMETERS),
    (
        "ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS",
        ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS,
    ),
];

const STRATEGIES: [(&str, SwOverhangStrategy); 4] = [
    ("SOFTCLIP", SwOverhangStrategy::SoftClip),
    ("INDEL", SwOverhangStrategy::Indel),
    ("LEADING_INDEL", SwOverhangStrategy::LeadingIndel),
    ("IGNORE", SwOverhangStrategy::Ignore),
];

/// The golden's row for the Java aligner: everything after the implementation column.
fn java_row<'a>(text: &'a str, kind: &str, key: &[&str]) -> &'a str {
    let prefix = format!("{kind}\t{}\t{JAVA}\t", key.join("\t"));
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("no {kind} row for {key:?}"))
}

#[test]
fn align_matches_the_golden_on_every_pair_parameter_set_and_strategy() {
    let text = golden();
    let mut compared = 0;
    for (label, reference, alternate) in PAIRS {
        for (parameter_name, parameters) in &PARAMETERS {
            for (strategy_name, strategy) in &STRATEGIES {
                let expected = java_row(&text, "align", &[label, parameter_name, strategy_name]);
                let ours = match SmithWatermanJavaAligner.align(
                    reference.as_bytes(),
                    alternate.as_bytes(),
                    parameters,
                    *strategy,
                ) {
                    Ok(alignment) => {
                        format!(
                            "{}\t{}",
                            alignment.alignment_offset,
                            alignment.cigar.to_text()
                        )
                    }
                    Err(error) => format!("-\t{}", exception(&error)),
                };
                assert_eq!(ours, expected, "{label} {parameter_name} {strategy_name}");
                compared += 1;
            }
        }
    }
    assert_eq!(compared, PAIRS.len() * 16);
}

#[test]
fn calculate_cigar_matches_the_golden_on_every_pair() {
    let text = golden();
    let mut compared = 0;
    for (label, reference, alternate) in PAIRS {
        for (parameter_name, parameters) in &PARAMETERS {
            for (strategy_name, strategy) in &STRATEGIES[..2] {
                let expected = java_row(&text, "cigar", &[label, parameter_name, strategy_name]);
                let ours = match calculate_cigar(
                    reference.as_bytes(),
                    alternate.as_bytes(),
                    &SmithWatermanJavaAligner,
                    parameters,
                    *strategy,
                ) {
                    Ok(Some(cigar)) => cigar.to_text(),
                    Ok(None) => "null".to_string(),
                    Err(error) => exception(&error),
                };
                assert_eq!(ours, expected, "{label} {parameter_name} {strategy_name}");
                compared += 1;
            }
        }
    }
    assert_eq!(compared, PAIRS.len() * 8);
}

#[test]
fn the_golden_pins_the_java_aligner() {
    let text = golden();
    assert!(text.contains(&format!("loaded\t{JAVA}\tyes\n")));
    // The Intel aligner never loads in the pinned image, so `FASTEST_AVAILABLE` falls back to the
    // Java one. A runner that one day loads it would add rows for it and fail the manifest first.
    assert!(text.contains("loaded\tSmithWatermanIntelAligner\tno:HardwareFeatureException\n"));
    assert!(text.contains(&format!("fastest\t{JAVA}\n")));
}
