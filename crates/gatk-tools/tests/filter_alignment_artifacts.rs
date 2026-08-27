//! Conformance for `FilterAlignmentArtifacts` against GATK 4.6.2.0, compared as the answer of
//! every `supportsVariant` case and every `findJointAlignments` case.
//!
//! Golden from `tools/readfilter-conformance/FilterAlignmentArtifactsDump.java`.
//!
//! The realignment itself is BWA's and the assembly that makes the unitigs is the Mutect2
//! assembler's; neither is measured or ported. What is compared here is the two rules the tool
//! owns, and the filter decision they feed.
//!
//! # What this suite is for
//!
//!  * **a SNP being matched by bases and an indel by cigar operator**;
//!  * **a variant position inside a deletion never supporting a SNP**;
//!  * **the running sum having to land on the offset exactly at a tolerance of zero**;
//!  * **an insertion never being supported by an `I` at that tolerance**, only by a clip;
//!  * **the tolerance loop freezing the sum once one element qualifies**;
//!  * **two unitigs joining only on the same strand and within the padding**;
//!  * **the best-scoring alignment of each unitig being the one kept**;
//!  * **and the filter needing BOTH differences to fall below their thresholds.**

use gatk_corpus as corpus;
use gatk_tools::filter_alignment_artifacts::{
    decide, find_joint_alignments, read_index_for_reference_coordinate, supports_variant,
    Alignment, Read, Variant, ALIGNMENT_ARTIFACT_FILTER_NAME,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_alignment_artifacts.txt.gz"),
    )
}

/// One `support` line of the golden.
fn support(text: &str, label: &str) -> bool {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("support\t{label}=")))
        .unwrap_or_else(|| panic!("the golden carries support/{label}"))
        == "true"
}

/// One `joint` line of the golden, as the sorted groups the dump wrote.
fn joint(text: &str, label: &str) -> Vec<Vec<String>> {
    let body = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("joint\t{label}=")))
        .unwrap_or_else(|| panic!("the golden carries joint/{label}"));
    if body.is_empty() {
        return Vec::new();
    }
    body.split(';')
        .map(|group| group.split(',').map(str::to_string).collect())
        .collect()
}

fn parse_cigar(text: &str) -> Vec<(char, i32)> {
    let mut out = Vec::new();
    let mut digits = String::new();
    for character in text.chars() {
        if character.is_ascii_digit() {
            digits.push(character);
        } else {
            out.push((character, digits.parse().expect("a length")));
            digits.clear();
        }
    }
    out
}

fn read(start: i32, cigar: &str, bases: &str) -> Read {
    Read {
        name: "r".to_string(),
        start,
        cigar: parse_cigar(cigar),
        bases: bases.as_bytes().to_vec(),
    }
}

fn variant(start: i32, reference: &str, alternate: &str) -> Variant {
    Variant {
        start,
        reference: reference.as_bytes().to_vec(),
        alternates: vec![alternate.as_bytes().to_vec()],
    }
}

/// The dump's fifteen cases, restated here exactly as it built them.
fn cases() -> Vec<(&'static str, Read, Variant, i32)> {
    let twenty = "ACGTACGTACGTACGTACGT";
    let twenty_four = "ACGTACGTACGTACGTACGTACGT";
    vec![
        (
            "snp-matching",
            read(1000, "20M", twenty),
            variant(1005, "A", "C"),
            0,
        ),
        (
            "snp-not-matching",
            read(1000, "20M", twenty),
            variant(1005, "A", "G"),
            0,
        ),
        (
            "snp-not-covered",
            read(1000, "3M", "ACG"),
            variant(1005, "A", "C"),
            0,
        ),
        (
            "snp-in-deletion",
            read(1000, "3M5D12M", "ACGTACGTACGTACG"),
            variant(1005, "A", "C"),
            0,
        ),
        (
            "del-exact-deletion",
            read(1000, "5M4D15M", twenty),
            variant(1005, "ACGTA", "A"),
            0,
        ),
        (
            "del-exact-clip",
            read(1000, "5M15S", twenty),
            variant(1005, "ACGTA", "A"),
            0,
        ),
        (
            "ins-exact-clip",
            read(1000, "5M15S", twenty),
            variant(1005, "A", "ACGTA"),
            0,
        ),
        (
            "ins-exact-insertion",
            read(1000, "5M4I15M", twenty_four),
            variant(1005, "A", "ACGTA"),
            0,
        ),
        (
            "ins-tolerant-insertion",
            read(1000, "5M4I15M", twenty_four),
            variant(1005, "A", "ACGTA"),
            100,
        ),
        (
            "ins-exact-deletion",
            read(1000, "5M4D15M", twenty),
            variant(1005, "A", "ACGTA"),
            0,
        ),
        (
            "del-exact-insertion",
            read(1000, "5M4I15M", twenty_four),
            variant(1005, "ACGTA", "A"),
            0,
        ),
        (
            "del-stepped-over",
            read(1000, "6M4D14M", twenty),
            variant(1005, "ACGTA", "A"),
            0,
        ),
        (
            "tolerance-0",
            read(1000, "12M4D8M", twenty),
            variant(1005, "ACGTA", "A"),
            0,
        ),
        (
            "tolerance-5",
            read(1000, "12M4D8M", twenty),
            variant(1005, "ACGTA", "A"),
            5,
        ),
        (
            "tolerance-100",
            read(1000, "12M4D8M", twenty),
            variant(1005, "ACGTA", "A"),
            100,
        ),
        (
            "tolerance-100-insertion",
            read(1000, "12M4D8M", twenty),
            variant(1005, "A", "ACGTA"),
            100,
        ),
    ]
}

fn alignment(
    reference_id: i32,
    reference_start: i32,
    reference_end: i32,
    score: i32,
    mismatches: i32,
    reverse_strand: bool,
) -> Alignment {
    Alignment {
        reference_id,
        reference_start,
        reference_end,
        score,
        mismatches,
        reverse_strand,
    }
}

fn describe(alignment: &Alignment) -> String {
    format!(
        "{}:{}-{}/{}/{}/{}",
        alignment.reference_id,
        alignment.reference_start,
        alignment.reference_end,
        alignment.score,
        alignment.mismatches,
        if alignment.reverse_strand { "-" } else { "+" }
    )
}

/// The dump's thirteen joint cases, restated as it built them.
fn joint_cases() -> Vec<(&'static str, Vec<Vec<Alignment>>, i32)> {
    vec![
        ("no-unitigs", vec![], 1000),
        (
            "one-unitig",
            vec![vec![
                alignment(0, 1000, 1100, 100, 0, false),
                alignment(0, 5000, 5100, 90, 2, false),
                alignment(1, 1000, 1100, 80, 4, false),
            ]],
            1000,
        ),
        (
            "two-same-strand",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 1050, 1150, 95, 1, false)],
            ],
            1000,
        ),
        (
            "two-opposite-strands",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 1050, 1150, 95, 1, true)],
            ],
            1000,
        ),
        (
            "far-apart-narrow",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 3000, 3100, 95, 1, false)],
            ],
            1000,
        ),
        (
            "far-apart-wide",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 3000, 3100, 95, 1, false)],
            ],
            100_000,
        ),
        (
            "best-score-kept",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![
                    alignment(0, 1050, 1150, 95, 1, false),
                    alignment(0, 1060, 1160, 120, 3, false),
                ],
            ],
            1000,
        ),
        (
            "two-loci",
            vec![
                vec![
                    alignment(0, 1000, 1100, 100, 0, false),
                    alignment(0, 50000, 50100, 98, 1, false),
                ],
                vec![
                    alignment(0, 1050, 1150, 95, 1, false),
                    alignment(0, 50050, 50150, 93, 2, false),
                ],
            ],
            1000,
        ),
        (
            "one-sided",
            vec![
                vec![
                    alignment(0, 1000, 1100, 100, 0, false),
                    alignment(0, 90000, 90100, 99, 0, false),
                ],
                vec![alignment(0, 1050, 1150, 95, 1, false)],
            ],
            1000,
        ),
        (
            "different-contigs",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(1, 1000, 1100, 95, 1, false)],
            ],
            1000,
        ),
        (
            "three-unitigs",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 1050, 1150, 95, 1, false)],
                vec![alignment(0, 1080, 1180, 90, 2, false)],
            ],
            1000,
        ),
        (
            "three-unitigs-one-missing",
            vec![
                vec![alignment(0, 1000, 1100, 100, 0, false)],
                vec![alignment(0, 1050, 1150, 95, 1, false)],
                vec![alignment(0, 90000, 90100, 90, 2, false)],
            ],
            1000,
        ),
    ]
}

/// The dump sorts both the groups and the alignments inside them, because the reference's own
/// order comes out of a HashSet over identity hashes.
fn sorted(groups: Vec<Vec<Alignment>>) -> Vec<Vec<String>> {
    let mut out: Vec<Vec<String>> = groups
        .into_iter()
        .map(|group| {
            let mut parts: Vec<String> = group.iter().map(describe).collect();
            parts.sort();
            parts
        })
        .collect();
    out.sort();
    out
}

#[test]
fn every_case_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, read, variant, tolerance) in cases() {
        assert_eq!(
            supports_variant(&read, &variant, tolerance),
            support(&text, label),
            "{label}"
        );
        compared += 1;
    }
    for (label, unitigs, max_fragment) in joint_cases() {
        assert_eq!(
            sorted(find_joint_alignments(&unitigs, max_fragment)),
            joint(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 16 + 12, "the cases the golden carries");
}

/// The bases are compared from the offset the coordinate maps to, and a coordinate that fell in a
/// deletion is refused before the bases are looked at.
#[test]
fn a_snp_is_matched_by_bases_unless_the_position_was_deleted() {
    let text = golden();
    assert!(support(&text, "snp-matching"));
    assert!(!support(&text, "snp-not-matching"));
    assert!(!support(&text, "snp-not-covered"));
    assert!(!support(&text, "snp-in-deletion"));

    // The deleted case really does find an offset: the coordinate is bracketed by the deletion,
    // which consumes no read bases, so the index is the last one BEFORE it and the operator is
    // `D`. Only that operator refuses it.
    let deleted = read(1000, "3M5D12M", "ACGTACGTACGTACG");
    let index = read_index_for_reference_coordinate(&deleted, 1005).expect("an index");
    assert_eq!(index.operator, 'D');
    assert_eq!(index.offset, 3, "the last base before the deletion");
    // A SNP whose alternate IS the base at that offset is still refused, so it is the operator and
    // nothing else that decides.
    assert_eq!(deleted.bases[index.offset], b'T');
    assert!(!supports_variant(&deleted, &variant(1005, "A", "T"), 0));
    // A read with no deletion at the same offset does support a SNP.
    assert!(supports_variant(
        &read(1000, "20M", "ACGTACGTACGTACGTACGT"),
        &variant(1005, "A", "C"),
        0
    ));
}

/// The sum has to land on the offset exactly, and an insertion's own bases push the offset past
/// where the sum can reach.
#[test]
fn an_insertion_is_never_supported_by_an_insertion_at_tolerance_zero() {
    let text = golden();
    assert!(support(&text, "del-exact-deletion"), "a deletion is");
    assert!(support(&text, "ins-exact-clip"), "and a clip supports one");
    assert!(
        !support(&text, "ins-exact-insertion"),
        "but an I cannot be reached"
    );

    // The reason: the coordinate maps PAST the inserted bases, so the offset is 9 while the sum
    // reaches 5 and then 9 only after adding the insertion's own length.
    let inserted = read(1000, "5M4I15M", "ACGTACGTACGTACGTACGTACGT");
    let index = read_index_for_reference_coordinate(&inserted, 1005).expect("an index");
    assert_eq!(index.offset, 9);
    // The deletion's coordinate, by contrast, maps to 5, which is where the sum lands.
    let deleted = read(1000, "5M4D15M", "ACGTACGTACGTACGTACGT");
    assert_eq!(
        read_index_for_reference_coordinate(&deleted, 1005)
            .expect("an index")
            .offset,
        5
    );
    // Widening the tolerance is what reaches the I, and it does so by freezing the sum at 0.
    assert!(support(&text, "ins-tolerant-insertion"));
}

/// A sum that steps over the offset never lands on it, and the loop stops advancing once one
/// element is within tolerance.
#[test]
fn the_tolerance_loop_freezes_the_sum() {
    let text = golden();
    assert!(
        !support(&text, "del-stepped-over"),
        "6M steps 0 to 6 over 5"
    );
    assert!(!support(&text, "tolerance-0"));
    assert!(support(&text, "tolerance-5"));
    assert!(support(&text, "tolerance-100"));
    // The freeze is not a licence to match anything: a D still cannot support an insertion at any
    // tolerance.
    assert!(!support(&text, "tolerance-100-insertion"));

    // The read at tolerance 5 has its indel at a sum of 12, seven away from the offset of 5. It
    // matches because the first element already qualifies and the sum never moves.
    let far = read(1000, "12M4D8M", "ACGTACGTACGTACGTACGT");
    assert!(!supports_variant(&far, &variant(1005, "ACGTA", "A"), 4));
    assert!(supports_variant(&far, &variant(1005, "ACGTA", "A"), 5));
}

/// Same strand, within the padding, and the best score of each unitig.
#[test]
fn two_unitigs_join_only_on_the_same_strand() {
    let text = golden();
    assert_eq!(joint(&text, "two-same-strand").len(), 1);
    assert!(joint(&text, "two-opposite-strands").is_empty());
    assert!(joint(&text, "different-contigs").is_empty());
    assert!(joint(&text, "far-apart-narrow").is_empty());
    assert_eq!(joint(&text, "far-apart-wide").len(), 1);
    // Three unitigs must ALL reach the locus.
    assert_eq!(joint(&text, "three-unitigs").len(), 1);
    assert!(joint(&text, "three-unitigs-one-missing").is_empty());
    // The higher-scoring of two candidates is the one kept, even though it is further away.
    assert_eq!(
        joint(&text, "best-score-kept"),
        vec![vec![
            "0:1000-1100/100/0/+".to_string(),
            "0:1060-1160/120/3/+".to_string()
        ]]
    );
    // One unitig makes each alignment its own group; none makes none.
    assert_eq!(joint(&text, "one-unitig").len(), 3);
    assert!(joint(&text, "no-unitigs").is_empty());
    // A locus only one unitig reaches is dropped, so two candidates give one group.
    assert_eq!(joint(&text, "one-sided").len(), 1);
    assert_eq!(joint(&text, "two-loci").len(), 2);
}

/// Another contig filters outright; otherwise BOTH differences must fall below their thresholds.
#[test]
fn the_filter_needs_both_differences() {
    let on_chr1 = vec![alignment(0, 1000, 1100, 100, 0, false)];
    let on_chr2 = vec![alignment(1, 1000, 1100, 100, 0, false)];

    // The best joint alignment on another contig, which is decided before any threshold.
    let elsewhere = decide(std::slice::from_ref(&on_chr2), 0, 100, 0.02, 0.02);
    assert!(elsewhere.filtered);
    assert_eq!(elsewhere.score_difference, None);
    assert_eq!(elsewhere.joint_alignment_count, 1);

    // One joint alignment on the right contig: nothing to compare, nothing filtered.
    let single = decide(std::slice::from_ref(&on_chr1), 0, 100, 0.02, 0.02);
    assert!(!single.filtered);
    assert_eq!(single.score_difference, None);

    // Two of them, close in score and in mismatches: a multimapping, so filtered.
    let close = vec![
        on_chr1.clone(),
        vec![alignment(0, 50000, 50100, 99, 0, false)],
    ];
    let decision = decide(&close, 0, 100, 0.02, 0.02);
    assert!(decision.filtered);
    assert_eq!(decision.score_difference, Some(1));
    assert_eq!(decision.joint_alignment_count, 2);

    // A clear score difference alone lifts the filter, even with the mismatches equal.
    let clear_score = vec![
        on_chr1.clone(),
        vec![alignment(0, 50000, 50100, 50, 0, false)],
    ];
    assert!(!decide(&clear_score, 0, 100, 0.02, 0.02).filtered);
    // And a clear MISMATCH difference alone lifts it too, with the scores one apart.
    let clear_mismatches = vec![
        on_chr1.clone(),
        vec![alignment(0, 50000, 50100, 99, 40, false)],
    ];
    let by_mismatches = decide(&clear_mismatches, 0, 100, 0.02, 0.02);
    assert!(!by_mismatches.filtered);
    assert_eq!(by_mismatches.score_difference, Some(1), "still reported");

    assert_eq!(ALIGNMENT_ARTIFACT_FILTER_NAME, "alignment_artifact");
}
