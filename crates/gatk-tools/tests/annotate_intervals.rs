//! Conformance for `AnnotateIntervals` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/AnnotateIntervalsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the BED off-by-one**, which is the only reason a track covering an interval exactly
//!    annotates eleven twelfths rather than one;
//!  * **the GC denominator**, which is the `ACGT` count and not the interval's length;
//!  * **`NaN` for an interval with no `ACGT` at all**, written as three characters;
//!  * **a missing score reading as one**;
//!  * **and the overlap check**, which accepts touching features and refuses overlapping ones.

use gatk_corpus as corpus;
use gatk_tools::annotate_intervals::{
    self, AnnotateError, BedFeature, GC_CONTENT, MAPPABILITY, SEGMENTAL_DUPLICATION_CONTENT,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/annotate_intervals.txt.gz"),
    )
}

/// The data rows of one case's table, with the SAM header and the column line dropped.
fn rows(text: &str, label: &str) -> Vec<String> {
    let prefix = format!("table\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {label}"))
        .split("\\n")
        .filter(|row| !row.is_empty())
        .map(|row| row.replace("\\t", "\t"))
        // Three `@` lines then the column names.
        .skip(4)
        .collect()
}

fn error(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"))
        .to_string()
}

/// `ReferenceQueryDump.FASTA`'s chr1, as the caching reader answers it: upper-cased, IUPAC to `N`.
const CHR1: &str = "ACGTACGTACGTACGTNNNNACGTACGTNNNNNNNNNNNACGT";

fn bases(start: i32, end: i32) -> &'static [u8] {
    &CHR1.as_bytes()[(start - 1) as usize..end as usize]
}

/// The dump's mappability track.
fn mappability() -> Vec<BedFeature> {
    vec![
        BedFeature {
            contig: "chr1".to_string(),
            start: 0,
            end: 12,
            score: 1.0,
        },
        BedFeature {
            contig: "chr1".to_string(),
            start: 12,
            end: 18,
            score: 0.5,
        },
        // No score column at all, which arrives as NaN and is read as one.
        BedFeature {
            contig: "chr1".to_string(),
            start: 24,
            end: 30,
            score: f64::NAN,
        },
    ]
}

fn overlapping_for(features: &[BedFeature], start: i32, end: i32) -> Vec<&BedFeature> {
    features
        .iter()
        .filter(|feature| feature.contig == "chr1" && feature.start < end && start <= feature.end)
        .collect()
}

#[test]
fn the_gc_rows_match_the_golden() {
    let text = golden();
    let expected = rows(&text, "gc-only");
    let mine: Vec<String> = [(1, 12), (13, 24), (25, 36)]
        .iter()
        .map(|(start, end)| {
            annotate_intervals::row(
                "chr1",
                *start,
                *end,
                &[annotate_intervals::gc_content(bases(*start, *end))],
            )
        })
        .collect();
    assert_eq!(mine, expected);
}

/// The IUPAC stretch is `ACGTNNNNNNNN` through the reader, so the denominator is four.
#[test]
fn the_gc_denominator_is_the_acgt_count() {
    assert_eq!(annotate_intervals::gc_content(b"ACGTNNNNNNNN"), 0.5);
    assert_eq!(annotate_intervals::gc_content(b"GGGG"), 1.0);
    assert_eq!(annotate_intervals::gc_content(b"AAAA"), 0.0);
}

#[test]
fn an_interval_with_no_acgt_is_nan() {
    let text = golden();
    let value = annotate_intervals::gc_content(bases(29, 32));
    assert!(value.is_nan());
    assert_eq!(
        annotate_intervals::row("chr1", 29, 32, &[value]),
        rows(&text, "all-n")[0]
    );
}

#[test]
fn the_single_base_rows_match_the_golden() {
    let text = golden();
    let expected = rows(&text, "single-base");
    let mine: Vec<String> = [1, 2]
        .iter()
        .map(|position| {
            annotate_intervals::row(
                "chr1",
                *position,
                *position,
                &[annotate_intervals::gc_content(bases(*position, *position))],
            )
        })
        .collect();
    assert_eq!(mine, expected);
}

/// The off-by-one: a track covering `0-12` over `1-12` annotates eleven twelfths.
#[test]
fn the_mappability_rows_match_the_golden() {
    let text = golden();
    let track = mappability();
    let expected = rows(&text, "mappability");
    let mine: Vec<String> = [(1, 12), (13, 24), (25, 30)]
        .iter()
        .map(|(start, end)| {
            let features = overlapping_for(&track, *start, *end);
            annotate_intervals::row(
                "chr1",
                *start,
                *end,
                &[
                    annotate_intervals::gc_content(bases(*start, *end)),
                    annotate_intervals::length_weighted_annotation(&features, *start, *end),
                ],
            )
        })
        .collect();
    assert_eq!(mine, expected);
    assert!(
        expected[0].ends_with("0.916667"),
        "eleven twelfths, not one: {}",
        expected[0]
    );
}

/// A missing score is one, which is why the third interval annotates above a half.
#[test]
fn a_missing_score_is_one() {
    let track = mappability();
    let scored =
        annotate_intervals::length_weighted_annotation(&overlapping_for(&track, 25, 30), 25, 30);
    assert!(scored > 0.5, "{scored}");
}

/// Two touching features are two after the merge; two overlapping ones are one.
#[test]
fn the_overlap_check_accepts_touching_and_refuses_overlapping() {
    let touching = vec![
        BedFeature {
            contig: "chr1".to_string(),
            start: 0,
            end: 12,
            score: 1.0,
        },
        BedFeature {
            contig: "chr1".to_string(),
            start: 12,
            end: 18,
            score: 1.0,
        },
    ];
    assert!(!annotate_intervals::track_has_overlaps(&touching));

    let overlapping = vec![
        BedFeature {
            contig: "chr1".to_string(),
            start: 0,
            end: 12,
            score: 1.0,
        },
        BedFeature {
            contig: "chr1".to_string(),
            start: 6,
            end: 18,
            score: 1.0,
        },
    ];
    assert!(annotate_intervals::track_has_overlaps(&overlapping));

    let text = golden();
    let refusal =
        AnnotateError::OverlappingTrack("/work/annotateintervals-dump/overlapping.bed".to_string());
    assert_eq!(
        format!("{}:{}", refusal.java_class(), refusal.message()),
        error(&text, "overlapping-track")
    );
}

/// The merging rule this tool insists on.
#[test]
fn the_default_merging_rule_is_refused() {
    let text = golden();
    let refusal = AnnotateError::MergingRuleNotOverlappingOnly;
    assert_eq!(
        format!("{}:{}", refusal.java_class(), refusal.message()),
        error(&text, "default-merging")
    );
}

/// The columns follow the annotators that ran, in the order they were added.
#[test]
fn the_columns_follow_the_annotators() {
    assert_eq!(
        annotate_intervals::columns(&[GC_CONTENT]),
        "CONTIG\tSTART\tEND\tGC_CONTENT"
    );
    assert_eq!(
        annotate_intervals::columns(&[GC_CONTENT, MAPPABILITY, SEGMENTAL_DUPLICATION_CONTENT]),
        "CONTIG\tSTART\tEND\tGC_CONTENT\tMAPPABILITY\tSEGMENTAL_DUPLICATION_CONTENT"
    );
    let text = golden();
    assert!(rows(&text, "both-tracks")[0].split('\t').count() == 6);
}
