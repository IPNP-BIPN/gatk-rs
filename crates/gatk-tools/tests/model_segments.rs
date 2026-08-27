//! Conformance for `ModelSegments` against GATK 4.6.2.0, compared as the segmentation each run
//! writes and the genotyping that decides it.
//!
//! Golden from `tools/readfilter-conformance/ModelSegmentsDump.java`, whose runs are all in the
//! multi-sample mode: it genotypes and segments and stops, so the one file it writes is the
//! interval list this suite reads.
//!
//! The kernel segmenter is not ported, so the changepoints themselves are taken FROM the golden.
//! What is checked is everything the port does around them.
//!
//! # What this suite is for
//!
//!  * **a changepoint being the last index of its segment**, so the golden's own intervals are
//!    rebuilt from the indices they imply;
//!  * **the window sizes being a set**, so one named twice segments as one named once;
//!  * **a floor above every site's total leaving no het at all**, so the run falls back on the
//!    copy ratios;
//!  * **the heterozygosity test at the default threshold keeping a 39/21 site**, which is where
//!    the allele-fraction step is and therefore where the default run cuts;
//!  * **a lowered threshold changing nothing** on sites this unambiguous;
//!  * **the cap counting segments rather than changepoints, and keeping the last break**;
//!  * **and the two samples whose copy-ratio intervals disagree being refused.**

use gatk_corpus as corpus;
use gatk_tools::model_segments::{
    filter_by_total_count, genotype_hets, homozygous_log_ratio, is_heterozygous,
    maximum_changepoints_per_chromosome, segments_from_changepoints, window_sizes, AllelicCount,
    DEFAULT_GENOTYPING_BASE_ERROR_RATE, DEFAULT_GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD,
    MISMATCHED_COPY_RATIO_INTERVALS_MESSAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/model_segments.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The payload of the one line whose first two fields are the kind and the name.
fn field(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .unwrap_or_else(|| panic!("{kind} {name} is in the golden"));
    unescape(&line[prefix.len()..])
}

/// One run's interval list, as the intervals it holds.
fn segmentation(text: &str, label: &str) -> Vec<(i32, i32)> {
    let payload = field(text, "out", label);
    let (_, content) = payload
        .split_once('=')
        .expect("the file name and its content");
    content
        .lines()
        .filter(|line| !line.starts_with('@') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse::<i32>().expect("a start"),
                columns[2].parse::<i32>().expect("an end"),
            )
        })
        .collect()
}

/// The copy-ratio points of a named input, which are the segmenter's data points.
fn points(text: &str, name: &str) -> Vec<(i32, i32)> {
    let prefix = format!("counts\t{name}=");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .expect("the input is in the golden");
    unescape(&line[prefix.len()..])
        .lines()
        .filter(|line| line.starts_with("chr1\t"))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse::<i32>().expect("a start"),
                columns[2].parse::<i32>().expect("an end"),
            )
        })
        .collect()
}

/// The allelic counts of a named input.
fn allelic_counts(text: &str, name: &str) -> Vec<AllelicCount> {
    let prefix = format!("counts\t{name}=");
    let line = text
        .lines()
        .find(|line| line.starts_with(&prefix))
        .expect("the input is in the golden");
    unescape(&line[prefix.len()..])
        .lines()
        .filter(|line| line.starts_with("chr1\t"))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            AllelicCount {
                position: columns[1].parse().expect("a position"),
                reference_count: columns[2].parse().expect("a reference count"),
                alternate_count: columns[3].parse().expect("an alternate count"),
            }
        })
        .collect()
}

/// The INTERIOR changepoint indices a segmentation implies: the index of the point each segment
/// ends on, less the last, which the segmenter never returns.
///
/// The last index is the closing one the reference appends itself. Handing it in as well would
/// have the reference append it a second time and then read one point past the end.
fn changepoints(points: &[(i32, i32)], segments: &[(i32, i32)]) -> Vec<usize> {
    let mut indices: Vec<usize> = segments
        .iter()
        .map(|(_, end)| {
            points
                .iter()
                .position(|point| point.1 == *end)
                .unwrap_or_else(|| panic!("a point ends at {end}"))
        })
        .collect();
    assert_eq!(
        indices.pop(),
        Some(points.len() - 1),
        "the run closes on the last point"
    );
    indices
}

/// Every run's intervals are rebuilt from the changepoints they imply, which is the loop.
#[test]
fn the_segments_are_rebuilt_from_their_changepoints() {
    let text = golden();
    let points = points(&text, "ratios-a");
    assert_eq!(points.len(), 120);
    for label in [
        "two-samples",
        "copy-ratios-only",
        "penalty-high",
        "penalty-low",
        "capped-segments",
        "one-window",
        "one-window-twice",
        "het-threshold-low",
        "minimum-count-high",
    ] {
        let expected = segmentation(&text, label);
        assert!(!expected.is_empty(), "{label} has a segmentation");
        let indices = changepoints(&points, &expected);
        assert_eq!(
            segments_from_changepoints(&points, &indices),
            expected,
            "{label}"
        );
    }
}

/// The point count is not a changepoint index, so the closing guard always fires: the interior
/// changepoints alone close the last segment on the last point.
#[test]
fn the_closing_index_is_appended() {
    let text = golden();
    let points = points(&text, "ratios-a");
    let expected = segmentation(&text, "two-samples");
    let indices = changepoints(&points, &expected);
    assert!(!indices.contains(&(points.len() - 1)));
    let built = segments_from_changepoints(&points, &indices);
    assert_eq!(built.len(), indices.len() + 1);
    assert_eq!(
        built.last().map(|segment| segment.1),
        Some(points[points.len() - 1].1)
    );
}

/// Naming one window size twice is naming it once, and the two runs agree.
#[test]
fn the_window_sizes_are_a_set() {
    let text = golden();
    assert_eq!(window_sizes(&[16, 16]), window_sizes(&[16]));
    assert_eq!(window_sizes(&[16, 16]), vec![16]);
    assert_eq!(
        segmentation(&text, "one-window-twice"),
        segmentation(&text, "one-window")
    );
}

/// A floor above every site's total leaves no het at all, and the run then segments on the copy
/// ratios alone.
#[test]
fn a_floor_above_every_total_leaves_no_het() {
    let text = golden();
    let counts = allelic_counts(&text, "counts-a");
    assert!(counts.iter().all(|count| count.total_read_count() < 1000));
    assert!(filter_by_total_count(&counts, 1000).is_empty());
    assert_eq!(filter_by_total_count(&counts, 0).len(), counts.len());
    assert_eq!(
        segmentation(&text, "minimum-count-high"),
        segmentation(&text, "copy-ratios-only")
    );
}

/// Both shapes of site are called heterozygous at the default threshold, which is why the default
/// run cuts where the allele fractions step rather than where they stop being het.
#[test]
fn the_default_threshold_keeps_both_shapes_of_site() {
    let text = golden();
    let points = points(&text, "ratios-a");
    let counts = allelic_counts(&text, "counts-a");
    let hets = genotype_hets(
        &counts,
        &points,
        0,
        DEFAULT_GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD,
        DEFAULT_GENOTYPING_BASE_ERROR_RATE,
    )
    .expect("the beta answered");
    assert_eq!(hets.len(), counts.len());
    let balanced = counts
        .iter()
        .find(|count| count.reference_count == 30)
        .copied()
        .expect("a balanced site");
    let skewed = counts
        .iter()
        .find(|count| count.reference_count == 39)
        .copied()
        .expect("a skewed site");
    for count in [balanced, skewed] {
        let ratio =
            homozygous_log_ratio(count, DEFAULT_GENOTYPING_BASE_ERROR_RATE).expect("a ratio");
        assert!(
            ratio < DEFAULT_GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD,
            "{ratio}"
        );
        // And at -30 as well, which is why the lowered-threshold run is the default one.
        assert!(ratio < -30.0, "{ratio}");
        assert!(is_heterozygous(count, -30.0, DEFAULT_GENOTYPING_BASE_ERROR_RATE).expect("a call"));
    }
    assert_eq!(
        segmentation(&text, "het-threshold-low"),
        segmentation(&text, "two-samples")
    );
}

/// The cap counts segments, and the break it keeps is the LAST of the uncapped run's breaks.
#[test]
fn the_cap_keeps_the_last_break() {
    let text = golden();
    let points = points(&text, "ratios-a");
    assert_eq!(maximum_changepoints_per_chromosome(2), 1);
    let capped = segmentation(&text, "capped-segments");
    assert_eq!(capped.len(), 2);
    let uncapped = segmentation(&text, "penalty-low");
    assert_eq!(uncapped.len(), 5);
    let capped_breaks = changepoints(&points, &capped);
    let uncapped_breaks = changepoints(&points, &uncapped);
    assert_eq!(capped_breaks.len(), 1);
    assert_eq!(uncapped_breaks.len(), 4);
    assert_eq!(
        capped_breaks[0],
        uncapped_breaks[uncapped_breaks.len() - 1],
        "the last break survives"
    );
    assert!(!capped_breaks.contains(&uncapped_breaks[0]));
}

/// Two samples whose copy-ratio intervals disagree are refused by name.
#[test]
fn mismatched_copy_ratio_intervals_are_refused() {
    let text = golden();
    let message = field(&text, "error", "mismatched-intervals");
    assert_eq!(
        message,
        format!("java.lang.IllegalArgumentException:{MISMATCHED_COPY_RATIO_INTERVALS_MESSAGE}")
    );
}

/// The allele fractions decide at the default penalty and the copy ratios do without them: the
/// two runs cut at different places, which is what the fixture was built to show.
#[test]
fn the_allele_fractions_dominate_at_the_default_penalty() {
    let text = golden();
    let with_counts = segmentation(&text, "two-samples");
    let without = segmentation(&text, "copy-ratios-only");
    assert_eq!(
        with_counts,
        vec![(1, 30000), (30001, 80000), (80001, 120000)]
    );
    assert_eq!(without, vec![(1, 40000), (40001, 120000)]);
    assert_eq!(segmentation(&text, "penalty-high"), vec![(1, 120000)]);
}
