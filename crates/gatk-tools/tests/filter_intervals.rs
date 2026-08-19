//! Conformance for `FilterIntervals` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FilterIntervalsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the solitary-interval rule**, which removes a contig's only survivor and can therefore turn
//!    a successful run into a refusal;
//!  * **the inclusive annotation bounds**;
//!  * **the strictly-greater count rules**, where one sample of two is not enough;
//!  * **and the interval-list output**, which is a format nothing else here writes.

use gatk_corpus as corpus;
use gatk_tools::filter_intervals::{
    self, FilterError, Interval, DEFAULT_EXTREME_MAXIMUM_PERCENTILE,
    DEFAULT_EXTREME_MINIMUM_PERCENTILE, DEFAULT_EXTREME_PERCENTAGE, DEFAULT_LOW_COUNT_PERCENTAGE,
    DEFAULT_LOW_COUNT_THRESHOLD, DEFAULT_MAXIMUM_GC_CONTENT, DEFAULT_MINIMUM_GC_CONTENT,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_intervals.txt.gz"),
    )
}

fn file(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .replace("\\t", "\t")
        .replace("\\n", "\n")
}

/// The dump's five intervals.
fn intervals() -> Vec<Interval> {
    (0..5)
        .map(|index| {
            let start = 1 + 100 * index;
            Interval {
                contig: "chr1".to_string(),
                start,
                end: start + 99,
            }
        })
        .collect()
}

fn contigs() -> Vec<String> {
    vec!["chr1".to_string(); 5]
}

const GC: [f64; 5] = [0.05, 0.1, 0.5, 0.9, 0.95];
const COUNTS_ONE: [f64; 5] = [5.0, 50.0, 100.0, 150.0, 5000.0];
const COUNTS_TWO: [f64; 5] = [5.0, 60.0, 110.0, 160.0, 6000.0];

const SEQUENCES: [(&str, i32); 1] = [("chr1", 1000)];

fn sequences() -> Vec<(String, i32)> {
    SEQUENCES
        .iter()
        .map(|(name, length)| (name.to_string(), *length))
        .collect()
}

/// The surviving intervals, written as the tool writes them.
fn survivors(mask: &[bool]) -> String {
    let kept: Vec<Interval> = intervals()
        .into_iter()
        .enumerate()
        .filter(|(index, _)| !mask[*index])
        .map(|(_, interval)| interval)
        .collect();
    filter_intervals::write(&sequences(), &kept)
}

/// The annotation path: the GC filter, then the solitary rule.
fn annotation_run(minimum: f64, maximum: f64) -> Result<String, FilterError> {
    let mut mask = vec![false; 5];
    filter_intervals::update_mask_by_annotation(&mut mask, &GC, minimum, maximum)?;
    filter_intervals::update_mask_by_solitary_intervals(&mut mask, &contigs())?;
    Ok(survivors(&mask))
}

/// The count path: low counts, then the percentile band, then the solitary rule.
fn count_run(
    counts: &[Vec<f64>],
    threshold: i32,
    minimum: f64,
    maximum: f64,
) -> Result<String, FilterError> {
    let mut mask = vec![false; 5];
    filter_intervals::update_mask_by_low_counts(
        &mut mask,
        counts,
        threshold,
        DEFAULT_LOW_COUNT_PERCENTAGE,
    )?;
    filter_intervals::update_mask_by_extreme_counts(
        &mut mask,
        counts,
        minimum,
        maximum,
        DEFAULT_EXTREME_PERCENTAGE,
    )?;
    filter_intervals::update_mask_by_solitary_intervals(&mut mask, &contigs())?;
    Ok(survivors(&mask))
}

#[test]
fn the_annotation_runs_match_the_golden() {
    let text = golden();
    assert_eq!(
        annotation_run(DEFAULT_MINIMUM_GC_CONTENT, DEFAULT_MAXIMUM_GC_CONTENT).expect("three pass"),
        file(&text, "list", "annotations")
    );
    assert_eq!(
        annotation_run(0.2, 0.92).expect("two pass"),
        file(&text, "list", "tight-gc")
    );
}

/// Bounds leaving exactly one interval leave none: the solitary rule takes it, and the count check
/// then refuses the run.
#[test]
fn a_solitary_survivor_is_removed_and_the_run_refused() {
    let text = golden();
    let error = annotation_run(0.2, 0.8).expect_err("one survivor is none");
    assert_eq!(error, FilterError::EverythingFiltered);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        file(&text, "error", "solitary").trim_end()
    );

    // The GC filter on its own leaves exactly one, which is what makes the rule visible.
    let mut mask = vec![false; 5];
    let passing =
        filter_intervals::update_mask_by_annotation(&mut mask, &GC, 0.2, 0.8).expect("one passes");
    assert_eq!(passing, 1);
}

#[test]
fn the_count_runs_match_the_golden() {
    let text = golden();
    let one = vec![COUNTS_ONE.to_vec()];
    let two = vec![COUNTS_ONE.to_vec(), COUNTS_TWO.to_vec()];

    assert_eq!(
        count_run(
            &one,
            DEFAULT_LOW_COUNT_THRESHOLD,
            DEFAULT_EXTREME_MINIMUM_PERCENTILE,
            DEFAULT_EXTREME_MAXIMUM_PERCENTILE
        )
        .expect("four pass"),
        file(&text, "list", "counts-one-sample")
    );
    assert_eq!(
        count_run(
            &two,
            DEFAULT_LOW_COUNT_THRESHOLD,
            DEFAULT_EXTREME_MINIMUM_PERCENTILE,
            DEFAULT_EXTREME_MAXIMUM_PERCENTILE
        )
        .expect("four pass"),
        file(&text, "list", "counts-two-samples")
    );
    assert_eq!(
        count_run(
            &one,
            120,
            DEFAULT_EXTREME_MINIMUM_PERCENTILE,
            DEFAULT_EXTREME_MAXIMUM_PERCENTILE
        )
        .expect("two pass"),
        file(&text, "list", "high-threshold")
    );
    assert_eq!(
        count_run(&one, DEFAULT_LOW_COUNT_THRESHOLD, 0.0, 100.0).expect("four pass"),
        file(&text, "list", "wide-percentiles")
    );
}

/// The low-count rule is strictly greater, so one sample of two below the threshold is not enough.
#[test]
fn one_sample_of_two_does_not_fail_an_interval() {
    let counts = vec![vec![5.0, 100.0], vec![100.0, 100.0]];
    let mut mask = vec![false; 2];
    filter_intervals::update_mask_by_low_counts(
        &mut mask,
        &counts,
        DEFAULT_LOW_COUNT_THRESHOLD,
        DEFAULT_LOW_COUNT_PERCENTAGE,
    )
    .expect("both pass");
    assert_eq!(mask, vec![false, false], "one of two is not more than half");

    // Both samples below it is.
    let counts = vec![vec![5.0, 100.0], vec![5.0, 100.0]];
    let mut mask = vec![false; 2];
    filter_intervals::update_mask_by_low_counts(
        &mut mask,
        &counts,
        DEFAULT_LOW_COUNT_THRESHOLD,
        DEFAULT_LOW_COUNT_PERCENTAGE,
    )
    .expect("one passes");
    assert_eq!(mask, vec![true, false]);
}

/// The two refusals that happen before any filtering.
#[test]
fn the_input_refusals_match_the_golden() {
    let text = golden();
    for (label, error) in [
        ("neither", FilterError::NoInputs),
        ("whole-contig", FilterError::EmptyIntersection),
    ] {
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            file(&text, "error", label).trim_end(),
            "{label}"
        );
    }
}

/// The output is a Picard interval list, whose last two columns carry no information.
#[test]
fn the_output_is_an_interval_list() {
    let text = golden();
    let written = file(&text, "list", "annotations");
    assert!(written.starts_with("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"));
    for line in written.lines().filter(|line| !line.starts_with('@')) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 5);
        assert_eq!(fields[3], "+");
        assert_eq!(fields[4], ".");
    }
}
