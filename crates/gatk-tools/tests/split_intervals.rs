//! Conformance for `SplitIntervals` against GATK 4.6.2.0, compared as the names of every file
//! written and the whole content of each.
//!
//! Golden from `tools/readfilter-conformance/SplitIntervalsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the file name's width comes from a logarithm of one less than the scatter count**, which
//!    at one is negative infinity, so `--interval-file-num-digits` is what saves the name;
//!  * **no intervals at all means the whole reference** filtered by `--min-contig-size`, and the
//!    filter is ignored the moment any `-L` is given;
//!  * **`--dont-mix-contigs` splits after the scatter**, so the files can outnumber the requested
//!    count, and the sublists come out in dictionary order;
//!  * **the merging rule is the common default**, which is ALL, so adjacent `-L` intervals are
//!    merged before anything is scattered;
//!  * **and the prefix and the extension are concatenated raw**.
//!
//! # What is compared
//!
//! Every byte of every file, sequence lines included, with the harness's mask on the `M5` and `UR`
//! fields. The intervals each run starts from are produced here by the ported interval argument
//! collection, so the merging rule is exercised rather than assumed.

use gatk_corpus as corpus;
use gatk_engine::interval::MergingRule;
use gatk_engine::interval_args::{load_intervals, SetRule};
use gatk_tools::preprocess_intervals::Sequence;
use gatk_tools::split_intervals::{split, Arguments, SplitError};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::interval::{Interval, IntervalList};

const CONTIG_LENGTH: i32 = 240;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/split_intervals.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn row(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(unescape)
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    header
}

fn sequences() -> Vec<(String, i32)> {
    vec![
        ("chr1".to_string(), CONTIG_LENGTH),
        ("chr2".to_string(), CONTIG_LENGTH),
    ]
}

/// The dictionary as the interval list carries it, with the harness's mask in place.
fn dictionary() -> Vec<Sequence> {
    sequences()
        .into_iter()
        .map(|(name, length)| Sequence {
            name,
            length,
            md5: Some("<masked>".to_string()),
            uri: Some("<masked>".to_string()),
        })
        .collect()
}

/// The `-L` arguments of a run, through the ported interval argument collection.
fn intervals(queries: &[&str], rule: MergingRule) -> Vec<Interval> {
    let header = header();
    let queries: Vec<String> = queries.iter().map(|q| (*q).to_string()).collect();
    let (loaded, _) = load_intervals(&queries, &header, SetRule::Union, rule, 0)
        .expect("intervals the dictionary allows");
    loaded
        .into_iter()
        .map(|interval| Interval::new(&interval.contig, interval.start, interval.end))
        .collect()
}

/// The interval list written for one shard, as the tool writes it.
///
/// The `@HD` line carries `SO:coordinate` for every mode but `INTERVAL_SUBDIVISION`, whose
/// `uniqued()` preprocessing clones the original header rather than stamping a sort order.
fn written(shard: &IntervalList, mode: gatk_engine::interval_list_scatter::ScatterMode) -> String {
    let intervals: Vec<gatk_tools::filter_intervals::Interval> = shard
        .intervals
        .iter()
        .map(|interval| gatk_tools::filter_intervals::Interval {
            contig: interval.contig.clone(),
            start: interval.start,
            end: interval.end,
        })
        .collect();
    gatk_tools::preprocess_intervals::write_list_with_sort_order(
        &dictionary(),
        &intervals,
        mode.stamps_sort_order().then_some("coordinate"),
    )
}

/// Every run of the dump: its intervals, if any, and its arguments.
fn run(label: &str) -> (Option<Vec<Interval>>, Arguments) {
    use gatk_engine::interval_list_scatter::ScatterMode::*;
    let base = Arguments::default();
    match label {
        "whole-genome-1" => (
            None,
            Arguments {
                scatter_count: 1,
                ..base
            },
        ),
        "whole-genome-3" => (
            None,
            Arguments {
                scatter_count: 3,
                ..base
            },
        ),
        "whole-genome-5" => (
            None,
            Arguments {
                scatter_count: 5,
                ..base
            },
        ),
        "min-contig-size" => (
            None,
            Arguments {
                scatter_count: 2,
                min_contig_size: 241,
                ..base
            },
        ),
        "min-contig-size-with-intervals" => (
            Some(intervals(&["chr1:1-100"], MergingRule::All)),
            Arguments {
                scatter_count: 2,
                min_contig_size: 241,
                ..base
            },
        ),
        "adjacent-merged" => (
            Some(intervals(&["chr1:1-50", "chr1:51-100"], MergingRule::All)),
            Arguments {
                scatter_count: 2,
                ..base
            },
        ),
        "adjacent-kept" => (
            Some(intervals(
                &["chr1:1-50", "chr1:51-100"],
                MergingRule::OverlappingOnly,
            )),
            Arguments {
                scatter_count: 2,
                ..base
            },
        ),
        "dont-mix-contigs" => (
            Some(intervals(
                &["chr1:1-100", "chr2:1-100", "chr2:201-240"],
                MergingRule::OverlappingOnly,
            )),
            Arguments {
                scatter_count: 2,
                dont_mix_contigs: true,
                ..base
            },
        ),
        "dont-mix-contigs-reversed" => (
            Some(intervals(
                &["chr2:1-100", "chr1:1-100"],
                MergingRule::OverlappingOnly,
            )),
            Arguments {
                scatter_count: 1,
                dont_mix_contigs: true,
                ..base
            },
        ),
        "named" => (
            Some(intervals(&["chr1:1-100"], MergingRule::All)),
            Arguments {
                scatter_count: 2,
                prefix: "shard-".to_string(),
                extension: ".list".to_string(),
                ..base
            },
        ),
        "one-digit" => (
            Some(intervals(&["chr1:1-100"], MergingRule::All)),
            Arguments {
                scatter_count: 2,
                num_digits: 1,
                ..base
            },
        ),
        "one-digit-twelve" => (
            Some(intervals(&["chr1:1-240"], MergingRule::All)),
            Arguments {
                scatter_count: 12,
                num_digits: 1,
                ..base
            },
        ),
        "one-digit-one" => (
            Some(intervals(&["chr1:1-100"], MergingRule::All)),
            Arguments {
                scatter_count: 1,
                num_digits: 1,
                ..base
            },
        ),
        "eight-digits" => (
            Some(intervals(&["chr1:1-100"], MergingRule::All)),
            Arguments {
                scatter_count: 2,
                num_digits: 8,
                ..base
            },
        ),
        "more-shards-than-intervals" => (
            Some(intervals(
                &["chr1:1-100", "chr1:201-240"],
                MergingRule::OverlappingOnly,
            )),
            Arguments {
                scatter_count: 10,
                subdivision_mode: BalancingWithoutIntervalSubdivision,
                ..base
            },
        ),
        label if label.starts_with("mode-") => (
            Some(intervals(
                &["chr1:1-100", "chr1:141-240", "chr2:1-50"],
                MergingRule::OverlappingOnly,
            )),
            Arguments {
                scatter_count: 3,
                subdivision_mode: match &label["mode-".len()..] {
                    "INTERVAL_SUBDIVISION" => IntervalSubdivision,
                    "BALANCING_WITHOUT_INTERVAL_SUBDIVISION" => BalancingWithoutIntervalSubdivision,
                    "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW" => {
                        BalancingWithoutIntervalSubdivisionWithOverflow
                    }
                    "INTERVAL_COUNT" => IntervalCount,
                    "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER" => {
                        IntervalCountWithDistributedRemainder
                    }
                    other => panic!("an unexpected mode: {other}"),
                },
                ..base
            },
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_run_writes_the_files_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "whole-genome-1",
        "whole-genome-3",
        "whole-genome-5",
        "min-contig-size",
        "min-contig-size-with-intervals",
        "mode-INTERVAL_SUBDIVISION",
        "mode-BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
        "mode-BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
        "mode-INTERVAL_COUNT",
        "mode-INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER",
        "adjacent-merged",
        "adjacent-kept",
        "dont-mix-contigs",
        "dont-mix-contigs-reversed",
        "named",
        "one-digit",
        "one-digit-twelve",
        "one-digit-one",
        "eight-digits",
        "more-shards-than-intervals",
    ] {
        let (given, arguments) = run(label);
        let shards =
            split(given.as_deref(), &sequences(), &arguments).expect("a run that finishes");
        let names: Vec<String> = shards.iter().map(|(name, _)| name.clone()).collect();
        assert_eq!(
            names.join(","),
            value(&text, "files", label),
            "{label}: the names"
        );
        for (name, shard) in &shards {
            assert_eq!(
                written(shard, arguments.subdivision_mode),
                value(&text, "list", &format!("{label},{name}")),
                "{label}/{name}"
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 20, "the golden's runs");
}

#[test]
fn the_two_refusals_carry_the_references_messages() {
    let text = golden();
    let given = intervals(&["chr1:1-100"], MergingRule::All);
    let zero = split(
        Some(&given),
        &sequences(),
        &Arguments {
            scatter_count: 0,
            ..Arguments::default()
        },
    )
    .expect_err("a scatter count of zero");
    assert_eq!(zero, SplitError::ScatterCountNotPositive);
    assert_eq!(
        format!("{}:{}", zero.java_class(), zero.message()),
        row(&text, "error", "zero-scatter-count").expect("the golden carries the refusal")
    );

    let digits = split(
        Some(&given),
        &sequences(),
        &Arguments {
            scatter_count: 2,
            num_digits: 0,
            ..Arguments::default()
        },
    )
    .expect_err("a zero number of digits");
    assert_eq!(digits, SplitError::DigitsOutOfRange(0));
    assert_eq!(
        format!("{}:{}", digits.java_class(), digits.message()),
        row(&text, "error", "zero-digits").expect("the golden carries the refusal")
    );
}
