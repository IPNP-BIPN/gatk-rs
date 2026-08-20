//! Conformance for the five interval-list scatter modes against Picard 3.4.0, as GATK 4.6.2.0
//! bundles it, compared as every ideal weight, every shard count and every shard.
//!
//! Golden from `tools/readfilter-conformance/IntervalListScatterDump.java`.
//!
//! # What this suite is for
//!
//!  * **only `INTERVAL_SUBDIVISION` uniques its input**, so the `overlapping` list comes out as
//!    two merged intervals under it and as four separate ones under the other four modes;
//!  * **the ideal weight is taken from the unique base count** even for the modes that never
//!    unique, and is floored at one;
//!  * **the no-subdivision modes raise it to the widest interval**, which is why `uneven` scattered
//!    ten ways comes back with far fewer than ten shards;
//!  * **the last shard takes everything left**, whatever its weight;
//!  * **the overflow mode reads the projection**, so its shards depend on how far into the scatter
//!    it is;
//!  * **the two interval-count modes put the remainder at opposite ends**;
//!  * **and the subdivision cut keeps the name and the strand** of the interval it cut.

use gatk_corpus as corpus;
use gatk_engine::interval_list_scatter::{scatter, ScatterError, ScatterMode};
use htsjdk_bam::interval::{Interval, IntervalList};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/interval_list_scatter.txt.gz"),
    )
}

fn named(contig: &str, start: i32, end: i32, negative: bool, name: &str) -> Interval {
    Interval::with_strand_and_name(contig, start, end, negative, Some(name))
}

/// The dump's seven interval lists, by name.
fn list(name: &str) -> IntervalList {
    let mut list = IntervalList::new(vec!["chr1".to_string(), "chr2".to_string()]);
    list.intervals = match name {
        "even" => (0..5)
            .map(|i| {
                named(
                    "chr1",
                    1 + i * 200,
                    100 + i * 200,
                    false,
                    &format!("even{i}"),
                )
            })
            .collect(),
        "uneven" => vec![
            named("chr1", 1, 1, false, "one"),
            named("chr1", 101, 110, false, "ten"),
            named("chr1", 201, 1200, false, "thousand"),
            named("chr1", 1301, 1350, false, "fifty"),
        ],
        "overlapping" => vec![
            named("chr1", 1, 100, false, "a"),
            named("chr1", 50, 150, false, "b"),
            named("chr1", 151, 200, false, "c"),
            named("chr1", 400, 500, false, "d"),
        ],
        "two-contigs" => vec![
            named("chr1", 1, 100, false, "first"),
            named("chr2", 1, 100, false, "second"),
            named("chr2", 201, 400, false, "third"),
        ],
        "unsorted" => vec![
            named("chr2", 1, 50, false, "z"),
            named("chr1", 301, 400, true, "y"),
            named("chr1", 1, 100, false, "x"),
        ],
        "single" => vec![named("chr1", 1, 1000, false, "only")],
        "one-base" => vec![named("chr1", 1, 1, false, "tiny")],
        other => panic!("{other} is in the golden but not configured here"),
    };
    list
}

fn mode(name: &str) -> ScatterMode {
    *ScatterMode::ALL
        .iter()
        .find(|mode| mode.name() == name)
        .unwrap_or_else(|| panic!("an unexpected mode: {name}"))
}

/// One shard, rendered the way the dump renders it.
fn render(shard: &IntervalList) -> String {
    shard
        .intervals
        .iter()
        .map(|interval| {
            format!(
                "{}:{}-{}:{}:{}",
                interval.contig,
                interval.start,
                interval.end,
                interval.strand(),
                interval.name.clone().unwrap_or_default()
            )
        })
        .collect::<Vec<String>>()
        .join(";")
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    let mut weights = 0;
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let (kind, rest) = line.split_once('\t').expect("a kind");
        if kind == "error" {
            let (label, message) = rest.split_once('\t').expect("a message");
            let count = if label == "zero-count" { 0 } else { -1 };
            let error = scatter(&list("even"), ScatterMode::IntervalSubdivision, count)
                .expect_err("a scatter count below one");
            assert_eq!(error, ScatterError::ScatterCountBelowOne);
            assert_eq!(
                format!("{}:{}", error.java_class(), error.message()),
                message
            );
            rows += 1;
            continue;
        }
        let (label, value) = rest.split_once('=').expect("a value");
        let parts: Vec<&str> = label.split(',').collect();
        let source = list(parts[0]);
        let scatterer = mode(parts[1]);
        let count: i32 = parts[2].parse().expect("a scatter count");
        let ours = match kind {
            "weight" => {
                weights += 1;
                let processed = scatterer.preprocess(&source);
                scatterer.ideal_split_weight(&processed, count).to_string()
            }
            "shards" => scatter(&source, scatterer, count)
                .expect("a scatter count of at least one")
                .len()
                .to_string(),
            "shard" => {
                let index: usize = parts[3].parse().expect("a shard index");
                let shards = scatter(&source, scatterer, count).expect("a scatter");
                render(&shards[index])
            }
            other => panic!("an unexpected row: {other}"),
        };
        assert_eq!(ours, value, "{kind} {label}");
        rows += 1;
    }
    assert_eq!(weights, 175, "the golden's combinations");
    assert_eq!(rows, 785, "the golden's row count");
}
