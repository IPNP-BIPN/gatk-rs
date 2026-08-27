//! Conformance for `DepthOfCoverage` against GATK 4.6.2.0, compared as the per-locus table of every
//! run that wrote one, the set of files every run wrote, and every refusal.
//!
//! Golden from `tools/readfilter-conformance/DepthOfCoverageDump.java`.
//!
//! The summary and interval tables are reported by the golden but not recomputed here: their
//! quantiles come from the partitioned data store, which is not ported. What is ported is which
//! bases are counted, how, and which files a set of arguments produces.
//!
//! # What this suite is for
//!
//!  * **every base of the interval being a row**, covered or not;
//!  * **the partition being the sample**, not the read group;
//!  * **`--min-base-quality` filtering a base rather than a read**, and a ceiling doing the same;
//!  * **the base-count breakdown's `A: C: G: T: N:` order**;
//!  * **each omission removing its own files**;
//!  * **and a base quality outside the byte range failing two different ways.**

use gatk_corpus as corpus;
use gatk_tools::depth_of_coverage::{
    base_counts, check_base_quality, per_locus, per_locus_header, per_locus_row, written_suffixes,
    CoverageError, Omissions, Read, BASE_COUNT_ORDER,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/depth_of_coverage.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn line(text: &str, kind: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
        .to_string()
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

fn samples() -> Vec<String> {
    vec!["sampleA".to_string(), "sampleB".to_string()]
}

/// The five reads the dump wrote, with the bases it generated for them.
fn reads() -> Vec<Read> {
    let of = |name: &str, sample: &str, start: i32, length: usize, quality: i32| {
        let mut bases = String::new();
        while bases.len() < length {
            bases.push_str("ACGT");
        }
        Read {
            name: name.to_string(),
            sample: sample.to_string(),
            contig: "chr1".to_string(),
            start,
            bases: bases.as_bytes()[..length].to_vec(),
            base_qualities: vec![quality; length],
        }
    };
    vec![
        of("a1", "sampleA", 1000, 10, 30),
        of("a2", "sampleA", 1000, 10, 30),
        of("b1", "sampleB", 1000, 5, 30),
        of("lowq", "sampleB", 1000, 10, 5),
        of("highq", "sampleB", 1000, 10, 60),
    ]
}

/// The per-locus table one run wrote, rendered by the port.
fn rendered(minimum: i32, maximum: i32, print_base_counts: bool) -> String {
    let samples = samples();
    let mut out = per_locus_header(&samples, print_base_counts);
    out.push('\n');
    for locus in per_locus(&reads(), &samples, "chr1", 995, 1015, minimum, maximum) {
        out.push_str(&per_locus_row(&locus, print_base_counts));
        out.push('\n');
    }
    out
}

#[test]
fn every_per_locus_table_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, minimum, maximum, print_base_counts) in [
        ("default", 0, 127, false),
        ("base-counts", 0, 127, true),
        ("min-baseq", 10, 127, false),
        ("max-baseq", 0, 40, false),
        ("omit-locus-table", 0, 127, false),
        ("omit-per-sample", 0, 127, false),
        ("omit-intervals", 0, 127, false),
    ] {
        assert_eq!(
            rendered(minimum, maximum, print_base_counts),
            unescape(&line(&text, "out", &format!("{label}.base"))),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 7, "the runs that wrote a per-locus table");
}

#[test]
fn every_run_wrote_the_files_its_arguments_ask_for() {
    let text = golden();
    let mut compared = 0;
    for (label, omissions) in [
        ("default", Omissions::default()),
        (
            "omit-locus-table",
            Omissions {
                locus_table: true,
                ..Omissions::default()
            },
        ),
        (
            "omit-per-base",
            Omissions {
                depth_output_at_each_base: true,
                ..Omissions::default()
            },
        ),
        (
            "omit-per-sample",
            Omissions {
                per_sample_statistics: true,
                ..Omissions::default()
            },
        ),
        (
            "omit-intervals",
            Omissions {
                interval_statistics: true,
                ..Omissions::default()
            },
        ),
    ] {
        assert_eq!(
            written_suffixes(&omissions).join(","),
            line(&text, "files", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the omissions the golden carries");
}

/// Every base of the interval, including the five before the reads and the five after.
#[test]
fn every_base_of_the_interval_is_a_row() {
    let table = per_locus(&reads(), &samples(), "chr1", 995, 1015, 0, 127);
    assert_eq!(table.len(), 21, "the interval's length, not the reads'");
    assert_eq!(table[0].position, 995);
    assert_eq!(table[0].total_depth(), 0, "before any read");
    assert_eq!(table.last().expect("a row").position, 1015);
    assert_eq!(table.last().expect("a row").total_depth(), 0, "after them");
    // The reads themselves cover 1000 to 1009.
    assert!(table.iter().any(|locus| locus.total_depth() > 0));
    assert_eq!(
        table.iter().filter(|locus| locus.total_depth() > 0).count(),
        10
    );
}

/// The two read groups of sampleA are ONE column of depth two.
#[test]
fn the_partition_is_the_sample() {
    let table = per_locus(&reads(), &samples(), "chr1", 1000, 1000, 0, 127);
    assert_eq!(
        table[0].depths.len(),
        2,
        "two samples, not three read groups"
    );
    assert_eq!(table[0].depths[0], 2, "sampleA's two read groups");
    assert_eq!(table[0].depths[1], 3, "sampleB's three reads");
    assert_eq!(table[0].total_depth(), 5);
    // The average divides by the SAMPLES, so it is 2.5 rather than 5 over three groups.
    assert!((table[0].average_depth() - 2.5).abs() < 1e-12);
}

/// A base, not a read, and from both directions.
#[test]
fn the_quality_bounds_filter_a_base() {
    // The floor removes the quality-5 read's bases and nothing else.
    let raised = per_locus(&reads(), &samples(), "chr1", 1000, 1000, 10, 127);
    assert_eq!(raised[0].depths[0], 2, "sampleA is untouched");
    assert_eq!(raised[0].depths[1], 2, "sampleB loses one");
    // The ceiling removes the quality-60 read's, which is a filter from above.
    let lowered = per_locus(&reads(), &samples(), "chr1", 1000, 1000, 0, 40);
    assert_eq!(lowered[0].depths[1], 2);
    assert!(base_counts(30, 0, 127));
    assert!(!base_counts(60, 0, 40), "the ceiling is inclusive");
    assert!(base_counts(40, 0, 40));
    assert!(!base_counts(5, 10, 127));
    assert!(base_counts(10, 10, 127), "and so is the floor");
}

/// A fixed order, each pair followed by a space.
#[test]
fn the_base_counts_are_a_fixed_order() {
    assert_eq!(BASE_COUNT_ORDER, b"ACGTN");
    let table = per_locus(&reads(), &samples(), "chr1", 1000, 1000, 0, 127);
    let row = per_locus_row(&table[0], true);
    assert!(row.contains("A:2 C:0 G:0 T:0 N:0 "), "{row}");
    assert!(row.ends_with(' '), "every row ends with the trailing space");
    // The header names a base-counts column per sample only when it was asked for.
    assert!(per_locus_header(&samples(), true).contains("sampleA_base_counts"));
    assert!(!per_locus_header(&samples(), false).contains("base_counts"));
}

/// Two arguments that cannot be used together, and two ways for a quality to be out of range.
#[test]
fn the_three_refusals() {
    let text = golden();

    let (class, message) = refusal(&text, "both-interval-arguments");
    assert_eq!(
        class,
        "org.broadinstitute.barclay.argparser.CommandLineException"
    );
    let produced = CoverageError::MutuallyExclusive {
        argument: "omit-interval-statistics".to_string(),
        other: "calculate-coverage-over-genes".to_string(),
    };
    assert!(
        message.starts_with(&produced.message()),
        "{message} against {}",
        produced.message()
    );

    let (class, message) = refusal(&text, "quality-too-high");
    assert_eq!(
        class,
        "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
    );
    let produced = check_base_quality("min-base-quality", 200).expect_err("out of the byte range");
    assert_eq!(
        produced,
        CoverageError::BadByte {
            argument: "min-base-quality".to_string(),
            value: "200".to_string()
        }
    );
    assert!(message.starts_with(&produced.message()), "{message}");

    let (class, message) = refusal(&text, "quality-negative");
    assert_eq!(
        class,
        "org.broadinstitute.barclay.argparser.CommandLineException$OutOfRangeArgumentValue"
    );
    let produced = check_base_quality("min-base-quality", -1).expect_err("below the range");
    assert_eq!(
        produced,
        CoverageError::OutOfRange {
            argument: "min-base-quality".to_string(),
            value: "-1".to_string(),
            minimum: 0.0,
            maximum: 127.0
        }
    );
    assert!(message.starts_with(&produced.message()), "{message}");

    // The two failures are different exceptions, which is what makes -1 and 200 different cases.
    assert_ne!(
        refusal(&text, "quality-too-high").0,
        refusal(&text, "quality-negative").0
    );
    // A quality inside the range is accepted.
    assert!(check_base_quality("min-base-quality", 10).is_ok());
    assert!(check_base_quality("min-base-quality", 127).is_ok());
}
