//! Conformance for the `VariantWalker` traversal against the oracle.
//!
//! Goldens from `tools/readfilter-conformance/VariantWalkerDump.java`, a real `VariantWalker` run
//! through the real command line.
//!
//! # What the golden settled, against the source
//!
//! `FeatureIntervalIterator.featureIsNovel` remembers **one** interval, so reading it suggests that
//! `-L chr1:100-160 -L chr1:500-600 -L chr1:100-160` hands the first records to `apply` twice. It
//! does not:
//!
//! ```text
//! count  gap-then-repeat  2      not 4
//! count  unsorted         4      and in sorted order, not the order given
//! ```
//!
//! The argument layer sorts and merges `-L` before the traversal ever sees it, so the
//! precondition the iterator documents is always satisfied through the command line and the
//! one-interval memory is unreachable from there. That is worth stating precisely: the behaviour is
//! real in the class and unreachable in the tool, which is a different claim from "it does not
//! happen".
//!
//! The intervals below are therefore the **merged** ones, as `IntervalArgumentCollection` produces
//! them, and [`gatk_engine::interval_args`] is what produces them in the port.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::variant_source::{traverse, Located};

/// One record of the fixture: everything the traversal and the dump look at.
struct Variant {
    contig: &'static str,
    start: i32,
    stop: i32,
    id: &'static str,
    alleles: &'static str,
    filters: &'static str,
}

impl Located for Variant {
    fn contig(&self) -> &str {
        self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn stop(&self) -> i32 {
        self.stop
    }
}

impl Variant {
    fn render(&self) -> String {
        format!(
            "{}:{}-{}|{}|{}|{}",
            self.contig, self.start, self.stop, self.id, self.alleles, self.filters
        )
    }
}

/// The dump's VCF, decoded. The stops are what the record decoder produces: `start + ref - 1`,
/// except for the record carrying `END`, whose stop is the END value.
const FIXTURE: &[Variant] = &[
    Variant {
        contig: "chr1",
        start: 100,
        stop: 100,
        id: ".",
        alleles: "A,T",
        filters: "PASS",
    },
    Variant {
        contig: "chr1",
        start: 150,
        stop: 153,
        id: "rs1",
        alleles: "ACGT,A",
        filters: "LowQual",
    },
    Variant {
        contig: "chr1",
        start: 200,
        stop: 200,
        id: ".",
        alleles: "G,C",
        filters: "unfiltered",
    },
    Variant {
        contig: "chr1",
        start: 201,
        stop: 201,
        id: ".",
        alleles: "T,G",
        filters: "PASS",
    },
    Variant {
        contig: "chr1",
        start: 300,
        stop: 400,
        id: ".",
        alleles: "A,<DEL>",
        filters: "PASS",
    },
    Variant {
        contig: "chr2",
        start: 100,
        stop: 100,
        id: ".",
        alleles: "C,G",
        filters: "PASS",
    },
];

/// Label, and the intervals **after** the argument layer has sorted, merged and subtracted. `None`
/// is an unrestricted traversal.
fn intervals(label: &str) -> Option<Vec<SimpleInterval>> {
    let interval = |contig: &str, start: i32, end: i32| {
        SimpleInterval::new(contig, start, end).expect("a valid interval")
    };
    match label {
        "all" => None,
        "one-interval" => Some(vec![interval("chr1", 100, 160)]),
        "boundary" => Some(vec![interval("chr1", 150, 200)]),
        "two-sorted" => Some(vec![interval("chr1", 100, 160), interval("chr1", 200, 210)]),
        // The two overlapping intervals are merged into one before the traversal.
        "two-overlapping" => Some(vec![interval("chr1", 100, 250)]),
        // The repeated interval collapses onto the first, which is why nothing is emitted twice.
        "gap-then-repeat" => Some(vec![interval("chr1", 100, 160), interval("chr1", 500, 600)]),
        // Given out of order, traversed in order.
        "unsorted" => Some(vec![interval("chr1", 100, 160), interval("chr1", 200, 210)]),
        "end-tail" => Some(vec![interval("chr1", 390, 410)]),
        "chr2" => Some(vec![interval("chr2", 1, 1000)]),
        "empty-contig-interval" => Some(vec![interval("chr2", 500, 600)]),
        // `-XL chr1` leaves the whole of chr2.
        "exclude" => Some(vec![interval("chr2", 1, 1000)]),
        "select-passing" => Some(vec![interval("chr1", 1, 1000)]),
        other => panic!("unknown label {other}"),
    }
}

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/variant_walker.txt.gz"),
    )
}

#[test]
fn every_traversal_matches_the_reference() {
    let text = golden();

    let mut labels: Vec<&str> = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("count\t") {
            let (label, _) = rest.split_once('\t').expect("a label and a count");
            labels.push(label);
        }
    }
    assert!(!labels.is_empty(), "the golden carries no traversals");

    for label in &labels {
        let restriction = intervals(label);
        let emitted = traverse(FIXTURE, restriction.as_deref());

        let expected_count: usize = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("count\t{label}\t")))
            .expect("a count row")
            .parse()
            .expect("a number");
        assert_eq!(
            emitted.len(),
            expected_count,
            "{label}: number of apply calls"
        );

        for (index, variant) in emitted.iter().enumerate() {
            let expected = text
                .lines()
                .find_map(|line| line.strip_prefix(&format!("apply\t{label}\t{index}\t")))
                .unwrap_or_else(|| panic!("{label}: no apply row {index}"));
            assert_eq!(variant.render(), expected, "{label}, apply {index}");
        }
    }

    println!("{} traversals identical", labels.len());
}

/// The two rows that correct what the source suggests: the argument layer sorts and merges, so the
/// iterator's one-interval memory is never exercised through a command line.
#[test]
fn the_argument_layer_makes_the_iterators_memory_unreachable() {
    let text = golden();
    let count = |label: &str| -> usize {
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("count\t{label}\t")))
            .expect("a count row")
            .parse()
            .expect("a number")
    };

    // Three intervals, the first and third identical, with a gap between them. A traversal that
    // took them as given would emit the first two records twice.
    assert_eq!(count("gap-then-repeat"), 2);
    // Given out of order, traversed in sorted order: the first apply is the earliest record.
    let first = text
        .lines()
        .find_map(|line| line.strip_prefix("apply\tunsorted\t0\t"))
        .expect("an apply row");
    assert!(first.starts_with("chr1:100-100"), "{first}");
    // Two overlapping intervals become one, so the records inside the overlap arrive once.
    assert_eq!(count("two-overlapping"), 4);
}
