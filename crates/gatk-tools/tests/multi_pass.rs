//! Conformance for the multi-pass walkers against the oracle.
//!
//! Golden from `tools/readfilter-conformance/MultiPassWalkerDump.java`, real
//! `MultiplePassVariantWalker`, `TwoPassVariantWalker` and `MultiplePassReadWalker` subclasses run
//! through the real command line.
//!
//! # What the golden settles
//!
//! ```text
//! filter  npass-1           0  1     one counter, one pass
//! filter  npass-3           0  3     the same counter, three passes: it accumulates
//! filter  read-three-pass   0  3     a new counter per pass...
//! filter  read-three-pass   2  3     ...so each reports only its own drops
//! filters-built  npass-3          1
//! filters-built  read-three-pass  3
//! ```
//!
//! and that `afterNthPass` runs after the **last** pass too, and that `TwoPassVariantWalker`
//! emits `afterFirstPass` between the passes and nothing at all after the second.
//!
//! # Where the fixtures come from
//!
//! The reads are `ReadWalkerDump.buildFixture`, which is the method `MultiPassWalkerDump` calls, so
//! the BAM and its index are taken from that suite's golden rather than rebuilt here: a fixture
//! rebuilt to match a description of a fixture is a second fixture. The variants are
//! `VariantWalkerDump.VCF`, whose six records are transcribed below because what this suite needs
//! from them is which one is filtered, and that is a property of the file rather than a coordinate
//! in the golden.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::counting::Counting;
use gatk_readfilter::with_header;
use gatk_tools::multi_pass::{
    two_pass_after_route, two_pass_apply_route, AfterPass, CountingVariantFilter,
    MultiplePassReadWalker, PassApply,
};
use htsjdk_bam::record::BamRecord;

/// One record of `VariantWalkerDump.VCF`: where it is, and whether a FILTER was applied to it.
struct Record {
    contig: &'static str,
    start: i64,
    /// `VariantContext.isFiltered()`: a non-empty FILTER column that is not `PASS` and not `.`.
    is_filtered: bool,
}

/// The six records of the fixture. Only `chr1:150` carries `LowQual`; `chr1:200` carries `.`,
/// which is *unfiltered* rather than filtered, and `PASSES_FILTERS` keeps it.
const RECORDS: [Record; 6] = [
    Record {
        contig: "chr1",
        start: 100,
        is_filtered: false,
    },
    Record {
        contig: "chr1",
        start: 150,
        is_filtered: true,
    },
    Record {
        contig: "chr1",
        start: 200,
        is_filtered: false,
    },
    Record {
        contig: "chr1",
        start: 201,
        is_filtered: false,
    },
    Record {
        contig: "chr1",
        start: 300,
        is_filtered: false,
    },
    Record {
        contig: "chr2",
        start: 100,
        is_filtered: false,
    },
];

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/multi_pass.txt.gz"),
    )
}

fn read_walker_golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_walker.txt.gz"),
    )
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

/// The rows a run of the port produces, in the dump's own shape.
struct Rows(Vec<String>);

impl Rows {
    fn new() -> Self {
        Rows(Vec::new())
    }

    fn events(&mut self, label: &str, events: &[String]) {
        for (index, event) in events.iter().enumerate() {
            self.0.push(format!("event\t{label}\t{index}\t{event}"));
        }
    }

    fn filters(&mut self, label: &str, counts: &[u64]) {
        for (index, count) in counts.iter().enumerate() {
            self.0.push(format!("filter\t{label}\t{index}\t{count}"));
        }
        self.0
            .push(format!("filters-built\t{label}\t{}", counts.len()));
        self.0.push(format!("summary\t{label}\tok"));
    }
}

/// `PASSES_FILTERS`, the predicate the probe walkers were given.
fn passes_filters(record: &&Record) -> bool {
    !record.is_filtered
}

/// A `MultiplePassVariantWalker` run, in the dump's event vocabulary.
///
/// `two_pass` routes through `TwoPassVariantWalker`'s two `final` methods, which is what makes the
/// missing `afterSecondPass` visible: the second pass ends with no event at all.
fn variant_run(passes: usize, two_pass: bool) -> (Vec<String>, Vec<u64>) {
    let records: Vec<&Record> = RECORDS.iter().collect();
    let mut filter = CountingVariantFilter::new(passes_filters);
    // Both callbacks append to one log, and the traversal holds both at once. A `RefCell` is the
    // shape the reference has anyway: `nthPassApply` and `afterNthPass` are two methods of one
    // walker writing to its own fields.
    let events: std::cell::RefCell<Vec<String>> = std::cell::RefCell::new(Vec::new());

    gatk_tools::multi_pass::traverse_multiple_pass(
        passes,
        &records,
        &mut filter,
        &mut |record, pass| {
            let where_ = format!("{}:{}", record.contig, record.start);
            if two_pass {
                match two_pass_apply_route(pass) {
                    PassApply::FirstPassApply => {
                        events.borrow_mut().push(format!("firstPassApply {where_}"))
                    }
                    PassApply::SecondPassApply => events
                        .borrow_mut()
                        .push(format!("secondPassApply {where_}")),
                    PassApply::Refused => panic!("the probe never runs a third pass"),
                }
            } else {
                events
                    .borrow_mut()
                    .push(format!("nthPassApply {pass} {where_}"));
            }
        },
        &mut |pass| {
            if two_pass {
                match two_pass_after_route(pass) {
                    AfterPass::AfterFirstPass => {
                        events.borrow_mut().push("afterFirstPass".to_string())
                    }
                    // The hole: nothing is emitted after the second pass.
                    AfterPass::Nothing => {}
                    AfterPass::Refused => panic!("the probe never runs a third pass"),
                }
            } else {
                events.borrow_mut().push(format!("afterNthPass {pass}"));
            }
        },
    );

    let mut events = events.into_inner();
    if two_pass {
        // `onTraversalSuccess`, which the probe overrides and which runs once after `traverse()`.
        events.push("onTraversalSuccess".to_string());
    }
    (events, vec![filter.filtered_count()])
}

#[test]
fn every_pass_is_the_one_the_reference_ran() {
    let text = golden();
    let expected: Vec<&str> = text
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();

    let mut rows = Rows::new();

    // The variant side. One counter for the run, so the counts are 1, 2 and 3 for one, two and
    // three passes over a file with one filtered record.
    for (label, passes, two_pass) in [
        ("two-pass", 2usize, true),
        ("npass-2", 2, false),
        ("npass-3", 3, false),
        ("npass-1", 1, false),
        ("npass-0", 0, false),
    ] {
        let (events, counts) = variant_run(passes, two_pass);
        rows.events(label, &events);
        rows.filters(label, &counts);
    }

    // The read side, over the fixture the ReadWalker suite carries.
    let dir = std::env::temp_dir().join(format!("gatk-rs-multipass-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let reads_golden = read_walker_golden();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    std::fs::write(&bam, corpus::decode_base64(field(&reads_golden, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&reads_golden, "bai"))).unwrap();
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    // A `MultiplePassReadWalker`'s default filter is `WellformedReadFilter` alone, and a new one
    // is built per pass after the first.
    let make_filter = {
        let header = header.clone();
        move || {
            let header = header.clone();
            Counting::leaf("WellformedReadFilter", move |read: &BamRecord| {
                with_header::wellformed(read, &header)
            })
        }
    };

    for (label, passes, interval_text) in [
        ("read-two-pass", 2usize, None),
        ("read-three-pass", 3, None),
        ("read-one-pass", 1, None),
        ("read-zero-pass", 0, None),
        ("read-two-pass-interval", 2, Some("chr1:100-200")),
    ] {
        let intervals: Vec<SimpleInterval> = interval_text
            .into_iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        let mut walker = MultiplePassReadWalker::new(make_filter.clone());
        let mut events: Vec<String> = Vec::new();
        for pass in 0..passes {
            walker
                .for_each_read(
                    &source,
                    &intervals,
                    false,
                    make_filter.clone(),
                    &mut |read: &BamRecord| {
                        events.push(format!("read {pass} {}", read.read_name));
                    },
                )
                .expect("the traversal runs");
            events.push(format!("endOfPass {pass}"));
        }
        let counts: Vec<u64> = walker
            .filters()
            .iter()
            .map(|filter| filter.filtered_count())
            .collect();
        rows.events(label, &events);
        rows.filters(label, &counts);
    }

    assert_eq!(rows.0.len(), expected.len(), "row count");
    for (produced, oracle) in rows.0.iter().zip(expected.iter()) {
        assert_eq!(produced, oracle);
    }
}

/// The routing table of `TwoPassVariantWalker`, including the case that is neither a callback nor
/// an error.
#[test]
fn the_second_pass_has_no_after_hook() {
    assert_eq!(two_pass_apply_route(0), PassApply::FirstPassApply);
    assert_eq!(two_pass_apply_route(1), PassApply::SecondPassApply);
    assert_eq!(two_pass_apply_route(2), PassApply::Refused);

    assert_eq!(two_pass_after_route(0), AfterPass::AfterFirstPass);
    assert_eq!(two_pass_after_route(1), AfterPass::Nothing);
    assert_eq!(two_pass_after_route(2), AfterPass::Refused);
}
