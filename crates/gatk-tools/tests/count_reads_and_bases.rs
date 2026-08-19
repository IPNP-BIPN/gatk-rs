//! Conformance for `CountReads` and `CountBases` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/CountReadsAndBasesDump.java`.
//!
//! # What this suite is for
//!
//! Both tools are an `apply` that increments a counter, so every row here is about what the engine
//! does before `apply`:
//!
//!  * **the default filters are not nothing**: 8 of the fixture's 11 records reach `apply`, and 11
//!    only with every default filter disabled;
//!  * **`CountBases` counts the sequence rather than the span**, which the varied fixture makes
//!    visible: 15 + 5 + 10 bases where the 10 spans twenty reference bases, and a fourth read whose
//!    sequence is empty contributes nothing while still being a read;
//!  * **an added filter is additional**, so `NotDuplicateReadFilter` takes 8 to 7;
//!  * **and the two tools see the same reads**, which is the claim that makes them one archetype.
//!
//! The eleven-record fixture travels in the `read-walker` golden as base64 and is read from there
//! rather than duplicated. The varied fixture is written here, because the reference built it
//! inside the dump and no golden carries its bytes: what has to agree is the records, and the port
//! writes them with its own writer.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::{self, IntervalArguments};
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::not_duplicate;
use gatk_tools::count_reads::{count_bases, count_reads, default_read_filter, output};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/count_reads_and_bases.txt.gz"),
    )
}

/// The `read-walker` golden, which carries the shared fixture's bytes.
fn read_walker_golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_walker.txt.gz"),
    )
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    let prefix = format!("{kind}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries a {kind} row"))
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

/// A directory nothing else writes into.
///
/// The tests run in parallel and two of them unpack the same fixture, so a name built from the
/// process alone is a race: one test opens the BAI while the other is still writing it. The counter
/// is what makes each call its own directory.
fn directory(name: &str) -> std::path::PathBuf {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "gatk-rs-countreads-{}-{}-{name}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    dir
}

/// The eleven-record fixture, unpacked from the read-walker golden.
fn shared_fixture() -> ReadsDataSource {
    let text = read_walker_golden();
    let dir = directory("shared");
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    std::fs::write(&bam, corpus::decode_base64(field(&text, "bam"))).expect("the fixture");
    std::fs::write(&bai, corpus::decode_base64(field(&text, "bai"))).expect("the index");
    ReadsDataSource::open(&bam, &bai).expect("the fixture opens")
}

/// One record of the varied fixture.
fn varied_record(name: &str, start: i32, cigar: Cigar, bases: &[u8]) -> BamRecord {
    let mut record = BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        cigar,
        read_bases: bases.to_vec(),
        base_qualities: vec![30; bases.len()],
        mate_reference_index: -1,
        ..Default::default()
    };
    record.tags.insert(
        htsjdk_bam::tag::Tag::new(b"RG"),
        htsjdk_bam::tag::TagValue::Str("rg1".into()),
    );
    record
}

fn cigar(elements: &[(Op, u32)]) -> Cigar {
    Cigar {
        elements: elements
            .iter()
            .map(|(op, length)| CigarElement {
                op: *op,
                length: *length,
            })
            .collect(),
    }
}

/// The dump's second fixture: reads of fifteen, five, ten and zero bases.
fn varied_fixture() -> ReadsDataSource {
    let mut header = SamHeader::new();
    header.sequences.push(SequenceRecord::new("chr1", 200));
    header.set_sort_order("coordinate");
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "sample1");
    group.attributes.set("PL", "ILLUMINA");
    header.read_groups.push(group);

    let records = vec![
        varied_record("v001", 10, cigar(&[(Op::M, 15)]), b"ACGTACGTACGTACG"),
        varied_record("v002", 40, cigar(&[(Op::M, 5)]), b"ACGTA"),
        // Five bases either side of a ten-base deletion: a span of twenty, a length of ten.
        varied_record(
            "v003",
            60,
            cigar(&[(Op::M, 5), (Op::D, 10), (Op::M, 5)]),
            b"ACGTACGTAC",
        ),
        // No sequence at all, which is a length of zero rather than a refusal.
        varied_record("v004", 100, Cigar::default(), b""),
    ];

    let dir = directory("varied");
    let bam = dir.join("varied.bam");
    let bai = dir.join("varied.bai");
    let writer = BamWriter::new(Vec::new(), &header)
        .expect("a writer")
        .with_index();
    let mut writer = writer;
    for record in &records {
        writer.write(record).expect("a record");
    }
    let (bytes, index) = writer.finish_with_index().expect("a complete file");
    std::fs::write(&bam, bytes).expect("the fixture");
    std::fs::write(&bai, index).expect("the index");
    ReadsDataSource::open(&bam, &bai).expect("the fixture opens")
}

fn intervals(source: &ReadsDataSource, queries: &[&str]) -> Vec<SimpleInterval> {
    if queries.is_empty() {
        return Vec::new();
    }
    let arguments = IntervalArguments {
        include: queries.iter().map(|q| q.to_string()).collect(),
        ..Default::default()
    };
    interval_args::parse_intervals(&arguments, source.header())
        .expect("the intervals parse")
        .intervals
}

#[test]
fn every_count_matches_the_golden() {
    let text = golden();
    let source = shared_fixture();
    let header = source.header().clone();

    // label, intervals, and whether the default filters run at all.
    let cases: Vec<(&str, Vec<&str>, bool, bool)> = vec![
        ("all", vec![], true, false),
        // A reference changes what `apply` is handed and not which reads reach it.
        ("all-withref", vec![], true, false),
        ("chr1", vec!["chr1"], true, false),
        ("chr2", vec!["chr2"], true, false),
        ("deletion", vec!["chr1:140-160"], true, false),
        ("no-duplicates", vec![], true, true),
        ("no-filters", vec![], false, false),
    ];

    for (label, queries, defaults, no_duplicates) in cases {
        let bounds = intervals(&source, &queries);
        let filter = |read: &BamRecord| {
            if !defaults {
                return true;
            }
            default_read_filter(read, &header) && (!no_duplicates || not_duplicate(read))
        };

        let reads = count_reads(&source, &bounds, &filter).expect("the traversal runs");
        let bases = count_bases(&source, &bounds, &filter).expect("the traversal runs");
        assert_eq!(output(reads), row(&text, "reads", label), "{label}: reads");
        assert_eq!(output(bases), row(&text, "bases", label), "{label}: bases");
    }
}

#[test]
fn the_varied_fixture_separates_the_two_tools() {
    let text = golden();
    let source = varied_fixture();
    let header = source.header().clone();

    let filtered = |read: &BamRecord| default_read_filter(read, &header);
    let everything = |_: &BamRecord| true;

    assert_eq!(
        output(count_reads(&source, &[], &filtered).expect("the traversal runs")),
        row(&text, "reads", "varied"),
    );
    assert_eq!(
        output(count_bases(&source, &[], &filtered).expect("the traversal runs")),
        row(&text, "bases", "varied"),
    );
    assert_eq!(
        output(count_reads(&source, &[], &everything).expect("the traversal runs")),
        row(&text, "reads", "varied-nofilters"),
    );
    // The read with no sequence is the only difference between the two counts, and it contributes
    // nothing, so the base count is the same with the filters off.
    assert_eq!(
        output(count_bases(&source, &[], &everything).expect("the traversal runs")),
        row(&text, "bases", "varied-nofilters"),
    );
}

/// `-L chrZ` is refused before any read is counted.
///
/// The reference calls it a `MalformedGenomeLoc`, the same class it uses for an interval that runs
/// off the end of a contig it does know: "Contig chrZ given as location, but this contig isn't
/// present in the Fasta sequence dictionary".
#[test]
fn an_unknown_contig_is_refused_rather_than_counted() {
    let text = golden();
    assert_eq!(
        row(&text, "error", "reads-unknown-contig"),
        "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc",
    );
    assert_eq!(
        row(&text, "error", "bases-unknown-contig"),
        "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc",
    );

    let source = shared_fixture();
    let arguments = IntervalArguments {
        include: vec!["chrZ".to_string()],
        ..Default::default()
    };
    interval_args::parse_intervals(&arguments, source.header())
        .expect_err("chrZ is in no dictionary");
}

/// The output is the number and nothing else: no newline, which `print` rather than `println` is.
#[test]
fn the_output_has_no_trailing_newline() {
    assert_eq!(output(80), "80");
    assert!(!output(80).ends_with('\n'));
    let text = golden();
    // The dump escapes the file's bytes, so a newline would have travelled as `\n`.
    assert!(
        !text
            .lines()
            .any(|line| line.starts_with("bases\t") && line.contains("\\n")),
        "no count row carries an escaped newline"
    );
}
