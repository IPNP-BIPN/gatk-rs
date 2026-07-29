//! The indexed query path, checked against a full scan of the same file.
//!
//! This is **not** conformance: the reference has not been asked anything here. The suite that
//! asks it is `readsdatasource` in `tools/conformance/manifest.json`, and it has no golden yet.
//!
//! What this does check is the half a golden would not isolate: that going through the `.bai`
//! returns exactly what filtering every record in the file returns. The filter is the same code
//! in both arms on purpose, so a wrong filter passes here and fails the oracle; what cannot pass
//! here is a wrong chunk list, a cursor that loses a record spanning two BGZF blocks, or a query
//! that reads a block it should have skipped. That is precisely the part where `noodles` parses
//! the index and this crate interprets it.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::{FilterState, IntervalFilter, ReadsDataSource};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;

const CONTIG_LENGTH: i32 = 100_000;

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2", "chr3"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    header
}

fn mapped(name: &str, reference_index: i32, start: i32, cigar: Vec<(Op, u32)>) -> BamRecord {
    let elements: Vec<CigarElement> = cigar
        .into_iter()
        .map(|(op, length)| CigarElement { op, length })
        .collect();
    let read_length: usize = elements
        .iter()
        .filter(|e| e.op.consumes_read_bases())
        .map(|e| e.length as usize)
        .sum();
    BamRecord {
        read_name: name.to_string(),
        flags: 0,
        reference_index,
        alignment_start: start,
        mapping_quality: 60,
        cigar: Cigar::new(elements),
        read_bases: b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"[..read_length].to_vec(),
        base_qualities: vec![30; read_length],
        ..BamRecord::default()
    }
}

fn unmapped_at_mate(name: &str, reference_index: i32, start: i32) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        // paired, first of pair, unmapped
        flags: 0x1 | 0x40 | 0x4,
        reference_index,
        alignment_start: start,
        mate_reference_index: reference_index,
        mate_alignment_start: start,
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![30; 10],
        ..BamRecord::default()
    }
}

fn unplaced(name: &str) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        flags: 0x4,
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![30; 10],
        ..BamRecord::default()
    }
}

/// The same shape as the oracle harness's fixture, built through htsjdk-rs's own writer and
/// indexer. Enough records to fill several BGZF blocks, so a record straddling a block boundary
/// is actually exercised rather than assumed.
fn fixture() -> Vec<BamRecord> {
    let mut records = vec![
        mapped("r001", 0, 100, vec![(Op::M, 10)]),
        mapped("r002", 0, 150, vec![(Op::M, 10)]),
        mapped("r003", 0, 195, vec![(Op::M, 10)]),
        mapped("r004", 0, 200, vec![(Op::M, 10)]),
        mapped("r005", 0, 250, vec![(Op::M, 10)]),
        mapped("r006", 0, 300, vec![(Op::M, 5), (Op::D, 10), (Op::M, 5)]),
        unmapped_at_mate("u001", 0, 300),
    ];
    // A filler run, coordinate-ordered, large enough to cross both a 16 kb linear-index window
    // and several BGZF blocks.
    for i in 0..2_000 {
        records.push(mapped(
            &format!("f{i:04}"),
            0,
            400 + i * 10,
            vec![(Op::M, 10)],
        ));
    }
    records.push(mapped("r009", 0, 99_995, vec![(Op::M, 6)]));
    records.push(mapped("r101", 1, 100, vec![(Op::M, 10)]));
    records.push(mapped("r102", 1, 5_000, vec![(Op::M, 10)]));
    // chr3 gets no reads: its linear index is empty.
    records.push(unplaced("x001"));
    records.push(unplaced("x002"));
    records.push(unplaced("x003"));
    records
}

fn write_fixture(dir: &std::path::Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let mut writer = BamWriter::new(Vec::new(), &header())
        .expect("header")
        .with_index();
    for record in fixture() {
        writer.write(&record).expect("write");
    }
    let (bam, bai) = writer.finish_with_index().expect("finish");

    let bam_path = dir.join("reads.bam");
    let bai_path = dir.join("reads.bai");
    std::fs::write(&bam_path, bam).unwrap();
    std::fs::write(&bai_path, bai).unwrap();
    (bam_path, bai_path)
}

/// The full-scan answer: every record in the file, through the same stateful filter.
fn by_full_scan(source: &ReadsDataSource, intervals: &[SimpleInterval]) -> Vec<String> {
    let converted: Vec<_> = intervals
        .iter()
        .map(|i| {
            gatk_engine::reads::convert_simple_interval_to_query_interval(i, source.header())
                .expect("the fixture declares every contig queried here")
        })
        .collect();
    let optimized = gatk_engine::reads::optimize_intervals(&converted);
    let mut filter = IntervalFilter::new(&optimized, false);
    let mut kept = Vec::new();
    for record in source.iter_all().expect("full scan") {
        match filter.compare_to_filter(&record) {
            FilterState::MatchesFilter => kept.push(record.read_name),
            FilterState::ContinueIteration => {}
            FilterState::StopIteration => break,
        }
    }
    kept
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval {
        contig: contig.to_string(),
        start,
        end,
    }
}

#[test]
fn the_indexed_query_returns_what_a_full_scan_returns() {
    let dir = std::env::temp_dir().join(format!("gatk-rs-readsquery-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (bam, bai) = write_fixture(&dir);
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");

    let queries: Vec<Vec<SimpleInterval>> = vec![
        vec![interval("chr1", 100, 100)],
        vec![interval("chr1", 105, 105)],
        vec![interval("chr1", 99, 99)],
        vec![interval("chr1", 150, 160)],
        vec![interval("chr1", 200, 200)],
        vec![interval("chr1", 300, 300)],
        vec![interval("chr1", 305, 305)],
        vec![interval("chr1", 1, 100_000)],
        vec![interval("chr1", 16_380, 16_395)],
        vec![interval("chr1", 16_384, 16_384)],
        vec![interval("chr1", 99_995, 100_000)],
        vec![interval("chr1", 100_001, 100_010)],
        vec![interval("chr2", 1, 1_000)],
        vec![interval("chr3", 1, 1_000)],
        // Abutting, so they merge into one.
        vec![interval("chr1", 100, 200), interval("chr1", 201, 300)],
        // The same interval twice.
        vec![interval("chr1", 100, 200), interval("chr1", 100, 200)],
        // Out of order, across contigs.
        vec![interval("chr2", 1, 1_000), interval("chr1", 100, 200)],
        vec![
            interval("chr1", 250, 260),
            interval("chr2", 1, 100_000),
            interval("chr1", 100, 200),
        ],
    ];

    let mut compared = 0;
    for intervals in &queries {
        let indexed: Vec<String> = source
            .query(intervals)
            .expect("query")
            .into_iter()
            .map(|r| r.read_name)
            .collect();
        let scanned = by_full_scan(&source, intervals);
        let label: Vec<String> = intervals
            .iter()
            .map(|i| format!("{}:{}-{}", i.contig, i.start, i.end))
            .collect();
        assert_eq!(indexed, scanned, "{}", label.join(" "));
        compared += 1;
    }

    // Named answers, so a change that breaks both arms the same way still fails.
    let names = |intervals: &[SimpleInterval]| -> Vec<String> {
        source
            .query(intervals)
            .unwrap()
            .into_iter()
            .map(|r| r.read_name)
            .collect()
    };
    assert_eq!(names(&[interval("chr1", 100, 100)]), ["r001"]);
    // r003 spans 195-204 and r004 starts at 200; r002 ends at 159.
    assert_eq!(names(&[interval("chr1", 200, 200)]), ["r003", "r004"]);
    // The unmapped read at its mate's coordinate is returned, and the deletion carries r006 here.
    assert_eq!(names(&[interval("chr1", 300, 300)]), ["r006", "u001"]);
    // Inside the deletion, so only r006: an unmapped read has no span to overlap with.
    assert_eq!(names(&[interval("chr1", 305, 305)]), ["r006"]);
    // Abutting intervals merge, so the read spanning the join comes back once.
    let abutting = names(&[interval("chr1", 100, 200), interval("chr1", 201, 300)]);
    assert_eq!(abutting.iter().filter(|n| *n == "r003").count(), 1);
    assert!(names(&[interval("chr3", 1, 1_000)]).is_empty());

    let unmapped: Vec<String> = source
        .query_unmapped()
        .expect("unmapped")
        .into_iter()
        .map(|r| r.read_name)
        .collect();
    assert_eq!(unmapped, ["x001", "x002", "x003"]);

    std::fs::remove_dir_all(&dir).ok();
    println!("{compared} interval queries, indexed and scanned answers identical");
}

#[test]
fn a_contig_the_reads_do_not_declare_is_an_error_and_not_an_empty_answer() {
    let dir = std::env::temp_dir().join(format!("gatk-rs-readsquery-x-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (bam, bai) = write_fixture(&dir);
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");

    let error = source.query(&[interval("chrX", 1, 100)]).unwrap_err();
    assert_eq!(
        error,
        gatk_engine::reads::ReadsError::ContigNotInDictionary("chrX".to_string())
    );

    std::fs::remove_dir_all(&dir).ok();
}
