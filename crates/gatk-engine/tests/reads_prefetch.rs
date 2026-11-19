//! The traversal and the cursor are handed their blocks by a thread that inflates ahead of them.
//! What they return must be what inflating each block as it was needed returned: the same
//! records, the same refusal at the same record after the ones before it, and a chunk that still
//! ends at the first record starting at its end, with blocks past the chunk read by the cursor's
//! own reader. The file is large enough for many handfuls of blocks.

use gatk_engine::reads::{Chunk, ReadsDataSource, ReadsError};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;
use htsjdk_bgzf::BgzfWriter;
use std::io::Write;

const RECORDS: usize = 40_000;

fn scratch(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("gatk-rs-prefetch-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The header's bytes and each record's, inflated.
fn stream() -> (Vec<u8>, Vec<Vec<u8>>) {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", 100_000_000));
    let mut writer = BamWriter::new(Vec::new(), &header).expect("header");
    for i in 0..RECORDS {
        let length = 20 + i % 30;
        writer
            .write(&BamRecord {
                read_name: format!("r{i:06}{}", "x".repeat(i % 7)),
                reference_index: 0,
                alignment_start: 1 + i as i32,
                mapping_quality: 60,
                cigar: Cigar::new(vec![CigarElement {
                    op: Op::M,
                    length: length as u32,
                }]),
                read_bases: b"ACGT".repeat(15)[..length].to_vec(),
                base_qualities: vec![30; length],
                ..BamRecord::default()
            })
            .expect("write");
    }
    let inflated = htsjdk_bgzf::read::decompress_all(&writer.finish().expect("finish")).unwrap();
    let start = BamReader::new(&inflated).unwrap().records_start();
    let mut records = Vec::new();
    let mut at = start;
    while at < inflated.len() {
        let length = i32::from_le_bytes(inflated[at..at + 4].try_into().unwrap()) as usize;
        records.push(inflated[at..at + 4 + length].to_vec());
        at += 4 + length;
    }
    (inflated[..start].to_vec(), records)
}

/// The BAM of `header` then `records`, and the virtual pointer each record starts at.
fn write(
    dir: &std::path::Path,
    header: &[u8],
    records: &[Vec<u8>],
) -> (std::path::PathBuf, Vec<u64>) {
    let mut writer = BgzfWriter::new(Vec::new());
    writer.write_all(header).unwrap();
    let mut pointers = Vec::new();
    for record in records {
        pointers.push(writer.file_pointer().unwrap());
        writer.write_all(record).unwrap();
    }
    let path = dir.join("reads.bam");
    std::fs::write(&path, writer.into_inner().unwrap()).unwrap();
    (path, pointers)
}

fn names(
    records: impl Iterator<Item = Result<BamRecord, ReadsError>>,
) -> Vec<Result<String, ReadsError>> {
    records.map(|r| r.map(|r| r.read_name)).collect()
}

/// `BamReader` over the whole inflated stream, which is what the traversal answers to.
fn whole_stream(header: &[u8], records: &[Vec<u8>]) -> Vec<Result<String, ReadsError>> {
    let inflated: Vec<u8> = header
        .iter()
        .chain(records.iter().flatten())
        .copied()
        .collect();
    BamReader::new(&inflated)
        .unwrap()
        .map(|r| {
            r.map(|r| r.read_name)
                .map_err(|e| ReadsError::Malformed(format!("{e:?}")))
        })
        .collect()
}

fn expected_names(records: &[Vec<u8>]) -> Vec<Result<String, ReadsError>> {
    records
        .iter()
        .map(|r| Ok(BamRecord::decode(r).unwrap().unwrap().0.read_name))
        .collect()
}

#[test]
fn a_record_that_does_not_decode_is_refused_after_the_ones_before_it() {
    let dir = scratch("refused");
    let (header, mut records) = stream();
    let damaged = 35_000;
    // A name longer than the record: the record's length still holds, its decoding fails.
    records[damaged][4 + 8] = 255;
    let refusal = format!("{:?}", BamRecord::decode(&records[damaged]).unwrap_err());
    let (path, pointers) = write(&dir, &header, &records);
    let source = ReadsDataSource::open_unindexed(&path).unwrap();

    let mut wanted = expected_names(&records[..damaged]);
    let traversal = names(source.records_all().unwrap());
    assert_eq!(traversal, whole_stream(&header, &records));
    wanted.push(Err(ReadsError::Malformed(refusal)));
    let chunk = Chunk {
        start: pointers[0],
        end: pointers[damaged + 10],
    };
    assert_eq!(names(source.read_chunk_iter(&chunk)), wanted);
}

#[test]
fn a_chunk_ends_at_the_first_record_starting_at_its_end() {
    let dir = scratch("chunk");
    let (header, records) = stream();
    let (path, pointers) = write(&dir, &header, &records);
    let source = ReadsDataSource::open_unindexed(&path).unwrap();
    for (from, to) in [
        (0, 2_500),
        (100, 33_333),
        (0, RECORDS - 1),
        (39_000, 39_001),
    ] {
        let chunk = Chunk {
            start: pointers[from],
            end: pointers[to],
        };
        assert_eq!(
            names(source.read_chunk_iter(&chunk)),
            expected_names(&records[from..to]),
            "{from}..{to}"
        );
    }
    assert_eq!(
        names(source.records_all().unwrap()),
        expected_names(&records)
    );
}

#[test]
fn a_stream_cut_inside_its_last_record_is_truncated_or_eof_as_before() {
    let dir = scratch("cut");
    let (header, mut records) = stream();
    let last = records.last_mut().unwrap();
    let length = last.len();
    last.truncate(length - 20);
    let (path, pointers) = write(&dir, &header, &records);
    let source = ReadsDataSource::open_unindexed(&path).unwrap();

    assert_eq!(
        names(source.records_all().unwrap()),
        whole_stream(&header, &records)
    );

    let mut wanted = expected_names(&records[..RECORDS - 1]);
    wanted.push(Err(ReadsError::PrematureEof {
        expected: (length - 4 - 32) as i64,
        received: length - 20 - 4 - 32,
        file: path.display().to_string(),
    }));
    let chunk = Chunk {
        start: pointers[0],
        end: u64::MAX,
    };
    assert_eq!(names(source.read_chunk_iter(&chunk)), wanted);
}
