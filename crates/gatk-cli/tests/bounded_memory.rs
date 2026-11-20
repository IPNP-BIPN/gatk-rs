//! A traversal's memory is bounded by a block and a read, not by the file.
//!
//! The conformance corpus is a few kilobytes, so nothing there notices a port that reads a BAM
//! whole: it does, and gets the same bytes. This guard runs `CountReads` over a BAM of several
//! gigabytes under an address-space limit several times smaller, which a source holding the file,
//! or a walker collecting its reads, cannot meet.
//!
//! The file is one header and one compressed unit of records, the unit repeated: BGZF blocks stand
//! alone, so the inflated stream is the unit's records over and over, and a file of gigabytes is
//! written at the speed of the disk rather than of a deflater.
//!
//! `SortSamSpark` sorts such a BAM by name, which every repeat of the unit reorders: its sorted runs
//! have to spill to `--tmp-dir` and merge back for the tool to fit.
//!
//! `ulimit -v` is `RLIMIT_AS`, which only Linux enforces, so the guard runs there. The tests take
//! turns, so that the disk holds one large file and its runs at a time.

#![cfg(target_os = "linux")]

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};
use std::io::Write;

/// The file's size: several times the limit below.
const FILE_BYTES: u64 = 3 << 30;
/// The sorted file's: twice the limit, and as much again in runs on the disk.
const SORTED_FILE_BYTES: u64 = 2 << 30;
/// The address space the tool runs in, in KiB, which is `ulimit -v`'s unit.
const LIMIT_KIB: u64 = 1 << 20;

/// A record that `WellformedReadFilter` keeps, with bases and qualities that hardly compress, so
/// the stream inflates to little more than the file.
fn record(index: u64, state: &mut u64) -> BamRecord {
    let mut next = || {
        // xorshift: a stream of bytes that deflate cannot shrink much, with no dependency.
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    };
    let mut tags = Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    BamRecord {
        tags,
        read_name: format!("r{index}"),
        reference_index: 0,
        alignment_start: 1000,
        mapping_quality: 60,
        cigar: Cigar::new(vec![CigarElement {
            op: Op::M,
            length: 150,
        }]),
        read_bases: (0..150).map(|_| b"ACGT"[(next() % 4) as usize]).collect(),
        base_qualities: (0..150).map(|_| (next() % 41) as u8).collect(),
        ..BamRecord::default()
    }
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", 1_000_000));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "sample1");
    header.read_groups.push(group);
    header
}

/// The BAM, of about `bytes`, and how many reads it holds.
fn write_bam(path: &std::path::Path, bytes: u64) -> u64 {
    let header = header();
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut payload = Vec::new();
    let reads_per_unit = 4096u64;
    for index in 0..reads_per_unit {
        record(index, &mut state)
            .encode_into(&mut payload)
            .expect("the record");
    }
    let mut unit = htsjdk_bgzf::BgzfWriter::new(Vec::new());
    unit.write_all(&payload).expect("the unit");
    let unit = unit.into_inner_without_terminator().expect("the unit");

    let mut file = std::io::BufWriter::new(std::fs::File::create(path).expect("the file"));
    file.write_all(&htsjdk_bam::writer::write_bam_header_block(&header).expect("the header"))
        .expect("the header");
    let repeats = bytes / unit.len() as u64 + 1;
    for _ in 0..repeats {
        file.write_all(&unit).expect("a unit");
    }
    file.write_all(&htsjdk_bgzf::EMPTY_GZIP_BLOCK)
        .expect("the terminator");
    file.flush().expect("the file");
    repeats * reads_per_unit
}

/// One test at a time.
static TURN: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// `gatk` with `arguments` (`$1`, `$2`, ... in them name `paths`) under the address-space limit.
fn run_limited(arguments: &str, paths: &[&std::path::Path]) -> std::process::Output {
    std::process::Command::new("sh")
        .arg("-c")
        .arg(format!("ulimit -v {LIMIT_KIB} && exec \"$0\" {arguments}"))
        .arg(env!("CARGO_BIN_EXE_gatk"))
        .args(paths)
        // Each glibc arena reserves 64 MiB of address space and a thread pool makes one per
        // thread, which is address space and not memory.
        .env("MALLOC_ARENA_MAX", "2")
        .env("RAYON_NUM_THREADS", "4")
        .output()
        .expect("the tool runs")
}

fn assert_success(output: &std::process::Output) {
    assert!(
        output.status.success(),
        "status {:?}\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn count_reads_walks_a_bam_larger_than_its_address_space() {
    let _turn = TURN.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let dir = std::env::temp_dir().join("gatk-cli-bounded-memory");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    let bam = dir.join("large.bam");
    let reads = write_bam(&bam, FILE_BYTES);
    assert!(std::fs::metadata(&bam).expect("the file").len() > 3 * (LIMIT_KIB << 10));

    let output = run_limited("CountReads --input \"$1\"", &[&bam]);
    let _ = std::fs::remove_dir_all(&dir);
    assert_success(&output);
    // `handleResult` prints "Tool returned:" and the count on the line after it.
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(stdout.lines().last(), Some(reads.to_string().as_str()));
}

#[test]
fn sort_sam_spark_spills_a_bam_larger_than_its_address_space() {
    let _turn = TURN.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let dir = std::env::temp_dir().join("gatk-cli-bounded-memory-sort");
    let _ = std::fs::remove_dir_all(&dir);
    let runs = dir.join("runs");
    std::fs::create_dir_all(&runs).expect("a scratch directory");
    let bam = dir.join("large.bam");
    write_bam(&bam, SORTED_FILE_BYTES);
    assert!(std::fs::metadata(&bam).expect("the file").len() > 2 * (LIMIT_KIB << 10));

    // One reducer: more would sample the keys, which the port refuses. The sorted file goes to
    // the null device, where it is written in place, so the disk holds only the input and the runs;
    // no `.sbi`, which would be written beside it.
    let output = run_limited(
        "SortSamSpark --input \"$1\" --output /dev/null --sort-order queryname \
         --num-reducers 1 --create-output-bam-splitting-index false --tmp-dir \"$2\"",
        &[&bam, &runs],
    );
    let left = std::fs::read_dir(&runs).map(Iterator::count).unwrap_or(0);
    let _ = std::fs::remove_dir_all(&dir);
    assert_success(&output);
    // The runs are unnamed files, gone with the sort.
    assert_eq!(left, 0);
}
