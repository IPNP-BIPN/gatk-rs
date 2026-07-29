//! Conformance for `CigarUtils`' clipping and `CigarBuilder` against GATK 4.6.2.0.
//!
//! The golden is produced by `tools/readfilter-conformance/CigarClipDump.java`: every cigar of a
//! sixteen-entry corpus clipped at every left and right boundary, with both clipping operators,
//! plus the alignment-start shift at every clip length, plus twelve element sequences fed to the
//! builder directly to reach its rewrites in isolation.
//!
//! `E` is the reference throwing. The builder validates rather than returning: a cigar that ends
//! up entirely soft-clipped, or whose sections are out of order, is a failure, and a port that
//! returned a cigar there would produce output where the reference produces none.

use gatk_corpus as corpus;
use gatk_engine::cigar_builder::CigarBuilder;
use gatk_engine::cigar_utils;
use htsjdk_bam::cigar::{CigarElement, Op};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cigar_clipping.txt.gz"),
    )
}

fn cigar(text: &str) -> htsjdk_bam::cigar::Cigar {
    htsjdk_bam::text_parse::parse_cigar(text).unwrap_or_else(|e| panic!("{text}: {e:?}"))
}

/// One cigar element from its text, `3M` to (3, M).
fn element(text: &str) -> CigarElement {
    let (length, op) = text.split_at(text.len() - 1);
    CigarElement {
        length: length.parse().unwrap(),
        op: match op {
            "M" => Op::M,
            "I" => Op::I,
            "D" => Op::D,
            "N" => Op::N,
            "S" => Op::S,
            "H" => Op::H,
            "P" => Op::P,
            "=" => Op::Eq,
            "X" => Op::X,
            other => panic!("unknown operator {other}"),
        },
    }
}

#[test]
fn every_clip_produces_the_cigar_the_reference_produces() {
    let text = golden();
    let mut clips = 0;
    let mut shifts = 0;
    let mut reverts = 0;
    let mut builds = 0;

    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts[0] {
            "clip" => {
                let (source, start, stop) = (parts[1], parts[2].parse().unwrap(), parts[3]);
                let stop: i32 = stop.parse().unwrap();
                let operator = match parts[4] {
                    "S" => Op::S,
                    "H" => Op::H,
                    other => panic!("unexpected clipping operator {other}"),
                };
                let ours = cigar_utils::clip_cigar(&cigar(source), start, stop, operator)
                    .map_or_else(|_| "E".to_string(), |c| c.to_text());
                assert_eq!(
                    ours, parts[5],
                    "{source} clipped {start}..{stop} with {}",
                    parts[4]
                );
                clips += 1;
            }
            "shift" => {
                let (source, clipped) = (parts[1], parts[2].parse().unwrap());
                assert_eq!(
                    cigar_utils::alignment_start_shift(&cigar(source), clipped).to_string(),
                    parts[3],
                    "{source} shift at {clipped}"
                );
                shifts += 1;
            }
            "revert" => {
                let source = parts[1];
                let ours = cigar_utils::revert_soft_clips(&cigar(source))
                    .map_or_else(|_| "E".to_string(), |c| c.to_text());
                assert_eq!(ours, parts[2], "{source} reverted");
                reverts += 1;
            }
            "build" => {
                let mut builder = CigarBuilder::default();
                let mut failed = false;
                for text in parts[1].split(',') {
                    if builder.add(element(text)).is_err() {
                        failed = true;
                        break;
                    }
                }
                let made = if failed {
                    Err(())
                } else {
                    builder.make(false).map_err(|_| ())
                };
                let ours = made
                    .as_ref()
                    .map_or_else(|_| "E".to_string(), |c| c.to_text());
                assert_eq!(ours, parts[2], "build {}", parts[1]);
                // The counters are read after make(), which is where the trailing removal happens.
                assert_eq!(
                    builder.leading_deletion_bases_removed().to_string(),
                    parts[3],
                    "build {}: leading deletion bases removed",
                    parts[1]
                );
                assert_eq!(
                    builder.trailing_deletion_bases_removed().to_string(),
                    parts[4],
                    "build {}: trailing deletion bases removed",
                    parts[1]
                );
                builds += 1;
            }
            _ => {}
        }
    }

    assert!(
        clips > 0 && shifts > 0 && reverts > 0 && builds > 0,
        "empty golden"
    );
    println!("{clips} clips, {shifts} shifts, {reverts} reverts, {builds} builds, all identical");
}
