//! Conformance for `SplitNCigarReads` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/SplitNCigarReadsDump.java`. The five input BAMs and
//! their indexes travel in full, base64, so the port opens the same bytes the reference did, and
//! each output BAM is compared byte for byte.
//!
//! # What this suite is for
//!
//!  * **a read with k N elements becomes k+1 reads**, each keeping every base and soft clipping the
//!    rest;
//!  * **a section beside a deletion is trimmed before the clip**, and `N-D-N` is refused by name;
//!  * **a cigar ending in N loses it, one beginning with N is passed through untouched**;
//!  * **the mapping quality transform rewrites 255 and only 255**, unless it is skipped;
//!  * **the default filter is ALLOW_ALL_READS**, so a malformed read is still written;
//!  * **an MC tag is rewritten to what the mate's cigar will become**;
//!  * **a secondary alignment is not split** unless asked for;
//!  * **an overhang across another read's splice is clipped** where the bases disagree;
//!  * **and the output is not in traversal order**, which is what the not-presorted writer is for.

use gatk_corpus as corpus;
use gatk_engine::overhang_fixing_manager::OverhangArguments;
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::sam_output::Options;
use gatk_tools::split_n_cigar_reads::{self, refactor_ndn_to_n, SplitArguments, SplitToolError};
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/split_n_cigar_reads.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn reference_bases(text: &str) -> String {
    rows(text, "reference")[0][0].to_string()
}

fn install_fixtures(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
    for row in rows(text, "fixtureindex") {
        if row[1] == "absent" {
            continue;
        }
        std::fs::write(
            dir.join(format!("{}.bai", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture index");
    }
}

fn parse_cigar(text: &str) -> Cigar {
    if text == "*" {
        return Cigar::new(Vec::new());
    }
    let mut elements = Vec::new();
    let mut length = 0u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
            continue;
        }
        let op = match byte {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            other => panic!("no cigar operator {}", other as char),
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    Cigar::new(elements)
}

/// The fixture and the arguments each labelled run used.
fn configuration(label: &str) -> (&str, SplitArguments) {
    let default = SplitArguments::default();
    match label {
        "splits" => ("splits", default),
        "qualities" => ("qualities", default),
        "qualities-skip-mq" => (
            "qualities",
            SplitArguments {
                skip_mq_transform: true,
                ..default
            },
        ),
        "pairs" => ("pairs", default),
        "pairs-secondary" => (
            "pairs",
            SplitArguments {
                process_secondary_alignments: true,
                ..default
            },
        ),
        "overhangs" => ("overhangs", default),
        "overhangs-not-fixed" => (
            "overhangs",
            SplitArguments {
                overhang: OverhangArguments {
                    do_not_fix_overhangs: true,
                    ..OverhangArguments::default()
                },
                ..default
            },
        ),
        "overhangs-strict" => (
            "overhangs",
            SplitArguments {
                overhang: OverhangArguments {
                    max_mismatches_in_overhang: 0,
                    ..OverhangArguments::default()
                },
                ..default
            },
        ),
        "ndn" => ("ndn", default),
        "ndn-refactored" => (
            "ndn",
            SplitArguments {
                refactor_ndn_cigar_reads: true,
                ..default
            },
        ),
        other => panic!("no run {other}"),
    }
}

/// One labelled run of the tool over its fixture.
fn run(
    dir: &std::path::Path,
    label: &str,
    chr1: &str,
) -> Result<(Vec<u8>, Option<Vec<u8>>), SplitToolError> {
    let (fixture, arguments) = configuration(label);
    let bam = dir.join(format!("{fixture}.bam"));
    let bai = dir.join(format!("{fixture}.bai"));
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");

    let command_line = command_line(label);
    let options = Options {
        command_line: &command_line,
        ..Options::default()
    };
    // `ALLOW_ALL_READS`, which is this tool's whole default filter.
    let filter = |_: &BamRecord| true;
    let mut reference = |contig: &str, start: i32, end: i32| -> Result<Vec<u8>, String> {
        if contig != "chr1" {
            return Err(format!("unknown contig {contig}"));
        }
        Ok(chr1.as_bytes()[(start - 1) as usize..end as usize].to_vec())
    };

    split_n_cigar_reads::split_n_cigar_reads(&source, &arguments, &options, &filter, &mut reference)
}

/// The `@PG` command line the reference recorded, taken from the golden rather than rebuilt.
///
/// Barclay expands every argument including the defaults, so the line is an input to the port and
/// not something it could invent.
fn command_line(label: &str) -> String {
    let text = golden();
    rows(&text, "commandline")
        .into_iter()
        .find(|row| row[0] == label)
        .map(|row| row[1].to_string())
        .unwrap_or_default()
}

#[test]
fn every_ndn_refactor_is_the_reference() {
    let text = golden();
    let cases = rows(&text, "ndn");
    assert_eq!(cases.len(), 7, "seven cigars through the refactor");
    for row in cases {
        assert_eq!(
            refactor_ndn_to_n(&parse_cigar(row[0])).to_text(),
            row[1],
            "ndn/{}",
            row[0]
        );
    }
}

#[test]
fn every_output_bam_is_the_reference_byte_for_byte() {
    let text = golden();
    let chr1 = reference_bases(&text);
    let dir = std::env::temp_dir().join(format!("gatk-rs-split-n-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let outputs = rows(&text, "output");
    assert!(outputs.len() >= 9, "nine runs finish");

    for row in &outputs {
        let label = row[0];
        let (bam, index) = run(&dir, label, &chr1).unwrap_or_else(|error| {
            panic!("{label} was refused: {}", error_message(&error));
        });
        assert_eq!(bam, corpus::decode_base64(row[1]), "output/{label}");

        let expected_index = rows(&text, "index")
            .into_iter()
            .find(|other| other[0] == label)
            .map(|other| other[1].to_string())
            .expect("every finished run dumps its index");
        match (index, expected_index.as_str()) {
            (None, "absent") => {}
            (Some(ours), expected) => {
                assert_eq!(ours, corpus::decode_base64(expected), "index/{label}")
            }
            (None, _) => panic!("{label} produced no index and the reference did"),
        }
    }
}

/// `N-D-N` is refused with the reference's own class and message.
#[test]
fn the_empty_section_is_refused_by_name() {
    let text = golden();
    let chr1 = reference_bases(&text);
    let dir = std::env::temp_dir().join(format!("gatk-rs-split-n-error-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 1, "one refusal");
    for row in errors {
        let error = run(&dir, row[0], &chr1).expect_err("this run is refused");
        assert_eq!(
            format!(
                "java.lang.IllegalArgumentException:{}",
                error_message(&error)
            ),
            row[1],
            "error/{}",
            row[0]
        );
    }
}

fn error_message(error: &SplitToolError) -> String {
    match error {
        SplitToolError::Split(split) => split.message(),
        SplitToolError::Reads(reads) => format!("{reads:?}"),
    }
}
