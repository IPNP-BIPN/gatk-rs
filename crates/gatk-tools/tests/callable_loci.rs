//! Conformance for `CallableLoci` against GATK 4.6.2.0, compared as the whole BED file and the
//! whole summary of every run.
//!
//! Golden from `tools/readfilter-conformance/CallableLociDump.java`.
//!
//! # What this suite is for
//!
//!  * **the run test never compares contigs**, so two `REF_N` stretches on different contigs whose
//!    coordinates run on come out as one line under the first contig's name;
//!  * **a deletion counts toward the QC depth whatever its base quality**, which the
//!    `high-base-quality` run isolates: four callable bases and nothing else;
//!  * **the poor-mapping-quality test comes before the depth tests**;
//!  * **the low count is `<=` and the passing count is `>=`**, so a run with the two thresholds
//!    equal counts a read at the threshold in both;
//!  * **`EXCESSIVE_COVERAGE` is tested on the raw depth**;
//!  * **and both files are fixed-width**: the BED is zero-based on its start, and the summary is a
//!    `%30s %d` table.

use gatk_corpus as corpus;
use gatk_engine::pileup::PileupElement;
use gatk_engine::read_pileup::pileup_from_reads;
use gatk_tools::callable_loci::{
    sample_refusal, state_at, write, Arguments, Element, OutputFormat,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/callable_loci.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

/// The harness's reference: chr1 upper-case ACGT to 60, lower-case to 120, sixty Ns to 180 and
/// upper-case again; chr2 all N but for an `AC` at 150 and 151.
fn reference_base(contig: &str, position: i32) -> u8 {
    match (contig, position) {
        ("chr1", 1..=60) => b"ACGT"[((position - 1) % 4) as usize],
        ("chr1", 61..=120) => b"acgt"[((position - 1) % 4) as usize],
        ("chr1", 121..=180) => b'N',
        ("chr1", _) => b"ACGT"[((position - 1) % 4) as usize],
        ("chr2", 150) => b'A',
        ("chr2", 151) => b'C',
        ("chr2", 61..=120) => b'n',
        ("chr2", _) => b'N',
        (other, _) => panic!("an unexpected contig: {other}"),
    }
}

fn read(name: &str, start: i32, cigar: &str, quality: u8, mapping_quality: u8) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length: usize = cigar
        .elements
        .iter()
        .filter(|element| {
            matches!(
                element.op,
                htsjdk_bam::cigar::Op::M | htsjdk_bam::cigar::Op::I
            )
        })
        .map(|element| element.length as usize)
        .sum();
    let mut tags = htsjdk_bam::tag::Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality,
        read_bases: vec![b'A'; length],
        base_qualities: vec![quality; length],
        cigar,
        tags,
        ..Default::default()
    }
}

/// The single-sample fixture, read for read.
fn fixture() -> Vec<BamRecord> {
    let mut reads = Vec::new();
    for i in 0..5 {
        reads.push(read(&format!("high{i}"), 1, "20M", 30, 60));
    }
    for i in 0..3 {
        reads.push(read(&format!("thin{i}"), 21, "20M", 30, 60));
    }
    for i in 0..5 {
        reads.push(read(&format!("low{i}"), 41, "20M", 30, 1));
    }
    for i in 0..2 {
        reads.push(read(&format!("mid{i}"), 41, "20M", 30, 20));
    }
    reads.push(read("del", 1, "4M4D12M", 5, 60));
    for i in 0..4 {
        reads.push(read(&format!("ns{i}"), 115, "20M", 30, 60));
    }
    reads
}

/// One pileup element, as the state depends on it.
fn element(element: &PileupElement<'_>) -> Element {
    Element {
        mapping_quality: i32::from(element.mapping_qual()),
        base_quality: i32::from(element.qual()),
        is_deletion: element.is_deletion(),
    }
}

/// The loci of one run, in traversal order, each with its state.
fn loci(
    reads: &[BamRecord],
    windows: &[(&str, i32, i32)],
    arguments: &Arguments,
) -> Vec<(String, i32, gatk_tools::callable_loci::State)> {
    let mut out = Vec::new();
    for (contig, from, to) in windows {
        for position in *from..=*to {
            // Only chr1 carries reads in the fixture, and a pileup on chr2 is empty.
            let pileup = if *contig == "chr1" {
                pileup_from_reads(contig, position, reads, |_| true, |_| true)
            } else {
                pileup_from_reads(contig, -1, reads, |_| true, |_| true)
            };
            let elements: Vec<Element> = pileup.elements.iter().map(element).collect();
            out.push((
                (*contig).to_string(),
                position,
                state_at(reference_base(contig, position), &elements, arguments),
            ));
        }
    }
    out
}

/// Every run of the dump: its intervals, its arguments and its output format.
fn run(label: &str) -> (Vec<(&'static str, i32, i32)>, Arguments, OutputFormat) {
    let base = Arguments::default();
    match label {
        "default" => (vec![("chr1", 1, 60)], base, OutputFormat::Bed),
        "n-run" => (vec![("chr1", 118, 125)], base, OutputFormat::Bed),
        "no-coverage" => (vec![("chr1", 200, 210)], base, OutputFormat::Bed),
        "state-per-base" => (vec![("chr1", 1, 40)], base, OutputFormat::StatePerBase),
        "contig-run-on" => (
            vec![("chr1", 171, 180), ("chr2", 181, 190)],
            base,
            OutputFormat::Bed,
        ),
        "contig-run-on-differing" => (
            vec![("chr1", 200, 210), ("chr2", 211, 220)],
            base,
            OutputFormat::Bed,
        ),
        "max-depth" => (
            vec![("chr1", 1, 60)],
            Arguments {
                max_depth: Some(3),
                ..base
            },
            OutputFormat::Bed,
        ),
        "min-depth-one" => (
            vec![("chr1", 1, 60)],
            Arguments {
                min_depth: 1,
                ..base
            },
            OutputFormat::Bed,
        ),
        "thresholds-equal" => (
            vec![("chr1", 1, 60)],
            Arguments {
                max_low_mapq: 20,
                min_mapping_quality: 20,
                ..base
            },
            OutputFormat::Bed,
        ),
        "poor-mapping-quality" => (
            vec![("chr1", 1, 60)],
            Arguments {
                min_depth_low_mapq: 1,
                max_low_mapq_fraction: 0.01,
                ..base
            },
            OutputFormat::Bed,
        ),
        "high-base-quality" => (
            vec![("chr1", 1, 60)],
            Arguments {
                min_base_quality: 60,
                min_depth: 1,
                ..base
            },
            OutputFormat::Bed,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let reads = fixture();
    let mut compared = 0;
    for label in [
        "default",
        "n-run",
        "no-coverage",
        "state-per-base",
        "contig-run-on",
        "contig-run-on-differing",
        "max-depth",
        "min-depth-one",
        "thresholds-equal",
        "poor-mapping-quality",
        "high-base-quality",
    ] {
        let (windows, arguments, format) = run(label);
        let (bed, summary) = write(&loci(&reads, &windows, &arguments), format);
        assert_eq!(bed, value(&text, "bed", label), "{label}: the BED");
        assert_eq!(
            summary,
            value(&text, "summary", label),
            "{label}: the summary"
        );
        compared += 1;
    }
    assert_eq!(compared, 11, "the golden's runs");
}

#[test]
fn the_cross_contig_run_is_one_line_under_the_first_contigs_name() {
    let text = golden();
    // Stated as itself, because it is the reference's bug and the port carries it on purpose.
    assert_eq!(
        value(&text, "bed", "contig-run-on"),
        "chr1\t170\t190\tREF_N\n"
    );
}

#[test]
fn two_samples_are_refused_by_name() {
    let text = golden();
    let refusal = sample_refusal(&["first".to_string(), "second".to_string()]);
    assert_eq!(
        refusal,
        text.lines()
            .find_map(|line| line.strip_prefix("error\ttwo-samples\t"))
            .expect("the golden carries the refusal")
    );
}
