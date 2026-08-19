//! Conformance for `FastaReferenceMaker` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FastaReferenceMakerDump.java`.
//!
//! # What this suite is for
//!
//!  * **the output sequences are numbered**, the contig going into the description, so `>1
//!    chr1:1-43` is a sequence called `1`;
//!  * **a gap starts a new sequence** and an abutting pair does not, which is
//!    `withinDistanceOf(interval, 1)` and not an interval count;
//!  * **a contig boundary is a gap**, so a run with no `-L` writes one sequence per contig;
//!  * **the bases are the caching reader's**, so the IUPAC stretch comes out as `N`s and the
//!    soft-masked one comes out upper-cased: this tool does not copy its input;
//!  * **and the index and dictionary are the writer's**, including the `M5` it emits by default.

use gatk_corpus as corpus;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::fasta_reference_maker::{self, MakerError, DEFAULT_LINE_WIDTH};
use htsjdk_bam::fasta_writer::FastaOutputs;
use std::io::Write;

/// `ReferenceQueryDump.FASTA`, which the dump writes and indexes.
const FASTA: &str = ">chr1 first contig\n\
                    ACGTACGTACGT\n\
                    acgtNNNNacgt\n\
                    ACGTRYKMSWBD\n\
                    HVNACGT\n\
                    >chr2\n\
                    TTTTGGGGCCCC\n\
                    AAAATTTTGGGG\n";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/fasta_reference_maker.txt.gz"),
    )
}

fn row(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|value| value.to_string())
}

/// The dump's escaping.
fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

/// A fresh reference, in a directory of this call's own.
fn reference() -> ReferenceFileSource {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "gatk-rs-fastamaker-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let fasta = dir.join("ref.fasta");
    let mut file = std::fs::File::create(&fasta).expect("the fixture");
    file.write_all(FASTA.as_bytes()).expect("the fixture");
    drop(file);
    let mut fai = std::fs::File::create(dir.join("ref.fasta.fai")).expect("the index");
    writeln!(fai, "chr1\t43\t19\t12\t13").expect("the index");
    writeln!(fai, "chr2\t24\t72\t12\t13").expect("the index");
    drop(fai);
    ReferenceFileSource::open(&fasta).expect("the reference opens")
}

fn arguments(include: &[&str]) -> IntervalArguments {
    IntervalArguments {
        include: include.iter().map(|q| q.to_string()).collect(),
        ..Default::default()
    }
}

fn check(text: &str, label: &str, outputs: &FastaOutputs) {
    assert_eq!(
        escape(&String::from_utf8_lossy(&outputs.fasta)),
        row(text, "fasta", label).expect("a fasta row"),
        "{label}: the FASTA"
    );
    assert_eq!(
        escape(&outputs.index),
        row(text, "fai", label).expect("a fai row"),
        "{label}: the index"
    );
    assert_eq!(
        escape(&outputs.dictionary),
        row(text, "dict", label).expect("a dict row"),
        "{label}: the dictionary"
    );
}

#[test]
fn every_written_reference_matches_the_golden() {
    let text = golden();

    for (label, include, width) in [
        ("all", vec![], DEFAULT_LINE_WIDTH),
        ("one-interval", vec!["chr1:1-12"], DEFAULT_LINE_WIDTH),
        ("gap", vec!["chr1:1-5", "chr1:7-12"], DEFAULT_LINE_WIDTH),
        (
            "abutting",
            vec!["chr1:1-5", "chr1:6-12"],
            DEFAULT_LINE_WIDTH,
        ),
        (
            "two-contigs",
            vec!["chr1:1-6", "chr2:1-6"],
            DEFAULT_LINE_WIDTH,
        ),
        ("masked", vec!["chr1:13-24"], DEFAULT_LINE_WIDTH),
        ("iupac", vec!["chr1:25-36"], DEFAULT_LINE_WIDTH),
        ("narrow-lines", vec!["chr1:1-12"], 5),
        ("one-base", vec!["chr1:7-7"], DEFAULT_LINE_WIDTH),
    ] {
        let mut source = reference();
        let outputs = fasta_reference_maker::run(&mut source, &arguments(&include), width)
            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
        check(&text, label, &outputs);
    }
}

/// A line width of zero is the writer's refusal, and it happens before the reference is read.
#[test]
fn a_line_width_of_zero_is_refused() {
    let text = golden();
    assert_eq!(
        row(&text, "error", "zero-width").expect("an error row"),
        "java.lang.IllegalArgumentException"
    );

    let mut source = reference();
    let error = fasta_reference_maker::run(&mut source, &arguments(&["chr1:1-12"]), 0)
        .expect_err("a width of zero");
    let MakerError::Writer(writer) = error else {
        panic!("the writer refuses it, not the traversal");
    };
    assert_eq!(writer.java_class(), "java.lang.IllegalArgumentException");
}

/// The tool does not copy its input: the same twelve bases are `ACGTRYKMSWBD` in the file and
/// `ACGTNNNNNNNN` in the output.
#[test]
fn the_iupac_stretch_is_written_as_ns() {
    let text = golden();
    assert!(
        row(&text, "fasta", "iupac")
            .expect("a fasta row")
            .contains("ACGTNNNNNNNN"),
        "the golden holds the flattened bases"
    );
    let mut source = reference();
    let outputs =
        fasta_reference_maker::run(&mut source, &arguments(&["chr1:25-36"]), DEFAULT_LINE_WIDTH)
            .expect("the run succeeds");
    assert!(String::from_utf8_lossy(&outputs.fasta).contains("ACGTNNNNNNNN"));
}
