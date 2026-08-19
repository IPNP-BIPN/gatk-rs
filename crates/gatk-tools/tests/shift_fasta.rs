//! Conformance for `ShiftFasta` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ShiftFastaDump.java`.
//!
//! # What this suite is for
//!
//!  * **the default offset is half the contig**, per contig and by integer division;
//!  * **the two halves are appended separately**, tail first, and the line break falls where the
//!    writer's running count puts it rather than at the join;
//!  * **a contig that is not shifted is dropped**, from all four outputs at once;
//!  * **the chain id counts across contigs**, so the second contig's records are 3 and 4;
//!  * **the two interval files differ by the contig's parity**, and by nothing else;
//!  * **and an offset list of the wrong length is refused** before anything is written.

use gatk_corpus as corpus;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::shift_fasta::{self, ShiftError, ShiftOutputs};
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

/// `FastaReferenceWriter.DEFAULT_BASES_PER_LINE`, which is the tool's default too.
const DEFAULT_LINE_WIDTH: usize = 60;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/shift_fasta.txt.gz"),
    )
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

fn reference() -> ReferenceFileSource {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let dir = std::env::temp_dir().join(format!(
        "gatk-rs-shiftfasta-{}-{}",
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

fn check(text: &str, label: &str, outputs: &ShiftOutputs) {
    assert_eq!(
        escape(&String::from_utf8_lossy(&outputs.reference.fasta)),
        row(text, "fasta", label),
        "{label}: the shifted FASTA"
    );
    assert_eq!(
        escape(&outputs.reference.index),
        row(text, "fai", label),
        "{label}: the index"
    );
    assert_eq!(
        escape(&outputs.reference.dictionary),
        row(text, "dict", label),
        "{label}: the dictionary"
    );
    assert_eq!(
        escape(&outputs.chain),
        row(text, "chain", label),
        "{label}: the chain file"
    );
    assert_eq!(
        escape(&outputs.intervals),
        row(text, "intervals", label),
        "{label}: the intervals"
    );
    assert_eq!(
        escape(&outputs.shifted_intervals),
        row(text, "shifted", label),
        "{label}: the shifted intervals"
    );
}

#[test]
fn every_shifted_reference_matches_the_golden() {
    let text = golden();

    for (label, offsets, width) in [
        ("halves", vec![], DEFAULT_LINE_WIDTH),
        ("explicit", vec![5, 7], DEFAULT_LINE_WIDTH),
        // A zero offset and an offset of the whole contig, both of which drop that contig.
        ("skip-first", vec![0, 7], DEFAULT_LINE_WIDTH),
        ("skip-whole", vec![43, 7], DEFAULT_LINE_WIDTH),
        ("shift-one", vec![1, 1], DEFAULT_LINE_WIDTH),
        ("narrow-lines", vec![], 7),
    ] {
        let mut source = reference();
        let outputs = shift_fasta::run(&mut source, &offsets, width)
            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
        check(&text, label, &outputs);
    }
}

/// The skipped contig is absent from every output, not present and unshifted in any of them.
#[test]
fn a_contig_that_is_not_shifted_is_in_no_output() {
    let mut source = reference();
    let outputs = shift_fasta::run(&mut source, &[0, 7], DEFAULT_LINE_WIDTH).expect("the run");
    let fasta = String::from_utf8_lossy(&outputs.reference.fasta).to_string();
    for (what, text) in [
        ("the FASTA", fasta.as_str()),
        ("the index", outputs.reference.index.as_str()),
        ("the dictionary", outputs.reference.dictionary.as_str()),
        ("the chain file", outputs.chain.as_str()),
        ("the intervals", outputs.intervals.as_str()),
    ] {
        assert!(!text.contains("chr1"), "chr1 is still in {what}");
        assert!(text.contains("chr2"), "chr2 is missing from {what}");
    }
}

/// An offset list of the wrong size is refused with the reference's class and message.
#[test]
fn an_offset_list_of_the_wrong_length_is_refused() {
    let text = golden();
    let mut source = reference();
    let error = shift_fasta::run(&mut source, &[5], DEFAULT_LINE_WIDTH).expect_err("one of two");
    assert_eq!(
        format!("{}:{}", error.java_class(), escape(&error.message())),
        row(&text, "error", "wrong-length")
    );
    assert!(matches!(error, ShiftError::BadOffsetList { .. }));
}
