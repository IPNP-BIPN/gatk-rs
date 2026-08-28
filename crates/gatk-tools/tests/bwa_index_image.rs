//! Conformance for `BwaMemIndexImageCreator` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/BwaMemIndexImageCreatorDump.java`, which does NOT
//! hold the image's bytes: building one reference twice in a single process gives two files of the
//! same length whose differing bytes are in-process pointers.
//! `docs/pointers-that-reach-the-output.md` writes down why neither masking nor freezing works.
//!
//! # What this suite is for
//!
//!  * **the default output appending `.img` to the whole name**;
//!  * **the size being a function of the reference**;
//!  * **the case of the bases not being one**;
//!  * **the contig name being one**;
//!  * **the image never being byte-stable**;
//!  * **and the refusal naming the file and the reason.**

use gatk_corpus as corpus;
use gatk_tools::bwa_mem_index_image_creator::{
    cannot_read_reference, default_output, IMAGE_BYTES_ARE_REPRODUCIBLE, IMAGE_EXTENSION,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bwa_index_image.txt.gz"),
    )
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

fn size(text: &str, case: &str) -> i64 {
    field(text, "size", case)
        .unwrap_or_else(|| panic!("size/{case}"))
        .parse()
        .expect("a number")
}

/// The default output is the input's whole name plus `.img`.
#[test]
fn the_default_output_appends_the_extension() {
    let text = golden();
    assert_eq!(
        field(&text, "wrote", "default-output").as_deref(),
        Some("default-output.fasta.img")
    );
    assert_eq!(field(&text, "wrote", "plain").as_deref(), Some("plain.img"));
    assert_eq!(default_output("reference.fasta"), "reference.fasta.img");
    assert_ne!(default_output("reference.fasta"), "reference.img");
    assert_eq!(IMAGE_EXTENSION, ".img");
}

/// The size is a function of the reference, and the golden says which parts of it.
#[test]
fn the_size_answers_to_the_reference() {
    let text = golden();
    // Five repeats of eight bases against twenty.
    assert_eq!(size(&text, "plain"), 1333);
    assert_eq!(size(&text, "longer"), 1551);
    // The same bases in lower case are the same size: they are upper-cased before indexing.
    assert_eq!(size(&text, "lower-case"), size(&text, "plain"));
    // The contig NAME is in the image, so a longer name is a longer file.
    assert_eq!(size(&text, "renamed-contig"), 1334);
    // A second contig, and a run of Ns, each add to it.
    assert_eq!(size(&text, "two-contigs"), 1465);
    assert_eq!(size(&text, "with-ns"), 1349);
    // One base still produces an image.
    assert_eq!(size(&text, "one-base"), 1291);
    // And the size is the same on this laptop and on the runner: it is the bytes that move.
    assert_eq!(size(&text, "plain"), size(&text, "plain-again"));
    assert_eq!(size(&text, "default-output"), size(&text, "plain"));
}

/// No two builds of one reference agree, which is why the bytes are not in the golden.
#[test]
fn the_image_is_never_byte_stable() {
    let text = golden();
    let cases = [
        "plain",
        "plain-again",
        "lower-case",
        "longer",
        "two-contigs",
        "renamed-contig",
        "with-ns",
        "one-base",
        "default-output",
    ];
    for case in cases {
        assert_eq!(
            field(&text, "stable", case).as_deref(),
            Some("false"),
            "{case}"
        );
    }
    // The golden holds no bytes at all for this tool.
    assert!(!text.contains("\nimage\t"));
    // And the port says so where a caller would look.
    assert_eq!(
        IMAGE_BYTES_ARE_REPRODUCIBLE,
        field(&text, "stable", "plain").as_deref() == Some("true")
    );
}

/// The refusal names the file and the reason the native side gave.
#[test]
fn a_missing_reference_is_refused_by_name() {
    let text = golden();
    let refusal = field(&text, "error", "missing-fasta").expect("a refusal");
    assert!(
        refusal
            .starts_with("org.broadinstitute.hellbender.utils.bwa.CouldNotReadReferenceException:"),
        "{refusal}"
    );
    assert!(
        refusal.contains("cannot read the reference file '<dir>/no-such.fasta'"),
        "{refusal}"
    );
    assert!(refusal.contains("input file unre"), "{refusal}");
    assert_eq!(
        cannot_read_reference("<dir>/no-such.fasta", "input file unreadable"),
        "cannot read the reference file '<dir>/no-such.fasta': input file unreadable"
    );
}
