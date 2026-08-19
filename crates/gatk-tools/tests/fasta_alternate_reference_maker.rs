//! Conformance for `FastaAlternateReferenceMaker` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FastaAlternateReferenceMakerDump.java`.
//!
//! # What this suite is for
//!
//!  * **a deletion drops the bases after it**, through a counter that crosses `apply` calls;
//!  * **an insertion emits the whole alternate**, so one locus contributes three bases;
//!  * **a filtered record is skipped** and the reference base is written instead;
//!  * **the first concrete alternate is the one used**, so `*` is passed over;
//!  * **the mask writes `N` and loses a tie unless it is given priority**;
//!  * **`--use-iupac-sample` reads the genotype**, not the alternate;
//!  * **and a sample homozygous for a spanning deletion crashes the reference**, which the port
//!    reproduces rather than papering over.
//!
//! The records are built here rather than parsed: the VCF text lives in the dump, the parse path
//! has its own suites, and what this tool consumes is the decoded records.

use gatk_corpus as corpus;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::fasta_alternate_reference_maker::{
    self, AlternateArguments, AlternateError, ArgumentError,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};
use std::io::Write;

const FASTA: &str = ">chr1 first contig\n\
                    ACGTACGTACGT\n\
                    acgtNNNNacgt\n\
                    ACGTRYKMSWBD\n\
                    HVNACGT\n\
                    >chr2\n\
                    TTTTGGGGCCCC\n\
                    AAAATTTTGGGG\n";

const DEFAULT_LINE_WIDTH: usize = 60;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/fasta_alternate_reference_maker.txt.gz"),
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
        "gatk-rs-fastaalt-{}-{}",
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

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::create(bases.as_bytes(), is_ref).expect("a valid allele")
}

fn genotypes(first: &[&str], second: &[&str]) -> Vec<Genotype> {
    vec![
        Genotype::new(
            "NA1",
            first.iter().map(|a| allele(a, false)).collect::<Vec<_>>(),
        ),
        Genotype::new(
            "NA2",
            second.iter().map(|a| allele(a, false)).collect::<Vec<_>>(),
        ),
    ]
}

/// The dump's six records, in its order.
fn variants() -> Vec<VariantContext> {
    let mut records = Vec::new();

    // 2 C>T, a plain SNP.
    let mut snp = VariantContext::new("chr1", 2, vec![allele("C", true), allele("T", false)]);
    snp.filters = Some(Vec::new());
    snp.genotypes = genotypes(&["C", "T"], &["T", "T"]);
    records.push(snp);

    // 5 A>AGG, a simple insertion.
    let mut insertion =
        VariantContext::new("chr1", 5, vec![allele("A", true), allele("AGG", false)]);
    insertion.filters = Some(Vec::new());
    insertion.genotypes = genotypes(&["A", "AGG"], &["AGG", "AGG"]);
    records.push(insertion);

    // 8 TAC>T, a simple deletion of the two bases after 8.
    let mut deletion =
        VariantContext::new("chr1", 8, vec![allele("TAC", true), allele("T", false)]);
    deletion.filters = Some(Vec::new());
    deletion.genotypes = genotypes(&["TAC", "T"], &["T", "T"]);
    records.push(deletion);

    // 15 N>A, filtered.
    let mut filtered = VariantContext::new("chr1", 15, vec![allele("N", true), allele("A", false)]);
    filtered.filters = Some(vec!["LowQual".to_string()]);
    filtered.genotypes = genotypes(&["N", "A"], &["A", "A"]);
    records.push(filtered);

    // 20 N>*,C, whose first alternate is the spanning deletion.
    let mut spanning = VariantContext::new(
        "chr1",
        20,
        vec![allele("N", true), allele("*", false), allele("C", false)],
    );
    spanning.filters = Some(Vec::new());
    spanning.genotypes = genotypes(&["N", "*"], &["*", "*"]);
    records.push(spanning);

    // 30 N>G, a het for NA1 and a hom var for NA2.
    let mut het = VariantContext::new("chr1", 30, vec![allele("N", true), allele("G", false)]);
    het.filters = Some(Vec::new());
    het.genotypes = genotypes(&["N", "G"], &["G", "G"]);
    records.push(het);

    records
}

/// The mask: one record at 3 and one at 30, where a call also sits.
fn mask() -> Vec<VariantContext> {
    let mut first = VariantContext::new("chr1", 3, vec![allele("G", true), allele("A", false)]);
    first.filters = Some(Vec::new());
    let mut second = VariantContext::new("chr1", 30, vec![allele("N", true), allele("T", false)]);
    second.filters = Some(Vec::new());
    vec![first, second]
}

fn arguments(include: &[&str]) -> IntervalArguments {
    IntervalArguments {
        include: include.iter().map(|q| q.to_string()).collect(),
        ..Default::default()
    }
}

const SAMPLES: [&str; 2] = ["NA1", "NA2"];

fn samples() -> Vec<String> {
    SAMPLES.iter().map(|name| name.to_string()).collect()
}

#[test]
fn every_written_reference_matches_the_golden() {
    let text = golden();
    let records = variants();
    let masked = mask();

    let cases: Vec<(&str, Vec<&str>, AlternateArguments)> = vec![
        ("plain", vec!["chr1"], AlternateArguments::default()),
        (
            "masked",
            vec!["chr1"],
            AlternateArguments {
                mask: Some(&masked),
                ..Default::default()
            },
        ),
        (
            "mask-priority",
            vec!["chr1"],
            AlternateArguments {
                mask: Some(&masked),
                mask_priority: true,
                ..Default::default()
            },
        ),
        (
            "iupac-het",
            vec!["chr1"],
            AlternateArguments {
                iupac_sample: Some("NA1".to_string()),
                ..Default::default()
            },
        ),
        // A window that starts inside the deletion, where the counter has nothing to carry.
        (
            "after-deletion",
            vec!["chr1:9-20"],
            AlternateArguments::default(),
        ),
    ];

    for (label, include, extra) in cases {
        let mut source = reference();
        let outputs = fasta_alternate_reference_maker::run(
            &mut source,
            &arguments(&include),
            DEFAULT_LINE_WIDTH,
            &records,
            &extra,
            &samples(),
        )
        .unwrap_or_else(|error| panic!("{label}: {error:?}"));
        assert_eq!(
            escape(&String::from_utf8_lossy(&outputs.fasta)),
            row(&text, "fasta", label),
            "{label}: the FASTA"
        );
        assert_eq!(
            escape(&outputs.index),
            row(&text, "fai", label),
            "{label}: the index"
        );
    }
}

/// The reference crashes here, and the port crashes the same way: a sample homozygous for a
/// spanning deletion writes a space, and the writer refuses it.
#[test]
fn a_hom_var_spanning_deletion_is_refused_by_the_writer() {
    let text = golden();
    assert_eq!(
        row(&text, "error", "iupac-homvar"),
        "java.lang.IllegalArgumentException:the input sequence contains invalid base calls like:  "
    );

    let mut source = reference();
    let error = fasta_alternate_reference_maker::run(
        &mut source,
        &arguments(&["chr1"]),
        DEFAULT_LINE_WIDTH,
        &variants(),
        &AlternateArguments {
            iupac_sample: Some("NA2".to_string()),
            ..Default::default()
        },
        &samples(),
    )
    .expect_err("the space is refused");
    let AlternateError::Maker(gatk_tools::fasta_reference_maker::MakerError::Writer(writer)) =
        error
    else {
        panic!("the writer refuses it");
    };
    assert_eq!(
        format!("{}:{}", writer.java_class(), writer.message()),
        "java.lang.IllegalArgumentException:the input sequence contains invalid base calls like:  "
    );
}

#[test]
fn both_argument_checks_match_the_golden() {
    let text = golden();
    let masked = mask();

    for (label, extra, expected) in [
        (
            "priority-without-mask",
            AlternateArguments {
                mask_priority: true,
                ..Default::default()
            },
            ArgumentError::PriorityWithoutMask,
        ),
        (
            "unknown-sample",
            AlternateArguments {
                mask: Some(&masked),
                iupac_sample: Some("NOBODY".to_string()),
                ..Default::default()
            },
            ArgumentError::UnknownIupacSample,
        ),
    ] {
        let mut source = reference();
        let error = fasta_alternate_reference_maker::run(
            &mut source,
            &arguments(&["chr1"]),
            DEFAULT_LINE_WIDTH,
            &variants(),
            &extra,
            &samples(),
        )
        .expect_err("an argument check");
        let AlternateError::Argument(argument) = error else {
            panic!("{label}: an argument check, not a traversal failure");
        };
        assert_eq!(argument, expected, "{label}");
        assert_eq!(
            format!("{}:{}", argument.java_class(), escape(&argument.message())),
            row(&text, "error", label),
            "{label}"
        );
    }
}
