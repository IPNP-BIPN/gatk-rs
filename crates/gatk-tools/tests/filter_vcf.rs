//! Conformance for `FilterVcf` against Picard 3.4.0, compared as the whole output file of every
//! run.
//!
//! Golden from `tools/readfilter-conformance/FilterVcfDump.java`, which carries each run's input as
//! well as its output.
//!
//! # What this suite is for
//!
//!  * **a passing genotype's `FT` never reaches the file**, so FT appears only on the records
//!    where something was filtered;
//!  * **a record whose genotypes were all filtered is `AllGtsFiltered`**;
//!  * **`LowQD` skips a record with no QD** however high the threshold, and `StrandBias` reaches
//!    one with no FS only when the threshold is negative;
//!  * **the allele-balance filter groups by the genotype's alleles**, so two samples with the same
//!    het call share one tally;
//!  * **and the four FILTER lines and the FT line are added to the header** whatever the file did.

use gatk_corpus as corpus;
use gatk_tools::filter_vcf::{filter, Thresholds};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/filter_vcf.txt.gz"),
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

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// The thresholds of each run.
fn thresholds(label: &str) -> Thresholds {
    let base = Thresholds::default();
    match label {
        "defaults" | "no-records" | "no-contigs" => base,
        "min-ab" => Thresholds {
            min_ab: 0.4,
            ..base
        },
        "max-fs" => Thresholds {
            max_fs: 0.5,
            ..base
        },
        "min-qd" => Thresholds {
            min_qd: 25.0,
            ..base
        },
        "min-gq" => Thresholds { min_gq: 50, ..base },
        "min-dp" => Thresholds { min_dp: 10, ..base },
        "negative-fs" => Thresholds {
            max_fs: -1.0,
            ..base
        },
        "everything" => Thresholds {
            min_ab: 0.45,
            max_fs: 0.5,
            min_qd: 25.0,
            min_gq: 50,
            min_dp: 10,
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_filtered_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "defaults",
        "min-ab",
        "max-fs",
        "min-qd",
        "min-gq",
        "min-dp",
        "negative-fs",
        "everything",
        "no-records",
    ] {
        let ours = filter(&value(&text, "input", label), &thresholds(label))
            .expect("a run the tool allows");
        assert_eq!(ours, value(&text, "filtered", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 9, "the golden's outputs");
}

#[test]
fn a_header_with_no_contigs_is_the_one_refusal() {
    let text = golden();
    let error = filter(
        &value(&text, "input", "no-contigs"),
        &thresholds("no-contigs"),
    )
    .expect_err("the dictionary refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-contigs")
    );
}

#[test]
fn all_gts_filtered_counts_only_variant_genotypes() {
    let input = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
        ##contig=<ID=chr1,length=1000>\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n\
        chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT:GQ:DP\t0/0:5:5\t0/0:5:5\n\
        chr1\t20\t.\tA\tG\t.\tPASS\t.\tGT:GQ:DP\t0/1:5:5\t0/0:50:50\n";
    let thresholds = Thresholds {
        min_gq: 30,
        ..Thresholds::default()
    };
    let output = filter(input, &thresholds).expect("filtered");
    let sites: Vec<&str> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(6).unwrap())
        .collect();
    // Every hom-ref genotype filtered is not a filtered site; the one het filtered is.
    assert_eq!(sites, ["PASS", "AllGtsFiltered"]);
}

/// The record-at-a-time run refuses what the whole-file read refuses, in its order: a record the
/// reader refuses, however late, beats a missing dictionary and an earlier record the encoder
/// refuses.
#[test]
fn a_late_read_refusal_beats_an_earlier_one() {
    let header = |contig: &str| {
        format!("##fileformat=VCFv4.2\n{contig}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    };
    let unknown_key = "chr1\t10\t.\tA\tG\t.\tPASS\tXX=1\n";
    let bad_position = "chr1\tten\t.\tA\tG\t.\tPASS\t.\n";
    let read_refusal = |text: &str| {
        let failure = htsjdk_vcf::reader::read_vcf(text).expect_err("a record the reader refuses");
        format!("{}:{}", failure.error.class(), failure.error.message())
    };
    let ours = |text: &str| {
        let error = filter(text, &Thresholds::default()).expect_err("refused");
        format!("{}:{}", error.java_class(), error.message())
    };

    let no_contigs = format!("{}{bad_position}", header(""));
    assert_eq!(ours(&no_contigs), read_refusal(&no_contigs));

    let contig = "##contig=<ID=chr1,length=100>\n";
    let encoder_then_reader = format!("{}{unknown_key}{bad_position}", header(contig));
    assert_eq!(
        ours(&encoder_then_reader),
        read_refusal(&encoder_then_reader)
    );

    let encoder_only = format!("{}{unknown_key}", header(contig));
    assert!(ours(&encoder_only).starts_with("java.lang.IllegalStateException:Key XX found"));
    let no_contigs_only = format!("{}{unknown_key}", header(""));
    assert!(ours(&no_contigs_only).starts_with("picard.PicardException:"));
}
