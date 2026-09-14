//! Conformance for `GtfToBed` against GATK 4.6.2.0, compared as the whole file it writes.
//!
//! Golden from `tools/readfilter-conformance/GtfToBedDump.java`.
//!
//! # What this suite is for
//!
//!  * **the one-based coordinates** of a file the tool calls a BED;
//!  * **a gene being as wide as its transcripts**, and `--use-basic-transcript` changing the gene
//!    rows because a transcript it drops never widens anything;
//!  * **`--sort-by-transcript` selecting rather than sorting**;
//!  * **the order**, which is the dictionary's contig index, then the start, then the key as a
//!    string;
//!  * **and the two refusals**, one from the traversal and one from the comparator.

use gatk_corpus as corpus;
use gatk_tools::gtf_to_bed::{run, Feature, GtfError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gtf_to_bed.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, prefix: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix))
            .unwrap_or_else(|| panic!("the golden carries {prefix}")),
    )
}

fn bed(text: &str, label: &str) -> String {
    field(text, &format!("bed\t{label}="))
}

fn refusal(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .expect("the golden carries the refusal")
        .to_string()
}

/// The contigs of a `.dict`, in the order it declares them.
fn dictionary(text: &str, name: &str) -> Vec<String> {
    field(text, &format!("dict\t{name}="))
        .lines()
        .filter(|line| line.starts_with("@SQ"))
        .map(|line| {
            line.split('\t')
                .find_map(|field| field.strip_prefix("SN:"))
                .expect("a sequence name")
                .to_string()
        })
        .collect()
}

/// Every gene and transcript line of the GTF; the exons the tool walks past are dropped here.
fn features(gtf: &str) -> Vec<Feature> {
    // The parser is the TOOL's, so the golden's own annotation is what covers it.
    gatk_tools::gtf_to_bed::parse_features(gtf)
}

fn check(text: &str, label: &str, by_transcript: bool, use_basic: bool) {
    let gtf = features(&field(text, "input\tannotation="));
    let dict = dictionary(text, "gtf");
    assert_eq!(
        run(&gtf, Some(&dict), by_transcript, use_basic).expect("an output"),
        bed(text, label),
        "{label}"
    );
}

#[test]
fn gene_rows_are_one_based_and_as_wide_as_their_transcripts() {
    let text = golden();
    check(&text, "genes", false, false);
    // The gene line said 100 200; its transcript reaches 50 250, and that is the row.
    assert!(bed(&text, "genes").starts_with("chr1\t50\t250\tbeta\n"));
}

#[test]
fn transcript_rows_carry_the_id_after_a_comma() {
    let text = golden();
    check(&text, "transcripts", true, false);
    assert!(bed(&text, "transcripts").contains("beta,TX_B1.1"));
}

#[test]
fn the_basic_flag_changes_the_gene_rows_too() {
    let text = golden();
    check(&text, "genes-basic-only", false, true);
    check(&text, "transcripts-basic-only", true, true);
    // gamma's only transcript is not basic, so with the flag on nothing widens it.
    assert!(bed(&text, "genes").contains("chr1\t300\t500\tgamma"));
    assert!(bed(&text, "genes-basic-only").contains("chr1\t300\t400\tgamma"));
}

#[test]
fn two_features_at_one_position_are_ordered_by_their_keys() {
    let text = golden();
    let genes = bed(&text, "genes");
    let rows: Vec<&str> = genes.lines().collect();
    // GENE_A.1 before GENE_C.1, both starting at 300.
    let alpha = rows
        .iter()
        .position(|row| row.ends_with("alpha"))
        .expect("alpha");
    let gamma = rows
        .iter()
        .position(|row| row.ends_with("gamma"))
        .expect("gamma");
    assert!(alpha < gamma);
}

#[test]
fn the_dictionary_is_required() {
    let text = golden();
    let gtf = features(&field(&text, "input\tannotation="));
    let error = run(&gtf, None, false, false).expect_err("a refusal");
    assert_eq!(error, GtfError::NoDictionary);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-dictionary")
    );
}

#[test]
fn a_contig_the_dictionary_does_not_know_is_refused_by_the_comparator() {
    let text = golden();
    let gtf = features(&field(&text, "input\tannotation="));
    let error =
        run(&gtf, Some(&dictionary(&text, "chr1only")), false, false).expect_err("a refusal");
    assert_eq!(
        error,
        GtfError::UnknownContig {
            contig: "chr2".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unknown-contig")
    );
}
