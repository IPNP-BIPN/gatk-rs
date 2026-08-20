//! Conformance for Picard's `UpdateVcfSequenceDictionary` against Picard 3.4.0, compared as the
//! whole output file of every run.
//!
//! Golden from `tools/readfilter-conformance/PicardUpdateVcfDictionaryDump.java`, which carries
//! each run's input and its dictionary.
//!
//! # What this suite is for
//!
//!  * **a contig the dictionary lacks is gone from the header while its records stay**, so the
//!    output declares fewer contigs than it uses;
//!  * **an empty dictionary is accepted**, leaving a file with no contig lines at all;
//!  * **the order is the dictionary's**, so reversing it reorders the header and nothing else;
//!  * **only `AS` survives into a contig line**, an `M5` and a `UR` having nowhere to go;
//!  * **and everything else in the header is untouched**, samples included.

use gatk_corpus as corpus;
use gatk_tools::picard_update_vcf_dictionary::{read_dictionary, update};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/picard_update_vcf_dictionary.txt.gz"),
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

#[test]
fn every_updated_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "both-contigs",
        "reversed",
        "missing-contig",
        "extra-contig",
        "with-attributes",
        "no-contigs-in",
        "empty-dictionary",
        "no-records",
    ] {
        let dictionary = read_dictionary(&value(&text, "dictionary", label));
        let ours =
            update(&value(&text, "input", label), &dictionary).expect("a run the tool allows");
        assert_eq!(ours, value(&text, "updated", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 8, "the golden's outputs");
}

#[test]
fn a_dictionary_with_no_sequences_leaves_a_header_with_no_contigs() {
    let text = golden();
    let dictionary = read_dictionary(&value(&text, "dictionary", "empty-dictionary"));
    assert!(dictionary.is_empty());
    let ours = update(&value(&text, "input", "empty-dictionary"), &dictionary).expect("a run");
    assert!(
        !ours.lines().any(|line| line.starts_with("##contig=")),
        "the output declares no contigs, and still carries records on two of them"
    );
    assert_eq!(
        ours.lines().filter(|line| line.starts_with("chr")).count(),
        2
    );
}
