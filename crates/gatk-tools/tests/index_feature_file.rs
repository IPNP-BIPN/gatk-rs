//! Conformance for `IndexFeatureFile` against GATK 4.6.2.0, compared as the whole index file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/IndexFeatureFileDump.java`, whose inputs and indexes
//! all travel as base64.
//!
//! # What this suite is for
//!
//!  * **the index bytes**, for the dynamic branch and the linear one, over a vcf and a bed;
//!  * **the same records under two names giving two different files**, `reads.vcf` a dynamic index
//!    and `reads.g.vcf` a linear one with a bin width of 128000;
//!  * **the default output**, appended to the whole name rather than replacing an extension;
//!  * **the refusal a block compressed input raises for an output that is not a `.tbi`**, whose
//!    message quotes the input;
//!  * **and the two refusals before any index is built**, an unreadable file and one in no
//!    supported format.
//!
//! # What the golden does not pin down here
//!
//! The tabix branch. Its bytes are in the golden for a later brick, and this port answers
//! `IndexKind::Tabix` for those names without writing anything. The suite asserts the choice and
//! the naming, not the bytes.
//!
//! The timestamp the reference embeds is zeroed by the dump, so the bytes compared here are the
//! reference's with those eight bytes cleared, which is what `Source::new` produces.

use gatk_corpus as corpus;
use gatk_tools::index_feature_file::{
    build, codec_for, default_output, index_kind, IndexKind, Refusal, Source,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/index_feature_file.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

fn bytes(text: &str, kind: &str, label: &str) -> Vec<u8> {
    corpus::decode_base64(&row(text, kind, label))
}

fn refusal(text: &str, label: &str) -> String {
    unescape(&row(text, "error", label))
}

/// The directory the dump ran in, which every index header quotes.
const DIR: &str = "/work/index-feature-file-dump";

/// The file each label was written under.
fn name(label: &str) -> &'static str {
    match label {
        "plain" => "reads.vcf",
        "gvcf" => "reads.g.vcf",
        "bed" => "regions.bed",
        "compressed" => "reads.vcf.gz",
        "unsorted" => "unsorted.vcf",
        "unknown" => "notes.txt",
        other => panic!("{other} is in the golden but not named here"),
    }
}

fn source(label: &str) -> Source {
    Source::new(&format!("{DIR}/{}", name(label)))
}

fn input(text: &str, label: &str) -> String {
    String::from_utf8(bytes(text, "input", label)).expect("a text fixture")
}

#[test]
fn every_index_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    // The three runs with no output given, plus the two explicit ones, which index the same file
    // and so must produce the same bytes.
    for (label, fixture) in [
        ("plain", "plain"),
        ("gvcf", "gvcf"),
        ("bed", "bed"),
        ("explicit", "plain"),
        ("explicit-tbi-name", "plain"),
    ] {
        let ours = build(&input(&text, fixture), &source(fixture), name(fixture))
            .expect("a run the tool allows");
        assert_eq!(ours, bytes(&text, "index", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's indexes");
}

/// The same records under two names are two different files, which is the whole point of the
/// extension rule.
#[test]
fn the_two_names_give_two_different_indexes() {
    let text = golden();
    assert_eq!(input(&text, "plain"), input(&text, "gvcf"));
    assert_ne!(
        bytes(&text, "index", "plain"),
        bytes(&text, "index", "gvcf")
    );
    assert_eq!(index_kind("reads.vcf"), IndexKind::Dynamic);
    assert_eq!(index_kind("reads.g.vcf"), IndexKind::Linear);
    assert_eq!(index_kind("reads.vcf.gz"), IndexKind::Tabix);
}

#[test]
fn the_default_output_is_appended_to_the_whole_name() {
    let text = golden();
    for label in ["plain", "gvcf", "bed", "compressed"] {
        assert_eq!(
            row(&text, "returned", label),
            format!("<dir>/{}", default_output(name(label))),
            "{label}"
        );
    }
    assert_eq!(default_output("reads.vcf"), "reads.vcf.idx");
    assert_eq!(default_output("reads.vcf.gz"), "reads.vcf.gz.tbi");
}

/// The port's own refusal is not one of the reference's, and says so.
#[test]
fn the_tabix_branch_refuses_in_the_ports_own_words() {
    let text = golden();
    // The reference writes a tabix index for this input, and the golden holds those bytes.
    assert!(!bytes(&text, "index", "compressed").is_empty());
    // The branch is chosen by the NAME, so the plain fixture's text under the compressed
    // fixture's name reaches it without this test having to decompress anything.
    let ours = build(
        &input(&text, "plain"),
        &source("compressed"),
        name("compressed"),
    )
    .expect_err("the port has no tabix writer");
    // No Java class, because the reference makes no such refusal: naming one of its exceptions
    // here would be a claim about the reference rather than about this port.
    assert_eq!(ours.java_class(), "");
    assert!(ours.message().contains("this port does not write yet"));
    // And it is NOT the refusal the reference makes for a mismatched extension, which is what a
    // covering array run against the binary caught: the two used to be the same variant, so a row
    // where the reference writes an index and the port cannot read as agreement.
    let theirs = Refusal::WrongIndexExtension {
        path: format!("<dir>/{}", name("compressed")),
    };
    assert_ne!(ours.message(), theirs.message());
}

#[test]
fn a_block_compressed_input_refuses_an_output_that_is_not_a_tbi() {
    let text = golden();
    let error = Refusal::WrongIndexExtension {
        path: format!("<dir>/{}", name("compressed")),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "compressed-wrong-extension")
    );
}

#[test]
fn the_refusals_before_any_index_match_the_golden() {
    let text = golden();
    // No codec for a name the manager does not know.
    assert_eq!(codec_for("notes.txt"), None);
    let error = Refusal::NoSuitableCodecs {
        path: "<dir>/notes.txt".to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unknown")
    );

    let error = Refusal::CouldNotReadInputFile {
        path: "<dir>/absent.vcf".to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "absent")
    );
}

/// A file whose last record is before its first is refused by the ordering check, and the port's
/// message is the reference's down to the two positions Tribble names.
#[test]
fn an_unsorted_file_is_refused_by_the_ordering_check() {
    let text = golden();
    let error = build(
        &input(&text, "unsorted"),
        &source("unsorted"),
        name("unsorted"),
    )
    .expect_err("the ordering refusal");
    let expected = refusal(&text, "unsorted");
    // The golden masks the dump's directory; the port has no reason to, so the mask is applied
    // here rather than invented in the port.
    let ours = format!("{}:{}", error.java_class(), error.message()).replace(DIR, "<dir>");
    assert_eq!(ours, expected);
}
