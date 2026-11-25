//! Conformance for `VcfFormatConverter` against Picard 3.4.0, compared as the whole output file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/VcfFormatConverterDump.java`, which carries the run's
//! input as well as its output.
//!
//! # What this suite is for
//!
//!  * **`REQUIRE_INDEX` defaults to true**, so a file with no index beside it is refused by the
//!    reader before a record is read;
//!  * **`CREATE_INDEX` also defaults to true** and is what refuses a file with no contig lines;
//!  * **the two are independent**, so dropping the requirement does not stop the indexing from
//!    refusing, and turning the indexing off accepts the same file;
//!  * **the header comes back in the writer's order**, whatever order the input wrote it in;
//!  * **and a conversion is a rewrite**, so every record passes through the decoder and the
//!    encoder.
//!
//! # The BCF output
//!
//! The golden records it as a digest and a length, which `BCF2Writer`'s bytes are compared to: a
//! record from the VCF codec, through `htsjdk_vcf::bcf2`, as `gatk-cli` writes a `.bcf`. The round
//! trip back from BCF returns the same text as the plain conversion.

use gatk_corpus as corpus;
use gatk_tools::vcf_format_converter::{convert, Arguments, Input};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/vcf_format_converter.txt.gz"),
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

/// The input of each run: which file, whether an index sits beside it, and the arguments.
fn run(label: &str) -> (&'static str, bool, Arguments) {
    match label {
        "vcf-to-vcf" | "to-gz" | "to-bcf" => ("indexed", true, Arguments::default()),
        "empty" => ("empty", true, Arguments::default()),
        "no-index" => ("bare", false, Arguments::default()),
        "no-index-allowed" => (
            "bare",
            false,
            Arguments {
                require_index: false,
                ..Arguments::default()
            },
        ),
        "no-contigs" => (
            "no-contigs",
            false,
            Arguments {
                require_index: false,
                ..Arguments::default()
            },
        ),
        "no-contigs-no-index" => (
            "no-contigs",
            false,
            Arguments {
                require_index: false,
                create_index: false,
            },
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn converted(
    text: &str,
    label: &str,
) -> Result<String, gatk_tools::vcf_format_converter::ConvertError> {
    let (file, indexed, arguments) = run(label);
    let input = value(text, "input", file);
    convert(
        &Input {
            path: &format!("<dir>/{file}.vcf"),
            text: &input,
            indexed,
        },
        &arguments,
    )
}

#[test]
fn every_converted_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "vcf-to-vcf",
        "no-index-allowed",
        "no-contigs-no-index",
        "empty",
    ] {
        let ours = converted(&text, label).expect("a run the tool allows");
        assert_eq!(ours, value(&text, "converted", label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 4, "the golden's text outputs");
}

/// `BCF2Writer` over the indexed input, whose bytes the golden digests.
#[test]
fn the_bcf_output_matches_the_golden_digest() {
    use htsjdk_vcf::bcf2::writer::{Bcf2Writer, Origin};
    use sha2::{Digest, Sha256};

    let text = golden();
    let file = htsjdk_vcf::reader::read_vcf(&value(&text, "input", "indexed")).expect("parses");
    let mut writer = Bcf2Writer::new(&file.header, false).expect("a header BCF2Writer takes");
    let mut bytes = writer.header_bytes();
    for record in &file.records {
        bytes.extend(writer.encode(record, Origin::VcfCodec).expect("encodes"));
    }
    let digest: String = Sha256::digest(&bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect();
    let golden = text
        .lines()
        .find_map(|line| line.strip_prefix("sha256\tto-bcf="))
        .expect("the golden carries the BCF's digest");
    assert_eq!(format!("{digest}\t{}", bytes.len()), golden);
}

/// The block compressed run and the round trip through BCF both decompress to the plain
/// conversion's text.
#[test]
fn the_other_formats_carry_the_same_text() {
    let text = golden();
    let ours = converted(&text, "vcf-to-vcf").expect("a run the tool allows");
    assert_eq!(ours, value(&text, "converted", "to-gz"));
    assert_eq!(ours, value(&text, "converted", "bcf-to-vcf"));
}

#[test]
fn the_two_refusals_match_the_golden() {
    let text = golden();
    for label in ["no-index", "no-contigs"] {
        let error = converted(&text, label).expect_err("a refusal");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }
}

/// The header of the input is deliberately out of the writer's order, and the conversion sorts it.
#[test]
fn the_header_comes_back_in_the_writers_order() {
    let text = golden();
    let keys = |file: &str| -> Vec<String> {
        file.lines()
            .filter(|line| line.starts_with("##"))
            .map(|line| {
                line.trim_start_matches('#')
                    .split('=')
                    .next()
                    .expect("a key")
                    .to_string()
            })
            .collect()
    };
    let input = keys(&value(&text, "input", "indexed"));
    let ours = keys(&converted(&text, "vcf-to-vcf").expect("a run the tool allows"));
    assert_eq!(
        input,
        vec!["fileformat", "INFO", "FILTER", "FORMAT", "ALT", "contig"]
    );
    assert_eq!(
        ours,
        vec!["fileformat", "ALT", "FILTER", "FORMAT", "INFO", "contig"]
    );
}
