//! Conformance for `ConvertCountsToDepthFile` against GATK 4.7.0.0, compared as the whole output
//! of every run: the text, or the BGZF bytes and their tabix index.
//!
//! Golden from `tools/readfilter-conformance/ConvertCountsToDepthFileDump.java`.
//!
//! # What this suite is for
//!
//!  * **the one-based counts written zero-based**, under the one sample the header names;
//!  * **`--sample-name` and `-n` winning over the read group's `SM`**;
//!  * **the block-compressed bytes at each level and the index beside them**;
//!  * **a header without a sample refused by the reader, wrapped twice**, whatever
//!    `--sample-name` says;
//!  * **and the two output-name refusals.**
//!
//! The refusals the engine makes before the tool runs (the parser's range, a missing or
//! wrong-typed input, `-L` without a dictionary) are the runner's, not this module's.

use gatk_corpus as corpus;
use gatk_tools::convert_counts_to_depth_file::{
    check_output, header_error_chain, read_counts, to_depth, write_block_compressed, write_text,
    BlockError, Counts,
};
use gatk_tools::sv_feature_codecs::Encoding;

/// Where the golden lands once a CI run has produced it.
const GOLDEN: &str = "tests/data/convert_counts_to_depth_file.txt.gz";

fn golden() -> Option<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(GOLDEN);
    // `golden-pending`: the first CI run on real x86-64 publishes the file this reads.
    path.exists().then(|| corpus::read_golden(&path))
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field<'a>(text: &'a str, kind: &str, label: &str) -> &'a str {
    let prefix = format!("{kind}\t{label}=");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
}

fn input(text: &str, label: &str) -> Counts {
    read_counts(&unescape(field(text, "input", label))).expect("a well-formed counts file")
}

/// `error` then every `cause`, as `class:message`.
fn refusals(text: &str, label: &str) -> Vec<String> {
    text.lines()
        .filter_map(|line| {
            let (kind, rest) = line.split_once('\t')?;
            let rest = rest.strip_prefix(label)?.strip_prefix('\t')?;
            (kind == "error" || kind == "cause").then(|| unescape(rest))
        })
        .collect()
}

fn converted(
    counts: &Counts,
    keep: impl Fn(&gatk_tools::print_read_counts::SimpleCount) -> bool,
) -> Vec<gatk_tools::condense_depth_evidence::DepthEvidence> {
    counts
        .records
        .iter()
        .filter(|count| keep(count))
        .map(to_depth)
        .collect()
}

#[test]
fn every_text_output_matches_the_golden() {
    let Some(text) = golden() else { return };
    for (label, sample) in [
        ("plain", None),
        ("sample-name", Some("delta")),
        ("sample-name-short", Some("epsilon")),
        ("plain-level-9", None),
        ("empty", None),
    ] {
        let counts = input(&text, label);
        let sample = sample.unwrap_or(&counts.sample);
        assert_eq!(
            write_text(sample, &converted(&counts, |_| true)),
            unescape(field(&text, "out", label)),
            "{label}"
        );
    }
    // -L chr1:150-160 -L chr1:180-250 -L chr2: the second count overlaps two and is written once.
    let counts = input(&text, "intervals");
    let kept = converted(&counts, |count| {
        (count.contig == "chr1" && count.start <= 250 && count.end >= 150) || count.contig == "chr2"
    });
    assert_eq!(
        write_text(&counts.sample, &kept),
        unescape(field(&text, "out", "intervals"))
    );
}

#[test]
fn every_block_compressed_output_and_index_match_the_golden() {
    let Some(text) = golden() else { return };
    // The index is deflated at htsjdk's own level, whatever --compression-level says: the dump
    // calls `instanceMain` directly, so `Main` never installs GATKConfig's two.
    let index_level = htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL;
    for (label, level) in [
        ("gz", 4),
        ("gz-level-1", 1),
        ("gz-level-9", 9),
        ("gz-empty", 4),
    ] {
        let counts = input(&text, label);
        let (data, index) = match write_block_compressed(
            &counts.sample,
            &converted(&counts, |_| true),
            &counts.header.sequences,
            level,
            index_level,
            htsjdk_bgzf::Deflater::Gkl,
        ) {
            Ok(written) => written,
            // Level 1 is ISA-L igzip, which only an x86-64 host reaches; CI compares it there.
            Err(BlockError::Deflate(_)) if !cfg!(target_arch = "x86_64") => continue,
            Err(error) => panic!("{label}: {error:?}"),
        };
        let plain = htsjdk_bgzf::read::decompress_all(&data).expect("BGZF");
        assert_eq!(
            String::from_utf8(plain).expect("text"),
            unescape(field(&text, "depth", label)),
            "{label}"
        );
        assert_eq!(
            data,
            corpus::decode_base64(field(&text, "gz", label)),
            "{label}"
        );
        assert_eq!(
            index,
            corpus::decode_base64(field(&text, "tbi", label)),
            "{label}"
        );
    }
}

#[test]
fn a_header_without_a_sample_is_refused_by_the_reader() {
    let Some(text) = golden() else { return };
    // With --sample-name or without: the reader refuses before the tool reads it.
    for label in ["no-sample", "no-sample-named", "no-read-groups"] {
        let inner = read_counts(&unescape(field(&text, "input", label)))
            .expect_err("a header refusal")
            .expect("of the header");
        let path = format!("<dir>/{label}/input.counts.tsv");
        let ours: Vec<String> = header_error_chain(&path, &inner)
            .into_iter()
            .map(|(class, message)| format!("{class}:{message}"))
            .collect();
        assert_eq!(ours, refusals(&text, label), "{label}");
    }
}

#[test]
fn the_two_output_refusals_match_the_golden() {
    let Some(text) = golden() else { return };
    let error = check_output(
        "<dir>/wrong-output-type/out.baf.txt",
        "file://<dir>/wrong-output-type/out.baf.txt",
    )
    .expect_err("another feature type");
    assert_eq!(
        vec![format!("{}:{}", error.java_class(), error.message())],
        refusals(&text, "wrong-output-type")
    );
    let error = check_output("<dir>/no-output-codec/out.txt", "").expect_err("no codec");
    assert_eq!(
        vec![format!("{}:{}", error.java_class(), error.message())],
        refusals(&text, "no-output-codec")
    );
    assert_eq!(
        check_output("x.RD.TXT.GZ", ""),
        Ok(Encoding::Text {
            block_compressed: true
        })
    );
}
