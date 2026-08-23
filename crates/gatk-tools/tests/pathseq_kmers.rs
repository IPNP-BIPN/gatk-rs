//! Conformance for `PathSeqBuildKmers` against GATK 4.6.2.0, compared as the whole k-mer set.
//!
//! Golden from `tools/readfilter-conformance/PathSeqBuildKmersDump.java`.
//!
//! # What this suite is for
//!
//!  * **the encoding**, two bits per base with the first base in the high bits;
//!  * **canonicalisation**, and the even size that `canonical` itself refuses;
//!  * **masking after canonicalisation**, and the entries it collapses;
//!  * **the restart a bad base costs**, and the lower case that costs nothing;
//!  * **the spacing**, which rewinds the valid-base counter rather than the position;
//!  * **and the empty set no run ever writes**.
//!
//! The containers are not ported, so the Bloom filter run is checked the only way it can be: every
//! k-mer the filter was asked about is one the set holds, and it answered yes to all of them.

use gatk_corpus as corpus;
use gatk_tools::pathseq_kmers::{build, get_mask, parse_mask, Kmer, KmerError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pathseq_kmers.txt.gz"),
    )
}

fn contigs(text: &str) -> Vec<Vec<u8>> {
    ["one", "two"]
        .iter()
        .map(|name| {
            text.lines()
                .find_map(|line| line.strip_prefix(&format!("fixture\t{name}=")))
                .unwrap_or_else(|| panic!("the golden carries contig {name}"))
                .as_bytes()
                .to_vec()
        })
        .collect()
}

/// The set of one run, as the dump printed it: the long and the bases it spells.
fn kmers(text: &str, label: &str) -> Vec<(u64, String)> {
    let prefix = format!("kmer\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let (value, bases) = rest.split_once('\t').expect("a kmer row");
            (value.parse().expect("a long"), bases.to_string())
        })
        .collect()
}

fn count(text: &str, label: &str) -> usize {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("count\t{label}=")))
        .expect("the golden carries the count")
        .parse()
        .expect("a count")
}

fn refusal(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .expect("the golden carries the refusal")
        .to_string()
}

fn check(text: &str, label: &str, kmer_size: usize, spacing: usize, mask: &str) {
    let set = build(&contigs(text), kmer_size, spacing, mask).expect("a set");
    let ours: Vec<(u64, String)> = set
        .iter()
        .map(|value| (*value, Kmer(*value).bases(kmer_size)))
        .collect();
    assert_eq!(ours, kmers(text, label), "{label}");
    assert_eq!(ours.len(), count(text, label), "{label}: the count");
}

#[test]
fn five_bases_with_no_mask() {
    let text = golden();
    check(&text, "k5", 5, 1, "");
}

#[test]
fn seven_bases() {
    let text = golden();
    check(&text, "k7", 7, 1, "");
}

#[test]
fn masking_the_middle_base() {
    let text = golden();
    check(&text, "k5-mask-middle", 5, 1, "2");
    // The mask clears the third base's bits, so every entry spells an A there.
    assert!(kmers(&text, "k5-mask-middle")
        .iter()
        .all(|(_, bases)| bases.as_bytes()[2] == b'A'));
}

#[test]
fn masking_two_positions_collapses_entries() {
    let text = golden();
    check(&text, "k5-mask-first-and-last", 5, 1, "0,4");
    // Fourteen entries where the unmasked set had twenty: the mask merged six.
    assert!(count(&text, "k5-mask-first-and-last") < count(&text, "k5"));
}

#[test]
fn the_spacing_rewinds_the_valid_base_counter() {
    let text = golden();
    check(&text, "k5-spacing-three", 5, 3, "");
    check(&text, "k5-spacing-five", 5, 5, "");
    assert!(count(&text, "k5-spacing-five") < count(&text, "k5-spacing-three"));
}

#[test]
fn an_even_size_is_refused_by_canonical_itself() {
    let text = golden();
    let error = build(&contigs(&text), 4, 1, "").expect_err("a refusal");
    assert_eq!(error, KmerError::EvenKmerSize);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "even-k")
    );
}

#[test]
fn a_mask_index_outside_the_kmer() {
    let text = golden();
    let error = build(&contigs(&text), 5, 1, "9").expect_err("a refusal");
    assert_eq!(
        error,
        KmerError::InvalidMaskIndex {
            index: "9".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "mask-out-of-range")
    );
    assert!(parse_mask("4", 5).is_ok());
    assert!(parse_mask("5", 5).is_err());
}

#[test]
fn a_size_longer_than_every_contig_writes_nothing() {
    let text = golden();
    let error = build(&contigs(&text), 31, 1, "").expect_err("a refusal");
    assert_eq!(error, KmerError::EmptySet);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "k-longer-than-reference")
    );
}

#[test]
fn the_bloom_filter_holds_the_same_kmers() {
    let text = golden();
    let set = build(&contigs(&text), 5, 1, "").expect("a set");
    let asked: Vec<(u64, bool)> = text
        .lines()
        .filter_map(|line| line.strip_prefix("bloom\tk5-bloom\tcontains-"))
        .map(|rest| {
            let (value, answer) = rest.split_once('=').expect("a contains row");
            (
                value.parse().expect("a long"),
                answer.parse().expect("a boolean"),
            )
        })
        .collect();
    assert_eq!(asked.len(), set.len());
    for (value, answer) in asked {
        assert!(set.contains(&value), "the set holds {value}");
        assert!(answer, "the filter answered yes for {value}");
    }
}

#[test]
fn the_mask_is_built_from_the_start_of_the_kmer() {
    // Position 0 clears the two highest bits of a five-base k-mer.
    assert_eq!(get_mask(&[0], 5), Kmer(!(3u64 << 8)));
    assert_eq!(get_mask(&[4], 5), Kmer(!3u64));
    assert_eq!(get_mask(&[], 5), Kmer(!0u64));
}
