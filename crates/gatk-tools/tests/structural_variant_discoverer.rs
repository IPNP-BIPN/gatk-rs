//! Conformance for `StructuralVariantDiscoverer` against GATK 4.6.2.0, compared as the variants the
//! queryname-sorted run called and the refusal the other one made.
//!
//! Golden from `tools/readfilter-conformance/StructuralVariantDiscovererDump.java`.
//!
//! The golden's reads carry no bases, and its contigs are `ACGT` repeated, which is what the
//! duplication's `SEQ_ALT_HAPLOTYPE` reads back. The REF column comes from a reference the golden
//! does not carry either, so every column but that one is compared.
//!
//! # What this suite is for
//!
//!  * **a reference gap being a deletion and an overlap a tandem duplication**;
//!  * **a strand flip producing nothing at all**;
//!  * **a lone alignment, a secondary one and an unmapped one producing nothing**;
//!  * **the contig name being carried onto the call**;
//!  * **and a coordinate-sorted input being refused.**

use gatk_corpus as corpus;
use gatk_tools::structural_variant_discoverer::{
    check_sort_order, discover, passes_default_read_filters, write_vcf, DiscoveryArguments,
    SortOrder, NOT_QUERYNAME_SORTED,
};
use gatk_tools::sv_contig_alignments::{Cigar, ContigRead, Dictionary, Interval, SvError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/structural_variant_discoverer.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

/// One line of the golden's reads: the record and its two filter flags.
struct Line {
    read: ContigRead,
    secondary: bool,
}

/// The reads the golden reports, in the queryname order the file was written in.
fn lines(text: &str) -> Vec<Line> {
    section(text, "bam", "reads")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("flags");
            let start: i32 = columns[3].parse().expect("a start");
            let cigar = Cigar::parse(columns[5]).expect("a cigar");
            Line {
                read: ContigRead {
                    name: columns[0].to_string(),
                    unmapped: flags & 0x4 != 0,
                    reverse_strand: flags & 0x10 != 0,
                    supplementary: flags & 0x800 != 0,
                    contig: columns[2].to_string(),
                    start,
                    end: start + cigar.reference_length() - 1,
                    mapping_quality: columns[4].parse().expect("a mapping quality"),
                    bases: "ACGT".repeat(cigar.read_length() as usize / 4).into_bytes(),
                    cigar,
                    nm: None,
                    alignment_score: None,
                },
                secondary: flags & 0x100 != 0,
            }
        })
        .collect()
}

/// The reads that reach `apply`, through the tool's two default filters.
fn applied(lines: &[Line]) -> Vec<ContigRead> {
    lines
        .iter()
        .filter(|line| passes_default_read_filters(line.read.unmapped, line.secondary))
        .map(|line| line.read.clone())
        .collect()
}

fn dictionary() -> Dictionary {
    Dictionary {
        names: vec!["chr1".to_string()],
        lengths: vec![100_000],
        assemblies: vec![None],
    }
}

/// The records the tool writes for these reads, as the VCF body with the REF column blanked.
fn called(reads: &[ContigRead]) -> Vec<String> {
    let dictionary = dictionary();
    let canonical = dictionary.names.clone();
    let mut reference = |_: &Interval| -> Result<Vec<u8>, SvError> { Ok(b"N".to_vec()) };
    let variants = discover(
        reads,
        &dictionary,
        &canonical,
        None,
        &DiscoveryArguments::default(),
        &mut reference,
    )
    .expect("the golden's reads are all interpreted");
    body(&write_vcf(&variants, &dictionary, &[]))
}

/// The record lines of a VCF, with the REF column blanked.
fn body(vcf: &str) -> Vec<String> {
    vcf.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let mut columns: Vec<&str> = line.split('\t').collect();
            columns[3] = "";
            columns.join("\t")
        })
        .collect()
}

fn of(all: &[Line], name: &str) -> Vec<ContigRead> {
    applied(all)
        .into_iter()
        .filter(|read| read.name == name)
        .collect()
}

#[test]
fn every_call_matches_the_golden() {
    let text = golden();
    let produced = called(&applied(&lines(&text)));
    assert_eq!(produced, body(&section(&text, "out", "default")));
    assert_eq!(produced.len(), 2, "two of the six contigs call anything");
}

/// A gap is a deletion, an overlap a tandem duplication.
#[test]
fn a_gap_is_a_deletion_and_an_overlap_a_duplication() {
    let text = golden();
    let all = lines(&text);
    let deletion = called(&of(&all, "ctg-del"));
    assert_eq!(deletion.len(), 1);
    let columns: Vec<&str> = deletion[0].split('\t').collect();
    assert_eq!(columns[1], "10099");
    assert_eq!(columns[2], "DEL_chr1_10099_10599");
    assert_eq!(columns[4], "<DEL>");
    let duplication = called(&of(&all, "ctg-overlap"));
    assert_eq!(duplication.len(), 1);
    let columns: Vec<&str> = duplication[0].split('\t').collect();
    assert_eq!(columns[4], "<DUP>");
    assert!(columns[2].starts_with("INS-DUPLICATION-TANDEM-EXPANSION_"));
}

/// A strand flip alone is not a signature: the reverse-strand piece of `ctg-inv` claims the same
/// contig bases as the forward one, so it is dropped as contained, and nothing is left to pair.
#[test]
fn a_strand_flip_produces_nothing() {
    let text = golden();
    let all = lines(&text);
    let inverted = of(&all, "ctg-inv");
    assert_eq!(inverted.len(), 2, "it really has two pieces");
    assert_ne!(
        inverted[0].reverse_strand, inverted[1].reverse_strand,
        "and they really are on different strands"
    );
    assert!(called(&inverted).is_empty());
    // The same geometry on ONE strand is a deletion, so the strand is the only difference.
    let same_strand: Vec<ContigRead> = inverted
        .iter()
        .map(|read| ContigRead {
            reverse_strand: false,
            ..read.clone()
        })
        .collect();
    let deletion = called(&same_strand);
    assert_eq!(deletion.len(), 1);
    assert!(deletion[0].contains("<DEL>"));
}

/// A lone alignment, a secondary one and an unmapped one all produce nothing.
#[test]
fn three_kinds_of_contig_produce_nothing() {
    let text = golden();
    let all = lines(&text);
    let single = of(&all, "ctg-single");
    assert_eq!(single.len(), 1, "one alignment, not two");
    assert!(called(&single).is_empty());
    // The filters are what remove the other two, before the tool sees them.
    assert!(of(&all, "ctg-secondary").is_empty());
    assert!(of(&all, "ctg-unmapped").is_empty());
    let secondary = all
        .iter()
        .find(|line| line.read.name == "ctg-secondary")
        .expect("the golden carries it");
    assert!(!passes_default_read_filters(
        secondary.read.unmapped,
        secondary.secondary
    ));
    let unmapped = all
        .iter()
        .find(|line| line.read.name == "ctg-unmapped")
        .expect("the golden carries it");
    assert!(!passes_default_read_filters(
        unmapped.read.unmapped,
        unmapped.secondary
    ));
    // And the deletion's two pieces pass them.
    assert_eq!(of(&all, "ctg-del").len(), 2);
}

/// Each call says which contig made it.
#[test]
fn the_contig_name_is_carried_onto_the_call() {
    let text = golden();
    let calls = called(&applied(&lines(&text)));
    assert!(calls[0].contains("CTG_NAMES=ctg-del;"));
    assert!(calls[1].contains("CTG_NAMES=ctg-overlap;"));
    let golden_body = section(&text, "out", "default");
    assert!(golden_body.contains("CTG_NAMES=ctg-del"));
    assert!(golden_body.contains("CTG_NAMES=ctg-overlap"));
}

/// The tool walks consecutive records of one name, so anything else is refused.
#[test]
fn a_coordinate_sorted_input_is_refused() {
    let text = golden();
    let (class, message) = refusal(&text, "coordinate-sorted");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException"
    );
    let produced = check_sort_order(SortOrder::Coordinate).expect_err("coordinate order");
    assert_eq!(produced.class, class);
    assert_eq!(produced.message, message);
    assert_eq!(NOT_QUERYNAME_SORTED, message);
    assert!(check_sort_order(SortOrder::Queryname).is_ok());
    // Unsorted is refused too, for the same reason.
    assert_eq!(
        check_sort_order(SortOrder::Unsorted).expect_err("unsorted"),
        produced
    );
}
