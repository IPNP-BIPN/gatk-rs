//! Conformance for `StructuralVariantDiscoverer` against GATK 4.6.2.0, compared as the variants the
//! queryname-sorted run called and the refusal the other one made.
//!
//! Golden from `tools/readfilter-conformance/StructuralVariantDiscovererDump.java`.
//!
//! The VCF's annotations beyond the contig names are in the golden and are not reproduced: they
//! come from the evidence classes, which are not ported.
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
    check_sort_order, discover, passes_read_filters, signature, Alignment, DiscovererError,
    Signature, SortOrder,
};

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

/// The reads the golden reports, in the queryname order the file was written in.
fn alignments(text: &str) -> Vec<Alignment> {
    section(text, "bam", "reads")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("flags");
            let start: i32 = columns[3].parse().expect("a start");
            // Every cigar in this fixture is `<n>M` with optional soft clips, so the reference span
            // is the M run alone.
            let matched: i32 = {
                let cigar = columns[5];
                let mut length = 0;
                let mut digits = String::new();
                for character in cigar.chars() {
                    if character.is_ascii_digit() {
                        digits.push(character);
                    } else {
                        if character == 'M' {
                            length += digits.parse::<i32>().expect("a length");
                        }
                        digits.clear();
                    }
                }
                length
            };
            Alignment {
                contig_name: columns[0].to_string(),
                reference: columns[2].to_string(),
                start,
                end: start + matched - 1,
                reverse_strand: flags & 0x10 != 0,
                supplementary: flags & 0x800 != 0,
                secondary: flags & 0x100 != 0,
                unmapped: flags & 0x4 != 0,
                mapping_quality: columns[4].parse().expect("a mapping quality"),
            }
        })
        .collect()
}

/// The variants one run wrote, as position, id and alternate.
fn measured(text: &str, label: &str) -> Vec<(i32, String, String)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse().expect("a position"),
                columns[2].to_string(),
                columns[4].to_string(),
            )
        })
        .collect()
}

#[test]
fn every_call_matches_the_golden() {
    let text = golden();
    let produced: Vec<(i32, String, String)> = discover(&alignments(&text))
        .into_iter()
        .map(|variant| (variant.start, variant.id, variant.alternate))
        .collect();
    assert_eq!(produced, measured(&text, "default"));
    assert_eq!(produced.len(), 2, "two of the six contigs call anything");
}

/// A gap is a deletion, an overlap a tandem duplication.
#[test]
fn a_gap_is_a_deletion_and_an_overlap_a_duplication() {
    let text = golden();
    let all = alignments(&text);
    let of = |name: &str| -> Vec<Alignment> {
        all.iter()
            .filter(|alignment| alignment.contig_name == name)
            .cloned()
            .collect()
    };
    let deletion = of("ctg-del");
    assert_eq!(deletion.len(), 2);
    assert_eq!(
        signature(&deletion),
        Signature::Deletion {
            start: 10099,
            end: 10599
        }
    );
    let duplication = of("ctg-overlap");
    assert_eq!(duplication.len(), 2);
    assert!(
        matches!(signature(&duplication), Signature::TandemDuplication { .. }),
        "{:?}",
        signature(&duplication)
    );

    // And the golden names them.
    let calls = measured(&text, "default");
    assert_eq!(calls[0].1, "DEL_chr1_10099_10599");
    assert_eq!(calls[0].2, "<DEL>");
    assert_eq!(calls[1].2, "<DUP>");
    assert!(calls[1].1.starts_with("INS-DUPLICATION-TANDEM-EXPANSION_"));
}

/// A strand flip alone is not a signature.
#[test]
fn a_strand_flip_produces_nothing() {
    let text = golden();
    let all = alignments(&text);
    let inverted: Vec<Alignment> = all
        .iter()
        .filter(|alignment| alignment.contig_name == "ctg-inv")
        .cloned()
        .collect();
    assert_eq!(inverted.len(), 2, "it really has two pieces");
    assert_ne!(
        inverted[0].reverse_strand, inverted[1].reverse_strand,
        "and they really are on different strands"
    );
    // The same geometry on ONE strand would have been a deletion, so the strand is the only
    // difference.
    let same_strand: Vec<Alignment> = inverted
        .iter()
        .map(|alignment| Alignment {
            reverse_strand: false,
            ..alignment.clone()
        })
        .collect();
    assert!(matches!(
        signature(&same_strand),
        Signature::Deletion { .. }
    ));
    assert_eq!(signature(&inverted), Signature::None);
    // And nothing at 20000 reaches the output.
    assert!(!measured(&text, "default")
        .iter()
        .any(|(start, ..)| (20000..21000).contains(start)));
}

/// A lone alignment, a secondary one and an unmapped one all produce nothing.
#[test]
fn three_kinds_of_contig_produce_nothing() {
    let text = golden();
    let all = alignments(&text);
    let of = |name: &str| -> Vec<Alignment> {
        all.iter()
            .filter(|alignment| alignment.contig_name == name)
            .cloned()
            .collect()
    };
    assert_eq!(signature(&of("ctg-single")), Signature::None);
    assert_eq!(of("ctg-single").len(), 1, "one alignment, not two");

    let secondary = of("ctg-secondary");
    assert_eq!(secondary.len(), 1);
    assert!(!passes_read_filters(&secondary[0]));
    assert_eq!(signature(&secondary), Signature::None);

    let unmapped = of("ctg-unmapped");
    assert_eq!(unmapped.len(), 1);
    assert!(!passes_read_filters(&unmapped[0]));
    assert_eq!(signature(&unmapped), Signature::None);

    // The two deletion pieces pass the filters, so the filters are what removes the others.
    for alignment in of("ctg-del") {
        assert!(passes_read_filters(&alignment));
    }
    // None of the three reaches the output.
    for start in [30000, 50000, 60000] {
        assert!(!measured(&text, "default")
            .iter()
            .any(|(at, ..)| (start..start + 1000).contains(at)));
    }
}

/// Each call says which contig made it.
#[test]
fn the_contig_name_is_carried_onto_the_call() {
    let text = golden();
    let calls = discover(&alignments(&text));
    assert_eq!(calls[0].contig_names, vec!["ctg-del".to_string()]);
    assert_eq!(calls[1].contig_names, vec!["ctg-overlap".to_string()]);
    // The golden carries the same names in CTG_NAMES.
    let body = section(&text, "out", "default");
    assert!(body.contains("CTG_NAMES=ctg-del"));
    assert!(body.contains("CTG_NAMES=ctg-overlap"));
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
    assert_eq!(produced, DiscovererError::NotQuerynameSorted);
    assert_eq!(produced.message(), message);
    assert!(check_sort_order(SortOrder::Queryname).is_ok());
    // Unsorted is refused too, for the same reason.
    assert_eq!(
        check_sort_order(SortOrder::Unsorted).expect_err("unsorted"),
        DiscovererError::NotQuerynameSorted
    );
}
