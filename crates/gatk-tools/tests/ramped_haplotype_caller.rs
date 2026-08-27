//! Conformance for `RampedHaplotypeCaller`'s ramps against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/RampedHaplotypeCallerDump.java`.
//!
//! The caller itself is not compared: this tool is `HaplotypeCaller` with a different engine, and
//! the engine is milestone G3's. What is compared is the ramp format and the orderings, which are
//! the tool's own.
//!
//! # What this suite is for
//!
//!  * **an entry being named by its region's coordinates, or bare at the root**;
//!  * **`info.json` being written last**;
//!  * **the haplotype table's fixed header, 1-or-0 reference column and `Double.toString` score**;
//!  * **a reference haplotype sorting last**;
//!  * **the score being compared by sign**;
//!  * **the read comparator starting from the strand**;
//!  * **a size mismatch being refused before anything is compared**;
//!  * **and the bam index path replacing EVERY `.bam`.**

use gatk_corpus as corpus;
use gatk_tools::ramped_haplotype_caller::{
    bam_index_path, compare_haplotypes, compare_reads, entry_name, haplotype_table,
    loc_filename_suffix, read_supp_name, verify_haplotypes, verify_reads, Haplotype, Interval,
    Read, VerificationError, HAPLOTYPE_HEADER,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/ramped_haplotype_caller.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn line(text: &str, kind: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
        .to_string()
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

fn region() -> Interval {
    Interval {
        contig: "chr1".to_string(),
        start: 1000,
        end: 2000,
    }
}

/// The three haplotypes the dump built, in the order it built them.
fn haplotypes() -> Vec<Haplotype> {
    let of = |bases: &str, is_reference: bool, score: f64| Haplotype {
        contig: "chr1".to_string(),
        start: 1000,
        end: 1007,
        is_reference,
        cigar: "8M".to_string(),
        bases: bases.to_string(),
        score,
        alignment_start_hap_wrt_ref: 0,
    };
    vec![
        of("ACGTACGT", false, -12.5),
        of("ACGTTCGT", true, 0.0),
        of("ACGTACGA", false, -12.5000000001),
    ]
}

/// The three reads the dump built. `common_text` is not in the golden, so it is set to what the
/// ordering needs: the names it sorts by are the ones the golden reports.
fn reads() -> Vec<Read> {
    let of = |name: &str, start: i32, reverse: bool| Read {
        reverse_strand: reverse,
        common_text: name.to_string(),
        bases: "ACGTACGT".to_string(),
        base_qualities: vec![30; 8],
        soft_start: start,
        soft_end: start + 7,
        start,
        end: start + 7,
        unclipped_start: start,
        unclipped_end: start + 7,
    };
    vec![
        of("r-reverse-early", 1000, true),
        of("r-forward-late", 5000, false),
        of("r-forward-early", 1000, false),
    ]
}

#[test]
fn every_value_matches_the_golden() {
    let text = golden();
    let mut compared = 0;

    // The entry names, in write order: the region's two, then the root's, then info.json last.
    assert_eq!(
        line(&text, "entry", "ramp"),
        [
            entry_name(Some(&region()), "reads.txt"),
            entry_name(Some(&region()), "haplotypes.csv"),
            entry_name(None, "root.txt"),
            entry_name(None, "info.json"),
        ]
        .join(",")
    );
    compared += 1;

    // The haplotype table, whole.
    assert_eq!(
        haplotype_table(&haplotypes()),
        unescape(&line(&text, "content", "chr1-1000-2000/haplotypes.csv"))
    );
    compared += 1;

    // The two orderings.
    let mut sorted = haplotypes();
    sorted.sort_by(|a, b| compare_haplotypes(a, b).cmp(&0));
    assert_eq!(
        sorted
            .iter()
            .map(|haplotype| haplotype.bases.clone())
            .collect::<Vec<String>>()
            .join(","),
        line(&text, "order", "haplotypes")
    );
    compared += 1;
    let mut sorted_reads = reads();
    sorted_reads.sort_by(|a, b| compare_reads(a, b).cmp(&0));
    assert_eq!(
        sorted_reads
            .iter()
            .map(|read| read.common_text.clone())
            .collect::<Vec<String>>()
            .join(","),
        line(&text, "order", "reads")
    );
    compared += 1;

    // The three signs.
    let haplotypes = haplotypes();
    assert_eq!(
        compare_haplotypes(&haplotypes[0], &haplotypes[2])
            .signum()
            .to_string(),
        line(&text, "compare", "score-hair")
    );
    assert_eq!(
        compare_haplotypes(&haplotypes[1], &haplotypes[0])
            .signum()
            .to_string(),
        line(&text, "compare", "reference-last")
    );
    let reads = reads();
    assert_eq!(
        compare_reads(&reads[0], &reads[1]).signum().to_string(),
        line(&text, "compare", "strand-first")
    );
    compared += 3;

    // The derived names and paths.
    assert_eq!(
        loc_filename_suffix(&region()),
        line(&text, "name", "loc-suffix")
    );
    assert_eq!(
        read_supp_name("readname", true),
        line(&text, "name", "supp-true")
    );
    assert_eq!(
        read_supp_name("readname", false),
        line(&text, "name", "supp-false")
    );
    assert_eq!(bam_index_path("/x/reads.bam"), line(&text, "path", "plain"));
    assert_eq!(
        bam_index_path("/x/run.bam.d/reads.bam"),
        line(&text, "path", "twice")
    );
    assert_eq!(bam_index_path("/x/reads.cram"), line(&text, "path", "none"));
    compared += 6;

    assert_eq!(compared, 13, "the values the golden carries");
}

/// The region's coordinates become a directory, and a root entry has no prefix.
#[test]
fn an_entry_is_named_by_its_region() {
    assert_eq!(loc_filename_suffix(&region()), "chr1-1000-2000");
    assert_eq!(
        entry_name(Some(&region()), "haplotypes.csv"),
        "chr1-1000-2000/haplotypes.csv"
    );
    assert_eq!(entry_name(None, "info.json"), "info.json");

    // And `info.json` is the last entry written, after every region entry.
    let text = golden();
    let names: Vec<&str> = line(&text, "entry", "ramp")
        .split(',')
        .map(str::to_string)
        .collect::<Vec<String>>()
        .leak()
        .iter()
        .map(String::as_str)
        .collect();
    assert_eq!(names.last(), Some(&"info.json"));
    assert!(names[0].starts_with("chr1-1000-2000/"));
}

/// A fixed header, a 1-or-0 reference column and `Double.toString` for the score.
#[test]
fn the_haplotype_table_has_a_fixed_shape() {
    let table = haplotype_table(&haplotypes());
    let mut lines = table.lines();
    assert_eq!(lines.next(), Some(HAPLOTYPE_HEADER));
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len(), 3);
    // The reference haplotype's column is 1 and the others' are 0, in INPUT order: the table is
    // written unsorted.
    assert!(rows[0].contains(",0,8M,ACGTACGT,"));
    assert!(rows[1].contains(",1,8M,ACGTTCGT,"));
    // `Double.toString`, so a whole number keeps its `.0` and a long one keeps every digit.
    assert!(rows[1].ends_with(",0.0,0"));
    assert!(rows[2].ends_with(",-12.5000000001,0"));
}

/// A reference haplotype sorts last, and the score is compared by sign.
#[test]
fn a_reference_haplotype_sorts_last() {
    let haplotypes = haplotypes();
    let mut sorted = haplotypes.clone();
    sorted.sort_by(|a, b| compare_haplotypes(a, b).cmp(&0));
    assert!(sorted.last().expect("a haplotype").is_reference);
    assert!(!sorted[0].is_reference);

    // The two non-reference ones differ only by 1e-10 in score, and are still ordered.
    assert_eq!(compare_haplotypes(&haplotypes[0], &haplotypes[2]), 1);
    assert_eq!(compare_haplotypes(&haplotypes[2], &haplotypes[0]), -1);
    // A difference of exactly zero falls through to the bases.
    let same_score = Haplotype {
        bases: "AAAAAAAA".to_string(),
        ..haplotypes[0].clone()
    };
    assert_eq!(compare_haplotypes(&haplotypes[0], &same_score), 1);
}

/// The strand is the first key, so a reverse read at 1000 sorts after a forward one at 5000.
#[test]
fn the_read_comparator_starts_from_the_strand() {
    let reads = reads();
    assert!(reads[0].reverse_strand);
    assert!(!reads[1].reverse_strand);
    assert!(reads[0].start < reads[1].start);
    assert_eq!(
        compare_reads(&reads[0], &reads[1]),
        1,
        "the reverse one is later"
    );

    let mut sorted = reads.clone();
    sorted.sort_by(|a, b| compare_reads(a, b).cmp(&0));
    assert!(sorted.last().expect("a read").reverse_strand);
    // Two forward reads then fall through to the common text.
    assert_eq!(sorted[0].common_text, "r-forward-early");
    assert_eq!(sorted[1].common_text, "r-forward-late");
}

/// The size is checked before anything is compared, and otherwise the index is named.
#[test]
fn the_verification_refuses_a_size_before_a_difference() {
    let text = golden();
    let haplotypes = haplotypes();
    assert!(verify_haplotypes(&haplotypes, &haplotypes.clone()).is_ok());
    assert_eq!(refusal(&text, "haplotypes-same").1, "accepted");

    let (class, message) = refusal(&text, "haplotypes-size");
    assert_eq!(class, "java.lang.RuntimeException");
    let produced = verify_haplotypes(&haplotypes, &haplotypes[..2]).expect_err("a size mismatch");
    assert_eq!(
        produced,
        VerificationError::HaplotypeSize { left: 3, right: 2 }
    );
    assert_eq!(produced.message(), message);

    let mut changed = haplotypes[..2].to_vec();
    changed.push(Haplotype {
        bases: "TTTTTTTT".to_string(),
        ..haplotypes[2].clone()
    });
    let (class, message) = refusal(&text, "haplotypes-different");
    assert_eq!(class, "java.lang.RuntimeException");
    let produced = verify_haplotypes(&haplotypes, &changed).expect_err("a difference");
    assert_eq!(produced, VerificationError::HaplotypeIndex { index: 0 });
    assert_eq!(produced.message(), message);

    let reads = reads();
    assert!(verify_reads(&reads, &reads.clone()).is_ok());
    let (_, message) = refusal(&text, "reads-size");
    let produced = verify_reads(&reads, &reads[..2]).expect_err("a size mismatch");
    assert_eq!(produced, VerificationError::ReadSize { left: 3, right: 2 });
    assert_eq!(produced.message(), message);

    // Each side is sorted first, so the order the two arrived in does not matter.
    let mut shuffled = haplotypes.clone();
    shuffled.reverse();
    assert!(verify_haplotypes(&haplotypes, &shuffled).is_ok());
}

/// It is `String.replace`, so every `.bam` goes, and a path with none is unchanged.
#[test]
fn the_bam_index_path_replaces_every_occurrence() {
    assert_eq!(bam_index_path("/x/reads.bam"), "/x/reads.bai");
    assert_eq!(
        bam_index_path("/x/run.bam.d/reads.bam"),
        "/x/run.bai.d/reads.bai",
        "the directory's .bam goes too"
    );
    assert_eq!(bam_index_path("/x/reads.cram"), "/x/reads.cram");
    // A suffix change would have left the directory alone, which is what makes this a quirk.
    assert_ne!(
        bam_index_path("/x/run.bam.d/reads.bam"),
        "/x/run.bam.d/reads.bai"
    );
}
