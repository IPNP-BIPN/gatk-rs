//! Conformance for `CreateHadoopBamSplittingIndex` against GATK 4.6.2.0, compared as the bytes of
//! the index it wrote.
//!
//! Golden from `tools/readfilter-conformance/CreateHadoopBamSplittingIndexDump.java`, which holds
//! each index as base64 and again as the fields read back out of it, so a mismatch says which
//! number moved rather than only that the bytes differ.
//!
//! # What this suite is for
//!
//!  * **the granularity counting records and not bytes**;
//!  * **the last entry being where the next record would have gone**;
//!  * **an empty BAM falling back to the file's length for it**;
//!  * **the default output appending `.sbi` rather than replacing an extension**;
//!  * **the `.bai` companion being named by replacing one**;
//!  * **only the `.bai` path caring about the sort order**;
//!  * **a granularity of nought or less being refused before anything is opened**;
//!  * **a file that is not a BAM being refused by its extension**;
//!  * **and the md5 and uuid being written as zeroes.**

use gatk_corpus as corpus;
use gatk_tools::create_hadoop_bam_splitting_index::{
    assert_granularity, assert_is_bam, bai_companion, default_output, make_file_pointer, offsets,
    write, DEFAULT_GRANULARITY, GRANULARITY_MESSAGE, SBI_MAGIC,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/splitting_index.txt.gz"),
    )
}

fn field(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

/// One case's `fields` row, split into its parts.
fn fields(text: &str, label: &str) -> Vec<String> {
    field(text, "fields", label)
        .unwrap_or_else(|| panic!("fields/{label}"))
        .split(',')
        .map(str::to_string)
        .collect()
}

fn number(parts: &[String], index: usize) -> u64 {
    parts[index].parse().expect("a number")
}

/// The index's own bytes, rebuilt from the numbers the reference wrote.
#[test]
fn the_bytes_are_the_ported_layout() {
    let text = golden();
    for label in [
        "granularity-default",
        "granularity-two",
        "granularity-one",
        "granularity-above-the-count",
        "default-output",
        "with-bai",
        "queryname-without-bai",
        "empty",
        "empty-with-bai",
    ] {
        let parts = fields(&text, label);
        assert_eq!(parts[0], "SBI\\1", "{label}");
        assert_eq!(parts[2], "0".repeat(32), "{label} md5");
        assert_eq!(parts[3], "0".repeat(32), "{label} uuid");
        let entries: Vec<u64> = parts[7..]
            .iter()
            .map(|p| p.parse().expect("a number"))
            .collect();
        let ours = write(
            number(&parts, 1),
            number(&parts, 4),
            number(&parts, 5),
            &entries,
        );
        let theirs = corpus::decode_base64(&field(&text, "index", label).expect("an index"));
        assert_eq!(ours, theirs, "{label}");
        assert_eq!(&ours[..4], SBI_MAGIC, "{label}");
        // The count field is the number of entries and not the number of records.
        assert_eq!(number(&parts, 6), entries.len() as u64, "{label}");
    }
}

/// The granularity counts records: five of them leave three offsets and a final entry at two, five
/// and a final entry at one, and one and a final entry above the count.
#[test]
fn the_granularity_counts_records() {
    let text = golden();
    let two = fields(&text, "granularity-two");
    let one = fields(&text, "granularity-one");
    let many = fields(&text, "granularity-above-the-count");
    let default = fields(&text, "granularity-default");
    assert_eq!((number(&two, 4), number(&two, 5)), (5, 2));
    assert_eq!(number(&two, 6), 4);
    assert_eq!(number(&one, 6), 6);
    assert_eq!(number(&many, 6), 2);
    // The default is htsjdk's, which is above this file's record count.
    assert_eq!(number(&default, 5), DEFAULT_GRANULARITY);
    assert_eq!(number(&default, 6), 2);
    // Which is what the port makes of the same offsets: the ones the granularity-one run holds,
    // whose last entry is the same next-start the other runs end on.
    let every: Vec<u64> = one[7..one.len() - 1]
        .iter()
        .map(|p| p.parse().expect("a number"))
        .collect();
    let next_start: u64 = one.last().expect("a last entry").parse().expect("a number");
    let expected: Vec<u64> = two[7..]
        .iter()
        .map(|p| p.parse().expect("a number"))
        .collect();
    assert_eq!(every.len(), 5);
    assert_eq!(offsets(&every, 2, next_start), expected);
    assert_eq!(offsets(&every, DEFAULT_GRANULARITY, next_start).len(), 2);
}

/// The last entry is where the next record would have gone, which is inside the last block.
#[test]
fn the_last_entry_is_where_the_next_record_would_have_gone() {
    let text = golden();
    let sorted = fields(&text, "granularity-two");
    let file_length = number(&sorted, 1);
    let last: u64 = sorted
        .last()
        .expect("a last entry")
        .parse()
        .expect("a number");
    // Not the file's length, and not past it either: the block address is below the file length.
    assert_ne!(last, make_file_pointer(file_length));
    assert!(last >> 16 < file_length, "{last} against {file_length}");
}

/// An empty BAM has no last record to ask, so the entry falls back to the file's length.
#[test]
fn an_empty_bam_falls_back_to_the_file_length() {
    let text = golden();
    let empty = fields(&text, "empty");
    let file_length = number(&empty, 1);
    assert_eq!(number(&empty, 4), 0, "no records");
    assert_eq!(number(&empty, 6), 1, "the final entry alone");
    let last: u64 = empty
        .last()
        .expect("a last entry")
        .parse()
        .expect("a number");
    assert_eq!(last, make_file_pointer(file_length));
    // And the index is still written: sixty-eight bytes of header and the one entry.
    assert_eq!(
        corpus::decode_base64(&field(&text, "index", "empty").expect("an index")).len(),
        76
    );
}

/// The default output appends, where the `.bai` companion replaces.
#[test]
fn the_default_output_appends_and_the_companion_replaces() {
    let text = golden();
    assert_eq!(
        field(&text, "wrote", "default-output").as_deref(),
        Some("sorted.bam.sbi")
    );
    assert_eq!(default_output("sorted.bam"), "sorted.bam.sbi");
    assert_ne!(default_output("sorted.bam"), "sorted.sbi");
    assert_eq!(bai_companion("sorted.bam.sbi"), "sorted.bam.bai");
    assert_eq!(bai_companion("with-bai.sbi"), "with-bai.bai");
    // The companion is written only when it is asked for.
    assert_ne!(field(&text, "bai", "with-bai").as_deref(), Some("absent"));
    assert_eq!(
        field(&text, "bai", "granularity-two").as_deref(),
        Some("absent")
    );
}

/// Only the `.bai` path reads the records, so only it refuses a file that is not sorted.
#[test]
fn only_the_bai_path_cares_about_the_sort_order() {
    let text = golden();
    let refusal = field(&text, "error", "queryname-with-bai").expect("a refusal");
    assert!(refusal.contains("UserException$BadInput"), "{refusal}");
    assert!(refusal.contains("Cannot create a .bai index for a file that isn't coordinate sorted."));
    // Without it the same file is indexed, and to offsets of its own.
    let queryname = fields(&text, "queryname-without-bai");
    assert_eq!(number(&queryname, 4), 5);
    assert_eq!(number(&queryname, 6), 4);
    assert!(field(&text, "error", "queryname-without-bai").is_none());
}

/// A granularity of nought or less is refused before anything is opened.
#[test]
fn a_granularity_of_nought_is_refused() {
    let text = golden();
    for label in ["granularity-zero", "granularity-negative"] {
        let refusal = field(&text, "error", label).unwrap_or_else(|| panic!("{label}"));
        assert!(
            refusal.contains("CommandLineException$BadArgumentValue"),
            "{refusal}"
        );
        assert!(refusal.contains(GRANULARITY_MESSAGE), "{refusal}");
        assert!(
            field(&text, "index", label).is_none(),
            "{label} wrote nothing"
        );
    }
    assert_eq!(assert_granularity(0), Err(GRANULARITY_MESSAGE.to_string()));
    assert_eq!(assert_granularity(-1), Err(GRANULARITY_MESSAGE.to_string()));
    assert_eq!(assert_granularity(1), Ok(()));
}

/// A file that is not a BAM is refused by its extension, which the message names.
#[test]
fn a_file_that_is_not_a_bam_is_refused_by_its_extension() {
    let text = golden();
    let refusal = field(&text, "error", "not-a-bam").expect("a refusal");
    assert!(refusal.contains(
        "A splitting index is only relevant for a bam file, but a file with extension sam was \
         specified."
    ));
    assert_eq!(
        assert_is_bam("plain.sam"),
        Err(
            "A splitting index is only relevant for a bam file, but a file with extension sam was \
         specified."
                .to_string()
        )
    );
    assert_eq!(assert_is_bam("sorted.bam"), Ok(()));
}
