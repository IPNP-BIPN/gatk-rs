//! Conformance for `AlignmentUtils.leftAlignIndels` and `normalizeAlleles` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/LeftAlignIndelsDump.java`. No files here: the
//! interface is a cigar, a reference, a read and a start in, and a cigar and two counts out, so
//! the golden is a table and every row is a call.
//!
//! # What this suite is for
//!
//! `LeftAlignIndels` is eight lines around this call, and four of the twenty-one rows are things
//! that would not survive a port written from the source alone:
//!
//!  * **two indels can cancel each other.** The cigar is walked right to left and an indel's
//!    ranges are only resolved at an alignment block, so `3M2I2M2D3M` over a homopolymer comes
//!    back as `10M`;
//!  * **a deletion that reaches the start of the cigar is dropped** and the reference bases it
//!    removed are reported instead. Seven of the twelve cases end this way, which makes moving the
//!    read the ordinary outcome rather than the corner;
//!  * **`normalizeAlleles` can return a negative start shift**, because trimming shared bases off
//!    the front moves the range right;
//!  * **soft and hard clips are not alignment blocks**, so an indel may not shift into one.
//!
//! Both refusals are `IllegalArgumentException` out of a util rather than `UserException` out of a
//! tool, and the golden carries their messages.

use gatk_corpus as corpus;
use gatk_engine::alignment_utils::{self as util, IndexRange};
use htsjdk_bam::text_parse::parse_cigar;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/left_align_indels.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// `[from,to)` as the dump prints it.
fn range(text: &str) -> IndexRange {
    let inner = text
        .strip_prefix('[')
        .and_then(|text| text.strip_suffix(')'))
        .expect("a half-open range");
    let (from, to) = inner.split_once(',').expect("two bounds");
    IndexRange::new(from.parse().unwrap(), to.parse().unwrap())
}

#[test]
fn every_left_aligned_cigar_is_the_reference_s() {
    let text = golden();
    let cases = rows(&text, "leftalign");
    assert_eq!(cases.len(), 12, "the golden lost cases");

    let mut dropped_a_leading_deletion = 0;
    for row in &cases {
        let (label, cigar, reference, read, start) = (row[0], row[1], row[2], row[3], row[4]);
        let result = util::left_align_indels(
            &parse_cigar(cigar).expect("a parsable cigar"),
            reference.as_bytes(),
            read.as_bytes(),
            start.parse().expect("a read start"),
        )
        .unwrap_or_else(|error| panic!("{label}: the reference did not refuse: {error:?}"));

        assert_eq!(result.cigar.to_text(), row[5], "{label}: the cigar");
        assert_eq!(
            result.leading_deletion_bases_removed.to_string(),
            row[6],
            "{label}: leading deletion bases removed"
        );
        assert_eq!(
            result.trailing_deletion_bases_removed.to_string(),
            row[7],
            "{label}: trailing deletion bases removed"
        );
        if result.leading_deletion_bases_removed > 0 {
            dropped_a_leading_deletion += 1;
        }
    }

    // The count matters: a deletion walking off the start is what makes the tool move the read,
    // and it is the ordinary outcome rather than the corner.
    assert_eq!(
        dropped_a_leading_deletion, 7,
        "seven of the twelve cases move the read"
    );
    println!("left-align-indels: {} cigars compared", cases.len());
}

/// Two indels with too few matching bases between them are merged, and the merge can be empty.
#[test]
fn two_indels_that_meet_are_merged() {
    let text = golden();
    let by_label = |label: &str| -> Vec<String> {
        rows(&text, "leftalign")
            .into_iter()
            .find(|row| row[0] == label)
            .map(|row| row.iter().map(|field| field.to_string()).collect())
            .unwrap_or_else(|| panic!("the golden lost {label}"))
    };

    let cancelled = by_label("insertion-then-deletion");
    assert_eq!(cancelled[1], "3M2I2M2D3M", "what went in");
    assert_eq!(cancelled[5], "10M", "what came out: no indel at all");
    assert_eq!(cancelled[6], "0", "and nothing was removed from the front");

    let merged = by_label("colliding-indels");
    assert_eq!(merged[1], "3M1D2M1D3M");
    assert_eq!(merged[5], "8M");
    assert_eq!(merged[6], "2", "two deletions became two removed bases");
}

#[test]
fn the_refusals_are_the_reference_s() {
    let text = golden();
    let refusals = rows(&text, "leftalignerror");
    assert_eq!(refusals.len(), 2);

    let inputs = [
        ("past-the-reference", "4M1D3M", "AAAA", "AAAAAAT"),
        (
            "cigar-misses-read-bases",
            "4M1D3M",
            "AAAAAAAAAAAA",
            "AAAAAAAAAAAA",
        ),
    ];
    for row in &refusals {
        let (label, expected) = (row[0], row[1]);
        let (_, cigar, reference, read) = inputs
            .iter()
            .find(|(name, _, _, _)| *name == label)
            .unwrap_or_else(|| panic!("{label} is in the golden but not configured here"));
        let error = util::left_align_indels(
            &parse_cigar(cigar).expect("a parsable cigar"),
            reference.as_bytes(),
            read.as_bytes(),
            0,
        )
        .expect_err("the reference refused");
        let (class, message) = expected.split_once(':').expect("a class and a message");
        assert_eq!(class, error.class(), "{label}: the class");
        assert_eq!(message, error.message(), "{label}: the message");
    }
}

/// The inner function, with the ranges it adjusts as a side effect.
#[test]
fn normalize_alleles_moves_the_ranges_the_reference_moved() {
    let text = golden();
    let cases = rows(&text, "normalize");
    assert_eq!(cases.len(), 7, "the golden lost cases");

    let mut negative = 0;
    for row in &cases {
        let (label, reference, read) = (row[0], row[1], row[2]);
        let mut bounds = [range(row[3]), range(row[4])];
        let max_shift: i32 = row[5].parse().expect("a max shift");
        let trim: bool = row[6].parse().expect("a trim flag");

        let (start_shift, end_shift) = util::normalize_alleles(
            &[reference.as_bytes(), read.as_bytes()],
            &mut bounds,
            max_shift,
            trim,
        )
        .unwrap_or_else(|error| panic!("{label}: the reference did not refuse: {error:?}"));

        assert_eq!(start_shift.to_string(), row[7], "{label}: the start shift");
        assert_eq!(end_shift.to_string(), row[8], "{label}: the end shift");
        assert_eq!(bounds[0], range(row[9]), "{label}: the reference range");
        assert_eq!(bounds[1], range(row[10]), "{label}: the read range");
        if start_shift < 0 {
            negative += 1;
        }
    }

    // The row the signed shift exists for. A port that typed it as unsigned passes every other
    // case here.
    assert_eq!(negative, 1, "one case shifts right");
}
