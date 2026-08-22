//! Conformance for `DownsampleByDuplicateSet` against GATK 4.6.2.0, compared as the reads every
//! run writes.
//!
//! Golden from `tools/readfilter-conformance/DownsampleByDuplicateSetDump.java`.
//!
//! # What this suite is for
//!
//!  * **the last duplicate set escaping every rejection rule**, measured three ways;
//!  * **a rejected set consuming no random draw**, so a molecule at the front changes nothing
//!    downstream;
//!  * **an odd set being rejected at the defaults**;
//!  * **the fixed seed of 142 and the draw being `nextDouble() < fraction`**;
//!  * **the sets being cut on the molecule number alone**;
//!  * **and the one refusal**, a molecule number that goes backwards.

use gatk_corpus as corpus;
use gatk_tools::downsample_by_duplicate_set::{reject_set, run, Arguments, DownsampleError, Read};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/downsample_by_duplicate_set.txt.gz"),
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

/// The reads of a SAM text, header dropped, as the walker sees them: name, molecule, strand.
///
/// An unmapped read is dropped here, because the walker's own read filters drop it before any set
/// is formed, and `MAPPED` is one of them.
fn reads(sam: &str) -> Vec<Read> {
    sam.lines()
        .filter(|line| !line.starts_with('@'))
        .filter(|line| {
            let flags: i32 = line
                .split('\t')
                .nth(1)
                .expect("a flag")
                .parse()
                .expect("a flag");
            flags & 4 == 0
        })
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let tag = fields
                .iter()
                .find_map(|field| field.strip_prefix("MI:Z:"))
                .expect("an MI tag");
            let (molecule, strand) = tag.split_once('/').expect("a molecule and a strand");
            Read {
                name: fields[0].to_string(),
                molecule: molecule.parse().expect("a molecule number"),
                strand: strand.to_string(),
            }
        })
        .collect()
}

/// The read names a run wrote, which is what the two sides are compared on.
fn names(reads: &[Read]) -> Vec<String> {
    reads.iter().map(|read| read.name.clone()).collect()
}

fn produced(
    text: &str,
    label: &str,
    arguments: &Arguments,
) -> Result<Vec<String>, DownsampleError> {
    Ok(names(&run(
        &reads(&value(text, "input", label)),
        arguments,
    )?))
}

fn expected(text: &str, label: &str) -> Vec<String> {
    names(&reads(&value(text, "sam", label)))
}

#[test]
fn every_downsampled_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, arguments) in [
        ("keep-all", Arguments::keeping(1.0)),
        ("keep-none", Arguments::keeping(0.0)),
        ("keep-half", Arguments::keeping(0.5)),
        ("keep-most", Arguments::keeping(0.95)),
        ("odd-set-in-front", Arguments::keeping(0.5)),
        ("odd-set-at-the-end", Arguments::keeping(1.0)),
        ("one-read-file", Arguments::keeping(1.0)),
        (
            "min-reads-4",
            Arguments {
                minimum_reads: 4,
                ..Arguments::keeping(1.0)
            },
        ),
        (
            "min-per-strand-2",
            Arguments {
                minimum_reads_per_strand: 2,
                ..Arguments::keeping(1.0)
            },
        ),
        ("unmapped-read-inside", Arguments::keeping(1.0)),
    ] {
        assert_eq!(
            produced(&text, label, &arguments).expect("a run that is not refused"),
            expected(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 10, "the golden's outputs");
}

/// `processLastReadSet` never asks `rejectSet`, so the trailing set is written whatever it is.
#[test]
fn the_last_set_escapes_every_rejection_rule() {
    let text = golden();

    // Ten two-read molecules under a minimum of four: only the last survives.
    let kept = produced(
        &text,
        "min-reads-4",
        &Arguments {
            minimum_reads: 4,
            ..Arguments::keeping(1.0)
        },
    )
    .expect("a run");
    assert_eq!(kept.len(), 2);
    assert_eq!(kept, expected(&text, "min-reads-4"));

    // The same under a per-strand minimum of two.
    let per_strand = produced(
        &text,
        "min-per-strand-2",
        &Arguments {
            minimum_reads_per_strand: 2,
            ..Arguments::keeping(1.0)
        },
    )
    .expect("a run");
    assert_eq!(per_strand.len(), 2);

    // A file of one read: odd, under the minimum, and written all the same.
    let alone = produced(&text, "one-read-file", &Arguments::keeping(1.0)).expect("a run");
    assert_eq!(alone.len(), 1);

    // And a three-read molecule at the end, which the rules would otherwise reject.
    let trailing = produced(&text, "odd-set-at-the-end", &Arguments::keeping(1.0)).expect("a run");
    assert_eq!(trailing.len(), 23);
}

/// The rejection happens before the draw, so a rejectable molecule at the front changes nothing
/// downstream: the two runs keep exactly the same molecules.
#[test]
fn a_rejected_set_consumes_no_draw() {
    let text = golden();
    let plain = produced(&text, "keep-half", &Arguments::keeping(0.5)).expect("a run");
    let with_odd = produced(&text, "odd-set-in-front", &Arguments::keeping(0.5)).expect("a run");
    // The reads are named by position in the file, so compare the molecules rather than the names.
    let molecules = |label: &str, kept: &[String]| {
        let all = reads(&value(&text, "input", label));
        kept.iter()
            .map(|name| {
                all.iter()
                    .find(|read| &read.name == name)
                    .expect("a read")
                    .molecule
            })
            .collect::<Vec<i32>>()
    };
    assert_eq!(
        molecules("keep-half", &plain),
        molecules("odd-set-in-front", &with_odd)
    );
    assert_eq!(
        molecules("keep-half", &plain),
        vec![1, 1, 2, 2, 3, 3, 6, 6, 7, 7, 9, 9]
    );
}

/// An odd set is rejected at the defaults, whatever the minimums say.
#[test]
fn an_odd_set_is_rejected_at_the_defaults() {
    let three: Vec<Read> = ["A", "A", "B"]
        .iter()
        .enumerate()
        .map(|(index, strand)| Read {
            name: format!("r{index}"),
            molecule: 0,
            strand: strand.to_string(),
        })
        .collect();
    assert!(reject_set(&three, &Arguments::keeping(1.0)));
    // Two reads pass every default.
    assert!(!reject_set(&three[..2], &Arguments::keeping(1.0)));
}

/// A fraction of one keeps everything and a fraction of zero keeps nothing, since `nextDouble()`
/// answers in `[0, 1)`.
#[test]
fn the_two_extreme_fractions_are_total() {
    let text = golden();
    assert_eq!(
        produced(&text, "keep-all", &Arguments::keeping(1.0))
            .expect("a run")
            .len(),
        20
    );
    assert!(produced(&text, "keep-none", &Arguments::keeping(0.0))
        .expect("a run")
        .is_empty());
}

#[test]
fn the_one_refusal_matches_the_golden() {
    let text = golden();
    let error = run(
        &reads(&value(&text, "input", "unsorted-molecule-ids")),
        &Arguments::keeping(1.0),
    )
    .expect_err("a refused run");
    assert_eq!(error, DownsampleError::NotSortedByMoleculeId);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unsorted-molecule-ids")
    );
}
