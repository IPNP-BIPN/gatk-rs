//! Conformance for the reference walker and `CountBasesInReference`, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/ReferenceWalkerDump.java`.
//!
//! # What this suite is for
//!
//!  * **an absent `-L` is the whole reference**, 67 loci over the two-contig fixture, which is the
//!    branch `IntervalWalker` cannot reach because its `requiresIntervals()` is true;
//!  * **`-XL` alone is legal here** and subtracts from that whole-reference list;
//!  * **the intervals are traversed in dictionary order**, not in the order they were given;
//!  * **every locus is one base**, so the apply count is the number of bases and not of intervals;
//!  * **the bases are the caching reader's**, upper-cased with every IUPAC code flattened to `N`,
//!    which is what makes the count table five rows over a fifteen-symbol FASTA;
//!  * **and `getReferenceWindow` is a method**: the widened walker is still called once per base
//!    and sees five, clipped at the contig start by its own `max(1, ...)`.
//!
//! The FASTA and its `.fai` are rebuilt here rather than carried: the dump writes
//! `ReferenceQueryDump.FASTA` to a temporary directory and indexes it, and the same bytes are
//! written here. The `apply` rows are what say the two agree.
//!
//! The golden also carries the tool's own standard output, because `onTraversalSuccess` prints the
//! table twice: to the file `-O` names and to stdout. Those bare rows are asserted beside the
//! `counts` row rather than skipped, since printing to both is the behaviour.

use std::io::Write;

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::count_bases_in_reference;
use gatk_tools::reference_walker::{self, TraversalError};

/// `ReferenceQueryDump.FASTA`, byte for byte.
const FASTA: &str = ">chr1 first contig\n\
                    ACGTACGTACGT\n\
                    acgtNNNNacgt\n\
                    ACGTRYKMSWBD\n\
                    HVNACGT\n\
                    >chr2\n\
                    TTTTGGGGCCCC\n\
                    AAAATTTTGGGG\n";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/reference_walker.txt.gz"),
    )
}

/// The FASTA and a `.fai` for it, written once per test into a temporary directory.
fn reference() -> ReferenceFileSource {
    let dir = std::env::temp_dir().join(format!(
        "reference-walker-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).expect("a temporary directory");
    let fasta = dir.join("ref.fasta");
    let mut file = std::fs::File::create(&fasta).expect("the fixture");
    file.write_all(FASTA.as_bytes()).expect("the fixture");
    drop(file);

    // The `.fai` the dump has htsjdk build: name, length, offset of the first base, bases per line
    // and bytes per line, one row per contig.
    let mut fai = std::fs::File::create(dir.join("ref.fasta.fai")).expect("the index");
    writeln!(fai, "chr1\t43\t19\t12\t13").expect("the index");
    writeln!(fai, "chr2\t24\t72\t12\t13").expect("the index");
    drop(fai);

    ReferenceFileSource::open(&fasta).expect("the reference opens")
}

fn arguments(include: &[&str], exclude: &[&str], padding: i32) -> IntervalArguments {
    IntervalArguments {
        include: include.iter().map(|s| s.to_string()).collect(),
        exclude: exclude.iter().map(|s| s.to_string()).collect(),
        padding,
        ..Default::default()
    }
}

/// The dump's traversals, in its order: label, arguments, and the window the walker uses.
type Case = (&'static str, IntervalArguments, bool);

fn cases() -> Vec<Case> {
    vec![
        ("all", arguments(&[], &[], 0), false),
        ("chr1-window", arguments(&["chr1:5-15"], &[], 0), false),
        ("chr2", arguments(&["chr2"], &[], 0), false),
        (
            "out-of-order",
            arguments(&["chr2:3-5", "chr1:1-3"], &[], 0),
            false,
        ),
        ("excluded-only", arguments(&[], &["chr1:5-40"], 0), false),
        ("padded", arguments(&["chr1:10-12"], &[], 2), false),
        (
            "abutting",
            arguments(&["chr1:1-3", "chr1:4-6"], &[], 0),
            false,
        ),
        ("past-the-end", arguments(&["chr2:20-40"], &[], 0), false),
        ("wide-window", arguments(&["chr1:1-4"], &[], 0), true),
    ]
}

/// The dump's `WideWalker.getReferenceWindow`: two bases either side, clipped at the contig start.
fn wide(locus: &SimpleInterval) -> SimpleInterval {
    SimpleInterval::new(&locus.contig, 1.max(locus.start - 2), locus.end + 2)
        .expect("the widened window is valid")
}

fn rows<'a>(text: &'a str, kind: &str, label: &str) -> Vec<&'a str> {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .collect()
}

#[test]
fn every_traversal_matches_the_golden() {
    let text = golden();
    let mut total = 0;

    for (label, args, widened) in cases() {
        let mut source = reference();
        let traversed = if widened {
            reference_walker::traverse(&mut source, &args, wide)
        } else {
            reference_walker::traverse(&mut source, &args, |locus: &SimpleInterval| locus.clone())
        };

        let expected_count: usize = rows(&text, "count", label)[0].parse().expect("a number");
        let expected_summary = rows(&text, "summary", label)[0];

        let applied = match traversed {
            Ok(applied) => {
                assert_eq!(expected_summary, "ok", "{label}: the reference refused");
                applied
            }
            Err(error) => {
                // The one measured refusal: an interval past the end of its contig is a malformed
                // genome location rather than a clip.
                assert_eq!(
                    expected_summary,
                    "E:org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc",
                    "{label}: {error:?}"
                );
                assert_eq!(expected_count, 0, "{label}: a refusal applies nothing");
                total += 1;
                continue;
            }
        };

        assert_eq!(applied.len(), expected_count, "{label}: the apply count");
        let expected = rows(&text, "apply", label);
        assert_eq!(expected.len(), applied.len(), "{label}: the apply rows");
        for (index, call) in applied.iter().enumerate() {
            let mine = format!(
                "{index}\t{}:{}-{}|{}",
                call.window.contig,
                call.window.start,
                call.window.end,
                String::from_utf8_lossy(&call.bases)
            );
            assert_eq!(mine, expected[index], "{label} apply {index}");
        }
        total += 1;
    }

    assert_eq!(total, 9, "every traversal in the golden");
}

#[test]
fn every_count_table_matches_the_golden() {
    let text = golden();

    for (label, args) in [
        ("all", arguments(&[], &[], 0)),
        ("chr1-window", arguments(&["chr1:5-15"], &[], 0)),
        ("iupac", arguments(&["chr1:25-36"], &[], 0)),
    ] {
        let mut source = reference();
        let counts = count_bases_in_reference::run(&mut source, &args).expect("the traversal runs");
        let expected = rows(&text, "counts", label)[0];
        // The dump escapes the file's text; the report is the same string unescaped.
        let mine = counts.report().replace('\n', "\\n");
        assert_eq!(mine, expected, "{label}: the count table");

        // And the same table on standard output, which the golden carries as bare rows.
        for line in counts.report().lines() {
            assert!(
                text.lines().any(|row| row == line),
                "{label}: the golden carries {line} on standard output"
            );
        }
    }
}

/// The refusal is the walker's, not the reference file's: nothing is read before it fires.
#[test]
fn an_interval_past_the_end_is_refused_before_any_base_is_read() {
    let mut source = reference();
    let error = reference_walker::traverse(
        &mut source,
        &arguments(&["chr2:20-40"], &[], 0),
        |locus: &SimpleInterval| locus.clone(),
    )
    .expect_err("past the end of chr2");
    assert!(
        matches!(error, TraversalError::Intervals(_)),
        "the interval parser refuses it, {error:?}"
    );
}
