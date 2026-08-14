//! Conformance for `CountFalsePositives` against GATK 4.6.2.0, compared as the whole output table
//! of every run and as the class and message of every refusal.
//!
//! Golden from `tools/readfilter-conformance/CountFalsePositivesDump.java`.
//!
//! # What this suite is for
//!
//!  * **everything that is not an indel is a SNP**, so seven passing records of seven shapes are
//!    counted 5 and 2;
//!  * **the id is the file name with one extension removed**, not a sample;
//!  * **the territory is the merged intervals in bases**;
//!  * **and the rates are per megabase, divided before they are scaled**, and keep their `.0`.
//!
//! # What is compared, and what is not
//!
//! The `no-intervals` row is asserted as a message rather than replayed: `requiresIntervals()`
//! makes `-L` a required argument, so that refusal is Barclay's and belongs to the argument suite.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_tools::count_false_positives::{
    count, id_from_path, table, target_territory, CountFalsePositivesError,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/count_false_positives.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The records of one input, decoded as far as the two branches look at them.
fn records(text: &str, label: &str) -> Vec<VariantContext> {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    );
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let mut alleles = vec![Allele::create(field[3].as_bytes(), true).expect("a reference")];
            if field[4] != "." {
                for alternate in field[4].split(',') {
                    alleles
                        .push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
                }
            }
            let mut variant =
                VariantContext::new(field[0], field[1].parse().expect("a position"), alleles);
            variant.stop = variant.start + field[3].len() as i64 - 1;
            variant.filters = Some(match field[6] {
                "PASS" | "." => Vec::new(),
                names => names.split(';').map(|name| name.to_string()).collect(),
            });
            variant
        })
        .collect()
}

fn interval(start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new("chr1", start, end).expect("a valid interval")
}

/// Every run that produced a table: which file it read, which records the intervals select, and
/// the merged intervals themselves.
struct Run {
    file: &'static str,
    id: &'static str,
    /// The merged `-L`, as `IntervalArgumentCollection` produces it.
    intervals: Vec<SimpleInterval>,
    /// The records the traversal hands to `apply`, as positions in the file.
    selected: Option<Vec<usize>>,
}

fn setup(label: &str) -> Run {
    let whole = |file: &'static str, id: &'static str| Run {
        file,
        id,
        intervals: vec![interval(1, 1000)],
        selected: None,
    };
    match label {
        "whole-contig" => whole("mixed", "mixed"),
        "all-filtered" => whole("all-filtered", "all-filtered"),
        "two-extensions" => whole("two-extensions.vcf", "two-extensions.vcf"),
        "overlapping-intervals" => Run {
            intervals: vec![interval(1, 400)],
            ..whole("one-snp", "one-snp")
        },
        "disjoint-intervals" => Run {
            intervals: vec![interval(1, 100), interval(900, 1000)],
            ..whole("one-snp", "one-snp")
        },
        "small-territory" => Run {
            intervals: vec![interval(98, 100)],
            ..whole("one-snp", "one-snp")
        },
        "integral-rate" => Run {
            intervals: vec![interval(1, 100)],
            ..whole("one-snp", "one-snp")
        },
        // `-L chr1:1-350` reaches the SNP at 100, the insertion at 200 and the deletion at 300.
        "interval-selects-some" => Run {
            intervals: vec![interval(1, 350)],
            selected: Some(vec![0, 1, 2]),
            ..whole("mixed", "mixed")
        },
        other => panic!("no setup for {other}"),
    }
}

#[test]
fn every_table_matches_the_golden_byte_for_byte() {
    let text = golden();
    let tables = rows(&text, "table");
    assert_eq!(tables.len(), 8, "eight of the ten runs write a table");

    for row in tables {
        let (label, expected) = (row[0], unescape(row[1]));
        let run = setup(label);
        let all = records(&text, run.file);
        let traversed: Vec<VariantContext> = match &run.selected {
            None => all,
            Some(indices) => indices.iter().map(|at| all[*at].clone()).collect(),
        };
        let counts = count(&traversed);
        let mine = table(run.id, counts, target_territory(&run.intervals));
        assert_eq!(mine, expected, "the table of {label}");
    }
}

#[test]
fn seven_passing_records_of_seven_shapes_are_counted_five_and_two() {
    let text = golden();
    let counts = count(&records(&text, "mixed"));
    assert_eq!(
        counts.snp, 5,
        "the MNP, the symbolic, the mixed and the no-variation are all here"
    );
    assert_eq!(counts.indel, 2);
}

#[test]
fn the_id_keeps_the_first_of_two_extensions() {
    let text = golden();
    let row = rows(&text, "table")
        .into_iter()
        .find(|row| row[0] == "two-extensions")
        .expect("the run over a name with two extensions");
    let written = unescape(row[1]);
    let id = written
        .lines()
        .nth(1)
        .expect("the record")
        .split('\t')
        .next()
        .expect("the id");
    assert_eq!(id, "two-extensions.vcf");
    assert_eq!(
        id_from_path("countfalsepositives-dump/two-extensions.vcf.vcf"),
        id
    );
}

#[test]
fn the_output_refusal_carries_the_references_class_and_words() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "output-is-a-directory")
        .expect("the output refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");

    let refused = CountFalsePositivesError::CouldNotOpenOutput {
        path: "countfalsepositives-dump/.".to_string(),
    };
    assert_eq!(refused.class(), class);
    assert_eq!(refused.message(), message);
}

/// The other refusal is the argument parser's, and it is asserted rather than replayed.
#[test]
fn the_missing_interval_is_refused_by_barclay() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "no-intervals")
        .expect("the missing -L");
    assert_eq!(
        row[1],
        "org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument:\
         Argument intervals was missing: Argument 'intervals' is required"
    );
}
