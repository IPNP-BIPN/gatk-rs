//! Conformance for `FlowPairHMMAlignReadsToHaplotypes` against GATK 4.6.2.0, compared as the whole
//! output file of every run.
//!
//! Golden from `tools/readfilter-conformance/FlowPairHMMAlignReadsToHaplotypesDump.java`.
//!
//! The alignment engine is not measured or ported: the scores are read off each run's own expanded
//! matrix and the concise file beside it is rebuilt from them.
//!
//! # What this suite is for
//!
//!  * **the expanded format being one column per haplotype in the FASTA's order**;
//!  * **the concise format's five columns**;
//!  * **the reference score being recorded only while the reference is the best so far**;
//!  * **which makes `Diff_from_ref` depend on the FASTA's order**;
//!  * **no reference haplotype, and one the FASTA does not name, being the same thing**;
//!  * **a read that matches nothing having an empty name, `-Infinity` and two `NaN`s**;
//!  * **the two engines disagreeing about that read alone**;
//!  * **and an unknown engine being a bare RuntimeException.**

use gatk_corpus as corpus;
use gatk_tools::flow_pairhmm_align_reads_to_haplotypes::{
    buffers, concise, concise_file, concise_row, expanded_file, expanded_header, is_known_engine,
    three_decimals, Haplotype, BUFFER_SIZE_LIMIT, CONCISE_HEADER, ENGINES, UNKNOWN_ENGINE_MESSAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/flow_pairhmm_align_reads_to_haplotypes.txt.gz"),
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

/// The haplotype names one of the golden's FASTAs declares, in its own order.
fn fasta_names(text: &str, name: &str) -> Vec<String> {
    section(text, "fasta", name)
        .lines()
        .filter(|line| line.starts_with('>'))
        .map(|line| line[1..].to_string())
        .collect()
}

fn haplotypes(text: &str, fasta: &str, reference: Option<&str>) -> Vec<Haplotype> {
    fasta_names(text, fasta)
        .iter()
        .map(|name| Haplotype::new(name, reference))
        .collect()
}

/// One run's expanded matrix, as its header and its rows.
fn expanded(text: &str, label: &str) -> (Vec<String>, Vec<(String, Vec<f64>)>) {
    let file = section(text, "out", label);
    let mut lines = file.lines();
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .split('\t')
        .skip(1)
        .map(str::to_string)
        .collect();
    let rows = lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[0].to_string(),
                columns[1..]
                    .iter()
                    .map(|value| value.parse().expect("a score"))
                    .collect(),
            )
        })
        .collect();
    (header, rows)
}

/// Every expanded file is what the port writes from its own numbers.
#[test]
fn every_expanded_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, fasta) in [
        ("expanded", "haplotypes"),
        ("expanded-ref", "haplotypes"),
        ("hmm-expanded", "haplotypes"),
        ("expanded-ref-last", "reordered"),
        ("expanded-tied", "tied"),
    ] {
        let haplotypes = haplotypes(&text, fasta, None);
        let (header, rows) = expanded(&text, label);
        assert_eq!(
            header,
            haplotypes
                .iter()
                .map(|h| h.name.clone())
                .collect::<Vec<_>>(),
            "{label}"
        );
        assert_eq!(
            expanded_file(&haplotypes, &rows),
            section(&text, "out", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the expanded runs the port reproduces");
}

/// The concise file's two difference columns, compared against a rebuild from the ROUNDED scores.
///
/// The expanded file is the only place the scores are available, and it carries them to three
/// decimals. A difference taken between two rounded numbers is not the rounding of the difference:
/// -0.210 less -6.199 is 5.989 where the true scores give 5.990. The name and the best score are
/// therefore compared exactly and the two differences to within one unit of the last decimal.
fn assert_concise_matches(produced: &str, written: &str, label: &str) {
    let rows = |file: &str| -> Vec<Vec<String>> {
        file.lines()
            .filter(|line| !line.is_empty())
            .map(|line| line.split('\t').map(str::to_string).collect())
            .collect()
    };
    let produced = rows(produced);
    let written = rows(written);
    assert_eq!(produced.len(), written.len(), "{label}");
    for (produced, written) in produced.iter().zip(written.iter()) {
        assert_eq!(produced[..3], written[..3], "{label}");
        for column in 3..5 {
            let (produced, written) = (&produced[column], &written[column]);
            // Rust parses `NaN` and `Infinity` as numbers, so the two words are compared as
            // words before any arithmetic is tried on them.
            let numeric = |value: &str| {
                value
                    .parse::<f64>()
                    .ok()
                    .filter(|value: &f64| value.is_finite())
            };
            match (numeric(produced), numeric(written)) {
                (Some(produced), Some(written)) => assert!(
                    (produced - written).abs() <= 0.0011,
                    "{label}: {produced} vs {written}"
                ),
                _ => assert_eq!(produced, written, "{label}"),
            }
        }
    }
}

/// Every concise file is what the port writes from the expanded matrix beside it.
#[test]
fn every_concise_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, from, fasta, reference) in [
        ("concise-ref", "expanded", "haplotypes", Some("hap_ref")),
        ("concise-no-ref", "expanded", "haplotypes", None),
        (
            "concise-unknown-ref",
            "expanded",
            "haplotypes",
            Some("nothere"),
        ),
        ("hmm-concise", "hmm-expanded", "haplotypes", Some("hap_ref")),
        (
            "concise-ref-last",
            "expanded-ref-last",
            "reordered",
            Some("hap_ref"),
        ),
    ] {
        let haplotypes = haplotypes(&text, fasta, reference);
        let (_, rows) = expanded(&text, from);
        assert_concise_matches(
            &concise_file(&haplotypes, &rows),
            &section(&text, "out", label),
            label,
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the concise runs the port reproduces");
}

/// The five columns, in that order.
#[test]
fn the_concise_format_is_five_columns() {
    let text = golden();
    let file = section(&text, "out", "concise-ref");
    assert_eq!(file.lines().next().expect("a header"), CONCISE_HEADER);
    assert_eq!(CONCISE_HEADER.split('\t').count(), 5);
    for line in file.lines().skip(1).filter(|line| !line.is_empty()) {
        assert_eq!(line.split('\t').count(), 5, "{line}");
    }
    // The expanded one is `Read` and the FASTA's own names.
    let haplotypes = haplotypes(&text, "haplotypes", None);
    assert_eq!(
        expanded_header(&haplotypes),
        section(&text, "out", "expanded")
            .lines()
            .next()
            .expect("a header")
    );
}

/// A reference that comes after a better haplotype is never recorded at all.
#[test]
fn the_reference_score_depends_on_the_fastas_order() {
    let text = golden();
    // The same haplotypes, the same reads, the same scores: only the FASTA's order differs.
    let (_, first) = expanded(&text, "expanded");
    let (_, last) = expanded(&text, "expanded-ref-last");
    let names = |rows: &[(String, Vec<f64>)]| -> Vec<String> {
        rows.iter().map(|(name, _)| name.clone()).collect()
    };
    assert_eq!(names(&first), names(&last));
    // With the reference first every read reports a real difference.
    let with_first = section(&text, "out", "concise-ref");
    let reference_column = |file: &str| -> Vec<String> {
        file.lines()
            .skip(1)
            .filter(|line| !line.is_empty())
            .map(|line| line.split('\t').nth(4).expect("a column").to_string())
            .collect()
    };
    assert_eq!(
        reference_column(&with_first),
        vec!["0.000", "5.990", "11.979", "17.969", "NaN"]
    );
    // With the reference last every read whose best is not the reference reports an infinity.
    let with_last = section(&text, "out", "concise-ref-last");
    assert_eq!(
        reference_column(&with_last),
        vec!["0.000", "Infinity", "Infinity", "Infinity", "NaN"]
    );
    // Which the port reproduces from the same numbers under the two orders.
    let ordered = haplotypes(&text, "haplotypes", Some("hap_ref"));
    let reordered = haplotypes(&text, "reordered", Some("hap_ref"));
    let scores = &first[1].1;
    assert!(concise(&ordered, scores).reference_score.is_finite());
    let shuffled = &last[1].1;
    assert_eq!(
        concise(&reordered, shuffled).reference_score,
        f64::NEG_INFINITY
    );
}

/// Both leave every haplotype non-reference.
#[test]
fn an_unnamed_reference_is_the_same_as_none() {
    let text = golden();
    assert_eq!(
        section(&text, "out", "concise-no-ref"),
        section(&text, "out", "concise-unknown-ref")
    );
    let none = haplotypes(&text, "haplotypes", None);
    let unknown = haplotypes(&text, "haplotypes", Some("nothere"));
    assert_eq!(none, unknown);
    assert!(none.iter().all(|haplotype| !haplotype.is_reference));
    // A name the FASTA does carry marks exactly one.
    let named = haplotypes(&text, "haplotypes", Some("hap_ref"));
    assert_eq!(named.iter().filter(|h| h.is_reference).count(), 1);
    assert_eq!(named[0].name, "hap_ref");
}

/// An empty name, a score of `-Infinity` and two `NaN`s.
#[test]
fn a_read_that_matches_nothing_has_no_best_haplotype() {
    let text = golden();
    let line = section(&text, "out", "concise-ref")
        .lines()
        .find(|line| line.starts_with("r-like-none"))
        .expect("its line")
        .to_string();
    assert_eq!(line, "r-like-none\t\t-Infinity\tNaN\tNaN");
    // Which the port produces from the row of negative infinities the expanded file carries.
    let (_, rows) = expanded(&text, "expanded");
    let scores = &rows
        .iter()
        .find(|(name, _)| name == "r-like-none")
        .expect("its row")
        .1;
    assert!(scores.iter().all(|score| *score == f64::NEG_INFINITY));
    let reduced = concise(&haplotypes(&text, "haplotypes", Some("hap_ref")), scores);
    assert_eq!(reduced.best_haplotype, "");
    assert_eq!(reduced.best_score, f64::NEG_INFINITY);
    assert!(reduced.difference_from_second().is_nan());
    assert!(reduced.difference_from_reference().is_nan());
    assert_eq!(concise_row("r-like-none", &reduced), line);
    // An infinity prints by name and a number to three decimals.
    assert_eq!(three_decimals(f64::NEG_INFINITY), "-Infinity");
    assert_eq!(three_decimals(f64::INFINITY), "Infinity");
    assert_eq!(three_decimals(f64::NAN), "NaN");
    assert_eq!(three_decimals(-0.2104), "-0.210");
}

/// FlowBased leaves it unmatched where FlowBasedHMM gives it the reference.
#[test]
fn the_two_engines_disagree_about_that_read_alone() {
    let text = golden();
    let best = |label: &str| -> Vec<String> {
        section(&text, "out", label)
            .lines()
            .skip(1)
            .filter(|line| !line.is_empty())
            .map(|line| line.split('\t').nth(1).expect("a column").to_string())
            .collect()
    };
    assert_eq!(
        best("concise-ref"),
        vec!["hap_ref", "hap_one", "hap_two", "hap_three", ""]
    );
    assert_eq!(
        best("hmm-concise"),
        vec!["hap_ref", "hap_one", "hap_two", "hap_three", "hap_ref"]
    );
    // The four reads that match a haplotype agree; only the fifth differs.
    assert_eq!(best("concise-ref")[..4], best("hmm-concise")[..4]);
    assert_ne!(best("concise-ref")[4], best("hmm-concise")[4]);
    // The numbers differ throughout, though.
    let (_, flow_based) = expanded(&text, "expanded");
    let (_, hmm) = expanded(&text, "hmm-expanded");
    assert_ne!(flow_based[0].1, hmm[0].1);
}

/// A printed difference of 0.000 is a rounding and not an equality.
#[test]
fn two_scores_that_print_the_same_need_not_be_tied() {
    let text = golden();
    let (_, rows) = expanded(&text, "expanded-tied");
    // The third read's two best columns print the same.
    let scores = &rows[2].1;
    assert_eq!(three_decimals(scores[1]), three_decimals(scores[2]));
    // And the concise line reports a difference of zero.
    let file = section(&text, "out", "concise-tied");
    let line = file
        .lines()
        .find(|line| line.starts_with("r-like-two"))
        .expect("its line");
    assert!(line.contains("\t0.000\t"), "{line}");
    // The best haplotype is the SECOND of the two, so the two scores were not tied at all: the
    // rounded numbers cannot say which won, which is why this run is left out of the rebuild.
    assert_eq!(line.split('\t').nth(1), Some("hap_right"));
    let rebuilt = concise(&haplotypes(&text, "tied", Some("hap_ref")), scores);
    assert_eq!(
        rebuilt.best_haplotype, "hap_left",
        "what the rounded numbers alone would say"
    );
    // The comparison is strict, so an EXACT tie keeps the first.
    let two = vec![
        Haplotype::new("first", None),
        Haplotype::new("second", None),
    ];
    let tied = concise(&two, &[-1.0, -1.0]);
    assert_eq!(tied.best_haplotype, "first");
    assert_eq!(tied.difference_from_second(), 0.0);
    // A later score that beats the best displaces it and pushes the old best down.
    let rising = concise(&two, &[-2.0, -1.0]);
    assert_eq!(rising.best_haplotype, "second");
    assert_eq!(rising.second_best_score, -2.0);
}

/// The concise file's differences are not the differences of the expanded file's numbers.
///
/// Both files round to three decimals, but the concise one rounds the difference where a reader
/// of the expanded one can only difference the roundings. The two disagree by a unit in the last
/// place often enough that the fixture shows it on three of its five rows.
#[test]
fn the_two_files_cannot_be_derived_from_each_other() {
    let text = golden();
    let (_, rows) = expanded(&text, "expanded");
    let haplotypes = haplotypes(&text, "haplotypes", Some("hap_ref"));
    let from_rounded = concise_file(&haplotypes, &rows);
    let written = section(&text, "out", "concise-ref");
    // They are not the same file.
    assert_ne!(from_rounded, written);
    // The first row differs in its second-best column: 5.989 against 5.990.
    let column = |file: &str, row: usize, at: usize| -> String {
        file.lines()
            .nth(row)
            .expect("a row")
            .split('\t')
            .nth(at)
            .expect("a column")
            .to_string()
    };
    assert_eq!(column(&from_rounded, 1, 3), "5.989");
    assert_eq!(column(&written, 1, 3), "5.990");
    // And the difference is exactly one unit in the last place, everywhere it appears.
    for row in 1..5 {
        for at in 3..5 {
            let produced = column(&from_rounded, row, at);
            let expected = column(&written, row, at);
            let numeric = |value: &str| {
                value
                    .parse::<f64>()
                    .ok()
                    .filter(|value: &f64| value.is_finite())
            };
            match (numeric(&produced), numeric(&expected)) {
                (Some(a), Some(b)) => assert!((a - b).abs() <= 0.0011, "{produced} vs {expected}"),
                _ => assert_eq!(produced, expected),
            }
        }
    }
    // The columns that are not differences agree exactly.
    for row in 1..6 {
        assert_eq!(column(&from_rounded, row, 1), column(&written, row, 1));
        assert_eq!(column(&from_rounded, row, 2), column(&written, row, 2));
    }
}

/// The engine name is checked against the two the tool has.
#[test]
fn an_unknown_engine_is_a_bare_runtime_exception() {
    let text = golden();
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix("error\tunknown-engine\t"))
        .expect("the golden carries error/unknown-engine");
    let (class, message) = row.split_once(':').expect("a class and a message");
    assert_eq!(class, "java.lang.RuntimeException");
    assert_eq!(message, UNKNOWN_ENGINE_MESSAGE);
    assert!(!is_known_engine("PairHMM"));
    for engine in ENGINES {
        assert!(is_known_engine(engine));
    }
    // The default is the first of the two, which is the engine every unlabelled run used.
    assert_eq!(ENGINES[0], "FlowBased");
}

/// The buffer decides when the matrix is computed, not what it holds.
#[test]
fn the_reads_are_scored_fifty_at_a_time() {
    let text = golden();
    assert_eq!(BUFFER_SIZE_LIMIT, 50);
    // The fixture has five reads, so one buffer, and its file has one header and five rows.
    let (_, rows) = expanded(&text, "expanded");
    assert_eq!(rows.len(), 5);
    assert_eq!(buffers(5), vec![5]);
    assert_eq!(buffers(0), vec![0]);
    // A count over the limit is split, and the remainder is flushed at the end.
    assert_eq!(buffers(120), vec![50, 50, 20]);
    // A count that is a multiple of the limit still gets a final flush, of nothing.
    assert_eq!(buffers(100), vec![50, 50, 0]);
}
