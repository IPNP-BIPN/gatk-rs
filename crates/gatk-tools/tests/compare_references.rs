//! Conformance for `CompareReferences` against GATK 4.6.2.0, compared as the whole table of every
//! run and as the analysis it printed.
//!
//! Golden from `tools/readfilter-conformance/CompareReferencesDump.java`, which carries every
//! fasta and every dictionary as well as both halves of each run's output.
//!
//! # What this suite is for
//!
//!  * **the table**, keyed by MD5, with `---` where a reference has no such sequence;
//!  * **the pair analysis**, including the pair that keeps both `DIFFER_IN_SEQUENCE_NAMES` and
//!    `DIFFER_IN_SEQUENCES_PRESENT` because the superset rule does not apply to it;
//!  * **the status order**, which is the enum's and not the order the flags were added;
//!  * **and the three MD5 modes**, one of which refuses and one of which trusts a dictionary that
//!    lies.
//!
//! The sequences' MD5s are read from the golden's own dictionaries rather than recomputed, which
//! is what makes the modes comparable here: `ALWAYS_RECALCULATE` is given the true digest and
//! `USE_DICT` the one the file carries.

use gatk_corpus as corpus;
use gatk_tools::compare_references::{
    build, compare_all, write_table, Md5Mode, Reference, Sequence, TableError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/compare_references.txt.gz"),
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

/// One reference, read out of the golden's own dictionary. `truth` is the dictionary whose M5 is
/// the real digest, which is what recalculating would find.
fn reference(text: &str, label: &str, dict: &str, truth: &str) -> Reference {
    let real: Vec<(String, String)> = sequences(&value(text, "dict", truth))
        .into_iter()
        .map(|(name, _, md5)| (name, md5.unwrap_or_default()))
        .collect();
    let sequences = sequences(&value(text, "dict", dict))
        .into_iter()
        .map(|(name, length, md5)| {
            let calculated = real
                .iter()
                .find(|(other, _)| *other == name)
                .map(|(_, digest)| digest.clone())
                .unwrap_or_default();
            Sequence {
                name,
                length,
                md5,
                calculated_md5: calculated,
            }
        })
        .collect();
    Reference {
        column: format!("{label}.fasta"),
        sequences,
    }
}

/// The `@SQ` lines of a dictionary: name, length and `M5` when it has one.
fn sequences(dict: &str) -> Vec<(String, i64, Option<String>)> {
    dict.lines()
        .filter(|line| line.starts_with("@SQ\t"))
        .map(|line| {
            let field = |tag: &str| {
                line.split('\t')
                    .find_map(|part| part.strip_prefix(&format!("{tag}:")))
                    .map(str::to_string)
            };
            (
                field("SN").expect("a name"),
                field("LN").and_then(|v| v.parse().ok()).expect("a length"),
                field("M5"),
            )
        })
        .collect()
}

/// The analysis as the reference prints it, banner included.
fn printed(pairs: &[gatk_tools::compare_references::Pair]) -> String {
    let mut out = String::from("*********************************************************\n");
    for pair in pairs {
        out.push_str(&pair.rendered());
        out.push('\n');
    }
    out
}

fn run(text: &str, labels: &[(&str, &str)], mode: Md5Mode) -> (String, String) {
    let references: Vec<Reference> = labels
        .iter()
        .map(|(label, dict)| reference(text, label, dict, label))
        .collect();
    let table = build(&references, mode).expect("a table the tool builds");
    let pairs = compare_all(&table, &references).expect("an analysis the tool allows");
    (write_table(&table), printed(&pairs))
}

#[test]
fn every_table_and_analysis_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, references) in [
        ("renamed", vec![("base", "base"), ("renamed", "renamed")]),
        ("altered", vec![("base", "base"), ("altered", "altered")]),
        ("superset", vec![("extra", "extra"), ("base", "base")]),
        ("subset", vec![("fewer", "fewer"), ("base", "base")]),
        (
            "three",
            vec![("base", "base"), ("renamed", "renamed"), ("extra", "extra")],
        ),
    ] {
        let (table, analysis) = run(&text, &references, Md5Mode::RecalculateIfMissing);
        assert_eq!(table, value(&text, "table", label), "{label}: the table");
        assert_eq!(
            analysis,
            value(&text, "stdout", label),
            "{label}: the analysis"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's runs");
}

/// The pair that both renames and omits keeps both flags: the superset rule is skipped when a
/// naming discrepancy was found.
#[test]
fn a_renaming_pair_is_never_a_superset() {
    let text = golden();
    let (_, analysis) = run(
        &text,
        &[("base", "base"), ("renamed", "renamed"), ("extra", "extra")],
        Md5Mode::RecalculateIfMissing,
    );
    // The third pair, renamed against extra.
    assert!(analysis.contains(
        "REFERENCE PAIR: renamed.fasta, extra.fasta\nStatus:\n\tDIFFER_IN_SEQUENCE_NAMES\n\tDIFFER_IN_SEQUENCES_PRESENT\n"
    ));
    // While the pair that only omits is a subset outright.
    assert!(analysis.contains("REFERENCE PAIR: base.fasta, extra.fasta\nStatus:\n\tSUBSET\n"));
}

/// A reference compared with itself is not a comparison at all: the two paths collapse into one
/// map entry, so no pair is generated and the tool walks off the empty list. The golden records
/// the `IndexOutOfBoundsException`, and `compare_all` answers the empty list that produces it.
#[test]
fn a_reference_compared_with_itself_produces_no_pair() {
    let text = golden();
    let base = reference(&text, "base", "base", "base");
    let table = build(std::slice::from_ref(&base), Md5Mode::RecalculateIfMissing)
        .expect("a table over one reference");
    let pairs = compare_all(&table, std::slice::from_ref(&base)).expect("no pairs");
    assert!(pairs.is_empty());
    assert_eq!(
        refusal(&text, "identical"),
        "java.lang.IndexOutOfBoundsException:Index 0 out of bounds for length 0"
    );
}

/// `USE_DICT` refuses a dictionary with no M5, and trusts one that lies.
#[test]
fn the_md5_modes_read_different_things() {
    let text = golden();

    let stripped = reference(&text, "no-md5", "no-md5-stripped", "no-md5");
    let base = reference(&text, "base", "base", "base");
    let error = build(&[base.clone(), stripped.clone()], Md5Mode::UseDict)
        .expect_err("the missing-M5 refusal");
    assert_eq!(
        error,
        TableError::MissingMd5 {
            sequence: "chr1".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "use-dict-missing")
    );

    // A dictionary whose M5 is a lie: trusted as written, so its sequences share no row with the
    // base's and the table has three of them.
    let lying = reference(&text, "wrong-md5", "wrong-md5-wrong", "wrong-md5");
    let (table, analysis) = {
        let references = vec![base.clone(), lying.clone()];
        let table = build(&references, Md5Mode::UseDict).expect("a table");
        let pairs = compare_all(&table, &references).expect("an analysis");
        (write_table(&table), printed(&pairs))
    };
    assert_eq!(table, value(&text, "table", "use-dict-wrong"));
    assert_eq!(analysis, value(&text, "stdout", "use-dict-wrong"));

    // And recalculating ignores the lie, so the two are an exact match again.
    let references = vec![base, lying];
    let table = build(&references, Md5Mode::AlwaysRecalculate).expect("a table");
    let pairs = compare_all(&table, &references).expect("an analysis");
    assert_eq!(
        write_table(&table),
        value(&text, "table", "always-recalculate")
    );
    assert_eq!(
        printed(&pairs),
        value(&text, "stdout", "always-recalculate")
    );
}
