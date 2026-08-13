//! Conformance for `SelectVariants`' sample selection against GATK 4.6.2.0, compared as the sample
//! columns of every output header, the class and message of every refusal, and the expression
//! matcher on its own.
//!
//! Golden from `tools/readfilter-conformance/SampleSelectionDump.java`.
//!
//! # What this suite is for
//!
//!  * **the expressions are unanchored `find()`**, so `s1` selects `xs10`;
//!  * **matching nothing selects everything**, by two different routes;
//!  * **the output order is sorted**, not the command line's;
//!  * **and exclusion beats inclusion**, emptying the set into a refusal.

use gatk_corpus as corpus;
use gatk_engine::java_regex::{filter_collection_by_expressions, Pattern};
use gatk_tools::select_variants::{
    create_sample_name_inclusion_list, SampleArguments, SampleSelectionError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_samples.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
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

fn input(text: &str) -> String {
    unescape(rows(text, "input").first().expect("an input")[1])
}

/// The samples as the file declares them, which is not their sorted order.
fn header_samples(text: &str) -> Vec<String> {
    input(text)
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header line")
        .split('\t')
        .skip(9)
        .map(|name| name.to_string())
        .collect()
}

/// The records of the input.
fn input_records(text: &str) -> Vec<String> {
    input(text)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.to_string())
        .collect()
}

/// The INFO column of a record, which is where a rewritten record gains AC, AF and AN.
fn info(record: &str) -> String {
    record
        .split('\t')
        .nth(7)
        .expect("an INFO column")
        .to_string()
}

/// A record's genotypes, keyed by the sample the given column order names.
fn genotypes(record: &str, samples: &[String]) -> Vec<(String, String)> {
    let mut pairs: Vec<(String, String)> = samples
        .iter()
        .cloned()
        .zip(record.split('\t').skip(9).map(|call| call.to_string()))
        .collect();
    pairs.sort();
    pairs
}

fn value(text: &str, kind: &str, label: &str) -> Option<String> {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .map(|row| unescape(row[1]))
}

fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// How the golden holds a refusal: the class, a colon, and whatever prefix the class adds to
/// `getMessage`. `UserException$BadInput` prepends "Bad input: ", which is the class's doing.
fn rendered(error: &SampleSelectionError) -> String {
    let prefix = match error {
        SampleSelectionError::SamplesNotInHeader(_) => "Bad input: ",
        _ => "",
    };
    format!("{}:{}{}", error.java_class(), prefix, error.message())
}

fn names(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

/// The arguments of each run, which the golden does not carry.
fn setup(run: &str) -> SampleArguments {
    let base = SampleArguments::default;
    match run {
        "all-samples" => base(),
        "one-name" => SampleArguments {
            sample_names: names(&["s1"]),
            ..base()
        },
        "two-names-reversed" => SampleArguments {
            sample_names: names(&["tumor", "s0"]),
            ..base()
        },
        "expression-substring" => SampleArguments {
            sample_expressions: names(&["s1"]),
            ..base()
        },
        "expression-anchored" => SampleArguments {
            sample_expressions: names(&["^s1$"]),
            ..base()
        },
        "expression-prefix" => SampleArguments {
            sample_expressions: names(&["^NA"]),
            ..base()
        },
        "expression-matches-nothing" => SampleArguments {
            sample_expressions: names(&["zzz"]),
            ..base()
        },
        "name-and-expression" => SampleArguments {
            sample_names: names(&["tumor"]),
            sample_expressions: names(&["^NA"]),
            ..base()
        },
        "missing-name" => SampleArguments {
            sample_names: names(&["ghost"]),
            ..base()
        },
        "missing-name-allowed" => SampleArguments {
            sample_names: names(&["ghost"]),
            allow_nonoverlapping_command_line_samples: true,
            ..base()
        },
        "missing-and-present-allowed" => SampleArguments {
            sample_names: names(&["ghost", "s1"]),
            allow_nonoverlapping_command_line_samples: true,
            ..base()
        },
        "two-missing-names" => SampleArguments {
            sample_names: names(&["zeta", "alpha"]),
            ..base()
        },
        "exclude-one" => SampleArguments {
            exclude_sample_names: names(&["s0"]),
            ..base()
        },
        "exclude-expression" => SampleArguments {
            exclude_sample_expressions: names(&["^s"]),
            ..base()
        },
        "exclude-what-was-included" => SampleArguments {
            sample_names: names(&["s0"]),
            exclude_sample_names: names(&["s0"]),
            ..base()
        },
        "exclude-some-of-what-was-included" => SampleArguments {
            sample_names: names(&["s0", "s1"]),
            exclude_sample_names: names(&["s0"]),
            ..base()
        },
        "exclude-everything" => SampleArguments {
            exclude_sample_names: names(&["tumor", "s1", "NA12891", "xs10", "s0", "NA12878"]),
            ..base()
        },
        "exclude-missing-name" => SampleArguments {
            exclude_sample_names: names(&["ghost"]),
            ..base()
        },
        "uncompilable-expression" => SampleArguments {
            sample_expressions: names(&["["]),
            ..base()
        },
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 19] = [
    "all-samples",
    "one-name",
    "two-names-reversed",
    "expression-substring",
    "expression-anchored",
    "expression-prefix",
    "expression-matches-nothing",
    "name-and-expression",
    "missing-name",
    "missing-name-allowed",
    "missing-and-present-allowed",
    "two-missing-names",
    "exclude-one",
    "exclude-expression",
    "exclude-what-was-included",
    "exclude-some-of-what-was-included",
    "exclude-everything",
    "exclude-missing-name",
    "uncompilable-expression",
];

/// The expression matcher on its own, over the samples in the order the file declares them.
const EXPRESSIONS: [(&str, &[&str]); 7] = [
    ("substring", &["s1"]),
    ("anchored", &["^s1$"]),
    ("prefix", &["^NA"]),
    ("nothing", &["zzz"]),
    ("exact-name", &["NA12878"]),
    ("two", &["^s", "tumor"]),
    ("everything", &["."]),
];

#[test]
fn every_expression_matches_what_the_reference_matched() {
    let text = golden();
    let samples = header_samples(&text);
    for (label, patterns) in EXPRESSIONS {
        let expressions = names(patterns);
        let matched = filter_collection_by_expressions(&samples, &expressions, false)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
        let expected = value(&text, "expressions", label).expect("a row");
        assert_eq!(matched.join(","), expected, "expressions/{label}");
    }
}

#[test]
fn every_selection_is_the_reference_s() {
    let text = golden();
    let samples = header_samples(&text);
    for run in RUNS {
        let result = create_sample_name_inclusion_list(&samples, &setup(run));
        match value(&text, "error", run) {
            Some(expected) => {
                let error = result.expect_err(run);
                assert_eq!(rendered(&error), expected, "error/{run}");
            }
            None => {
                let selection = result.unwrap_or_else(|error| panic!("{run}: {}", error.message()));
                let expected = value(&text, "samples", run).expect("a sample row");
                assert_eq!(selection.samples.join(","), expected, "samples/{run}");
            }
        }
    }
}

/// The two ways an empty selection becomes the whole cohort, neither of them an empty output.
#[test]
fn matching_nothing_selects_everything() {
    let text = golden();
    let samples = header_samples(&text);
    let mut sorted = samples.clone();
    sorted.sort();

    for run in ["expression-matches-nothing", "missing-name-allowed"] {
        let selection = create_sample_name_inclusion_list(&samples, &setup(run)).expect(run);
        assert_eq!(selection.samples, sorted, "samples/{run}");
        assert!(selection.no_samples_specified, "flag/{run}");
        assert_eq!(
            value(&text, "samples", run).expect("a row"),
            sorted.join(",")
        );
        // And no record was rewritten: same INFO, without the AC, AF and AN a subset gains, and
        // the same call for every sample.
        let ours = written(&text, run);
        let theirs = input_records(&text);
        assert_eq!(ours.len(), theirs.len(), "records/{run}");
        for (out, input) in ours.iter().zip(theirs.iter()) {
            assert_eq!(info(out), info(input), "info/{run}");
            assert_eq!(
                genotypes(out, &sorted),
                genotypes(input, &samples),
                "genotypes/{run}"
            );
        }
    }
}

/// The columns are the header's order and the header's order is sorted, so every run permutes them
/// even when it selects everything and rewrites nothing.
#[test]
fn the_columns_come_out_sorted_even_when_no_record_is_rewritten() {
    let text = golden();
    let samples = header_samples(&text);
    let mut sorted = samples.clone();
    sorted.sort();
    assert_ne!(
        samples, sorted,
        "the input is written out of order on purpose"
    );

    let selection =
        create_sample_name_inclusion_list(&samples, &setup("all-samples")).expect("all");
    assert!(selection.no_samples_specified);
    assert_eq!(selection.samples, sorted);

    let first = written(&text, "all-samples")[0].clone();
    let input_first = input_records(&text)[0].clone();
    // Same values, different columns: the calls follow the sample names rather than their places.
    assert_eq!(
        first.split('\t').skip(9).collect::<Vec<_>>(),
        vec!["0/1:60", "0/1:30", "0/0:50", "0/0:60", "0/1:20", "1/1:40"]
    );
    assert_eq!(
        genotypes(&first, &sorted),
        genotypes(&input_first, &samples)
    );
}

/// The unanchored search, which is the difference between `find()` and `matches()`.
#[test]
fn an_expression_reaches_further_than_it_reads() {
    let text = golden();
    let samples = header_samples(&text);
    assert!(samples.contains(&"xs10".to_string()));

    let selection =
        create_sample_name_inclusion_list(&samples, &setup("expression-substring")).expect("kept");
    assert_eq!(
        selection.samples,
        vec!["s1".to_string(), "xs10".to_string()]
    );
    assert_eq!(
        value(&text, "samples", "expression-substring").expect("a row"),
        "s1,xs10"
    );

    // The same pattern anchored keeps one, which is what the argument was meant to say.
    let anchored =
        create_sample_name_inclusion_list(&samples, &setup("expression-anchored")).expect("kept");
    assert_eq!(anchored.samples, vec!["s1".to_string()]);
}

/// The command line's order is not the output's.
#[test]
fn the_output_order_is_sorted() {
    let text = golden();
    let samples = header_samples(&text);
    let selection =
        create_sample_name_inclusion_list(&samples, &setup("two-names-reversed")).expect("kept");
    assert_eq!(
        selection.samples,
        vec!["s0".to_string(), "tumor".to_string()]
    );
    assert_eq!(
        value(&text, "samples", "two-names-reversed").expect("a row"),
        "s0,tumor"
    );
}

/// Both routes to an emptied set are the same refusal, and the missing-name one keeps its order.
#[test]
fn the_refusals_are_the_reference_s() {
    let text = golden();
    let samples = header_samples(&text);

    for run in ["exclude-what-was-included", "exclude-everything"] {
        let error = create_sample_name_inclusion_list(&samples, &setup(run)).unwrap_err();
        assert_eq!(error, SampleSelectionError::AllExcluded, "{run}");
        assert_eq!(rendered(&error), value(&text, "error", run).expect("a row"));
    }

    // The names are listed in the order they were given, while everything else here is sorted.
    let error =
        create_sample_name_inclusion_list(&samples, &setup("two-missing-names")).unwrap_err();
    assert!(error.message().contains("\n\nzeta,alpha\n\n"));
    assert_eq!(
        rendered(&error),
        value(&text, "error", "two-missing-names").expect("a row")
    );
}

/// The regex engine's own refusal, reaching the user unwrapped.
#[test]
fn an_uncompilable_expression_is_the_pattern_compilers_refusal() {
    let text = golden();
    let expected = value(&text, "error", "uncompilable-expression").expect("a row");
    assert!(expected.starts_with("java.util.regex.PatternSyntaxException:"));

    let error = Pattern::compile("[").expect_err("refused");
    assert_eq!(error.description, "Unclosed character class");
    assert_eq!(error.index, 0);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        expected
    );
}
