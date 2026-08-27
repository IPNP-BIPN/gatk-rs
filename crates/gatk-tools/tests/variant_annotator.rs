//! Conformance for `VariantAnnotator` against GATK 4.6.2.0, compared as the annotations of every
//! record of every run.
//!
//! Golden from `tools/readfilter-conformance/VariantAnnotatorDump.java`.
//!
//! The read-based annotations are not measured or ported: this suite is about the expression
//! machinery, which copies from one VCF to another.
//!
//! # What this suite is for
//!
//!  * **the key being `<tag>.<field>`**, so one file under two tags is two annotations;
//!  * **a per-allele field not crossing to a different alternate**, whatever the arguments say;
//!  * **a scalar field crossing unless allele concordance is asked for**;
//!  * **`--comp` adding a bare flag**;
//!  * **an unknown field being silent and an unknown tag being refused**;
//!  * **and every annotation the input carried surviving.**

use gatk_corpus as corpus;
use gatk_tools::variant_annotator::{
    annotate, check_expressions, comparison_annotation, value_of, Annotation, AnnotatorError,
    Arity, Expression, Record,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/variant_annotator.txt.gz"),
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

/// One VCF the golden carries, read as its records.
fn records(text: &str, name: &str) -> Vec<Record> {
    section(text, "vcf", name)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Record {
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                id: columns[2].to_string(),
                reference: columns[3].to_string(),
                alternates: columns[4].split(',').map(str::to_string).collect(),
                filters: columns[6].split(';').map(str::to_string).collect(),
                attributes: if columns[7] == "." {
                    Vec::new()
                } else {
                    columns[7]
                        .split(';')
                        .filter_map(|part| part.split_once('='))
                        .map(|(k, v)| (k.to_string(), v.to_string()))
                        .collect()
                },
            }
        })
        .collect()
}

/// The INFO column one run wrote, per position, as the writer wrote it.
fn measured(text: &str, label: &str) -> Vec<(i32, String)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse().expect("a position"),
                columns[7].to_string(),
            )
        })
        .collect()
}

/// `AC` and `AF` are Number=A; everything else in this fixture is scalar.
fn arity(field: &str) -> Arity {
    match field {
        "AC" | "AF" => Arity::PerAllele,
        _ => Arity::Scalar,
    }
}

/// The INFO column the port produces for one input record, rendered the way the writer renders it:
/// the record's own attributes and the new annotations, all sorted by key.
fn rendered(input: &Record, annotations: &[Annotation]) -> String {
    let mut parts: Vec<String> = input
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    for annotation in annotations {
        parts.push(match &annotation.value {
            Some(value) => format!("{}={value}", annotation.key),
            None => annotation.key.clone(),
        });
    }
    if parts.is_empty() {
        return ".".to_string();
    }
    parts.sort();
    parts.join(";")
}

/// label, expressions, allele concordance.
fn runs() -> Vec<(&'static str, Vec<&'static str>, bool)> {
    vec![
        ("one-expression", vec!["res.AC"], false),
        ("scalar-expression", vec!["res.NOTE"], false),
        ("scalar-concordance", vec!["res.NOTE"], true),
        ("two-expressions", vec!["res.AC", "res.AF"], false),
        (
            "id-alt-filter",
            vec!["res.ID", "res.ALT", "res.FILTER"],
            false,
        ),
        ("allele-concordance", vec!["res.AC"], true),
        ("unknown-field", vec!["res.MISSING"], false),
        ("no-resource", vec![], false),
    ]
}

#[test]
fn every_annotation_matches_the_golden() {
    let text = golden();
    let inputs = records(&text, "input");
    let resource = records(&text, "resource");
    let tagged: Vec<(String, &Record)> = resource
        .iter()
        .map(|record| ("res".to_string(), record))
        .collect();
    let mut compared = 0;
    for (label, expressions, concordance) in runs() {
        let parsed: Vec<Expression> = expressions
            .iter()
            .map(|text| Expression::parse(text).expect("an expression"))
            .collect();
        let produced: Vec<(i32, String)> = inputs
            .iter()
            .map(|input| {
                let annotations = annotate(input, &tagged, &parsed, arity, concordance);
                (input.position, rendered(input, &annotations))
            })
            .collect();
        assert_eq!(produced, measured(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 8, "the runs the port reproduces");
}

/// The key is the tag and the field, so the same file under two tags is two annotations.
#[test]
fn the_key_is_the_tag_and_the_field() {
    let text = golden();
    let inputs = records(&text, "input");
    let resource = records(&text, "resource");
    let two_tags: Vec<(String, &Record)> = resource
        .iter()
        .flat_map(|record| [("res".to_string(), record), ("other".to_string(), record)])
        .collect();
    let expressions = [
        Expression::parse("res.AC").expect("an expression"),
        Expression::parse("other.AC").expect("an expression"),
    ];
    let produced: Vec<(i32, String)> = inputs
        .iter()
        .map(|input| {
            let annotations = annotate(input, &two_tags, &expressions, arity, false);
            (input.position, rendered(input, &annotations))
        })
        .collect();
    assert_eq!(produced, measured(&text, "two-tags"));
    assert_eq!(expressions[0].key(), "res.AC");
    assert_eq!(expressions[1].key(), "other.AC");
}

/// It is withheld whatever the arguments say; a scalar one is withheld only when asked.
#[test]
fn a_per_allele_field_cannot_cross_to_a_different_alternate() {
    let text = golden();
    let at = |label: &str, position: i32| {
        measured(&text, label)
            .into_iter()
            .find(|(at, _)| *at == position)
            .expect("a record")
            .1
    };
    // Position 2000 is where the resource's alternate differs from the input's.
    assert_eq!(at("one-expression", 2000), ".", "AC never crosses");
    assert_eq!(
        at("allele-concordance", 2000),
        ".",
        "and the flag changes nothing"
    );
    assert_eq!(at("one-expression", 2000), at("allele-concordance", 2000));
    // The scalar one crosses by default and is withheld when concordance is asked for.
    assert_eq!(at("scalar-expression", 2000), "res.NOTE=second");
    assert_eq!(at("scalar-concordance", 2000), ".");
    // Where the alternates agree, both cross under either setting.
    assert!(at("one-expression", 3000).contains("res.AC=9"));
    assert!(at("scalar-concordance", 3000).contains("res.NOTE=third"));
    // ID, ALT and FILTER cross the discordant site too.
    assert!(at("id-alt-filter", 2000).contains("res.ID=rs2"));
}

/// Read off the record rather than out of its INFO column.
#[test]
fn three_fields_are_not_info_attributes() {
    let text = golden();
    let resource = records(&text, "resource");
    let first = &resource[0];
    assert_eq!(value_of(first, "ID").as_deref(), Some("rs1"));
    assert_eq!(value_of(first, "ALT").as_deref(), Some("C"));
    assert_eq!(value_of(first, "FILTER").as_deref(), Some("PASS"));
    assert_eq!(value_of(first, "AC").as_deref(), Some("5"));
    // A field the resource does not carry yields nothing, which is not a refusal.
    assert_eq!(value_of(first, "MISSING"), None);
    for (_, info) in measured(&text, "unknown-field") {
        assert!(!info.contains("res.MISSING"), "{info}");
    }
    // And the run that asked for it produced exactly what the run with no resource did.
    assert_eq!(
        measured(&text, "unknown-field"),
        measured(&text, "no-resource")
    );
}

/// A bare key, not a key and a value.
#[test]
fn a_comparison_adds_a_bare_flag() {
    let text = golden();
    let inputs = records(&text, "input");
    let resource = records(&text, "resource");
    let produced: Vec<(i32, String)> = inputs
        .iter()
        .map(|input| {
            let comparison = resource
                .iter()
                .find(|record| record.position == input.position);
            let annotations: Vec<Annotation> = comparison_annotation(input, "cmp", comparison)
                .into_iter()
                .collect();
            (input.position, rendered(input, &annotations))
        })
        .collect();
    assert_eq!(produced, measured(&text, "comparison"));
    // The flag has no value at all.
    let flag = comparison_annotation(&inputs[0], "cmp", Some(&resource[0])).expect("a flag");
    assert_eq!(flag.key, "cmp");
    assert_eq!(flag.value, None);
    // And it is subject to the same allele test as a scalar expression.
    assert!(comparison_annotation(&inputs[1], "cmp", Some(&resource[1])).is_none());
}

/// An unknown field is silent; an unknown tag is not.
#[test]
fn an_unknown_tag_is_refused() {
    let text = golden();
    let (class, message) = refusal(&text, "unknown-tag");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    let expressions = [Expression::parse("nothere.AC").expect("an expression")];
    let produced = check_expressions(&expressions, &["res".to_string()]).expect_err("no such tag");
    assert_eq!(
        produced,
        AnnotatorError::UnknownResource {
            expression: "nothere.AC".to_string()
        }
    );
    assert_eq!(produced.message(), message);
    // A known tag with an unknown FIELD passes this check and stays silent later.
    assert!(check_expressions(
        &[Expression::parse("res.MISSING").expect("an expression")],
        &["res".to_string()]
    )
    .is_ok());
}

/// The input's own annotations are added to, never replaced.
#[test]
fn the_inputs_annotations_survive() {
    let text = golden();
    // The first site carries NOTE=kept, which every run keeps.
    for (label, ..) in runs() {
        let first = measured(&text, label)
            .into_iter()
            .find(|(at, _)| *at == 1000)
            .expect("the first record")
            .1;
        assert!(first.contains("NOTE=kept"), "{label}: {first}");
    }
    // A site the resource does not mention is left alone entirely.
    for (label, ..) in runs() {
        let last = measured(&text, label)
            .into_iter()
            .find(|(at, _)| *at == 4000)
            .expect("the last record")
            .1;
        assert_eq!(last, ".", "{label}");
    }
}
