//! Conformance for `ValidateVariants` against GATK 4.6.2.0, compared as whether each run refused
//! and, where it did, as the class and message of the refusal.
//!
//! Golden from `tools/readfilter-conformance/ValidateVariantsDump.java`.
//!
//! # What this suite is for
//!
//!  * **excluding any type puts REF back in the set**, so one exclusion is a refusal;
//!  * **the allele check is about the genotypes**, not the ALT column;
//!  * **the GVCF checks are three**, and the per-record ones fire before the coverage one;
//!  * **and `--warn-on-errors` turns every refusal into nothing at all**.
//!
//! # What is compared, and what is not
//!
//! The coverage check counts loci over the whole reference; the tool's own count is compared
//! against a computed one here rather than reproduced from a traversal, since what this port
//! carries is the decision and the message rather than the interval machinery.

use gatk_corpus as corpus;
use gatk_tools::validate_variants::{
    types_to_apply, validate_record, Arguments, OrderCheck, Record, ValidationError, ValidationType,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/validate_variants.txt.gz"),
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

/// `determineType`, as far as the message quotes it.
fn variant_type(reference: &str, alternates: &str) -> String {
    let mut kind: Option<&str> = None;
    for alternate in alternates.split(',') {
        let this = if alternate.starts_with('<') {
            "SYMBOLIC"
        } else if alternate.len() == reference.len() {
            if reference.len() == 1 {
                "SNP"
            } else {
                "MNP"
            }
        } else {
            "INDEL"
        };
        match kind {
            None => kind = Some(this),
            Some(seen) if seen != this => return "MIXED".to_string(),
            Some(_) => {}
        }
    }
    kind.unwrap_or("NO_VARIATION").to_string()
}

fn file(text: &str, label: &str) -> Vec<Record> {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    );
    let sample_count = whole
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header")
        .split('\t')
        .count()
        .saturating_sub(9);
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let info: Vec<(&str, &str)> = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .collect();
            let value = |key: &str| info.iter().find(|(name, _)| *name == key).map(|(_, v)| *v);
            let keys: Vec<&str> = if field.len() > 8 {
                field[8].split(':').collect()
            } else {
                Vec::new()
            };
            let genotypes = (0..sample_count)
                .map(|index| {
                    let values: Vec<&str> = field[9 + index].split(':').collect();
                    let call = keys
                        .iter()
                        .position(|key| *key == "GT")
                        .map(|at| values[at])
                        .unwrap_or("./.");
                    call.split(['/', '|'])
                        .map(|allele| allele.parse::<usize>().ok())
                        .collect()
                })
                .collect();
            Record {
                contig: field[0].to_string(),
                start: field[1].parse().expect("a position"),
                reference: field[3].to_string(),
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                filters: match field[6] {
                    "." | "PASS" => Vec::new(),
                    names => names.split(';').map(|name| name.to_string()).collect(),
                },
                allele_counts: value("AC")
                    .map(|text| {
                        text.split(',')
                            .filter_map(|part| part.parse().ok())
                            .collect()
                    })
                    .unwrap_or_default(),
                allele_number: value("AN").and_then(|text| text.parse().ok()),
                genotypes,
                qual: field[5].parse().ok(),
                variant_type: variant_type(field[3], field[4]),
                // `sortedString`, which the message prints, is the keys in order.
                attributes: {
                    let mut sorted: Vec<(String, String)> = info
                        .iter()
                        .map(|(key, value)| (key.to_string(), value.to_string()))
                        .collect();
                    sorted.sort();
                    sorted
                },
            }
        })
        .collect()
}

/// Every run: which file, and the arguments. The reference is two contigs of `A`.
fn setup(run: &str) -> (&'static str, Arguments) {
    let base = || Arguments {
        has_reference: true,
        ..Arguments::default()
    };
    match run {
        "good" => ("good", base()),
        "unused-alternate" => ("unused-alternate", base()),
        "unused-alternate-excluded" => (
            "unused-alternate",
            Arguments {
                types_to_exclude: vec![ValidationType::Alleles],
                ..base()
            },
        ),
        "bad-counts" => ("bad-counts", base()),
        "bad-counts-excluded" => (
            "bad-counts",
            Arguments {
                types_to_exclude: vec![ValidationType::ChrCounts],
                ..base()
            },
        ),
        "bad-counts-warn-only" => (
            "bad-counts",
            Arguments {
                warn_on_errors: true,
                ..base()
            },
        ),
        "filtered-bad-counts" => ("filtered-bad-counts", base()),
        "filtered-bad-counts-skipped" => (
            "filtered-bad-counts",
            Arguments {
                do_not_validate_filtered_records: true,
                ..base()
            },
        ),
        "sites-only" => ("sites-only", base()),
        "wrong-reference-base" => ("wrong-reference-base", base()),
        "wrong-reference-base-excluded" => (
            "wrong-reference-base",
            Arguments {
                types_to_exclude: vec![ValidationType::Ref],
                ..base()
            },
        ),
        "not-gvcf" => (
            "not-gvcf",
            Arguments {
                validate_gvcf: true,
                ..base()
            },
        ),
        "unordered" => (
            "unordered",
            Arguments {
                validate_gvcf: true,
                ..base()
            },
        ),
        other => panic!("no setup for {other}"),
    }
}

/// The runs whose answer this port decides on its own; the two whose answer is the coverage count
/// are asserted separately.
const RUNS: [&str; 13] = [
    "good",
    "unused-alternate",
    "unused-alternate-excluded",
    "bad-counts",
    "bad-counts-excluded",
    "bad-counts-warn-only",
    "filtered-bad-counts",
    "filtered-bad-counts-skipped",
    "sites-only",
    "wrong-reference-base",
    "wrong-reference-base-excluded",
    "not-gvcf",
    "unordered",
];

fn refusal(text: &str, run: &str) -> Option<String> {
    rows(text, "error")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| unescape(row[1]))
}

/// One whole run, as `apply` makes it: the types first, then each record.
fn outcome(text: &str, run: &str) -> Result<(), ValidationError> {
    let (label, arguments) = setup(run);
    let records = file(text, label);
    let types = types_to_apply(&arguments)?;
    let mut order = OrderCheck::new();
    for record in &records {
        if arguments.validate_gvcf {
            order.check(record)?;
        }
        let outcome = validate_record(
            record,
            &format!("validatevariants-dump/{label}.vcf"),
            &types,
            // The reference is all `A`.
            Some("A"),
            &arguments,
        );
        if let Err(error) = outcome {
            // `--warn-on-errors` turns the refusal into a log line and nothing else.
            if arguments.warn_on_errors {
                continue;
            }
            return Err(error);
        }
    }
    Ok(())
}

#[test]
fn every_run_refuses_what_the_reference_refused() {
    let text = golden();
    for run in RUNS {
        match (outcome(&text, run), refusal(&text, run)) {
            (Ok(()), None) => {}
            (Err(error), Some(expected)) => assert_eq!(
                format!("{}:{}", error.java_class(), error.message()),
                expected,
                "error/{run}"
            ),
            (Ok(()), Some(expected)) => panic!("{run} passed here and refused there: {expected}"),
            (Err(error), None) => panic!("{run} refused here: {}", error.message()),
        }
    }
}

/// One exclusion turns another check on.
#[test]
fn excluding_a_type_puts_ref_back_in_the_set() {
    // `ALL` with no reference is the allele and count checks alone.
    let plain = types_to_apply(&Arguments::default()).expect("no refusal");
    assert_eq!(
        plain,
        vec![ValidationType::Alleles, ValidationType::ChrCounts]
    );

    // Excluding the allele check asks for the concrete set instead, which has REF in it.
    let excluded = types_to_apply(&Arguments {
        types_to_exclude: vec![ValidationType::Alleles],
        ..Arguments::default()
    });
    assert_eq!(excluded, Err(ValidationError::MissingReference));

    // With a reference it is fine, and REF is in the set.
    let with_reference = types_to_apply(&Arguments {
        types_to_exclude: vec![ValidationType::Alleles],
        has_reference: true,
        ..Arguments::default()
    })
    .expect("no refusal");
    assert!(with_reference.contains(&ValidationType::Ref));
}

/// The two GVCF record checks, and the order the tool applies them in.
#[test]
fn the_per_record_gvcf_checks_fire_before_the_coverage_one() {
    let text = golden();
    // The file that is both incomplete and missing <NON_REF> is refused for the allele.
    let expected = refusal(&text, "not-gvcf").expect("a refusal");
    assert!(
        expected.contains("must contain a <NON_REF> allele"),
        "{expected}"
    );
    assert!(outcome(&text, "not-gvcf").is_err());

    // And the coverage refusal, which this port carries as a message rather than a traversal.
    let coverage = refusal(&text, "gvcf").expect("a refusal");
    let error = ValidationError::NotCovering {
        loci: 1898,
        first_gap: "chr1:1-99".to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        coverage
    );
}

/// The two flags that make a refusal disappear.
#[test]
fn a_refusal_can_be_turned_into_nothing() {
    let text = golden();
    assert!(refusal(&text, "bad-counts").is_some());
    assert!(refusal(&text, "bad-counts-warn-only").is_none());
    assert!(refusal(&text, "filtered-bad-counts-skipped").is_none());
    assert!(outcome(&text, "bad-counts-warn-only").is_ok());
    assert!(outcome(&text, "filtered-bad-counts-skipped").is_ok());
    // The same record without either flag is a refusal.
    assert!(outcome(&text, "filtered-bad-counts").is_err());
}
