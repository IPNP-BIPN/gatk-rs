//! Conformance for `SelectVariants`' `--discordance` and `--concordance` against GATK 4.6.2.0,
//! compared as which records each run wrote.
//!
//! Golden from `tools/readfilter-conformance/SelectVariantsConcordanceDump.java`.
//!
//! # What this suite is for
//!
//!  * **without a sample, both flags are about presence alone**;
//!  * **with one, they compare allele sets**, so `0/1` matches `1/0` and `1/1` matches `1`;
//!  * **a filtered genotype never matches anything**, another filtered one included;
//!  * **and `--exclude-filtered` moves no record**, which is what makes that clause dead code.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{Allele, Variant};
use gatk_tools::select_variants::{
    create_sample_name_inclusion_list, is_concordant, is_discordant, ComparisonArguments, Record,
    SampleArguments,
};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_concordance.txt.gz"),
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

/// One of the two input files, by its label.
fn file(text: &str, label: &str) -> (Vec<String>, Vec<Record>) {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    );
    let samples: Vec<String> = whole
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header")
        .split('\t')
        .skip(9)
        .map(|name| name.to_string())
        .collect();
    let records = whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let mut alleles = vec![Allele::new(field[3].as_bytes(), true)];
            for alternate in field[4].split(',') {
                alleles.push(Allele::new(alternate.as_bytes(), false));
            }
            let keys: Vec<&str> = field[8].split(':').collect();
            let genotypes = (0..samples.len())
                .map(|index| {
                    let values: Vec<&str> = field[9 + index].split(':').collect();
                    let by_key: BTreeMap<&str, &str> =
                        keys.iter().copied().zip(values.iter().copied()).collect();
                    let call = by_key.get("GT").copied().unwrap_or("./.");
                    Genotype {
                        alleles: call
                            .split(['/', '|'])
                            .map(|allele| allele.parse::<usize>().ok())
                            .collect(),
                        pl: None,
                        gq: by_key.get("GQ").and_then(|gq| gq.parse().ok()),
                        ad: None,
                        dp: None,
                        attributes: by_key
                            .get("FT")
                            .filter(|ft| **ft != "." && **ft != "PASS")
                            .map(|ft| vec![("FT".to_string(), ft.to_string())])
                            .unwrap_or_default(),
                    }
                })
                .collect();
            Record {
                variant: Variant {
                    contig: field[0].to_string(),
                    start: field[1].parse().expect("a position"),
                    stop: field[1].parse().expect("a position"),
                    alleles,
                    genotypes,
                    attributes: Vec::new(),
                },
                samples: samples.clone(),
            }
        })
        .collect();
    (samples, records)
}

fn names(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

fn setup(run: &str) -> (SampleArguments, ComparisonArguments) {
    let none = SampleArguments::default;
    let one = || SampleArguments {
        sample_names: names(&["s0"]),
        ..none()
    };
    let both = || SampleArguments {
        sample_names: names(&["s0", "s1"]),
        ..none()
    };
    let discordance = ComparisonArguments {
        discordance_only: true,
        ..ComparisonArguments::default()
    };
    let concordance = ComparisonArguments {
        concordance_only: true,
        ..ComparisonArguments::default()
    };
    match run {
        "no-comparison" => (none(), ComparisonArguments::default()),
        "discordance" => (none(), discordance),
        "concordance" => (none(), concordance),
        "discordance-one-sample" => (one(), discordance),
        "concordance-one-sample" => (one(), concordance),
        "discordance-both-samples" => (both(), discordance),
        "concordance-both-samples" => (both(), concordance),
        "discordance-exclude-filtered" => (
            one(),
            ComparisonArguments {
                exclude_filtered: true,
                ..discordance
            },
        ),
        "concordance-exclude-filtered" => (
            one(),
            ComparisonArguments {
                exclude_filtered: true,
                ..concordance
            },
        ),
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 9] = [
    "no-comparison",
    "discordance",
    "concordance",
    "discordance-one-sample",
    "concordance-one-sample",
    "discordance-both-samples",
    "concordance-both-samples",
    "discordance-exclude-filtered",
    "concordance-exclude-filtered",
];

fn kept(text: &str, run: &str) -> Vec<String> {
    rows(text, "kept")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| {
            row[1]
                .split(',')
                .filter(|at| !at.is_empty())
                .map(|at| at.to_string())
                .collect()
        })
        .unwrap_or_default()
}

/// One run, as `applyFirstRoundOfFiltering` decides it.
fn survivors(text: &str, run: &str) -> Vec<String> {
    let (samples, records) = file(text, "records");
    let (_, comparison) = file(text, "comparison");
    let (sample_arguments, comparison_arguments) = setup(run);
    let selection =
        create_sample_name_inclusion_list(&samples, &sample_arguments).expect("a selection");

    records
        .iter()
        .filter(|record| {
            let others: Vec<Record> = comparison
                .iter()
                .filter(|other| {
                    other.variant.contig == record.variant.contig
                        && other.variant.start == record.variant.start
                })
                .cloned()
                .collect();
            if comparison_arguments.discordance_only
                && !is_discordant(record, &others, &selection, &comparison_arguments)
            {
                return false;
            }
            if comparison_arguments.concordance_only
                && !is_concordant(record, &others, &selection, &comparison_arguments)
            {
                return false;
            }
            true
        })
        .map(|record| record.variant.start.to_string())
        .collect()
}

#[test]
fn every_run_keeps_what_the_reference_kept() {
    let text = golden();
    for run in RUNS {
        assert_eq!(survivors(&text, run), kept(&text, run), "kept/{run}");
    }
}

/// Without a sample, neither flag looks at a genotype.
#[test]
fn without_a_sample_the_flags_are_about_presence_alone() {
    let text = golden();
    // The comparison file has nothing at 400 and something everywhere else.
    assert_eq!(survivors(&text, "discordance"), vec!["400".to_string()]);
    let concordant = survivors(&text, "concordance");
    assert_eq!(concordant.len(), 8);
    assert!(!concordant.contains(&"400".to_string()));
    // Record 300's calls disagree with the other file's, and it is still concordant.
    assert!(concordant.contains(&"300".to_string()));
}

/// The clause that cannot fire, and the flag that cannot matter.
#[test]
fn a_filtered_genotype_never_matches_and_the_flag_changes_nothing() {
    let text = golden();
    // 700 is filtered on both sides, and it is discordant all the same.
    assert!(survivors(&text, "discordance-one-sample").contains(&"700".to_string()));
    // The flag whose only clause is about two filtered genotypes moves no record.
    assert_eq!(
        survivors(&text, "discordance-exclude-filtered"),
        survivors(&text, "discordance-one-sample")
    );
    assert_eq!(
        survivors(&text, "concordance-exclude-filtered"),
        survivors(&text, "concordance-one-sample")
    );
}

/// Allele sets rather than genotypes, and the hom-ref that is never discordant.
#[test]
fn the_comparison_is_of_allele_sets() {
    let text = golden();
    let discordant = survivors(&text, "discordance-one-sample");
    let concordant = survivors(&text, "concordance-one-sample");

    // 200 is `0/1` against `1/0`, and 900 is `1/1` against a haploid `1`: both concordant.
    assert!(concordant.contains(&"200".to_string()));
    assert!(concordant.contains(&"900".to_string()));
    // 500's selected sample is hom-ref, so discordance skips it entirely.
    assert!(!discordant.contains(&"500".to_string()));
    // Selecting the other sample as well makes the same record discordant.
    assert!(survivors(&text, "discordance-both-samples").contains(&"500".to_string()));
    // And they are not negations: 900 is neither discordant nor, with both samples, concordant.
    assert!(!discordant.contains(&"900".to_string()));
    assert!(!survivors(&text, "concordance-both-samples").contains(&"900".to_string()));
}
