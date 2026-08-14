//! Conformance for `SelectVariants`' record filters against GATK 4.6.2.0, compared as which
//! records each run wrote and as the class and message of every refusal.
//!
//! Golden from `tools/readfilter-conformance/SelectVariantsFiltersDump.java`.
//!
//! # What this suite is for
//!
//!  * **the filtered-genotype fraction is an integer division**, so its argument does nothing;
//!  * **the no-call fraction, one line below, casts and works**;
//!  * **`--invert-select` inverts each expression** rather than their disjunction;
//!  * **and an expression over a per-allele annotation is a refusal**, not a false.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{Allele, Variant};
use gatk_tools::select_variants::{
    create_sample_name_inclusion_list, keeps_after_subset, keeps_before_subset, subset_record,
    AlleleRestriction, FilterArguments, FilterRecord, Record, SampleArguments, SelectError,
    SubsetArguments, VariantType,
};
use std::collections::HashMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_filters.txt.gz"),
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

/// The input file as the two views the port takes of it.
fn input(text: &str) -> (Vec<String>, Vec<(Record, FilterRecord)>) {
    let whole = unescape(rows(text, "input").first().expect("an input")[1]);
    let samples: Vec<String> = whole
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header")
        .split('\t')
        .skip(9)
        .map(|name| name.to_string())
        .collect();

    // The header decides which INFO fields decode to a LIST rather than a scalar, and that is what
    // makes an expression over one of them a refusal: `Number=A` is one value per alternate.
    let per_allele: Vec<String> = whole
        .lines()
        .filter(|line| line.starts_with("##INFO=<") && line.contains("Number=A,"))
        .filter_map(|line| {
            line.split_once("ID=")
                .and_then(|(_, rest)| rest.split(',').next())
                .map(|id| id.to_string())
        })
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
            let mut genotypes = Vec::new();
            let mut genotype_fields = Vec::new();
            for index in 0..samples.len() {
                let values: Vec<&str> = field[9 + index].split(':').collect();
                let by_key: HashMap<String, String> = keys
                    .iter()
                    .zip(values.iter())
                    .map(|(key, value)| (key.to_string(), value.to_string()))
                    .collect();
                let call = by_key
                    .get("GT")
                    .cloned()
                    .unwrap_or_else(|| "./.".to_string());
                genotypes.push(Genotype {
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
                        .filter(|ft| *ft != "." && *ft != "PASS")
                        .map(|ft| vec![("FT".to_string(), ft.clone())])
                        .unwrap_or_default(),
                });
                genotype_fields.push(by_key);
            }
            let mut attributes = Vec::new();
            let mut info = HashMap::new();
            for entry in field[7].split(';') {
                if let Some((key, value)) = entry.split_once('=') {
                    attributes.push((key.to_string(), value.to_string()));
                    // A per-allele field is an ArrayList once decoded, and an expression sees its
                    // `toString`: `[2]`, or `[1, 1]` for two alternates. Nothing numeric compares
                    // against that, which is the refusal the golden holds.
                    let rendered = if per_allele.iter().any(|id| id == key) {
                        format!("[{}]", value.split(',').collect::<Vec<_>>().join(", "))
                    } else {
                        value.to_string()
                    };
                    info.insert(key.to_string(), rendered);
                }
            }
            let filters = match field[6] {
                "." | "PASS" => Vec::new(),
                names => names.split(';').map(|name| name.to_string()).collect(),
            };
            (
                Record {
                    variant: Variant {
                        contig: field[0].to_string(),
                        start: field[1].parse().expect("a position"),
                        stop: field[1].parse::<i32>().expect("a position") + field[3].len() as i32
                            - 1,
                        alleles,
                        genotypes,
                        attributes,
                    },
                    samples: samples.clone(),
                },
                FilterRecord {
                    id: field[2].to_string(),
                    filters,
                    info,
                    genotype_fields,
                },
            )
        })
        .collect();
    (samples, records)
}

fn names(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

/// The arguments of each run, which the golden does not carry.
fn setup(run: &str) -> (SampleArguments, FilterArguments) {
    let base = FilterArguments::default;
    let no_samples = SampleArguments::default;
    let filters = |arguments: FilterArguments| (no_samples(), arguments);
    match run {
        "no-filter" => filters(base()),
        "type-snp" => filters(FilterArguments {
            types_to_include: vec![VariantType::Snp],
            ..base()
        }),
        "type-indel" => filters(FilterArguments {
            types_to_include: vec![VariantType::Indel],
            ..base()
        }),
        "type-exclude-snp" => filters(FilterArguments {
            types_to_exclude: vec![VariantType::Snp],
            ..base()
        }),
        "type-include-and-exclude" => filters(FilterArguments {
            types_to_include: vec![VariantType::Snp, VariantType::Indel],
            types_to_exclude: vec![VariantType::Snp],
            ..base()
        }),
        "biallelic" => filters(FilterArguments {
            allele_restriction: AlleleRestriction::Biallelic,
            ..base()
        }),
        "multiallelic" => filters(FilterArguments {
            allele_restriction: AlleleRestriction::Multiallelic,
            ..base()
        }),
        "max-indel-size" => filters(FilterArguments {
            max_indel_size: 10,
            ..base()
        }),
        "min-indel-size" => filters(FilterArguments {
            min_indel_size: 5,
            ..base()
        }),
        "keep-ids" => filters(FilterArguments {
            keep_ids: names(&["rs100"]),
            ..base()
        }),
        "exclude-ids" => filters(FilterArguments {
            exclude_ids: names(&["rs100"]),
            ..base()
        }),
        "exclude-filtered" => filters(FilterArguments {
            exclude_filtered: true,
            ..base()
        }),
        "max-filtered-genotypes" => filters(FilterArguments {
            max_filtered_genotypes: 1,
            ..base()
        }),
        "max-fraction-filtered-genotypes" => filters(FilterArguments {
            max_fraction_filtered_genotypes: 0.1,
            ..base()
        }),
        "max-nocall-number" => filters(FilterArguments {
            max_nocall_number: 1,
            ..base()
        }),
        "max-nocall-fraction" => filters(FilterArguments {
            max_nocall_fraction: 0.1,
            ..base()
        }),
        "exclude-non-variants" => filters(FilterArguments {
            exclude_non_variants: true,
            ..base()
        }),
        "exclude-non-variants-subset" => (
            SampleArguments {
                sample_names: names(&["s0", "s1"]),
                ..no_samples()
            },
            FilterArguments {
                exclude_non_variants: true,
                ..base()
            },
        ),
        "select-one" => filters(FilterArguments {
            select_expressions: names(&["QD > 10.0"]),
            ..base()
        }),
        "select-inverted" => filters(FilterArguments {
            select_expressions: names(&["QD > 10.0"]),
            invert_select: true,
            ..base()
        }),
        "select-genotype" => filters(FilterArguments {
            select_genotype_expressions: names(&["GQ > 55"]),
            ..base()
        }),
        "select-two" => filters(FilterArguments {
            select_expressions: names(&["QD > 10.0", "AC > 1"]),
            ..base()
        }),
        "select-two-inverted" => filters(FilterArguments {
            select_expressions: names(&["QD > 10.0", "AC > 1"]),
            invert_select: true,
            ..base()
        }),
        "select-ac-after-subset" => (
            SampleArguments {
                sample_names: names(&["s0"]),
                ..no_samples()
            },
            FilterArguments {
                select_expressions: names(&["AC > 1"]),
                ..base()
            },
        ),
        other => panic!("no setup for {other}"),
    }
}

/// The runs whose whole answer is the list of records written.
const RUNS: [&str; 21] = [
    "no-filter",
    "type-snp",
    "type-indel",
    "type-exclude-snp",
    "type-include-and-exclude",
    "biallelic",
    "multiallelic",
    "max-indel-size",
    "min-indel-size",
    "keep-ids",
    "exclude-ids",
    "exclude-filtered",
    "max-filtered-genotypes",
    "max-fraction-filtered-genotypes",
    "max-nocall-number",
    "max-nocall-fraction",
    "exclude-non-variants",
    "exclude-non-variants-subset",
    "select-one",
    "select-inverted",
    "select-genotype",
];

fn kept(text: &str, run: &str) -> Vec<String> {
    rows(text, "kept")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| {
            row[1]
                .split(',')
                .filter(|position| !position.is_empty())
                .map(|position| position.to_string())
                .collect()
        })
        .unwrap_or_default()
}

fn refusal(text: &str, run: &str) -> Option<String> {
    rows(text, "error")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| unescape(row[1]))
}

/// Both stages, as `apply` runs them, for one run.
fn survivors(
    records: &[(Record, FilterRecord)],
    header_samples: &[String],
    run: &str,
) -> Result<Vec<String>, SelectError> {
    let (sample_arguments, filter_arguments) = setup(run);
    let selection =
        create_sample_name_inclusion_list(header_samples, &sample_arguments).expect("a selection");
    let mut written = Vec::new();
    for (record, filter_record) in records {
        if !keeps_before_subset(record, filter_record, &filter_arguments, &selection)? {
            continue;
        }
        let subset =
            subset_record(record, &selection, &SubsetArguments::default()).expect("a subset");
        if !keeps_after_subset(&subset, filter_record, &filter_arguments)? {
            continue;
        }
        written.push(record.variant.start.to_string());
    }
    Ok(written)
}

#[test]
fn every_run_writes_the_records_the_reference_wrote() {
    let text = golden();
    let (samples, records) = input(&text);
    for run in RUNS {
        let ours = survivors(&records, &samples, run)
            .unwrap_or_else(|error| panic!("{run}: {}", error.message()));
        assert_eq!(ours, kept(&text, run), "kept/{run}");
    }
}

/// The refusals, which the golden holds instead of a record list.
#[test]
fn an_expression_over_a_per_allele_annotation_is_a_refusal() {
    let text = golden();
    let (samples, records) = input(&text);

    for run in [
        "select-two",
        "select-two-inverted",
        "select-ac-after-subset",
    ] {
        let expected = refusal(&text, run).unwrap_or_else(|| panic!("no refusal for {run}"));
        let error = survivors(&records, &samples, run).expect_err(run);
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            expected,
            "error/{run}"
        );
    }

    // The unparseable one never reaches a record: it is the argument parser's refusal.
    let expected = refusal(&text, "select-unparseable").expect("a refusal");
    let error = SelectError::Unparseable {
        index: 0,
        text: "QD >".to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        expected
    );
}

/// One line apart, one casts and one does not.
#[test]
fn the_filtered_fraction_is_an_integer_division_and_the_nocall_fraction_is_not() {
    let text = golden();
    let (samples, records) = input(&text);

    // Record 800 has two filtered genotypes out of four, which is a fraction of 0.5.
    let by_count = survivors(&records, &samples, "max-filtered-genotypes").expect("kept");
    assert!(
        !by_count.contains(&"800".to_string()),
        "the count gate drops it"
    );

    // 0.5 is above 0.1, and yet the fraction gate keeps it, because 2 / 4 is 0 in int arithmetic.
    let by_fraction =
        survivors(&records, &samples, "max-fraction-filtered-genotypes").expect("kept");
    assert!(
        by_fraction.contains(&"800".to_string()),
        "the fraction gate does nothing"
    );
    assert_eq!(by_fraction, kept(&text, "no-filter"));

    // The no-call fraction, written one line below in the reference, casts: 1/4 is above 0.1.
    let no_calls = survivors(&records, &samples, "max-nocall-fraction").expect("kept");
    assert!(!no_calls.contains(&"800".to_string()));
    assert!(!no_calls.contains(&"900".to_string()));
}

/// Inverting one expression keeps exactly the records it rejected, which is not the complement of
/// a disjunction of several.
#[test]
fn invert_select_inverts_each_expression_and_not_their_disjunction() {
    let text = golden();
    let (samples, records) = input(&text);
    let plain = survivors(&records, &samples, "select-one").expect("kept");
    let inverted = survivors(&records, &samples, "select-inverted").expect("kept");
    let all = kept(&text, "no-filter");

    assert_eq!(plain, kept(&text, "select-one"));
    assert_eq!(inverted, kept(&text, "select-inverted"));
    // Here, with one expression, the two are complements of each other.
    let mut union = plain.clone();
    union.extend(inverted.iter().cloned());
    union.sort();
    assert_eq!(union, all);
}

/// The spanning deletion, and the record only the unselected samples vary at.
#[test]
fn exclude_non_variants_runs_after_the_subset() {
    let text = golden();
    let (samples, records) = input(&text);

    let whole = survivors(&records, &samples, "exclude-non-variants").expect("kept");
    // 700 is `A -> *`, which is not a variant even though a genotype calls it.
    assert!(!whole.contains(&"700".to_string()));
    // 600 is carried by s2 and s3, which this run keeps.
    assert!(whole.contains(&"600".to_string()));

    let subset = survivors(&records, &samples, "exclude-non-variants-subset").expect("kept");
    // The same record, with only s0 and s1 selected, is no longer polymorphic.
    assert!(!subset.contains(&"600".to_string()));
    assert_eq!(subset, kept(&text, "exclude-non-variants-subset"));
}
