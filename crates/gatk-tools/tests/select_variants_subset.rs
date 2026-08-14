//! Conformance for `SelectVariants`' record subsetting against GATK 4.6.2.0, compared as the
//! alleles, the INFO attributes and every genotype field of every record of every run.
//!
//! Golden from `tools/readfilter-conformance/SelectVariantsSubsetDump.java`.
//!
//! # What this suite is for
//!
//!  * **the INFO DP is replaced by a sum over the kept, unfiltered genotypes**;
//!  * **a filtered genotype is invisible to AN and AC** as well as to that sum;
//!  * **MLEAC and MLEAF are stripped** whenever the record is rewritten;
//!  * **and `--remove-unused-alternates` rewrites AD and PL** with the allele list.
//!
//! # What is compared, and what is not
//!
//! The golden holds whole VCF lines; this compares the alleles, the INFO map and the genotype
//! fields by name. The FORMAT column's own order is the writer's doing (`GT` first, then the union
//! of the keys, sorted) and htsjdk-rs has suites over that; reproducing it here would be testing
//! the encoder rather than the tool. `FT` is compared where the reference wrote it, and its absence
//! is read as `PASS`, which is what an unfiltered genotype means.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{Allele, Variant};
use gatk_tools::select_variants::{
    create_sample_name_inclusion_list, subset_record, Record, SampleArguments, SubsetArguments,
};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_subset.txt.gz"),
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

/// One VCF record, as much of it as this port models.
#[derive(Debug, Clone, PartialEq)]
struct Line {
    reference: String,
    alternates: Vec<String>,
    info: BTreeMap<String, String>,
    /// Per sample, the FORMAT keys and their values.
    calls: Vec<BTreeMap<String, String>>,
}

fn parse_line(line: &str, samples: &[String]) -> Line {
    let field: Vec<&str> = line.split('\t').collect();
    let alternates = if field[4] == "." {
        Vec::new()
    } else {
        field[4].split(',').map(|alt| alt.to_string()).collect()
    };
    let mut info = BTreeMap::new();
    if field[7] != "." {
        for entry in field[7].split(';') {
            match entry.split_once('=') {
                Some((key, value)) => info.insert(key.to_string(), value.to_string()),
                None => info.insert(entry.to_string(), String::new()),
            };
        }
    }
    let keys: Vec<&str> = field[8].split(':').collect();
    let calls = (0..samples.len())
        .map(|index| {
            let values: Vec<&str> = field[9 + index].split(':').collect();
            keys.iter()
                .zip(values.iter())
                .filter(|(_, value)| **value != ".")
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect()
        })
        .collect();
    Line {
        reference: field[3].to_string(),
        alternates,
        info,
        calls,
    }
}

/// The input file, parsed into the port's model.
fn input(text: &str) -> (Vec<String>, Vec<Record>) {
    let whole = unescape(rows(text, "input").first().expect("an input")[1]);
    let header = whole
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header");
    let samples: Vec<String> = header.split('\t').skip(9).map(|s| s.to_string()).collect();

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
                        pl: numbers(by_key.get("PL").copied()),
                        gq: by_key.get("GQ").and_then(|gq| gq.parse().ok()),
                        ad: numbers(by_key.get("AD").copied()),
                        dp: by_key.get("DP").and_then(|dp| dp.parse().ok()),
                        attributes: by_key
                            .get("FT")
                            .filter(|ft| **ft != ".")
                            .map(|ft| vec![("FT".to_string(), ft.to_string())])
                            .unwrap_or_default(),
                    }
                })
                .collect();
            let mut attributes = Vec::new();
            for entry in field[7].split(';') {
                if let Some((key, value)) = entry.split_once('=') {
                    attributes.push((key.to_string(), value.to_string()));
                }
            }
            Record {
                variant: Variant {
                    contig: field[0].to_string(),
                    start: field[1].parse().expect("a position"),
                    stop: field[1].parse::<i32>().expect("a position") + field[3].len() as i32 - 1,
                    alleles,
                    genotypes,
                    attributes,
                },
                samples: samples.clone(),
            }
        })
        .collect();
    (samples, records)
}

fn numbers(value: Option<&str>) -> Option<Vec<i32>> {
    let value = value?;
    if value == "." {
        return None;
    }
    value
        .split(',')
        .map(|part| part.parse().ok())
        .collect::<Option<Vec<i32>>>()
}

/// What the port produced, in the same shape as a parsed golden line.
fn rendered(record: &Record) -> Line {
    let calls = record
        .variant
        .genotypes
        .iter()
        .map(|genotype| {
            let mut call: BTreeMap<String, String> = BTreeMap::new();
            call.insert(
                "GT".to_string(),
                genotype
                    .alleles
                    .iter()
                    .map(|allele| match allele {
                        Some(index) => index.to_string(),
                        None => ".".to_string(),
                    })
                    .collect::<Vec<_>>()
                    .join("/"),
            );
            if let Some(gq) = genotype.gq {
                call.insert("GQ".to_string(), gq.to_string());
            }
            if let Some(dp) = genotype.dp {
                call.insert("DP".to_string(), dp.to_string());
            }
            if let Some(ad) = &genotype.ad {
                call.insert("AD".to_string(), join(ad));
            }
            if let Some(pl) = &genotype.pl {
                call.insert("PL".to_string(), join(pl));
            }
            for (key, value) in &genotype.attributes {
                call.insert(key.clone(), value.clone());
            }
            call
        })
        .collect();
    Line {
        reference: String::from_utf8(record.variant.alleles[0].bases.clone()).expect("bases"),
        alternates: record.variant.alleles[1..]
            .iter()
            .map(|allele| String::from_utf8(allele.bases.clone()).expect("bases"))
            .collect(),
        info: record
            .variant
            .attributes
            .iter()
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect(),
        calls,
    }
}

fn join(values: &[i32]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// The arguments of each run, which the golden does not carry.
fn setup(run: &str) -> (SampleArguments, SubsetArguments) {
    let names =
        |values: &[&str]| -> Vec<String> { values.iter().map(|value| value.to_string()).collect() };
    let two = SampleArguments {
        sample_names: names(&["s0", "s1"]),
        ..SampleArguments::default()
    };
    match run {
        "all-samples" => (SampleArguments::default(), SubsetArguments::default()),
        "subset-two" => (two, SubsetArguments::default()),
        "subset-two-remove-unused" => (
            two,
            SubsetArguments {
                remove_unused_alternates: true,
                ..SubsetArguments::default()
            },
        ),
        "subset-two-preserve" => (
            two,
            SubsetArguments {
                remove_unused_alternates: true,
                preserve_alleles: true,
                ..SubsetArguments::default()
            },
        ),
        "subset-two-keep-original" => (
            two,
            SubsetArguments {
                keep_original_chr_counts: true,
                keep_original_depth: true,
                ..SubsetArguments::default()
            },
        ),
        "subset-two-remove-unused-keep-original" => (
            two,
            SubsetArguments {
                remove_unused_alternates: true,
                keep_original_chr_counts: true,
                ..SubsetArguments::default()
            },
        ),
        "subset-one" => (
            SampleArguments {
                sample_names: names(&["s2"]),
                ..SampleArguments::default()
            },
            SubsetArguments::default(),
        ),
        "subset-carriers" => (
            SampleArguments {
                sample_names: names(&["s2", "s3"]),
                ..SampleArguments::default()
            },
            SubsetArguments::default(),
        ),
        "all-remove-unused" => (
            SampleArguments::default(),
            SubsetArguments {
                remove_unused_alternates: true,
                ..SubsetArguments::default()
            },
        ),
        "exclude-one" => (
            SampleArguments {
                exclude_sample_names: names(&["s3"]),
                ..SampleArguments::default()
            },
            SubsetArguments::default(),
        ),
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 10] = [
    "all-samples",
    "subset-two",
    "subset-two-remove-unused",
    "subset-two-preserve",
    "subset-two-keep-original",
    "subset-two-remove-unused-keep-original",
    "subset-one",
    "subset-carriers",
    "all-remove-unused",
    "exclude-one",
];

/// The reference writes no `FT` for a record whose genotypes are all unfiltered, and `PASS` for the
/// unfiltered genotypes of a record where one is. The port carries the field as it was read, so an
/// absent `FT` and a `PASS` are the same answer here.
fn without_passing_filters(mut call: BTreeMap<String, String>) -> BTreeMap<String, String> {
    if call.get("FT").map(String::as_str) == Some("PASS") {
        call.remove("FT");
    }
    call
}

fn compare(ours: &Line, theirs: &Line, what: &str) {
    assert_eq!(ours.reference, theirs.reference, "REF/{what}");
    assert_eq!(ours.alternates, theirs.alternates, "ALT/{what}");
    assert_eq!(ours.info, theirs.info, "INFO/{what}");
    assert_eq!(ours.calls.len(), theirs.calls.len(), "samples/{what}");
    for (index, (ours, theirs)) in ours.calls.iter().zip(theirs.calls.iter()).enumerate() {
        assert_eq!(
            without_passing_filters(ours.clone()),
            without_passing_filters(theirs.clone()),
            "call {index}/{what}"
        );
    }
}

#[test]
fn every_subset_is_the_reference_s() {
    let text = golden();
    let (header_samples, records) = input(&text);

    for run in RUNS {
        let (sample_arguments, subset_arguments) = setup(run);
        let selection = create_sample_name_inclusion_list(&header_samples, &sample_arguments)
            .unwrap_or_else(|error| panic!("{run}: {}", error.message()));

        let expected: Vec<String> = rows(&text, "vcfline")
            .into_iter()
            .filter(|row| row[0] == run)
            .map(|row| unescape(row[1]))
            .collect();
        assert_eq!(expected.len(), records.len(), "records/{run}");

        for (record, line) in records.iter().zip(expected.iter()) {
            let ours = subset_record(record, &selection, &subset_arguments)
                .unwrap_or_else(|error| panic!("{run}: {error:?}"));
            let theirs = parse_line(line, &ours.samples);
            compare(
                &rendered(&ours),
                &theirs,
                &format!("{run}@{}", record.variant.start),
            );
        }
    }
}

/// The DP a subset record carries is not the DP it arrived with.
#[test]
fn the_info_depth_becomes_a_sum_over_the_kept_unfiltered_genotypes() {
    let text = golden();
    let (header_samples, records) = input(&text);
    let (sample_arguments, subset_arguments) = setup("subset-two");
    let selection =
        create_sample_name_inclusion_list(&header_samples, &sample_arguments).expect("selected");

    // Every record arrived saying DP=200.
    for record in &records {
        assert_eq!(
            record
                .variant
                .attributes
                .iter()
                .find(|(key, _)| key == "DP")
                .map(|(_, value)| value.as_str()),
            Some("200")
        );
    }

    let depth = |record: &Record| -> String {
        subset_record(record, &selection, &subset_arguments)
            .expect("subset")
            .variant
            .attributes
            .iter()
            .find(|(key, _)| key == "DP")
            .map(|(_, value)| value.clone())
            .expect("a DP")
    };

    // 10 + 20 for the first record, and 10 alone for the third, whose second genotype is filtered.
    assert_eq!(depth(&records[0]), "30");
    assert_eq!(depth(&records[2]), "10");
    // The last record's genotypes have no DP at all, so the record keeps its own.
    assert_eq!(depth(&records[4]), "200");
}

/// A filtered genotype is not a called one, so it counts towards neither AN nor AC.
#[test]
fn a_filtered_genotype_is_invisible_to_the_chromosome_counts() {
    let text = golden();
    let (header_samples, records) = input(&text);
    let selection = create_sample_name_inclusion_list(
        &header_samples,
        &SampleArguments {
            sample_names: vec!["s0".to_string(), "s1".to_string()],
            ..SampleArguments::default()
        },
    )
    .expect("selected");

    let ours = subset_record(&records[2], &selection, &SubsetArguments::default()).expect("subset");
    let info: BTreeMap<String, String> = ours.variant.attributes.iter().cloned().collect();
    // Two kept genotypes, four chromosomes, but the second is filtered: AN is 2 and not 4.
    assert_eq!(info.get("AN").map(String::as_str), Some("2"));
    assert_eq!(info.get("AC").map(String::as_str), Some("1"));
    assert_eq!(info.get("AF").map(String::as_str), Some("0.500"));
}

/// The tags that describe a calling no longer in the record.
#[test]
fn the_maximum_likelihood_tags_are_stripped_and_the_counts_are_replaced() {
    let text = golden();
    let (header_samples, records) = input(&text);
    let before: BTreeMap<String, String> = records[0].variant.attributes.iter().cloned().collect();
    assert!(before.contains_key("MLEAC") && before.contains_key("AC"));

    let selection = create_sample_name_inclusion_list(
        &header_samples,
        &SampleArguments {
            sample_names: vec!["s0".to_string(), "s1".to_string()],
            ..SampleArguments::default()
        },
    )
    .expect("selected");
    let ours = subset_record(&records[0], &selection, &SubsetArguments::default()).expect("subset");
    let after: BTreeMap<String, String> = ours.variant.attributes.iter().cloned().collect();
    assert!(!after.contains_key("MLEAC"), "MLEAC survived");
    assert!(!after.contains_key("MLEAF"), "MLEAF survived");
    // AC is replaced rather than stripped, and the alternate nobody kept now counts zero.
    assert_eq!(after.get("AC").map(String::as_str), Some("0"));
    assert_eq!(after.get("AF").map(String::as_str), Some("0.00"));
}

/// A whole-cohort run with no flag returns the record itself.
#[test]
fn selecting_everything_rewrites_nothing() {
    let text = golden();
    let (header_samples, records) = input(&text);
    let selection = create_sample_name_inclusion_list(&header_samples, &SampleArguments::default())
        .expect("selected");
    for record in &records {
        let ours =
            subset_record(record, &selection, &SubsetArguments::default()).expect("unchanged");
        assert_eq!(&ours, record);
    }
}
