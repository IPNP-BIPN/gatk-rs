//! Conformance for `SelectVariants`' output stage against GATK 4.6.2.0, compared as the order the
//! records come out in and as what each of them carries.
//!
//! Golden from `tools/readfilter-conformance/SelectVariantsOutputDump.java`.
//!
//! # What this suite is for
//!
//!  * **trimming moves a record right and the queue puts the file back in order**;
//!  * **the drain is `<=`**, so a record trimmed onto a later record's start goes first;
//!  * **`--set-filtered-gt-to-nocall` recomputes the counts** and keeps the FT;
//!  * **and dropping an annotation removes what the tool itself wrote**.

use gatk_corpus as corpus;
use gatk_engine::subset_alleles::Genotype;
use gatk_engine::variant_context_utils::{Allele, Variant};
use gatk_tools::select_variants::{
    create_sample_name_inclusion_list, drop_annotations, set_filtered_genotypes_to_no_call,
    subset_record, OutputArguments, PendingWriter, Record, SampleArguments, SubsetArguments,
};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/select_variants_output.txt.gz"),
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

fn input(text: &str) -> (Vec<String>, Vec<Record>) {
    let whole = unescape(rows(text, "input").first().expect("an input")[1]);
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
                        dp: by_key.get("DP").and_then(|dp| dp.parse().ok()),
                        // FT and XX are extended attributes; the rest have their own fields.
                        attributes: ["FT", "XX"]
                            .iter()
                            .filter_map(|key| {
                                by_key
                                    .get(key)
                                    .filter(|value| **value != ".")
                                    .map(|value| (key.to_string(), value.to_string()))
                            })
                            .collect(),
                    }
                })
                .collect();
            let attributes = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect();
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

fn names(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

fn setup(run: &str) -> (SampleArguments, SubsetArguments, OutputArguments) {
    let two = || SampleArguments {
        sample_names: names(&["s0", "s1"]),
        ..SampleArguments::default()
    };
    let plain = SubsetArguments::default;
    match run {
        "all-samples" => (
            SampleArguments::default(),
            plain(),
            OutputArguments::default(),
        ),
        "subset" => (two(), plain(), OutputArguments::default()),
        "subset-preserve" => (
            two(),
            SubsetArguments {
                preserve_alleles: true,
                ..plain()
            },
            OutputArguments::default(),
        ),
        "set-filtered-to-nocall" => (
            SampleArguments::default(),
            plain(),
            OutputArguments {
                set_filtered_genotypes_to_no_call: true,
                ..OutputArguments::default()
            },
        ),
        "set-filtered-to-nocall-subset" => (
            two(),
            plain(),
            OutputArguments {
                set_filtered_genotypes_to_no_call: true,
                ..OutputArguments::default()
            },
        ),
        "drop-info" => (
            SampleArguments::default(),
            plain(),
            OutputArguments {
                info_annotations_to_drop: names(&["QD"]),
                ..OutputArguments::default()
            },
        ),
        "drop-info-recomputed" => (
            two(),
            plain(),
            OutputArguments {
                info_annotations_to_drop: names(&["AC"]),
                ..OutputArguments::default()
            },
        ),
        "drop-genotype" => (
            SampleArguments::default(),
            plain(),
            OutputArguments {
                genotype_annotations_to_drop: names(&["XX"]),
                ..OutputArguments::default()
            },
        ),
        "drop-both" => (
            SampleArguments::default(),
            plain(),
            OutputArguments {
                info_annotations_to_drop: names(&["QD"]),
                genotype_annotations_to_drop: names(&["XX"]),
                ..OutputArguments::default()
            },
        ),
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 9] = [
    "all-samples",
    "subset",
    "subset-preserve",
    "set-filtered-to-nocall",
    "set-filtered-to-nocall-subset",
    "drop-info",
    "drop-info-recomputed",
    "drop-genotype",
    "drop-both",
];

/// One whole run of `apply`, ending with the drain `onTraversalSuccess` does.
fn written(records: &[Record], header_samples: &[String], run: &str) -> Vec<Record> {
    let (sample_arguments, subset_arguments, output_arguments) = setup(run);
    let selection =
        create_sample_name_inclusion_list(header_samples, &sample_arguments).expect("a selection");
    let mut writer = PendingWriter::new();
    let mut out = Vec::new();
    for record in records {
        out.extend(writer.drain_before(&record.variant.contig, record.variant.start));
        let mut result = subset_record(record, &selection, &subset_arguments).expect("a subset");
        if output_arguments.set_filtered_genotypes_to_no_call {
            set_filtered_genotypes_to_no_call(&mut result);
        }
        drop_annotations(&mut result, &output_arguments);
        writer.add(result);
    }
    out.extend(writer.drain());
    out
}

fn order(text: &str, run: &str) -> Vec<String> {
    rows(text, "order")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| row[1].split(',').map(|at| at.to_string()).collect())
        .unwrap_or_default()
}

fn lines(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

#[test]
fn every_run_writes_in_the_reference_s_order() {
    let text = golden();
    let (samples, records) = input(&text);
    for run in RUNS {
        let ours: Vec<String> = written(&records, &samples, run)
            .iter()
            .map(|record| format!("{}:{}", record.variant.contig, record.variant.start))
            .collect();
        assert_eq!(ours, order(&text, run), "order/{run}");
    }
}

/// Trimming moves the first record onto the third's start, and the queue writes it first.
#[test]
fn a_trimmed_record_overtakes_the_one_that_followed_it() {
    let text = golden();
    let (samples, records) = input(&text);

    // Untouched, the file is in the order it was read.
    assert_eq!(
        order(&text, "all-samples"),
        vec!["chr1:100", "chr1:101", "chr1:104", "chr1:200", "chr2:50", "chr2:51"]
    );
    // Subsetting trims, and the record at 100 lands on 104, after the record at 101.
    assert_eq!(
        order(&text, "subset"),
        vec!["chr1:101", "chr1:104", "chr1:104", "chr1:200", "chr2:51", "chr2:54"]
    );
    // --preserve-alleles puts the original order back, which is how we know the order is the
    // trimming's doing and not the queue's.
    assert_eq!(order(&text, "subset-preserve"), order(&text, "all-samples"));

    let ours = written(&records, &samples, "subset");
    // The first of the two records at 104 is the one that was trimmed: its alternate is `AC`.
    assert_eq!(ours[1].variant.start, 104);
    assert_eq!(ours[1].variant.alleles[1].bases, b"AC");
    assert_eq!(ours[2].variant.start, 104);
    assert_eq!(ours[2].variant.alleles[1].bases, b"G");
}

/// The replacement recomputes the counts, which is visible as an AF the input never had.
#[test]
fn no_calling_a_filtered_genotype_recomputes_the_counts() {
    let text = golden();
    let (samples, records) = input(&text);

    let before: Vec<&str> = records[3]
        .variant
        .attributes
        .iter()
        .map(|(key, _)| key.as_str())
        .collect();
    assert!(!before.contains(&"AF"), "the input has no AF");

    let ours = written(&records, &samples, "set-filtered-to-nocall");
    let record = ours.iter().find(|r| r.variant.start == 200).expect("200");
    let info: BTreeMap<&str, &str> = record
        .variant
        .attributes
        .iter()
        .map(|(key, value)| (key.as_str(), value.as_str()))
        .collect();
    assert_eq!(info.get("AN"), Some(&"2"));
    assert_eq!(info.get("AC"), Some(&"1"));
    assert_eq!(info.get("AF"), Some(&"0.500"));

    // The two replaced genotypes are no-calls that still carry the filter that no-called them.
    assert!(record.variant.genotypes[0]
        .alleles
        .iter()
        .all(Option::is_none));
    assert!(record.variant.genotypes[0]
        .attributes
        .iter()
        .any(|(key, value)| key == "FT" && value == "LowGQ"));

    // The golden's own line for that record, which carries both halves of this.
    let theirs = lines(&text, "set-filtered-to-nocall")
        .into_iter()
        .find(|line| line.split('\t').nth(1) == Some("200"))
        .expect("the record at 200");
    assert!(theirs.contains("AC=1;AF=0.500;AN=2"), "{theirs}");
    assert!(theirs.contains("./.:10:LowGQ"), "{theirs}");
}

/// Dropping an annotation removes what the tool would have written, not what it read.
#[test]
fn dropping_an_annotation_reaches_the_recomputed_one() {
    let text = golden();
    let (samples, records) = input(&text);

    let ours = written(&records, &samples, "drop-info-recomputed");
    for record in &ours {
        let keys: Vec<&str> = record
            .variant
            .attributes
            .iter()
            .map(|(key, _)| key.as_str())
            .collect();
        assert!(!keys.contains(&"AC"), "AC survived the drop");
        // The recomputed AF is still there, which is what shows AC was the recomputed one too.
        assert!(keys.contains(&"AF"));
    }

    // A genotype annotation is dropped from every sample, and the fields that are not extended
    // attributes stay.
    let dropped = written(&records, &samples, "drop-genotype");
    for record in &dropped {
        for genotype in &record.variant.genotypes {
            assert!(!genotype.attributes.iter().any(|(key, _)| key == "XX"));
            assert!(genotype.gq.is_some());
            assert!(genotype.dp.is_some());
        }
    }
}
