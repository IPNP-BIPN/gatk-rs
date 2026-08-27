//! Conformance for `CombineGVCFs` against GATK 4.6.2.0, compared as the records of every run.
//!
//! Golden from `tools/readfilter-conformance/CombineGVCFsDump.java`.
//!
//! The likelihoods the merge expands to the joint allele set are in the golden and are not
//! reproduced: they come from the genotype machinery, which is not ported. What is compared is
//! which records the output has, which sample carries data on each, and what each contributes.
//!
//! # What this suite is for
//!
//!  * **the records being the union of every input's edges**;
//!  * **a variant in one sample cutting the others' blocks**;
//!  * **each sample keeping its own quality**;
//!  * **a finished sample keeping its column and losing its fields**;
//!  * **the two band arguments, and base-pair resolution winning over the grid**;
//!  * **and the same file twice being refused.**

use gatk_corpus as corpus;
use gatk_tools::combine_gvcfs::{
    boundaries, check_inputs, combine, BandArguments, CombineError, Record,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/combine_gvcfs.txt.gz"),
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

fn field(genotype: &str, format: &str, key: &str) -> Option<i32> {
    let keys: Vec<&str> = format.split(':').collect();
    let values: Vec<&str> = genotype.split(':').collect();
    keys.iter()
        .position(|name| *name == key)
        .and_then(|at| values.get(at))
        .and_then(|text| text.parse().ok())
}

/// One input GVCF, read as its records.
fn input(text: &str, sample: &str) -> Vec<Record> {
    section(text, "vcf", sample)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[1].parse().expect("a position");
            let end = columns[7]
                .split(';')
                .find_map(|part| part.strip_prefix("END="))
                .map(|value| value.parse().expect("an end"))
                .unwrap_or(start);
            Record {
                sample: sample.to_string(),
                start,
                end,
                alternates: columns[4]
                    .split(',')
                    .filter(|allele| *allele != "<NON_REF>")
                    .map(str::to_string)
                    .collect(),
                genotype_quality: field(columns[9], columns[8], "GQ").expect("a quality"),
            }
        })
        .collect()
}

/// One run's records, as the span and each sample's quality.
fn measured(text: &str, label: &str) -> Vec<(i32, i32, Vec<Option<i32>>)> {
    let body = section(text, "out", label);
    let mut lines = body.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let sample_count = header.len() - 9;
    lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[1].parse().expect("a position");
            let end = columns[7]
                .split(';')
                .find_map(|part| part.strip_prefix("END="))
                .map(|value| value.parse().expect("an end"))
                .unwrap_or(start);
            (
                start,
                end,
                (0..sample_count)
                    .map(|index| field(columns[9 + index], columns[8], "GQ"))
                    .collect(),
            )
        })
        .collect()
}

fn produced(
    records: &[Record],
    samples: &[String],
    arguments: &BandArguments,
) -> Vec<(i32, i32, Vec<Option<i32>>)> {
    combine(records, samples, arguments)
        .into_iter()
        .map(|record| (record.start, record.end, record.qualities))
        .collect()
}

fn samples(names: &[&str]) -> Vec<String> {
    names.iter().map(|name| name.to_string()).collect()
}

fn all_records(text: &str, names: &[&str]) -> Vec<Record> {
    names.iter().flat_map(|name| input(text, name)).collect()
}

#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let three = samples(&["s1", "s2", "s3"]);
    let mut compared = 0;
    for (label, names, arguments) in [
        (
            "three-samples",
            &["s1", "s2", "s3"][..],
            BandArguments::default(),
        ),
        ("two-samples", &["s1", "s2"][..], BandArguments::default()),
        ("one-sample", &["s1"][..], BandArguments::default()),
        (
            "base-pair-resolution",
            &["s1", "s2", "s3"][..],
            BandArguments {
                base_pair_resolution: true,
                ..BandArguments::default()
            },
        ),
        (
            "break-bands-100",
            &["s1", "s2", "s3"][..],
            BandArguments {
                break_bands_at_multiples_of: 100,
                ..BandArguments::default()
            },
        ),
        (
            "break-bands-50",
            &["s1", "s2", "s3"][..],
            BandArguments {
                break_bands_at_multiples_of: 50,
                ..BandArguments::default()
            },
        ),
        (
            "both-band-arguments",
            &["s1", "s2", "s3"][..],
            BandArguments {
                base_pair_resolution: true,
                break_bands_at_multiples_of: 100,
            },
        ),
        // Calling the genotypes changes the GT and nothing this suite compares.
        (
            "call-genotypes",
            &["s1", "s2", "s3"][..],
            BandArguments::default(),
        ),
    ] {
        let records = all_records(&text, names);
        let names = samples(names);
        assert_eq!(
            produced(&records, &names, &arguments),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    let _ = three;
    assert_eq!(compared, 8, "the runs the port reproduces");
}

/// A block no input broke is broken wherever any other input has an edge.
#[test]
fn the_records_are_the_union_of_every_inputs_edges() {
    let text = golden();
    let records = all_records(&text, &["s1", "s2", "s3"]);
    // s1's first block runs 1000 to 1199 unbroken.
    let first = records
        .iter()
        .find(|record| record.sample == "s1" && record.start == 1000)
        .expect("s1's first block");
    assert_eq!(first.end, 1199);
    // The output breaks it at 1100, where s2 has an edge, and at 1151, where s3 ends.
    let out = measured(&text, "three-samples");
    assert_eq!(out[0].0, 1000);
    assert_eq!(out[0].1, 1099);
    assert_eq!(out[1].0, 1100);
    assert_eq!(out[1].1, 1150);
    assert_eq!(out[2].0, 1151);
    assert_eq!(out[2].1, 1199);
    // And those are exactly the boundaries the port computes.
    assert_eq!(
        boundaries(&records, &BandArguments::default())[..3],
        [1000, 1100, 1151]
    );
    // Two samples alone give fewer edges, because s3's are gone.
    let two = measured(&text, "two-samples");
    assert!(!two.iter().any(|(start, ..)| *start == 1151));
}

/// It cuts the others' blocks at its own base, and every sample is written there.
#[test]
fn a_variant_cuts_the_other_samples_blocks() {
    let text = golden();
    let out = measured(&text, "three-samples");
    let at_variant = out
        .iter()
        .find(|(start, ..)| *start == 1200)
        .expect("s1's variant");
    assert_eq!(at_variant.1, 1200, "one base");
    // s2 has a block there and is written with its own quality.
    assert_eq!(at_variant.2[1], Some(50));
    // The block s2 had across it is cut either side.
    assert!(out
        .iter()
        .any(|(start, end, _)| *start == 1151 && *end == 1199));
    assert!(out
        .iter()
        .any(|(start, end, _)| *start == 1201 && *end == 1299));
    // And the other sample's variant does the same at 1300.
    let other = out
        .iter()
        .find(|(start, ..)| *start == 1300)
        .expect("s2's variant");
    assert_eq!(other.2[0], Some(40), "s1's own quality there");
}

/// Each sample's own quality, on every record.
#[test]
fn each_sample_keeps_its_own_quality() {
    let text = golden();
    let out = measured(&text, "three-samples");
    assert_eq!(out[0].2, vec![Some(40), Some(30), Some(20)]);
    assert_eq!(out[1].2, vec![Some(40), Some(50), Some(20)]);
    // Nothing is summarised: the three differ on the same record.
    assert_eq!(out[0].2.iter().flatten().count(), 3);
    let unique: std::collections::BTreeSet<i32> = out[0].2.iter().flatten().copied().collect();
    assert_eq!(unique.len(), 3);
}

/// It keeps its column and loses its fields, rather than being padded or dropped.
#[test]
fn a_finished_sample_keeps_its_column() {
    let text = golden();
    let out = measured(&text, "three-samples");
    // s3's input ends at 1150.
    let records = input(&text, "s3");
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].end, 1150);
    // It carries a quality up to there and nothing after.
    assert_eq!(out[1].2[2], Some(20), "the record ending at 1150");
    for (start, _, qualities) in &out {
        if *start > 1150 {
            assert_eq!(qualities[2], None, "at {start}");
        }
    }
    // The column is still there: every record has three entries.
    assert!(out.iter().all(|(_, _, qualities)| qualities.len() == 3));
    // And the writer renders it as a bare no-call.
    assert!(section(&text, "out", "three-samples").contains("\t./.\n"));
}

/// The grid cuts wherever the data does not, and base-pair resolution wins over it.
#[test]
fn the_two_band_arguments() {
    let text = golden();
    // A grid of 100 cuts 1301-1400 at 1400, which no input asked for.
    let plain = measured(&text, "three-samples");
    assert!(plain
        .iter()
        .any(|(start, end, _)| *start == 1301 && *end == 1400));
    let grid = measured(&text, "break-bands-100");
    assert!(grid
        .iter()
        .any(|(start, end, _)| *start == 1301 && *end == 1399));
    assert!(grid
        .iter()
        .any(|(start, end, _)| *start == 1400 && *end == 1400));
    // A finer grid cuts more.
    assert!(measured(&text, "break-bands-50").len() > grid.len());
    // Base-pair resolution is one record per base.
    let bases = measured(&text, "base-pair-resolution");
    assert!(bases.iter().all(|(start, end, _)| start == end));
    // The two together are the base-pair run, not a refusal and not the grid.
    assert_eq!(measured(&text, "both-band-arguments"), bases);
    assert_eq!(
        BandArguments {
            base_pair_resolution: true,
            break_bands_at_multiples_of: 100
        }
        .effective_grid(),
        Some(1)
    );
}

/// Caught by the feature-input check rather than by anything about samples.
#[test]
fn the_same_file_twice_is_refused() {
    let text = golden();
    let (class, message) = refusal(&text, "duplicate-sample");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    let produced = check_inputs(&["a.g.vcf".to_string(), "a.g.vcf".to_string()])
        .expect_err("the same file twice");
    assert_eq!(
        produced,
        CombineError::DuplicateInput {
            path: "a.g.vcf".to_string()
        }
    );
    assert!(
        message.starts_with("Bad input: Feature inputs must be unique"),
        "{message}"
    );
    // Two different files are accepted.
    assert!(check_inputs(&["a.g.vcf".to_string(), "b.g.vcf".to_string()]).is_ok());
}
