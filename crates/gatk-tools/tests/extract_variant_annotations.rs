//! Conformance for `ExtractVariantAnnotations` against GATK 4.6.2.0, compared as the rows of both
//! matrices and the sites-only VCF of every run.
//!
//! Golden from `tools/readfilter-conformance/ExtractVariantAnnotationsDump.java`.
//!
//! Reading a VCF and writing HDF5 are not measured or ported, and neither is the random stream:
//! the reservoir's indices are derived from what the golden kept rather than generated.
//!
//! # What this suite is for
//!
//!  * **the mode deciding which variant types are kept**;
//!  * **a label coming only from a tag whose value is the string `true`**;
//!  * **`snp` being reserved, because the matrix carries one of its own**;
//!  * **an unlabelled row going to a SEPARATE file, or nowhere**;
//!  * **the reservoir not being in genomic order**;
//!  * **the columns being sorted by name**;
//!  * **every absence being the same NaN**;
//!  * **the default matching strategy being the loosest one**;
//!  * **only the minimal representation reconciling a padded allele**;
//!  * **an allele-specific annotation switching the whole run**;
//!  * **and the alternate array being flat.**

use gatk_corpus as corpus;
use gatk_tools::extract_variant_annotations::{
    alternate_alleles, attributes, check_resource_labels, context_type, datasets,
    decode_annotation, extract, is_allele_in_list, label_column, passes_filters, reference_alleles,
    snp_column, sorted_annotation_names, sorted_labels, variant_type, Arguments, ContextType,
    MatchingStrategy, Record, Reservoir, Resource, Row, VariantType, ALLELES_ALT_PATH,
    ALLELES_REF_PATH, DEFAULT_MATCHING_STRATEGY, SNP_LABEL,
};
use std::collections::BTreeSet;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/extract_variant_annotations.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
        .map(unescape)
}

fn section(text: &str, kind: &str, name: &str) -> String {
    field(text, kind, name).unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

/// The records of one VCF the golden carries.
fn records(text: &str, name: &str) -> Vec<Record> {
    section(text, "vcf", name)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Record {
                contig: columns[0].to_string(),
                start: columns[1].parse().expect("a position"),
                reference: columns[3].to_string(),
                alternates: columns[4].split(',').map(str::to_string).collect(),
                filters: if columns[6] == "PASS" || columns[6] == "." {
                    Vec::new()
                } else {
                    columns[6].split(';').map(str::to_string).collect()
                },
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

/// The intervals one matrix carries, as `contig:start-end`.
fn intervals(text: &str, name: &str) -> Vec<String> {
    section(text, "intervals", name)
        .split(',')
        .map(str::to_string)
        .collect()
}

/// The numbers one matrix carries, row by row, rendered the way `Double.toString` renders them.
fn rows(text: &str, name: &str) -> Vec<Vec<String>> {
    section(text, "rows", name)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// One label's column, as the golden wrote it.
fn labels(text: &str, name: &str, key: &str) -> Vec<bool> {
    section(text, "labels", name)
        .lines()
        .find_map(|line| line.strip_prefix(&format!("{key}: ")))
        .unwrap_or_else(|| panic!("the golden carries labels/{name}/{key}"))
        .split(',')
        .map(|value| value == "true")
        .collect()
}

/// `Double.toString`, as far as this fixture's numbers need it: every value here is a NaN or a
/// double with one decimal.
fn java_double(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    format!("{value:?}")
}

fn rendered(rows: &[Row]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|row| row.annotations.iter().map(|v| java_double(*v)).collect())
        .collect()
}

fn placed(rows: &[Row]) -> Vec<String> {
    rows.iter()
        .map(|row| format!("{}:{}-{}", row.contig, row.start, row.end))
        .collect()
}

fn resource(name: &str, tags: &[(&str, &str)], records: Vec<Record>) -> Resource {
    Resource {
        name: name.to_string(),
        tags: attributes(tags),
        records,
        // The fixture's resources carry genotypes and are polymorphic in them.
        has_genotypes: true,
        is_polymorphic: true,
    }
}

fn training(text: &str) -> Resource {
    resource("train", &[("training", "true")], records(text, "training"))
}

fn calibration(text: &str) -> Resource {
    resource(
        "cal",
        &[("calibration", "true")],
        records(text, "calibration"),
    )
}

fn modes(kinds: &[VariantType]) -> BTreeSet<VariantType> {
    kinds.iter().copied().collect()
}

/// A run with no reservoir: the labelled matrix alone.
fn run(text: &str, resources: &[Resource], arguments: Arguments, requested: &[&str]) -> Vec<Row> {
    let names: Vec<String> = requested.iter().map(|n| n.to_string()).collect();
    let mut never = |_: usize| unreachable!("no reservoir was asked for");
    extract(
        &records(text, "input"),
        resources,
        &arguments,
        &names,
        0,
        &mut never,
    )
    .labeled
}

/// Every run without a reservoir, compared as its intervals, its columns and its numbers.
#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let train = training(&text);
    let cal = calibration(&text);
    let cases: Vec<(&str, Vec<Resource>, Arguments, Vec<&str>)> = vec![
        (
            "snp",
            vec![train.clone()],
            Arguments {
                modes: modes(&[VariantType::Snp]),
                ..Arguments::default()
            },
            vec!["QD", "MQ"],
        ),
        (
            "indel",
            vec![train.clone()],
            Arguments {
                modes: modes(&[VariantType::Indel]),
                ..Arguments::default()
            },
            vec!["QD", "MQ"],
        ),
        (
            "both-modes",
            vec![train.clone()],
            Arguments::default(),
            vec!["QD", "MQ"],
        ),
        (
            "two-labels",
            vec![train.clone(), cal.clone()],
            Arguments::default(),
            // Asked for out of order, which the columns do not keep.
            vec!["MQ", "FS", "QD"],
        ),
        (
            "match-start_position",
            vec![train.clone()],
            Arguments {
                strategy: MatchingStrategy::StartPosition,
                ..Arguments::default()
            },
            vec!["QD"],
        ),
        (
            "match-start_position_and_given_representation",
            vec![train.clone()],
            Arguments {
                strategy: MatchingStrategy::StartPositionAndGivenRepresentation,
                ..Arguments::default()
            },
            vec!["QD"],
        ),
        (
            "match-start_position_and_minimal_representation",
            vec![train.clone()],
            Arguments {
                strategy: MatchingStrategy::StartPositionAndMinimalRepresentation,
                ..Arguments::default()
            },
            vec!["QD"],
        ),
        (
            "omit-alleles",
            vec![train.clone()],
            Arguments {
                modes: modes(&[VariantType::Snp]),
                ..Arguments::default()
            },
            vec!["QD"],
        ),
    ];
    let mut compared = 0;
    for (label, resources, arguments, requested) in cases {
        let produced = run(&text, &resources, arguments, &requested);
        assert_eq!(placed(&produced), intervals(&text, label), "{label}");
        assert_eq!(rendered(&produced), rows(&text, label), "{label}");
        assert_eq!(
            sorted_annotation_names(&requested.iter().map(|n| n.to_string()).collect::<Vec<_>>())
                .join(","),
            section(&text, "names", label),
            "{label}"
        );
        assert_eq!(
            snp_column(&produced),
            labels(&text, label, SNP_LABEL),
            "{label}"
        );
        assert_eq!(
            label_column(&produced, "training"),
            labels(&text, label, "training"),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the runs the port reproduces");
}

/// A record of the wrong type is dropped from every file, not written with no label.
#[test]
fn the_mode_decides_which_types_are_kept() {
    let text = golden();
    let input = records(&text, "input");
    // The fixture: SNPs at 1000, 2000, 5000, 6000 and 7000, indels at 3000 and 4000.
    let kinds: Vec<(i32, VariantType)> = input
        .iter()
        .map(|record| (record.start, variant_type(record).expect("a type")))
        .collect();
    assert_eq!(
        kinds,
        vec![
            (1000, VariantType::Snp),
            (2000, VariantType::Snp),
            (3000, VariantType::Indel),
            (4000, VariantType::Indel),
            (5000, VariantType::Snp),
            (6000, VariantType::Snp),
            (7000, VariantType::Snp),
        ]
    );
    // The SNP run and the INDEL run partition the both-modes run between them.
    let snp = intervals(&text, "snp");
    let indel = intervals(&text, "indel");
    let both = intervals(&text, "both-modes");
    assert_eq!(snp.len() + indel.len(), both.len());
    for place in snp.iter().chain(indel.iter()) {
        assert!(both.contains(place), "{place}");
    }
    // The multiallelic record is a SNP, because both its alternates are.
    assert_eq!(context_type(&input[6]), ContextType::Snp);
    assert_eq!(input[6].alternates, vec!["C", "G"]);
}

/// A tag written `training=false` labels nothing, so the run extracts nothing.
#[test]
fn a_label_comes_only_from_a_true_tag() {
    let text = golden();
    let false_tag = resource(
        "train",
        &[("training", "false")],
        records(&text, "training"),
    );
    assert!(false_tag.labels().is_empty());
    let produced = run(
        &text,
        &[false_tag],
        Arguments {
            modes: modes(&[VariantType::Snp]),
            ..Arguments::default()
        },
        &["QD"],
    );
    assert!(produced.is_empty());
    // Which is the run that wrote no matrix at all, though it still wrote a VCF.
    assert_eq!(
        field(&text, "none", "false-tag").as_deref(),
        Some("no annotations hdf5")
    );
    assert!(field(&text, "out", "false-tag").is_some());
    // The true tag does label.
    assert_eq!(
        training(&text).labels().into_iter().collect::<Vec<_>>(),
        vec!["training".to_string()]
    );
}

/// Because the matrix carries a `snp` column of its own for every run.
#[test]
fn the_snp_label_is_reserved() {
    let text = golden();
    let (class, message) = refusal(&text, "reserved-label");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    let reserved: BTreeSet<String> = [SNP_LABEL.to_string()].into_iter().collect();
    assert_eq!(
        check_resource_labels(&reserved).expect_err("reserved"),
        message
    );
    // A resource label that is not the reserved one is fine.
    let ordinary: BTreeSet<String> = ["training".to_string()].into_iter().collect();
    assert!(check_resource_labels(&ordinary).is_ok());
    // And the column is there whatever the resources say, sorting before the others.
    assert_eq!(sorted_labels(&ordinary), vec!["snp", "training"]);
    assert_eq!(
        labels(&text, "snp", SNP_LABEL),
        vec![true, true],
        "both rows of the SNP run"
    );
    assert_eq!(labels(&text, "indel", SNP_LABEL), vec![false]);
}

/// Never into the labelled one, and only when a reservoir was asked for at all.
#[test]
fn an_unlabeled_row_goes_to_a_separate_file() {
    let text = golden();
    let mut never = |_: usize| 0;
    let produced = extract(
        &records(&text, "input"),
        &[training(&text)],
        &Arguments::default(),
        &["QD".to_string(), "MQ".to_string()],
        10,
        &mut never,
    );
    assert_eq!(placed(&produced.labeled), intervals(&text, "unlabeled-all"));
    assert_eq!(rendered(&produced.labeled), rows(&text, "unlabeled-all"));
    let reservoir = produced.unlabeled.as_ref().expect("a reservoir");
    assert_eq!(
        placed(&reservoir.rows),
        intervals(&text, "unlabeled-all.unlabeled")
    );
    assert_eq!(
        rendered(&reservoir.rows),
        rows(&text, "unlabeled-all.unlabeled")
    );
    // No row is in both files.
    for place in placed(&reservoir.rows) {
        assert!(!placed(&produced.labeled).contains(&place), "{place}");
    }
    // Every unlabelled row's `training` column is false, and the `snp` column still varies.
    assert_eq!(
        label_column(&reservoir.rows, "training"),
        labels(&text, "unlabeled-all.unlabeled", "training")
    );
    assert_eq!(
        snp_column(&reservoir.rows),
        labels(&text, "unlabeled-all.unlabeled", SNP_LABEL)
    );
    // A run with no reservoir writes no such file, and drops those rows entirely.
    assert_eq!(
        field(&text, "none", "both-modes.unlabeled").as_deref(),
        Some("no annotations hdf5")
    );
    assert!(produced.writes_labeled_matrix());
    assert!(produced.writes_unlabeled_matrix());
}

/// The reservoir fills in order and is then overwritten in place.
#[test]
fn the_reservoir_is_not_in_genomic_order() {
    let text = golden();
    // Four unlabelled records survive when every filter is ignored: 2000, 4000, 5000, 6000.
    let arguments = Arguments {
        ignore_all_filters: true,
        ..Arguments::default()
    };
    let candidates = {
        let mut all = |_: usize| 0;
        let produced = extract(
            &records(&text, "input"),
            &[training(&text)],
            &arguments,
            &["QD".to_string()],
            10,
            &mut all,
        );
        placed(&produced.unlabeled.expect("a reservoir").rows)
    };
    assert_eq!(
        candidates,
        vec![
            "chr1:2000-2000",
            "chr1:4000-4002",
            "chr1:5000-5000",
            "chr1:6000-6000"
        ]
    );
    // Seeds 0, 1 and 100 all keep 2000 and 6000, which Algorithm R reaches by sending the third
    // record and then the fourth to slot one.
    for label in [
        "unlabeled-two-seed-0",
        "unlabeled-two-seed-1",
        "unlabeled-two-seed-100",
    ] {
        let mut indices = [1usize, 1].into_iter();
        let mut random = |_: usize| indices.next().expect("two draws");
        let produced = extract(
            &records(&text, "input"),
            &[training(&text)],
            &arguments,
            &["QD".to_string()],
            2,
            &mut random,
        );
        let reservoir = produced.unlabeled.expect("a reservoir");
        assert_eq!(
            placed(&reservoir.rows),
            intervals(&text, &format!("{label}.unlabeled")),
            "{label}"
        );
    }
    // Seed 42 keeps a different pair, and keeps it out of genomic order.
    let mut indices = [1usize, 0].into_iter();
    let mut random = |_: usize| indices.next().expect("two draws");
    let produced = extract(
        &records(&text, "input"),
        &[training(&text)],
        &arguments,
        &["QD".to_string()],
        2,
        &mut random,
    );
    let reservoir = produced.unlabeled.expect("a reservoir");
    assert_eq!(
        placed(&reservoir.rows),
        intervals(&text, "unlabeled-two-seed-42.unlabeled")
    );
    assert_eq!(
        placed(&reservoir.rows),
        vec!["chr1:6000-6000", "chr1:5000-5000"],
        "the later position first"
    );
    // The draw is consulted only once the reservoir is full.
    let mut small = Reservoir::new(2);
    small.offer(vec![produced.labeled[0].clone()], 99);
    small.offer(vec![produced.labeled[1].clone()], 99);
    assert_eq!(small.rows.len(), 2);
    assert_eq!(small.seen, 2);
}

/// Sorted by name, whatever order they were asked for in, and deduplicated.
#[test]
fn the_columns_are_sorted_by_name() {
    let text = golden();
    let requested: Vec<String> = ["MQ", "FS", "QD"].iter().map(|n| n.to_string()).collect();
    assert_eq!(sorted_annotation_names(&requested), vec!["FS", "MQ", "QD"]);
    assert_eq!(section(&text, "names", "two-labels"), "FS,MQ,QD");
    // The SNP run asked for QD and then MQ and got them the other way round.
    assert_eq!(section(&text, "names", "snp"), "MQ,QD");
    let repeated: Vec<String> = ["QD", "QD", "MQ"].iter().map(|n| n.to_string()).collect();
    assert_eq!(sorted_annotation_names(&repeated), vec!["MQ", "QD"]);
}

/// Missing, unparseable and infinite are one absence, not three.
#[test]
fn every_absence_is_the_same_nan() {
    let text = golden();
    let input = records(&text, "input");
    // The record at 6000 carries no QD at all and an MQ of `Infinity`.
    let sixth = input.iter().find(|r| r.start == 6000).expect("a record");
    assert_eq!(sixth.attribute("QD"), None);
    assert_eq!(sixth.attribute("MQ"), Some("Infinity"));
    assert!(decode_annotation(sixth, None, "QD", false).is_nan());
    assert!(decode_annotation(sixth, None, "MQ", false).is_nan());
    // Which is the row the golden wrote as two NaNs side by side.
    let unlabeled = rows(&text, "unlabeled-all.unlabeled");
    assert_eq!(unlabeled.last().expect("a row"), &vec!["NaN", "NaN"]);
    // A value that does not parse is the same absence.
    let odd = Record {
        attributes: attributes(&[("QD", "high")]),
        ..sixth.clone()
    };
    assert!(decode_annotation(&odd, None, "QD", false).is_nan());
    // A negative infinity too, so no infinity of either sign reaches the matrix.
    let negative = Record {
        attributes: attributes(&[("QD", "-Infinity")]),
        ..sixth.clone()
    };
    assert!(decode_annotation(&negative, None, "QD", false).is_nan());
    // And an ordinary value survives.
    let first = input.iter().find(|r| r.start == 1000).expect("a record");
    assert_eq!(decode_annotation(first, None, "QD", false), 1.5);
}

/// The loosest one: the start position and the variant class, and nothing about the alleles.
#[test]
fn the_default_matching_strategy_is_the_loosest() {
    let text = golden();
    assert_eq!(DEFAULT_MATCHING_STRATEGY, MatchingStrategy::StartPosition);
    // The default run and the explicit START_POSITION run kept the same rows.
    assert_eq!(
        intervals(&text, "both-modes"),
        intervals(&text, "match-start_position")
    );
    // At 7000 the input is `A>C,G` and the resource is `A>T`: no allele in common, and the
    // loosest strategy labels it anyway.
    assert!(intervals(&text, "match-start_position").contains(&"chr1:7000-7000".to_string()));
    assert!(
        !intervals(&text, "match-start_position_and_given_representation")
            .contains(&"chr1:7000-7000".to_string())
    );
    assert!(
        !intervals(&text, "match-start_position_and_minimal_representation")
            .contains(&"chr1:7000-7000".to_string())
    );
}

/// The same insertion written two ways, which only one strategy recognises.
#[test]
fn only_the_minimal_representation_reconciles_a_padded_allele() {
    let text = golden();
    let resource_records = records(&text, "training");
    // The resource writes the insertion at 3000 with a padding base.
    let padded = resource_records
        .iter()
        .find(|r| r.start == 3000)
        .expect("a record");
    assert_eq!(padded.reference, "CG");
    assert_eq!(padded.alternates, vec!["CATG"]);
    assert!(is_allele_in_list(
        "C",
        "CAT",
        &padded.reference,
        &padded.alternates
    ));
    // As written they share no alternate, so the middle strategy does not match.
    assert!(!padded.alternates.contains(&"CAT".to_string()));
    let given = intervals(&text, "match-start_position_and_given_representation");
    let minimal = intervals(&text, "match-start_position_and_minimal_representation");
    assert!(!given.contains(&"chr1:3000-3000".to_string()));
    assert!(minimal.contains(&"chr1:3000-3000".to_string()));
    // A genuinely different event is not reconciled, however it is padded.
    assert!(!is_allele_in_list("C", "CAT", "CA", &["CAAT".to_string()]));
}

/// One row per alternate, and the strictest strategy whatever was asked for.
#[test]
fn an_allele_specific_annotation_switches_the_whole_run() {
    let text = golden();
    let arguments = Arguments {
        allele_specific: true,
        ..Arguments::default()
    };
    let mut all = |_: usize| 0;
    let produced = extract(
        &records(&text, "input"),
        &[training(&text)],
        &arguments,
        &["AS_QD".to_string()],
        10,
        &mut all,
    );
    assert_eq!(
        placed(&produced.labeled),
        intervals(&text, "allele-specific")
    );
    assert_eq!(rendered(&produced.labeled), rows(&text, "allele-specific"));
    let reservoir = produced.unlabeled.as_ref().expect("a reservoir");
    assert_eq!(
        placed(&reservoir.rows),
        intervals(&text, "allele-specific.unlabeled")
    );
    assert_eq!(
        rendered(&reservoir.rows),
        rows(&text, "allele-specific.unlabeled")
    );
    // The multiallelic record became TWO rows at the same position, one per alternate.
    let at_seven: Vec<&Row> = reservoir.rows.iter().filter(|r| r.start == 7000).collect();
    assert_eq!(at_seven.len(), 2);
    assert_eq!(at_seven[0].alternates, vec!["C"]);
    assert_eq!(at_seven[1].alternates, vec!["G"]);
    assert_eq!(at_seven[0].annotations, vec![7.5]);
    assert_eq!(at_seven[1].annotations, vec![8.5]);
    // And it fell OUT of the labelled matrix, because the run now uses the strictest strategy.
    assert!(!intervals(&text, "allele-specific").contains(&"chr1:7000-7000".to_string()));
    assert!(intervals(&text, "both-modes").contains(&"chr1:7000-7000".to_string()));
}

/// A multiallelic record contributes both alternates, so the two arrays cannot be read together.
#[test]
fn the_alternate_array_is_flat() {
    let text = golden();
    let produced = run(
        &text,
        &[training(&text)],
        Arguments {
            modes: modes(&[VariantType::Snp]),
            ..Arguments::default()
        },
        &["QD", "MQ"],
    );
    assert_eq!(
        reference_alleles(&produced).join(","),
        section(&text, "alleles", "snp/alleles/ref")
    );
    assert_eq!(
        alternate_alleles(&produced).join(","),
        section(&text, "alleles", "snp/alleles/alt")
    );
    assert_eq!(reference_alleles(&produced).len(), 2);
    assert_eq!(alternate_alleles(&produced).len(), 3, "one row gave two");
    // --omit-alleles-in-hdf5 takes both datasets out and leaves everything else.
    assert!(!datasets(true)[ALLELES_REF_PATH]);
    assert!(!datasets(true)[ALLELES_ALT_PATH]);
    assert!(datasets(false)[ALLELES_REF_PATH]);
    assert!(datasets(false)[ALLELES_ALT_PATH]);
    assert_eq!(
        field(&text, "none", "omit-alleles/alleles/ref").as_deref(),
        Some("no alleles")
    );
    assert_eq!(
        field(&text, "none", "omit-alleles/alleles/alt").as_deref(),
        Some("no alleles")
    );
    // The rest of the matrix is unchanged by it.
    assert_eq!(intervals(&text, "omit-alleles"), intervals(&text, "snp"));
}

/// One filter argument names a filter and the other takes them all.
#[test]
fn a_filtered_record_is_dropped() {
    let text = golden();
    let input = records(&text, "input");
    let filtered = input.iter().find(|r| r.start == 5000).expect("a record");
    assert_eq!(filtered.filters, vec!["LOW"]);
    assert!(!passes_filters(filtered, &Arguments::default()));
    let named = Arguments {
        ignored_filters: ["LOW".to_string()].into_iter().collect(),
        ..Arguments::default()
    };
    assert!(passes_filters(filtered, &named));
    let all = Arguments {
        ignore_all_filters: true,
        ..Arguments::default()
    };
    assert!(passes_filters(filtered, &all));
    // Both arguments let it into the reservoir, and the golden agrees they let in the same rows.
    assert_eq!(
        intervals(&text, "ignore-filter.unlabeled"),
        intervals(&text, "ignore-all-filters.unlabeled")
    );
    assert!(intervals(&text, "ignore-filter.unlabeled").contains(&"chr1:5000-5000".to_string()));
    // The default run left it out.
    assert!(!intervals(&text, "unlabeled-all.unlabeled").contains(&"chr1:5000-5000".to_string()));
    // A record with two filters needs both named, not one.
    let two = Record {
        filters: vec!["LOW".to_string(), "BAD".to_string()],
        ..filtered.clone()
    };
    assert!(!passes_filters(&two, &named));
    assert!(passes_filters(&two, &all));
}

/// A run that keeps nothing writes no matrix and still writes its VCF.
#[test]
fn a_run_that_extracts_nothing_still_writes_a_vcf() {
    let text = golden();
    let produced = run(&text, &[], Arguments::default(), &["QD"]);
    assert!(produced.is_empty());
    assert_eq!(
        field(&text, "none", "extracts-nothing").as_deref(),
        Some("no annotations hdf5")
    );
    // The VCF is there, and it is a header alone.
    let vcf = section(&text, "out", "extracts-nothing");
    assert_eq!(vcf.lines().count(), 1);
    assert!(vcf.starts_with("#CHROM"));
}
