//! Conformance for the four BQSR covariates against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/CovariatesDump.java`. A recalibration table is indexed
//! by these keys, so the whole key matrix of every read travels: shape (event type) x (read
//! position) x (covariate), printed cell by cell.
//!
//! # What this suite is for
//!
//!  * **the key matrix is reused between reads and does not leak**, which the golden shows by
//!    carrying the same corpus twice, once with the shared cache BQSR uses and once with a fresh
//!    one per read;
//!  * **the read group is identified by PU and not by ID**, and a missing group has three different
//!    ends: -1 from `keyFromValue`, a null dereference from a read, and `missing key 99` from
//!    `formatKey`;
//!  * **the context is of the read's own bases**, clipped to N at the low-quality tail and
//!    reverse-complemented on the negative strand, with the length packed into the key's low four
//!    bits;
//!  * **the cycle's sign is the low bit**, and all four strand/pair corners are in the corpus;
//!  * **indel cycle keys are -1 within four bases of either end**;
//!  * **a read with bases and no qualities is an exception**, not a covariate value.

use gatk_corpus as corpus;
use gatk_engine::covariates::{
    context_from_key, cycle_from_key, key_from_context, key_from_cycle, ContextCovariate,
    CovariateError, CovariateKind, PerReadCovariateMatrix, QualityScoreCovariate,
    ReadGroupCovariate, RecalibrationArguments, StandardCovariateList, CUSHION_FOR_INDELS,
    MISSING_READ_GROUP_KEY, UNKNOWN_OR_ERROR_CONTEXT_CODE,
};
use gatk_engine::recal_datum::EventType;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/covariates.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn constant(text: &str, name: &str) -> String {
    rows(text, "const")
        .into_iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("the golden has no constant {name}"))[1]
        .to_string()
}

/// The Java exception the reference threw where the port answers this error.
fn java_exception(error: &CovariateError) -> &'static str {
    match error {
        CovariateError::NoReadGroupInHeader => "NullPointerException",
        CovariateError::Clip(_) => "ArrayIndexOutOfBoundsException",
        CovariateError::CycleTooBig { .. } => "UserException",
        CovariateError::MissingKey(_) => "IllegalStateException",
        CovariateError::NegativeContextKey => "GATKException",
        CovariateError::ContextSizeTooBig { .. } => "BadArgumentValue",
        CovariateError::ContextSizeNotPositive { .. } => "CommandLineException",
    }
}

fn event_from(name: &str) -> EventType {
    EventType::from_representation(name).unwrap_or_else(|| panic!("no event type {name}"))
}

/// The keys the port computes for one read, or the error it stops with.
fn matrix_for(
    covariates: &StandardCovariateList,
    header: &SamHeader,
    record: &BamRecord,
    record_indel_values: bool,
) -> Result<PerReadCovariateMatrix, CovariateError> {
    let mut matrix = PerReadCovariateMatrix::new(record.read_bases.len(), covariates.size());
    covariates.populate_per_read_covariate_matrix(
        record,
        header,
        &mut matrix,
        record_indel_values,
    )?;
    Ok(matrix)
}

#[test]
fn the_constants_are_the_references() {
    let text = golden();
    assert_eq!(
        constant(&text, "MISSING_READ_GROUP_KEY"),
        MISSING_READ_GROUP_KEY.to_string()
    );
    assert_eq!(
        constant(&text, "UNKNOWN_OR_ERROR_CONTEXT_CODE"),
        UNKNOWN_OR_ERROR_CONTEXT_CODE.to_string()
    );
    assert_eq!(
        constant(&text, "CUSHION_FOR_INDELS"),
        CUSHION_FOR_INDELS.to_string()
    );
    let arguments = RecalibrationArguments::default();
    assert_eq!(
        constant(&text, "MISMATCHES_CONTEXT_SIZE"),
        arguments.mismatches_context_size.to_string()
    );
    assert_eq!(
        constant(&text, "INDELS_CONTEXT_SIZE"),
        arguments.indels_context_size.to_string()
    );
    assert_eq!(
        constant(&text, "MAXIMUM_CYCLE_VALUE"),
        arguments.maximum_cycle_value.to_string()
    );
    assert_eq!(
        constant(&text, "LOW_QUAL_TAIL"),
        arguments.low_qual_tail.to_string()
    );
    assert_eq!(constant(&text, "size"), "4");
    assert_eq!(constant(&text, "numberOfSpecialCovariates"), "2");
}

/// The read group's identifier is its platform unit, which is what makes the table's keys what they
/// are.
#[test]
fn the_read_group_is_identified_by_its_platform_unit() {
    let text = golden();
    let header = corpus::header(&text);
    let ids = rows(&text, "rgids")[0][0];
    assert_eq!(ReadGroupCovariate::read_group_ids(&header).join(","), ids);

    for row in rows(&text, "rgidentifier") {
        let (id, platform_unit, expected) = (row[0], row[1], row[2]);
        let group = header.read_groups.iter().find(|group| group.id == id);
        match group {
            Some(group) => {
                assert_eq!(ReadGroupCovariate::read_group_identifier(group), expected);
                assert_eq!(
                    group.attributes.get("PU").unwrap_or("null"),
                    platform_unit,
                    "{id}: PU"
                );
            }
            // The dump also builds a group the header does not hold, to show the fallback to ID.
            None => {
                let bare = htsjdk_bam::header::ReadGroup::new(id);
                assert_eq!(platform_unit, "null");
                assert_eq!(ReadGroupCovariate::read_group_identifier(&bare), expected);
            }
        }
    }
}

#[test]
fn the_list_is_the_four_in_order() {
    let text = golden();
    let header = corpus::header(&text);
    let covariates =
        StandardCovariateList::from_header(&RecalibrationArguments::default(), &header).unwrap();

    assert_eq!(covariates.covariate_names(), rows(&text, "names")[0][0]);
    for row in rows(&text, "covariate") {
        let index: usize = row[0].parse().unwrap();
        let kind = covariates.kinds()[index];
        assert_eq!(kind.class_name(), row[1], "index {index}: class name");
        assert_eq!(kind.parsed_name(), row[2], "index {index}: parsed name");
        assert_eq!(
            covariates.maximum_key_value(kind).to_string(),
            row[3],
            "index {index}: maximum key"
        );
        assert_eq!(
            covariates.index_by_class(kind).to_string(),
            row[4],
            "index {index}: index by class"
        );
    }
    for row in rows(&text, "byname") {
        let found = covariates
            .covariate_by_parsed_name(row[0])
            .map(|kind| kind.class_name())
            .unwrap_or("null");
        assert_eq!(found, row[1], "byname {}", row[0]);
    }
    assert_eq!(
        covariates.additional_covariates(),
        [CovariateKind::Context, CovariateKind::Cycle]
    );
    assert_eq!(covariates.number_of_special_covariates(), 2);
}

/// The whole key matrix of every read, cell by cell, for both indel settings.
#[test]
fn every_key_of_every_read_is_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let records = corpus::records(&text);
    let covariates =
        StandardCovariateList::from_header(&RecalibrationArguments::default(), &header).unwrap();

    let mut compared = 0;
    let mut refused = 0;
    for (label, record_indel_values) in [("fresh", true), ("fresh-no-indels", false)] {
        for (index, record) in records.iter().enumerate() {
            let expected: Vec<Vec<&str>> = rows(&text, "matrix")
                .into_iter()
                .filter(|row| row[0] == label && row[1] == index.to_string())
                .collect();
            assert!(!expected.is_empty(), "{label}: no rows for read {index}");

            match matrix_for(&covariates, &header, record, record_indel_values) {
                Ok(matrix) => {
                    for row in expected {
                        assert!(
                            row[2] != "-",
                            "{label}/{index}: the reference stopped here and the port did not: {}",
                            row[4]
                        );
                        let event = event_from(row[2]);
                        let covariate: usize = row[3].parse().unwrap();
                        let ours: Vec<String> = matrix
                            .matrix_for_error_model(event)
                            .iter()
                            .map(|position| position[covariate].to_string())
                            .collect();
                        assert_eq!(
                            ours.join(","),
                            row[4],
                            "{label}/{index}/{event:?}/{covariate}"
                        );
                        compared += 1;
                    }
                }
                Err(error) => {
                    // The reference stopped too, with the same exception, and said so in one row.
                    assert_eq!(
                        expected.len(),
                        1,
                        "{label}/{index}: the port stopped and the reference did not"
                    );
                    let expected_exception = expected[0][4]
                        .strip_prefix("E:")
                        .and_then(|rest| rest.split_once(':'))
                        .map(|(kind, _)| kind)
                        .unwrap_or_else(|| {
                            panic!("{label}/{index}: the port stopped with {error:?} and the reference did not")
                        });
                    assert_eq!(
                        java_exception(&error),
                        expected_exception,
                        "{label}/{index}"
                    );
                    refused += 1;
                }
            }
        }
    }
    println!("covariates: {compared} key rows compared, {refused} reads refused");
    // Three reads whose RG names a group the header does not declare, and one with no qualities,
    // in each of the two indel settings.
    assert_eq!(refused, 8);
}

/// The shared cache and the fresh one agree, which is the measurement that lets this port allocate
/// per read.
#[test]
fn the_shared_key_cache_does_not_leak() {
    let text = golden();
    let mut shared = Vec::new();
    let mut fresh = Vec::new();
    for row in rows(&text, "matrix") {
        match row[0] {
            "shared" => shared.push(row[1..].join("\t")),
            "fresh" => fresh.push(row[1..].join("\t")),
            _ => {}
        }
    }
    assert!(!shared.is_empty());
    assert_eq!(
        shared, fresh,
        "the reference's cache carried something over"
    );
}

/// The bases the context covariate sees: N over the low-quality tail, then reverse-complemented.
#[test]
fn the_clipped_stranded_bases_are_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let records = corpus::records(&text);
    for row in rows(&text, "clipped") {
        let index: usize = row[0].parse().unwrap();
        let low_qual_tail: u8 = row[1].parse().unwrap();
        let arguments = RecalibrationArguments {
            low_qual_tail,
            ..RecalibrationArguments::default()
        };
        let covariate = ContextCovariate::new(&arguments).unwrap();
        match covariate.stranded_clipped_bytes(&records[index], Some(&header)) {
            Ok(bases) => assert_eq!(
                String::from_utf8(bases).unwrap(),
                row[2],
                "read {index} at {low_qual_tail}"
            ),
            Err(_) => assert!(
                row[2].starts_with("E:"),
                "read {index} at {low_qual_tail}: the port stopped and the reference did not"
            ),
        }
    }
}

#[test]
fn the_context_encoding_is_the_reference() {
    let text = golden();
    for row in rows(&text, "context") {
        let (dna, key, back) = (row[0], row[1], row[2]);
        assert_eq!(key_from_context(dna.as_bytes()).to_string(), key, "{dna}");
        match context_from_key(key.parse().unwrap()) {
            Ok(decoded) => assert_eq!(decoded, back, "{dna}: back"),
            Err(_) => assert_eq!(back, "E", "{dna}: back"),
        }
    }
    for row in rows(&text, "contextfromkey") {
        let key: i32 = row[0].parse().unwrap();
        // The dump prints an empty decoding as an empty field, which splits to one element.
        let expected = row.get(1).copied().unwrap_or("");
        assert_eq!(context_from_key(key).unwrap(), expected, "key {key}");
    }
}

#[test]
fn the_cycle_encoding_is_the_reference() {
    let text = golden();
    for row in rows(&text, "cycle") {
        let cycle: i32 = row[0].parse().unwrap();
        let max_cycle: i32 = row[1].parse().unwrap();
        let key = key_from_cycle(cycle, max_cycle).unwrap();
        assert_eq!(key.to_string(), row[2], "cycle {cycle}");
        assert_eq!(cycle_from_key(key).to_string(), row[3], "cycle {cycle}");
    }
    for row in rows(&text, "cyclefromkey") {
        let key: i32 = row[0].parse().unwrap();
        assert_eq!(cycle_from_key(key).to_string(), row[1], "key {key}");
    }
}

/// `formatKey` and `keyFromValue`, which are how a recalibration report is written and read.
#[test]
fn formatting_and_parsing_keys_is_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let arguments = RecalibrationArguments::default();
    let covariates = StandardCovariateList::from_header(&arguments, &header).unwrap();
    let quality = QualityScoreCovariate;

    for row in rows(&text, "format") {
        let (covariate, key, expected) = (row[0], row[1].parse::<i32>().unwrap(), row[2]);
        let ours = match covariate {
            "ReadGroup" => covariates.read_group.format_key(key).unwrap().to_string(),
            "QualityScore" => quality.format_key(key),
            "Context" => covariates
                .context
                .format_key(key)
                .unwrap()
                .unwrap_or_else(|| "null".to_string()),
            "Cycle" => covariates.cycle.format_key(key),
            other => panic!("no covariate {other}"),
        };
        assert_eq!(ours, expected, "{covariate} formatKey({key})");
    }

    for row in rows(&text, "fromvalue") {
        let (covariate, value, expected) = (row[0], row[1], row[2]);
        let ours = match covariate {
            "ReadGroup" => covariates.read_group.key_from_value(value),
            // The three branches of QualityScoreCovariate.keyFromValue all end at the same number,
            // which is why the port has one function and the golden has three rows.
            "QualityScore" => value
                .split_once(':')
                .map(|(_, number)| number.parse::<i32>().unwrap())
                .unwrap(),
            "Context" => covariates.context.key_from_value(value),
            "Cycle" => covariates
                .cycle
                .key_from_value(
                    value
                        .split_once(':')
                        .map(|(_, number)| number.parse::<i32>().unwrap())
                        .unwrap(),
                )
                .unwrap(),
            other => panic!("no covariate {other}"),
        };
        assert_eq!(
            ours.to_string(),
            expected,
            "{covariate} keyFromValue({value})"
        );
    }
}

/// Every argument the covariates refuse, worded as the reference words it.
#[test]
fn the_refusals_are_worded_like_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let message = |what: &str| -> String {
        rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == what)
            .unwrap_or_else(|| panic!("no error row {what}"))[2]
            .to_string()
    };

    for size in [14, 0, -1] {
        let arguments = RecalibrationArguments {
            mismatches_context_size: size,
            ..RecalibrationArguments::default()
        };
        assert_eq!(
            ContextCovariate::new(&arguments).unwrap_err().message(),
            message(&format!("mismatches-context-size@{size}")),
            "mismatches context size {size}"
        );
    }
    for size in [14, 0] {
        let arguments = RecalibrationArguments {
            indels_context_size: size,
            ..RecalibrationArguments::default()
        };
        assert_eq!(
            ContextCovariate::new(&arguments).unwrap_err().message(),
            message(&format!("indels-context-size@{size}")),
            "indels context size {size}"
        );
    }

    let covariates =
        StandardCovariateList::from_header(&RecalibrationArguments::default(), &header).unwrap();
    assert_eq!(
        covariates.read_group.format_key(99).unwrap_err().message(),
        message("readgroup-format-unknown")
    );
    assert_eq!(
        context_from_key(-1).unwrap_err().message(),
        message("context-from-negative")
    );
    for cycle in [501, -501] {
        assert_eq!(
            key_from_cycle(cycle, 500).unwrap_err().message(),
            message(&format!("keyFromCycle@{cycle}")),
            "cycle {cycle}"
        );
    }
    // A maximum of three against a six-base read, which is where a run meets the refusal.
    let tiny = RecalibrationArguments {
        maximum_cycle_value: 3,
        ..RecalibrationArguments::default()
    };
    let tiny_list = StandardCovariateList::from_header(&tiny, &header).unwrap();
    let mut record = corpus::records(&text)[0].clone();
    record.read_bases = b"ACGTAC".to_vec();
    record.base_qualities = vec![40; 6];
    let mut matrix = PerReadCovariateMatrix::new(6, tiny_list.size());
    let error = tiny_list
        .populate_per_read_covariate_matrix(&record, &header, &mut matrix, true)
        .unwrap_err();
    assert_eq!(error.message(), message("cycle-past-maximum"));
}
