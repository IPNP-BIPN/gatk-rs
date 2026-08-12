//! Conformance for `BQSRReadTransformer` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/BqsrTransformerDump.java`. This is the thing ApplyBQSR
//! is, and every other suite in the recalibration machinery is an input to it.
//!
//! # What this suite is for
//!
//!  * **the estimate is `y3 + y4 - y2`**, and a null datum contributes nothing rather than a zero
//!    delta;
//!  * **the datum cache makes it order-dependent**, measured directly on one datum asked twice;
//!  * **the rounding is `fastRound`**, so the double below a half rounds twice;
//!  * **a quality below `preserveQLessThan` is left alone entirely**, not even quantized;
//!  * **`--allow-missing-read-group` covers a narrower case than its name suggests**;
//!  * **the static mapping sorts its argument list** and rounds in probability space.

use std::cell::RefCell;
use std::rc::Rc;

use gatk_corpus as corpus;
use gatk_engine::bqsr_transformer::{
    bounded_integer_qual, construct_static_quantized_mapping,
    hierarchical_bayesian_quality_estimate, ApplyBqsrArguments, BqsrReadTransformer,
};
use gatk_engine::covariates::{RecalibrationArguments, StandardCovariateList};
use gatk_engine::qual_quantizer::{QualQuantizer, QuantizationInfo, MIN_USABLE_Q_SCORE};
use gatk_engine::recal_datum::RecalDatum;
use gatk_engine::recalibration_tables::{RecalibrationTables, SharedDatum};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::TagValue;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/bqsr_transformer.txt.gz"),
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

fn datum(observations: i64, mismatches: f64, quality: i8) -> SharedDatum {
    Rc::new(RefCell::new(
        RecalDatum::new(observations, mismatches, quality).unwrap(),
    ))
}

fn bits(value: f64) -> String {
    format!("{:x}", value.to_bits())
}

/// The same table the dump built: a datum in every place the transformer looks.
fn build_tables(covariates: &StandardCovariateList) -> RecalibrationTables {
    let mut tables = RecalibrationTables::new(covariates).unwrap();
    tables
        .read_group_table_mut()
        .put(datum(100_000, 1000.0, 30), &[0, 0])
        .unwrap();
    for quality in [2i32, 5, 6, 20, 30, 40] {
        tables
            .quality_score_table_mut()
            .put(datum(10_000, 50.0, quality as i8), &[0, quality, 0])
            .unwrap();
        for key in 0..260i32 {
            tables.all_tables[2]
                .put(
                    datum(1000, 5.0 + (key % 7) as f64, quality as i8),
                    &[0, quality, key, 0],
                )
                .unwrap();
        }
        for key in 0..24i32 {
            tables.all_tables[3]
                .put(
                    datum(1000, 3.0 + (key % 5) as f64, quality as i8),
                    &[0, quality, key, 0],
                )
                .unwrap();
        }
    }
    tables
}

/// `new QuantizationInfo(tables, QUANTIZING_LEVELS)`: the empirical quality histogram of the quality
/// score table, quantized to sixteen levels.
fn quantization_info(tables: &RecalibrationTables) -> QuantizationInfo {
    let mut histogram = vec![0i64; 94];
    for (_, datum) in tables.quality_score_table().all_leaves() {
        let empirical = gatk_engine::math_utils::fast_round(datum.borrow_mut().empirical_quality());
        histogram[empirical as usize] += datum.borrow().num_observations();
    }
    let quantizer = QualQuantizer::new(&histogram, 16, MIN_USABLE_Q_SCORE).unwrap();
    QuantizationInfo::new(quantizer.original_to_quantized_map, histogram)
}

/// The four reads the dump built, taken from the corpus it printed alongside them.
fn read(text: &str, name: &str) -> BamRecord {
    corpus::records(text)
        .into_iter()
        .find(|record| record.read_name == name)
        .unwrap_or_else(|| panic!("the golden has no read {name}"))
}

/// One run of the transformer, over the same reads and arguments the dump used.
fn run(
    label: &str,
    header: &SamHeader,
    reads: &[BamRecord],
    arguments: &ApplyBqsrArguments,
    read_groups: &[String],
) -> Vec<(String, String)> {
    let covariates =
        StandardCovariateList::new(&RecalibrationArguments::default(), read_groups).unwrap();
    let mut tables = build_tables(&covariates);
    let mut quantization = quantization_info(&tables);
    let mut transformer = BqsrReadTransformer::new(
        header,
        &mut tables,
        &mut quantization,
        &covariates,
        arguments,
    )
    .unwrap();

    let _ = label;
    reads
        .iter()
        .map(|record| {
            let outcome = match transformer.apply(record) {
                Ok(out) => out
                    .base_qualities
                    .iter()
                    .map(|quality| quality.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                Err(error) => format!("E:GATKException:{}", error.message()),
            };
            (record.read_name.clone(), outcome)
        })
        .collect()
}

#[test]
fn the_static_mapping_is_the_reference() {
    let text = golden();
    for row in rows(&text, "static") {
        let label = row[0];
        let quals: Vec<i32> = label
            .split('@')
            .next()
            .unwrap()
            .trim_matches(|c| c == '[' || c == ']')
            .split(", ")
            .filter(|piece| !piece.is_empty())
            .map(|piece| piece.parse().unwrap())
            .collect();
        let round_down = label.ends_with("@down");
        let mut mutable = quals.clone();
        let mapping = construct_static_quantized_mapping(&mut mutable, round_down);
        // The reference's mapping is a `byte[]`, so it prints signed and the identity map wraps.
        let ours: Vec<String> = mapping
            .iter()
            .map(|value| (*value as i8).to_string())
            .collect();
        assert_eq!(ours.join(","), row[1], "{label}");

        // The call sorted the caller's own list.
        let expected_sorted = rows(&text, "sorted")
            .into_iter()
            .find(|sorted| sorted[0] == label)
            .unwrap()[1]
            .to_string();
        let sorted_text = format!(
            "[{}]",
            mutable
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        assert_eq!(sorted_text, expected_sorted, "{label}: sorted in place");
    }
}

#[test]
fn the_rounding_is_the_reference() {
    let text = golden();
    for row in rows(&text, "round") {
        let value: f64 = row[0].parse().unwrap();
        assert_eq!(
            gatk_engine::math_utils::fast_round(value).to_string(),
            row[1],
            "fastRound({value})"
        );
        assert_eq!(
            bounded_integer_qual(value).to_string(),
            row[2],
            "boundQual({value})"
        );
    }
}

#[test]
fn the_hierarchical_estimate_is_the_reference() {
    let text = golden();
    let expected = |label: &str| -> String {
        rows(&text, "estimate")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no estimate row {label}"))[1]
            .to_string()
    };
    let prior = 25.0;

    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            None,
            None,
            &[None, None]
        )),
        expected("all-null")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            Some(&datum(10000, 100.0, 30)),
            None,
            &[None, None]
        )),
        expected("read-group-only")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            Some(&datum(10000, 100.0, 30)),
            Some(&datum(5000, 20.0, 30)),
            &[None, None]
        )),
        expected("two-covariates")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            Some(&datum(10000, 100.0, 30)),
            Some(&datum(5000, 20.0, 30)),
            &[Some(datum(1000, 1.0, 30)), None]
        )),
        expected("one-special")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            Some(&datum(10000, 100.0, 30)),
            Some(&datum(5000, 20.0, 30)),
            &[Some(datum(1000, 1.0, 30)), Some(datum(2000, 200.0, 30))]
        )),
        expected("two-specials")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            prior,
            None,
            None,
            &[Some(datum(1000, 1.0, 30)), None]
        )),
        expected("special-only")
    );

    // The cache, which is what makes a run order-dependent.
    let shared = datum(1000, 1.0, 30);
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            25.0,
            None,
            Some(&shared),
            &[]
        )),
        expected("cached-first-25")
    );
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            45.0,
            None,
            Some(&shared),
            &[]
        )),
        expected("cached-then-45")
    );
    let fresh = datum(1000, 1.0, 30);
    assert_eq!(
        bits(hierarchical_bayesian_quality_estimate(
            45.0,
            None,
            Some(&fresh),
            &[]
        )),
        expected("fresh-45")
    );
}

/// Every run of the transformer, over the same reads and arguments the dump used.
#[test]
fn every_recalibrated_read_is_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let all_groups: Vec<String> = vec![
        "unit-rg1".to_string(),
        "unit-rg2".to_string(),
        "unit-rg3".to_string(),
    ];
    let one_group = vec!["unit-rg1".to_string()];

    let first = read(&text, "first");
    let second = read(&text, "second");
    let low = read(&text, "low");
    let unknown = read(&text, "unknown");

    let expected = |label: &str| -> Vec<(String, String)> {
        rows(&text, "apply")
            .into_iter()
            .filter(|row| row[0] == label)
            .map(|row| (row[1].to_string(), row[2].to_string()))
            .collect()
    };

    let cases: Vec<(&str, Vec<BamRecord>, ApplyBqsrArguments, &[String])> = vec![
        (
            "in-order",
            vec![first.clone(), second.clone(), low.clone()],
            ApplyBqsrArguments::default(),
            &all_groups,
        ),
        (
            "reversed",
            vec![second.clone(), first.clone(), low.clone()],
            ApplyBqsrArguments::default(),
            &all_groups,
        ),
        (
            "quantized-4",
            vec![first.clone(), second.clone()],
            ApplyBqsrArguments {
                quantization_levels: 4,
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "no-quantization",
            vec![first.clone(), second.clone()],
            ApplyBqsrArguments::default(),
            &all_groups,
        ),
        (
            "static-quals",
            vec![first.clone(), second.clone()],
            ApplyBqsrArguments {
                static_quantization_quals: vec![10, 20, 30, 40],
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "global-prior",
            vec![first.clone(), second.clone()],
            ApplyBqsrArguments {
                global_qscore_prior: 20.0,
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "preserve-31",
            vec![first.clone(), second.clone(), low.clone()],
            ApplyBqsrArguments {
                preserve_qscores_less_than: 31,
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "emit-original",
            vec![first.clone()],
            ApplyBqsrArguments {
                emit_original_quals: true,
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "no-datum-for-read-group",
            vec![unknown.clone()],
            ApplyBqsrArguments::default(),
            &all_groups,
        ),
        (
            "no-datum-for-read-group-allowed",
            vec![unknown.clone()],
            ApplyBqsrArguments {
                allow_missing_read_groups: true,
                ..ApplyBqsrArguments::default()
            },
            &all_groups,
        ),
        (
            "covariate-missing-read-group",
            vec![first.clone(), unknown.clone()],
            ApplyBqsrArguments::default(),
            &one_group,
        ),
        (
            "covariate-missing-read-group-allowed",
            vec![first.clone(), unknown.clone()],
            ApplyBqsrArguments {
                allow_missing_read_groups: true,
                ..ApplyBqsrArguments::default()
            },
            &one_group,
        ),
    ];

    let mut compared = 0;
    for (label, reads, arguments, groups) in cases {
        let ours = run(label, &header, &reads, &arguments, groups);
        let theirs = expected(label);
        assert_eq!(ours.len(), theirs.len(), "{label}: read count");
        for ((name, ours), (their_name, theirs)) in ours.iter().zip(&theirs) {
            assert_eq!(name, their_name, "{label}: read order");
            assert_eq!(ours, theirs, "{label}/{name}");
            compared += 1;
        }
    }
    println!("bqsr-transformer: {compared} reads recalibrated");
}

/// The `OQ` tag, which is written only when it is not already there.
#[test]
fn the_original_qualities_tag_is_the_reference() {
    let text = golden();
    let header = corpus::header(&text);
    let groups: Vec<String> = vec![
        "unit-rg1".to_string(),
        "unit-rg2".to_string(),
        "unit-rg3".to_string(),
    ];
    let first = read(&text, "first");

    let covariates =
        StandardCovariateList::new(&RecalibrationArguments::default(), &groups).unwrap();
    let mut tables = build_tables(&covariates);
    let mut quantization = quantization_info(&tables);
    let arguments = ApplyBqsrArguments {
        emit_original_quals: true,
        ..ApplyBqsrArguments::default()
    };
    let mut transformer = BqsrReadTransformer::new(
        &header,
        &mut tables,
        &mut quantization,
        &covariates,
        &arguments,
    )
    .unwrap();
    let out = transformer.apply(&first).unwrap();
    let tag = out
        .tags
        .iter()
        .find(|(tag, _)| tag.name() == *b"OQ")
        .map(|(_, value)| match value {
            TagValue::Str(text) => text.clone(),
            other => format!("{other:?}"),
        })
        .unwrap();
    let expected = rows(&text, "oqtag")
        .into_iter()
        .find(|row| row[0] == "emit-original")
        .unwrap()[2]
        .to_string();
    assert_eq!(tag, expected);
}
