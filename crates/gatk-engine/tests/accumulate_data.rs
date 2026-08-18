//! Conformance for `accumulateData` and the pass schedule against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/AccumulateDataDump.java`.
//!
//! # What this suite is for
//!
//!  * **four passes, not two**, three of which accumulate and only two of which learn parameters;
//!  * **a record whose only alternate is `<NON_REF>` is skipped**, as is one with no alternate;
//!  * **a record with no `TLOD` is a refusal**, not a skip;
//!  * **and the three accumulators do not move together**: a triallelic record whose alternates are
//!    both obvious artifacts accumulates no clustering data, two probabilities and two obvious
//!    artifacts.
//!
//! # What the golden pins here, and what it does not
//!
//! The `accumulated` rows are **counts**, not values: the clustering model's data size, the
//! threshold calculator's probability count, and the obvious-artifact count. Those counts name which
//! branch each alternate took, and this test supplies probabilities on that branch. What it does not
//! do is recompute the probabilities the eighteen real filters produced for these records, which
//! needs the engine assembled around them; that is its own slice. So this compares the control flow
//! exactly and the arithmetic not at all, which is what the golden holds.

use gatk_corpus as corpus;
use gatk_engine::accumulate_data::{
    accumulate_data, action_after_pass, action_for_pass, AccumulationAllele, AfterPassAction,
    PassAction, NUMBER_OF_LEARNING_PASSES, NUMBER_OF_PASSES,
};
use gatk_engine::somatic_clustering_model::{
    AlternateAllele, PriorArguments, SomaticClusteringModel,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/accumulate_data.txt.gz"),
    )
}

fn rows() -> Vec<(String, String, String)> {
    golden()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut fields = line.splitn(3, '\t');
            (
                fields.next().expect("a kind").to_string(),
                fields.next().expect("a label").to_string(),
                fields.next().expect("a payload").to_string(),
            )
        })
        .collect()
}

fn substitution() -> AccumulationAllele {
    AccumulationAllele {
        allele: AlternateAllele {
            length: 1,
            symbolic: false,
        },
        non_ref: false,
    }
}

fn non_ref() -> AccumulationAllele {
    AccumulationAllele {
        allele: AlternateAllele {
            length: 0,
            symbolic: true,
        },
        non_ref: true,
    }
}

/// `<DEL>`: symbolic, and **not** `<NON_REF>`, so the guard does not skip it.
fn symbolic_deletion() -> AccumulationAllele {
    AccumulationAllele {
        allele: AlternateAllele {
            length: 0,
            symbolic: true,
        },
        non_ref: false,
    }
}

/// One record's arguments. The probabilities are on the branch the golden's counts name.
struct Case {
    alternates: Vec<AccumulationAllele>,
    depths: Vec<i32>,
    tumor_log_odds: Option<Vec<f64>>,
    artifact: Vec<f64>,
    non_somatic: Vec<f64>,
    combined: Vec<f64>,
}

fn case(label: &str) -> Case {
    match label {
        "real-alternate" | "three-records" => Case {
            alternates: vec![substitution()],
            depths: vec![80, 20],
            tumor_log_odds: Some(vec![46.05]),
            artifact: vec![0.0],
            non_somatic: vec![0.0],
            combined: vec![0.0],
        },
        // Both alternates were obvious artifacts: no data, two probabilities, two counted.
        "two-alternates" => Case {
            alternates: vec![substitution(), substitution()],
            depths: vec![80, 20, 5],
            tumor_log_odds: Some(vec![46.05, 13.8]),
            artifact: vec![0.95, 0.95],
            non_somatic: vec![0.0, 0.0],
            combined: vec![0.95, 0.95],
        },
        "only-non-ref" => Case {
            alternates: vec![non_ref()],
            depths: vec![80, 20],
            tumor_log_odds: Some(vec![46.05]),
            artifact: vec![0.0],
            non_somatic: vec![0.0],
            combined: vec![0.0],
        },
        "alternate-and-non-ref" => Case {
            alternates: vec![substitution(), non_ref()],
            depths: vec![80, 20, 0],
            tumor_log_odds: Some(vec![46.05, 0.0]),
            artifact: vec![0.0, 0.0],
            non_somatic: vec![0.0, 0.0],
            // `ErrorProbabilities` removed the symbolic allele's entry before the threshold saw it.
            combined: vec![0.0],
        },
        // Not skipped by the guard, and the symbolic removal left nothing to accumulate.
        "symbolic-not-non-ref" => Case {
            alternates: vec![symbolic_deletion()],
            depths: vec![80, 20],
            tumor_log_odds: Some(vec![46.05]),
            artifact: vec![0.0],
            non_somatic: vec![0.0],
            combined: Vec::new(),
        },
        "reference-only" => Case {
            alternates: Vec::new(),
            depths: vec![80],
            tumor_log_odds: Some(Vec::new()),
            artifact: Vec::new(),
            non_somatic: Vec::new(),
            combined: Vec::new(),
        },
        "no-tlod" => Case {
            alternates: vec![substitution()],
            depths: vec![80, 20],
            tumor_log_odds: None,
            artifact: vec![0.0],
            non_somatic: vec![0.0],
            combined: vec![0.0],
        },
        other => panic!("no case named {other}"),
    }
}

/// The three counts one label produces, as the dump prints them.
fn counts(label: &str) -> Result<String, String> {
    let mut model = SomaticClusteringModel::new(PriorArguments::new(), None);
    let mut probabilities = Vec::new();
    // `three-records` puts the same record through one engine three times.
    let times = if label == "three-records" { 3 } else { 1 };
    for _ in 0..times {
        let case = case(label);
        let mut depths = case.depths.clone();
        accumulate_data(
            &mut model,
            &mut probabilities,
            &case.alternates,
            &mut depths,
            case.tumor_log_odds.as_deref(),
            &case.artifact,
            &case.non_somatic,
            &case.combined,
            1,
        )
        .map_err(|error| format!("{}:{}", error.class(), error.message()))?;
    }
    Ok(format!(
        "{},{},{}",
        model.accumulated(),
        probabilities.len(),
        model.obvious_artifact_count()
    ))
}

#[test]
fn every_row_matches_the_golden() {
    let rows = rows();
    assert_eq!(rows.len(), 12, "the golden's row count");
    for (kind, label, payload) in &rows {
        match kind.as_str() {
            "passes" => {
                let ours = match label.as_str() {
                    "numberOfPasses" => NUMBER_OF_PASSES,
                    "learningPasses" => NUMBER_OF_LEARNING_PASSES,
                    // The passes that accumulate: 0 through NUMBER_OF_LEARNING_PASSES inclusive.
                    "accumulatingPasses" => (0..NUMBER_OF_PASSES)
                        .filter(|pass| action_for_pass(*pass) == Some(PassAction::Accumulate))
                        .count() as i32,
                    // The pass that applies and writes.
                    "applyingPass" => (0..NUMBER_OF_PASSES)
                        .find(|pass| action_for_pass(*pass) == Some(PassAction::ApplyAndWrite))
                        .expect("one pass applies"),
                    other => panic!("no schedule value named {other}"),
                };
                assert_eq!(ours.to_string(), *payload, "passes {label}");
            }
            "accumulated" => assert_eq!(counts(label).expect("accumulated"), *payload, "{label}"),
            "error" => assert_eq!(counts(label).expect_err("refused"), *payload, "{label}"),
            other => panic!("no row kind {other}"),
        }
    }
    // The schedule the golden's four numbers summarise, in full.
    assert_eq!(
        action_after_pass(2),
        Some(AfterPassAction::LearnThresholdOnly)
    );
}
