//! Conformance for `SomaticClusteringModel` as constructed against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/SomaticClusteringModelDump.java`.
//!
//! # What this suite is for
//!
//!  * **the prior getter mutates the map it reads**, so the same lengths asked in two orders are
//!    asked of two models here, exactly as the dump asks them;
//!  * **only a SNV gets `LOG_ONE_THIRD`**;
//!  * **the mitochondrial defaults are guarded by `==` against the ordinary default**;
//!  * **`record` zeroes the caller's array in place**, printed before and after;
//!  * **two thresholds that look alike do different things**, one counting and one silent;
//!  * **and the initial weights do not sum to one**, which makes
//!    `logLikelihoodGivenSomatic(0, 0)` positive.
//!
//! # Every row, including the five that once waited on the EM iteration
//!
//! The `artifactprior <label>-learned` rows are the dump reading back what `record` accumulated by
//! running `learnAndClearAccumulatedData`. That call is the EM iteration and its quantile
//! initialisation, which the sibling `somatic-clustering-learn` slice ported, so the five rows are
//! reproduced here rather than deferred.
//!
//! # Where the exponential shows
//!
//! Everything on this path that carries `exp` -- `logSumExp` under `logLikelihoodGivenSomatic` and
//! `probabilityOfSequencingError`, and the round of EM behind the `-learned` rows -- is bounded
//! rather than equal, `exp` having no bit-exact transcription under htsjdk-rs decision 0014.
//!
//! Only two rows need the allowance and both need **90 ulps**: `one-alt-learned` and
//! `symbolic-alt-learned`, whose data carry a TLOD of 6 where the sequencing-error probability is
//! large enough to matter. One ulp of `exp` there moves a somatic probability, which moves the
//! variant count, which the prior takes the log of. The three `-learned` rows beside them are
//! bit-identical, as are every `loglike` and `seqerror` row. `EXP_BOUNDED` names the two with the
//! size each needs, so neither a third row nor a larger divergence can appear quietly.

use gatk_corpus as corpus;
use gatk_engine::allele_fraction_cluster::Datum;
use gatk_engine::somatic_clustering_model::{
    indel_length, AlternateAllele, PriorArguments, SomaticClusteringModel,
};
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/somatic_clustering_model.txt.gz"),
    )
}

/// The rows allowed to differ because they carry `exp`, with the allowance each needs.
const EXP_BOUNDED: [(&str, i64); 2] = [("one-alt-learned", 90), ("symbolic-alt-learned", 90)];

fn double(text: &str) -> f64 {
    text.parse().unwrap_or_else(|_| panic!("a double: {text}"))
}

fn snp() -> AlternateAllele {
    AlternateAllele {
        length: 1,
        symbolic: false,
    }
}

/// An alternate of the indel length the label names, against a one-base reference for an insertion
/// and a longer one for a deletion, matching the records the dump builds.
fn record_shape(length: i32) -> (i32, AlternateAllele) {
    if length >= 0 {
        (
            1,
            AlternateAllele {
                length: 1 + length,
                symbolic: false,
            },
        )
    } else {
        (1 - length, snp())
    }
}

/// A model built the way the dump builds its default one.
fn model() -> SomaticClusteringModel {
    SomaticClusteringModel::new(PriorArguments::new(), None)
}

fn prior_line(model: &mut SomaticClusteringModel, label: &str, length: i32) -> String {
    let (reference_length, alternate) = record_shape(length);
    let value = model.log_prior_of_somatic_variant(indel_length(reference_length, alternate));
    format!("prior\t{label}\t{}", java_double_to_string(value))
}

/// Every `prior` and `artifactprior` row, in the dump's order, since the models are stateful.
fn prior_rows() -> Vec<String> {
    let mut rows = Vec::new();
    let mut in_window = model();
    for length in [0, 1, -1, 2, -2, 10, -10] {
        rows.push(prior_line(
            &mut in_window,
            &format!("in-window-{length}"),
            length,
        ));
    }
    let mut ascending = model();
    for length in [11, -11, 50, -50] {
        rows.push(prior_line(
            &mut ascending,
            &format!("ascending-{length}"),
            length,
        ));
    }
    let mut descending = model();
    for length in [-50, 50, -11, 11] {
        rows.push(prior_line(
            &mut descending,
            &format!("descending-{length}"),
            length,
        ));
    }
    let mut repeated = model();
    rows.push(prior_line(&mut repeated, "repeated-first", 11));
    rows.push(prior_line(&mut repeated, "repeated-second", 11));
    rows.push(prior_line(&mut repeated, "repeated-in-window", 0));

    let mut tuned_arguments = PriorArguments::new();
    tuned_arguments.log_snv_prior = -5.0;
    tuned_arguments.log_indel_prior = -4.0;
    tuned_arguments.initial_log_prior_of_variant_versus_artifact = -1.0;
    let mut tuned = SomaticClusteringModel::new(tuned_arguments, None);
    rows.push(prior_line(&mut tuned, "tuned-snv", 0));
    rows.push(prior_line(&mut tuned, "tuned-indel", 3));
    rows.push(format!(
        "artifactprior\ttuned\t{}",
        java_double_to_string(tuned.log_prior_of_variant_versus_artifact())
    ));
    rows.push(format!(
        "artifactprior\tdefault\t{}",
        java_double_to_string(model().log_prior_of_variant_versus_artifact())
    ));

    let mut mito_arguments = PriorArguments::new();
    mito_arguments.mitochondria = true;
    let mut mito = SomaticClusteringModel::new(mito_arguments, None);
    rows.push(prior_line(&mut mito, "mito-snv", 0));
    rows.push(prior_line(&mut mito, "mito-indel", 3));
    let mut overridden_arguments = PriorArguments::new();
    overridden_arguments.mitochondria = true;
    overridden_arguments.log_snv_prior = -5.0;
    let mut overridden = SomaticClusteringModel::new(overridden_arguments, None);
    rows.push(prior_line(&mut overridden, "mito-snv-overridden", 0));
    rows.push(prior_line(&mut overridden, "mito-indel-untouched", 3));
    rows
}

/// The `ads` rows, which are the arrays `record` was handed, before and after.
fn ads_rows() -> Vec<String> {
    let mut rows = Vec::new();
    let mut push = |label: &str,
                    ads: &mut Vec<i32>,
                    odds: &[f64],
                    artifact: &[f64],
                    non_somatic: &[f64],
                    alternates: &[AlternateAllele],
                    reference_length: i32| {
        rows.push(format!("ads\t{label}-before\t{ads:?}"));
        let mut model = SomaticClusteringModel::new(PriorArguments::new(), Some(1000.0));
        if model
            .record(
                ads,
                odds,
                artifact,
                non_somatic,
                alternates,
                reference_length,
            )
            .is_ok()
        {
            rows.push(format!("ads\t{label}-after\t{ads:?}"));
            // What the model accumulated, read back through the prior one round of EM rewrites: a
            // datum that was dropped cannot move it.
            model
                .learn_and_clear_accumulated_data()
                .expect("the initial shapes stay in range");
            rows.push(format!(
                "artifactprior\t{label}-learned\t{}",
                java_double_to_string(model.log_prior_of_variant_versus_artifact())
            ));
        }
    };
    push(
        "one-alt",
        &mut vec![80, 20],
        &[6.0],
        &[0.0],
        &[0.0],
        &[snp()],
        1,
    );
    push(
        "symbolic-alt",
        &mut vec![80, 20, 5],
        &[6.0, 6.0],
        &[0.0, 0.0],
        &[0.0, 0.0],
        &[
            snp(),
            AlternateAllele {
                length: 0,
                symbolic: true,
            },
        ],
        1,
    );
    push(
        "obvious-artifact",
        &mut vec![80, 20],
        &[6.0],
        &[0.95],
        &[0.0],
        &[snp()],
        1,
    );
    push(
        "obvious-non-somatic",
        &mut vec![80, 20],
        &[6.0],
        &[0.0],
        &[0.95],
        &[snp()],
        1,
    );
    push(
        "at-threshold",
        &mut vec![80, 20],
        &[6.0],
        &[0.9],
        &[0.9],
        &[snp()],
        1,
    );
    push(
        "short-array",
        &mut vec![80],
        &[6.0],
        &[0.0],
        &[0.0],
        &[snp()],
        1,
    );
    rows.push(
        "error\tshort-array\tjava.lang.IllegalArgumentException:tumorADs must have one entry per \
         allele including the ref allele"
            .to_string(),
    );
    rows
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let lines: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();

    // The stateful rows are reproduced as whole sequences, since each answer depends on the ones
    // before it.
    let mut ours: Vec<String> = prior_rows();
    let fresh = model();
    for (total, alt) in [
        (10, 0),
        (10, 1),
        (10, 5),
        (10, 10),
        (100, 3),
        (100, 50),
        (1000, 7),
        (0, 0),
    ] {
        let value = fresh
            .log_likelihood_given_somatic(total, alt)
            .expect("finite");
        ours.push(format!(
            "loglike\t{total},{alt}\t{}",
            java_double_to_string(value)
        ));
    }
    // The sequencing-error rows share one model with the likelihood rows above, and the last of them
    // inserts an indel length into its prior map on the way through.
    let mut sequencing = fresh;
    let error_row =
        |model: &mut SomaticClusteringModel, odds: f64, total: i32, alt: i32, length: i32| {
            let datum = Datum::new(odds, 0.0, 0.0, alt, total, length);
            let value = model
                .probability_of_sequencing_error(&datum)
                .expect("finite");
            format!(
                "seqerror\t{}-{total},{alt}-{length}\t{}",
                java_double_to_string(odds),
                java_double_to_string(value)
            )
        };
    for odds in [0.0, 5.0, 20.0, -5.0] {
        for (total, alt) in [(10, 5), (100, 3)] {
            ours.push(error_row(&mut sequencing, odds, total, alt, 0));
        }
    }
    ours.push(error_row(&mut sequencing, 5.0, 10, 5, 3));
    ours.push(error_row(&mut sequencing, 5.0, 10, 5, -3));
    ours.push(error_row(&mut sequencing, 5.0, 10, 5, 40));

    ours.extend(ads_rows());

    for (label, reference_length, alternate) in [
        ("snp", 1, snp()),
        (
            "insertion",
            1,
            AlternateAllele {
                length: 4,
                symbolic: false,
            },
        ),
        ("deletion", 4, snp()),
        (
            "symbolic",
            1,
            AlternateAllele {
                length: 0,
                symbolic: true,
            },
        ),
    ] {
        ours.push(format!(
            "indellength\t{label}\t{}",
            indel_length(reference_length, alternate)
        ));
    }

    let expected = &lines;
    assert_eq!(lines.len(), 66, "the golden's row count");

    for (mine, theirs) in ours.iter().zip(expected) {
        if mine == theirs {
            continue;
        }
        let label = theirs.split('\t').nth(1).unwrap_or_default();
        let allowed = EXP_BOUNDED
            .iter()
            .find(|(name, _)| *name == label)
            .map(|(_, ulps)| *ulps)
            .unwrap_or_else(|| panic!("{label}: {mine} against {theirs}"));
        let value: f64 = double(mine.rsplit('\t').next().expect("a value"));
        let reference: f64 = double(theirs.rsplit('\t').next().expect("a value"));
        let ulps = ((value.to_bits() as i64) - (reference.to_bits() as i64)).abs();
        assert!(
            ulps <= allowed,
            "{label}: {mine} against {theirs}, {ulps} ulps"
        );
    }
    assert_eq!(ours.len(), expected.len(), "every row is accounted for");
}
