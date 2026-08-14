//! Conformance for `GermlineFilter.germlineProbability` against GATK 4.6.2.0, compared as the
//! probability of every call in the dump's four sweeps.
//!
//! Golden from `tools/readfilter-conformance/GermlineProbabilityDump.java`.
//!
//! # What this suite is for
//!
//!  * **the population frequency decides both ends**, `0.0` at zero and `1.0` at one;
//!  * **the answer is the first entry of the normalisation**, not the second;
//!  * **the somatic prior enters both sides**, so a prior of one answers `0.0`;
//!  * **the hom-alt hypothesis is switched off by a negative infinity**;
//!  * **and an infinite normal log odds is refused where a negative infinity is zero**.
//!
//! Thirty-one of the thirty-four values are bit-identical; three sit one ulp away, for the same reason
//! the sibling `mutect-engine-arithmetic` suite has one: every entry of the normalisation goes
//! through `exp`, whose bit-exact transcription is withdrawn under htsjdk-rs decision 0014. The three
//! are named in `ONE_ULP_APART`, and the comparison fails if a fourth appears **or** if one of those
//! three starts agreeing.

use gatk_corpus as corpus;
use gatk_engine::germline_filter::germline_probability;
use gatk_engine::tsv_table::java_double_to_string;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/germline_probability.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn expected(text: &str, label: &str) -> String {
    rows(text, "germline")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no probability {label}"))[1]
        .to_string()
}

const NORMAL_LOG_ODDS: f64 = 5.0;
const HET_VS_SOMATIC: f64 = 1.0;
const HOM_ALT_VS_SOMATIC: f64 = 0.0;
const POPULATION_AF: f64 = 0.001;
const LOG_PRIOR_SOMATIC: f64 = -13.0;

/// The dump's label for a double, which is `Double.toString` after string concatenation.
fn label(prefix: &str, value: f64) -> String {
    format!("{prefix}{}", java_double_to_string(value))
}

fn ours(normal: f64, het: f64, hom_alt: f64, af: f64, prior: f64) -> String {
    germline_probability(normal, het, hom_alt, af, prior)
        .map(java_double_to_string)
        .unwrap_or_else(|_| "refused".to_string())
}

/// The inputs on which the port is one ulp from the reference rather than equal to it.
///
/// Every entry of the normalisation goes through `exp`, whose bit-exact transcription is withdrawn
/// under htsjdk-rs decision 0014, so a divergence in the last place is expected on some inputs and
/// is bounded rather than absent. These are the ones in this suite, named so a third cannot appear
/// without failing the test.
const ONE_ULP_APART: [&str; 3] = ["af-1.0E-8", "normal--10.0", "normal--50.0"];

/// Compare, allowing exactly the named inputs to sit one ulp away.
fn same(ours: &str, theirs: &str, label: &str) {
    if ours == theirs {
        assert!(
            !ONE_ULP_APART.contains(&label),
            "{label} agrees now: take it out of ONE_ULP_APART"
        );
        return;
    }
    assert!(
        ONE_ULP_APART.contains(&label),
        "{label}: {ours} against {theirs}"
    );
    let mine: f64 = ours.parse().expect("a double");
    let reference: f64 = theirs.parse().expect("a double");
    let ulps = ((mine.to_bits() as i64) - (reference.to_bits() as i64)).abs();
    assert!(ulps <= 1, "{label}: {ours} against {theirs}, {ulps} ulps");
}

#[test]
fn every_probability_matches_the_golden() {
    let text = golden();
    same(
        &ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        ),
        &expected(&text, "baseline"),
        "baseline",
    );

    for af in [0.0, 1.0e-8, 1.0e-4, 0.001, 0.1, 0.5, 0.9, 0.999, 1.0] {
        let key = label("af-", af);
        same(
            &ours(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                af,
                LOG_PRIOR_SOMATIC,
            ),
            &expected(&text, &key),
            &key,
        );
    }

    for odds in [-50.0, -10.0, -1.0, 0.0, 1.0, 10.0, 50.0, f64::NEG_INFINITY] {
        let key = label("normal-", odds);
        same(
            &ours(
                odds,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC,
            ),
            &expected(&text, &key),
            &key,
        );
    }

    for het in [-10.0, -1.0, 0.0, 1.0, 10.0] {
        let key = label("het-", het);
        same(
            &ours(
                NORMAL_LOG_ODDS,
                het,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC,
            ),
            &expected(&text, &key),
            &key,
        );
    }

    for prior in [-50.0, -13.0, -1.0, -1.0e-9, 0.0, f64::NEG_INFINITY] {
        let key = label("prior-", prior);
        same(
            &ours(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                prior,
            ),
            &expected(&text, &key),
            &key,
        );
    }
}

#[test]
fn the_population_frequency_decides_both_ends() {
    let text = golden();
    assert_eq!(expected(&text, "af-0.0"), "0.0");
    assert_eq!(expected(&text, "af-1.0"), "1.0");
    // And it rises in between.
    let values: Vec<f64> = ["af-1.0E-8", "af-1.0E-4", "af-0.1", "af-0.9"]
        .iter()
        .map(|label| expected(&text, label).parse().expect("a double"))
        .collect();
    assert!(
        values.windows(2).all(|pair| pair[0] < pair[1]),
        "{values:?}"
    );
}

#[test]
fn the_somatic_prior_enters_both_sides() {
    let text = golden();
    // Certainly somatic, whatever the odds say.
    assert_eq!(expected(&text, "prior-0.0"), "0.0");
    // Certainly not somatic.
    assert_eq!(expected(&text, "prior--Infinity"), "1.0");
    assert_eq!(
        ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            POPULATION_AF,
            0.0
        ),
        "0.0"
    );
}

#[test]
fn the_hom_alt_hypothesis_is_switched_off_by_a_value() {
    let text = golden();
    same(
        &ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            f64::NEG_INFINITY,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        ),
        &expected(&text, "homalt-off"),
        "homalt-off",
    );
    same(
        &ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            f64::NEG_INFINITY,
            0.5,
            LOG_PRIOR_SOMATIC,
        ),
        &expected(&text, "homalt-off-common"),
        "homalt-off-common",
    );
    // Switching it off can only lower the germline probability, at either frequency.
    let rare_on: f64 = expected(&text, "baseline").parse().expect("a double");
    let rare_off: f64 = expected(&text, "homalt-off").parse().expect("a double");
    let common_on: f64 = expected(&text, "af-0.5").parse().expect("a double");
    let common_off: f64 = expected(&text, "homalt-off-common")
        .parse()
        .expect("a double");
    assert!(rare_off < rare_on);
    assert!(common_off < common_on);
}

#[test]
fn an_infinite_normal_log_odds_is_refused() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "normal-Infinity")
        .expect("a refusal");
    let (class, message) = row[1].split_once(':').expect("class and message");
    let error = germline_probability(
        f64::INFINITY,
        HET_VS_SOMATIC,
        HOM_ALT_VS_SOMATIC,
        POPULATION_AF,
        LOG_PRIOR_SOMATIC,
    )
    .expect_err("not finite");
    assert_eq!(error.class(), class);
    assert_eq!(error.message(), message);
    // While the other infinity is simply zero.
    assert_eq!(expected(&text, "normal--Infinity"), "0.0");
}

#[test]
fn both_corners_normalise_rather_than_refusing() {
    let text = golden();
    assert_eq!(
        ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            0.0,
            0.0
        ),
        expected(&text, "af-zero-and-certain-somatic")
    );
    assert_eq!(
        ours(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            1.0,
            f64::NEG_INFINITY
        ),
        expected(&text, "af-one-and-impossible-somatic")
    );
    assert_eq!(expected(&text, "af-zero-and-certain-somatic"), "0.0");
    assert_eq!(expected(&text, "af-one-and-impossible-somatic"), "1.0");
}
