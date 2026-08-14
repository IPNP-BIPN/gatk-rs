//! `GermlineFilter.germlineProbability`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.mutect.filtering.GermlineFilter` (GATK 4.6.2.0).
//!
//! The probability that an allele Mutect called somatic is really germline, from five numbers: the
//! normal's log odds, the two log odds of germline against somatic, the population allele frequency
//! and the log prior that the site is somatic.
//!
//! # The population frequency is three priors at once
//!
//! ```java
//! final double logPriorGermlineHet = Math.log(2*populationAF*(1-populationAF));
//! final double logPriorGermlineHomAlt = Math.log( MathUtils.square(populationAF));
//! final double logPriorNotGermline = Math.log(MathUtils.square(1 - populationAF));
//! ```
//!
//! One number decides all three, so a frequency of zero makes both germline hypotheses impossible
//! and the answer `0.0`, while a frequency of one makes the somatic hypothesis impossible and the
//! answer `1.0`. Neither is a special case in the code: both fall out of `log(0)`.
//!
//! # The answer is the first entry, not the second
//!
//! ```java
//! return NaturalLogUtils.normalizeLog(new double[] {logProbGermline, logProbSomatic}, false, true)[0];
//! ```
//!
//! [`crate::mutect_engine::posterior_probability_of_error`] has the same shape and returns index
//! **1**. These two functions answer opposite questions, and taking the wrong index turns a germline
//! filter into a somatic one that never fires.
//!
//! # The somatic prior enters twice
//!
//! Once as itself on the somatic side, and once through `log1mexp` on **both** germline sides. A
//! prior of one therefore answers `0.0` whatever the odds say, and a prior of negative infinity
//! answers `1.0`.

use crate::natural_log_utils::{
    log1mexp, log_sum_exp, normalize_from_log_to_linear_space, NonFiniteSum,
};

/// `germlineProbability`.
///
/// `log_odds_of_germline_hom_alt_vs_somatic` is negative infinity when the caller judged the allele
/// fraction too low for a germline hom alt: the hypothesis is switched off by a value rather than by
/// a flag.
pub fn germline_probability(
    normal_log_odds: f64,
    log_odds_of_germline_het_vs_somatic: f64,
    log_odds_of_germline_hom_alt_vs_somatic: f64,
    population_af: f64,
    log_prior_somatic: f64,
) -> Result<f64, NonFiniteSum> {
    let log_prior_not_somatic = log1mexp(log_prior_somatic);
    let log_prior_germline_het = (2.0 * population_af * (1.0 - population_af)).ln();
    let log_prior_germline_hom_alt = (population_af * population_af).ln();
    let log_prior_not_germline = ((1.0 - population_af) * (1.0 - population_af)).ln();

    // Unnormalized, and the normal's odds reach both germline hypotheses and neither somatic one.
    let log_prob_germline_het = log_prior_germline_het
        + log_odds_of_germline_het_vs_somatic
        + normal_log_odds
        + log_prior_not_somatic;
    let log_prob_germline_hom_alt = log_prior_germline_hom_alt
        + log_odds_of_germline_hom_alt_vs_somatic
        + normal_log_odds
        + log_prior_not_somatic;
    let log_prob_germline = log_sum_exp(&[log_prob_germline_het, log_prob_germline_hom_alt])?;
    let log_prob_somatic = log_prior_not_germline + log_prior_somatic;

    let normalized = normalize_from_log_to_linear_space(&[log_prob_germline, log_prob_somatic])?;
    // The FIRST entry: germline.
    Ok(normalized[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    const NORMAL_LOG_ODDS: f64 = 5.0;
    const HET_VS_SOMATIC: f64 = 1.0;
    const HOM_ALT_VS_SOMATIC: f64 = 0.0;
    const POPULATION_AF: f64 = 0.001;
    const LOG_PRIOR_SOMATIC: f64 = -13.0;

    fn baseline(af: f64, prior: f64) -> f64 {
        germline_probability(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            af,
            prior,
        )
        .expect("finite")
    }

    #[test]
    fn the_population_frequency_decides_both_ends() {
        // Nothing can be germline.
        assert_eq!(baseline(0.0, LOG_PRIOR_SOMATIC), 0.0);
        // Nothing can be somatic.
        assert_eq!(baseline(1.0, LOG_PRIOR_SOMATIC), 1.0);
        // And in between it rises with the frequency.
        assert!(baseline(1.0e-8, LOG_PRIOR_SOMATIC) < baseline(1.0e-4, LOG_PRIOR_SOMATIC));
        assert!(baseline(1.0e-4, LOG_PRIOR_SOMATIC) < baseline(0.5, LOG_PRIOR_SOMATIC));
    }

    #[test]
    fn the_somatic_prior_enters_twice() {
        // A prior of one: certainly somatic, whatever the odds.
        assert_eq!(baseline(POPULATION_AF, 0.0), 0.0);
        // A prior of zero: certainly not.
        assert_eq!(baseline(POPULATION_AF, f64::NEG_INFINITY), 1.0);
        // And it moves the answer everywhere in between.
        assert!(baseline(POPULATION_AF, -1.0) < baseline(POPULATION_AF, -13.0));
    }

    #[test]
    fn the_normal_odds_move_it_monotonically_and_an_infinity_is_refused() {
        let mut previous = 0.0;
        for odds in [-50.0, -10.0, -1.0, 0.0, 1.0, 10.0, 50.0] {
            let value = germline_probability(
                odds,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC,
            )
            .expect("finite");
            assert!(value >= previous, "{odds}");
            previous = value;
        }
        // Negative infinity is simply zero.
        assert_eq!(
            germline_probability(
                f64::NEG_INFINITY,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                POPULATION_AF,
                LOG_PRIOR_SOMATIC
            )
            .expect("finite"),
            0.0
        );
        // Positive infinity is a refusal from the normaliser.
        let error = germline_probability(
            f64::INFINITY,
            HET_VS_SOMATIC,
            HOM_ALT_VS_SOMATIC,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        )
        .expect_err("not finite");
        assert_eq!(
            error.message(),
            "logValues must be non-infinite and non-NAN"
        );
    }

    #[test]
    fn the_hom_alt_hypothesis_is_switched_off_by_a_value() {
        let on = baseline(POPULATION_AF, LOG_PRIOR_SOMATIC);
        let off = germline_probability(
            NORMAL_LOG_ODDS,
            HET_VS_SOMATIC,
            f64::NEG_INFINITY,
            POPULATION_AF,
            LOG_PRIOR_SOMATIC,
        )
        .expect("finite");
        // At a rare allele the hom alt hypothesis is worth almost nothing.
        assert!(off < on);
        assert!(on - off < 1.0e-6);
    }

    #[test]
    fn both_corners_normalise_rather_than_refusing() {
        // No germline hypothesis and a certain somatic one.
        assert_eq!(
            germline_probability(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                0.0,
                0.0
            )
            .expect("one side survives"),
            0.0
        );
        // No somatic hypothesis and an impossible somatic prior.
        assert_eq!(
            germline_probability(
                NORMAL_LOG_ODDS,
                HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC,
                1.0,
                f64::NEG_INFINITY
            )
            .expect("one side survives"),
            1.0
        );
    }
}
