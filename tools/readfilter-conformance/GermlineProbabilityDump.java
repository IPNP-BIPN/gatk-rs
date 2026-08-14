/*
 * GermlineFilter.germlineProbability, taken from the reference.
 *
 * The probability that an allele Mutect called somatic is really germline, from five numbers: the
 * normal's log odds, the two log odds of germline against somatic, the population allele frequency,
 * and the log prior that the site is somatic. Six behaviours this is built to catch.
 *
 *   - THE POPULATION FREQUENCY IS THREE PRIORS AT ONCE, `log(2f(1-f))` for a het, `log(f^2)` for a
 *     hom alt and `log((1-f)^2)` for neither, so a frequency of zero makes both germline hypotheses
 *     impossible and a frequency of one makes the somatic hypothesis impossible;
 *   - THE ANSWER IS THE FIRST ENTRY of the normalisation, the germline one, where the engine's own
 *     posterior returns the second: the two functions look alike and answer opposite questions;
 *   - THE SOMATIC PRIOR ENTERS TWICE, once as itself on the somatic side and once through
 *     `log1mexp` on both germline sides, so a prior of one drives the germline probability to zero
 *     whatever the odds say;
 *   - THE HOM ALT HYPOTHESIS IS SWITCHED OFF BY A NEGATIVE INFINITY rather than by a flag, which is
 *     what the caller passes when the allele fraction is too low for a germline hom alt;
 *   - A FREQUENCY OF ZERO AGAINST A SOMATIC PRIOR OF ONE LEAVES NOTHING TO NORMALISE, and what
 *     comes back is measured rather than assumed;
 *   - AND THE NORMAL'S LOG ODDS ARE ADDED TO BOTH GERMLINE HYPOTHESES, so they move the answer
 *     without ever touching the somatic side.
 *
 * Output:
 *
 *     germline\t<label>\t<the probability>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GermlineProbabilityDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.GermlineFilter;

public class GermlineProbabilityDump {

    /** The baseline every sweep varies one axis of. */
    static final double NORMAL_LOG_ODDS = 5.0;
    static final double HET_VS_SOMATIC = 1.0;
    static final double HOM_ALT_VS_SOMATIC = 0.0;
    static final double POPULATION_AF = 0.001;
    static final double LOG_PRIOR_SOMATIC = -13.0;

    public static void main(final String[] args) {
        System.out.println("# GermlineProbabilityDump: five numbers in, one probability out");

        probability("baseline", NORMAL_LOG_ODDS, HET_VS_SOMATIC, HOM_ALT_VS_SOMATIC, POPULATION_AF,
                LOG_PRIOR_SOMATIC);

        // The population frequency, from impossible to certain.
        for (final double af : new double[] {0.0, 1.0e-8, 1.0e-4, 0.001, 0.1, 0.5, 0.9, 0.999, 1.0}) {
            probability("af-" + af, NORMAL_LOG_ODDS, HET_VS_SOMATIC, HOM_ALT_VS_SOMATIC, af,
                    LOG_PRIOR_SOMATIC);
        }

        // The normal's log odds, which reach both germline hypotheses and neither somatic one.
        for (final double odds : new double[] {-50.0, -10.0, -1.0, 0.0, 1.0, 10.0, 50.0,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY}) {
            probability("normal-" + odds, odds, HET_VS_SOMATIC, HOM_ALT_VS_SOMATIC, POPULATION_AF,
                    LOG_PRIOR_SOMATIC);
        }

        // The het odds.
        for (final double het : new double[] {-10.0, -1.0, 0.0, 1.0, 10.0}) {
            probability("het-" + het, NORMAL_LOG_ODDS, het, HOM_ALT_VS_SOMATIC, POPULATION_AF,
                    LOG_PRIOR_SOMATIC);
        }

        // The hom alt hypothesis, switched off the way the caller switches it off.
        probability("homalt-off", NORMAL_LOG_ODDS, HET_VS_SOMATIC, Double.NEGATIVE_INFINITY,
                POPULATION_AF, LOG_PRIOR_SOMATIC);
        probability("homalt-off-common", NORMAL_LOG_ODDS, HET_VS_SOMATIC, Double.NEGATIVE_INFINITY,
                0.5, LOG_PRIOR_SOMATIC);

        // The somatic prior, from impossible to certain.
        for (final double prior : new double[] {-50.0, -13.0, -1.0, -1.0e-9, 0.0,
                Double.NEGATIVE_INFINITY}) {
            probability("prior-" + prior, NORMAL_LOG_ODDS, HET_VS_SOMATIC, HOM_ALT_VS_SOMATIC,
                    POPULATION_AF, prior);
        }

        // Nothing left to normalise: no germline hypothesis and no somatic one either.
        probability("af-zero-and-certain-somatic", NORMAL_LOG_ODDS, HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC, 0.0, 0.0);
        probability("af-one-and-impossible-somatic", NORMAL_LOG_ODDS, HET_VS_SOMATIC,
                HOM_ALT_VS_SOMATIC, 1.0, Double.NEGATIVE_INFINITY);
    }

    static void probability(final String label, final double normalLogOdds,
                            final double logOddsOfGermlineHetVsSomatic,
                            final double logOddsOfGermlineHomAltVsSomatic,
                            final double populationAF, final double logPriorSomatic) {
        try {
            System.out.printf("germline\t%s\t%s%n", label, Double.toString(
                    GermlineFilter.germlineProbability(normalLogOdds, logOddsOfGermlineHetVsSomatic,
                            logOddsOfGermlineHomAltVsSomatic, populationAF, logPriorSomatic)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
