/*
 * Mutect2FilteringEngine's static arithmetic, taken from the reference.
 *
 * The three public statics every filter in the engine goes through, measured on their own before the
 * engine that calls them exists. Five behaviours this is built to catch.
 *
 *   - THE POSTERIOR IS A TWO-ENTRY NORMALISATION, `[logOdds + logPrior, log1mexp(logPrior)]` put
 *     through `normalizeFromLogToLinearSpace`, of which the SECOND entry is returned. It is
 *     therefore the probability of ERROR, and it falls as the odds of being real rise;
 *   - A PRIOR OF ONE MAKES THE ERROR IMPOSSIBLE, `log1mexp(0)` being negative infinity, and a prior
 *     of zero makes it certain whatever the odds say;
 *   - INFINITE ODDS DO NOT ALWAYS GIVE A CLEAN ANSWER: the sum inside the normalisation can be
 *     `-inf + -inf` or `inf - inf`, and what comes back is measured rather than assumed;
 *   - `roundFinitePrecisionErrors` IS A CLAMP, `max(min(p, 1), 0)`, which is not a rounding at all:
 *     it moves a probability of 1.0000000001 to 1 and leaves NaN exactly where it is, because
 *     Java's `Math.min` and `Math.max` propagate it;
 *   - AND `getTumorLogOdds` CONVERTS FROM LOG10 TO NATURAL LOG, which is why every threshold
 *     compared against it is in units the annotation is not written in. A missing annotation is
 *     null rather than an empty array.
 *
 * Output:
 *
 *     posterior\t<log odds>,<log prior>\t<the probability of error>
 *     round\t<label>\t<the clamped value>
 *     tumorlogodds\t<label>\t<the converted array, or null>
 *
 * Usage: MutectEngineArithmeticDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.util.Arrays;
import java.util.List;

public class MutectEngineArithmeticDump {

    public static void main(final String[] args) {
        System.out.println("# MutectEngineArithmeticDump: two entries normalised, one clamp, one conversion");

        final double[] logOdds = {-10.0, -1.0, 0.0, 1.0, 10.0, 100.0,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY};
        // A prior of one is log 0; a prior of zero is log of negative infinity.
        final double[] logPriors = {0.0, -0.001, -1.0, -10.0, -100.0, Double.NEGATIVE_INFINITY};

        for (final double odds : logOdds) {
            for (final double prior : logPriors) {
                posterior(odds, prior);
            }
        }

        // The clamp, on both sides and on the values that are not numbers.
        for (final double probability : new double[] {-0.5, -1.0e-12, 0.0, 0.5, 1.0,
                1.0 + 1.0e-10, 2.0, Double.NaN, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY}) {
            System.out.printf("round\t%s\t%s%n", Double.toString(probability),
                    Double.toString(Mutect2FilteringEngine.roundFinitePrecisionErrors(probability)));
        }

        // The conversion, and the missing annotation.
        tumorLogOdds("one", List.of(6.0));
        tumorLogOdds("two", List.of(6.0, 4.0));
        tumorLogOdds("zero-and-negative", List.of(0.0, -3.0));
        tumorLogOdds("absent", null);

        System.out.printf("round\tEPSILON\t%s%n", Double.toString(Mutect2FilteringEngine.EPSILON));
        System.out.printf("round\tMIN_REPORTABLE_ERROR_PROBABILITY\t%s%n",
                Double.toString(Mutect2FilteringEngine.MIN_REPORTABLE_ERROR_PROBABILITY));
    }

    static void posterior(final double logOdds, final double logPrior) {
        final String label = Double.toString(logOdds) + "," + Double.toString(logPrior);
        try {
            System.out.printf("posterior\t%s\t%s%n", label, Double.toString(
                    Mutect2FilteringEngine.posteriorProbabilityOfError(logOdds, logPrior)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void tumorLogOdds(final String label, final List<Double> values) {
        final VariantContextBuilder builder = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C));
        if (values != null) {
            builder.attribute("TLOD", values);
        }
        final VariantContext vc = builder.make();
        final double[] converted = Mutect2FilteringEngine.getTumorLogOdds(vc);
        System.out.printf("tumorlogodds\t%s\t%s%n", label,
                converted == null ? "null" : Arrays.toString(converted));
    }
}
