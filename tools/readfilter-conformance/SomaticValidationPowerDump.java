/*
 * The power a validation pileup has, taken from the reference.
 *
 * `ValidateBasicSomaticShortMutations` asks: given what the discovery pileup saw, how likely is a
 * validation pileup of this depth to show the variant at all? The answer is one minus a
 * beta-binomial cumulative probability, and the count that counts as "showing it" is a binomial
 * quantile floored at two.
 *
 * The beta-binomial's `logProbability` is already pinned by the `beta-binomial` suite and ported.
 * What is measured here is the layer above it, which that suite does not reach.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE CUMULATIVE PROBABILITY IS A PLAIN LOOP over `probability(i)`, summed with `+=` and NOT
 *     compensated, so it is not the `DoubleStream.sum` every other accumulation in this port has
 *     had to be;
 *   - EVERY TERM GOES THROUGH exp(log), because `probability(k)` is `Math.exp(logProbability(k))`.
 *     Under decision 0014 that is the platform's `exp`, so these rows carry the same one-ulp
 *     exposure the clustering model's do;
 *   - `k > n` SHORT-CIRCUITS to negative infinity before any of the three log terms is evaluated,
 *     and the cumulative probability therefore stops growing past the trial count;
 *   - A NEGATIVE k IS REFUSED by both `logProbability` and `cumulativeProbability` rather than
 *     answering zero, which is what makes `minCount - 1` a precondition in the power calculation;
 *   - THE MINIMUM COUNT IS A BINOMIAL QUANTILE FLOORED AT TWO. `Math.max(inverseCumulative(0.99),
 *     2)`, so a shallow pileup and a clean one both answer two;
 *   - AND THE POWER IS `1 - cumulativeProbability(minCount - 1)`, computed on a distribution whose
 *     shapes are `discoveryAlt + 1` and `discoveryTotal - discoveryAlt + 1` -- never zero, which is
 *     what keeps the constructor's own refusals out of reach.
 *
 * Output:
 *
 *     cumulative\t<alpha bits>,<beta bits>,<n>,<k>=<bits>
 *     moments\t<alpha bits>,<beta bits>,<n>=<mean bits>,<variance bits>
 *     power\t<validationTotal>,<discoveryAlt>,<discoveryTotal>,<minCount>=<bits>
 *     mincount\t<validationTotal>,<maxSignalRatio bits>=<count>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SomaticValidationPowerDump
 */

import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.PowerCalculationUtils;

public class SomaticValidationPowerDump {

    /** Shapes a real discovery pileup produces, which are `altCount + 1` and `refCount + 1`. */
    static final double[][] SHAPES = {
            {1.0, 1.0}, {2.0, 9.0}, {6.0, 6.0}, {11.0, 1.0}, {0.5, 0.5}, {31.0, 71.0},
    };

    /** Validation depths. */
    static final int[] TRIALS = {0, 1, 5, 20};

    public static void main(final String[] args) {
        System.out.println("# SomaticValidationPowerDump: the power a validation pileup has");

        for (final double[] shape : SHAPES) {
            for (final int n : TRIALS) {
                final BetaBinomialDistribution distribution =
                        new BetaBinomialDistribution(null, shape[0], shape[1], n);
                System.out.printf("moments\t%016x,%016x,%d=%016x,%016x%n",
                        Double.doubleToRawLongBits(shape[0]), Double.doubleToRawLongBits(shape[1]), n,
                        Double.doubleToRawLongBits(distribution.getNumericalMean()),
                        Double.doubleToRawLongBits(distribution.getNumericalVariance()));
                // One past the trial count, where the short circuit shows.
                for (int k = 0; k <= n + 1; k++) {
                    System.out.printf("cumulative\t%016x,%016x,%d,%d=%016x%n",
                            Double.doubleToRawLongBits(shape[0]),
                            Double.doubleToRawLongBits(shape[1]), n, k,
                            Double.doubleToRawLongBits(distribution.cumulativeProbability(k)));
                }
            }
        }

        for (final int validationTotal : new int[] {0, 1, 10, 50}) {
            for (final int[] discovery : new int[][] {{1, 10}, {5, 10}, {10, 10}, {30, 100}}) {
                for (final int minCount : new int[] {2, 3, 5}) {
                    System.out.printf("power\t%d,%d,%d,%d=%016x%n",
                            validationTotal, discovery[0], discovery[1], minCount,
                            Double.doubleToRawLongBits(PowerCalculationUtils.calculatePower(
                                    validationTotal, discovery[0], discovery[1], minCount)));
                }
            }
        }

        for (final int validationTotal : new int[] {0, 1, 10, 50, 317}) {
            for (final double ratio : new double[] {0.0, 0.01, 0.05, 0.2, 1.0}) {
                System.out.printf("mincount\t%d,%016x=%d%n", validationTotal,
                        Double.doubleToRawLongBits(ratio),
                        PowerCalculationUtils.calculateMinCountForSignal(validationTotal, ratio));
            }
        }

        error("negative-k-cumulative",
                () -> new BetaBinomialDistribution(null, 1.0, 1.0, 5).cumulativeProbability(-1));
        error("ratio-above-one", () -> (double) PowerCalculationUtils.calculateMinCountForSignal(10, 1.5));
        error("negative-total", () -> (double) PowerCalculationUtils.calculateMinCountForSignal(-1, 0.1));
    }

    interface Body {
        double run();
    }

    static void error(final String label, final Body body) {
        try {
            System.out.printf("unexpected\t%s\t%016x%n", label,
                    Double.doubleToRawLongBits(body.run()));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }
}
