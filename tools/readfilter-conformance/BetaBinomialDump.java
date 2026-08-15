/*
 * The beta-binomial likelihood Mutect's clusters are built on, taken from the reference.
 *
 * `BetaBinomialDistribution.logProbability` is three commons-math calls:
 *
 *     binomialCoefficientLog(n, k) + logBeta(k + alpha, n - k + beta) - logBeta(alpha, beta)
 *
 * and it is what `SomaticClusteringModel.logLikelihoodGivenSomatic` sums over its clusters. Five
 * behaviours this is built to catch.
 *
 *   - THE THREE PIECES ARE MEASURED SEPARATELY as well as together, because a port has a choice at
 *     `logBeta`: commons-math computes it by its own expansion, and the obvious identity
 *     `logGamma(p) + logGamma(q) - logGamma(p + q)` is a different sequence of roundings. The
 *     golden is what says whether the two agree on the values the clusters actually use;
 *   - THE FLAT BETA IS THE BACKGROUND CLUSTER, `alpha = beta = 1`, where `logBeta(1, 1)` is
 *     NEGATIVE zero and the beta-binomial is uniform over the counts — though not bit-identically
 *     so: at `n = 10` the answer is `...707` at `k = 0`, `...716` at `k = 1` and `...700` at
 *     `k = 5`, the cancellation between the coefficient and the two beta terms being approximate;
 *   - THE HIGH-AF CLUSTER IS `alpha = 10, beta = 1`, so its likelihood rises steeply with the
 *     alternate count and the two clusters disagree most at the ends;
 *   - k GREATER THAN n IS NEGATIVE INFINITY by an explicit branch rather than by arithmetic;
 *   - AND A COUNT OF ZERO IS NOT A SPECIAL CASE: `binomialCoefficientLog(n, 0)` is zero and the
 *     rest carries the answer.
 *
 * Output:
 *
 *     logbeta\t<p>,<q>\t<value>
 *     binomlog\t<n>,<k>\t<value>
 *     betabinom\t<alpha>,<beta>,<n>,<k>\t<value>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: BetaBinomialDump
 */

import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;

public class BetaBinomialDump {

    /** The two shapes `SomaticClusteringModel` starts with: the flat background and the high-AF cluster. */
    static final double[][] CLUSTERS = {{1.0, 1.0}, {10.0, 1.0}};

    public static void main(final String[] args) {
        System.out.println("# BetaBinomialDump: three commons-math calls, and the likelihood they make");

        final double[] shapes = {0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0};
        for (final double p : shapes) {
            for (final double q : shapes) {
                logBeta(p, q);
            }
        }
        // The shapes a fuzzy binomial cluster is built from: a mean with a small standard deviation.
        logBeta(0.5, 0.5);
        logBeta(4999.5, 4999.5);
        logBeta(1.0e-6, 1.0);

        for (final int n : new int[] {0, 1, 10, 100, 1000}) {
            for (final int k : counts(n)) {
                binomialCoefficientLog(n, k);
            }
        }

        for (final double[] cluster : CLUSTERS) {
            for (final int n : new int[] {1, 10, 100}) {
                for (final int k : counts(n)) {
                    betaBinomial(cluster[0], cluster[1], n, k);
                }
            }
        }
        // A count past the total, which is an explicit branch rather than arithmetic.
        betaBinomial(1.0, 1.0, 10, 11);
        // And a negative count, which the argument check refuses.
        betaBinomial(1.0, 1.0, 10, -1);
    }

    /** The counts of interest for a total, without the repeats a small total would produce. */
    static java.util.List<Integer> counts(final int n) {
        final java.util.Set<Integer> distinct = new java.util.LinkedHashSet<>();
        for (final int k : new int[] {0, 1, n / 2, n}) {
            if (k >= 0 && k <= n) {
                distinct.add(k);
            }
        }
        return new java.util.ArrayList<>(distinct);
    }

    static void logBeta(final double p, final double q) {
        System.out.printf("logbeta\t%s,%s\t%s%n", Double.toString(p), Double.toString(q),
                Double.toString(Beta.logBeta(p, q)));
    }

    static void binomialCoefficientLog(final int n, final int k) {
        System.out.printf("binomlog\t%d,%d\t%s%n", n, k,
                Double.toString(CombinatoricsUtils.binomialCoefficientLog(n, k)));
    }

    static void betaBinomial(final double alpha, final double beta, final int n, final int k) {
        final String label = Double.toString(alpha) + "," + Double.toString(beta) + "," + n + "," + k;
        try {
            System.out.printf("betabinom\t%s\t%s%n", label, Double.toString(
                    new BetaBinomialDistribution(null, alpha, beta, n).logProbability(k)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
