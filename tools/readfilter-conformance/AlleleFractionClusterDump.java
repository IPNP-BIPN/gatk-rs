/*
 * The two allele-fraction clusters, taken from the reference.
 *
 * `SomaticClusteringModel` carries a `BinomialCluster` per allele fraction it tracks and a
 * `BetaBinomialCluster` beside them, and asks each for two numbers: a plain beta-binomial
 * `logLikelihood(total, alt)` and a `correctedLogLikelihood(datum)`, which is not a likelihood but a
 * TLOD corrected for the cluster's shape. Five behaviours this is built to catch.
 *
 *   - THE CORRECTION IS FOUR DIRICHLET NORMALISATIONS, `g(a, b) - g(a + alt, b + ref) - g(1, 1) +
 *     g(1 + alt, 1 + ref)`, where `g` is `logGamma(sum) - sum(logGamma)`. The original shape is
 *     always FLAT_BETA, so two of the four terms are constants -- and they do not cancel in doubles;
 *   - `BinomialCluster` HAS NO BINOMIAL IN IT: its constructor turns a mean into a beta shape with a
 *     fixed standard-deviation-over-mean of 0.01, CLAMPING THE MEAN AT `1 - 0.01` first, so a
 *     cluster asked for a mean of one is not given one and its two shapes are not symmetric;
 *   - `alphaPlusBeta` IS `((1 - mean) / (mean * 0.0001)) - 1`, which is enormous at a small mean:
 *     these are the shapes that reach logBeta's `a >= 10` branch;
 *   - `Datum`'s NON-SEQUENCING-ERROR PROBABILITY IS A COMBINATION computed in the constructor,
 *     `1 - (1 - artifactProb) * (1 - nonSomaticProb)`, and the datum keeps no field for the second;
 *   - AND `BetaDistributionShape` REFUSES IN TWO DIFFERENT CLASSES, `ParamUtils` on alpha and
 *     `Utils.validateArg` on beta, with two differently-worded messages.
 *
 * `learn` is out of scope: ten epochs of gradient ascent over `Gamma.digamma` belong with the
 * model's fitting.
 *
 * Output:
 *
 *     dirichlet\t<parameters>\t<value>
 *     fuzzy\t<mean>\t<alpha>,<beta>
 *     loglike\t<label>\t<value>
 *     corrected\t<label>\t<value>
 *     datum\t<label>\t<non-sequencing-error probability>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AlleleFractionClusterDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.SomaticLikelihoodsEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.AlleleFractionCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.BetaBinomialCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.BinomialCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;

public class AlleleFractionClusterDump {

    /** The means a model tracks, from a subclone to a clonal variant, and the two ends. */
    static final double[] MEANS = {0.001, 0.01, 0.1, 0.25, 0.5, 0.9, 0.99, 1.0, 2.0};

    /** Counts spanning a shallow site, a deep one, and a site with no alternate reads at all. */
    static final int[][] COUNTS = {{10, 0}, {10, 1}, {10, 5}, {10, 10}, {100, 3}, {100, 50}, {1000, 7}};

    public static void main(final String[] args) {
        System.out.println("# AlleleFractionClusterDump: two clusters, and the correction under both");

        // The normalisation the correction is built from, on its own.
        dirichlet(1.0, 1.0);
        dirichlet(1.0);
        dirichlet(0.5);
        dirichlet(10.0, 1.0);
        dirichlet(1.0, 10.0);
        dirichlet(1.0e-6, 1.0);
        dirichlet(1.0, 2.0, 3.0);
        dirichlet(9999.0, 99.0);
        dirichlet(0.0, 1.0);

        // The mean-to-shape map, including the clamp at the top end.
        for (final double mean : MEANS) {
            fuzzy(mean);
        }

        for (final double mean : MEANS) {
            final AlleleFractionCluster cluster = binomialCluster(mean);
            if (cluster == null) {
                continue;
            }
            for (final int[] counts : COUNTS) {
                logLikelihood("binomial-" + Double.toString(mean), cluster, counts[0], counts[1]);
                corrected("binomial-" + Double.toString(mean), cluster, datum(counts[0], counts[1]));
            }
        }

        // The beta-binomial cluster, at the flat shape and at two the model can reach.
        final double[][] shapes = {{1.0, 1.0}, {10.0, 1.0}, {1.0, 10.0}, {0.5, 0.5}};
        for (final double[] shape : shapes) {
            final String label = "betabinomial-" + Double.toString(shape[0]) + "," + Double.toString(shape[1]);
            final AlleleFractionCluster cluster =
                    new BetaBinomialCluster(new BetaDistributionShape(shape[0], shape[1]));
            for (final int[] counts : COUNTS) {
                logLikelihood(label, cluster, counts[0], counts[1]);
                corrected(label, cluster, datum(counts[0], counts[1]));
            }
        }

        // The TLOD is carried through the correction unchanged, so a datum's odds move the answer.
        final AlleleFractionCluster flat = new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA);
        for (final double odds : new double[] {0.0, 5.0, -5.0, Double.NEGATIVE_INFINITY, Double.NaN}) {
            corrected("tlod-" + Double.toString(odds), flat,
                    new Datum(odds, 0.0, 0.0, 5, 10, 0));
        }

        // The combination in `Datum`'s constructor, which is not either input.
        datumProbability("both-zero", 0.0, 0.0);
        datumProbability("artifact-only", 0.3, 0.0);
        datumProbability("non-somatic-only", 0.0, 0.3);
        datumProbability("both", 0.3, 0.3);
        datumProbability("artifact-certain", 1.0, 0.5);
        datumProbability("tiny", 1.0e-10, 1.0e-10);

        // The two refusals, which are two different exception classes.
        shape("alpha-zero", 0.0, 1.0);
        shape("beta-zero", 1.0, 0.0);
        shape("alpha-negative", -1.0, 1.0);
        shape("beta-nan", 1.0, Double.NaN);
    }

    static void dirichlet(final double... parameters) {
        final StringBuilder label = new StringBuilder();
        for (final double parameter : parameters) {
            if (label.length() > 0) {
                label.append(',');
            }
            label.append(Double.toString(parameter));
        }
        System.out.printf("dirichlet\t%s\t%s%n", label,
                Double.toString(SomaticLikelihoodsEngine.logDirichletNormalization(parameters)));
    }

    /** `BinomialCluster`'s mean-to-shape map, read back through the cluster's own toString. */
    static void fuzzy(final double mean) {
        final BetaDistributionShape shape = fuzzyBinomial(mean);
        if (shape == null) {
            return;
        }
        System.out.printf("fuzzy\t%s\t%s,%s%n", Double.toString(mean),
                Double.toString(shape.getAlpha()), Double.toString(shape.getBeta()));
    }

    /**
     * `BinomialCluster.getFuzzyBinomial`, which is private, reproduced here so its two shapes can be
     * printed. The cluster built beside it is what every likelihood below actually goes through, so a
     * disagreement between the two would show up as a wrong likelihood rather than being hidden.
     */
    static BetaDistributionShape fuzzyBinomial(final double unboundedMean) {
        final double mean = Math.min(unboundedMean, 1 - 0.01);
        final double alphaPlusBeta = ((1 - mean) / (mean * 0.01 * 0.01)) - 1;
        try {
            return new BetaDistributionShape(mean * alphaPlusBeta, alphaPlusBeta - mean * alphaPlusBeta);
        } catch (final Exception e) {
            System.out.printf("error\tfuzzy-%s\t%s:%s%n", Double.toString(unboundedMean),
                    e.getClass().getName(), ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return null;
        }
    }

    static AlleleFractionCluster binomialCluster(final double mean) {
        try {
            return new BinomialCluster(mean);
        } catch (final Exception e) {
            System.out.printf("error\tbinomial-%s\t%s:%s%n", Double.toString(mean),
                    e.getClass().getName(), ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return null;
        }
    }

    static void logLikelihood(final String prefix, final AlleleFractionCluster cluster,
                              final int total, final int alt) {
        final String label = prefix + "-" + total + "," + alt;
        try {
            System.out.printf("loglike\t%s\t%s%n", label,
                    Double.toString(cluster.logLikelihood(total, alt)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void corrected(final String prefix, final AlleleFractionCluster cluster, final Datum datum) {
        final String label = prefix + "-" + datum.getTotalCount() + "," + datum.getAltCount();
        try {
            System.out.printf("corrected\t%s\t%s%n", label,
                    Double.toString(cluster.correctedLogLikelihood(datum)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void datumProbability(final String label, final double artifact, final double nonSomatic) {
        System.out.printf("datum\t%s\t%s%n", label,
                Double.toString(new Datum(0.0, artifact, nonSomatic, 5, 10, 0)
                        .getNonSequencingErrorProb()));
    }

    static void shape(final String label, final double alpha, final double beta) {
        try {
            final BetaDistributionShape shape = new BetaDistributionShape(alpha, beta);
            System.out.printf("fuzzy\t%s\t%s,%s%n", label, Double.toString(shape.getAlpha()),
                    Double.toString(shape.getBeta()));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static Datum datum(final int total, final int alt) {
        return new Datum(0.0, 0.0, 0.0, alt, total, 0);
    }
}
