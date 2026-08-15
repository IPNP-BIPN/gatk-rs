/*
 * The EM iteration and the quantile initialisation, taken from the reference.
 *
 * What `SomaticClusteringModel.learnAndClearAccumulatedData` does with the data `record` collected,
 * and the two `learn` methods under it. Six behaviours this is built to catch.
 *
 *   - `BetaBinomialCluster.learn` IS GRADIENT ASCENT WITH FLOORS, ten epochs over the data at a rate
 *     of 0.01, each step `Math.max(alpha + rate * gradient * responsibility, 1.0)` and the same for
 *     beta at 0.5. THE FLOOR IS APPLIED EVERY STEP, not at the end, so a shape that would have gone
 *     below it and come back cannot;
 *   - THE UPDATE IS SEQUENTIAL WITHIN AN EPOCH: alpha and beta are read and written per datum, so
 *     the same data in a different order learn a different shape. The same data are therefore
 *     learned in two orders here;
 *   - `BinomialCluster.learn` ADDS 0.0001 TO BOTH SUMS before dividing, which is not a tie-breaker:
 *     it is what keeps a cluster with no responsibility at all from dividing zero by zero;
 *   - `initializeClusters` SPLITS PEAKS OFF THE BACKGROUND UNTIL THE BIC STOPS IMPROVING, at most
 *     five times, and REFUSES A PEAK BELOW THE 0.1 QUANTILE. What comes out is visible in
 *     `clusteringMetadata`, whose weights are formatted `%.4f` and whose cluster descriptions are
 *     `%.2f` and `%.3f`, so a cluster that moved in the last digits reads as unchanged;
 *   - THE PRIOR MAP IS REWRITTEN FROM THE CALLABLE-SITE COUNT, floored at 1.0e-8 for a SNV and
 *     1.0e-9 for an indel, so a model told how many sites were callable ends with different priors
 *     from one that was not -- and the lengths the getter inserted earlier are rewritten too;
 *   - AND THE QUANTILE RESPONSIBILITIES ARE BINOMIAL DENSITIES, `binomialProbability(n, k, f)` times
 *     `n + 1`, which is commons-math's `BinomialDistribution.probability` and its saddle point.
 *
 * Output:
 *
 *     learn\t<label>\t<alpha>,<beta>
 *     probe\t<label>\t<logLikelihood at the probe counts, at full precision>
 *     metadata\t<label>\t<key>=<value>
 *     binomprob\t<n>,<k>,<f>\t<value>
 *     prior\t<label>\t<log prior of a somatic variant>
 *     artifactprior\t<label>\t<log prior of variant versus artifact>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SomaticClusteringLearnDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.AlleleFractionCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.BetaBinomialCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.BinomialCluster;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.readorientation.BetaDistributionShape;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SomaticClusteringLearnDump {

    /** The counts a learned shape is probed at, which is what shows a move the toString rounds away. */
    static final int[][] PROBES = {{10, 5}, {100, 3}, {100, 50}};

    public static void main(final String[] args) {
        System.out.println("# SomaticClusteringLearnDump: what the model learns, and what it rounds away");

        // The binomial density the quantile responsibilities are built from.
        for (final int n : new int[] {10, 100}) {
            for (final int k : new int[] {0, 1, 5, 10}) {
                if (k <= n) {
                    for (final double f : new double[] {0.0, 0.01, 0.1, 0.5, 0.99, 1.0}) {
                        binomialProbability(n, k, f);
                    }
                }
            }
        }

        // `learn` on its own, over a data set with a clear allele fraction.
        final List<Datum> clonal = data(new int[][] {{100, 50}, {100, 48}, {100, 52}, {100, 51}});
        learned("betabinomial-flat-clonal", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                clonal, ones(clonal.size()));
        learned("betabinomial-highaf-clonal", new BetaBinomialCluster(new BetaDistributionShape(10, 1)),
                clonal, ones(clonal.size()));
        learned("binomial-clonal", new BinomialCluster(0.5), clonal, ones(clonal.size()));

        // The same data in the opposite order, which a sequential update does not answer the same.
        final List<Datum> reversed = new ArrayList<>(clonal);
        Collections.reverse(reversed);
        learned("betabinomial-flat-reversed", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                reversed, ones(reversed.size()));
        learned("binomial-reversed", new BinomialCluster(0.5), reversed, ones(reversed.size()));

        // Half the responsibility, which scales the step but not the floor.
        learned("betabinomial-flat-halved", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                clonal, filled(clonal.size(), 0.5));
        learned("binomial-halved", new BinomialCluster(0.5), clonal, filled(clonal.size(), 0.5));

        // No responsibility at all, where the beta-binomial cannot move and the binomial divides
        // 0.0001 by 0.0001.
        learned("betabinomial-flat-zero", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                clonal, filled(clonal.size(), 0.0));
        learned("binomial-zero", new BinomialCluster(0.5), clonal, filled(clonal.size(), 0.0));

        // No data at all.
        learned("betabinomial-flat-empty", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                List.of(), new double[0]);
        learned("binomial-empty", new BinomialCluster(0.5), List.of(), new double[0]);

        // A subclonal data set, low enough that the floors bite.
        final List<Datum> subclonal = data(new int[][] {{100, 5}, {100, 7}, {100, 4}, {100, 6}});
        learned("betabinomial-flat-subclonal", new BetaBinomialCluster(BetaDistributionShape.FLAT_BETA),
                subclonal, ones(subclonal.size()));
        learned("binomial-subclonal", new BinomialCluster(0.5), subclonal, ones(subclonal.size()));

        // The whole of `learnAndClearAccumulatedData`, over four data sets.
        model("clonal", clonal, 10000.0);
        model("subclonal", subclonal, 10000.0);
        model("bimodal", data(new int[][] {{100, 45}, {100, 48}, {100, 50}, {100, 14}, {100, 16},
                {100, 15}, {100, 47}, {100, 13}}), 10000.0);
        model("one-datum", data(new int[][] {{100, 50}}), 10000.0);
        model("no-data", List.of(), 10000.0);
        // The same clonal data with no callable-site count, which switches off the empirical priors.
        model("clonal-no-callable-sites", clonal, Double.NaN);
        // And with a count below one, which the constructor treats as none at all.
        model("clonal-zero-callable-sites", clonal, 0.0);

        // Learning twice: the second round starts from what the first left, and the clusters are not
        // initialised again.
        twice("twice", clonal, subclonal, 10000.0);
    }

    static void binomialProbability(final int n, final int k, final double f) {
        System.out.printf("binomprob\t%d,%d,%s\t%s%n", n, k, Double.toString(f),
                Double.toString(MathUtils.binomialProbability(n, k, f)));
    }

    /** `learn` on one cluster, printed as its own description and as three likelihoods. */
    static void learned(final String label, final AlleleFractionCluster cluster,
                        final List<Datum> data, final double[] responsibilities) {
        try {
            cluster.learn(data, responsibilities);
            System.out.printf("learn\t%s\t%s%n", label, cluster.toString());
            for (final int[] probe : PROBES) {
                System.out.printf("probe\t%s-%d,%d\t%s%n", label, probe[0], probe[1],
                        Double.toString(cluster.logLikelihood(probe[0], probe[1])));
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** One model, given data through `record` and then told to learn. */
    static void model(final String label, final List<Datum> data, final double callableSites) {
        try {
            final SomaticClusteringModel model = build(callableSites);
            record(model, data);
            model.learnAndClearAccumulatedData();
            metadata(label, model);
            priors(label, model);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** Two rounds of learning, with new data between them. */
    static void twice(final String label, final List<Datum> first, final List<Datum> second,
                      final double callableSites) {
        try {
            final SomaticClusteringModel model = build(callableSites);
            record(model, first);
            model.learnAndClearAccumulatedData();
            metadata(label + "-first", model);
            priors(label + "-first", model);
            record(model, second);
            model.learnAndClearAccumulatedData();
            metadata(label + "-second", model);
            priors(label + "-second", model);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static SomaticClusteringModel build(final double callableSites) {
        final List<MutectStats> stats = Double.isNaN(callableSites) ? List.of()
                : List.of(new MutectStats("callable", callableSites));
        return new SomaticClusteringModel(new M2FiltersArgumentCollection(), stats);
    }

    /** Feed the data in one datum at a time, which is what a record of one alternate does. */
    static void record(final SomaticClusteringModel model, final List<Datum> data) {
        for (final Datum datum : data) {
            final VariantContext vc = new VariantContextBuilder("dump", "chr1", 100, 100,
                    List.of(Allele.REF_A, Allele.ALT_C)).make();
            model.record(new int[] {datum.getTotalCount() - datum.getAltCount(), datum.getAltCount()},
                    new double[] {datum.getTumorLogOdds()},
                    List.of(datum.getArtifactProb()), List.of(0.0), vc);
        }
    }

    static void metadata(final String label, final SomaticClusteringModel model) {
        for (final Pair<String, String> entry : model.clusteringMetadata()) {
            System.out.printf("metadata\t%s\t%s=%s%n", label,
                    ReferenceQueryDump.escape(entry.getLeft()),
                    ReferenceQueryDump.escape(entry.getRight()));
        }
    }

    /** The priors after learning, at three lengths inside the window and one outside it. */
    static void priors(final String label, final SomaticClusteringModel model) {
        for (final int length : new int[] {0, 1, -1, 40}) {
            final VariantContext vc = length == 0
                    ? snp()
                    : length > 0
                            ? new VariantContextBuilder("dump", "chr1", 100, 100,
                                    List.of(Allele.create("A", true),
                                            Allele.create("A" + "C".repeat(length), false))).make()
                            : new VariantContextBuilder("dump", "chr1", 100, 100 - length,
                                    List.of(Allele.create("A" + "C".repeat(-length), true),
                                            Allele.create("A", false))).make();
            System.out.printf("prior\t%s-%d\t%s%n", label, length,
                    Double.toString(model.getLogPriorOfSomaticVariant(vc, 0)));
        }
        System.out.printf("artifactprior\t%s\t%s%n", label,
                Double.toString(model.getLogPriorOfVariantVersusArtifact()));
    }

    static VariantContext snp() {
        return new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C)).make();
    }

    /** Data at a TLOD high enough to be somatic, one datum per total and alternate count. */
    static List<Datum> data(final int[][] counts) {
        final List<Datum> data = new ArrayList<>();
        for (final int[] pair : counts) {
            data.add(new Datum(20.0, 0.0, 0.0, pair[1], pair[0], 0));
        }
        return data;
    }

    static double[] ones(final int size) {
        return filled(size, 1.0);
    }

    static double[] filled(final int size, final double value) {
        final double[] array = new double[size];
        java.util.Arrays.fill(array, value);
        return array;
    }
}
