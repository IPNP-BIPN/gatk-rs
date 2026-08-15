/*
 * SomaticClusteringModel as constructed, taken from the reference.
 *
 * The state the constructor puts the model in, and the four read-only answers it gives from that
 * state, before any learning has happened. Five behaviours this is built to catch.
 *
 *   - `getLogPriorOfSomaticVariant` MUTATES THE MAP IT READS. An indel length outside the +/-10
 *     window is not in `logVariantPriors`, so the getter INSERTS THE CURRENT MINIMUM for it and
 *     returns that. Asking about a length once changes what the map holds, so the same question
 *     asked in a different order can be answered differently -- which is why the same lengths are
 *     asked twice here, in two models built the same way;
 *   - A SNV'S PRIOR HAS `LOG_ONE_THIRD` ADDED TO IT AND AN INDEL'S DOES NOT, so the two are not one
 *     prior with a different number;
 *   - `record` MUTATES THE CALLER'S ARRAY: it zeroes `tumorADs` for symbolic alleles in place before
 *     summing, so the array the engine passed in comes back changed. The array is printed before and
 *     after;
 *   - TWO THRESHOLDS THAT LOOK ALIKE DO DIFFERENT THINGS: an artifact probability above 0.9
 *     increments the obvious-artifact count and drops the datum, while a non-somatic probability
 *     above 0.9 drops it silently. Neither is visible directly, so what is measured is the count of
 *     data the model accumulated, read back through `learnAndClearAccumulatedData` being able to run;
 *   - AND THE INITIAL WEIGHTS ARE `log1p(0.01)` AND `log(0.01)`, which are not a normalised pair:
 *     the background weight is log(1.01) and not log(0.99). Their effect is visible in
 *     `logLikelihoodGivenSomatic`, a logSumExp over the two initial clusters and the first thing
 *     here that carries `exp`.
 *
 * The EM iteration, the quantile initialisation and everything that learns are out of scope.
 *
 * Output:
 *
 *     prior\t<label>\t<log prior of a somatic variant>
 *     artifactprior\t<label>\t<log prior of variant versus artifact>
 *     loglike\t<total>,<alt>\t<logLikelihoodGivenSomatic>
 *     seqerror\t<label>\t<probabilityOfSequencingError>
 *     ads\t<label>\t<the caller's array, after record returned>
 *     indellength\t<label>\t<indelLength(vc, altIndex)>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SomaticClusteringModelDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.MutectStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.Datum;
import org.broadinstitute.hellbender.tools.walkers.mutect.clustering.SomaticClusteringModel;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;

import java.util.Arrays;
import java.util.List;

public class SomaticClusteringModelDump {

    public static void main(final String[] args) {
        System.out.println("# SomaticClusteringModelDump: the model as constructed, before it learns");

        // The prior for every length in the window, and the two just outside it.
        final SomaticClusteringModel inWindow = model();
        for (final int length : new int[] {0, 1, -1, 2, -2, 10, -10}) {
            prior("in-window-" + length, inWindow, length);
        }

        // Outside the window, where the getter inserts. Two models, asked the same lengths in
        // opposite orders, which is what shows the insertion changing the later answer.
        final SomaticClusteringModel ascending = model();
        for (final int length : new int[] {11, -11, 50, -50}) {
            prior("ascending-" + length, ascending, length);
        }
        final SomaticClusteringModel descending = model();
        for (final int length : new int[] {-50, 50, -11, 11}) {
            prior("descending-" + length, descending, length);
        }
        // And the same length twice, which the second time is a plain read.
        final SomaticClusteringModel repeated = model();
        prior("repeated-first", repeated, 11);
        prior("repeated-second", repeated, 11);
        prior("repeated-in-window", repeated, 0);

        // The prior a non-default argument collection produces.
        final M2FiltersArgumentCollection loud = new M2FiltersArgumentCollection();
        loud.logSNVPrior = -5.0;
        loud.logIndelPrior = -4.0;
        loud.initialLogPriorOfVariantVersusArtifact = -1.0;
        final SomaticClusteringModel tuned = new SomaticClusteringModel(loud, List.of());
        prior("tuned-snv", tuned, 0);
        prior("tuned-indel", tuned, 3);
        System.out.printf("artifactprior\ttuned\t%s%n",
                Double.toString(tuned.getLogPriorOfVariantVersusArtifact()));
        System.out.printf("artifactprior\tdefault\t%s%n",
                Double.toString(model().getLogPriorOfVariantVersusArtifact()));

        // The mitochondrial defaults, which only apply while the field still holds the default: the
        // getter compares the field with the default by ==, so a mitochondrial run that was told the
        // default value explicitly gets the mitochondrial number and one told anything else does not.
        final M2FiltersArgumentCollection mito = new M2FiltersArgumentCollection();
        mito.mitochondria = true;
        final SomaticClusteringModel mitoModel = new SomaticClusteringModel(mito, List.of());
        prior("mito-snv", mitoModel, 0);
        prior("mito-indel", mitoModel, 3);
        final M2FiltersArgumentCollection mitoSet = new M2FiltersArgumentCollection();
        mitoSet.mitochondria = true;
        mitoSet.logSNVPrior = -5.0;
        final SomaticClusteringModel mitoSetModel = new SomaticClusteringModel(mitoSet, List.of());
        prior("mito-snv-overridden", mitoSetModel, 0);
        prior("mito-indel-untouched", mitoSetModel, 3);

        // The emission likelihood over the two initial clusters.
        final SomaticClusteringModel fresh = model();
        for (final int[] counts : new int[][] {{10, 0}, {10, 1}, {10, 5}, {10, 10}, {100, 3},
                {100, 50}, {1000, 7}, {0, 0}}) {
            logLikelihood(fresh, counts[0], counts[1]);
        }

        // The sequencing-error posterior, which is the same weighted sum through the corrected
        // likelihood and then a prior.
        for (final double odds : new double[] {0.0, 5.0, 20.0, -5.0}) {
            for (final int[] counts : new int[][] {{10, 5}, {100, 3}}) {
                sequencingError(fresh, odds, counts[0], counts[1], 0);
            }
        }
        // The same counts at an indel length, which reaches a different prior.
        sequencingError(fresh, 5.0, 10, 5, 3);
        sequencingError(fresh, 5.0, 10, 5, -3);
        // And at a length outside the window, which inserts on the way through.
        sequencingError(fresh, 5.0, 10, 5, 40);

        // `record`, which changes the array it was given.
        recorded("one-alt", new int[] {80, 20},
                new double[] {6.0}, List.of(0.0), List.of(0.0),
                vc(Allele.REF_A, Allele.ALT_C));
        recorded("symbolic-alt", new int[] {80, 20, 5},
                new double[] {6.0, 6.0}, List.of(0.0, 0.0), List.of(0.0, 0.0),
                vc(Allele.REF_A, Allele.ALT_C, Allele.create("<NON_REF>", false)));
        recorded("obvious-artifact", new int[] {80, 20},
                new double[] {6.0}, List.of(0.95), List.of(0.0),
                vc(Allele.REF_A, Allele.ALT_C));
        recorded("obvious-non-somatic", new int[] {80, 20},
                new double[] {6.0}, List.of(0.0), List.of(0.95),
                vc(Allele.REF_A, Allele.ALT_C));
        // Exactly at the threshold, which is `>` and therefore keeps the datum.
        recorded("at-threshold", new int[] {80, 20},
                new double[] {6.0}, List.of(0.9), List.of(0.9),
                vc(Allele.REF_A, Allele.ALT_C));
        // An array whose length does not match the record's allele count.
        recorded("short-array", new int[] {80},
                new double[] {6.0}, List.of(0.0), List.of(0.0),
                vc(Allele.REF_A, Allele.ALT_C));

        // The indel length a record reports per alternate allele.
        indelLength("snp", vc(Allele.REF_A, Allele.ALT_C), 0);
        indelLength("insertion", vc(Allele.create("A", true), Allele.create("ACCC", false)), 0);
        indelLength("deletion", vc(Allele.create("ACCC", true), Allele.create("A", false)), 0);
        indelLength("symbolic", vc(Allele.REF_A, Allele.create("<DEL>", false)), 0);
    }

    static SomaticClusteringModel model() {
        return new SomaticClusteringModel(new M2FiltersArgumentCollection(), List.of());
    }

    /** A model told how many callable sites there were, which is what switches on the empirical priors. */
    static SomaticClusteringModel modelWithStats(final double callableSites) {
        return new SomaticClusteringModel(new M2FiltersArgumentCollection(),
                List.of(new MutectStats("callable", callableSites)));
    }

    static void prior(final String label, final SomaticClusteringModel model, final int indelLength) {
        final VariantContext vc = indelLength == 0 ? vc(Allele.REF_A, Allele.ALT_C)
                : indelLength > 0 ? vc(Allele.create("A", true), Allele.create("A" + "C".repeat(indelLength), false))
                : vc(Allele.create("A" + "C".repeat(-indelLength), true), Allele.create("A", false));
        try {
            System.out.printf("prior\t%s\t%s%n", label,
                    Double.toString(model.getLogPriorOfSomaticVariant(vc, 0)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void logLikelihood(final SomaticClusteringModel model, final int total, final int alt) {
        final String label = total + "," + alt;
        try {
            System.out.printf("loglike\t%s\t%s%n", label,
                    Double.toString(model.logLikelihoodGivenSomatic(total, alt)));
        } catch (final Exception e) {
            System.out.printf("error\tloglike-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void sequencingError(final SomaticClusteringModel model, final double odds,
                                final int total, final int alt, final int indelLength) {
        final String label = Double.toString(odds) + "-" + total + "," + alt + "-" + indelLength;
        try {
            System.out.printf("seqerror\t%s\t%s%n", label, Double.toString(
                    model.probabilityOfSequencingError(new Datum(odds, 0.0, 0.0, alt, total, indelLength))));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void recorded(final String label, final int[] tumorADs, final double[] tumorLogOdds,
                         final List<Double> artifactProbabilities,
                         final List<Double> nonSomaticProbabilities, final VariantContext vc) {
        final SomaticClusteringModel model = modelWithStats(1000.0);
        System.out.printf("ads\t%s-before\t%s%n", label, Arrays.toString(tumorADs));
        try {
            model.record(tumorADs, tumorLogOdds, artifactProbabilities, nonSomaticProbabilities, vc);
            System.out.printf("ads\t%s-after\t%s%n", label, Arrays.toString(tumorADs));
            // What the model accumulated, read back through the prior the EM iteration rewrites: a
            // datum that was dropped cannot move it.
            model.learnAndClearAccumulatedData();
            System.out.printf("artifactprior\t%s-learned\t%s%n", label,
                    Double.toString(model.getLogPriorOfVariantVersusArtifact()));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void indelLength(final String label, final VariantContext vc, final int altIndex) {
        try {
            System.out.printf("indellength\t%s\t%d%n", label,
                    SomaticClusteringModel.indelLength(vc, altIndex));
        } catch (final Exception e) {
            System.out.printf("error\tindellength-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static VariantContext vc(final Allele... alleles) {
        final Allele reference = alleles[0];
        return new VariantContextBuilder("dump", "chr1", 100,
                100 + reference.length() - 1, List.of(alleles)).make();
    }
}
