/*
 * AlleleSubsettingUtils.subsetAlleles, taken from the reference.
 *
 * What splitVariantContextToBiallelics does to a genotype that carries likelihoods, which is the
 * one thing the split port refuses for want of a measurement. The index machinery underneath is
 * already measured by the genotype-index suite; this is what the reference does with it.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE NEW PLs ARE THE OLD ONES PERMUTED AND THEN RESCALED, `scaleLogSpaceArrayForNumericalStability`
 *     subtracting the maximum, so the smallest phred of the subset is 0 whatever it was before;
 *   - AND A PL ARRAY OF THE WRONG LENGTH IS DROPPED ENTIRELY, `originalLikelihoods.length ==
 *     expectedNumLikelihoods ? ... : null`, so a genotype whose PLs do not match its ploidy and
 *     allele count loses them without a word;
 *   - THE GQ IS RECOMPUTED FROM THE NEW PLs, and only KEPT when the subset is the reference alone;
 *   - AND A GENOTYPE WITH A GQ AND NO PLs KEEPS ITS GQ, through a different branch;
 *   - A HOM-REF OR NO-CALL GENOTYPE WITH GQ 0 IS TURNED INTO A NO-CALL before the method is even
 *     consulted, AND IF ITS DP IS ALSO 0 EVERYTHING IS CLEARED: no PL, no DP, no AD, no GQ, no
 *     attributes;
 *   - BEST_MATCH_TO_ORIGINAL ALSO NO-CALLS ON GQ 0 with a leading PL of 0, which is a second and
 *     different test of the same idea;
 *   - USE_PLS_TO_ASSIGN CALLS THE MOST LIKELY GENOTYPE, and no-calls only when that is the
 *     reference AND the log10 GQ is above SUM_GL_THRESH_NOCALL, which is a GQ near zero: a
 *     CONFIDENT hom-ref is called hom-ref, an UNINFORMATIVE one is no-called. The threshold reads
 *     like a confidence test and is the opposite of one;
 *   - THE AD IS PERMUTED, NOT SUMMED, so an allele that is dropped takes its depth with it;
 *   - POSTERIORS AND PRIORS ARE REMOVED WHATEVER HAPPENS, being invalid for a new allele list,
 *     while any other attribute is kept;
 *   - AND SET_TO_NO_CALL EMPTIES THE CALL AND KEEPS EVERYTHING ELSE, where
 *     SET_TO_NO_CALL_NO_ANNOTATIONS keeps only the depth: two names one word apart, two different
 *     genotypes out.
 *
 * Output:
 *
 *     in\t<label>\t<the genotype>
 *     out\t<label>\t<method>\t<kept alleles>\t<the genotype>
 *     error\t<label>\t<exception class>:<message>
 *
 * A genotype is printed as `alleles|PL|GQ|AD|DP|attributes`, each empty when absent.
 *
 * Usage: SubsetAllelesDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;

import java.util.ArrayList;
import java.util.List;

public class SubsetAllelesDump {

    static final List<Allele> THREE = List.of(
            Allele.create("A", true), Allele.create("C", false), Allele.create("G", false));

    public static void main(final String[] args) {
        System.out.println("# SubsetAllelesDump: subsetAlleles, from the reference");

        // A het with likelihoods, subset each way.
        final Genotype het = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(1)))
                .PL(new int[] {50, 0, 60, 40, 30, 70}).GQ(50).DP(30).AD(new int[] {10, 12, 8}).make();
        run("het-keep-first", het, new int[] {0, 1}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        run("het-keep-second", het, new int[] {0, 2}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        run("het-keep-ref-only", het, new int[] {0}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        run("het-keep-all", het, new int[] {0, 1, 2}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        run("het-swapped", het, new int[] {0, 2, 1}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        // The same genotype through the other methods.
        run("het-use-pls", het, new int[] {0, 1}, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);
        run("het-no-call", het, new int[] {0, 1}, GenotypeAssignmentMethod.SET_TO_NO_CALL);
        run("het-no-annotations", het, new int[] {0, 1},
                GenotypeAssignmentMethod.SET_TO_NO_CALL_NO_ANNOTATIONS);

        // A CONFIDENT hom-ref, which USE_PLS_TO_ASSIGN calls hom-ref: the no-call is for the
        // uninformative one, the threshold being on a GQ near zero.
        final Genotype homRef = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(0)))
                .PL(new int[] {0, 60, 90, 60, 90, 90}).GQ(60).DP(30).make();
        run("hom-ref-use-pls", homRef, new int[] {0, 1}, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN);
        run("hom-ref-best-match", homRef, new int[] {0, 1},
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);

        // A hom-ref with GQ 0, which is no-called before the method is consulted, and the same with
        // DP 0, which clears everything.
        final Genotype noData = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(0)))
                .PL(new int[] {0, 0, 0, 0, 0, 0}).GQ(0).DP(3).AD(new int[] {1, 1, 1}).make();
        run("gq-zero", noData, new int[] {0, 1}, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
        final Genotype trulyNoData = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(0)))
                .PL(new int[] {0, 0, 0, 0, 0, 0}).GQ(0).DP(0).AD(new int[] {0, 0, 0}).make();
        run("gq-and-dp-zero", trulyNoData, new int[] {0, 1},
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);

        // A genotype whose PL array is the wrong length for its ploidy and alleles.
        final Genotype wrongLength = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(1)))
                .PL(new int[] {10, 0, 20}).GQ(10).make();
        run("wrong-length-pls", wrongLength, new int[] {0, 1},
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);

        // A genotype with a GQ and no likelihoods at all, which keeps its GQ.
        final Genotype noLikelihoods = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(1)))
                .GQ(40).DP(20).AD(new int[] {5, 6, 7}).make();
        run("no-likelihoods", noLikelihoods, new int[] {0, 1},
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);

        // Posteriors and priors, which are removed whatever happens.
        final Genotype withPosteriors = new GenotypeBuilder("s", List.of(THREE.get(0), THREE.get(1)))
                .PL(new int[] {50, 0, 60, 40, 30, 70}).GQ(50)
                .attribute("GP", new double[] {1.0, 2.0, 3.0, 4.0, 5.0, 6.0})
                .attribute("PG", new int[] {0, 1, 2, 3, 4, 5})
                .attribute("DPX", 7)
                .make();
        run("posteriors-removed", withPosteriors, new int[] {0, 1},
                GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL);
    }

    static void run(final String label, final Genotype genotype, final int[] keep,
                    final GenotypeAssignmentMethod method) {
        System.out.printf("in\t%s\t%s%n", label, render(genotype));
        final List<Allele> kept = new ArrayList<>();
        final List<String> keptNames = new ArrayList<>();
        for (final int index : keep) {
            kept.add(THREE.get(index));
            keptNames.add(String.valueOf(index));
        }
        final GenotypesContext result;
        try {
            result = AlleleSubsettingUtils.subsetAlleles(
                    GenotypesContext.create(new ArrayList<>(List.of(genotype))), 2, THREE, kept,
                    null, method);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        for (final Genotype out : result) {
            System.out.printf("out\t%s\t%s\t%s\t%s%n", label, method, String.join(",", keptNames),
                    render(out));
        }
    }

    /** `alleles|PL|GQ|AD|DP|attributes`, each empty when the genotype does not have it. */
    static String render(final Genotype genotype) {
        final List<String> alleles = new ArrayList<>();
        for (final Allele allele : genotype.getAlleles()) {
            alleles.add(allele.getDisplayString());
        }
        return String.join("/", alleles)
                + "|" + (genotype.hasPL() ? ints(genotype.getPL()) : "")
                + "|" + (genotype.hasGQ() ? String.valueOf(genotype.getGQ()) : "")
                + "|" + (genotype.hasAD() ? ints(genotype.getAD()) : "")
                + "|" + (genotype.hasDP() ? String.valueOf(genotype.getDP()) : "")
                + "|" + attributes(genotype);
    }

    static String attributes(final Genotype genotype) {
        final List<String> out = new ArrayList<>();
        genotype.getExtendedAttributes().keySet().stream().sorted().forEach(key ->
                out.add(key + "=" + render(genotype.getExtendedAttribute(key))));
        return String.join(";", out);
    }

    static String render(final Object value) {
        if (value instanceof int[]) {
            return ints((int[]) value);
        }
        if (value instanceof double[]) {
            final List<String> out = new ArrayList<>();
            for (final double element : (double[]) value) {
                out.add(String.valueOf(element));
            }
            return String.join(",", out);
        }
        if (value instanceof Object[]) {
            final List<String> out = new ArrayList<>();
            for (final Object element : (Object[]) value) {
                out.add(String.valueOf(element));
            }
            return String.join(",", out);
        }
        return String.valueOf(value);
    }

    static String ints(final int[] values) {
        final List<String> out = new ArrayList<>();
        for (final int value : values) {
            out.add(String.valueOf(value));
        }
        return String.join(",", out);
    }
}
