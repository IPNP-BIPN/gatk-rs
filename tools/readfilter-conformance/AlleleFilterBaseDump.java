/*
 * The allele filter's shared machinery and the tumour-evidence filter, taken from the reference.
 *
 * `Mutect2AlleleFilter` and `Mutect2Filter` are what every per-allele filter stands on, and
 * `TumorEvidenceFilter` is the thinnest thing built on them. Six behaviours this is built to catch.
 *
 *   - `combineDataByAllele` ZIPS TWO ITERATORS AND STOPS AT THE SHORTER ONE. A genotype whose
 *     per-allele list is shorter than the record's allele count contributes to the first alleles and
 *     NOTHING to the rest, silently; a longer one has its tail dropped. No exception, no padding;
 *   - THE MAP IS KEYED BY Allele AND MERGED WITH `(a, b) -> a`, so a record carrying the same allele
 *     twice collapses to one key and the second one's data lands in the first one's list;
 *   - A MISSING REQUIRED ANNOTATION IS AN EMPTY LIST, NOT A ZERO, and an empty list is what
 *     ErrorProbabilities drops entirely rather than counting;
 *   - `weightedMedianPosteriorProbability` SORTS THE CALLER'S LIST IN PLACE and returns the first
 *     posterior at which twice the cumulative alt depth reaches the total. An empty list is 0, and
 *     the comparison is `>=`, so an even split returns the lower half's last element;
 *   - `sumADsOverSamples` INDEXES EVERY GENOTYPE'S AD BY THE RECORD'S ALLELE COUNT, so a genotype
 *     with a shorter AD array is an ArrayIndexOutOfBoundsException rather than a skip;
 *   - AND `TumorEvidenceFilter` IS `probabilityOfSequencingError` PER ALLELE over a Datum built with
 *     both probabilities zero, so it reads the clustering model without contributing to it.
 *
 * Output:
 *
 *     byallele\t<label>\t<allele>=<the list gathered for it>
 *     altbyallele\t<label>\t<allele>=<the list gathered for it>
 *     median\t<label>\t<the weighted median posterior>
 *     sorted\t<label>\t<the caller's list after the call, which the call reordered>
 *     ads\t<label>\t<sumADsOverSamples, as an array>
 *     evidence\t<label>\t<one error probability per alternate allele>
 *     name\t<label>\t<filterName>,<errorType>,<phredScaledPosteriorAnnotationName>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AlleleFilterBaseDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2AlleleFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.TumorEvidenceFilter;

import java.io.File;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class AlleleFilterBaseDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# AlleleFilterBaseDump: the base class every per-allele filter stands on");

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "T2", "N1"));
        final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(
                new M2FiltersArgumentCollection(), header, new File("no-such-stats-file.tsv"));

        // Two tumour samples and a normal, each carrying one value per allele.
        final VariantContext triallelic = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G))
                .attribute("TLOD", List.of(20.0, 6.0))
                .genotypes(List.of(
                        genotype("T1", new int[] {80, 20, 5}, List.of(1, 2, 3)),
                        genotype("T2", new int[] {70, 30, 10}, List.of(4, 5, 6)),
                        genotype("N1", new int[] {99, 1, 0}, List.of(7, 8, 9))))
                .make();
        byAllele("three-alleles-three-values", triallelic, engine);
        altByAllele("three-alleles-three-values", triallelic, engine);

        // A genotype whose list is SHORT: the zip stops, and the later alleles get nothing from it.
        final VariantContext ragged = new VariantContextBuilder(triallelic)
                .genotypes(List.of(
                        genotype("T1", new int[] {80, 20, 5}, List.of(1)),
                        genotype("T2", new int[] {70, 30, 10}, List.of(4, 5, 6)),
                        genotype("N1", new int[] {99, 1, 0}, List.of(7, 8, 9))))
                .make();
        byAllele("short-list", ragged, engine);
        altByAllele("short-list", ragged, engine);

        // A genotype whose list is LONG: the tail is dropped.
        final VariantContext overlong = new VariantContextBuilder(triallelic)
                .genotypes(List.of(
                        genotype("T1", new int[] {80, 20, 5}, List.of(1, 2, 3, 4, 5)),
                        genotype("N1", new int[] {99, 1, 0}, List.of(7, 8, 9))))
                .make();
        byAllele("long-list", overlong, engine);

        // No genotype passes the precondition at all.
        byAllele("no-tumor", new VariantContextBuilder(triallelic)
                .genotypes(List.of(genotype("N1", new int[] {99, 1, 0}, List.of(7, 8, 9)))).make(), engine);

        // The weighted median, over the shapes that decide which element it lands on.
        median("even-split", List.of(pair(10, 0.1), pair(10, 0.9)));
        median("one-dominant", List.of(pair(1, 0.1), pair(99, 0.9)));
        median("out-of-order", List.of(pair(5, 0.9), pair(5, 0.2), pair(5, 0.5)));
        median("all-equal", List.of(pair(4, 0.5), pair(4, 0.5), pair(4, 0.5)));
        median("zero-depth", List.of(pair(0, 0.3), pair(0, 0.7)));
        median("single", List.of(pair(7, 0.42)));
        median("empty", List.of());

        // sumADsOverSamples, over the three flag combinations and a short AD array.
        ads("tumor-only", triallelic, engine, true, false);
        ads("normal-only", triallelic, engine, false, true);
        ads("both", triallelic, engine, true, true);
        ads("neither", triallelic, engine, false, false);
        ads("short-ad", new VariantContextBuilder(triallelic)
                .genotypes(List.of(genotype("T1", new int[] {80, 20}, List.of(1, 2, 3)),
                        genotype("N1", new int[] {99, 1, 0}, List.of(7, 8, 9)))).make(),
                engine, true, false);

        // The filter itself: its identity, and its probabilities on records that differ only in TLOD.
        final TumorEvidenceFilter filter = new TumorEvidenceFilter();
        System.out.printf("name\ttumor-evidence\t%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"));
        evidence("strong", filter, engine, triallelic);
        evidence("weak", filter, engine,
                new VariantContextBuilder(triallelic).attribute("TLOD", List.of(1.0, 0.5)).make());
        evidence("negative", filter, engine,
                new VariantContextBuilder(triallelic).attribute("TLOD", List.of(-3.0, -0.1)).make());
        // No TLOD at all: an empty list rather than a probability.
        evidence("no-tlod", filter, engine,
                new VariantContextBuilder(triallelic).rmAttribute("TLOD").make());
        // A TLOD list shorter than the alternate alleles.
        evidence("short-tlod", filter, engine,
                new VariantContextBuilder(triallelic).attribute("TLOD", List.of(20.0)).make());
    }

    static Genotype genotype(final String sample, final int[] ad, final List<Integer> perAllele) {
        return new GenotypeBuilder(sample, List.of(Allele.REF_A, Allele.ALT_C))
                .AD(ad).attribute("VALUES", perAllele).make();
    }

    static ImmutablePair<Integer, Double> pair(final int depth, final double posterior) {
        return ImmutablePair.of(depth, posterior);
    }

    /** `getDataByAllele`, whose map includes the reference allele. */
    static void byAllele(final String label, final VariantContext vc,
                         final Mutect2FilteringEngine engine) {
        try {
            final Map<Allele, List<Integer>> data = Mutect2AlleleFilter.getDataByAllele(vc,
                    engine::isTumor, g -> values(g), engine);
            print("byallele", label, data);
        } catch (final Exception e) {
            System.out.printf("error\tbyallele-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** `getAltDataByAllele`, whose map does not. */
    static void altByAllele(final String label, final VariantContext vc,
                            final Mutect2FilteringEngine engine) {
        try {
            final Map<Allele, List<Integer>> data = Mutect2AlleleFilter.getAltDataByAllele(vc,
                    engine::isTumor, g -> values(g), engine);
            print("altbyallele", label, data);
        } catch (final Exception e) {
            System.out.printf("error\taltbyallele-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    @SuppressWarnings("unchecked")
    static List<Integer> values(final Genotype g) {
        return (List<Integer>) g.getExtendedAttribute("VALUES");
    }

    static void print(final String kind, final String label, final Map<Allele, List<Integer>> data) {
        for (final Map.Entry<Allele, List<Integer>> entry : data.entrySet()) {
            System.out.printf("%s\t%s\t%s=%s%n", kind, label,
                    entry.getKey().getDisplayString(), entry.getValue());
        }
    }

    /** `weightedMedianPosteriorProbability`, which is protected and static, and sorts in place. */
    static void median(final String label, final List<ImmutablePair<Integer, Double>> input) {
        final List<ImmutablePair<Integer, Double>> mutable = new ArrayList<>(input);
        try {
            final Method method = Class.forName(
                    "org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter")
                    .getDeclaredMethod("weightedMedianPosteriorProbability", List.class);
            method.setAccessible(true);
            final double answer = (double) method.invoke(null, mutable);
            System.out.printf("median\t%s\t%s%n", label, Double.toString(answer));
            System.out.printf("sorted\t%s\t%s%n", label, mutable);
        } catch (final Exception e) {
            final Throwable cause = e.getCause() == null ? e : e.getCause();
            System.out.printf("error\tmedian-%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
        }
    }

    static void ads(final String label, final VariantContext vc, final Mutect2FilteringEngine engine,
                    final boolean includeTumor, final boolean includeNormal) {
        try {
            System.out.printf("ads\t%s\t%s%n", label,
                    Arrays.toString(engine.sumADsOverSamples(vc, includeTumor, includeNormal)));
        } catch (final Exception e) {
            System.out.printf("error\tads-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void evidence(final String label, final TumorEvidenceFilter filter,
                         final Mutect2FilteringEngine engine, final VariantContext vc) {
        try {
            System.out.printf("evidence\t%s\t%s%n", label,
                    filter.errorProbabilities(vc, engine, null));
        } catch (final Exception e) {
            System.out.printf("error\tevidence-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
