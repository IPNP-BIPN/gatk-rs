/*
 * `NormalArtifactFilter`, and the commons-math it reaches through, taken from the reference.
 *
 * The second of the ten filters the `filter-mutect-calls` golden needs, and the first that reaches
 * `Beta.regularizedBeta` -- a modified Lentz continued fraction that nothing here has ported.
 * Seven behaviours this is built to catch.
 *
 *   - THE FILTER READS ONE INDEX TWO DIFFERENT WAYS. `indexOfMaxTumorLod` is the ALTERNATE allele's
 *     index, so the depths are read at `indexOfMaxTumorLod + 1`, the AD array counting the
 *     reference, and the normal-artifact log odds at `indexOfMaxTumorLod`, NALOD not counting it;
 *   - THE RATIO GATE COMPARES ALLELE FRACTIONS, AND A TUMOUR WITH NO DEPTH IS NaN RATHER THAN A
 *     REFUSAL. `normalDepth == 0 ? 0 : ...` guards the normal side alone; the tumour side is a bare
 *     division, so an empty tumour makes the gate `x < NaN`, which is false, and the record falls
 *     through to the arithmetic instead of returning zero;
 *   - A MISSING MBQ IS NOT THE IMPUTED 30, IT IS AN IndexOutOfBoundsException.
 *     `getAttributeAsIntList(key, 30)` answers an EMPTY LIST for an absent key -- the default
 *     applies to a null or "." ELEMENT of a present list -- and the filter calls `.get(0)` on it;
 *   - `cumulativeProbability(normalAltDepth - 1)` IS EVALUATED AT -1 when the normal has no alt
 *     read, which commons-math answers `0.0` rather than refusing, so the p-value is exactly 1 and
 *     the filter falls through to the posterior. At the other end `x >= trials` answers `1.0`, so
 *     the p-value is 0, below any threshold, and the filter returns a hard `1.0`;
 *   - `Mutect2VariantFilter.errorProbabilities` ANSWERS 0.0 PER ALLELE FOR A MISSING ANNOTATION,
 *     where the per-allele class answers an EMPTY LIST. The two base classes disagree, and only the
 *     empty list is dropped by `ErrorProbabilities` rather than counted;
 *   - `regularizedBeta`'s SECOND CONJUNCT IS IMPLIED BY ITS FIRST IN THE REALS:
 *     `1 - x <= (b + 1) / (2 + b + a)` rearranges to `x >= (a + 1) / (2 + b + a)`, which the strict
 *     `x > (a + 1) / (2 + b + a)` already asserts. The two can only disagree by a rounding, so the
 *     grid walks the boundary itself;
 *   - AND THE CONTINUED FRACTION'S EPSILON IS `Beta`'s 1E-14, NOT `ContinuedFraction`'s OWN 10E-9,
 *     with `|v| <= 1e-50` as its zero test rather than `v == 0.0`.
 *
 * Output:
 *
 *     regbeta\t<label>\t<x>,<a>,<b>=<Beta.regularizedBeta>
 *     cdf\t<label>\t<trials>,<p>,<x>=<BinomialDistribution.cumulativeProbability>
 *     qualerr\t<label>\t<qual>=<QualityUtils.qualToErrorProb>
 *     prior\tlog-prior-variant-versus-artifact\t<the clustering model's prior>
 *     posterior\t<label>\t<negative log odds>=<posteriorProbabilityOfNormalArtifact>
 *     name\tnormal-artifact\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     prob\t<label>\t<one error probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: NormalArtifactFilterDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.special.Beta;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NormalArtifactFilter;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class NormalArtifactFilterDump {

    /** `M2FiltersArgumentCollection.DEFAULT_NORMAL_P_VALUE_THRESHOLD`, which is private. */
    static final double DEFAULT_NORMAL_P_VALUE_THRESHOLD = 0.001;

    public static void main(final String[] args) throws Exception {
        System.out.println("# NormalArtifactFilterDump: the normal-artifact filter and its binomial p-value");

        // `Beta.regularizedBeta`, over both arrangements of the symmetry branch and the guards.
        // The direct branch: x at or below (a + 1) / (2 + b + a).
        regbeta("direct-small", 0.01, 2.0, 5.0);
        regbeta("direct-mid", 0.2, 2.0, 5.0);
        regbeta("direct-equal", 1.0 / 3.0, 2.0, 5.0);
        // The symmetry branch, which recurses on 1 - x with a and b swapped.
        regbeta("symmetric-half", 0.5, 2.0, 2.0);
        regbeta("symmetric-large", 0.99, 2.0, 5.0);
        regbeta("symmetric-just-over", 0.4, 2.0, 5.0);
        // The boundary itself, where the two conjuncts can only disagree by a rounding.
        regbeta("boundary-exact", 3.0 / 9.0, 2.0, 5.0);
        regbeta("boundary-below", Math.nextDown(3.0 / 9.0), 2.0, 5.0);
        regbeta("boundary-above", Math.nextUp(3.0 / 9.0), 2.0, 5.0);
        // The shapes a binomial CDF actually asks for: a small p, a large a, a b of one.
        regbeta("binomial-like", 0.001, 1.0, 100.0);
        regbeta("binomial-deep", 0.001, 11.0, 90.0);
        regbeta("b-is-one", 0.25, 4.0, 1.0);
        regbeta("a-is-one", 0.25, 1.0, 4.0);
        regbeta("both-large", 0.5, 500.0, 500.0);
        regbeta("both-tiny", 0.5, 0.001, 0.001);
        // The endpoints, where the logs go to infinity rather than the guards catching them.
        regbeta("x-zero", 0.0, 2.0, 5.0);
        regbeta("x-one", 1.0, 2.0, 5.0);
        // The guards, every one of which answers NaN rather than throwing.
        regbeta("x-below-zero", -0.5, 2.0, 5.0);
        regbeta("x-above-one", 1.5, 2.0, 5.0);
        regbeta("a-zero", 0.5, 0.0, 5.0);
        regbeta("b-negative", 0.5, 2.0, -1.0);
        regbeta("x-nan", Double.NaN, 2.0, 5.0);

        // `BinomialDistribution.cumulativeProbability`, over the three arms and the shapes the
        // filter reaches. The constructor takes a null RandomGenerator, which it accepts.
        cdf("no-alt-read", 100, QualityUtils.qualToErrorProb(30.0), -1);
        cdf("one-alt-read", 100, QualityUtils.qualToErrorProb(30.0), 0);
        cdf("ten-alt-reads", 100, QualityUtils.qualToErrorProb(30.0), 9);
        cdf("every-read-alt", 100, QualityUtils.qualToErrorProb(30.0), 99);
        cdf("past-the-trials", 100, QualityUtils.qualToErrorProb(30.0), 100);
        cdf("no-trials", 0, QualityUtils.qualToErrorProb(30.0), 0);
        cdf("low-quality", 100, QualityUtils.qualToErrorProb(2.0), 9);
        cdf("high-quality", 100, QualityUtils.qualToErrorProb(60.0), 9);
        cdf("shallow-normal", 10, QualityUtils.qualToErrorProb(30.0), 4);
        cdf("half", 20, 0.5, 9);
        cdf("p-zero", 20, 0.0, 9);
        cdf("p-one", 20, 1.0, 9);

        // `QualityUtils.qualToErrorProb(double)`, which the int median base quality widens into.
        for (final int qual : new int[] {0, 2, 20, 30, 60, 93}) {
            System.out.printf("qualerr\tq%d\t%d=%s%n", qual, qual,
                    Double.toString(QualityUtils.qualToErrorProb((double) qual)));
        }

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(
                new M2FiltersArgumentCollection(), header, new File("no-such-stats-file.tsv"));

        System.out.printf("prior\tlog-prior-variant-versus-artifact\t%s%n",
                Double.toString(engine.getSomaticClusteringModel().getLogPriorOfVariantVersusArtifact()));
        // `posteriorProbabilityOfNormalArtifact`, which is the prior above against one log odds.
        for (final double nalod : new double[] {-10.0, -2.0, -0.5, 0.0, 0.5, 2.0, 10.0}) {
            final double negativeLogOdds = -(nalod * Math.log(10.0));
            System.out.printf("posterior\tnalod%s\t%s=%s%n", Double.toString(nalod),
                    Double.toString(negativeLogOdds),
                    Double.toString(engine.posteriorProbabilityOfNormalArtifact(negativeLogOdds)));
        }

        final NormalArtifactFilter filter = new NormalArtifactFilter(DEFAULT_NORMAL_P_VALUE_THRESHOLD);
        System.out.printf("name\tnormal-artifact\t%s,%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"), "NALOD,TLOD");

        // A triallelic record whose normal carries enough of the alternate allele to pass the gate.
        final VariantContext base = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G))
                .attribute("TLOD", List.of(20.0, 6.0))
                .attribute("NALOD", List.of(2.0, 0.5))
                .attribute("MBQ", List.of(30, 30, 30))
                .genotypes(List.of(
                        genotype("T1", new int[] {80, 20, 5}),
                        genotype("N1", new int[] {90, 10, 2})))
                .make();
        prob("normal-carries-the-allele", filter, engine, base);

        // The second alternate allele wins the tumour log odds: the depths move by two, the NALOD
        // by one.
        prob("second-allele", filter, engine,
                new VariantContextBuilder(base).attribute("TLOD", List.of(6.0, 20.0)).make());

        // The ratio gate: a normal with no alternate read at all is a fraction of zero.
        prob("normal-is-clean", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}),
                        genotype("N1", new int[] {90, 0, 0}))).make());

        // A normal whose fraction is exactly a tenth of the tumour's, which the strict `<` keeps.
        prob("normal-at-the-ratio", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}),
                        genotype("N1", new int[] {80, 2, 0}))).make());

        // No normal sample at all: the normal depth is zero, so the fraction is the guarded zero.
        prob("no-normal-sample", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}))).make());

        // No TUMOUR depth: the gate becomes `0 < NaN`, which is false, and the record falls through.
        prob("no-tumor-depth", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {0, 0, 0}),
                        genotype("N1", new int[] {90, 10, 2}))).make());

        // And with no tumour depth and a clean normal, the CDF is asked for -1.
        prob("no-tumor-depth-clean-normal", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {0, 0, 0}),
                        genotype("N1", new int[] {90, 0, 0}))).make());

        // Every read in the normal supports the allele: `x >= trials`, so the p-value is zero.
        prob("normal-is-all-alt", filter, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}),
                        genotype("N1", new int[] {0, 20, 0}))).make());

        // The median base quality decides the p-value, and the threshold turns it into a hard 1.0.
        prob("low-median-base-quality", filter, engine,
                new VariantContextBuilder(base).attribute("MBQ", List.of(2, 2, 2)).make());
        prob("high-median-base-quality", filter, engine,
                new VariantContextBuilder(base).attribute("MBQ", List.of(60, 60, 60)).make());
        // A missing MBQ is the empty list, not the imputed 30.
        prob("no-mbq", filter, engine, new VariantContextBuilder(base).rmAttribute("MBQ").make());

        // A NALOD saying the normal looks like an artifact, and one saying it does not. Both carry
        // the low median base quality, because at MBQ 30 the p-value fires first and the posterior
        // never reaches the answer.
        prob("negative-nalod", filter, engine, new VariantContextBuilder(base)
                .attribute("MBQ", List.of(2, 2, 2)).attribute("NALOD", List.of(-5.0, -1.0)).make());
        prob("large-nalod", filter, engine, new VariantContextBuilder(base)
                .attribute("MBQ", List.of(2, 2, 2)).attribute("NALOD", List.of(20.0, 20.0)).make());

        // The missing required annotations, which this base class answers with 0.0 per allele
        // rather than the empty list the per-allele base class answers.
        prob("no-nalod", filter, engine, new VariantContextBuilder(base).rmAttribute("NALOD").make());
        prob("no-tlod", filter, engine, new VariantContextBuilder(base).rmAttribute("TLOD").make());

        // A biallelic record, where the alternate index is the only one there is.
        prob("biallelic", filter, engine, new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C))
                .attribute("TLOD", List.of(20.0))
                .attribute("NALOD", List.of(2.0))
                .attribute("MBQ", List.of(30, 30))
                .genotypes(List.of(genotype("T1", new int[] {80, 20}),
                        genotype("N1", new int[] {90, 10})))
                .make());
    }

    static Genotype genotype(final String sample, final int[] ad) {
        return new GenotypeBuilder(sample, List.of(Allele.REF_A, Allele.ALT_C)).AD(ad).make();
    }

    static void regbeta(final String label, final double x, final double a, final double b) {
        try {
            System.out.printf("regbeta\t%s\t%s,%s,%s=%s%n", label, Double.toString(x),
                    Double.toString(a), Double.toString(b),
                    Double.toString(Beta.regularizedBeta(x, a, b)));
        } catch (final Exception e) {
            System.out.printf("error\tregbeta-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void cdf(final String label, final int trials, final double p, final int x) {
        try {
            System.out.printf("cdf\t%s\t%d,%s,%d=%s%n", label, trials, Double.toString(p), x,
                    Double.toString(new BinomialDistribution(null, trials, p).cumulativeProbability(x)));
        } catch (final Exception e) {
            System.out.printf("error\tcdf-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void prob(final String label, final NormalArtifactFilter filter,
                     final Mutect2FilteringEngine engine, final VariantContext vc) {
        try {
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine, null));
        } catch (final Exception e) {
            System.out.printf("error\tprob-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
