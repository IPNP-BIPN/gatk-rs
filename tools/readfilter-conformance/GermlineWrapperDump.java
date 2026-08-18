/*
 * `GermlineFilter.calculateErrorProbability` and the two helpers it calls, taken from the reference.
 *
 * The core `germlineProbability` is measured by the `germline-probability` suite; this is the layer
 * that reads the record and decides what to hand it. Six behaviours it is built to catch.
 *
 *   - TWO EARLY RETURNS BRACKET THE POPULATION FREQUENCY: below 1e-10 the answer is 0 and above
 *     1 - 1e-10 it is 1, and the frequency is `Math.pow(10, -POPAF[maxLodIndex])`, so a POPAF of 0
 *     means an allele fixed in the population and answers a hard one;
 *   - `maxLodIndex` IS CHOSEN FROM `TLOD` AND THEN INDEXES THREE ARRAYS with three conventions:
 *     `POPAF` per alternate, the depth array with `+ 1`, and the weighted allele fractions per
 *     alternate;
 *   - A TOTAL DEPTH OF ZERO ANSWERS 0, per a comment about GGA mode;
 *   - THE GERMLINE HOM-ALT HYPOTHESIS IS SWITCHED OFF BY A VALUE, NOT A FLAG: an alternate allele
 *     fraction below 0.9 sets its log odds to negative infinity;
 *   - THE NORMAL'S LOG ODDS ENTER NEGATED, and a record with no `NLOD` uses 0 rather than skipping
 *     the normal;
 *   - AND WITH NO TUMOUR SEGMENTATION TABLE every sample's minor allele fraction is 0.5, so the
 *     helper is a depth-weighted average of one repeated constant whose denominator is the
 *     tumour-only depth sum.
 *
 * Output:
 *
 *     prob\t<label>\t<one probability per alternate allele>
 *     weightedaf\t<label>\t<weightedAverageOfTumorAFs>
 *     maf\t<label>\t<computeMinorAlleleFraction>
 *     name\tgermline\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GermlineWrapperDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.GermlineFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class GermlineWrapperDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT_C = Allele.create("C", false);
    static final Allele ALT_G = Allele.create("G", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# GermlineWrapperDump: the germline filter's wrapper");

        final GermlineFilter filter = new GermlineFilter(List.of());
        System.out.printf("name\tgermline\t%s,%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"), "TLOD,POPAF");

        // A biallelic record with a rare population allele and a low tumour allele fraction.
        prob("rare-allele", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2}));
        // A common population allele.
        prob("common-allele", record(List.of(REF, ALT_C), List.of(20.0), List.of(1.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2}));
        // Fixed in the population: `POPAF` of zero is a frequency of one.
        prob("fixed-in-the-population", record(List.of(REF, ALT_C), List.of(20.0), List.of(0.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2}));
        // Vanishingly rare: below the epsilon, the answer is a hard zero.
        prob("below-the-epsilon", record(List.of(REF, ALT_C), List.of(20.0), List.of(400.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2}));

        // The hom-alt hypothesis, switched off below an allele fraction of 0.9 and on above it.
        prob("hom-alt-off", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {10, 80}, new double[] {0.89}));
        prob("hom-alt-on", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {10, 80}, new double[] {0.9}));

        // The normal's log odds, present and absent.
        prob("normal-says-germline", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-8.0), new int[] {80, 20}, new double[] {0.2}));
        prob("normal-says-somatic", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(8.0), new int[] {80, 20}, new double[] {0.2}));
        prob("no-nlod", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                null, new int[] {80, 20}, new double[] {0.2}));

        // Two alternates, where the second wins the tumour log odds.
        prob("second-allele-wins", record(List.of(REF, ALT_C, ALT_G), List.of(6.0, 20.0),
                List.of(6.0, 2.0), List.of(-4.0, -4.0), new int[] {60, 20, 20},
                new double[] {0.2, 0.2}));

        // No depth at all.
        prob("no-depth", record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {0, 0}, new double[] {0.0}));

        // The missing required annotations.
        prob("no-tlod", stripped(List.of(REF, ALT_C), "TLOD"));
        prob("no-popaf", stripped(List.of(REF, ALT_C), "POPAF"));

        // The two helpers, on the record every case above is a variation of.
        final VariantContext base = record(List.of(REF, ALT_C), List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2});
        helpers("one-tumour", base);
        // Two tumour samples with different depths and fractions: both helpers weight by depth.
        helpers("two-tumours", twoTumours());
    }

    static Genotype tumour(final String sample, final int[] ad, final double[] alleleFractions) {
        return new GenotypeBuilder(sample, List.of(REF, REF)).AD(ad)
                .attribute("AF", alleleFractions).make();
    }

    static VariantContext record(final List<Allele> alleles, final List<Double> tumorLog10Odds,
                                 final List<Double> populationAf, final List<Double> normalLog10Odds,
                                 final int[] tumorDepths, final double[] alleleFractions) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", 100, 100, alleles)
                        .attribute("TLOD", tumorLog10Odds)
                        .attribute("POPAF", populationAf)
                        .genotypes(List.of(tumour("T1", tumorDepths, alleleFractions),
                                new GenotypeBuilder("N1", List.of(REF, REF))
                                        .AD(new int[] {90, 1}).make()));
        if (normalLog10Odds != null) {
            builder.attribute("NLOD", normalLog10Odds);
        }
        return builder.make();
    }

    static VariantContext stripped(final List<Allele> alleles, final String key) {
        return new VariantContextBuilder(record(alleles, List.of(20.0), List.of(6.0),
                List.of(-4.0), new int[] {80, 20}, new double[] {0.2})).rmAttribute(key).make();
    }

    static VariantContext twoTumours() {
        return new VariantContextBuilder("dump", "chr1", 100, 100, List.of(REF, ALT_C))
                .attribute("TLOD", List.of(20.0))
                .attribute("POPAF", List.of(6.0))
                .attribute("NLOD", List.of(-4.0))
                .genotypes(List.of(tumour("T1", new int[] {80, 20}, new double[] {0.2}),
                        tumour("T2", new int[] {20, 30}, new double[] {0.6}),
                        new GenotypeBuilder("N1", List.of(REF, REF)).AD(new int[] {90, 1}).make()))
                .make();
    }

    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "T2", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static void prob(final String label, final VariantContext vc) {
        final GermlineFilter filter = new GermlineFilter(List.of());
        try {
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine(), null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** `weightedAverageOfTumorAFs`, which is public, and the private minor allele fraction. */
    static void helpers(final String label, final VariantContext vc) throws Exception {
        final Mutect2FilteringEngine engine = engine();
        System.out.printf("weightedaf\t%s\t%s%n", label,
                Arrays.toString(engine.weightedAverageOfTumorAFs(vc)));
        final int[] alleleCounts = engine.sumADsOverSamples(vc, true, false);
        final Method method = GermlineFilter.class.getDeclaredMethod("computeMinorAlleleFraction",
                VariantContext.class, Mutect2FilteringEngine.class, int[].class);
        method.setAccessible(true);
        final double maf = (double) method.invoke(new GermlineFilter(List.of()), vc, engine, alleleCounts);
        System.out.printf("maf\t%s\t%s%n", label, Double.toString(maf));
    }
}
