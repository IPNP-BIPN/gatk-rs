/*
 * `PolymeraseSlippageFilter`, taken from the reference.
 *
 * The seventh of the ten filters the `filter-mutect-calls` golden needs, and the first to reach
 * `Beta.regularizedBeta` for its own likelihood rather than through a binomial CDF. Six behaviours
 * this is built to catch.
 *
 *   - THE BETA IS ASKED ABOUT `ADs[1]`, NOT ABOUT THE ALTERNATE COUNT.
 *     `regularizedBeta(slippageRate, ADs[1] + 1, ADs[0] + 1)` reads the depth of the FIRST
 *     alternate allele, while the somatic likelihood beside it uses `depth` and `depth - ADs[0]`,
 *     which sums EVERY alternate. On a triallelic record the two halves of one log-odds are about
 *     different alleles;
 *   - THE PRIOR IS `getLogPriorOfSomaticVariant(vc, 0)`, WITH THE INDEX HARD-CODED, so whichever
 *     allele slipped, the prior comes from alternate zero's indel length;
 *   - `getLogPriorOfSomaticVariant` INSERTS BEFORE IT READS, so the dump calls the filter twice on
 *     one record with one engine, to pin whether two identical calls are two identical answers;
 *   - `RPA` IS PARSED WITH `Integer.parseInt(String.valueOf(o))`, so a non-integer entry is a
 *     `NumberFormatException` out of a filter rather than a skip. `rpa.length < 2` is the only
 *     length guard there is;
 *   - THE GATE IS `ru.length() * rpa[0] >= minSlippageLength && Math.abs(rpa[0] - rpa[1]) == 1`.
 *     An empty `RU` makes the base count zero, and exactly one slip is required in either
 *     direction: a two-repeat contraction is not filtered at all;
 *   - AND THE `MaxCountExceededException` FALLBACK CANNOT BE REACHED, `regularizedBeta(x, a, b)`
 *     passing `Integer.MAX_VALUE` iterations.
 *
 * Output:
 *
 *     default\t<argument>\t<value>
 *     name\tslippage\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     prob\t<label>\t<one error probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PolymeraseSlippageFilterDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.PolymeraseSlippageFilter;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class PolymeraseSlippageFilterDump {

    static final Allele REF = Allele.create("AA", true);
    static final Allele DELETION = Allele.create("A", false);
    static final Allele INSERTION = Allele.create("AAA", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# PolymeraseSlippageFilterDump: the short-tandem-repeat slippage filter");

        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        System.out.printf("default\tminSlippageLength\t%d%n", arguments.minSlippageLength);
        System.out.printf("default\tslippageRate\t%s%n", Double.toString(arguments.slippageRate));

        final PolymeraseSlippageFilter filter =
                new PolymeraseSlippageFilter(arguments.minSlippageLength, arguments.slippageRate);
        System.out.printf("name\tslippage\t%s,%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"), "RPA,RU");

        // A ten-base reference repeat contracted by one: the gate opens.
        final VariantContext base = new VariantContextBuilder("dump", "chr1", 100, 101,
                List.of(REF, DELETION))
                .attribute("RPA", List.of(10, 9))
                .attribute("RU", "A")
                .genotypes(List.of(genotype("T1", new int[] {80, 20}), genotype("N1", new int[] {90, 1})))
                .make();
        prob("contracted-by-one", base);

        // The same record twice through ONE engine: does the model's insert change the answer?
        twice("called-twice", base);

        // An expansion rather than a contraction: `numPCRSlips` is negative and `abs` accepts it.
        prob("expanded-by-one", new VariantContextBuilder(base)
                .attribute("RPA", List.of(9, 10)).make());

        // Two slips, which the `== 1` refuses.
        prob("contracted-by-two", new VariantContextBuilder(base)
                .attribute("RPA", List.of(10, 8)).make());

        // A repeat unit longer than one base multiplies the base count.
        prob("two-base-repeat-unit", new VariantContextBuilder(base)
                .attribute("RPA", List.of(5, 4)).attribute("RU", "AT").make());

        // Below the minimum slippage length, and exactly at it.
        prob("below-the-minimum", new VariantContextBuilder(base)
                .attribute("RPA", List.of(7, 6)).make());
        prob("at-the-minimum", new VariantContextBuilder(base)
                .attribute("RPA", List.of(8, 7)).make());

        // An empty repeat unit: the base count is zero however long the repeat is.
        prob("empty-repeat-unit", new VariantContextBuilder(base).attribute("RU", "").make());

        // One entry in RPA, which is the only length guard the filter has.
        prob("one-entry-in-rpa", new VariantContextBuilder(base)
                .attribute("RPA", List.of(10)).make());

        // A non-integer RPA entry, parsed with Integer.parseInt.
        prob("non-integer-rpa", new VariantContextBuilder(base)
                .attribute("RPA", List.of("ten", "nine")).make());
        prob("decimal-rpa", new VariantContextBuilder(base)
                .attribute("RPA", List.of(10.0, 9.0)).make());

        // The missing required annotations, which this base class answers with 0.0 per allele.
        prob("no-rpa", new VariantContextBuilder(base).rmAttribute("RPA").make());
        prob("no-ru", new VariantContextBuilder(base).rmAttribute("RU").make());

        // The depths the beta reads: no alternate reads, every read alternate, and a deep site.
        prob("no-alt-reads", withTumor(base, new int[] {100, 0}));
        prob("all-alt-reads", withTumor(base, new int[] {0, 100}));
        prob("shallow", withTumor(base, new int[] {4, 2}));
        prob("deep", withTumor(base, new int[] {800, 200}));

        // A triallelic record: the beta reads ADs[1] alone while the likelihood sums both
        // alternates, and the prior is alternate ZERO's indel length either way.
        final VariantContext triallelic = new VariantContextBuilder("dump", "chr1", 100, 101,
                List.of(REF, DELETION, INSERTION))
                .attribute("RPA", List.of(10, 9, 11))
                .attribute("RU", "A")
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 40}),
                        genotype("N1", new int[] {90, 1, 0})))
                .make();
        prob("triallelic", triallelic);

        // A slippage rate at each end of its range.
        prob("slippage-rate-one", base, 1.0);
        prob("slippage-rate-zero", base, 0.0);
        prob("slippage-rate-tiny", base, 1.0e-12);
    }

    static Genotype genotype(final String sample, final int[] ad) {
        return new GenotypeBuilder(sample, List.of(REF, DELETION)).AD(ad).make();
    }

    static VariantContext withTumor(final VariantContext vc, final int[] ad) {
        return new VariantContextBuilder(vc)
                .genotypes(List.of(genotype("T1", ad), genotype("N1", new int[] {90, 1}))).make();
    }

    /** A fresh engine per case, so one case's clustering model cannot reach the next. */
    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static void prob(final String label, final VariantContext vc) {
        prob(label, vc, new M2FiltersArgumentCollection().slippageRate);
    }

    static void prob(final String label, final VariantContext vc, final double slippageRate) {
        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        final PolymeraseSlippageFilter filter =
                new PolymeraseSlippageFilter(arguments.minSlippageLength, slippageRate);
        try {
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine(), null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** The same record twice through one engine, whose clustering model the first call touched. */
    static void twice(final String label, final VariantContext vc) {
        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        final PolymeraseSlippageFilter filter =
                new PolymeraseSlippageFilter(arguments.minSlippageLength, arguments.slippageRate);
        final Mutect2FilteringEngine engine = engine();
        try {
            System.out.printf("prob\t%s-first\t%s%n", label, filter.errorProbabilities(vc, engine, null));
            System.out.printf("prob\t%s-second\t%s%n", label, filter.errorProbabilities(vc, engine, null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
