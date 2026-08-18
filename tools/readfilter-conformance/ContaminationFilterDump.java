/*
 * `ContaminationFilter`, taken from the reference.
 *
 * The eighth of the ten filters the `filter-mutect-calls` golden needs, and one of only two in the
 * `NON_SOMATIC` error type. Six behaviours this is built to catch.
 *
 *   - THIS FILTER CAN ANSWER NaN. `depthsAndPosteriorsPerAllele` is sized from the record's
 *     alleles, but every loop that fills it runs over `alleleFrequencies.length`, which is
 *     `POPAF`'s. An allele the annotation does not cover keeps an EMPTY list, and an empty list is
 *     `Double.NaN`, which `roundFinitePrecisionErrors` passes through unchanged because
 *     `Math.min`/`Math.max` propagate NaN;
 *   - A `POPAF` LONGER THAN THE RECORD IS AN ArrayIndexOutOfBoundsException, the same loop indexing
 *     `altADs`, which is sized from the genotype's `AD`;
 *   - THE CONTAMINATION IS CLAMPED ON BOTH SIDES: `Math.max(0, Math.min(fromFile, 1 - EPSILON))`,
 *     where `EPSILON` is an INSTANCE field of `1.0e-10` and the comment says the clamp exists "to
 *     handle file with contamination == 1";
 *   - THE TWO CONTAMINATION HYPOTHESES ARE COMPARED BY MAXIMUM AND LOGGED ONCE:
 *     `log(max(singleContaminant, manyContaminant))` picks the larger rather than summing, so which
 *     hypothesis wins can change with the depth alone;
 *   - THE PRIOR'S ALLELE INDEX IS THE LOOP'S, unlike `PolymeraseSlippageFilter`, which hard-codes
 *     zero in the same call: each alternate gets its own indel length;
 *   - AND THE ANSWER IS A WEIGHTED MEDIAN ACROSS SAMPLES, weighted by each sample's alternate depth,
 *     so two tumour samples do not average.
 *
 * The dump passes no contamination tables, so every sample takes `defaultContamination`; a table
 * naming a sample twice is an `IllegalStateException` out of `Collectors.toMap`, which needs files
 * and is out of scope here.
 *
 * Output:
 *
 *     default\tcontaminationEstimate\t<value>
 *     name\tcontamination\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     prob\t<label>\t<one error probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ContaminationFilterDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ContaminationFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class ContaminationFilterDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# ContaminationFilterDump: the cross-sample contamination filter");

        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        System.out.printf("default\tcontaminationEstimate\t%s%n",
                Double.toString(arguments.contaminationEstimate));

        final ContaminationFilter filter = new ContaminationFilter(List.of(), 0.05);
        System.out.printf("name\tcontamination\t%s,%s,%s,%s%n", filter.filterName(), filter.errorType(),
                filter.phredScaledPosteriorAnnotationName().orElse("none"), "POPAF");

        // A triallelic record with one tumour and one normal, and a population allele frequency for
        // each alternate: POPAF is -log10, so 2.0 is a frequency of 0.01.
        final VariantContext base = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G))
                .attribute("POPAF", List.of(2.0, 3.0))
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}),
                        genotype("N1", new int[] {90, 1, 0})))
                .make();

        // The contamination estimate, across its range and both ends of its clamp.
        prob("contamination-default", base, arguments.contaminationEstimate);
        prob("contamination-five-percent", base, 0.05);
        prob("contamination-half", base, 0.5);
        prob("contamination-one", base, 1.0);
        prob("contamination-above-one", base, 2.0);
        prob("contamination-negative", base, -0.5);

        // The population frequencies: common, vanishingly rare, and certain.
        prob("common-allele", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(0.5, 0.5)).make(), 0.05);
        prob("rare-allele", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(9.0, 9.0)).make(), 0.05);
        prob("frequency-of-one", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(0.0, 0.0)).make(), 0.05);
        prob("infinite-popaf", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY))
                .make(), 0.05);

        // A POPAF shorter than the record's alternates: the second allele keeps an empty list.
        prob("short-popaf", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(2.0)).make(), 0.05);
        // And one longer than the genotype's AD.
        prob("long-popaf", new VariantContextBuilder(base)
                .attribute("POPAF", List.of(2.0, 3.0, 4.0)).make(), 0.05);
        // No POPAF at all: the required-annotation check answers an empty list.
        prob("no-popaf", new VariantContextBuilder(base).rmAttribute("POPAF").make(), 0.05);

        // Only a normal sample: every allele's list is empty, so every allele is NaN.
        prob("normal-only", new VariantContextBuilder(base)
                .genotypes(List.of(genotype("N1", new int[] {90, 1, 0}))).make(), 0.05);

        // Two tumour samples, whose posteriors are combined by a depth-weighted median rather than
        // averaged.
        prob("two-tumours", new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}),
                        genotype("T2", new int[] {40, 1, 30}),
                        genotype("N1", new int[] {90, 1, 0}))).make(), 0.05);

        // The depths the binomials read.
        prob("no-alt-reads", withTumor(base, new int[] {100, 0, 0}), 0.05);
        prob("all-alt-reads", withTumor(base, new int[] {0, 100, 0}), 0.05);
        prob("shallow", withTumor(base, new int[] {4, 2, 1}), 0.05);
        prob("deep", withTumor(base, new int[] {800, 200, 50}), 0.05);
        // A depth at which the many-contaminant hypothesis overtakes the single one.
        prob("one-alt-read-deep", withTumor(base, new int[] {999, 1, 0}), 0.05);
    }

    static Genotype genotype(final String sample, final int[] ad) {
        return new GenotypeBuilder(sample, List.of(Allele.REF_A, Allele.ALT_C)).AD(ad).make();
    }

    static VariantContext withTumor(final VariantContext vc, final int[] ad) {
        return new VariantContextBuilder(vc)
                .genotypes(List.of(genotype("T1", ad), genotype("N1", new int[] {90, 1, 0}))).make();
    }

    /** A fresh engine per case, so one case's clustering model cannot reach the next. */
    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "T2", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static void prob(final String label, final VariantContext vc, final double contamination) {
        final ContaminationFilter filter = new ContaminationFilter(List.of(), contamination);
        try {
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine(), null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
