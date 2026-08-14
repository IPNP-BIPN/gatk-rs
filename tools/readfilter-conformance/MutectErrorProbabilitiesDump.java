/*
 * ErrorProbabilities, taken from the reference.
 *
 * What every filter's answer becomes once the engine has them all: one probability per alternate
 * allele. Five behaviours this is built to catch.
 *
 *   - A FILTER THAT ANSWERS AN EMPTY LIST IS DROPPED ENTIRELY rather than counted as zero, which is
 *     how a switched-off filter differs from one that fired and found nothing;
 *   - THE FILTERS OF ONE ERROR TYPE ARE COMBINED BY THEIR MAXIMUM, not by a product: the worst
 *     filter of a type decides that type's probability for each allele;
 *   - THE TYPES ARE THEN COMBINED AS INDEPENDENT, `1 - prod(1 - p)`, and the result goes through
 *     `roundFinitePrecisionErrors`. Every probability reachable here is 0 or 1, so what this dump
 *     pins is the shape rather than the arithmetic: a fractional probability needs either the
 *     somatic clustering model or a contamination table;
 *   - TWO FILTERS THAT ANSWER DIFFERENT NUMBERS OF ALLELES ARE A REFUSAL, `transpose` validating
 *     that every list is the same size, which a record whose annotation is shorter than its allele
 *     list produces;
 *   - AND A SYMBOLIC ALTERNATE ALLELE IS REMOVED FROM EVERY LIST before any of that.
 *
 * Only ARTIFACT-type filters are used here: every NON_SOMATIC filter needs either the somatic
 * clustering model or a contamination table, so the combination ACROSS types is not exercised by
 * this dump and needs one of its own.
 *
 * Output:
 *
 *     combined\t<label>\t<one probability per alternate allele>
 *     artifact\t<label>\t<the ARTIFACT type's probability per allele>
 *     filters\t<label>\t<how many filters survived>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MutectErrorProbabilitiesDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.BaseQualityFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ClusteredEventsFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorProbabilities;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NormalArtifactFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ReadPositionFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrictStrandBiasFilter;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class MutectErrorProbabilitiesDump {

    public static void main(final String[] args) {
        System.out.println("# MutectErrorProbabilitiesDump: many filters in, one probability per allele");

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(
                new M2FiltersArgumentCollection(), header, new File("no-such-stats-file.tsv"));

        // Three hard filters of the same error type, and one that answers an empty list.
        final List<Mutect2Filter> filters = List.of(
                new BaseQualityFilter(20),
                new ReadPositionFilter(5),
                new ClusteredEventsFilter(2, 2),
                new StrictStrandBiasFilter(0));

        // One alternate allele, failing one filter of the three.
        probabilities("one-allele-one-failure", filters, engine, snp(
                entry("MBQ", List.of(30, 10)),
                entry("MPOS", List.of(20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));

        // The same allele failing two of them: the maximum is still one, not more.
        probabilities("one-allele-two-failures", filters, engine, snp(
                entry("MBQ", List.of(30, 10)),
                entry("MPOS", List.of(2)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));

        // Nothing failing at all.
        probabilities("one-allele-no-failure", filters, engine, snp(
                entry("MBQ", List.of(30, 35)),
                entry("MPOS", List.of(20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));

        // Two alternates, one of them failing: the site filter's answer is copied to both.
        probabilities("two-alleles", filters, engine, triallelic(
                entry("MBQ", List.of(30, 10, 35)),
                entry("MPOS", List.of(2, 20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));

        // A site filter that fires: every allele takes it.
        probabilities("two-alleles-site-filter", filters, engine, triallelic(
                entry("MBQ", List.of(30, 35, 35)),
                entry("MPOS", List.of(20, 20)),
                entry("ECNT", 9),
                entry("ECNTH", List.of(1))));

        // An annotation shorter than the allele list: two filters answer different lengths.
        probabilities("ragged-lists", filters, engine, triallelic(
                entry("MBQ", List.of(30, 10)),
                entry("MPOS", List.of(2, 20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));

        // A filter that consults the engine rather than an annotation alone. Its own answer would be
        // a posterior, but a normal whose allele fraction matches the tumour's trips the pileup
        // p-value short circuit first and the answer is 1.0: every probability reachable without the
        // clustering model is a bound rather than a fraction.
        final List<Mutect2Filter> withNormalArtifact = List.of(
                new BaseQualityFilter(20),
                new ReadPositionFilter(5),
                new NormalArtifactFilter(0.001));
        final VariantContext artifactish = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C))
                .attribute("MBQ", List.of(30, 35))
                .attribute("MPOS", List.of(20))
                .attribute("TLOD", List.of(6.0))
                .attribute("NALOD", List.of(-0.5))
                .genotypes(List.of(
                        new GenotypeBuilder("T1", List.of(Allele.REF_A, Allele.ALT_C))
                                .AD(new int[] {80, 20}).make(),
                        new GenotypeBuilder("N1", List.of(Allele.REF_A, Allele.ALT_C))
                                .AD(new int[] {80, 20}).make()))
                .make();
        probabilities("normal-artifact-pvalue", withNormalArtifact, engine, artifactish);

        // The same record with a filter that fires outright beside it: the maximum wins.
        probabilities("normal-artifact-pvalue-with-hard-failure", withNormalArtifact, engine,
                new VariantContextBuilder(artifactish).attribute("MPOS", List.of(2)).make());

        // A symbolic alternate, which is removed from every list before the combination.
        probabilities("symbolic-allele", filters, engine, symbolic(
                entry("MBQ", List.of(30, 10, 35)),
                entry("MPOS", List.of(2, 20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1))));
    }

    static void probabilities(final String label, final List<Mutect2Filter> filters,
                              final Mutect2FilteringEngine engine, final VariantContext vc) {
        try {
            final ErrorProbabilities probabilities =
                    new ErrorProbabilities(filters, vc, engine, null);
            System.out.printf("combined\t%s\t%s%n", label,
                    probabilities.getCombinedErrorProbabilities());
            System.out.printf("artifact\t%s\t%s%n", label,
                    probabilities.getTechnicalArtifactProbabilities());
            System.out.printf("filters\t%s\t%d%n", label,
                    probabilities.getProbabilitiesByFilter().size());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    @SafeVarargs
    static VariantContext snp(final Map.Entry<String, Object>... attributes) {
        return build(List.of(Allele.REF_A, Allele.ALT_C), attributes);
    }

    @SafeVarargs
    static VariantContext triallelic(final Map.Entry<String, Object>... attributes) {
        return build(List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), attributes);
    }

    @SafeVarargs
    static VariantContext symbolic(final Map.Entry<String, Object>... attributes) {
        return build(List.of(Allele.REF_A, Allele.ALT_C, Allele.create("<NON_REF>", false)),
                attributes);
    }

    @SafeVarargs
    static VariantContext build(final List<Allele> alleles,
                                final Map.Entry<String, Object>... attributes) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", 100, 100, alleles);
        for (final Map.Entry<String, Object> attribute : attributes) {
            builder.attribute(attribute.getKey(), attribute.getValue());
        }
        return builder.make();
    }

    static Map.Entry<String, Object> entry(final String key, final Object value) {
        return Map.entry(key, value);
    }
}
