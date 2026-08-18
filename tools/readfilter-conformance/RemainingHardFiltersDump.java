/*
 * The four remaining hard filters, taken from the reference.
 *
 * `DuplicatedAltReadFilter`, `NRatioFilter`, `MinAlleleFractionFilter` and `PanelOfNormalsFilter`:
 * the rest of the `HardFilter`/`HardAlleleFilter` family, none of which does any arithmetic. Six
 * behaviours this is built to catch.
 *
 *   - ALL FOUR ARE BUILT UNCONDITIONALLY AND NONE OF THEM CAN FIRE WITH THE DEFAULT ARGUMENTS.
 *     `DEFAULT_MIN_UNIQUE_ALT_READS = 0` against `count <= uniqueAltReadCount`,
 *     `DEFAULT_MAX_N_RATIO = Double.POSITIVE_INFINITY` against `ratio >= maxNRatio`, and
 *     `DEFAULT_MIN_AF = 0` against `max < minAf`. `PanelOfNormalsFilter` needs an annotation that
 *     only exists when Mutect2 was given a panel;
 *   - `NRatioFilter`'S COMMENT CLAIMS A GUARD THE CODE DOES NOT HAVE. "if there is no NCount
 *     annotation or the altCount is 0, don't apply the filter", but only the `altCount == 0` arm is
 *     written: a missing `NCount` is `getAttributeAsInt(key, 0)`, a zero rather than a skip;
 *   - `MinAlleleFractionFilter` REQUIRES NO ANNOTATION AT ALL, so it is evaluated on every record,
 *     and an allele with no data is `orElse(1.0)`: it answers "not an artifact" from an absence;
 *   - AND ITS `.filter(entry -> !vc.getReference().equals(entry.getKey()))` IS DEAD CODE, because
 *     `getAltDataByAllele` keys its map on the alternate alleles alone;
 *   - `PanelOfNormalsFilter` READS PRESENCE, NOT VALUE: `hasAttribute` is true for `PON=false`;
 *   - AND `DuplicatedAltReadFilter`'S LIST LENGTH COMES FROM THE ANNOTATION, NOT THE RECORD, so a
 *     short list answers for fewer alleles than the record has and an absent one is the empty list.
 *
 * Output:
 *
 *     name\t<filter>\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     default\t<argument>\t<value>
 *     filter\t<label>\t<one probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: RemainingHardFiltersDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.DuplicatedAltReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MinAlleleFractionFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NRatioFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.PanelOfNormalsFilter;

import java.io.File;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class RemainingHardFiltersDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# RemainingHardFiltersDump: the rest of the hard filter family");

        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        System.out.printf("default\tuniqueAltReadCount\t%d%n", arguments.uniqueAltReadCount);
        System.out.printf("default\tnRatio\t%s%n", Double.toString(arguments.nRatio));
        System.out.printf("default\tminAf\t%s%n", Double.toString(arguments.minAf));

        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(arguments, header,
                new File("no-such-stats-file.tsv"));

        final Mutect2Filter duplicated = new DuplicatedAltReadFilter(arguments.uniqueAltReadCount);
        final Mutect2Filter duplicatedTuned = new DuplicatedAltReadFilter(2);
        final Mutect2Filter nRatio = new NRatioFilter(arguments.nRatio);
        final Mutect2Filter nRatioTuned = new NRatioFilter(0.5);
        final Mutect2Filter minAf = new MinAlleleFractionFilter(arguments.minAf);
        final Mutect2Filter minAfTuned = new MinAlleleFractionFilter(0.1);
        final Mutect2Filter pon = new PanelOfNormalsFilter();
        for (final Mutect2Filter filter : List.of(duplicated, nRatio, minAf, pon)) {
            System.out.printf("name\t%s\t%s,%s,%s%n", filter.getClass().getSimpleName(),
                    filter.filterName(), filter.errorType(),
                    filter.phredScaledPosteriorAnnotationName().orElse("none"));
        }

        // A triallelic record carrying every annotation the four read, with a tumour and a normal.
        final VariantContext base = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G))
                .attribute("AS_UNIQ_ALT_READ_COUNT", List.of(3, 1))
                .attribute("NCount", 4)
                .genotypes(List.of(
                        genotype("T1", new int[] {80, 20, 5}, List.of(0.2, 0.05)),
                        genotype("N1", new int[] {90, 1, 0}, List.of(0.01, 0.0))))
                .make();

        // DuplicatedAltReadFilter: the default threshold of 0 against counts of 3 and 1.
        filter("duplicate-default", duplicated, engine, base);
        filter("duplicate-threshold-two", duplicatedTuned, engine, base);
        // A list shorter than the alternate alleles, and one longer.
        filter("duplicate-short-list", duplicatedTuned, engine,
                new VariantContextBuilder(base).attribute("AS_UNIQ_ALT_READ_COUNT", List.of(1)).make());
        filter("duplicate-long-list", duplicatedTuned, engine,
                new VariantContextBuilder(base)
                        .attribute("AS_UNIQ_ALT_READ_COUNT", List.of(1, 1, 1, 1)).make());
        // No annotation: the required-annotation check answers an empty list.
        filter("duplicate-no-annotation", duplicatedTuned, engine,
                new VariantContextBuilder(base).rmAttribute("AS_UNIQ_ALT_READ_COUNT").make());

        // NRatioFilter: 4 Ns against 26 alternate reads over both samples.
        filter("nratio-default", nRatio, engine, base);
        filter("nratio-threshold-half", nRatioTuned, engine, base);
        // The alternate count is zero, which is the one guard the filter has.
        filter("nratio-no-alt-reads", nRatioTuned, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 0, 0}, List.of(0.0, 0.0)),
                        genotype("N1", new int[] {90, 0, 0}, List.of(0.0, 0.0)))).make());
        // No NCount: the required-annotation check answers 0.0 before the default can be used.
        filter("nratio-no-annotation", nRatioTuned, engine,
                new VariantContextBuilder(base).rmAttribute("NCount").make());
        // A ratio exactly at the threshold, which `>=` calls an artifact.
        filter("nratio-at-the-threshold", nRatioTuned, engine,
                new VariantContextBuilder(base).attribute("NCount", 13).make());

        // MinAlleleFractionFilter: it requires no annotation, so every record reaches it.
        filter("minaf-default", minAf, engine, base);
        filter("minaf-threshold-tenth", minAfTuned, engine, base);
        // No AF on any genotype at all: every allele is `orElse(1.0)`.
        filter("minaf-no-allele-fraction", minAfTuned, engine, new VariantContextBuilder(base)
                .genotypes(List.of(
                        new GenotypeBuilder("T1", List.of(Allele.REF_A, Allele.ALT_C))
                                .AD(new int[] {80, 20, 5}).make())).make());
        // The normal carries a low fraction and the tumour a high one: only the tumour is read.
        filter("minaf-normal-is-low", minAfTuned, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}, List.of(0.5, 0.5)),
                        genotype("N1", new int[] {90, 1, 0}, List.of(0.001, 0.001)))).make());
        // An AF list as long as the RECORD rather than as long as the alternates: the zip stops at
        // the map's two entries and the last value is dropped.
        filter("minaf-full-length-list", minAfTuned, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}, List.of(0.9, 0.05, 0.9))))
                .make());
        // A fraction exactly at the threshold, which the strict `<` keeps.
        filter("minaf-at-the-threshold", minAfTuned, engine, new VariantContextBuilder(base)
                .genotypes(List.of(genotype("T1", new int[] {80, 20, 5}, List.of(0.1, 0.1)))).make());

        // PanelOfNormalsFilter: presence, not value.
        filter("pon-absent", pon, engine, base);
        filter("pon-present", pon, engine,
                new VariantContextBuilder(base).attribute("PON", true).make());
        filter("pon-false", pon, engine,
                new VariantContextBuilder(base).attribute("PON", false).make());
        filter("pon-empty-string", pon, engine,
                new VariantContextBuilder(base).attribute("PON", "").make());
    }

    static Genotype genotype(final String sample, final int[] ad, final List<Double> alleleFractions) {
        return new GenotypeBuilder(sample, List.of(Allele.REF_A, Allele.ALT_C))
                .AD(ad).attribute("AF", alleleFractions).make();
    }

    static void filter(final String label, final Mutect2Filter filter,
                       final Mutect2FilteringEngine engine, final VariantContext vc) {
        try {
            System.out.printf("filter\t%s\t%s%n", label, filter.errorProbabilities(vc, engine, null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
