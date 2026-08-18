/*
 * The whole filtering engine over one record: every filter's answer, the per-type maxima, the
 * combined probability, and what `applyFilters` makes of them. Taken from the reference.
 *
 * Every piece is measured on its own elsewhere; this is the thing that holds them. Five behaviours
 * it is built to catch.
 *
 *   - AN EMPTY LIST, A ZERO AND A NaN ARE THREE DIFFERENT ANSWERS. `ErrorProbabilities` DROPS the
 *     filters that answered an empty list, counts the zeroes, and combines the NaN;
 *   - THE PER-TYPE MAXIMUM COMES FIRST AND THE INDEPENDENCE PRODUCT SECOND, so a filter can only be
 *     masked by another of its own error type;
 *   - THE SYMBOLIC REMOVAL HAPPENS INSIDE `ErrorProbabilities`, once per filter, and the engine must
 *     not do it again;
 *   - WHICH FILTERS ANSWERED AT ALL IS A PROPERTY OF THE RECORD'S ANNOTATIONS, not of the mode: a
 *     bare record leaves most of the list unevaluated;
 *   - AND THE COMBINED PROBABILITY IS `1 - prod(1 - p)` ROUNDED, which is not the maximum.
 *
 * Output:
 *
 *     filter\t<label>\t<class>=<probabilities>
 *     type\t<label>\t<error type>=<per-allele maxima>
 *     combined\t<label>\t<per-allele>
 *     applied\t<label>\t<FILTER column>|<AS_FilterStatus>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: EngineAssemblyDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorProbabilities;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorType;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.lang.reflect.Field;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class EngineAssemblyDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT_C = Allele.create("C", false);
    static final Allele ALT_G = Allele.create("G", false);
    static final Allele NON_REF = Allele.create("<NON_REF>", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# EngineAssemblyDump: the eighteen filters over one record");

        // A biallelic record carrying every annotation the eighteen filters read.
        run("fully-annotated", false, annotated(List.of(REF, ALT_C), 1));
        // The same record in mitochondrial mode, where six filters are not built.
        run("fully-annotated-mitochondria", true, annotated(List.of(REF, ALT_C), 1));
        // Two alternates, so the transposes and the per-allele maxima have something to do.
        run("two-alternates", false, annotated(List.of(REF, ALT_C, ALT_G), 2));
        // A record with a symbolic alternate beside a real one.
        run("with-non-ref", false, annotated(List.of(REF, ALT_C, NON_REF), 2));
        // A bare record: most of the list answers an empty list and is dropped.
        run("bare", false, new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(REF, ALT_C)).genotypes(genotypes(2)).make());
        // A record in the panel of normals, which is the one filter that reads presence.
        run("in-panel-of-normals", false,
                new VariantContextBuilder(annotated(List.of(REF, ALT_C), 1))
                        .attribute("PON", true).make());
        // Weak evidence: a TLOD that says the alternate is barely there.
        run("weak-evidence", false, new VariantContextBuilder(annotated(List.of(REF, ALT_C), 1))
                .attribute("TLOD", List.of(0.5)).make());
        // A poor alternate base quality, which the base-quality filter fires on.
        run("poor-base-quality", false, new VariantContextBuilder(annotated(List.of(REF, ALT_C), 1))
                .attribute("MBQ", List.of(30, 2)).make());
    }

    /** One tumour and one normal, hom-ref so the genotypes are valid for any allele list. */
    static List<Genotype> genotypes(final int alleleCount) {
        final int[] tumorDepths = new int[alleleCount];
        final int[] normalDepths = new int[alleleCount];
        tumorDepths[0] = 80;
        normalDepths[0] = 90;
        for (int i = 1; i < alleleCount; i++) {
            tumorDepths[i] = 30 - 10 * i;
            normalDepths[i] = 1;
        }
        final GenotypeBuilder tumor = new GenotypeBuilder("T1", List.of(REF, REF)).AD(tumorDepths);
        final double[] alleleFractions = new double[alleleCount - 1];
        for (int i = 0; i < alleleFractions.length; i++) {
            alleleFractions[i] = 0.2 / (i + 1);
        }
        tumor.attribute("AF", alleleFractions);
        tumor.attribute("PGT", "0|1");
        tumor.attribute("PID", "100_A_C");
        return List.of(tumor.make(),
                new GenotypeBuilder("N1", List.of(REF, REF)).AD(normalDepths).make());
    }

    /** Every annotation the eighteen filters read, sized for the record's alternates. */
    static VariantContext annotated(final List<Allele> alleles, final int alternateCount) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", 100, 100, alleles)
                        .genotypes(genotypes(alleles.size()));
        // Per-alternate lists.
        builder.attribute("TLOD", repeatDouble(20.0, alternateCount));
        builder.attribute("NALOD", repeatDouble(2.0, alternateCount));
        builder.attribute("POPAF", repeatDouble(6.0, alternateCount));
        // Lists that carry the reference first.
        builder.attribute("MBQ", repeatInt(30, alternateCount + 1));
        builder.attribute("MMQ", repeatInt(60, alternateCount + 1));
        builder.attribute("MFRL", repeatInt(300, alternateCount + 1));
        // `MPOS` has no reference entry.
        builder.attribute("MPOS", repeatInt(25, alternateCount));
        builder.attribute("AS_UNIQ_ALT_READ_COUNT", repeatInt(8, alternateCount));
        // The strand table has one entry per allele, reference included.
        final StringBuilder table = new StringBuilder("40,40");
        for (int i = 0; i < alternateCount; i++) {
            table.append("|10,10");
        }
        builder.attribute("AS_SB_TABLE", table.toString());
        builder.attribute("NCount", 0);
        builder.attribute("ECNT", 1);
        builder.attribute("ECNTH", 1);
        builder.attribute("RPA", repeatInt(10, alternateCount + 1));
        builder.attribute("RU", "A");
        return builder.make();
    }

    static List<Double> repeatDouble(final double value, final int count) {
        final Double[] values = new Double[count];
        java.util.Arrays.fill(values, value);
        return List.of(values);
    }

    static List<Integer> repeatInt(final int value, final int count) {
        final Integer[] values = new Integer[count];
        java.util.Arrays.fill(values, value);
        return List.of(values);
    }

    @SuppressWarnings("unchecked")
    static void run(final String label, final boolean mitochondria, final VariantContext vc)
            throws Exception {
        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        arguments.mitochondria = mitochondria;
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine =
                new Mutect2FilteringEngine(arguments, header, new File("no-such-stats-file.tsv"));
        final Field field = Mutect2FilteringEngine.class.getDeclaredField("filters");
        field.setAccessible(true);
        final List<Mutect2Filter> filters = (List<Mutect2Filter>) field.get(engine);

        try {
            final ErrorProbabilities probabilities = new ErrorProbabilities(filters, vc, engine, null);
            // Every filter that answered, in the engine's construction order.
            for (final Map.Entry<Mutect2Filter, List<Double>> entry :
                    probabilities.getProbabilitiesByFilter().entrySet()) {
                System.out.printf("filter\t%s\t%s=%s%n", label,
                        entry.getKey().getClass().getSimpleName(), entry.getValue());
            }
            for (final ErrorType type : ErrorType.values()) {
                final List<Double> perAllele = type == ErrorType.ARTIFACT
                        ? probabilities.getTechnicalArtifactProbabilities()
                        : (type == ErrorType.NON_SOMATIC ? probabilities.getNonSomaticProbabilities() : null);
                if (perAllele != null) {
                    System.out.printf("type\t%s\t%s=%s%n", label, type, perAllele);
                }
            }
            System.out.printf("combined\t%s\t%s%n", label,
                    probabilities.getCombinedErrorProbabilities());

            final VariantContext filtered = engine.applyFiltersAndAccumulateOutputStats(vc, null);
            final String column = filtered.isNotFiltered() ? "PASS"
                    : String.join(";", new TreeSet<>(filtered.getFilters()));
            System.out.printf("applied\t%s\t%s|%s%n", label, column,
                    filtered.getAttributeAsString("AS_FilterStatus", "(none)"));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
