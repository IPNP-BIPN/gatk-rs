/*
 * `Mutect2FilteringEngine.accumulateData` and the pass schedule `FilterMutectCalls` drives it with,
 * taken from the reference.
 *
 * The last computation between the ported pieces and an end-to-end run; what remains after it is
 * file I/O. Four behaviours this is built to catch.
 *
 *   - THE TOOL TRAVERSES THE VCF FOUR TIMES, NOT TWO. `numberOfPasses()` is
 *     `NUMBER_OF_LEARNING_PASSES + 2`. Passes 0, 1 and 2 all accumulate; only 0 and 1 learn
 *     parameters afterwards. Pass 2 accumulates a whole traversal and uses it ONLY to relearn the
 *     threshold, the filters' parameters being deliberately frozen. Pass 3 applies and writes;
 *   - A RECORD WHOSE ONLY ALTERNATE IS `<NON_REF>` IS SKIPPED ENTIRELY, contributing to neither the
 *     clustering model nor the threshold. A record with no alternate is skipped by the same test,
 *     and a symbolic alternate that is NOT `<NON_REF>` is not skipped;
 *   - THE FILTERS ACCUMULATE BEFORE THE CLUSTERING MODEL RECORDS, and `record` mutates the tumour
 *     depth array while `getTumorLogOdds` mutates the `TLOD` array in place;
 *   - AND THE THRESHOLD CALCULATOR IS FED THE COMBINED PROBABILITIES, one per alternate allele.
 *
 * Output:
 *
 *     passes\t<name>\t<value>
 *     accumulated\t<label>\t<clustering data>,<threshold probabilities>
 *     mutated\t<label>\t<the caller's array after the call>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AccumulateDataDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class AccumulateDataDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT_C = Allele.create("C", false);
    static final Allele NON_REF = Allele.create("<NON_REF>", false);
    static final Allele DELETION_SYMBOL = Allele.create("<DEL>", false);

    public static void main(final String[] args) throws Exception {
        System.out.println("# AccumulateDataDump: what a non-final pass records");

        // The pass schedule, which is a protected method on the tool.
        final Method passes = FilterMutectCalls.class.getDeclaredMethod("numberOfPasses");
        passes.setAccessible(true);
        System.out.printf("passes\tnumberOfPasses\t%d%n", (int) passes.invoke(new FilterMutectCalls()));
        System.out.println("passes\tlearningPasses\t2");
        System.out.println("passes\taccumulatingPasses\t3");
        System.out.println("passes\tapplyingPass\t3");

        // A record with a real alternate, which is accumulated.
        accumulate("real-alternate", record(List.of(REF, ALT_C), new int[] {80, 20}, List.of(20.0)));
        // Two alternates: two combined probabilities reach the threshold calculator.
        accumulate("two-alternates", record(List.of(REF, ALT_C, Allele.create("G", false)),
                new int[] {80, 20, 5}, List.of(20.0, 6.0)));
        // Only `<NON_REF>`: skipped before anything is touched.
        accumulate("only-non-ref", record(List.of(REF, NON_REF), new int[] {80, 20}, List.of(20.0)));
        // A real alternate beside `<NON_REF>`, which is not skipped.
        accumulate("alternate-and-non-ref", record(List.of(REF, ALT_C, NON_REF),
                new int[] {80, 20, 0}, List.of(20.0, 0.0)));
        // A symbolic alternate that is not `<NON_REF>`, which is also not skipped.
        accumulate("symbolic-not-non-ref",
                record(List.of(REF, DELETION_SYMBOL), new int[] {80, 20}, List.of(20.0)));
        // No alternate at all.
        accumulate("reference-only", record(List.of(REF), new int[] {80}, List.of()));
        // No TLOD at all.
        accumulate("no-tlod", record(List.of(REF, ALT_C), new int[] {80, 20}, null));

        // Three records in a row, to show the accumulators growing rather than resetting.
        final Mutect2FilteringEngine engine = engine();
        for (int i = 0; i < 3; i++) {
            invokeAccumulate(engine, record(List.of(REF, ALT_C), new int[] {80, 20}, List.of(20.0)));
        }
        System.out.printf("accumulated\tthree-records\t%s%n", counts(engine));
    }

    /** Hom-ref calls, so the genotype's alleles are in every record whatever its alternates are. */
    static Genotype genotype(final String sample, final int[] ad) {
        return new GenotypeBuilder(sample, List.of(REF, REF)).AD(ad).make();
    }

    static VariantContext record(final List<Allele> alleles, final int[] ad,
                                 final List<Double> tumorLog10Odds) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", 100, 100, alleles)
                        .genotypes(List.of(genotype("T1", ad),
                                new GenotypeBuilder("N1", List.of(REF, REF))
                                        .AD(new int[] {90, 1}).make()));
        if (tumorLog10Odds != null) {
            builder.attribute("TLOD", tumorLog10Odds);
        }
        return builder.make();
    }

    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static void invokeAccumulate(final Mutect2FilteringEngine engine, final VariantContext vc) {
        engine.accumulateData(vc, null);
    }

    /** One record through a fresh engine, with the accumulators read before and after. */
    static void accumulate(final String label, final VariantContext vc) throws Exception {
        final Mutect2FilteringEngine engine = engine();
        try {
            invokeAccumulate(engine, vc);
            System.out.printf("accumulated\t%s\t%s%n", label, counts(engine));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /**
     * The clustering model's data, the threshold calculator's probabilities, and the count of
     * alternates the model skipped as obvious artifacts, by reflection.
     */
    static String counts(final Mutect2FilteringEngine engine) throws Exception {
        return listSize(engine, "somaticClusteringModel", "data") + ","
                + listSize(engine, "thresholdCalculator", "errorProbabilities") + ","
                + obviousArtifacts(engine);
    }

    static int obviousArtifacts(final Mutect2FilteringEngine engine) throws Exception {
        final Field ownerField =
                Mutect2FilteringEngine.class.getDeclaredField("somaticClusteringModel");
        ownerField.setAccessible(true);
        final Object model = ownerField.get(engine);
        final Field field = model.getClass().getDeclaredField("obviousArtifactCount");
        field.setAccessible(true);
        return ((org.apache.commons.lang3.mutable.MutableInt) field.get(model)).intValue();
    }

    static int listSize(final Mutect2FilteringEngine engine, final String owner, final String name)
            throws Exception {
        final Field ownerField = Mutect2FilteringEngine.class.getDeclaredField(owner);
        ownerField.setAccessible(true);
        final Object holder = ownerField.get(engine);
        final Field field = holder.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return ((List<?>) field.get(holder)).size();
    }
}
