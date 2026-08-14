/*
 * Mutect's engine-free hard filters, taken from the reference.
 *
 * Six of `FilterMutectCalls`' filters read one INFO annotation and answer yes or no without
 * consulting the filtering engine at all, which is what makes them measurable on their own. Seven
 * behaviours this is built to catch.
 *
 *   - ONLY A LONG INSERTION IS JUDGED BY THE REFERENCE'S MAPPING QUALITY. `MappingQualityFilter`
 *     substitutes the ref annotation for any allele whose indel length reaches the long-indel size,
 *     because an indel that maps uniquely still maps badly, and `getIndelLengths()` is the ALT
 *     length minus the REF length: negative for a deletion. A seven-base deletion is therefore
 *     judged on its own poor mapping quality while a seven-base insertion is rescued;
 *   - AND THE `remove(0)` IT DOES THAT WITH TOUCHES NOTHING, `getAttributeAsIntList` having built a
 *     list of its own: the record's annotation is the same before and after;
 *   - A NEGATIVE MEDIAN READ POSITION IS NEVER AN ARTIFACT. `ReadPositionFilter` guards with
 *     `readPos > -1` and the reference's own comment calls it a bug, linking GATK issue 5492;
 *   - MPOS HAS NO REFERENCE ENTRY while MBQ and MMQ do, so one filter reads the whole list and the
 *     other two skip the first element;
 *   - THE FRAGMENT LENGTH FILTER LOOKS AT ONE ALLELE ONLY, `get(1) - get(0)`, and being a site-level
 *     filter its single answer is then copied to every alternate allele;
 *   - STRICT STRAND BIAS ANSWERS AN EMPTY LIST rather than a list of falses when it is switched off
 *     or the strand table is missing, which is not the same length as the allele list;
 *   - AND THE MULTIALLELIC FILTER COUNTS ALLELES OVER A HARD-CODED LOD OF 5.0, separate from the
 *     allele-count threshold it was given.
 *
 * Output:
 *
 *     filter\t<filter>-<record>\t<one boolean per alternate allele, comma separated>
 *     error\t<filter>-<record>\t<exception class>:<message>
 *
 * Usage: MutectHardFiltersDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.BaseQualityFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ClusteredEventsFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FragmentLengthFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MappingQualityFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MultiallelicFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ReadPositionFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrictStrandBiasFilter;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class MutectHardFiltersDump {

    public static void main(final String[] args) {
        System.out.println("# MutectHardFiltersDump: one annotation in, one boolean per allele out");

        final Map<String, VariantContext> records = new LinkedHashMap<>();

        // A biallelic SNP whose alternate is poor on every annotation.
        records.put("poor-snp", snp(
                entry("MBQ", List.of(30, 10)),
                entry("MMQ", List.of(60, 20)),
                entry("MFRL", List.of(300, 380)),
                entry("MPOS", List.of(2)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1)),
                entry("TLOD", List.of(6.0)),
                entry("AS_SB_TABLE", "5,5|0,7")));

        // The same record with everything comfortably good.
        records.put("good-snp", snp(
                entry("MBQ", List.of(30, 35)),
                entry("MMQ", List.of(60, 60)),
                entry("MFRL", List.of(300, 305)),
                entry("MPOS", List.of(20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1)),
                entry("TLOD", List.of(6.0)),
                entry("AS_SB_TABLE", "5,5|4,3")));

        // Three alternates, two of them over the multiallelic LOD and one under it.
        records.put("triallelic", triallelic(
                entry("MBQ", List.of(30, 10, 35)),
                entry("MMQ", List.of(60, 20, 60)),
                entry("MFRL", List.of(300, 380, 302)),
                entry("MPOS", List.of(2, 20)),
                entry("ECNT", 3),
                entry("ECNTH", List.of(1, 4)),
                entry("TLOD", List.of(6.0, 4.0)),
                entry("AS_SB_TABLE", "5,5|0,7|3,4")));

        // A long deletion, whose own mapping quality is poor and whose reference's is not.
        records.put("long-deletion", deletion(
                entry("MBQ", List.of(30, 35)),
                entry("MMQ", List.of(60, 20)),
                entry("MFRL", List.of(300, 305)),
                entry("MPOS", List.of(20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1)),
                entry("TLOD", List.of(6.0)),
                entry("AS_SB_TABLE", "5,5|4,3")));

        // A long INSERTION, whose indel length is positive and therefore reaches the long-indel size.
        records.put("long-insertion", insertion(
                entry("MBQ", List.of(30, 35)),
                entry("MMQ", List.of(60, 20)),
                entry("MFRL", List.of(300, 305)),
                entry("MPOS", List.of(20)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1)),
                entry("TLOD", List.of(6.0)),
                entry("AS_SB_TABLE", "5,5|4,3")));

        // A median read position of -1, which the reference's own comment calls a bug.
        records.put("negative-position", snp(
                entry("MBQ", List.of(30, 35)),
                entry("MMQ", List.of(60, 60)),
                entry("MFRL", List.of(300, 305)),
                entry("MPOS", List.of(-1)),
                entry("ECNT", 1),
                entry("ECNTH", List.of(1)),
                entry("TLOD", List.of(6.0)),
                entry("AS_SB_TABLE", "5,5|4,3")));

        // Nothing at all: every filter meets a missing key.
        records.put("no-annotations", snp());

        for (final Map.Entry<String, VariantContext> record : records.entrySet()) {
            final String name = record.getKey();
            final VariantContext vc = record.getValue();
            alleleFilter("base-quality", name, new BaseQualityFilter(20), vc);
            alleleFilter("mapping-quality", name, new MappingQualityFilter(30, 5), vc);
            alleleFilter("read-position", name, new ReadPositionFilter(5), vc);
            alleleFilter("strict-strand", name, new StrictStrandBiasFilter(1), vc);
            siteFilter("fragment-length", name, new FragmentLengthFilter(50), vc);
            siteFilter("clustered-events", name, new ClusteredEventsFilter(2, 2), vc);
            siteFilter("multiallelic", name, new MultiallelicFilter(1), vc);
        }

        // Switched off, the strand filter answers an empty list whatever the record holds.
        alleleFilter("strict-strand-off", "poor-snp", new StrictStrandBiasFilter(0),
                records.get("poor-snp"));

        // The list the mapping-quality filter was handed, before and after.
        final List<Integer> handed = new ArrayList<>(List.of(60, 20));
        final VariantContext mutated = new VariantContextBuilder("dump", "chr1", 100, 100,
                List.of(Allele.REF_A, Allele.ALT_C))
                .attribute("MMQ", handed).make();
        System.out.printf("filter\tmutation-before\t%s%n", handed);
        new MappingQualityFilter(30, 5).areAllelesArtifacts(mutated, null, null);
        System.out.printf("filter\tmutation-after\t%s%n", handed);
    }

    static void alleleFilter(final String filter, final String record,
                             final org.broadinstitute.hellbender.tools.walkers.mutect.filtering.HardAlleleFilter instance,
                             final VariantContext vc) {
        try {
            System.out.printf("filter\t%s-%s\t%s%n", filter, record,
                    instance.areAllelesArtifacts(vc, null, null));
        } catch (final Exception e) {
            System.out.printf("error\t%s-%s\t%s:%s%n", filter, record, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void siteFilter(final String filter, final String record,
                           final org.broadinstitute.hellbender.tools.walkers.mutect.filtering.HardFilter instance,
                           final VariantContext vc) {
        try {
            System.out.printf("filter\t%s-%s\t%s%n", filter, record,
                    instance.isArtifact(vc, null));
        } catch (final Exception e) {
            System.out.printf("error\t%s-%s\t%s:%s%n", filter, record, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    @SafeVarargs
    static VariantContext snp(final Map.Entry<String, Object>... attributes) {
        return build(100, 100, List.of(Allele.REF_A, Allele.ALT_C), attributes);
    }

    @SafeVarargs
    static VariantContext triallelic(final Map.Entry<String, Object>... attributes) {
        return build(100, 100, List.of(Allele.REF_A, Allele.ALT_C, Allele.ALT_G), attributes);
    }

    /** A deletion of seven bases, whose indel length is NEGATIVE seven. */
    @SafeVarargs
    static VariantContext deletion(final Map.Entry<String, Object>... attributes) {
        return build(100, 107, List.of(Allele.create("ATTTTTTT", true), Allele.create("A", false)),
                attributes);
    }

    /** An insertion of seven bases, whose indel length is positive seven. */
    @SafeVarargs
    static VariantContext insertion(final Map.Entry<String, Object>... attributes) {
        return build(100, 100, List.of(Allele.create("A", true), Allele.create("ATTTTTTT", false)),
                attributes);
    }

    @SafeVarargs
    static VariantContext build(final int start, final int end, final List<Allele> alleles,
                                final Map.Entry<String, Object>... attributes) {
        final VariantContextBuilder builder =
                new VariantContextBuilder("dump", "chr1", start, end, alleles);
        for (final Map.Entry<String, Object> attribute : attributes) {
            builder.attribute(attribute.getKey(), attribute.getValue());
        }
        return builder.make();
    }

    static Map.Entry<String, Object> entry(final String key, final Object value) {
        return Map.entry(key, value);
    }
}
