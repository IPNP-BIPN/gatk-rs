/*
 * Coverage, MappingQualityZero and CountNs, taken from the reference through their own interface.
 *
 * These are the three annotations whose whole input is the likelihood matrix. Their arithmetic is
 * three counts; what separates them is the guard, the Java type, and what they count over, and no
 * two of the three agree on any of it.
 *
 *   - Coverage returns nothing for an EMPTY matrix, MappingQualityZero writes a zero for it, and
 *     CountNs writes a zero too, because only Coverage tests evidenceCount();
 *   - MappingQualityZero returns nothing at a NON-VARIANT site and the other two do not test the
 *     site at all;
 *   - Coverage and MappingQualityZero go through String.format("%d", ...), so their values are
 *     Strings; CountNs puts a long into an ImmutableMap, so its value is a Long. All three render
 *     identically and are three different objects to a consumer;
 *   - CountNs compares against the byte 'N', upper case only, and reaches the base through
 *     ReadUtils.getReadBaseAtReferenceCoordinate, whose bounds test uses the alignment span while
 *     its index uses the soft start: a read whose N is inside a soft clip does not count.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[<class>];...    (empty for an empty map)
 *
 * Usage: LikelihoodAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.CountNs;
import org.broadinstitute.hellbender.tools.walkers.annotator.Coverage;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.MappingQualityZero;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class LikelihoodAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();

    /** The variant sits at 105, ten bases into a read that starts at 100. */
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# LikelihoodAnnotationDump: Coverage, MappingQualityZero and CountNs");

        final InfoFieldAnnotation coverage = new Coverage();
        final InfoFieldAnnotation mq0 = new MappingQualityZero();
        final InfoFieldAnnotation countNs = new CountNs();

        final VariantContext variant = variantSite();
        final VariantContext monomorphic = referenceOnlySite();

        // No likelihoods at all: every one of the three returns an empty map.
        emit("null-likelihoods", coverage, mq0, countNs, variant, null);

        // An empty matrix: Coverage says nothing, the other two write a zero.
        emit("empty-matrix", coverage, mq0, countNs, variant, likelihoods(List.of()));

        // A sample present but carrying no evidence, which is the same evidence count by a
        // different route.
        emit("sample-without-evidence", coverage, mq0, countNs, variant,
                likelihoodsWithEmptySample());

        // Ordinary reads, none of them MAPQ 0, none of them carrying an N.
        emit("three-plain-reads", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTACGTAC", 100, "10M"),
                read("r1", 60, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"))));

        // Two of the three at MAPQ 0.
        emit("two-mapq-zero", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 0, "ACGTACGTAC", 100, "10M"),
                read("r1", 0, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"))));

        // An N exactly at the variant's start, one read out of three.
        emit("one-n-at-start", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTAN GTAC".replace(" ", ""), 100, "10M"),
                read("r1", 60, "ACGTACGTAC", 100, "10M"),
                read("r2", 60, "ACGTACGTAC", 100, "10M"))));

        // A lower-case n, which is not the byte 'N'.
        emit("lower-case-n", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTAnGTAC".replace(" ", ""), 100, "10M"))));

        // An N one base away from the variant, which is not counted.
        emit("n-beside-the-start", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTNCGTAC", 100, "10M"))));

        // An N inside a soft clip covering the variant's position: the bounds test uses the
        // alignment span, so the base cannot be reached and the read does not count.
        emit("n-in-soft-clip", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "NNNNNNACGT", 106, "6S4M"))));

        // A read that does not span the variant at all.
        emit("read-past-the-variant", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTACGTAC", 200, "10M"))));

        // A deletion over the variant's position: the operator consumes no read base.
        emit("deletion-over-start", coverage, mq0, countNs, variant, likelihoods(List.of(
                read("r0", 60, "ACGTACGTAC", 100, "5M3D5M"))));

        // Two samples, to show that all three count across the whole matrix.
        emit("two-samples", coverage, mq0, countNs, variant, twoSampleLikelihoods());

        // A non-variant site: MappingQualityZero says nothing, the other two still count.
        emit("monomorphic-site", coverage, mq0, countNs, monomorphic, likelihoods(List.of(
                read("r0", 0, "ACGTANGTAC", 100, "10M"))));
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static VariantContext variantSite() {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF, ALT1)).make();
    }

    static VariantContext referenceOnlySite() {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF)).make();
    }

    static GATKRead read(final String name, final int mappingQuality, final String bases,
                         final int start, final String cigar) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(mappingQuality);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoods(final List<GATKRead> reads) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        return new AlleleLikelihoods<>(new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(REF, ALT1), bySample);
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoodsWithEmptySample() {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>());
        bySample.put("s2", new ArrayList<>());
        return new AlleleLikelihoods<>(new IndexedSampleList(List.of("s1", "s2")),
                new IndexedAlleleList<>(REF, ALT1), bySample);
    }

    static AlleleLikelihoods<GATKRead, Allele> twoSampleLikelihoods() {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(List.of(
                read("a0", 0, "ACGTANGTAC", 100, "10M"),
                read("a1", 60, "ACGTACGTAC", 100, "10M"))));
        bySample.put("s2", new ArrayList<>(List.of(
                read("b0", 0, "ACGTACGTAC", 100, "10M"))));
        return new AlleleLikelihoods<>(new IndexedSampleList(List.of("s1", "s2")),
                new IndexedAlleleList<>(REF, ALT1), bySample);
    }

    static void emit(final String label, final InfoFieldAnnotation coverage,
                     final InfoFieldAnnotation mq0, final InfoFieldAnnotation countNs,
                     final VariantContext vc, final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        one("Coverage", label, coverage, vc, likelihoods);
        one("MappingQualityZero", label, mq0, vc, likelihoods);
        one("CountNs", label, countNs, vc, likelihoods);
    }

    static void one(final String name, final String label, final InfoFieldAnnotation annotation,
                    final VariantContext vc,
                    final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        try {
            final Map<String, Object> result = annotation.annotate(null, vc, likelihoods);
            if (result == null) {
                System.out.printf("anno\t%s\t%s\tnull%n", name, label);
                return;
            }
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                final Object value = entry.getValue();
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), value,
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\t%s\tE:%s:%s%n", name, label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }
}
