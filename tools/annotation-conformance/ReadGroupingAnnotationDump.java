/*
 * UniqueAltReadCount, BaseQualityHistogram and ReferenceBases, taken from the reference.
 *
 * Three annotations whose output SHAPE is more interesting than their arithmetic:
 *
 *   - AS_UNIQ_ALT_READ_COUNT counts distinct (start, fragmentLength) pairs, not reads, so a
 *     hundred PCR duplicates of one fragment count once and two genuinely different fragments
 *     that share both numbers also count once. The value is a String joined with "|", because it
 *     is an allele-specific raw annotation;
 *   - BQHIST is a flat list: one entry per distinct quality, then one count per allele OF THE
 *     MATRIX in matrix order, so its length depends on the matrix rather than on the variant;
 *   - REF_BASES is twenty-one bases taken out of whatever window the reference context carries,
 *     with the discard clamped at zero and the right end padded with N, so the string has a fixed
 *     length and is not necessarily centred on the variant.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[<class>]
 *     refbases\t<label>\t<window start>\t<bases>\t<variant start>\t<result>
 *
 * Usage: ReadGroupingAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQualityHistogram;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReferenceBases;
// The engine has a class of the same name for the bases themselves, so it is named in full below.
import org.broadinstitute.hellbender.tools.walkers.annotator.UniqueAltReadCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class ReadGroupingAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# ReadGroupingAnnotationDump: AS_UNIQ_ALT_READ_COUNT, BQHIST, REF_BASES");

        emit("empty-matrix", site(false), matrix(List.of(), false));
        emit("three-distinct-fragments", site(false), fragments(new int[][] {
                {100, 300}, {101, 300}, {102, 300}}, 3));
        emit("three-duplicates", site(false), fragments(new int[][] {
                {100, 300}, {100, 300}, {100, 300}}, 3));
        emit("same-start-different-length", site(false), fragments(new int[][] {
                {100, 300}, {100, 301}}, 2));
        emit("mixed", site(false), fragments(new int[][] {
                {100, 300}, {100, 300}, {101, 300}, {101, 300}, {102, 400}}, 5));
        emit("ref-and-alt", site(false), mixedAlleles());
        emit("two-alternates", site(true), threeAlleleMatrix());
        emit("varied-qualities", site(false), variedQualities());
        emit("mapq-zero-dropped", site(false), mapqZero());

        // REF_BASES on its own, over windows that are centred, off-centre and short.
        final String bases = "ACGTACGTACGTACGTACGTACGTACGTAC";
        refBases("centred", 95, bases, VARIANT_START);
        refBases("window-starts-at-variant", VARIANT_START, bases, VARIANT_START);
        refBases("window-starts-late", 104, bases, VARIANT_START);
        refBases("short-window", 100, "ACGTACGTACGT", VARIANT_START);
        refBases("one-base-window", VARIANT_START, "A", VARIANT_START);
        refBases("exactly-twenty-one", 95, "ACGTACGTACGTACGTACGTA", VARIANT_START);
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static VariantContext site(final boolean triallelic) {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(triallelic ? List.of(REF, ALT, ALT2) : List.of(REF, ALT)).make();
    }

    static GATKRead read(final String name, final int start, final int fragmentLength,
                         final int baseQuality, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setMappingQuality(mappingQuality);
        record.setInferredInsertSize(fragmentLength);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /** Every read supports the alternate, with the given (start, fragmentLength) pairs. */
    static AlleleLikelihoods<GATKRead, Allele> fragments(final int[][] pairs, final int unused) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < pairs.length; i++) {
            reads.add(read("r" + i, pairs[i][0], pairs[i][1], 30, 60));
        }
        return matrix(reads, false);
    }

    static AlleleLikelihoods<GATKRead, Allele> mixedAlleles() {
        final List<GATKRead> reads = List.of(
                read("r0", 100, 300, 30, 60),
                read("r1", 100, 300, 31, 60),
                read("a0", 101, 300, 32, 60),
                read("a1", 102, 300, 33, 60));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 2;
            m.set(0, e, isRef ? -1 : -10);
            m.set(1, e, isRef ? -10 : -1);
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> threeAlleleMatrix() {
        final List<GATKRead> reads = List.of(
                read("r0", 100, 300, 30, 60),
                read("a0", 101, 300, 31, 60),
                read("b0", 102, 300, 32, 60));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(List.of(REF, ALT, ALT2)), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < 3; a++) {
            for (int e = 0; e < reads.size(); e++) {
                m.set(a, e, a == e ? -1 : -10);
            }
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> variedQualities() {
        final List<GATKRead> reads = List.of(
                read("r0", 100, 300, 20, 60),
                read("r1", 100, 301, 30, 60),
                read("a0", 101, 300, 20, 60),
                read("a1", 101, 301, 40, 60));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 2;
            m.set(0, e, isRef ? -1 : -10);
            m.set(1, e, isRef ? -10 : -1);
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> mapqZero() {
        final List<GATKRead> reads = List.of(
                read("r0", 100, 300, 30, 0),
                read("a0", 101, 300, 31, 0),
                read("a1", 102, 300, 32, 60));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 1;
            m.set(0, e, isRef ? -1 : -10);
            m.set(1, e, isRef ? -10 : -1);
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> matrix(final List<GATKRead> reads,
                                                       final boolean unused) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            m.set(0, e, -10);
            m.set(1, e, -1);
        }
        return likelihoods;
    }

    static void refBases(final String label, final int windowStart, final String bases,
                         final int variantStart) {
        final SimpleInterval window = new SimpleInterval("chr1", windowStart,
                windowStart + bases.length() - 1);
        final org.broadinstitute.hellbender.utils.reference.ReferenceBases sequence =
                new org.broadinstitute.hellbender.utils.reference.ReferenceBases(
                        bases.getBytes(), window);
        final ReferenceContext context = new ReferenceContext(
                new ReferenceMemorySource(sequence, HEADER.getSequenceDictionary()), window);
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(variantStart)
                .stop(variantStart).alleles(List.of(REF, ALT)).make();
        try {
            System.out.printf("refbases\t%s\t%d\t%s\t%d\t%s%n", label, windowStart, bases,
                    variantStart, ReferenceBases.annotate(context, vc));
        } catch (final Exception | AssertionError e) {
            System.out.printf("refbases\t%s\t%d\t%s\t%d\tE:%s%n", label, windowStart, bases,
                    variantStart, e.getClass().getName());
        }
    }

    static void emit(final String label, final VariantContext vc,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        one("UniqueAltReadCount", label, new UniqueAltReadCount(), vc, likelihoods);
        one("BaseQualityHistogram", label, new BaseQualityHistogram(), vc, likelihoods);
    }

    static void one(final String name, final String label, final InfoFieldAnnotation annotation,
                    final VariantContext vc,
                    final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        try {
            final Map<String, Object> result = annotation.annotate(null, vc, likelihoods);
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                final Object value = entry.getValue();
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), value,
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
        }
    }
}
