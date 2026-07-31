/*
 * BaseQRankSum, MQRankSum, ReadPosRankSum and ClippingRankSum, taken from the reference through
 * their own interface.
 *
 * A Mann-Whitney U test of the ALTERNATE reads against the REFERENCE reads, reported as a Z score.
 * Four things decide what a consumer sees:
 *
 *   - the alternate goes FIRST. Swapping the two arrays flips the sign of every value this family
 *     reports, and nothing downstream would notice;
 *   - the value is a STRING, formatted with String.format("%.3f", z), so it rounds half-up on the
 *     decimal expansion rather than however a Double would render;
 *   - a NaN Z score is not written at all: the key is absent from the record rather than present
 *     with a placeholder, so a site with no alternate reads has no MQRankSum field;
 *   - three filters run, and ReadPosRankSum overrides one of them to add a soft-clip test against
 *     vc.getEnd() + 1, so the four members do not see the same reads.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[<class>];...    (empty for an empty map)
 *
 * Usage: RankSumAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.ClippingRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.MappingQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosRankSumTest;
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

public class RankSumAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# RankSumAnnotationDump: BaseQRankSum, MQRankSum, ReadPosRankSum, ClippingRankSum");

        final VariantContext vc = site(true);
        final VariantContext noGenotypes = site(false);

        emit("null-likelihoods", vc, null);
        emit("no-genotypes", noGenotypes, twoSided(6, 6, 30, 20, 60, 40));
        emit("empty-matrix", vc, likelihoods(List.of(), List.of()));

        // Below the exact/normal boundary on both sides, so the permutation test runs.
        emit("three-and-three", vc, twoSided(3, 3, 30, 20, 60, 40));
        emit("nine-and-nine", vc, twoSided(9, 9, 30, 20, 60, 40));
        // At the boundary, where the normal approximation takes over.
        emit("ten-and-nine", vc, twoSided(10, 9, 30, 20, 60, 40));
        emit("twelve-and-twelve", vc, twoSided(12, 12, 30, 20, 60, 40));

        // No difference at all between the two groups, which is where the tie handling shows.
        emit("identical-groups", vc, twoSided(12, 12, 30, 30, 60, 60));
        // Only reference reads, so one series is empty and the Z is NaN.
        emit("ref-only", vc, twoSided(12, 0, 30, 20, 60, 40));
        emit("alt-only", vc, twoSided(0, 12, 30, 20, 60, 40));

        // Mapping quality 0 and 255, which isUsableRead drops before any value is taken.
        emit("mapq-zero-and-unavailable", vc, mixedMappingQualities());

        // Reads that do not span the variant, where the four members disagree about what to do.
        emit("reads-past-the-variant", vc, pastTheVariant());
        // Hard clips, which only ClippingRankSum counts and which ReadPosRankSum measures around.
        emit("hard-clipped", vc, hardClipped());
        // A read opening with an insertion at exactly vc.getEnd() + 1.
        emit("leading-insertion", vc, leadingInsertion());
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static VariantContext site(final boolean withGenotypes) {
        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr("chr1").start(VARIANT_START).stop(VARIANT_START).alleles(List.of(REF, ALT));
        if (withGenotypes) {
            final Genotype genotype = new GenotypeBuilder("s1", List.of(REF, ALT)).make();
            builder.genotypes(List.of(genotype));
        }
        return builder.make();
    }

    static GATKRead read(final String name, final int mappingQuality, final int baseQuality,
                         final int start, final String cigar) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        final int length = readLength(cigar);
        final byte[] bases = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setMappingQuality(mappingQuality);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static int readLength(final String cigar) {
        int total = 0;
        int number = 0;
        for (final char c : cigar.toCharArray()) {
            if (Character.isDigit(c)) {
                number = number * 10 + (c - '0');
            } else {
                if (c == 'M' || c == 'I' || c == 'S') {
                    total += number;
                }
                number = 0;
            }
        }
        return total;
    }

    /** `refCount` reads supporting the reference and `altCount` supporting the alternate. */
    static AlleleLikelihoods<GATKRead, Allele> twoSided(final int refCount, final int altCount,
                                                        final int refBaseQuality,
                                                        final int altBaseQuality,
                                                        final int refMappingQuality,
                                                        final int altMappingQuality) {
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < refCount; i++) {
            reads.add(read("r" + i, refMappingQuality - i, refBaseQuality - i, 100, "10M"));
        }
        for (int i = 0; i < altCount; i++) {
            reads.add(read("a" + i, altMappingQuality - i, altBaseQuality - i, 100, "10M"));
        }
        final double[][] values = new double[2][reads.size()];
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < refCount;
            values[0][e] = isRef ? -1 : -10;
            values[1][e] = isRef ? -10 : -1;
        }
        return likelihoodsFor(reads, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> mixedMappingQualities() {
        final List<GATKRead> reads = List.of(
                read("r0", 0, 30, 100, "10M"),
                read("r1", 255, 30, 100, "10M"),
                read("r2", 60, 30, 100, "10M"),
                read("a0", 0, 20, 100, "10M"),
                read("a1", 255, 20, 100, "10M"),
                read("a2", 40, 20, 100, "10M"));
        final double[][] values = new double[2][reads.size()];
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 3;
            values[0][e] = isRef ? -1 : -10;
            values[1][e] = isRef ? -10 : -1;
        }
        return likelihoodsFor(reads, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> pastTheVariant() {
        final List<GATKRead> reads = List.of(
                read("r0", 60, 30, 100, "10M"),
                read("r1", 60, 30, 200, "10M"),
                read("a0", 40, 20, 100, "10M"),
                read("a1", 40, 20, 200, "10M"));
        final double[][] values = new double[2][reads.size()];
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 2;
            values[0][e] = isRef ? -1 : -10;
            values[1][e] = isRef ? -10 : -1;
        }
        return likelihoodsFor(reads, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> hardClipped() {
        final List<GATKRead> reads = List.of(
                read("r0", 60, 30, 100, "10M"),
                read("r1", 60, 30, 100, "3H10M"),
                read("a0", 40, 20, 100, "5H10M5H"),
                read("a1", 40, 20, 100, "10M2H"));
        final double[][] values = new double[2][reads.size()];
        for (int e = 0; e < reads.size(); e++) {
            final boolean isRef = e < 2;
            values[0][e] = isRef ? -1 : -10;
            values[1][e] = isRef ? -10 : -1;
        }
        return likelihoodsFor(reads, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> leadingInsertion() {
        final List<GATKRead> reads = List.of(
                read("r0", 60, 30, 100, "10M"),
                read("a0", 40, 20, VARIANT_START + 1, "3I7M"));
        final double[][] values = {{-1, -10}, {-10, -1}};
        return likelihoodsFor(reads, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoods(final List<GATKRead> reads,
                                                            final List<Allele> unusedAlleles) {
        return likelihoodsFor(reads, new double[2][reads.size()]);
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoodsFor(final List<GATKRead> reads,
                                                               final double[][] values) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
        for (int a = 0; a < 2; a++) {
            for (int e = 0; e < reads.size(); e++) {
                matrix.set(a, e, values[a][e]);
            }
        }
        return likelihoods;
    }

    static void emit(final String label, final VariantContext vc,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        one("BaseQualityRankSumTest", label, new BaseQualityRankSumTest(), vc, likelihoods);
        one("MappingQualityRankSumTest", label, new MappingQualityRankSumTest(), vc, likelihoods);
        one("ReadPosRankSumTest", label, new ReadPosRankSumTest(), vc, likelihoods);
        one("ClippingRankSumTest", label, new ClippingRankSumTest(), vc, likelihoods);
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
            System.out.printf("anno\t%s\t%s\tE:%s:%s%n", name, label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }
}
