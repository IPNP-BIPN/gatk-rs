/*
 * BaseQuality, MappingQuality, ReadPosition and FragmentLength, taken from the reference through
 * their own interface.
 *
 * MBQ, MMQ, MPOS and MFRL: the median, per allele, of one number taken from each read that best
 * supports it. One shared parent, four values, and four different answers to the same question
 * about an allele no read supports:
 *
 *     BaseQuality       empty -> 0
 *     MappingQuality    empty -> 60     "we don't want a GGA mode allele with no reads to
 *     ReadPosition      empty -> 50      prejudice us against a site"
 *     FragmentLength    empty -> 0
 *
 * Three filters run before a read contributes, and they are not the same filter: isInformative
 * is the likelihood confidence against the log10 threshold, isUsableRead is a mapping quality
 * that is neither 0 nor 255, and getValueForRead may decline for its own reason. So MappingQuality
 * never sees a read of quality 0, and its median can only be 0 by way of the empty case, which is
 * 60.
 *
 * The value under the key is a Java int[], so it prints through Arrays.toString and carries the
 * class [I. And includeRefAllele() is false in the parent and overridden to true in three of the
 * four, so MPOS reports one number per ALTERNATE allele where the others report one per allele.
 *
 * Output:
 *
 *     anno\t<annotation>\t<label>\t<key>=<value>[<class>];...    (empty for an empty map)
 *
 * Usage: PerAlleleAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.BaseQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.FragmentLength;
import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.MappingQuality;
import org.broadinstitute.hellbender.tools.walkers.annotator.ReadPosition;
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

public class PerAlleleAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();

    /** The variant sits at 105, six bases into a read that starts at 100. */
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# PerAlleleAnnotationDump: MBQ, MMQ, MPOS and MFRL");

        final VariantContext variant = site(List.of(REF, ALT1));
        final VariantContext triallelic = site(List.of(REF, ALT1, ALT2));

        // No likelihoods at all.
        emit("null-likelihoods", variant, null);

        // An empty matrix: every allele takes its annotation's value for no reads, and the four
        // disagree about what that is.
        emit("empty-matrix", variant, likelihoods(List.of(), List.of(REF, ALT1)));

        // One read per allele, each clearly supporting its own, with distinct qualities so the
        // median is identifiable.
        emit("one-read-each", variant, twoReadsOneEach());

        // Three reads on the alternate, so the median is over an odd count.
        emit("three-on-alt", variant, supporting(ALT1, List.of(
                read("r0", 20, 30, 100, "10M", 300),
                read("r1", 40, 35, 100, "10M", 400),
                read("r2", 60, 40, 100, "10M", 500))));

        // Four reads, so the median interpolates and the rounding shows.
        emit("four-on-alt", variant, supporting(ALT1, List.of(
                read("r0", 20, 30, 100, "10M", 300),
                read("r1", 30, 31, 100, "10M", 301),
                read("r2", 40, 32, 100, "10M", 302),
                read("r3", 50, 33, 100, "10M", 303))));

        // Two reads whose median lands on a half, which is where FastMath.round and Math.round
        // are one apart in general and where the arithmetic one rounds up.
        emit("median-on-a-half", variant, supporting(ALT1, List.of(
                read("r0", 20, 30, 100, "10M", 300),
                read("r1", 21, 31, 100, "10M", 301))));

        // A read at mapping quality 0, which isUsableRead drops before any value is taken.
        emit("mapq-zero", variant, supporting(ALT1, List.of(
                read("r0", 0, 30, 100, "10M", 300),
                read("r1", 40, 35, 100, "10M", 400))));

        // A read at mapping quality 255, the unavailable value, dropped by the same filter.
        emit("mapq-unavailable", variant, supporting(ALT1, List.of(
                read("r0", 255, 30, 100, "10M", 300),
                read("r1", 40, 35, 100, "10M", 400))));

        // An uninformative read: the two alleles are within the log10 threshold of each other.
        emit("uninformative", variant, uninformative());

        // A read that does not span the variant: MBQ and MPOS decline, MMQ and MFRL do not.
        emit("read-past-the-variant", variant, supporting(ALT1, List.of(
                read("r0", 40, 35, 200, "10M", 400))));

        // A deletion over the variant's position: no read base there, so MBQ declines while MPOS
        // still has an index.
        emit("deletion-over-start", variant, supporting(ALT1, List.of(
                read("r0", 40, 35, 100, "5M3D5M", 400))));

        // Hard clips, which MPOS counts as if the bases were still there.
        emit("hard-clipped", variant, supporting(ALT1, List.of(
                read("r0", 40, 35, 100, "3H10M5H", 400))));

        // A read whose soft clip covers the variant: the base is unreachable for MBQ.
        emit("soft-clipped-over-start", variant, supporting(ALT1, List.of(
                read("r0", 40, 35, 106, "6S4M", 400))));

        // A negative fragment length, which MFRL takes the absolute value of.
        emit("negative-fragment-length", variant, supporting(ALT1, List.of(
                read("r0", 40, 35, 100, "10M", -400),
                read("r1", 40, 35, 100, "10M", -300))));

        // Three alleles, one of which the matrix holds and no read supports, and one the variant
        // holds and the matrix does not.
        emit("three-alleles", triallelic, likelihoodsFor(
                List.of(read("r0", 40, 35, 100, "10M", 400)),
                List.of(REF, ALT1, ALT2), new double[][] {{-10}, {-1}, {-10}}));
        emit("allele-missing-from-matrix", triallelic,
                likelihoodsFor(List.of(read("r0", 40, 35, 100, "10M", 400)),
                        List.of(REF, ALT1), new double[][] {{-10}, {-1}}));

        // Two samples, to show the traversal covers the whole matrix.
        emit("two-samples", variant, twoSamples());
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static VariantContext site(final List<Allele> alleles) {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(alleles).make();
    }

    static GATKRead read(final String name, final int mappingQuality, final int baseQuality,
                         final int start, final String cigar, final int fragmentLength) {
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
        record.setInferredInsertSize(fragmentLength);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /** The read length a cigar implies, so the bases and qualities are always the right size. */
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

    /** A matrix where every read strongly supports one named allele. */
    static AlleleLikelihoods<GATKRead, Allele> supporting(final Allele allele,
                                                          final List<GATKRead> reads) {
        final List<Allele> alleles = List.of(REF, ALT1);
        final double[][] values = new double[alleles.size()][reads.size()];
        for (int a = 0; a < alleles.size(); a++) {
            for (int e = 0; e < reads.size(); e++) {
                values[a][e] = alleles.get(a).equals(allele) ? -1 : -10;
            }
        }
        return likelihoodsFor(reads, alleles, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> twoReadsOneEach() {
        final List<GATKRead> reads = List.of(
                read("ref0", 30, 20, 100, "10M", 200),
                read("alt0", 50, 40, 100, "10M", 600));
        return likelihoodsFor(reads, List.of(REF, ALT1),
                new double[][] {{-1, -10}, {-10, -1}});
    }

    /** A read whose two likelihoods are within the log10 informative threshold of each other. */
    static AlleleLikelihoods<GATKRead, Allele> uninformative() {
        final List<GATKRead> reads = List.of(read("r0", 40, 35, 100, "10M", 400));
        return likelihoodsFor(reads, List.of(REF, ALT1), new double[][] {{-1.0}, {-1.1}});
    }

    static AlleleLikelihoods<GATKRead, Allele> twoSamples() {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(List.of(read("a0", 40, 30, 100, "10M", 300))));
        bySample.put("s2", new ArrayList<>(List.of(read("b0", 50, 40, 100, "10M", 500))));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1", "s2")),
                new IndexedAlleleList<>(REF, ALT1), bySample);
        for (int s = 0; s < 2; s++) {
            final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(s);
            matrix.set(0, 0, -10);
            matrix.set(1, 0, -1);
        }
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoods(final List<GATKRead> reads,
                                                           final List<Allele> alleles) {
        final double[][] values = new double[alleles.size()][reads.size()];
        return likelihoodsFor(reads, alleles, values);
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoodsFor(final List<GATKRead> reads,
                                                               final List<Allele> alleles,
                                                               final double[][] values) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(alleles), bySample);
        final LikelihoodMatrix<GATKRead, Allele> matrix = likelihoods.sampleMatrix(0);
        for (int a = 0; a < alleles.size(); a++) {
            for (int e = 0; e < reads.size(); e++) {
                matrix.set(a, e, values[a][e]);
            }
        }
        return likelihoods;
    }

    static void emit(final String label, final VariantContext vc,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        one("BaseQuality", label, new BaseQuality(), vc, likelihoods);
        one("MappingQuality", label, new MappingQuality(), vc, likelihoods);
        one("ReadPosition", label, new ReadPosition(), vc, likelihoods);
        one("FragmentLength", label, new FragmentLength(), vc, likelihoods);
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
                // An int[] prints as its identity hash through %s, which would make the golden
                // depend on the run. Arrays.toString is what makes it a value.
                final String shown = value instanceof int[] ? Arrays.toString((int[]) value)
                        : String.valueOf(value);
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), shown,
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s\t%s%n", name, label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\t%s\tE:%s:%s%n", name, label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }
}
