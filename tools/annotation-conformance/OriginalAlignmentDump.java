/*
 * OriginalAlignment, taken from the reference through its own interface.
 *
 * It counts the reads supporting the BEST alternate allele whose original alignment was on a
 * different contig, which is how a Mutect2 call on the mitochondrion is checked against a NuMT.
 * Four of its decisions live somewhere other than in the class:
 *
 *   - the allele it counts for is picked by TLOD through getTumorLogOdds and MathUtils
 *     .maxElementIndex, so a tie goes to the earliest alternate, and a TLOD of "." is not missing
 *     but -1, which after the log conversion is an ordinary number that can win the maximum;
 *   - the contig comparison is a string field: getOAContig splits the OA tag on "," and takes
 *     field zero. An unmapped read's tag is "*,0,*,*,0,0;", so its original contig is "*", which
 *     differs from every real contig and therefore counts;
 *   - the filter takes BestAllele.isInformative(), which compares the confidence against the
 *     log10 threshold whatever base the matrix is in;
 *   - unlike Coverage, MappingQualityZero and CountNs, this one calls Utils.nonNull(likelihoods)
 *     rather than treating a null as nothing to say.
 *
 * Output:
 *
 *     anno\t<label>\t<key>=<value>[<class>]    (empty for an empty map, E:... for a throw)
 *     oacontig\t<label>\t<getOAContig>
 *
 * Usage: OriginalAlignmentDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.AddOriginalAlignmentTags;
import org.broadinstitute.hellbender.tools.walkers.annotator.OriginalAlignment;
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

public class OriginalAlignmentDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;
    /** The likelihoods that make a read informative for the first alternate allele. */
    static final double[][] FOR_ALT1 = {{-5.0}, {0.0}, {-5.0}};
    /** ... and for the second. */
    static final double[][] FOR_ALT2 = {{-5.0}, {-5.0}, {0.0}};
    /** ... and a matrix where nothing is informative, every allele being equally likely. */
    static final double[][] UNINFORMATIVE = {{0.0}, {0.0}, {0.0}};

    public static void main(final String[] args) {
        System.out.println("# OriginalAlignmentDump: OCM, and the tag it reads");

        final OriginalAlignment annotation = new OriginalAlignment();
        final org.broadinstitute.hellbender.engine.ReferenceContext reference = referenceContext();

        // No TLOD at all: the annotation logs once and says nothing.
        emit("no-tlod", annotation, reference, site(null), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT1));

        // One TLOD, one read whose original contig differs.
        emit("one-mismatch", annotation, reference, site("10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT1));

        // The same read, whose original contig is the current one: not counted.
        emit("same-contig", annotation, reference, site("10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr1,100,+,10M,60,0;")), FOR_ALT1));

        // No OA tag: not counted, whatever else is true.
        emit("no-oa-tag", annotation, reference, site("10.0"), matrix(List.of(
                read("r0", 60)), FOR_ALT1));

        // An unmapped read's tag, whose contig field is "*".
        emit("unmapped-oa", annotation, reference, site("10.0"), matrix(List.of(
                withOa(read("r0", 60), "*,0,*,*,0,0;")), FOR_ALT1));

        // Informative for the OTHER alternate allele, so not counted for the best one.
        emit("informative-for-alt2", annotation, reference, site("10.0,1.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT2));

        // Two TLODs, the second larger: the allele counted for is the second alternate.
        emit("second-alt-wins", annotation, reference, site("1.0,10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT2));

        // A tie between the two TLODs: maxElementIndex gives it to the first.
        emit("tlod-tie", annotation, reference, site("10.0,10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT1));

        // A TLOD of ".", which is -1 and not missing.
        emit("tlod-missing", annotation, reference, site("."), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), FOR_ALT1));

        // Uninformative likelihoods: the read supports nothing well enough to be counted.
        emit("uninformative", annotation, reference, site("10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;")), UNINFORMATIVE));

        // Several reads, two of them on another contig.
        emit("three-reads", annotation, reference, site("10.0"), matrix(List.of(
                withOa(read("r0", 60), "chr2,100,+,10M,60,0;"),
                withOa(read("r1", 60), "chr1,100,+,10M,60,0;"),
                withOa(read("r2", 60), "chr3,100,+,10M,60,0;")),
                new double[][] {{-5.0, -5.0, -5.0}, {0.0, 0.0, 0.0}, {-5.0, -5.0, -5.0}}));

        // Null likelihoods, which this annotation refuses rather than guarding.
        emit("null-likelihoods", annotation, reference, site("10.0"), null);

        // getOAContig on its own, including the shapes a caller might not expect.
        oaContig("ordinary", "chr2,100,+,10M,60,0;");
        oaContig("unmapped", "*,0,*,*,0,0;");
        oaContig("contig-with-underscore", "chr_2,100,+,10M,60,0;");
        oaContig("no-comma", "chr2");
        oaContig("empty", "");
        oaContig("leading-comma", ",100,+,10M,60,0;");
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000),
                new SAMSequenceRecord("chr3", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static org.broadinstitute.hellbender.engine.ReferenceContext referenceContext() {
        return new org.broadinstitute.hellbender.engine.ReferenceContext(null, new SimpleInterval("chr1", START, START));
    }

    static VariantContext site(final String tlod) {
        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr("chr1").start(START).stop(START).alleles(List.of(REF, ALT1, ALT2));
        if (tlod != null) {
            final Map<String, Object> attributes = new LinkedHashMap<>();
            attributes.put("TLOD", tlod);
            builder.attributes(attributes);
        }
        return builder.make();
    }

    static GATKRead read(final String name, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes());
        final byte[] qualities = new byte[10];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(mappingQuality);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static GATKRead withOa(final GATKRead read, final String oa) {
        read.setAttribute(AddOriginalAlignmentTags.OA_TAG_NAME, oa);
        return read;
    }

    static AlleleLikelihoods<GATKRead, Allele> matrix(final List<GATKRead> reads,
                                                      final double[][] values) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(REF, ALT1, ALT2), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < 3; a++) {
            for (int e = 0; e < reads.size(); e++) {
                m.set(a, e, values[a][e % values[a].length]);
            }
        }
        return likelihoods;
    }

    static void emit(final String label, final OriginalAlignment annotation,
                     final org.broadinstitute.hellbender.engine.ReferenceContext reference, final VariantContext vc,
                     final AlleleLikelihoods<GATKRead, Allele> likelihoods) {
        try {
            final Map<String, Object> result = annotation.annotate(reference, vc, likelihoods);
            final StringJoiner joiner = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry : result.entrySet()) {
                final Object value = entry.getValue();
                joiner.add(String.format("%s=%s[%s]", entry.getKey(), value,
                        value == null ? "null" : value.getClass().getName()));
            }
            System.out.printf("anno\t%s\t%s%n", label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static void oaContig(final String label, final String oa) {
        try {
            final GATKRead read = withOa(read("probe", 60), oa);
            System.out.printf("oacontig\t%s\t%s%n", label,
                    AddOriginalAlignmentTags.getOAContig(read));
        } catch (final Exception | AssertionError e) {
            System.out.printf("oacontig\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }
}
