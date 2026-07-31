/*
 * OrientationBiasReadCounts and FragmentDepthPerAlleleBySample, and the fragment grouping under
 * them, taken from the reference.
 *
 *   - groupEvidence SUMS the log likelihoods of a group rather than averaging them, so a pair whose
 *     reads each support an allele at -1 gives the fragment a -2, and the informativeness threshold
 *     is easier to clear after grouping than before;
 *   - the new evidence order is a HashMap's over the read names, because groupingBy builds one;
 *   - the pair's orientation, usability and base quality are all read off the FIRST read of the
 *     fragment, which is the first in the sample's own order, so which read is consulted depends on
 *     the order the reads were added;
 *   - isF2R1 is isReverseStrand() == isFirstOfPair(), and an unpaired read has isFirstOfPair()
 *     false, so a forward unpaired read lands in F2R1 and a reverse one in F1R2;
 *   - the counts are keyed on the MATRIX's alleles and read back by the VARIANT's, so an allele the
 *     variant declares and the matrix does not is a NullPointerException.
 *
 * Output:
 *
 *     frag\t<label>\t<sample>\t<start>-<end>:<readNames>;...
 *     grouped\t<label>\t<allele>=<bits>,<bits>,...;...
 *     f1r2\t<label>\t<key>=<value>[<class>];...
 *     fad\t<label>\t<key>=<value>[<class>];...
 *
 * Usage: FragmentCountsDump
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

import org.broadinstitute.hellbender.tools.walkers.annotator.FragmentDepthPerAlleleBySample;
import org.broadinstitute.hellbender.tools.walkers.annotator.OrientationBiasReadCounts;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class FragmentCountsDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# FragmentCountsDump: F1R2, F2R1 and FAD");

        for (final String label : new String[] {
                "paired-f1r2", "paired-f2r1", "unpaired-forward", "unpaired-reverse",
                "mixed", "low-mapping-quality", "unavailable-mapping-quality",
                "low-base-quality", "singleton-and-pair", "three-reads-one-name",
                "second-read-low-quality", "empty"}) {
            fragments(label);
            grouped(label);
            f1r2(label);
            fad(label);
        }
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    /** name, flags, mapping quality, base quality, best allele index. */
    static Object[][] composition(final String label) {
        switch (label) {
            case "paired-f1r2": return new Object[][] {
                    {"p", 0x41, 60, 30, 0}, {"p", 0x81, 60, 30, 0},
                    {"q", 0x51, 60, 30, 1}, {"q", 0x91, 60, 30, 1}};
            case "paired-f2r1": return new Object[][] {
                    {"p", 0x51, 60, 30, 0}, {"p", 0x91, 60, 30, 0},
                    {"q", 0x41, 60, 30, 1}, {"q", 0x81, 60, 30, 1}};
            case "unpaired-forward": return new Object[][] {
                    {"a", 0, 60, 30, 0}, {"b", 0, 60, 30, 1}};
            case "unpaired-reverse": return new Object[][] {
                    {"a", 0x10, 60, 30, 0}, {"b", 0x10, 60, 30, 1}};
            case "mixed": return new Object[][] {
                    {"p", 0x41, 60, 30, 0}, {"p", 0x81, 60, 30, 0},
                    {"q", 0x51, 60, 30, 1}, {"r", 0x41, 60, 30, 1},
                    {"s", 0x91, 60, 30, 0}};
            case "low-mapping-quality": return new Object[][] {
                    {"a", 0x41, 0, 30, 0}, {"b", 0x41, 60, 30, 1}};
            case "unavailable-mapping-quality": return new Object[][] {
                    {"a", 0x41, 255, 30, 0}, {"b", 0x41, 60, 30, 1}};
            case "low-base-quality": return new Object[][] {
                    {"a", 0x41, 60, 5, 0}, {"b", 0x41, 60, 30, 1}};
            case "singleton-and-pair": return new Object[][] {
                    {"p", 0x41, 60, 30, 0}, {"p", 0x81, 60, 30, 0},
                    {"a", 0x41, 60, 30, 1}};
            case "three-reads-one-name": return new Object[][] {
                    {"p", 0x41, 60, 30, 0}, {"p", 0x81, 60, 30, 0},
                    {"p", 0x841, 60, 30, 0}, {"b", 0x41, 60, 30, 1}};
            case "second-read-low-quality": return new Object[][] {
                    {"p", 0x41, 60, 30, 0}, {"p", 0x81, 60, 2, 0},
                    {"b", 0x41, 60, 30, 1}};
            case "empty": return new Object[0][];
            default: throw new IllegalArgumentException(label);
        }
    }

    static AlleleLikelihoods<GATKRead, Allele> readLikelihoods(final String label) {
        final Object[][] composition = composition(label);
        final List<GATKRead> reads = new ArrayList<>();
        final List<Integer> best = new ArrayList<>();
        for (final Object[] row : composition) {
            reads.add(read((String) row[0], (Integer) row[1], (Integer) row[2], (Integer) row[3]));
            best.add((Integer) row[4]);
        }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(REF, ALT), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            for (int a = 0; a < 2; a++) {
                m.set(a, e, a == best.get(e) ? -1 : -10);
            }
        }
        return likelihoods;
    }

    static AlleleLikelihoods<Fragment, Allele> fragmentLikelihoods(final String label) {
        return readLikelihoods(label).groupEvidence(GATKRead::getName, Fragment::createAndAvoidFailure);
    }

    static void fragments(final String label) {
        final AlleleLikelihoods<Fragment, Allele> likelihoods = fragmentLikelihoods(label);
        final StringJoiner joiner = new StringJoiner(";");
        for (final Fragment fragment : likelihoods.sampleEvidence(0)) {
            final StringJoiner names = new StringJoiner(",");
            for (final GATKRead read : fragment.getReads()) {
                names.add(read.getName() + "/" + read.getFlags());
            }
            joiner.add(fragment.getStart() + "-" + fragment.getEnd() + ":" + names);
        }
        System.out.printf("frag\t%s\ts1\t%s%n", label, joiner);
    }

    static void grouped(final String label) {
        final AlleleLikelihoods<Fragment, Allele> likelihoods = fragmentLikelihoods(label);
        final LikelihoodMatrix<Fragment, Allele> m = likelihoods.sampleMatrix(0);
        final StringJoiner joiner = new StringJoiner(";");
        for (int a = 0; a < m.numberOfAlleles(); a++) {
            final StringJoiner values = new StringJoiner(",");
            for (int e = 0; e < m.evidenceCount(); e++) {
                values.add(Long.toString(Double.doubleToRawLongBits(m.get(a, e))));
            }
            joiner.add(m.getAllele(a).getDisplayString() + "=" + values);
        }
        System.out.printf("grouped\t%s\t%s%n", label, joiner);
    }

    static VariantContext variantContext() {
        return new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(List.of(REF, ALT))
                .genotypes(List.of(new GenotypeBuilder("s1", List.of(REF, ALT)).make())).make();
    }

    static void f1r2(final String label) {
        final VariantContext vc = variantContext();
        final Genotype g = vc.getGenotype("s1");
        final GenotypeBuilder gb = new GenotypeBuilder(g);
        try {
            new OrientationBiasReadCounts().annotate(null, null, vc, g, gb,
                    readLikelihoods(label), fragmentLikelihoods(label), null);
            emit("f1r2", label, gb.make());
        } catch (final Exception | AssertionError e) {
            System.out.printf("f1r2\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void fad(final String label) {
        final VariantContext vc = variantContext();
        final Genotype g = vc.getGenotype("s1");
        final GenotypeBuilder gb = new GenotypeBuilder(g);
        try {
            new FragmentDepthPerAlleleBySample().annotate(null, null, vc, g, gb,
                    readLikelihoods(label), fragmentLikelihoods(label), null);
            emit("fad", label, gb.make());
        } catch (final Exception | AssertionError e) {
            System.out.printf("fad\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void emit(final String kind, final String label, final Genotype g) {
        final StringJoiner joiner = new StringJoiner(";");
        final List<String> keys = new ArrayList<>(g.getExtendedAttributes().keySet());
        keys.sort(String::compareTo);
        for (final String key : keys) {
            final Object value = g.getExtendedAttribute(key);
            joiner.add(String.format("%s=%s[%s]", key,
                    value instanceof int[] ? Arrays.toString((int[]) value) : value,
                    value.getClass().getName()));
        }
        System.out.printf("%s\t%s\t%s%n", kind, label, joiner);
    }

    static GATKRead read(final String name, final int flags, final int mappingQuality,
                         final int baseQuality) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setMappingQuality(mappingQuality);
        record.setFlags(flags | record.getFlags());
        return new SAMRecordToGATKReadAdapter(record);
    }
}
