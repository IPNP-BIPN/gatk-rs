/*
 * AS_FisherStrand, AS_StrandOddsRatio and the StrandBiasUtils machinery under them, taken from the
 * reference.
 *
 *   - AS_SB_TABLE writes an entry for the reference too, and puts the delimiter BETWEEN entries, so
 *     it does not start with one. The rank-sum family skips the reference's value but keeps its
 *     slot, so theirs does. Two allele-specific families, two conventions for the same delimiter;
 *   - a sample contributes nothing unless its whole table holds MORE THAN two informative reads,
 *     so a sample with exactly two is dropped even for the allele that had both;
 *   - AS_StrandOddsRatio computes a value for the reference allele against itself and then never
 *     prints it; AS_FisherStrand filters the reference out first. The two differ in that one line;
 *   - AS_FS floors its p-value at 1e-320, which is a subnormal, before phred-scaling it;
 *   - a reduced string skips an alternate the data does not carry rather than writing the missing
 *     value, so the field can come out shorter than the alternate count.
 *
 * Output:
 *
 *     sbraw\t<label>\t<key>=<value>[<class>];...
 *     sbdirect\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     sbcombine\t<label>\t<key>=<value>[<class>];...
 *     sbfinal\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     phred\t<errorRate>\t<bits>
 *
 * Usage: AlleleSpecificStrandBiasDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_FisherStrand;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandBiasTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_StrandOddsRatio;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.QualityUtils;
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

public class AlleleSpecificStrandBiasDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# AlleleSpecificStrandBiasDump: AS_FS, AS_SOR and AS_SB_TABLE");

        // The raw table, over matrices with different strand compositions.
        for (final String label : new String[] {
                "balanced", "skewed", "ref-only", "alt-only", "one-read", "two-reads",
                "three-reads", "two-samples-each-small", "two-samples-one-large", "multiallelic",
                "all-forward", "all-reverse", "empty", "null-likelihoods"}) {
            raw(label);
            direct("AS_FisherStrand", label);
            direct("AS_StrandOddsRatio", label);
        }

        // Combining and finalising raw strings.
        for (final String label : new String[] {
                "one-source", "two-sources", "three-alleles", "empty-entry", "zero-entry",
                "bracketed", "spaced", "wrong-count", "extreme"}) {
            combineAndFinalize(label);
        }

        // The phred scale, whose floor is the one AS_FS clamps against.
        for (final double rate : new double[] {
                1.0, 0.5, 0.1, 1.0E-10, 1.0E-300, 1.0E-320, 4.9E-324, 0.0}) {
            phred(rate);
        }
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    /** Per sample, the forward/reverse counts for the reference and for each alternate. */
    static int[][][] composition(final String label) {
        switch (label) {
            case "balanced": return new int[][][] {{{5, 5}, {5, 5}}};
            case "skewed": return new int[][][] {{{10, 0}, {0, 10}}};
            case "ref-only": return new int[][][] {{{5, 5}, {0, 0}}};
            case "alt-only": return new int[][][] {{{0, 0}, {5, 5}}};
            case "one-read": return new int[][][] {{{1, 0}, {0, 0}}};
            case "two-reads": return new int[][][] {{{1, 1}, {0, 0}}};
            case "three-reads": return new int[][][] {{{2, 1}, {0, 0}}};
            case "two-samples-each-small": return new int[][][] {{{1, 1}, {0, 0}}, {{1, 1}, {0, 0}}};
            case "two-samples-one-large": return new int[][][] {{{1, 1}, {0, 0}}, {{5, 5}, {5, 5}}};
            case "multiallelic": return new int[][][] {{{6, 6}, {3, 3}, {2, 2}}};
            case "all-forward": return new int[][][] {{{8, 0}, {6, 0}}};
            case "all-reverse": return new int[][][] {{{0, 8}, {0, 6}}};
            case "empty": return new int[][][] {{{0, 0}, {0, 0}}};
            case "null-likelihoods": return null;
            default: throw new IllegalArgumentException(label);
        }
    }

    static List<Allele> allelesFor(final String label) {
        return "multiallelic".equals(label) ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
    }

    static VariantContext variantContext(final String label) {
        final List<Allele> alleles = allelesFor(label);
        final int[][][] composition = composition(label);
        final int samples = composition == null ? 1 : composition.length;
        final List<htsjdk.variant.variantcontext.Genotype> genotypes = new ArrayList<>();
        for (int s = 0; s < samples; s++) {
            genotypes.add(new GenotypeBuilder("s" + s, List.of(REF, ALT)).make());
        }
        return new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(alleles).genotypes(genotypes).make();
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoods(final String label) {
        final int[][][] composition = composition(label);
        if (composition == null) {
            return null;
        }
        final List<Allele> alleles = allelesFor(label);
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        final List<String> samples = new ArrayList<>();
        final List<List<Integer>> bestPerSample = new ArrayList<>();
        for (int s = 0; s < composition.length; s++) {
            final List<GATKRead> reads = new ArrayList<>();
            final List<Integer> best = new ArrayList<>();
            for (int a = 0; a < composition[s].length; a++) {
                for (int strand = 0; strand < 2; strand++) {
                    for (int i = 0; i < composition[s][a][strand]; i++) {
                        reads.add(read("s" + s + "a" + a + "d" + strand + "i" + i, strand == 1));
                        best.add(a);
                    }
                }
            }
            samples.add("s" + s);
            bySample.put("s" + s, reads);
            bestPerSample.add(best);
        }
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(samples), new IndexedAlleleList<>(alleles), bySample);
        for (int s = 0; s < composition.length; s++) {
            final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(s);
            final List<Integer> best = bestPerSample.get(s);
            for (int e = 0; e < best.size(); e++) {
                for (int a = 0; a < alleles.size(); a++) {
                    m.set(a, e, a == best.get(e) ? -1 : -10);
                }
            }
        }
        return likelihoods;
    }

    static void raw(final String label) {
        final VariantContext vc = variantContext(label);
        try {
            emitMap("sbraw", null, label,
                    new AS_FisherStrand().annotateRawData(null, vc, likelihoods(label)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("sbraw\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static void direct(final String name, final String label) {
        final AS_StrandBiasTest annotation = "AS_FisherStrand".equals(name)
                ? new AS_FisherStrand() : new AS_StrandOddsRatio();
        final VariantContext vc = variantContext(label);
        try {
            emitMap("sbdirect", name, label, annotation.annotate(null, vc, likelihoods(label)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("sbdirect\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
        }
    }

    static String[] rawStrings(final String label) {
        switch (label) {
            case "one-source": return new String[] {"10,8|3,4"};
            case "two-sources": return new String[] {"10,8|3,4", "2,3|1,1"};
            case "three-alleles": return new String[] {"10,8|3,4|2,2"};
            case "empty-entry": return new String[] {"10,8|"};
            case "zero-entry": return new String[] {"10,8|0,0"};
            case "bracketed": return new String[] {"[10,8|3, 4]"};
            case "spaced": return new String[] {"10, 8|3, 4"};
            case "wrong-count": return new String[] {"10,8"};
            case "extreme": return new String[] {"4000,1|1,4000"};
            default: throw new IllegalArgumentException(label);
        }
    }

    static void combineAndFinalize(final String label) {
        final List<Allele> alleles = "three-alleles".equals(label)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final AS_FisherStrand fs = new AS_FisherStrand();
        final List<ReducibleAnnotationData<?>> data = new ArrayList<>();
        for (final String raw : rawStrings(label)) {
            data.add(new AlleleSpecificAnnotationData<List<Integer>>(alleles, raw));
        }
        String combined = null;
        try {
            final Map<String, Object> result = fs.combineRawData(alleles, data);
            combined = result.values().iterator().next().toString();
            emitMap("sbcombine", null, label, result);
        } catch (final Exception | AssertionError e) {
            System.out.printf("sbcombine\t%s\tE:%s%n", label, e.getClass().getName());
        }
        if (combined == null) {
            return;
        }
        final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START).stop(START)
                .alleles(alleles).attribute("AS_SB_TABLE", combined).make();
        for (final String name : new String[] {"AS_FisherStrand", "AS_StrandOddsRatio"}) {
            final AS_StrandBiasTest annotation = "AS_FisherStrand".equals(name)
                    ? new AS_FisherStrand() : new AS_StrandOddsRatio();
            try {
                emitMap("sbfinal", name, label, annotation.finalizeRawData(vc, vc));
            } catch (final Exception | AssertionError e) {
                System.out.printf("sbfinal\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
            }
        }
    }

    static void phred(final double errorRate) {
        try {
            // Raw bits, because this is the value AS_FS writes after a %.3f and the last ulp is
            // what a format string would hide.
            System.out.printf("phred\t%s\t%d%n", Double.toString(errorRate),
                    Double.doubleToRawLongBits(QualityUtils.phredScaleErrorRate(errorRate)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("phred\t%s\tE:%s%n", Double.toString(errorRate),
                    e.getClass().getName());
        }
    }

    static void emitMap(final String kind, final String name, final String label,
                        final Map<String, Object> result) {
        final StringJoiner joiner = new StringJoiner(";");
        if (result != null) {
            final List<String> keys = new ArrayList<>(result.keySet());
            keys.sort(String::compareTo);
            for (final String key : keys) {
                final Object value = result.get(key);
                joiner.add(String.format("%s=%s[%s]", key, value, value.getClass().getName()));
            }
        }
        if (name == null) {
            System.out.printf("%s\t%s\t%s%n", kind, label, joiner);
        } else {
            System.out.printf("%s\t%s\t%s\t%s%n", kind, name, label, joiner);
        }
    }

    static GATKRead read(final String name, final boolean reverse) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        record.setReadNegativeStrandFlag(reverse);
        return new SAMRecordToGATKReadAdapter(record);
    }
}
