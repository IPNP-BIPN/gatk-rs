/*
 * AS_BaseQualityRankSumTest, AS_MappingQualityRankSumTest and AS_ReadPosRankSumTest, plus the
 * Histogram and CompressedDataList they carry through a gVCF, taken from the reference.
 *
 *   - the direct annotate() path of an allele-specific annotation is NOT allele-specific: it pools
 *     every alternate's reads into one series and reports one Z score, exactly as its
 *     non-allele-specific parent does, under a different key;
 *   - the raw string starts with its delimiter, because the reference allele is skipped but its
 *     slot is not, and the parser reads the slots positionally;
 *   - a site with no unambiguous reference read produces a raw string of bare delimiters;
 *   - each Z score is stored as a one-entry histogram binned to a tenth by a FLOOR, so a negative
 *     score bins away from zero and the binning is not symmetric about it;
 *   - Histogram.median walks the bins, uses the lower of the two middle positions for an even
 *     count, and answers null when empty;
 *   - an empty Histogram prints as the four characters NaN, not as the empty string.
 *
 * Output:
 *
 *     hist\t<label>\t<toString>\t<median or null>
 *     cdl\t<label>\t<toString>\t<iteration>
 *     bin\t<value>\t<count in bin 0 or null>
 *     as\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     asraw\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     ascombine\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *     asfinal\t<annotation>\t<label>\t<key>=<value>[<class>];...
 *
 * Usage: AlleleSpecificRankSumDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.InfoFieldAnnotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_BaseQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_MappingQualityRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_RankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AS_ReadPosRankSumTest;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.AlleleSpecificAnnotationData;
import org.broadinstitute.hellbender.tools.walkers.annotator.allelespecific.ReducibleAnnotationData;
import org.broadinstitute.hellbender.utils.CompressedDataList;
import org.broadinstitute.hellbender.utils.Histogram;
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

public class AlleleSpecificRankSumDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int START = 105;

    public static void main(final String[] args) {
        System.out.println("# AlleleSpecificRankSumDump: AS_BaseQRankSum, AS_MQRankSum, AS_ReadPosRankSum");

        // Histogram, on its own.
        histogram("empty");
        histogram("one-value");
        histogram("two-values");
        histogram("three-values");
        histogram("four-values");
        histogram("negative");
        histogram("straddling-zero");
        histogram("on-a-boundary");
        histogram("just-below-a-boundary");
        histogram("repeated");
        histogram("wide-bins");

        // CompressedDataList, on its own.
        compressed("empty");
        compressed("ascending");
        compressed("descending");
        compressed("repeated");
        compressed("negative");

        // The binning rule, one value at a time.
        for (final double value : new double[] {
                0.0, -0.0, 0.05, 0.099, 0.1, -0.001, -0.0009, -0.05, -0.1, 1.0, -1.23, -1.27,
                1.23, 1.27, 0.049, -0.049}) {
            bin(value);
        }

        // The three annotations, over one sample as the raw path requires.
        emit("separated");
        emit("overlapping");
        emit("ref-only");
        emit("alt-only");
        emit("multiallelic");
        emit("single-read");
        emit("null-likelihoods");
        emit("no-genotypes");
        emit("two-samples");

        // Combining and finalising raw strings.
        combineAndFinalize("one-source", new String[] {"|-1.3,1|0.4,1"});
        combineAndFinalize("two-sources", new String[] {"|-1.3,1|0.4,1", "|-1.3,1|0.5,2"});
        combineAndFinalize("even-count", new String[] {"|1.0,1", "|2.0,1"});
        combineAndFinalize("odd-count", new String[] {"|1.0,1", "|2.0,2"});
        combineAndFinalize("empty-slots", new String[] {"||"});
        combineAndFinalize("bracketed", new String[] {"[|-1.3,1|0.4,1]"});
        combineAndFinalize("missing-alt", new String[] {"|-1.3,1"});
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static double[] histogramValues(final String label) {
        switch (label) {
            case "empty": return new double[0];
            case "one-value": return new double[] {0.5};
            case "two-values": return new double[] {0.1, 0.2};
            case "three-values": return new double[] {0.1, 0.2, 0.3};
            case "four-values": return new double[] {0.1, 0.2, 0.3, 0.4};
            case "negative": return new double[] {-1.5, -2.5, -3.5};
            case "straddling-zero": return new double[] {-0.2, -0.1, 0.1, 0.2};
            case "on-a-boundary": return new double[] {0.1, 0.2, 0.30000000000000004};
            case "just-below-a-boundary": return new double[] {0.0999, 0.1999};
            case "repeated": return new double[] {0.5, 0.5, 0.5, 0.5, 0.5};
            case "wide-bins": return new double[] {1.0, 2.0, 3.0};
            default: throw new IllegalArgumentException(label);
        }
    }

    static void histogram(final String label) {
        // "wide-bins" is the only one that does not use the default bin size.
        final Histogram h = "wide-bins".equals(label) ? new Histogram(0.01) : new Histogram();
        try {
            for (final double value : histogramValues(label)) {
                h.add(value);
            }
            final Double median = h.median();
            System.out.printf("hist\t%s\t%s\t%s%n", label, h,
                    median == null ? "null" : Double.toString(median));
        } catch (final Exception | AssertionError e) {
            System.out.printf("hist\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }

    static int[] compressedValues(final String label) {
        switch (label) {
            case "empty": return new int[0];
            case "ascending": return new int[] {1, 2, 3};
            case "descending": return new int[] {3, 2, 1};
            case "repeated": return new int[] {2, 2, 2, 5, 5};
            case "negative": return new int[] {-3, 0, 3};
            default: throw new IllegalArgumentException(label);
        }
    }

    static void compressed(final String label) {
        final CompressedDataList<Integer> list = new CompressedDataList<>();
        for (final int value : compressedValues(label)) {
            list.add(value);
        }
        final StringJoiner joiner = new StringJoiner(",");
        for (final Integer value : list) {
            joiner.add(value.toString());
        }
        System.out.printf("cdl\t%s\t%s\t%s%n", label, list, joiner);
    }

    static void bin(final double value) {
        final Histogram h = new Histogram();
        try {
            h.add(value);
            // The bin the value landed in, reported by asking for the count at bin zero and at the
            // value itself, plus the rendering, which names the bin's centre.
            System.out.printf("bin\t%s\t%s\t%s%n", Double.toString(value), h,
                    h.get(0.0) == null ? "null" : h.get(0.0).toString());
        } catch (final Exception | AssertionError e) {
            System.out.printf("bin\t%s\tE:%s%n", Double.toString(value), e.getClass().getName());
        }
    }

    static VariantContext variantContext(final String label) {
        final List<Allele> alleles = "multiallelic".equals(label)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final VariantContextBuilder builder = new VariantContextBuilder().chr("chr1")
                .start(START).stop(START).alleles(alleles);
        if ("no-genotypes".equals(label)) {
            return builder.make();
        }
        if ("two-samples".equals(label)) {
            return builder.genotypes(List.of(
                    new GenotypeBuilder("s1", List.of(REF, ALT)).make(),
                    new GenotypeBuilder("s2", List.of(REF, ALT)).make())).make();
        }
        return builder.genotypes(List.of(new GenotypeBuilder("s1", List.of(REF, ALT)).make())).make();
    }

    static AlleleLikelihoods<GATKRead, Allele> likelihoods(final String label) {
        if ("null-likelihoods".equals(label)) {
            return null;
        }
        final int count = "single-read".equals(label) ? 1 : 12;
        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            reads.add(read("r" + i, 20 + i * 3, 30 + i));
        }
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", reads);
        final List<Allele> alleles = "multiallelic".equals(label)
                ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")),
                new IndexedAlleleList<>(alleles), bySample);
        final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(0);
        for (int e = 0; e < reads.size(); e++) {
            final int best;
            if ("ref-only".equals(label)) {
                best = 0;
            } else if ("alt-only".equals(label)) {
                best = 1;
            } else if ("multiallelic".equals(label)) {
                best = e % 3;
            } else {
                best = e < reads.size() / 2 ? 0 : 1;
            }
            for (int a = 0; a < alleles.size(); a++) {
                final double strong = "overlapping".equals(label) ? -1 - (e * 0.1) : -1;
                final double weak = "overlapping".equals(label) ? -5 - (e * 0.1) : -10;
                m.set(a, e, a == best ? strong : weak);
            }
        }
        return likelihoods;
    }

    static void emit(final String label) {
        one("AS_BaseQualityRankSumTest", label, new AS_BaseQualityRankSumTest());
        one("AS_MappingQualityRankSumTest", label, new AS_MappingQualityRankSumTest());
        one("AS_ReadPosRankSumTest", label, new AS_ReadPosRankSumTest());
    }

    static void one(final String name, final String label, final AS_RankSumTest annotation) {
        final VariantContext vc = variantContext(label);
        final AlleleLikelihoods<GATKRead, Allele> likelihoods = likelihoods(label);
        try {
            emitMap("as", name, label, ((InfoFieldAnnotation) annotation).annotate(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("as\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
        }
        try {
            emitMap("asraw", name, label, annotation.annotateRawData(null, vc, likelihoods));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asraw\t%s\t%s\tE:%s%n", name, label, e.getClass().getName());
        }
    }

    static void combineAndFinalize(final String label, final String[] rawStrings) {
        final List<Allele> alleles = rawStrings[0].replace("[", "").replace("]", "")
                .split("\\|", -1).length > 2 ? List.of(REF, ALT, ALT2) : List.of(REF, ALT);
        final AS_BaseQualityRankSumTest annotation = new AS_BaseQualityRankSumTest();
        final List<ReducibleAnnotationData<?>> data = new ArrayList<>();
        for (final String raw : rawStrings) {
            data.add(new AlleleSpecificAnnotationData<Histogram>(alleles, raw));
        }
        String combined = null;
        try {
            final Map<String, Object> result = annotation.combineRawData(alleles, data);
            combined = result.values().iterator().next().toString();
            emitMap("ascombine", "AS_BaseQualityRankSumTest", label, result);
        } catch (final Exception | AssertionError e) {
            System.out.printf("ascombine\t%s\t%s\tE:%s%n", "AS_BaseQualityRankSumTest", label,
                    e.getClass().getName());
        }
        if (combined == null) {
            return;
        }
        try {
            final VariantContext vc = new VariantContextBuilder().chr("chr1").start(START)
                    .stop(START).alleles(alleles)
                    .attribute("AS_RAW_BaseQRankSum", combined).make();
            emitMap("asfinal", "AS_BaseQualityRankSumTest", label,
                    annotation.finalizeRawData(vc, vc));
        } catch (final Exception | AssertionError e) {
            System.out.printf("asfinal\t%s\t%s\tE:%s%n", "AS_BaseQualityRankSumTest", label,
                    e.getClass().getName());
        }
    }

    static void emitMap(final String kind, final String name, final String label,
                        final Map<String, Object> result) {
        final StringJoiner joiner = new StringJoiner(";");
        if (result != null) {
            final List<String> keys = new ArrayList<>(result.keySet());
            // The map is a HashMap, so sort for a stable dump; the port compares as a set.
            keys.sort(String::compareTo);
            for (final String key : keys) {
                final Object value = result.get(key);
                joiner.add(String.format("%s=%s[%s]", key, value, value.getClass().getName()));
            }
        }
        System.out.printf("%s\t%s\t%s\t%s%n", kind, name, label, joiner);
    }

    static GATKRead read(final String name, final int mappingQuality, final int baseQuality) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(100);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) Math.min(baseQuality, 60));
        record.setBaseQualities(qualities);
        record.setMappingQuality(Math.min(mappingQuality, 60));
        return new SAMRecordToGATKReadAdapter(record);
    }
}
