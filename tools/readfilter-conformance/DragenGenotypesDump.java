/*
 * The DRAGEN genotyping models HaplotypeCaller switches on with --apply-bqd and --apply-frd
 * (G4.1, #1330): DRAGENGenotypesModel.calculateLikelihoods over hand-built read likelihoods, with
 * GenotypeLikelihoodCalculatorDRAGEN's base quality dropout (BQD) and foreign read detection (FRD)
 * models, FRDBQDUtils' homopolymer adjustments, DRAGENMappingQualityReadTransformer's mapping
 * quality table, and the DragstrReferenceAnalyzer the indel prior is read through.
 *
 * Each case is computed three ways, BQD alone, FRD alone and both, and each sample's log10
 * genotype likelihoods are printed as raw bits followed by the PLs they store as.
 *
 * Output:
 *
 *     reference\t<bases>\t<start>
 *     mq\t<mapping quality>\t<score bits>             (or mq\t<mq>\terror\t<class>: <message>)
 *     xq\t<label>\t<mapping quality after the transformer>
 *     str\t<label>\t<bases>\t<start>\t<end>\t<max period>\t<position>\t<period>\t<repeat length>
 *     strerror\t<label>\t<class>: <message>
 *     case\t<label>\t<alleles>\t<ploidy>\t<variant>\t<padding>\t<offset>\t<str>\t<max effective depth>\t<samples>
 *     read\t<label>\t<sample>\t<name>\t<start>\t<cigar>\t<reverse>\t<mq>\t<quals,...>
 *     lk\t<label>\t<sample index>\t<allele index>\t<values,...>
 *     filtered\t<label>\t<sample index>\t<names,...>
 *     gl\t<label>\t<BQD|FRD|BOTH>\t<sample index>\t<bits,...>\t<PLs,...>
 *     error\t<label>\t<BQD|FRD|BOTH>\t<class>: <message>
 *
 * The reference is REF at chr1:1000; reads are 30 bases whose names say their spec.
 *
 * Usage: DragenGenotypesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.tools.walkers.genotyper.DRAGENGenotypesModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingData;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypingLikelihoods;
import org.broadinstitute.hellbender.tools.walkers.genotyper.HomogeneousPloidyModel;
import org.broadinstitute.hellbender.transformers.DRAGENMappingQualityReadTransformer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.dragstr.DragstrReferenceAnalyzer;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class DragenGenotypesDump {

    // The G at offset 30 sits between four A's and five, so an A error allele there carries the
    // full homopolymer adjustment on both strands; a CA repeat runs from offset 43.
    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTAAAAGAAAAATCGATCGCACACACACAGTCCATGACTTGA";
    static final int START = 1000;
    static final int LENGTH = 30;

    /** start:cigar:strand(F|R):mq:support:strength[:X], X marking a read the HMM disqualified. */
    record Read(String sample, int start, String cigar, boolean reverse, int mq, int support, int strength, boolean filtered) {
        static Read of(final String sample, final String spec) {
            final String[] f = spec.split(":");
            return new Read(sample, Integer.parseInt(f[0]), f[1], f[2].equals("R"), Integer.parseInt(f[3]),
                    Integer.parseInt(f[4]), Integer.parseInt(f[5]), f.length > 6 && f[6].equals("X"));
        }
    }

    public static void main(final String[] args) {
        System.out.println("# DragenGenotypesDump: DRAGENGenotypesModel, GenotypeLikelihoodCalculatorDRAGEN, FRDBQDUtils");
        System.out.println("reference\t" + REF + "\t" + START);
        mappingQualityTable();
        transformer();
        strAnalyzer();

        final List<String> mixed = List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:1:6", "1012:30M:F:50:0:6", "1018:30M:R:60:1:5",
                "1020:30M:F:60:0:6", "1005:30M:R:20:1:4", "1022:30M:F:60:1:6", "1009:30M:R:3:0:6");
        run("snp-het", "G,A", 2, 1030, 2, false, 0, Map.of("sA", mixed));
        run("snp-het-depth-cap", "G,A", 2, 1030, 2, false, 3, Map.of("sA", mixed));
        run("snp-het-depth-cap-high", "G,A", 2, 1030, 2, false, 100, Map.of("sA", mixed));
        run("snp-het-padding-0", "G,A", 2, 1030, 0, false, 0, Map.of("sA", mixed));

        // Alternate support confined to the reads' ends: what BQD exists to catch.
        run("bqd-dropout", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1001:30M:F:60:1:2", "1002:30M:F:60:1:2", "1010:30M:F:60:0:6", "1015:30M:F:60:0:6",
                "1020:30M:F:60:0:6", "1025:30M:F:60:0:6", "1030:30M:R:60:0:6", "1029:30M:R:60:0:6",
                "1028:30M:R:60:1:2", "1027:30M:R:60:0:6")));
        run("bqd-error-c", "G,C", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1001:30M:F:60:1:2", "1002:30M:F:60:1:2", "1010:30M:F:60:0:6", "1015:30M:F:60:0:6",
                "1030:30M:R:60:0:6", "1028:30M:R:60:1:2")));
        run("bqd-forward-only", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1001:30M:F:60:1:2", "1003:30M:F:60:1:3", "1010:30M:F:60:0:6", "1015:30M:F:60:0:6",
                "1020:30M:F:60:0:6")));
        run("bqd-ties", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1010:30M:F:60:0:6", "1010:30M:F:60:1:3", "1010:30M:R:60:0:6", "1010:30M:R:60:1:3",
                "1012:30M:F:60:0:6", "1012:30M:R:60:1:3")));

        // Alternate support confined to reads of poor mapping quality: what FRD exists to catch.
        run("frd-foreign", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:0:6", "1012:30M:F:60:0:6", "1018:30M:R:60:0:6",
                "1020:30M:F:5:1:6", "1005:30M:R:10:1:6", "1022:30M:F:2:1:6", "1009:30M:R:0:1:6")));
        run("frd-foreign-depth-cap", "G,A", 2, 1030, 2, false, 4, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:0:6", "1012:30M:F:60:0:6", "1018:30M:R:60:0:6",
                "1020:30M:F:5:1:6", "1005:30M:R:10:1:6", "1022:30M:F:2:1:6", "1009:30M:R:0:1:6")));
        run("frd-extended-mq", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:250:0:6", "1016:30M:R:120:0:6", "1012:30M:F:80:0:6", "1018:30M:R:200:0:6",
                "1020:30M:F:45:1:6", "1005:30M:R:15:1:6")));

        // Clips: the feather ends count soft and hard clipped bases differently.
        run("clipped", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1017:5S25M:F:60:1:2", "1012:25M5S:F:60:0:6", "1020:3H27M:F:60:1:3", "1008:27M3H:R:60:0:6",
                "1006:2S28M:R:60:1:2", "1014:28M2S:R:60:0:6", "1013:30M:F:60:0:6")));
        // A read deleting the variant, one whose soft clip alone covers it, one ending before it.
        run("not-covering", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1020:10M2D20M:F:60:0:6", "1031:5S25M:R:60:0:6", "990:30M:F:60:0:6", "1010:30M:F:60:1:3",
                "1015:30M:R:60:1:3", "1012:30M:F:60:0:6")));
        // Reads the HMM disqualified are sorted with the others but carry no likelihood.
        run("disqualified", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:1:6", "1003:30M:F:60:1:30:X", "1004:30M:R:20:0:30:X",
                "1020:30M:F:60:0:6", "1018:30M:R:60:1:5", "1011:30M:F:60:0:30:X")));
        run("all-disqualified-strand", "G,A", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:F:60:1:6", "1003:30M:R:60:1:30:X", "1004:30M:R:20:0:30:X")));

        // Indels: BQD leaves them alone, FRD reads the indel prior.
        run("insertion-flat-prior", "G,GA", 2, 1030, 2, false, 0, Map.of("sA", mixed));
        run("insertion-dragstr", "G,GA", 2, 1030, 2, true, 0, Map.of("sA", mixed));
        run("deletion-dragstr-repeat", "AC,A", 2, 1044, 2, true, 0, Map.of("sA", List.of(
                "1025:30M:F:60:0:6", "1026:30M:R:60:1:6", "1022:30M:F:60:0:6", "1028:30M:R:30:1:5",
                "1030:30M:F:8:1:6", "1015:30M:R:60:0:6")));
        run("snp-dragstr", "G,A", 2, 1030, 2, true, 0, Map.of("sA", mixed));

        // More alleles and other ploidies.
        run("triallelic", "G,A,T", 2, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:1:6", "1012:30M:F:50:2:6", "1018:30M:R:60:1:5",
                "1020:30M:F:4:2:6", "1005:30M:R:20:1:4", "1022:30M:F:60:2:6", "1009:30M:R:3:0:6")));
        run("triallelic-indel", "G,A,GA", 2, 1030, 2, true, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:1:6", "1012:30M:F:50:2:6", "1018:30M:R:60:1:5",
                "1020:30M:F:4:2:6", "1005:30M:R:20:1:4")));
        run("haploid", "G,A", 1, 1030, 2, false, 0, Map.of("sA", mixed));
        run("triploid", "G,A", 3, 1030, 2, false, 5, Map.of("sA", mixed));
        run("triploid-triallelic", "G,A,T", 3, 1030, 2, false, 0, Map.of("sA", List.of(
                "1015:30M:F:60:0:6", "1016:30M:R:60:1:6", "1012:30M:F:50:2:6", "1018:30M:R:60:1:5",
                "1020:30M:F:4:2:6", "1005:30M:R:20:1:4")));

        final Map<String, List<String>> two = new LinkedHashMap<>();
        two.put("sA", mixed);
        two.put("sB", List.of("1001:30M:F:60:1:2", "1002:30M:F:60:1:2", "1010:30M:F:60:0:6", "1015:30M:R:9:1:6",
                "1020:30M:R:60:0:6"));
        run("two-samples", "G,A", 2, 1030, 2, false, 0, two);
        final Map<String, List<String>> empty = new LinkedHashMap<>();
        empty.put("sA", mixed);
        empty.put("sB", List.of());
        run("empty-sample", "G,A", 2, 1030, 2, false, 0, empty);
    }

    static String bits(final double d) {
        return Long.toHexString(Double.doubleToRawLongBits(d));
    }

    static void mappingQualityTable() {
        for (int mq = -1; mq <= 257; mq++) {
            try {
                System.out.println("mq\t" + mq + "\t" + bits(DRAGENMappingQualityReadTransformer.mapMappingQualityToPhredLikelihoodScore(mq)));
            } catch (final Exception e) {
                System.out.println("mq\t" + mq + "\terror\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
    }

    static void transformer() {
        final SAMFileHeader header = header(List.of("sA"));
        final DRAGENMappingQualityReadTransformer t = new DRAGENMappingQualityReadTransformer();
        final GATKRead tagged = ArtificialReadUtils.createArtificialRead(header, "tagged", 0, 1000, 10);
        tagged.setMappingQuality(60);
        tagged.setAttribute("XQ", 187);
        System.out.println("xq\ttagged\t" + t.apply(tagged).getMappingQuality());
        final GATKRead plain = ArtificialReadUtils.createArtificialRead(header, "plain", 0, 1000, 10);
        plain.setMappingQuality(42);
        System.out.println("xq\tplain\t" + t.apply(plain).getMappingQuality());
        final GATKRead zero = ArtificialReadUtils.createArtificialRead(header, "zero", 0, 1000, 10);
        zero.setMappingQuality(60);
        zero.setAttribute("XQ", 0);
        System.out.println("xq\tzero\t" + t.apply(zero).getMappingQuality());
    }

    static void strAnalyzer() {
        final String[][] specs = {
                {"ref", REF, "0", String.valueOf(REF.length()), "8"},
                {"ref-window", REF, "20", "50", "8"},
                {"ref-period-2", REF, "38", "52", "2"},
                {"ref-period-1", REF, "20", "40", "1"},
                {"trinucleotide", "TTCAGCAGCAGCAGCAGTTACGACGACGTT", "0", "30", "8"},
                {"run-to-end", "ACGTAAAAAAAA", "2", "8", "8"},
                {"short", "ACAC", "0", "4", "8"},
                {"empty-window", REF, "10", "10", "8"},
        };
        for (final String[] s : specs) {
            final byte[] bases = s[1].getBytes(StandardCharsets.US_ASCII);
            final int start = Integer.parseInt(s[2]);
            final int end = Integer.parseInt(s[3]);
            final int maxPeriod = Integer.parseInt(s[4]);
            try {
                final DragstrReferenceAnalyzer a = DragstrReferenceAnalyzer.of(bases, start, end, maxPeriod);
                for (int p = start; p < end; p++) {
                    System.out.println("str\t" + s[0] + "\t" + s[1] + "\t" + start + "\t" + end + "\t" + maxPeriod + "\t" + p
                            + "\t" + a.period(p) + "\t" + a.repeatLength(p));
                }
            } catch (final Exception e) {
                System.out.println("strerror\t" + s[0] + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
        for (final int[] bad : new int[][]{{-1, 5}, {5, 4}, {0, REF.length() + 1}}) {
            try {
                DragstrReferenceAnalyzer.of(REF.getBytes(StandardCharsets.US_ASCII), bad[0], bad[1], 8);
                System.out.println("strerror\tbad-" + bad[0] + "-" + bad[1] + "\tnone");
            } catch (final Exception e) {
                System.out.println("strerror\tbad-" + bad[0] + "-" + bad[1] + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
    }

    static SAMFileHeader header(final List<String> samples) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 5000))));
        for (final String s : samples) {
            final SAMReadGroupRecord g = new SAMReadGroupRecord("rg" + s);
            g.setSample(s);
            header.addReadGroup(g);
        }
        return header;
    }

    /** Base qualities that move along the read, so each read's quality at the variant differs. */
    static byte[] quals(final int seed, final int length) {
        final byte[] q = new byte[length];
        for (int i = 0; i < length; i++) {
            q[i] = (byte) (12 + ((i * 7 + seed * 11) % 29));
        }
        return q;
    }

    static int readLength(final String cigar) {
        int total = 0;
        int n = 0;
        for (final char c : cigar.toCharArray()) {
            if (Character.isDigit(c)) {
                n = n * 10 + (c - '0');
            } else {
                if (c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X') {
                    total += n;
                }
                n = 0;
            }
        }
        return total;
    }

    static void run(final String label, final String alleleSpec, final int ploidy, final int variant, final int padding,
                    final boolean str, final int maxEffectiveDepth, final Map<String, List<String>> readsBySample) {
        final List<String> samples = new ArrayList<>(readsBySample.keySet());
        samples.sort(null);
        final List<Allele> alleles = new ArrayList<>();
        final String[] names = alleleSpec.split(",");
        for (int i = 0; i < names.length; i++) {
            alleles.add(Allele.create(names[i], i == 0));
        }
        final int offset = variant - START;
        System.out.println("case\t" + label + "\t" + alleleSpec + "\t" + ploidy + "\t" + variant + "\t" + padding + "\t" + offset
                + "\t" + str + "\t" + maxEffectiveDepth + "\t" + String.join(",", samples));
        final SAMFileHeader header = header(samples);
        final Map<String, List<GATKRead>> evidence = new LinkedHashMap<>();
        final Map<String, List<Read>> specs = new LinkedHashMap<>();
        int count = 0;
        for (final String s : samples) {
            evidence.put(s, new ArrayList<>());
            specs.put(s, new ArrayList<>());
            for (final String spec : readsBySample.get(s)) {
                final Read r = Read.of(s, spec);
                final String name = s + "r" + count++;
                final int length = readLength(r.cigar());
                final byte[] bases = new byte[length];
                Arrays.fill(bases, (byte) 'A');
                final byte[] q = quals(count, length);
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, r.start(), bases, q, r.cigar());
                read.setIsReverseStrand(r.reverse());
                read.setMappingQuality(r.mq());
                read.setReadGroup("rg" + s);
                System.out.println("read\t" + label + "\t" + s + "\t" + name + "\t" + r.start() + "\t" + r.cigar() + "\t" + r.reverse()
                        + "\t" + r.mq() + "\t" + Arrays.toString(q).replaceAll("[\\[\\] ]", ""));
                evidence.get(s).add(read);
                specs.get(s).add(r);
            }
        }
        for (final String config : List.of("BQD", "FRD", "BOTH")) {
            try {
                final AlleleLikelihoods<GATKRead, Allele> likelihoods = new AlleleLikelihoods<>(new IndexedSampleList(samples),
                        new IndexedAlleleList<>(alleles), evidence.entrySet().stream()
                        .collect(Collectors.toMap(Map.Entry::getKey, e -> new ArrayList<>(e.getValue()), (a, b) -> a, LinkedHashMap::new)));
                for (int s = 0; s < samples.size(); s++) {
                    final LikelihoodMatrix<GATKRead, Allele> m = likelihoods.sampleMatrix(s);
                    final List<Read> sampleSpecs = specs.get(samples.get(s));
                    for (int a = 0; a < alleles.size(); a++) {
                        final StringBuilder b = new StringBuilder();
                        for (int r = 0; r < sampleSpecs.size(); r++) {
                            final Read spec = sampleSpecs.get(r);
                            final double v = spec.support() == a ? -0.5 - 0.125 * (r % 3) - (spec.filtered() ? spec.strength() : 0)
                                    : -0.5 - spec.strength() - 0.25 * ((a + r) % 4);
                            m.set(a, r, v);
                            b.append(r == 0 ? "" : ",").append(bits(v));
                        }
                        if (config.equals("BQD")) {
                            System.out.println("lk\t" + label + "\t" + s + "\t" + a + "\t" + b);
                        }
                    }
                }
                likelihoods.filterPoorlyModeledEvidence(read -> -20.0);
                if (config.equals("BQD")) {
                    for (int s = 0; s < samples.size(); s++) {
                        System.out.println("filtered\t" + label + "\t" + s + "\t"
                                + likelihoods.filteredSampleEvidence(s).stream().map(GATKRead::getName).collect(Collectors.joining(",")));
                    }
                }
                likelihoods.setVariantCallingSubsetUsed(new SimpleInterval("chr1", variant - padding, variant + padding));
                final byte[] ref = REF.getBytes(StandardCharsets.US_ASCII);
                final DragstrReferenceAnalyzer dragstrs = str
                        ? DragstrReferenceAnalyzer.of(ref, 0, ref.length, DragstrParams.DEFAULT.maximumPeriod()) : null;
                final DRAGENGenotypesModel model = new DRAGENGenotypesModel(!config.equals("FRD"), !config.equals("BQD"),
                        padding, maxEffectiveDepth, DragstrParams.DEFAULT);
                final GenotypingLikelihoods<Allele> gls = model.calculateLikelihoods(new IndexedAlleleList<>(alleles),
                        new GenotypingData<>(new HomogeneousPloidyModel(new IndexedSampleList(samples), ploidy), likelihoods),
                        ref, offset, dragstrs);
                for (int s = 0; s < samples.size(); s++) {
                    final GenotypeLikelihoods gl = gls.sampleLikelihoods(s);
                    final String values = Arrays.stream(gl.getAsVector()).mapToObj(DragenGenotypesDump::bits).collect(Collectors.joining(","));
                    final String pls = Arrays.stream(gl.getAsPLs()).mapToObj(Integer::toString).collect(Collectors.joining(","));
                    System.out.println("gl\t" + label + "\t" + config + "\t" + s + "\t" + values + "\t" + pls);
                }
            } catch (final Exception e) {
                System.out.println("error\t" + label + "\t" + config + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            }
        }
    }
}
