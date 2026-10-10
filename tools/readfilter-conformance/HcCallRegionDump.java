/*
 * HaplotypeCallerEngine.callRegion: a whole active region called, as HaplotypeCaller calls it with
 * its defaults (no reference confidence), through a real engine over a temporary FASTA and the
 * default annotations (StandardAnnotation and StandardHCAnnotation).
 *
 * The region is assembled, its variation events trimmed to, its read stubs and non-passing reads
 * (short, mapping quality under 20, mate on another contig) removed, the reads scored against the
 * trimmed haplotypes, realigned to their best haplotype, genotyped, annotated (the filtered reads
 * added to the annotation matrix with no likelihood) and phased. The random generator is reset
 * before each region.
 *
 * Output:
 *
 *     reference\t<chr1 bases>
 *     case\t<label>\t<active start>\t<active end>\t<padding>\t<samples,...>
 *     read\t<label>\t<name>\t<sample>\t<start>\t<cigar>\t<bases>\t<quals phred+33>\t<flags>\t<mapq>
 *     call\t<label>\t<VCF line>
 *     calls\t<label>\t<count>
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: HcCallRegionDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.tools.walkers.annotator.Annotation;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.Utils;
import java.util.Collections;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.SortedSet;
import java.util.stream.Collectors;

public class HcCallRegionDump {

    record Spec(String sample, int start, String cigar, String bases, String quals, int flags, int mapq) {
    }

    static String ref;
    static Path fasta;
    static final String[] ANNOTATIONS = {"BaseQualityRankSumTest", "ChromosomeCounts", "Coverage", "DepthPerAlleleBySample",
            "DepthPerSampleHC", "ExcessHet", "FisherStrand", "InbreedingCoeff", "MappingQualityRankSumTest", "QualByDepth",
            "RMSMappingQuality", "ReadPosRankSumTest", "StrandOddsRatio"};

    public static void main(final String[] args) throws Exception {
        System.out.println("# HcCallRegionDump: HaplotypeCallerEngine.callRegion");
        final Random random = new Random(20261016L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 1500; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        ref = sb.toString();
        System.out.println("reference\t" + ref);
        final Path dir = Files.createTempDirectory("hccall");
        fasta = dir.resolve("ref.fasta");
        try (final FastaReferenceWriter w = new FastaReferenceWriterBuilder().setFastaFile(fasta)
                .setMakeFaiOutput(true).setMakeDictOutput(true).build()) {
            w.addSequence(new ReferenceSequence("chr1", 0, ref.getBytes(StandardCharsets.US_ASCII)));
        }
        final List<String> one = List.of("s1");
        run("het-snp", 700, 800, 100, one, het("s1", snp(ref, 749), 0));
        run("hom-snp", 700, 800, 100, one, het("s1", snp(ref, 749), 100));
        run("deletion", 700, 800, 100, one, het("s1", ref.substring(0, 750) + ref.substring(754), 0));
        run("insertion", 700, 800, 100, one, het("s1", ref.substring(0, 750) + "TTG" + ref.substring(750), 0));
        run("two-snps-cis", 700, 800, 100, one, het("s1", snp(snp(ref, 729), 769), 0));
        run("no-variation", 700, 800, 100, one, het("s1", ref, 0));
        run("soft-clipped", 700, 800, 100, one, softClipped("s1", snp(ref, 749)));
        run("low-mapq-reads", 700, 800, 100, one, lowMapq(het("s1", snp(ref, 749), 0)));
        run("two-samples", 700, 800, 100, List.of("s1", "s2"),
                concat(het("s1", snp(ref, 749), 0), het("s2", snp(ref, 749), 100)));
        run("near-contig-start", 50, 150, 100, one, hetFrom("s1", snp(ref, 99), 0, 1, 230));
    }

    static List<Spec> concat(final List<Spec> a, final List<Spec> b) {
        final List<Spec> all = new ArrayList<>(a);
        all.addAll(b);
        return all;
    }

    static List<Spec> lowMapq(final List<Spec> in) {
        final List<Spec> out = new ArrayList<>();
        for (int i = 0; i < in.size(); i++) {
            final Spec s = in.get(i);
            out.add(new Spec(s.sample(), s.start(), s.cigar(), s.bases(), s.quals(), s.flags(), i % 4 == 1 ? 15 : 60));
        }
        return out;
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    static List<Spec> het(final String sample, final String alt, final int altHomPercent) {
        return hetFrom(sample, alt, altHomPercent, 600, 900);
    }

    /**
     * Reads of 60 bases every 4 from `from`, alternating reference and alternate haplotype (all
     * alternate when `altHomPercent` is 100), with qualities that vary and some low tails. The
     * alternate's reads are aligned by position only while it has the reference's length;
     * otherwise the read carries the indel in its CIGAR, cut at base 750.
     */
    static List<Spec> hetFrom(final String sample, final String alt, final int altHomPercent, final int from, final int to) {
        final List<Spec> out = new ArrayList<>();
        final int delta = alt.length() - ref.length();
        int i = 0;
        for (int start = from; start + 60 <= to + 1 && start + 60 <= ref.length(); start += 4) {
            final boolean useAlt = altHomPercent == 100 || i % 2 == 1;
            String bases;
            String cigar = "60M";
            final int crossing = 751 - start;
            final boolean indelInside = crossing >= 5 && crossing + Math.max(delta, 0) <= 55;
            if (useAlt && delta != 0 && start <= 751 && start + 60 > 751 && !indelInside) {
                // An alternate read whose indel would sit at an end: the reference read instead.
                bases = ref.substring(start - 1, start - 1 + 60);
            } else if (!useAlt || delta == 0 || start + 60 <= 751 || start > 751) {
                final String src = useAlt ? alt : ref;
                final int offset = (useAlt && start > 751) ? delta : 0;
                bases = src.substring(start - 1 + offset, start - 1 + offset + 60);
            } else {
                final int left = 751 - start;
                if (delta > 0) {
                    bases = alt.substring(start - 1, start - 1 + 60);
                    cigar = left + "M" + delta + "I" + (60 - left - delta) + "M";
                } else {
                    bases = alt.substring(start - 1, start - 1 + 60);
                    cigar = left + "M" + (-delta) + "D" + (60 - left) + "M";
                }
            }
            final StringBuilder q = new StringBuilder();
            for (int k = 0; k < 60; k++) {
                q.append((char) (33 + ((i % 5 == 0 && k >= 56) ? 8 : 25 + (k * 3 + i) % 15)));
            }
            out.add(new Spec(sample, start, cigar, bases, q.toString(), i % 3 == 0 ? 16 : 0, 60));
            i++;
        }
        return out;
    }

    static List<Spec> softClipped(final String sample, final String alt) {
        final List<Spec> out = new ArrayList<>(het(sample, alt, 0));
        for (int start = 735; start <= 745; start += 2) {
            final String bases = "GGGGG" + alt.substring(start - 1, start - 1 + 55);
            out.add(new Spec(sample, start, "5S55M", bases, "I".repeat(60), 0, 60));
        }
        return out;
    }

    static void run(final String label, final int activeStart, final int activeEnd, final int padding, final List<String> samples,
                    final List<Spec> unsorted) {
        System.out.println("case\t" + label + "\t" + activeStart + "\t" + activeEnd + "\t" + padding + "\t" + String.join(",", samples));
        try {
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", ref.length()))));
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            for (final String sample : samples) {
                final SAMReadGroupRecord g = new SAMReadGroupRecord("rg" + sample);
                g.setSample(sample);
                header.addReadGroup(g);
            }
            final List<Spec> specs = new ArrayList<>(unsorted);
            specs.sort(Comparator.comparingInt(Spec::start));
            final AssemblyRegion region = new AssemblyRegion(new SimpleInterval("chr1", activeStart, activeEnd), padding, header);
            int n = 0;
            for (final Spec s : specs) {
                final byte[] q = s.quals().getBytes(StandardCharsets.US_ASCII);
                for (int k = 0; k < q.length; k++) {
                    q[k] -= 33;
                }
                final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "r" + n, 0, s.start(),
                        s.bases().getBytes(StandardCharsets.US_ASCII), q, s.cigar());
                read.setReadGroup("rg" + s.sample());
                read.setMappingQuality(s.mapq());
                read.setIsReverseStrand((s.flags() & 16) != 0);
                if (!region.getPaddedSpan().overlaps(read)) {
                    continue;
                }
                System.out.println("read\t" + label + "\tr" + n + "\t" + s.sample() + "\t" + s.start() + "\t" + s.cigar() + "\t" + s.bases()
                        + "\t" + s.quals() + "\t" + s.flags() + "\t" + s.mapq());
                region.add(read);
                n++;
            }
            final List<Annotation> annotations = new ArrayList<>();
            for (final String name : ANNOTATIONS) {
                annotations.add((Annotation) Class.forName("org.broadinstitute.hellbender.tools.walkers.annotator." + name)
                        .getDeclaredConstructor().newInstance());
            }
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            final HaplotypeCallerEngine engine = new HaplotypeCallerEngine(hcArgs,
                    new org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection(), false, false, header,
                    new CachingIndexedFastaSequenceFile(fasta),
                    new VariantAnnotatorEngine(annotations, null, Collections.emptyList(), false, false));
            Utils.resetRandomGenerator();
            final List<VariantContext> calls = engine.callRegion(region, new FeatureContext(),
                    new ReferenceContext(new ReferenceFileSource(fasta), region.getPaddedSpan()));
            final VCFEncoder encoder = new VCFEncoder(new VCFHeader(Collections.emptySet(), samples), true, false);
            for (final VariantContext vc : calls) {
                System.out.println("call\t" + label + "\t" + encoder.encode(vc));
            }
            System.out.println("calls\t" + label + "\t" + calls.size());
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }
}
