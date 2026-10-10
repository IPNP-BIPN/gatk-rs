/*
 * HaplotypeCallerEngine.isActive: the probability that a locus holds a variant, which the activity
 * profile turns into assembly regions.
 *
 * Every locus `LocusIteratorByState` yields over a case's reads (deletions and Ns kept, no
 * downsampling, as `AssemblyRegionIterator` builds it) goes through the engine's own `isActive`:
 * each sample's pileup becomes reference-versus-any genotype likelihoods
 * (`ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny`, bases at or under the minimum
 * quality skipped, a deletion scored at `--reference-model-deletion-quality`), stored as PLs;
 * one sample answers `calculateSingleSampleRefVsAnyActiveStateProfileValue`, several the QUAL of
 * the active-region genotyping engine (calling confidence 4) as a probability. Alternate bases
 * next to a soft clip average the read's high-quality soft-clipped bases, and an average above 6
 * marks the state HIGH_QUALITY_SOFT_CLIPS.
 *
 * Output:
 *
 *     reference\t<bases of chr1>
 *     case\t<label>\t<samples,...>\t<ploidy>\t<heterozygosity bits>\t<min base quality>\t<ref model deletion quality>
 *     read\t<label>\t<name>\t<sample>\t<start>\t<cigar>\t<bases>\t<quals phred+33>
 *     active\t<label>\t<pos>\t<depth>\t<prob bits>\t<type>\t<result value bits>\t<original active prob bits>
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: HcIsActiveDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class HcIsActiveDump {

    record Spec(String sample, int start, String cigar, String bases, String quals) {
    }

    static String ref;
    static Path fasta;

    public static void main(final String[] args) throws Exception {
        System.out.println("# HcIsActiveDump: HaplotypeCallerEngine.isActive");
        final Random random = new Random(20261015L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 300; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        ref = sb.toString();
        System.out.println("reference\t" + ref);
        final Path dir = Files.createTempDirectory("isactive");
        fasta = dir.resolve("ref.fasta");
        try (final FastaReferenceWriter w = new FastaReferenceWriterBuilder().setFastaFile(fasta)
                .setMakeFaiOutput(true).setMakeDictOutput(true).build()) {
            w.addSequence(new ReferenceSequence("chr1", 0, ref.getBytes(StandardCharsets.US_ASCII)));
        }

        final List<String> one = List.of("sA");
        final List<String> two = List.of("sA", "sB");
        run("het-snp", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'x', 20));
        run("hom-snp", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'x', 0));
        run("rare-snp", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'x', 80));
        run("low-quality-alt", one, 2, 0.001, 10, 30, het("sA", 100, 10, 'x', 20));
        run("quality-at-threshold", one, 2, 0.001, 10, 30, het("sA", 100, 11, 'x', 20));
        run("deletion", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'd', 20));
        run("insertion", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'i', 20));
        run("soft-clips", one, 2, 0.001, 10, 30, softClipped("sA", 100));
        run("ns", one, 2, 0.001, 10, 30, het("sA", 100, 0, 'n', 20));
        run("two-samples", two, 2, 0.001, 10, 30, concat(het("sA", 100, 0, 'x', 20), het("sB", 102, 0, 'x', 0)));
        run("two-samples-one-covered", two, 2, 0.001, 10, 30, het("sA", 100, 0, 'x', 20));
        run("two-samples-ref", two, 2, 0.001, 10, 30, concat(het("sA", 100, 0, 'x', 100), het("sB", 102, 0, 'x', 100)));
        run("ploidy-1", one, 1, 0.001, 10, 30, het("sA", 100, 0, 'x', 20));
        run("ploidy-3", one, 3, 0.001, 10, 30, het("sA", 100, 0, 'x', 20));
        run("ploidy-3-two-samples", two, 3, 0.001, 10, 30, concat(het("sA", 100, 0, 'x', 30), het("sB", 102, 0, 'x', 60)));
        run("heterozygosity", one, 2, 0.01, 10, 30, het("sA", 100, 0, 'x', 80));
        run("min-base-quality-20", one, 2, 0.001, 20, 30, het("sA", 100, 15, 'x', 20));
        run("deletion-quality-10", one, 2, 0.001, 10, 10, het("sA", 100, 0, 'd', 20));
    }

    static List<Spec> concat(final List<Spec> a, final List<Spec> b) {
        final List<Spec> all = new ArrayList<>(a);
        all.addAll(b);
        return all;
    }

    /**
     * Reads of 40 bases every 3 from 60 to 135, both haplotypes alternating; the alternate one
     * carries at 150 a substitution (`x`), a 3-base deletion (`d`), a 2-base insertion (`i`) or an
     * N (`n`). `altPercentStep` makes every read whose index modulo 100 is at or above it a
     * reference read (0: every read is alternate). `altQual` > 0 gives the alternate base that
     * quality.
     */
    static List<Spec> het(final String sample, final int salt, final int altQual, final char kind, final int refFraction) {
        final List<Spec> out = new ArrayList<>();
        int n = salt;
        for (int start = 60; start <= 135; start += 3) {
            for (int copy = 0; copy < 2; copy++) {
                n += 37;
                final boolean alt = (n % 100) >= refFraction;
                String bases = ref.substring(start - 1, start - 1 + 40);
                String cigar = "40M";
                final char[] q = "I".repeat(40).toCharArray();
                final int at = 150 - start;
                if (alt && at >= 1 && at < 36) {
                    final char c = bases.charAt(at);
                    switch (kind) {
                        case 'x' -> {
                            bases = bases.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + bases.substring(at + 1);
                            if (altQual > 0) {
                                q[at] = (char) (altQual + 33);
                            }
                        }
                        case 'n' -> bases = bases.substring(0, at) + 'N' + bases.substring(at + 1);
                        case 'd' -> {
                            bases = ref.substring(start - 1, start - 1 + at + 1) + ref.substring(start - 1 + at + 4, start - 1 + 43);
                            cigar = (at + 1) + "M3D" + (40 - at - 1) + "M";
                        }
                        case 'i' -> {
                            bases = ref.substring(start - 1, start - 1 + at + 1) + "TT" + ref.substring(start - 1 + at + 1, start - 1 + 38);
                            cigar = (at + 1) + "M2I" + (40 - at - 3) + "M";
                        }
                        default -> throw new IllegalArgumentException();
                    }
                }
                out.add(new Spec(sample, start, cigar, bases, new String(q)));
            }
        }
        return out;
    }

    /** Reads soft-clipped after an alternate base: high-quality clips on some, low on others. */
    static List<Spec> softClipped(final String sample, final int salt) {
        final List<Spec> out = new ArrayList<>();
        int n = salt;
        for (int start = 112; start <= 150; start += 3) {
            n += 37;
            final int aligned = 150 - start + 1;
            String bases = ref.substring(start - 1, start - 1 + aligned);
            final char c = bases.charAt(aligned - 1);
            bases = bases.substring(0, aligned - 1) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4);
            final String clip = "GATTACAGATTACA".substring(0, 10);
            final String q = "I".repeat(aligned) + ((n % 3 == 0) ? "#".repeat(10) : "I".repeat(10));
            out.add(new Spec(sample, start, aligned + "M10S", bases + clip, q));
            out.add(new Spec(sample, start, "40M", ref.substring(start - 1, start - 1 + 40), "I".repeat(40)));
        }
        return out;
    }

    static void run(final String label, final List<String> samples, final int ploidy, final double heterozygosity,
                    final int minBaseQuality, final int refModelDelQual, final List<Spec> specs) throws Exception {
        System.out.println("case\t" + label + "\t" + String.join(",", samples) + "\t" + ploidy + "\t"
                + Long.toHexString(Double.doubleToRawLongBits(heterozygosity)) + "\t" + minBaseQuality + "\t" + refModelDelQual);
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", ref.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String s : samples) {
            final SAMReadGroupRecord g = new SAMReadGroupRecord("rg" + s);
            g.setSample(s);
            header.addReadGroup(g);
        }
        final List<GATKRead> reads = new ArrayList<>();
        int i = 0;
        final List<Spec> sorted = new ArrayList<>(specs);
        sorted.sort(Comparator.comparingInt(Spec::start));
        for (final Spec s : sorted) {
            final String name = "r" + i++;
            final byte[] q = s.quals().getBytes(StandardCharsets.US_ASCII);
            for (int j = 0; j < q.length; j++) {
                q[j] -= 33;
            }
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, s.start(),
                    s.bases().getBytes(StandardCharsets.US_ASCII), q, s.cigar());
            read.setReadGroup("rg" + s.sample());
            read.setMappingQuality(60);
            reads.add(read);
            System.out.println("read\t" + label + "\t" + name + "\t" + s.sample() + "\t" + s.start() + "\t" + s.cigar() + "\t" + s.bases() + "\t" + s.quals());
        }
        try {
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            hcArgs.standardArgs.genotypeArgs.samplePloidy = ploidy;
            hcArgs.standardArgs.genotypeArgs.snpHeterozygosity = heterozygosity;
            hcArgs.minBaseQualityScore = (byte) minBaseQuality;
            hcArgs.refModelDelQual = (byte) refModelDelQual;
            final HaplotypeCallerEngine engine = new HaplotypeCallerEngine(hcArgs, new AssemblyRegionArgumentCollection(), false, false,
                    header, new CachingIndexedFastaSequenceFile(fasta),
                    new VariantAnnotatorEngine(Collections.emptyList(), null, Collections.emptyList(), false, false));
            final ReferenceFileSource source = new ReferenceFileSource(fasta);
            final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(), DownsamplingMethod.NONE,
                    ReadUtils.getSamplesFromHeader(header), header, true);
            while (libs.hasNext()) {
                final AlignmentContext context = libs.next();
                final SimpleInterval loc = new SimpleInterval(context.getContig(), (int) context.getPosition(), (int) context.getPosition());
                final ActivityProfileState state = engine.isActive(context, new ReferenceContext(source, loc), new FeatureContext());
                System.out.println("active\t" + label + "\t" + context.getPosition() + "\t" + context.getBasePileup().size()
                        + "\t" + Long.toHexString(Double.doubleToRawLongBits(state.isActiveProb()))
                        + "\t" + state.getResultState()
                        + "\t" + (state.getResultValue() == null ? "null" : Long.toHexString(Double.doubleToRawLongBits(state.getResultValue().doubleValue())))
                        + "\t" + Long.toHexString(Double.doubleToRawLongBits(state.getOriginalActiveProb())));
            }
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }
}
