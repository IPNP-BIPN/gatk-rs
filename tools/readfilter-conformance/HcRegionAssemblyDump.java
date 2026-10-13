/*
 * HaplotypeCaller's region assembly up to genotyping: AssemblyBasedCallerUtils.assembleReads (the
 * region finalized as HaplotypeCaller asks, the reference padded by 500, the reference haplotype,
 * the assembler of HaplotypeCaller's argument collection), AssemblyResultSet.getVariationEvents,
 * AssemblyRegionTrimmer.trim and AssemblyResultSet.trimTo (Haplotype.trim to the trimmed padded
 * span, duplicates collapsed, the reference rebuilt, sorted by length and bases).
 *
 * Output:
 *
 *     reference\t<chr1 bases>
 *     case\t<label>\t<active start>\t<active end>\t<padding>
 *     read\t<label>\t<name>\t<sample>\t<start>\t<cigar>\t<bases>\t<quals phred+33>\t<flags>\t<mate start>\t<tlen>
 *     finalized\t<label>\t<name>\t<start>\t<cigar>\t<bases>\t<quals>      (the region's reads after finalizeRegion)
 *     assembled\t<label>\t<index>\t<ref>\t<bases>\t<cigar>\t<alignment start>\t<location>
 *     event\t<label>\t<start>\t<ref>\t<alt>
 *     trim\t<label>\t<variant span>\t<padded span>
 *     trimmed\t<label>\t<index>\t<ref>\t<bases>\t<cigar>\t<alignment start>\t<location>\t<kmer>
 *     set\t<label>\t<haplotype count>\t<variation present>\t<genotyping span>\t<genotyping padded span>\t<reads>
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: HcRegionAssemblyDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
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

public class HcRegionAssemblyDump {

    record Spec(String sample, int start, String cigar, String bases, String quals, int flags, int mateStart, int tlen) {
    }

    static String ref;
    static Path fasta;

    public static void main(final String[] args) throws Exception {
        System.out.println("# HcRegionAssemblyDump: assembleReads, getVariationEvents, trim and trimTo");
        final Random random = new Random(20261016L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 1500; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        ref = sb.toString();
        System.out.println("reference\t" + ref);
        final Path dir = Files.createTempDirectory("hcasm");
        fasta = dir.resolve("ref.fasta");
        try (final FastaReferenceWriter w = new FastaReferenceWriterBuilder().setFastaFile(fasta)
                .setMakeFaiOutput(true).setMakeDictOutput(true).build()) {
            w.addSequence(new ReferenceSequence("chr1", 0, ref.getBytes(StandardCharsets.US_ASCII)));
        }

        run("het-snp", 700, 800, 100, het(snp(ref, 749), 0));
        run("hom-snp", 700, 800, 100, het(snp(ref, 749), 100));
        run("deletion", 700, 800, 100, het(ref.substring(0, 750) + ref.substring(754), 0));
        run("insertion", 700, 800, 100, het(ref.substring(0, 750) + "TTG" + ref.substring(750), 0));
        run("two-snps", 700, 800, 100, het(snp(snp(ref, 729), 769), 0));
        run("snp-outside-active", 700, 800, 100, het(snp(ref, 660), 0));
        run("no-variation", 700, 800, 100, het(ref, 0));
        run("near-contig-start", 50, 150, 100, hetFrom(snp(ref, 99), 0, 1, 230));
        run("near-contig-end", 1350, 1450, 100, hetFrom(snp(ref, 1399), 0, 1250, 1450));
        run("soft-clipped", 700, 800, 100, softClipped(snp(ref, 749)));
        run("small-padding", 700, 800, 25, het(snp(ref, 749), 0));
    }

    static String snp(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    static List<Spec> het(final String alt, final int altHomPercent) {
        return hetFrom(alt, altHomPercent, 600, 900);
    }

    /**
     * Reads of 60 bases every 4 from `from`, alternating reference and alternate haplotype (all
     * alternate when `altHomPercent` is 100), with qualities that vary and some low tails. The
     * alternate's reads are aligned by position only while it has the reference's length;
     * otherwise the read carries the indel in its CIGAR, cut at base 750.
     */
    static List<Spec> hetFrom(final String alt, final int altHomPercent, final int from, final int to) {
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
            out.add(new Spec("s1", start, cigar, bases, q.toString(), i % 3 == 0 ? 16 : 0, 0, 0));
            i++;
        }
        return out;
    }

    static List<Spec> softClipped(final String alt) {
        final List<Spec> out = new ArrayList<>(het(alt, 0));
        for (int start = 735; start <= 745; start += 2) {
            final String bases = "GGGGG" + alt.substring(start - 1, start - 1 + 55);
            out.add(new Spec("s1", start, "5S55M", bases, "I".repeat(60), 0, 0, 0));
        }
        return out;
    }

    static String show(final Haplotype h) {
        return h.isReference() + "\t" + h.getBaseString() + "\t" + h.getCigar() + "\t" + h.getAlignmentStartHapwrtRef() + "\t" + h.getGenomeLocation();
    }

    static void run(final String label, final int activeStart, final int activeEnd, final int padding, final List<Spec> unsorted) {
        System.out.println("case\t" + label + "\t" + activeStart + "\t" + activeEnd + "\t" + padding);
        try {
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", ref.length()))));
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            final SAMReadGroupRecord g = new SAMReadGroupRecord("rg1");
            g.setSample("s1");
            header.addReadGroup(g);
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
                read.setReadGroup("rg1");
                read.setMappingQuality(60);
                read.setIsReverseStrand((s.flags() & 16) != 0);
                if (!region.getPaddedSpan().overlaps(read)) {
                    continue;
                }
                System.out.println("read\t" + label + "\tr" + n + "\t" + s.sample() + "\t" + s.start() + "\t" + s.cigar() + "\t" + s.bases()
                        + "\t" + s.quals() + "\t" + s.flags() + "\t" + s.mateStart() + "\t" + s.tlen());
                region.add(read);
                n++;
            }
            final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
            final CachingIndexedFastaSequenceFile reader = new CachingIndexedFastaSequenceFile(fasta);
            final AssemblyResultSet untrimmed = AssemblyBasedCallerUtils.assembleReads(region, hcArgs, header,
                    new IndexedSampleList(List.of("s1")), LogManager.getLogger("dump"), reader, hcArgs.createReadThreadingAssembler(),
                    SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA), !hcArgs.doNotCorrectOverlappingBaseQualities,
                    hcArgs.fbargs, false);
            for (final GATKRead r : region.getReads()) {
                System.out.println("finalized\t" + label + "\t" + r.getName() + "\t" + r.getStart() + "\t" + r.getCigar() + "\t"
                        + r.getBasesString() + "\t" + qualString(r.getBaseQualities()));
            }
            final List<Haplotype> assembled = untrimmed.getHaplotypeList();
            for (int i = 0; i < assembled.size(); i++) {
                System.out.println("assembled\t" + label + "\t" + i + "\t" + show(assembled.get(i)));
            }
            final SortedSet<Event> events = untrimmed.getVariationEvents(0);
            for (final Event e : events) {
                System.out.println("event\t" + label + "\t" + e.getStart() + "\t" + e.refAllele().getBaseString() + "\t" + e.altAllele().getBaseString());
            }
            final AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer(new org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection(), header.getSequenceDictionary());
            final AssemblyRegionTrimmer.Result result = trimmer.trim(region, events,
                    new ReferenceContext(new ReferenceFileSource(fasta), region.getPaddedSpan()));
            if (!result.isVariationPresent()) {
                System.out.println("trim\t" + label + "\tnull\tnull");
                return;
            }
            final AssemblyRegion variantRegion = result.getVariantRegion();
            System.out.println("trim\t" + label + "\t" + variantRegion.getSpan() + "\t" + variantRegion.getPaddedSpan());
            final AssemblyResultSet trimmed = untrimmed.trimTo(variantRegion);
            final List<Haplotype> haps = trimmed.getHaplotypeList();
            for (int i = 0; i < haps.size(); i++) {
                System.out.println("trimmed\t" + label + "\t" + i + "\t" + show(haps.get(i)) + "\t" + haps.get(i).getKmerSize());
            }
            final AssemblyRegion forGenotyping = trimmed.getRegionForGenotyping();
            System.out.println("set\t" + label + "\t" + trimmed.getHaplotypeCount() + "\t" + trimmed.isVariationPresent() + "\t"
                    + forGenotyping.getSpan() + "\t" + forGenotyping.getPaddedSpan() + "\t"
                    + forGenotyping.getReads().stream().map(GATKRead::getName).collect(Collectors.joining(",")));
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }

    static String qualString(final byte[] q) {
        final StringBuilder b = new StringBuilder();
        for (final byte x : q) {
            b.append((char) (x + 33));
        }
        return b.toString();
    }
}
