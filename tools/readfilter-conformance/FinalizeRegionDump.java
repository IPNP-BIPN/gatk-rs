/*
 * AssemblyBasedCallerUtils.finalizeRegion: the reads of an assembly region as the assembler and the
 * genotyper see them.
 *
 * Each read loses its soft clips (hard-clipped, or reverted with `os`/`oe` when soft-clipped bases
 * are used and the fragment size is well defined), its low-quality ends (soft- or hard-clipped,
 * the threshold 6 under error correction), its adaptor, and what lies outside the padded span;
 * an empty or unaligned remainder is dropped. The rest is sorted by coordinate, and the overlapping
 * mates of a fragment have the qualities of their shared bases capped at 20, or zeroed where the
 * bases differ, per sample.
 *
 * Output:
 *
 *     contig\t<bases>
 *     read\t<name>\t<flags>\t<start>\t<cigar>\t<mate start>\t<fragment length>\t<sample>\t<bases>\t<quals>
 *     case\t<label>\t<errorCorrect>\t<dontUseSoftClipped>\t<minTail>\t<correctOverlap>\t<softClipLowQualEnds>\t<override>\t<track>\t<prefinalized>
 *     kept\t<label>\t<list>\t<index>\t<name>\t<start>\t<cigar>\t<bases>\t<quals>\t<os>\t<oe>
 *     region\t<label>\t<reads>\t<hard-clipped reads>\t<finalized>
 *     error\t<label>\t<exception class>: <message>
 *
 * `list` is `reads` or `hardclipped`. Qualities are Phred + 33.
 *
 * Usage: FinalizeRegionDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class FinalizeRegionDump {

    record Spec(String name, int flags, int start, String cigar, int mateStart, int tlen, String sample, String bases, byte[] quals) {
    }

    static String contig;
    static SAMFileHeader header;
    static final List<Spec> SPECS = new ArrayList<>();

    public static void main(final String[] args) {
        System.out.println("# FinalizeRegionDump: AssemblyBasedCallerUtils.finalizeRegion");
        final Random random = new Random(20261012L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 600; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        contig = sb.toString();
        System.out.println("contig\t" + contig);
        header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", contig.length()))));
        for (final String s : List.of("s1", "s2")) {
            final SAMReadGroupRecord g = new SAMReadGroupRecord("rg-" + s);
            g.setSample(s);
            header.addReadGroup(g);
        }

        // Flags: 1 paired, 2 proper, 16 reverse, 32 mate reverse, 64 first, 128 second.
        final int fwdFirst = 1 | 2 | 32 | 64;
        final int revSecond = 1 | 2 | 16 | 128;
        spec("r-outside-left", 0, 120, "50M", 0, 0, "s1", ref(120, 50), q(50, 30));
        spec("r-plain", 0, 160, "50M", 0, 0, "s1", ref(160, 50), q(50, 30));
        spec("r-softclip-left", 0, 170, "5S45M", 0, 0, "s1", "TTTTT" + ref(170, 45), q(50, 30));
        spec("r-softclip-right", 0, 180, "45M5S", 0, 0, "s1", ref(180, 45) + "GGGGG", q(50, 30));
        spec("r-lowqual-tail", 0, 185, "50M", 0, 0, "s1", ref(185, 50), tail(q(50, 30), 6, 5));
        spec("r-lowqual-head", 0, 190, "50M", 0, 0, "s1", ref(190, 50), head(q(50, 30), 4, 3));
        spec("r-lowqual-tail-8", 0, 192, "50M", 0, 0, "s1", ref(192, 50), tail(q(50, 30), 5, 8));
        spec("r-all-low", 0, 195, "50M", 0, 0, "s1", ref(195, 50), q(50, 4));
        spec("pairA", fwdFirst, 200, "50M", 230, 80, "s1", ref(200, 50), q(50, 35));
        spec("r-deletion", 0, 205, "20M5D30M", 0, 0, "s1", ref(205, 20) + ref(230, 30), q(50, 30));
        spec("pairB", fwdFirst, 210, "50M", 225, 65, "s1", ref(210, 50), q(50, 35));
        spec("pairS2", fwdFirst, 212, "50M", 222, 60, "s2", ref(212, 50), q(50, 33));
        spec("r-softclip-both", 0, 215, "4S40M6S", 0, 0, "s2", "AAAA" + ref(215, 40) + "CCCCCC", q(50, 30));
        spec("pairS2", revSecond, 222, "50M", 212, -60, "s2", ref(222, 50), q(50, 31));
        spec("pairB", revSecond, 225, "50M", 210, -65, "s1", mutate(ref(225, 50), 10), q(50, 32));
        spec("pairA", revSecond, 230, "50M", 200, -80, "s1", ref(230, 50), q(50, 25));
        spec("pairC", fwdFirst, 240, "50M", 245, 40, "s1", ref(240, 50), q(50, 30));
        spec("pairC", revSecond, 245, "50M", 240, -40, "s1", ref(245, 50), q(50, 30));
        spec("pairD", fwdFirst | 0, 270, "5S45M", 290, 70, "s1", "GGGGG" + ref(270, 45), q(50, 30));
        spec("pairD", revSecond, 290, "45M5S", 270, -70, "s1", ref(290, 45) + "TTTTT", q(50, 30));
        spec("r-softclip-lowqual", 0, 300, "6S44M", 0, 0, "s1", "ACGTAC" + ref(300, 44), head(q(50, 30), 8, 2));
        spec("r-near-end", 0, 330, "50M", 0, 0, "s2", ref(330, 50), q(50, 30));
        spec("r-after-padding", 0, 360, "50M", 0, 0, "s1", ref(360, 50), q(50, 30));
        for (final Spec s : SPECS) {
            System.out.println("read\t" + s.name + "\t" + build(s).convertToSAMRecord(header).getFlags() + "\t" + s.start + "\t" + s.cigar + "\t" + s.mateStart + "\t" + s.tlen
                    + "\t" + s.sample + "\t" + s.bases + "\t" + phred(s.quals));
        }

        run("hc", false, false, 10, true, false, false, false, false);
        run("dont-use-soft-clipped", false, true, 10, true, false, false, false, false);
        run("override-fragment-check", false, false, 10, true, false, true, false, false);
        run("error-correct", true, false, 10, true, false, false, false, false);
        run("min-tail-20", false, false, 20, true, false, false, false, false);
        run("soft-clip-low-qual-ends", false, false, 10, true, true, false, false, false);
        run("no-overlap-correction", false, false, 10, false, false, false, false, false);
        run("track-hardclipped", false, false, 10, true, false, false, true, false);
        run("track-hardclipped-soft-ends", false, false, 10, true, true, true, true, false);
        run("already-finalized", false, false, 10, true, false, false, false, true);
    }

    static String ref(final int start, final int length) {
        return contig.substring(start - 1, start - 1 + length);
    }

    static String mutate(final String s, final int at) {
        final char c = s.charAt(at);
        return s.substring(0, at) + "ACGT".charAt(("ACGT".indexOf(c) + 1) % 4) + s.substring(at + 1);
    }

    static byte[] q(final int n, final int value) {
        final byte[] b = new byte[n];
        java.util.Arrays.fill(b, (byte) value);
        return b;
    }

    static byte[] tail(final byte[] q, final int n, final int value) {
        for (int i = q.length - n; i < q.length; i++) {
            q[i] = (byte) value;
        }
        return q;
    }

    static byte[] head(final byte[] q, final int n, final int value) {
        for (int i = 0; i < n; i++) {
            q[i] = (byte) value;
        }
        return q;
    }

    static String phred(final byte[] q) {
        final StringBuilder b = new StringBuilder();
        for (final byte x : q) {
            b.append((char) (x + 33));
        }
        return b.toString();
    }

    static void spec(final String name, final int flags, final int start, final String cigar, final int mateStart, final int tlen,
                     final String sample, final String bases, final byte[] quals) {
        SPECS.add(new Spec(name, flags, start, cigar, mateStart, tlen, sample, bases, quals));
    }

    static GATKRead build(final Spec s) {
        final GATKRead r = ArtificialReadUtils.createArtificialRead(header, s.name, 0, s.start,
                s.bases.getBytes(StandardCharsets.US_ASCII), s.quals.clone(), s.cigar);
        r.setReadGroup("rg-" + s.sample);
        if ((s.flags & 1) != 0) {
            r.setIsPaired(true);
            r.setIsProperlyPaired((s.flags & 2) != 0);
            r.setIsReverseStrand((s.flags & 16) != 0);
            r.setMateIsReverseStrand((s.flags & 32) != 0);
            if ((s.flags & 128) != 0) {
                r.setIsSecondOfPair();
            } else {
                r.setIsFirstOfPair();
            }
            r.setMatePosition("chr1", s.mateStart);
            r.setFragmentLength(s.tlen);
        }
        return r;
    }

    static void run(final String label, final boolean errorCorrect, final boolean dontUseSoft, final int minTail, final boolean correctOverlap,
                    final boolean softClipLowQualEnds, final boolean override, final boolean track, final boolean prefinalized) {
        System.out.println("case\t" + label + "\t" + errorCorrect + "\t" + dontUseSoft + "\t" + minTail + "\t" + correctOverlap + "\t"
                + softClipLowQualEnds + "\t" + override + "\t" + track + "\t" + prefinalized);
        try {
            final AssemblyRegion region = new AssemblyRegion(new SimpleInterval("chr1", 200, 300), 50, header);
            for (final Spec s : SPECS) {
                final GATKRead r = build(s);
                if (region.getPaddedSpan().overlaps(r)) {
                    region.add(r);
                }
            }
            if (prefinalized) {
                region.setFinalized(true);
            }
            AssemblyBasedCallerUtils.finalizeRegion(region, errorCorrect, dontUseSoft, (byte) minTail, header,
                    new IndexedSampleList(List.of("s1", "s2")), correctOverlap, softClipLowQualEnds, override, track);
            print(label, "reads", region.getReads());
            print(label, "hardclipped", region.getHardClippedPileupReads());
            System.out.println("region\t" + label + "\t" + region.getReads().size() + "\t" + region.getHardClippedPileupReads().size() + "\t" + region.isFinalized());
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }

    static void print(final String label, final String list, final List<GATKRead> reads) {
        for (int i = 0; i < reads.size(); i++) {
            final GATKRead r = reads.get(i);
            System.out.println("kept\t" + label + "\t" + list + "\t" + i + "\t" + r.getName() + "\t" + r.getStart() + "\t" + r.getCigar()
                    + "\t" + r.getBasesString() + "\t" + phred(r.getBaseQualities())
                    + "\t" + r.getAttributeAsInteger("os") + "\t" + r.getAttributeAsInteger("oe"));
        }
    }
}
