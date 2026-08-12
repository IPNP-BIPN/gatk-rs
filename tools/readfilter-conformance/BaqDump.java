/*
 * BAQ, taken from the reference.
 *
 * Base Alignment Quality: a hidden Markov model that caps each base's quality by how confidently the
 * aligner could have placed it. BaseRecalibrator runs it over every read before counting, so the
 * recalibration table depends on these bytes and they are decided by a forward-backward pass in
 * doubles.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE EMISSION TABLE IS 256 x 256 x 94 AND MOST OF IT IS 1.0. Only the sixteen (ACGT x acgt)
 *     pairs are ever filled, so an N against an A has an emission probability of ONE, not of a
 *     mismatch. That is not a special case in the algorithm: it is what the uninitialised table
 *     holds;
 *   - AND ITS QUALITY IS FLOORED, NOT THE BASE'S. `qual2prob[q < minBaseQual ? minBaseQual : q]`
 *     raises a quality below the minimum for the emission only, leaving the read's own byte alone;
 *   - THE BAND IS A DIAGONAL WINDOW AND THE INDEXING IS BIT-PACKED. `set_u(b, i, k)` is
 *     `(k + 1 - max(i - b, 0)) * 3`, three states per cell, and the state array carries the aligned
 *     position shifted left by two with the indel flag in the low bits;
 *   - THE FORWARD PASS RESCALES BY A RUNNING SUM and the backward pass does not, so the two are not
 *     symmetric and the posterior is a ratio of differently scaled numbers;
 *   - THE REFERENCE WINDOW IS NOT THE READ'S SPAN. It is the alignment start minus half the band
 *     width minus the FIRST insertion's offset, to the end plus half the band width plus the LAST
 *     insertion's offset, clamped at one;
 *   - A READ WHOSE WINDOW RUNS PAST THE CONTIG GETS NO BAQ AT ALL, a null rather than a shorter
 *     window;
 *   - AN INSERTION OR A SOFT CLIP KEEPS ITS RAW QUALITY, because the capping loop copies rawQuals
 *     over bq for those elements, and the `case S:` falls through into `case I:` after moving the
 *     reference, which the reference's own comment questions;
 *   - AND THE BQ TAG IS THE DIFFERENCE, NOT THE VALUE: encodeBQTag writes
 *     `(char)(quality - baq + 64)`, so a read whose BAQ equals its quality carries a tag of all `@`.
 *
 * Output:
 *
 *     const\t<name>\t<value>
 *     qual2prob\t<q>\t<bits>\t<decimal>
 *     epsilon\t<ref>\t<read>\t<qual>\t<bits>\t<decimal>
 *     hmm\t<label>\t<return>\t<comma separated bq>\t<comma separated state>
 *     window\t<label>\t<contig>:<start>-<stop>
 *     baq\t<label>\t<comma separated bq or null>
 *     tag\t<label>\t<the BQ tag>
 *     fromtag\t<label>\t<comma separated quals>
 *     error\t<what>\t<exception>\t<message>
 *
 * Usage: BaqDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.baq.BAQ;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;

public class BaqDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# BaqDump: BAQ, the hidden Markov model BaseRecalibrator caps qualities with");

        constants();
        emissions();
        hmm();
        reads();
        tags();
    }

    static void constants() {
        final BAQ baq = new BAQ();
        System.out.printf("const\tDEFAULT_GOP\t%s%n", BAQ.DEFAULT_GOP);
        System.out.printf("const\tDEFAULT_BANDWIDTH\t%d%n", BAQ.DEFAULT_BANDWIDTH);
        System.out.printf("const\tBAQ_TAG\t%s%n", BAQ.BAQ_TAG);
        System.out.printf("const\tminBaseQual\t%d%n", baq.getMinBaseQual());
        System.out.printf("const\tgapOpenProb\t%s%n", baq.getGapOpenProb());
        System.out.printf("const\tgapExtensionProb\t%s%n", baq.getGapExtensionProb());
        System.out.printf("const\tbandWidth\t%d%n", baq.getBandWidth());
    }

    /**
     * The quality-to-probability cache and a slice of the emission table, including the pairs the
     * table never fills.
     */
    static void emissions() throws Exception {
        final Field field = BAQ.class.getDeclaredField("qual2prob");
        field.setAccessible(true);
        final double[] qual2prob = (double[]) field.get(null);
        for (final int q : new int[] {0, 1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 93, 255}) {
            System.out.printf("qual2prob\t%d\t%s\t%s%n", q,
                    Long.toHexString(Double.doubleToRawLongBits(qual2prob[q])), qual2prob[q]);
        }

        final BAQ baq = new BAQ();
        final java.lang.reflect.Method epsilon =
                BAQ.class.getDeclaredMethod("calcEpsilon", byte.class, byte.class, byte.class);
        epsilon.setAccessible(true);
        // Matching and mismatching bases, in both cases, at qualities on both sides of the floor,
        // plus the pairs the table leaves at one.
        for (final char ref : new char[] {'A', 'a', 'C', 'N', '-'}) {
            for (final char read : new char[] {'A', 'a', 'C', 'N'}) {
                for (final int q : new int[] {0, 3, 4, 5, 20, 93}) {
                    final double value = (double) epsilon.invoke(baq, (byte) ref, (byte) read,
                            (byte) q);
                    System.out.printf("epsilon\t%c\t%c\t%d\t%s\t%s%n", ref, read, q,
                            Long.toHexString(Double.doubleToRawLongBits(value)), value);
                }
            }
        }
    }

    /** The model itself, over sequences shaped to exercise each of its branches. */
    static void hmm() {
        final BAQ baq = new BAQ();
        // A reference with no internal repeat, so a matching read has ONE place to sit and the
        // model can be confident. A repeating one cannot tell the placements apart and floors every
        // base at minBaseQual, which the `repeat` case below shows on purpose.
        final String unique = "GATTACAGGCTCTAGCAT";
        final String[][] cases = {
                // An exact match against a unique reference, which is what confidence looks like.
                {"exact", unique, "TTACAGGC", "IIIIIIII"},
                // One mismatch in the middle of the same placement.
                {"mismatch", unique, "TTACTGGC", "IIIIIIII"},
                // A read that could sit anywhere, which is what BAQ exists to catch.
                {"repeat", "ACACACACACACACACACAC", "ACACAC", "IIIIII"},
                // Low qualities, so the floor at minBaseQual matters.
                {"low-quality", unique, "TTACAGGC", "!!!!!!!!"},
                // Qualities on both sides of the floor in one read.
                {"mixed-quality", unique, "TTACAGGC", "!\"#$%&\'("},
                // An N in the read, whose emission probability is the uninitialised one.
                {"n-in-read", unique, "TTANAGGC", "IIIIIIII"},
                // An N in the reference.
                {"n-in-ref", "GATTACANGCTCTAGCAT", "TTACAGGC", "IIIIIIII"},
                // A one-base query, which is the shortest the model accepts.
                {"one-base", unique, "A", "I"},
                // A query as long as the reference, so the band reaches both ends.
                {"full-length", "ACGTACGTAC", "ACGTACGTAC", "IIIIIIIIII"},
                // A query longer than the reference, which the band cannot cover.
                {"query-longer-than-ref", "ACGT", "ACGTACGT", "IIIIIIII"},
        };
        for (final String[] one : cases) {
            final byte[] ref = one[1].getBytes(StandardCharsets.UTF_8);
            final byte[] query = one[2].getBytes(StandardCharsets.UTF_8);
            final byte[] quals = one[3].getBytes(StandardCharsets.UTF_8);
            for (int i = 0; i < quals.length; i++) {
                quals[i] -= 33;
            }
            final int[] state = new int[query.length];
            final byte[] bq = new byte[query.length];
            final int returned = baq.hmm_glocal(ref, query, 0, query.length, quals, state, bq);
            System.out.printf("hmm\t%s\t%d\t%s\t%s%n", one[0], returned, bytes(bq), ints(state));
        }

        // The two state decoders, over the values the array carries.
        for (final int value : new int[] {0, 1, 2, 3, 4, 5, 8, 13, 400, 401}) {
            System.out.printf("state\t%d\t%b\t%d%n", value, BAQ.stateIsIndel(value),
                    BAQ.stateAlignedPosition(value));
        }
    }

    /** The reference window and the per-read calculation, over cigars that move it. */
    static void reads() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                java.util.List.of(new SAMSequenceRecord("chr1", 100))));

        final String[][] cases = {
                {"aligned", "10M", "50"},
                {"leading-insertion", "3I7M", "50"},
                {"trailing-insertion", "7M3I", "50"},
                {"both-insertions", "2I6M2I", "50"},
                {"soft-clipped", "3S7M", "50"},
                {"deletion", "4M2D6M", "50"},
                // At the very start, where the window clamps at one.
                {"at-contig-start", "10M", "1"},
                // At the very end, where the window runs past the contig and BAQ answers null.
                {"at-contig-end", "10M", "95"},
                // A cigar with an N, which BAQ refuses outright.
                {"skipped-region", "4M2N4M", "50"},
                // Entirely soft-clipped, which is the "completely clipped away" branch.
                {"all-soft-clipped", "10S", "50"},
        };
        final String reference = repeat("ACGT", 25);
        for (final String[] one : cases) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(one[0]);
            record.setReferenceName("chr1");
            record.setAlignmentStart(Integer.parseInt(one[2]));
            record.setCigarString(one[1]);
            final int length = 10;
            record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
            final byte[] quals = new byte[length];
            java.util.Arrays.fill(quals, (byte) 40);
            record.setBaseQualities(quals);
            record.setMappingQuality(60);
            final GATKRead read = new SAMRecordToGATKReadAdapter(record);

            try {
                System.out.printf("window\t%s\t%s%n", one[0],
                        BAQ.getReferenceWindowForRead(read, BAQ.DEFAULT_BANDWIDTH));
            } catch (final Exception e) {
                System.out.printf("error\twindow@%s\t%s\t%s%n", one[0],
                        e.getClass().getSimpleName(), e.getMessage());
            }

            // The offset form, which takes the reference bases directly rather than a data source.
            try {
                final BAQ.BAQCalculationResult result = new BAQ().calcBAQFromHMM(read,
                        reference.getBytes(StandardCharsets.UTF_8), 0);
                System.out.printf("baq\t%s\t%s%n", one[0],
                        result == null ? "null" : bytes(result.bq));
            } catch (final Exception e) {
                System.out.printf("error\tbaq@%s\t%s\t%s%n", one[0], e.getClass().getSimpleName(),
                        e.getMessage());
            }
        }
    }

    /** The BQ tag, which carries the difference between the quality and the BAQ. */
    static void tags() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                java.util.List.of(new SAMSequenceRecord("chr1", 100))));

        final byte[][] cases = {
                // Equal to the qualities, so every character is the zero of the encoding.
                {40, 40, 40, 40},
                // One below, one far below.
                {40, 39, 30, 0},
                // Above the quality, which the encoding does not forbid.
                {41, 45, 40, 40},
        };
        int index = 0;
        for (final byte[] baq : cases) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName("tag" + index);
            record.setReferenceName("chr1");
            record.setAlignmentStart(10);
            record.setCigarString(baq.length + "M");
            record.setReadBases("ACGTACGTAC".substring(0, baq.length)
                    .getBytes(StandardCharsets.UTF_8));
            final byte[] quals = new byte[baq.length];
            java.util.Arrays.fill(quals, (byte) 40);
            record.setBaseQualities(quals);
            final GATKRead read = new SAMRecordToGATKReadAdapter(record);

            System.out.printf("tag\t%d\t%s%n", index, BAQ.encodeBQTag(read, baq));
            BAQ.addBAQTag(read, baq);
            System.out.printf("hastag\t%d\t%b%n", index, BAQ.hasBAQTag(read));
            System.out.printf("fromtag\t%d\t%s%n", index,
                    bytes(BAQ.calcBAQFromTag(read, false, false)));
            // Overwriting the read's own qualities, which is what the caller usually asks for.
            System.out.printf("fromtag\t%d-overwrite\t%s%n", index,
                    bytes(BAQ.calcBAQFromTag(read, true, false)));
            index++;
        }

        // A read with no tag at all, asked for one both ways.
        final SAMRecord bare = new SAMRecord(header);
        bare.setReadName("bare");
        bare.setReferenceName("chr1");
        bare.setAlignmentStart(10);
        bare.setCigarString("4M");
        bare.setReadBases("ACGT".getBytes(StandardCharsets.UTF_8));
        bare.setBaseQualities(new byte[] {40, 40, 40, 40});
        final GATKRead read = new SAMRecordToGATKReadAdapter(bare);
        System.out.printf("hastag\tbare\t%b%n", BAQ.hasBAQTag(read));
        System.out.printf("fromtag\tbare-raw\t%s%n", bytes(BAQ.calcBAQFromTag(read, false, true)));
        try {
            System.out.printf("fromtag\tbare-strict\t%s%n",
                    bytes(BAQ.calcBAQFromTag(read, false, false)));
        } catch (final Exception e) {
            System.out.printf("error\tbare-strict\t%s\t%s%n", e.getClass().getSimpleName(),
                    e.getMessage());
        }
    }

    static String repeat(final String unit, final int times) {
        final StringBuilder out = new StringBuilder();
        for (int i = 0; i < times; i++) {
            out.append(unit);
        }
        return out.toString();
    }

    static String bytes(final byte[] values) {
        if (values == null) {
            return "null";
        }
        final StringBuilder out = new StringBuilder();
        for (final byte value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }

    static String ints(final int[] values) {
        final StringBuilder out = new StringBuilder();
        for (final int value : values) {
            if (out.length() != 0) {
                out.append(',');
            }
            out.append(value);
        }
        return out.toString();
    }
}
