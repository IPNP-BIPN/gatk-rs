/*
 * ReadUtils' reference-to-read coordinate mapping, taken from the reference.
 *
 * This is the arithmetic every walker, every annotation and every clipping operation stands on: a
 * reference position in, a read offset out. An off-by-one here does not fail loudly, it reads the
 * neighbouring base, and every number computed from that base is wrong by a plausible amount.
 *
 * Output, over the same corpus as the decision matrix:
 *
 *     soft\t<index>\t<softStart>\t<softEnd>\t<lastInsertionOffset>\t<basesReverseComplement>
 *     coord\t<index>\t<refCoord>\t<readIndex>\t<operator>\t<base>\t<quality>\t<inside>
 *
 * A field is `E` where the reference throws rather than answering, and `.` where it answers
 * "absent" (an empty Optional, or a null operator). Those are different outcomes: the first stops
 * a tool, the second is a decision the caller acts on.
 *
 * Usage: ReadCoordinateDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.List;
import java.util.Optional;

public class ReadCoordinateDump {

    /** How far outside a read's own span the probe walks, so the not-found answers are covered. */
    static final int MARGIN = 3;

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = ReadFilterDump.header();
        final List<SAMRecord> corpus = ReadFilterDump.corpus(header);

        System.out.println("# ReadCoordinateDump: ReadUtils coordinate mapping, from the reference");
        ReadFilterDump.printCorpus(header, corpus);

        for (int i = 0; i < corpus.size(); i++) {
            final GATKRead read = new SAMRecordToGATKReadAdapter(corpus.get(i));

            System.out.printf("soft\t%d\t%s\t%s\t%s\t%s%n",
                    i,
                    call(() -> String.valueOf(ReadUtils.getSoftStart(read))),
                    call(() -> String.valueOf(ReadUtils.getSoftEnd(read))),
                    call(() -> String.valueOf(ReadUtils.getLastInsertionOffset(read))),
                    call(() -> ReadUtils.getBasesReverseComplement(read)));

            final int from = read.getStart() - MARGIN;
            final int to = read.getEnd() + MARGIN;
            for (int refCoord = from; refCoord <= to; refCoord++) {
                final int coord = refCoord;
                System.out.printf("coord\t%d\t%d\t%s\t%s\t%s\t%s\t%s%n",
                        i, coord,
                        call(() -> String.valueOf(
                                ReadUtils.getReadIndexForReferenceCoordinate(read, coord).getLeft())),
                        call(() -> {
                            final Object op =
                                ReadUtils.getReadIndexForReferenceCoordinate(read, coord).getRight();
                            return op == null ? "." : op.toString();
                        }),
                        call(() -> optionalByte(
                                ReadUtils.getReadBaseAtReferenceCoordinate(read, coord))),
                        call(() -> optionalByte(
                                ReadUtils.getReadBaseQualityAtReferenceCoordinate(read, coord))),
                        call(() -> ReadUtils.isInsideRead(read, coord) ? "1" : "0"));
            }
        }
    }

    static String optionalByte(final Optional<Byte> value) {
        return value.map(b -> String.valueOf((int) (byte) b)).orElse(".");
    }

    /**
     * The value, or `E` where the reference throws.
     *
     * Several of these do throw on this corpus: a read with no qualities indexes an empty array,
     * and a read with no cigar elements asks for the last one. Recording that as its own outcome
     * keeps a crash from being ported as a plausible number.
     */
    static String call(final java.util.function.Supplier<String> supplier) {
        try {
            return supplier.get();
        } catch (final Exception | AssertionError e) {
            return "E";
        }
    }
}
