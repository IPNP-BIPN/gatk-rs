/*
 * ReadClipper's output, taken from the reference.
 *
 * Clipping rewrites the read: its bases, its qualities, its cigar, its alignment start and, for a
 * flow read group, some of its tags. This dump prints the whole clipped read for every entry point
 * over the same corpus as the read filters, so the port is compared against the read the reference
 * produced rather than against one field of it.
 *
 * Output, one row per (read, operation):
 *
 *     clipped\t<index>\t<operation>\t<start>|<flags>|<cigar>|<bases>|<quals>|<mapq>|<tags>
 *
 * `E` where the reference throws: clipping the middle of a read, and soft-clipping an unmapped
 * one, are exceptions rather than reads, and a port that returned something there would emit a
 * read where the reference stops.
 *
 * Usage: ReadClipperDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.utils.clipping.ClippingRepresentation;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

public class ReadClipperDump {

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = ReadFilterDump.header();
        final List<SAMRecord> corpus = ReadFilterDump.corpus(header);

        System.out.println("# ReadClipperDump: ReadClipper's output, from the reference");
        ReadFilterDump.printCorpus(header, corpus);

        for (int i = 0; i < corpus.size(); i++) {
            final int index = i;
            final SAMRecord record = corpus.get(i);
            final GATKRead read = new SAMRecordToGATKReadAdapter(record);
            final int start = read.getStart();
            final int end = read.getEnd();

            // Reference coordinates chosen around the read's own span, so each entry point gets a
            // clip that lands inside it, one that lands on its edge and one that misses entirely.
            final List<Integer> coordinates = new ArrayList<>();
            for (final int offset : new int[] {-2, 0, 1, 3, 5}) {
                coordinates.add(start + offset);
                coordinates.add(end - offset);
            }

            for (final int coord : coordinates) {
                emit(index, "leftTail@" + coord, () ->
                        ReadClipper.hardClipByReferenceCoordinatesLeftTail(copy(record), coord));
                emit(index, "rightTail@" + coord, () ->
                        ReadClipper.hardClipByReferenceCoordinatesRightTail(copy(record), coord));
            }

            emit(index, "toRegion@" + (start + 1) + "," + (end - 1), () ->
                    ReadClipper.hardClipToRegion(copy(record), start + 1, end - 1));
            emit(index, "toRegion@" + (start - 5) + "," + (start - 2), () ->
                    ReadClipper.hardClipToRegion(copy(record), start - 5, start - 2));
            emit(index, "revertSoftClipped", () ->
                    ReadClipper.revertSoftClippedBases(copy(record)));

            // The soft-clip family. Its coordinates are the *unclipped* span, because
            // softClipToRegionIncludingClippedBases tests the read's reach with its clipped bases
            // counted, so a read can be inside a region its aligned part misses.
            emit(index, "softToRegion@" + (start + 1) + "," + (end - 1), () ->
                    ReadClipper.softClipToRegionIncludingClippedBases(
                            copy(record), start + 1, end - 1));
            emit(index, "softToRegion@" + (start - 5) + "," + (start - 2), () ->
                    ReadClipper.softClipToRegionIncludingClippedBases(
                            copy(record), start - 5, start - 2));
            emit(index, "softToRegion@" + start + "," + end, () ->
                    ReadClipper.softClipToRegionIncludingClippedBases(copy(record), start, end));
            emit(index, "softBothEnds@" + start + "," + end, () ->
                    ReadClipper.softClipBothEndsByReferenceCoordinates(copy(record), start, end));
            emit(index, "softBothEnds@" + (start + 1) + "," + (end - 1), () ->
                    ReadClipper.softClipBothEndsByReferenceCoordinates(
                            copy(record), start + 1, end - 1));
            // left == right: clipping both ends at one coordinate clips everything.
            emit(index, "softBothEnds@" + start + "," + start, () ->
                    ReadClipper.softClipBothEndsByReferenceCoordinates(copy(record), start, start));

            // Read coordinates, not reference ones: the whole read, one end, and the middle, which
            // is refused.
            final int length = read.getLength();
            emit(index, "softByRead@0," + (length - 1), () ->
                    ReadClipper.softClipByReadCoordinates(copy(record), 0, length - 1));
            emit(index, "softByRead@0,2", () ->
                    ReadClipper.softClipByReadCoordinates(copy(record), 0, 2));
            emit(index, "softByRead@" + (length - 3) + "," + (length - 1), () ->
                    ReadClipper.softClipByReadCoordinates(copy(record), length - 3, length - 1));
            emit(index, "softByRead@1,2", () ->
                    ReadClipper.softClipByReadCoordinates(copy(record), 1, 2));

            // Hard-clipping the soft clips away, with and without extra aligned bases.
            for (final int extra : new int[] {0, 1, 3}) {
                emit(index, "hardClipSoftClipped@" + extra, () ->
                        ReadClipper.hardClipSoftClippedBases(copy(record), extra));
            }
            emit(index, "hardClipAdaptor", () ->
                    ReadClipper.hardClipAdaptorSequence(copy(record)));
            for (final byte lowQual : new byte[] {0, 20, 30}) {
                emit(index, "hardClipLowQual@" + lowQual, () ->
                        ReadClipper.hardClipLowQualEnds(copy(record), lowQual));
                emit(index, "softClipLowQual@" + lowQual, () ->
                        ReadClipper.clipLowQualEnds(copy(record), lowQual,
                                ClippingRepresentation.SOFTCLIP_BASES));
                emit(index, "writeNsLowQual@" + lowQual, () ->
                        ReadClipper.clipLowQualEnds(copy(record), lowQual,
                                ClippingRepresentation.WRITE_NS));
                emit(index, "writeQ0sLowQual@" + lowQual, () ->
                        ReadClipper.clipLowQualEnds(copy(record), lowQual,
                                ClippingRepresentation.WRITE_Q0S));
                emit(index, "writeNsQ0sLowQual@" + lowQual, () ->
                        ReadClipper.clipLowQualEnds(copy(record), lowQual,
                                ClippingRepresentation.WRITE_NS_Q0S));
            }
        }
    }

    /**
     * A fresh adapter over a fresh SAMRecord for every call.
     *
     * ReadClipper copies before it mutates, but the corpus record itself must not carry a change
     * from one operation into the next: the rows would then describe a read nobody asked about.
     */
    static GATKRead copy(final SAMRecord record) {
        return new SAMRecordToGATKReadAdapter(record.deepCopy());
    }

    static void emit(final int index, final String operation, final Supplier<GATKRead> op) {
        String rendered;
        try {
            rendered = render(op.get());
        } catch (final Exception | AssertionError e) {
            rendered = "E";
        }
        System.out.printf("clipped\t%d\t%s\t%s%n", index, operation, rendered);
    }

    /** The whole read, in the order the port reads it. */
    static String render(final GATKRead read) {
        final StringBuilder quals = new StringBuilder();
        for (final byte q : read.getBaseQualities()) {
            if (quals.length() != 0) quals.append(',');
            quals.append(q);
        }
        final StringBuilder tags = new StringBuilder();
        for (final SAMRecord.SAMTagAndValue tag
                : ((SAMRecordToGATKReadAdapter) read).getEncapsulatedSamRecord().getAttributes()) {
            if (tags.length() != 0) tags.append(';');
            tags.append(tag.tag).append('=').append(ReadFilterDump.tagValue(tag.value));
        }
        return String.join("|",
                String.valueOf(read.getStart()),
                String.valueOf(((SAMRecordToGATKReadAdapter) read)
                        .getEncapsulatedSamRecord().getFlags()),
                read.getCigar().toString(),
                new String(read.getBases()),
                quals.toString(),
                String.valueOf(read.getMappingQuality()),
                tags.toString());
    }
}
