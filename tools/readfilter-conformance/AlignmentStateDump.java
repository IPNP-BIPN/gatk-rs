/*
 * Every stop an AlignmentStateMachine makes, taken from the reference.
 *
 * This is the bottom of every locus-based tool. LocusIteratorByState runs one machine per read and
 * merges their positions into a pileup, so every depth, allele count and pileup annotation rests on
 * the offsets recorded here. Comparing them one read at a time is the cheapest place to catch a
 * divergence that would otherwise surface as a wrong number in a caller.
 *
 * One row per stepForwardOnGenome() call, so the cigars that stop in surprising places are visible
 * rather than summarised:
 *
 *   - a deletion returns once per reference base it spans, with the read offset frozen;
 *   - I, S, H and P are consumed whole inside one step, so no row ever lands on them and the
 *     caller sees the base after an insertion instead;
 *   - a read ending on an insertion advances the genome offset one past its end anyway;
 *   - a cigar starting or ending with a deletion is a malformed read and throws, while a
 *     zero-length element is skipped in silence.
 *
 * Output:
 *
 *     step\t<cigar>\t<n>\t<operator>|<read offset>|<genome offset>|<genome position>|<cigar element offset>|<offset into element>|<left edge>|<right edge>
 *     end\t<cigar>\t<ok|E:class>\t<number of steps>
 *
 * The reads are built here rather than read from a file: a cigar is the whole input, and a fixture
 * BAM would hide it behind an encoder.
 *
 * Usage: AlignmentStateDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import htsjdk.samtools.CigarOperator;

public class AlignmentStateDump {

    static final int CONTIG_LENGTH = 200;

    /**
     * The cigars probed, chosen for the machine's corners rather than for realism. The read's
     * length follows from the cigar, so every one of these is a complete input.
     */
    static final String[] CIGARS = {
        "10M",
        "1M",
        // Deletions: one stop per reference base, in the middle, and next to an insertion.
        "5M3D5M",
        "5M1D5M",
        "5M3I5M",
        "5M3D3I5M",
        "5M3I3D5M",
        // Skips, which look like deletions and mean something else.
        "5M10N5M",
        // Clipping on either side, and both at once.
        "3S7M",
        "7M3S",
        "3S4M3S",
        "3H7M",
        "7M3H",
        "3H2S5M2S3H",
        // Reads that end on something other than a match.
        "5M3I",
        "5M3S",
        "5M3H",
        // Reads that start on something other than a match.
        "3I5M",
        // Padding, which consumes neither the read nor the reference.
        "5M3P5M",
        // The two malformed cigars the machine refuses.
        "3D5M",
        "5M3D",
        // Zero-length elements, which are skipped rather than refused.
        "0M10M",
        "5M0D5M",
        "5M0I5M",
        // Sequence match and mismatch operators, which behave as M.
        "5=5X",
        // A long deletion, so the per-base repetition is unmistakable.
        "2M10D2M",
    };

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                java.util.List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));

        System.out.println("# AlignmentStateDump: every stop an AlignmentStateMachine makes");
        for (final String cigar : CIGARS) {
            dump(header, cigar);
        }
    }

    static void dump(final SAMFileHeader header, final String cigarText) {
        int steps = 0;
        String outcome = "ok";
        try {
            // Inside the try: a zero-length cigar element may be refused by the record itself
            // rather than by the machine, and which of the two refuses is worth recording.
            final GATKRead read = makeRead(header, cigarText);
            final AlignmentStateMachine machine = new AlignmentStateMachine(read);
            while (true) {
                final CigarOperator op = machine.stepForwardOnGenome();
                System.out.printf("step\t%s\t%d\t%s|%d|%d|%d|%d|%d|%s|%s%n",
                        cigarText,
                        steps,
                        op == null ? "null" : op.toString(),
                        machine.getReadOffset(),
                        machine.getGenomeOffset(),
                        machine.getGenomePosition(),
                        machine.getCurrentCigarElementOffset(),
                        machine.getOffsetIntoCurrentCigarElement(),
                        machine.isLeftEdge(),
                        machine.isRightEdge());
                steps++;
                if (op == null) {
                    break;
                }
            }
        } catch (final Exception e) {
            // The class, not the message: the message quotes the read, and the class is what
            // distinguishes the two malformed cigars from an ordinary end of read.
            outcome = "E:" + e.getClass().getName();
        }
        System.out.printf("end\t%s\t%s\t%d%n", cigarText, outcome, steps);
    }

    /** A read at chr1:101 whose bases and qualities are as long as the cigar requires. */
    static GATKRead makeRead(final SAMFileHeader header, final String cigarText) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("read-" + cigarText);
        record.setReferenceName("chr1");
        record.setAlignmentStart(101);
        record.setCigar(TextCigarCodec.decode(cigarText));
        final int length = record.getCigar().getReadLength();
        final byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        for (int i = 0; i < length; i++) {
            bases[i] = "ACGT".getBytes()[i % 4];
            quals[i] = (byte) 30;
        }
        record.setReadBases(bases);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        return new SAMRecordToGATKReadAdapter(record);
    }
}
