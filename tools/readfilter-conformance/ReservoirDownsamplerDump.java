/*
 * Which reads a ReservoirDownsampler keeps, taken from the reference.
 *
 * This decides what a deep pileup contains, and four of its behaviours change the answer without
 * being inherent to reservoir sampling:
 *
 *   - the draw happens even when the slot is discarded. nextInt(totalReadsSeen) is called for every
 *     read past the target, and only then is the slot compared against the target. So a read that
 *     is thrown away still advances the shared random stream, and every later draw anywhere in the
 *     run depends on how many reads were seen here. That is why the probe reports the generator's
 *     position afterwards, not only the reads it kept;
 *   - the bound grows with the count, so the same read at the same place draws differently
 *     depending on how many preceded it;
 *   - replacement is in place, so the reservoir's order is the order slots were last written rather
 *     than the reads' order;
 *   - setNonRandomReplacementMode swaps the draw for Math.abs(name.hashCode()) % totalReadsSeen,
 *     and Math.abs is wrong for Integer.MIN_VALUE, so a name hashing to it yields a negative slot.
 *
 * The generator is reset before each case, because Utils.getRandomGenerator() is a single static
 * stream and a case that inherited the previous one's position would measure the order of the
 * cases rather than the downsampler.
 *
 * Output:
 *
 *     keep\t<label>\t<comma-separated read names, in reservoir order>
 *     stats\t<label>\t<size>\t<discarded>\t<next int from the shared generator>
 *     error\t<label>\t<class>
 *
 * Usage: ReservoirDownsamplerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.ReservoirDownsampler;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class ReservoirDownsamplerDump {

    static final int CONTIG_LENGTH = 1000;

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));

        System.out.println("# ReservoirDownsamplerDump: which reads the reservoir keeps");

        // Fewer reads than the target: nothing is drawn and nothing is discarded.
        probe(header, "under-target", 3, 10, false);
        // Exactly the target: still no draw, because the test is <=.
        probe(header, "at-target", 10, 10, false);
        // One over: exactly one draw, whose slot may or may not land inside the reservoir.
        probe(header, "one-over", 11, 10, false);
        // Well over, so the growing bound and the in-place replacement both show.
        probe(header, "many-over", 50, 10, false);
        probe(header, "very-many-over", 500, 10, false);
        // A target of one, where every read past the first competes for the same slot.
        probe(header, "target-one", 20, 1, false);
        // The non-random mode, which takes no draw at all and leaves the stream where it was.
        probe(header, "nonrandom", 50, 10, true);
        // A target of zero, which the constructor refuses.
        zeroTarget();
    }

    static void probe(final SAMFileHeader header, final String label, final int readCount,
                      final int target, final boolean nonRandom) {
        // The shared stream is reset per case, so each case measures the downsampler rather than
        // the order the cases happen to run in.
        Utils.resetRandomGenerator();

        final ReservoirDownsampler downsampler = new ReservoirDownsampler(target);
        downsampler.setNonRandomReplacementMode(nonRandom);

        final List<GATKRead> reads = new ArrayList<>();
        for (int i = 0; i < readCount; i++) {
            reads.add(makeRead(header, String.format("r%03d", i), 101 + i));
        }
        for (final GATKRead read : reads) {
            downsampler.submit(read);
        }
        downsampler.signalEndOfInput();

        final int discarded = downsampler.getNumberOfDiscardedItems();
        final List<GATKRead> kept = downsampler.consumeFinalizedItems();
        final StringJoiner names = new StringJoiner(",");
        for (final GATKRead read : kept) {
            names.add(read.getName());
        }
        System.out.printf("keep\t%s\t%s%n", label, names);
        // The generator's next value after the case, which is what shows that discarded reads still
        // consumed draws: two cases that kept the same reads can leave the stream in different
        // places.
        System.out.printf("stats\t%s\t%d\t%d\t%d%n",
                label, kept.size(), discarded, Utils.getRandomGenerator().nextInt());
    }

    static void zeroTarget() {
        try {
            new ReservoirDownsampler(0);
            System.out.printf("error\t%s\t%s%n", "zero-target", "none");
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s%n", "zero-target", e.getClass().getName());
        }
    }

    static GATKRead makeRead(final SAMFileHeader header, final String name, final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigar(TextCigarCodec.decode("10M"));
        final byte[] bases = new byte[10];
        final byte[] quals = new byte[10];
        for (int i = 0; i < 10; i++) {
            bases[i] = "ACGT".getBytes()[i % 4];
            quals[i] = (byte) 30;
        }
        record.setReadBases(bases);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        return new SAMRecordToGATKReadAdapter(record);
    }
}
