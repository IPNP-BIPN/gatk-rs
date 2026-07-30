/*
 * The iterators between a shard of reads and an activity profile, taken from the reference.
 *
 * LocusIteratorByState yields nothing at all for a locus no read covers: a position with zero depth
 * is silently absent from its output. The assembly-region traversal cannot work that way, because
 * the activity profile has to see the gap in order to close a region over it. So AssemblyRegionIterator
 * wraps the locus iterator in IntervalAlignmentContextIterator, whose whole job is to walk every
 * base of every interval and manufacture an empty AlignmentContext wherever the wrapped iterator
 * has none. The comment upstream is unusually direct about the stakes: "This is critical for
 * reproducing GATK 3.x behavior!"
 *
 * Three things are dumped, smallest first, because the third is built out of the first two.
 *
 *   - ShardedIntervalIterator, whose arithmetic is absolute rather than relative: it computes shard
 *     indices over the interval's LENGTH with IntervalUtils.shardIndex and converts each index back
 *     to coordinates, so the last shard is truncated against the interval's end. At shard size 1,
 *     the only size IntervalLocusIterator ever uses, the arithmetic is invisible;
 *   - IntervalLocusIterator, which is that iterator at size 1;
 *   - IntervalAlignmentContextIterator, over a corpus with deliberate coverage gaps at both ends
 *     and in the middle, and over intervals that start before the reads, end after them, sit
 *     entirely inside a gap, and cross a contig boundary.
 *
 * The pileup is dumped as its size and its read names, not as an object: what is being measured is
 * which locus got which pileup, and whether an empty one was manufactured there.
 *
 * Output:
 *
 *     shard\t<label>\t<n>\t<contig>:<start>-<end>
 *     locus\t<label>\t<n>\t<contig>:<start>-<end>
 *     ctx\t<label>\t<n>\t<contig>:<pos>\t<size>\t<read names>
 *     count\t<label>\t<number of items>
 *     error\t<label>\t<class>\t<message>
 *
 * Usage: IntervalIteratorsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.iterators.IntervalLocusIterator;
import org.broadinstitute.hellbender.utils.iterators.ShardedIntervalIterator;
import org.broadinstitute.hellbender.utils.locusiterator.IntervalAlignmentContextIterator;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.StringJoiner;

public class IntervalIteratorsDump {

    static final int CHR1 = 300;
    static final int CHR2 = 200;

    /**
     * Reads with deliberate gaps: nothing before 101, a hole at 111-119, and nothing after 140.
     * The holes are what the wrapper exists for.
     */
    static final String[][] READS = {
        {"r1", "rg1", "10M", "101"},
        {"r2", "rg1", "10M", "105"},
        {"r3", "rg1", "10M", "120"},
        {"r4", "rg2", "10M", "131"},
    };

    public static void main(final String[] args) {
        final SAMFileHeader header = header();

        System.out.println("# IntervalIteratorsDump: the iterators between a shard and a profile");

        // ShardedIntervalIterator. Sizes either side of the interval length, and a size that does
        // not divide it, which is where the truncation of the last shard shows.
        for (final int size : new int[] {1, 2, 3, 4, 11, 100}) {
            shards("shard-10-" + size, size, new SimpleInterval("chr1", 10, 20));
        }
        // An interval of one base, and two intervals in a row, so the reset between them is visible.
        shards("shard-one-base", 3, new SimpleInterval("chr1", 10, 10));
        shards("shard-two-intervals", 4,
                new SimpleInterval("chr1", 10, 20), new SimpleInterval("chr2", 5, 9));
        // A shard size the class refuses.
        shards("shard-zero", 0, new SimpleInterval("chr1", 10, 20));
        shards("shard-negative", -1, new SimpleInterval("chr1", 10, 20));
        // No intervals at all.
        shards("shard-empty", 5);

        // IntervalLocusIterator, which is the same class at size 1.
        loci("loci-one", new SimpleInterval("chr1", 10, 14));
        loci("loci-two-contigs",
                new SimpleInterval("chr1", 10, 12), new SimpleInterval("chr2", 5, 6));
        loci("loci-adjacent",
                new SimpleInterval("chr1", 10, 12), new SimpleInterval("chr1", 13, 14));
        loci("loci-empty");

        // IntervalAlignmentContextIterator. Each case names a different relationship between the
        // intervals asked for and the coverage that exists.
        contexts(header, "ctx-covered", new SimpleInterval("chr1", 101, 110));
        // Starts before any read, so the leading loci are manufactured.
        contexts(header, "ctx-leading-gap", new SimpleInterval("chr1", 95, 106));
        // Ends after every read, so the trailing loci are manufactured.
        contexts(header, "ctx-trailing-gap", new SimpleInterval("chr1", 135, 145));
        // Spans the hole in the middle, which is the case the wrapper exists for.
        contexts(header, "ctx-interior-gap", new SimpleInterval("chr1", 108, 122));
        // Entirely inside the hole: every locus is manufactured and none is covered.
        contexts(header, "ctx-all-gap", new SimpleInterval("chr1", 112, 118));
        // Two intervals with a large uncovered stretch between them, which is where the
        // uncovered branch has to decide whether to advance the wrapped iterator.
        contexts(header, "ctx-two-intervals",
                new SimpleInterval("chr1", 101, 103), new SimpleInterval("chr1", 130, 133));
        // A contig with no reads at all.
        contexts(header, "ctx-other-contig", new SimpleInterval("chr2", 10, 14));
        // Both contigs in one traversal, in dictionary order.
        contexts(header, "ctx-both-contigs",
                new SimpleInterval("chr1", 138, 141), new SimpleInterval("chr2", 10, 12));
        // The whole of chr1, so the leading gap, the interior gap and the trailing gap are all in
        // one traversal.
        contexts(header, "ctx-whole-contig", new SimpleInterval("chr1", 1, CHR1));
        contexts(header, "ctx-empty");
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1), new SAMSequenceRecord("chr2", CHR2))));
        header.addReadGroup(readGroup("rg1", "sampleA"));
        header.addReadGroup(readGroup("rg2", "sampleB"));
        return header;
    }

    static SAMReadGroupRecord readGroup(final String id, final String sample) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        group.setSample(sample);
        return group;
    }

    static void shards(final String label, final int size, final SimpleInterval... intervals) {
        try {
            final ShardedIntervalIterator iterator =
                    new ShardedIntervalIterator(List.of(intervals).iterator(), size);
            int index = 0;
            while (iterator.hasNext()) {
                System.out.printf("shard\t%s\t%d\t%s%n", label, index++, iterator.next());
            }
            System.out.printf("count\t%s\t%d%n", label, index);
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static void loci(final String label, final SimpleInterval... intervals) {
        try {
            final IntervalLocusIterator iterator =
                    new IntervalLocusIterator(List.of(intervals).iterator());
            int index = 0;
            while (iterator.hasNext()) {
                System.out.printf("locus\t%s\t%d\t%s%n", label, index++, iterator.next());
            }
            System.out.printf("count\t%s\t%d%n", label, index);
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static void contexts(final SAMFileHeader header, final String label,
                         final SimpleInterval... intervals) {
        try {
            final List<GATKRead> reads = new ArrayList<>();
            for (final String[] spec : READS) {
                reads.add(read(header, spec[0], spec[1], spec[2], Integer.parseInt(spec[3])));
            }

            final LocusIteratorByState libs = new LocusIteratorByState(reads.iterator(),
                    DownsamplingMethod.NONE, ReadUtils.getSamplesFromHeader(header), header, true);
            final IntervalLocusIterator locusIterator =
                    new IntervalLocusIterator(List.of(intervals).iterator());
            final Iterator<AlignmentContext> iterator = new IntervalAlignmentContextIterator(
                    libs, locusIterator, header.getSequenceDictionary());

            int index = 0;
            while (iterator.hasNext()) {
                final AlignmentContext context = iterator.next();
                final StringJoiner names = new StringJoiner(",");
                for (final PileupElement element : context.getBasePileup()) {
                    names.add(element.getRead().getName());
                }
                System.out.printf("ctx\t%s\t%d\t%s:%d\t%d\t%s%n", label, index++,
                        context.getContig(), context.getPosition(),
                        context.getBasePileup().size(), names);
            }
            System.out.printf("count\t%s\t%d%n", label, index);
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static GATKRead read(final SAMFileHeader header, final String name, final String group,
                         final String cigar, final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigar(TextCigarCodec.decode(cigar));
        record.setMappingQuality(60);
        record.setAttribute("RG", group);
        final int length = record.getCigar().getReadLength();
        final byte[] bases = new byte[length];
        final byte[] quals = new byte[length];
        for (int i = 0; i < length; i++) {
            bases[i] = "ACGT".getBytes()[i % 4];
            quals[i] = 30;
        }
        record.setReadBases(bases);
        record.setBaseQualities(quals);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
