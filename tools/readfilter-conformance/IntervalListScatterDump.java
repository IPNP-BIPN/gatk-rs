/*
 * The five ways an interval list is scattered, taken from the reference.
 *
 * `SplitIntervals` and every scatter-gather pipeline built on it divide an interval list into
 * shards through `picard.util.IntervalList.IntervalListScatterMode`. The five modes disagree about
 * what a shard is, and the disagreements are not small.
 *
 * Nine behaviours this is built to catch.
 *
 *   - ONLY INTERVAL_SUBDIVISION UNIQUES ITS INPUT. Its `preprocessIntervalList` is `uniqued()`,
 *     which merges abutting and overlapping intervals and drops the names; every other mode's is
 *     `sorted()`, which keeps them apart. The same overlapping input therefore scatters into
 *     different intervals, not merely into different shards;
 *   - THE IDEAL WEIGHT IS FLOORED AT ONE and computed with `floorDiv` on the UNIQUE base count,
 *     even for the modes that never unique their list, so a list with overlaps has an ideal
 *     smaller than its own weight;
 *   - THE NO-SUBDIVISION MODES RAISE THE IDEAL TO THE WIDEST INTERVAL, which means a scatter count
 *     larger than the list can honour comes back with fewer shards than asked for;
 *   - THE LAST SHARD TAKES EVERYTHING LEFT. The loop stops offering intervals once it has returned
 *     `scatterCount - 1` lists and flushes the queue into the last one, whatever its weight;
 *   - THE PROJECTED REMAINING SIZE IS A DOUBLE DIVISION BY `scatterCount - intervalsReturned`,
 *     which is what the overflow mode compares against, so its decision depends on how many shards
 *     have already been returned rather than on the list alone;
 *   - INTERVAL_COUNT TAKES WHILE THERE IS ROOM, `idealSplitWeight - currentSize > 0`, so a
 *     remainder lands in the last shard;
 *   - INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER COMPARES THE PROJECTION AGAINST THE CURRENT SIZE
 *     instead, `projectedSizeOfRemaining > currentSize`, which spreads the remainder over the
 *     early shards rather than the last;
 *   - A SCATTER COUNT LARGER THAN THE WEIGHT still produces shards, because the ideal is floored
 *     at one and the queue empties;
 *   - AND THE SUBDIVISION CUT KEEPS THE NAME AND THE STRAND of the interval it cut, so a shard
 *     boundary inside a named interval leaves two intervals of that name.
 *
 * Output:
 *
 *     weight\t<list>,<mode>,<count>=<ideal split weight>
 *     shards\t<list>,<mode>,<count>=<number of shards>
 *     shard\t<list>,<mode>,<count>,<index>=<contig:start-end:strand:name;...>
 *
 * The intervals of a shard are separated by semicolons rather than by pipes, because `uniqued()`
 * joins the names of the intervals it merged with a pipe and a shard row would otherwise be
 * ambiguous.
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: IntervalListScatterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;

import java.util.ArrayList;
import java.util.List;

public class IntervalListScatterDump {

    static final int CONTIG_LENGTH = 10000;

    /** The interval lists this dump scatters, by name. */
    static IntervalList list(final String name) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : new String[] {"chr1", "chr2"}) {
            dictionary.addSequence(new SAMSequenceRecord(contig, CONTIG_LENGTH));
        }
        header.setSequenceDictionary(dictionary);
        final IntervalList list = new IntervalList(header);
        switch (name) {
            case "even":
                // Five intervals of a hundred bases each, which every mode can divide evenly.
                for (int i = 0; i < 5; i++) {
                    list.add(new Interval("chr1", 1 + i * 200, 100 + i * 200, false, "even" + i));
                }
                break;
            case "uneven":
                // One interval far wider than the rest, which is what raises the ideal weight of
                // the no-subdivision modes.
                list.add(new Interval("chr1", 1, 1, false, "one"));
                list.add(new Interval("chr1", 101, 110, false, "ten"));
                list.add(new Interval("chr1", 201, 1200, false, "thousand"));
                list.add(new Interval("chr1", 1301, 1350, false, "fifty"));
                break;
            case "overlapping":
                // Two overlaps and one abuttal, which only the uniquing mode merges.
                list.add(new Interval("chr1", 1, 100, false, "a"));
                list.add(new Interval("chr1", 50, 150, false, "b"));
                list.add(new Interval("chr1", 151, 200, false, "c"));
                list.add(new Interval("chr1", 400, 500, false, "d"));
                break;
            case "two-contigs":
                list.add(new Interval("chr1", 1, 100, false, "first"));
                list.add(new Interval("chr2", 1, 100, false, "second"));
                list.add(new Interval("chr2", 201, 400, false, "third"));
                break;
            case "unsorted":
                // Out of order and on the second contig first, so the sort is visible.
                list.add(new Interval("chr2", 1, 50, false, "z"));
                list.add(new Interval("chr1", 301, 400, true, "y"));
                list.add(new Interval("chr1", 1, 100, false, "x"));
                break;
            case "single":
                list.add(new Interval("chr1", 1, 1000, false, "only"));
                break;
            case "one-base":
                list.add(new Interval("chr1", 1, 1, false, "tiny"));
                break;
            default:
                throw new IllegalArgumentException("no such list: " + name);
        }
        return list;
    }

    public static void main(final String[] args) {
        System.out.println("# IntervalListScatterDump: the five ways an interval list is scattered");

        final String[] names = {"even", "uneven", "overlapping", "two-contigs", "unsorted",
                "single", "one-base"};
        final int[] counts = {1, 2, 3, 5, 10};

        for (final String name : names) {
            for (final IntervalListScatterMode mode : IntervalListScatterMode.values()) {
                for (final int count : counts) {
                    scatter(name, mode, count);
                }
            }
        }

        // The two refusals a caller can reach: a scatter count of zero and a negative one, both
        // through the ideal weight's own division.
        error("zero-count", "even", IntervalListScatterMode.INTERVAL_SUBDIVISION, 0);
        error("negative-count", "even", IntervalListScatterMode.INTERVAL_SUBDIVISION, -1);
    }

    static void scatter(final String name, final IntervalListScatterMode mode, final int count) {
        final IntervalList list = list(name);
        final IntervalListScatterer scatterer = mode.make();
        final String label = name + "," + mode.name() + "," + count;
        System.out.printf("weight\t%s=%d%n", label,
                scatterer.deduceIdealSplitWeight(scatterer.preprocessIntervalList(list), count));
        final List<IntervalList> shards = scatterer.scatter(list, count);
        System.out.printf("shards\t%s=%d%n", label, shards.size());
        for (int i = 0; i < shards.size(); i++) {
            System.out.printf("shard\t%s,%d=%s%n", label, i, render(shards.get(i)));
        }
    }

    /** One shard, as `contig:start-end:strand:name` per interval. */
    static String render(final IntervalList shard) {
        final List<String> parts = new ArrayList<>();
        for (final Interval interval : shard) {
            parts.add(String.format("%s:%d-%d:%s:%s", interval.getContig(), interval.getStart(),
                    interval.getEnd(), interval.isNegativeStrand() ? "-" : "+", interval.getName()));
        }
        return String.join(";", parts);
    }

    static void error(final String label, final String name, final IntervalListScatterMode mode,
                      final int count) {
        try {
            final List<IntervalList> shards = mode.make().scatter(list(name), count);
            System.out.printf("unexpected\t%s\t%d%n", label, shards.size());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }
}
