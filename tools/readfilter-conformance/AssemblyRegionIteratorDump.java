/*
 * The assembly-region traversal itself, taken from the reference.
 *
 * This is where the activity profile, the assembly region and the locus iterators meet: loci go in,
 * assembly regions come out, and every region carries the reads overlapping its padded span.
 * HaplotypeCaller and Mutect2 are this loop plus an assembler, so a region boundary that moves by
 * one base here changes their calls.
 *
 * Four behaviours are measured, and three of them are invisible from any tool's output.
 *
 *   - a shard has TWO interval lists and they are built the same way. MultiIntervalLocalReadShard
 *     passes both through IntervalUtils.getIntervalsWithFlanks, the padded one with the padding and
 *     the unpadded one with ZERO. Passing zero is not a no-op: the function sorts and merges with
 *     IntervalMergingRule.ALL whatever the padding is. So getIntervals() is already sorted and
 *     already has its adjacent intervals joined, and the reads are queried over one list while the
 *     loci are walked over the other;
 *   - the profile is popped BEFORE the current pileup is added, and the force flag is
 *     "this locus does not continue the profile". The upstream comment says "Ordering matters
 *     here". Adding first would make the profile contiguous with the new pileup, the force would
 *     never fire, and a region would never be closed at a gap between two intervals;
 *   - a region popped from the profile is not ready. It waits in a queue until the loci have
 *     advanced past the end of its PADDED span, because until then the reads that belong in it have
 *     not been read;
 *   - reads are carried forward from the previous region as well as taken from the cache, so a read
 *     spanning a boundary is in both regions, and a read is only left in the cache when its START
 *     is past the region's padded end.
 *
 * The evaluator is a probe, not a caller: activity is a declared function of the pileup depth, so
 * the regions are decided by the traversal rather than by anything statistical.
 *
 * Output:
 *
 *     shard\t<label>\t<intervals>\t<paddedIntervals>
 *     args\t<label>\t<ok|E:class:message>
 *     region\t<label>\t<n>\t<span>|<paddedSpan>|<isActive>|<nReads>
 *     rread\t<label>\t<n>\t<read names, comma-separated>
 *     count\t<label>\t<number of regions>
 *     error\t<label>\t<class>\t<message>
 *
 * Usage: AssemblyRegionIteratorDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class AssemblyRegionIteratorDump {

    static final int CHR1 = 1000;
    static final int CHR2 = 500;

    /**
     * Reads laid out so the depth, and therefore the activity, rises and falls: a covered block at
     * 101-140, a hole, a second covered block at 201-230, and a lone read at 400.
     */
    static final String[][] READS = {
        {"a1", "rg1", "20M", "101"},
        {"a2", "rg1", "20M", "111"},
        {"a3", "rg2", "20M", "121"},
        {"b1", "rg1", "20M", "201"},
        {"b2", "rg2", "20M", "211"},
        {"c1", "rg1", "20M", "400"},
    };

    /** A concrete argument collection, since the base class is abstract in name only here. */
    public static class Args extends AssemblyRegionArgumentCollection {
        private static final long serialVersionUID = 1L;
        @Override protected int defaultMinAssemblyRegionSize() { return 5; }
        @Override protected int defaultMaxAssemblyRegionSize() { return 50; }
        @Override protected int defaultAssemblyRegionPadding() { return 10; }
        @Override protected int defaultMaxReadsPerAlignmentStart() { return 50; }
        @Override protected double defaultActiveProbThreshold() { return 0.002; }
        @Override protected int defaultMaxProbPropagationDistance() { return 50; }
    }

    public static void main(final String[] args) {
        final SAMFileHeader header = header();

        System.out.println("# AssemblyRegionIteratorDump: the assembly-region traversal");

        // The two interval lists of a shard, over inputs that are unsorted, adjacent, overlapping
        // and far apart, so the sorting, the merging at zero padding and the merging after padding
        // are each visible on their own.
        shard("shard-one", 10, new SimpleInterval("chr1", 100, 200));
        shard("shard-adjacent", 10,
                new SimpleInterval("chr1", 100, 200), new SimpleInterval("chr1", 201, 300));
        shard("shard-overlapping", 10,
                new SimpleInterval("chr1", 100, 220), new SimpleInterval("chr1", 200, 300));
        shard("shard-unsorted", 10,
                new SimpleInterval("chr1", 300, 400), new SimpleInterval("chr1", 100, 200));
        // Fifteen bases apart: separate without padding, merged with ten bases of it on each side.
        shard("shard-merged-by-padding", 10,
                new SimpleInterval("chr1", 100, 200), new SimpleInterval("chr1", 216, 300));
        // Padding that runs off the front of the contig.
        shard("shard-off-contig", 10, new SimpleInterval("chr1", 5, 20));
        shard("shard-zero-padding", 0,
                new SimpleInterval("chr1", 100, 200), new SimpleInterval("chr1", 201, 300));
        shard("shard-negative-padding", -1, new SimpleInterval("chr1", 100, 200));

        // The argument collection's validation, in the order it checks.
        validate("args-default", 5, 50, 10, 50, 20, 75);
        validate("args-zero-min", 0, 50, 10, 50, 20, 75);
        validate("args-zero-max", 5, 0, 10, 50, 20, 75);
        validate("args-min-above-max", 60, 50, 10, 50, 20, 75);
        validate("args-negative-padding", 5, 50, -1, 50, 20, 75);
        validate("args-negative-max-reads", 5, 50, 10, -1, 20, 75);
        validate("args-negative-snp-padding", 5, 50, 10, 50, -1, 75);
        validate("args-negative-indel-padding", 5, 50, 10, 50, 20, -1);
        // Two things wrong at once, to fix which one is reported.
        validate("args-two-wrong", 0, 50, -1, 50, 20, 75);

        // The traversal. The threshold is what decides how much of the covered block is active, so
        // it is varied rather than fixed.
        traverse(header, "trav-whole-contig", 10, 5, 50,
                new SimpleInterval("chr1", 1, CHR1));
        traverse(header, "trav-covered-only", 10, 5, 50,
                new SimpleInterval("chr1", 101, 140));
        // Two intervals with a gap between them, which is where the force conversion fires.
        traverse(header, "trav-two-intervals", 10, 5, 50,
                new SimpleInterval("chr1", 101, 140), new SimpleInterval("chr1", 201, 230));
        // The same two intervals close enough to merge once padded, so the gap disappears.
        traverse(header, "trav-merged-intervals", 40, 5, 50,
                new SimpleInterval("chr1", 101, 140), new SimpleInterval("chr1", 201, 230));
        // A maximum region size below the active stretch, which forces a cut at a local minimum.
        traverse(header, "trav-small-max", 10, 5, 20,
                new SimpleInterval("chr1", 1, CHR1));
        // A minimum region size above the whole interval.
        traverse(header, "trav-large-min", 10, 400, 500,
                new SimpleInterval("chr1", 101, 140));
        // No padding at all, so the region's two spans coincide.
        traverse(header, "trav-no-padding", 0, 5, 50,
                new SimpleInterval("chr1", 101, 140));
        // An interval on a contig with no reads: every locus is an empty pileup.
        traverse(header, "trav-no-reads", 10, 5, 50,
                new SimpleInterval("chr2", 100, 140));
        // The lone read at 400, whose region has to be closed by the final force conversion.
        traverse(header, "trav-lone-read", 10, 5, 50,
                new SimpleInterval("chr1", 395, 425));
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1), new SAMSequenceRecord("chr2", CHR2))));
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        rg1.setSample("sampleA");
        header.addReadGroup(rg1);
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        rg2.setSample("sampleB");
        header.addReadGroup(rg2);
        return header;
    }

    /**
     * The shard's two interval lists, taken straight from IntervalUtils rather than through a
     * ReadsDataSource: the shard's constructor is exactly these two calls, and building a data
     * source here would add a BAM to the fixture without adding an observation.
     */
    static void shard(final String label, final int padding, final SimpleInterval... intervals) {
        try {
            final SAMSequenceDictionary dictionary = header().getSequenceDictionary();
            if (padding < 0) {
                // The shard's own precondition, which IntervalUtils does not make.
                throw new IllegalArgumentException("intervalPadding must be >= 0");
            }
            final List<SimpleInterval> plain =
                    IntervalUtils.getIntervalsWithFlanks(List.of(intervals), 0, dictionary);
            final List<SimpleInterval> padded =
                    IntervalUtils.getIntervalsWithFlanks(List.of(intervals), padding, dictionary);
            System.out.printf("shard\t%s\t%s\t%s%n", label, join(plain), join(padded));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static String join(final List<SimpleInterval> intervals) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final SimpleInterval interval : intervals) {
            joiner.add(interval.toString());
        }
        return joiner.toString();
    }

    static void validate(final String label, final int min, final int max, final int padding,
                         final int maxReads, final int snpPadding, final int indelPadding) {
        final Args args = new Args();
        args.minAssemblyRegionSize = min;
        args.maxAssemblyRegionSize = max;
        args.assemblyRegionPadding = padding;
        args.maxReadsPerAlignmentStart = maxReads;
        args.snpPaddingForGenotyping = snpPadding;
        args.indelPaddingForGenotyping = indelPadding;
        try {
            args.validate();
            System.out.printf("args\t%s\tok%n", label);
        } catch (final Exception | AssertionError e) {
            System.out.printf("args\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /**
     * Run the traversal by hand rather than through AssemblyRegionWalker: the walker adds a command
     * line, a BAM on disk and a reference, none of which changes which regions come out. The
     * iterator is given the same shard, header, arguments and evaluator the walker would give it.
     */
    static void traverse(final SAMFileHeader header, final String label, final int padding,
                         final int minSize, final int maxSize, final SimpleInterval... intervals) {
        try {
            final Args args = new Args();
            args.assemblyRegionPadding = padding;
            args.minAssemblyRegionSize = minSize;
            args.maxAssemblyRegionSize = maxSize;

            final List<GATKRead> reads = new ArrayList<>();
            for (final String[] spec : READS) {
                reads.add(read(header, spec[0], spec[1], spec[2], Integer.parseInt(spec[3])));
            }

            final List<AssemblyRegion> regions = ProbeTraversal.run(
                    header, reads, List.of(intervals), args);

            for (int i = 0; i < regions.size(); i++) {
                final AssemblyRegion region = regions.get(i);
                System.out.printf("region\t%s\t%d\t%s|%s|%b|%d%n", label, i, region.getSpan(),
                        region.getPaddedSpan(), region.isActive(), region.size());
                final StringJoiner names = new StringJoiner(",");
                for (final GATKRead read : region.getReads()) {
                    names.add(read.getName());
                }
                System.out.printf("rread\t%s\t%d\t%s%n", label, i, names);
            }
            System.out.printf("count\t%s\t%d%n", label, regions.size());
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /** The probe evaluator: activity is a declared function of the pileup depth. */
    static double activity(final AlignmentContext pileup) {
        final int depth = pileup.getBasePileup().size();
        return depth >= 2 ? 1.0 : 0.0;
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

    /** Kept apart so the traversal's own imports stay next to it. */
    static class ProbeTraversal {
        static List<AssemblyRegion> run(final SAMFileHeader header, final List<GATKRead> reads,
                                        final List<SimpleInterval> intervals,
                                        final AssemblyRegionArgumentCollection args) {
            final org.broadinstitute.hellbender.engine.MultiIntervalShard<GATKRead> shard =
                    new InMemoryShard(header, reads, intervals, args.assemblyRegionPadding);
            final org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator evaluator =
                    new org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator() {
                        @Override
                        public ActivityProfileState isActive(final AlignmentContext locusPileup,
                                                             final ReferenceContext referenceContext,
                                                             final FeatureContext featureContext) {
                            return new ActivityProfileState(
                                    new SimpleInterval(locusPileup.getContig(),
                                            (int) locusPileup.getPosition(),
                                            (int) locusPileup.getPosition()),
                                    activity(locusPileup));
                        }
                    };

            final List<AssemblyRegion> out = new ArrayList<>();
            final java.util.Iterator<AssemblyRegion> iterator =
                    new org.broadinstitute.hellbender.engine.AssemblyRegionIterator(
                            shard, header, null, null, evaluator, args, false);
            while (iterator.hasNext()) {
                out.add(iterator.next());
            }
            return out;
        }
    }

    /**
     * A MultiIntervalShard backed by a list rather than by a BAM, with the shard's own interval
     * arithmetic reproduced: getIntervals is the input at zero padding and getPaddedIntervals is
     * the input at the assembly-region padding, both through getIntervalsWithFlanks.
     */
    static class InMemoryShard
            implements org.broadinstitute.hellbender.engine.MultiIntervalShard<GATKRead> {
        final List<GATKRead> reads;
        final List<SimpleInterval> intervals;
        final List<SimpleInterval> paddedIntervals;

        InMemoryShard(final SAMFileHeader header, final List<GATKRead> reads,
                      final List<SimpleInterval> intervals, final int padding) {
            final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
            this.intervals = IntervalUtils.getIntervalsWithFlanks(intervals, 0, dictionary);
            this.paddedIntervals =
                    IntervalUtils.getIntervalsWithFlanks(intervals, padding, dictionary);
            final List<GATKRead> kept = new ArrayList<>();
            for (final GATKRead read : reads) {
                for (final SimpleInterval interval : this.paddedIntervals) {
                    if (interval.overlaps(read)) {
                        kept.add(read);
                        break;
                    }
                }
            }
            this.reads = kept;
        }

        @Override
        public List<SimpleInterval> getIntervals() {
            return intervals;
        }

        @Override
        public List<SimpleInterval> getPaddedIntervals() {
            return paddedIntervals;
        }

        @Override
        public java.util.Iterator<GATKRead> iterator() {
            return reads.iterator();
        }
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
