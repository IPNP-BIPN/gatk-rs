/*
 * The assembly region itself: its two spans, the reads it holds, and what trimming does to both.
 *
 * ActivityProfileDump settled where a region starts and stops. This one is what the region *is*
 * once it exists, which is the object HaplotypeCaller assembles: a primary span it calls variants
 * in, a padded span it assembles over, and the reads overlapping the padded one.
 *
 * Four behaviours are measured here because none of them is what the class documentation says.
 *
 *   - the padded span of the (span, padding) constructor goes through
 *     IntervalUtils.trimIntervalToContig, which returns NULL rather than throwing when the interval
 *     cannot be placed on the contig at all. The null then reaches Utils.nonNull inside the other
 *     constructor, so the failure is reported as a null padded span and never mentions padding;
 *   - trim(span, padding) does NOT produce the region its own javadoc describes. The javadoc worked
 *     an example: active 100-200 padded 50-250, trimmed to 150-225, "here we represent the assembly
 *     region as a region from 150-200 with 25 bp of padding". The code expands the REQUESTED span by
 *     the requested padding and intersects that with the old padded span, so the padding is not
 *     recomputed to fit and the answer is a different interval. The golden is the arbiter;
 *   - trim re-clips every read to the new padded span, drops the ones left empty or no longer
 *     overlapping, and then SORTS what survives with ReadCoordinateComparator. A region built by
 *     adding reads in one order can therefore come out of trim in another, and the comparator is
 *     part of the region's observable output rather than an internal detail;
 *   - trim carries the reads across and does NOT carry the hard-clipped pileup reads: they are
 *     dropped, because the new region is constructed empty and only addAll is called.
 *
 * add() is the other half. It validates in an order that matters: it builds a SimpleInterval from
 * the read BEFORE testing the overlap, so an unmapped read at position 0 fails on the interval's own
 * validation and the message never says the read is unmapped.
 *
 * Output:
 *
 *     ctor\t<label>\t<activeSpan>|<paddedSpan>|<isActive>|<nReads>
 *     add\t<label>\t<n>\t<ok|E:class:message>
 *     region\t<label>\t<activeSpan>|<paddedSpan>|<isActive>|<nReads>
 *     read\t<label>\t<n>\t<name>|<start>|<cigar>|<bases>
 *     trim\t<label>\t<activeSpan>|<paddedSpan>|<isActive>|<nReads>
 *     cmp\t<i>\t<j>\t<sign>
 *     error\t<label>\t<class>\t<message>
 *
 * Usage: AssemblyRegionDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;

public class AssemblyRegionDump {

    static SAMFileHeader header;
    static List<SAMRecord> corpus;

    public static void main(final String[] args) {
        header = ReadFilterDump.header();
        corpus = ReadFilterDump.corpus(header);

        System.out.println("# AssemblyRegionDump: the region itself, its spans, its reads, its trim");
        ReadFilterDump.printCorpus(header, corpus);

        // The comparator first, because trim's output order depends on it and nothing else in the
        // port has measured it. Every ordered pair, so the antisymmetry is measured too: a
        // comparator that is not antisymmetric sorts differently under a different algorithm.
        for (int i = 0; i < corpus.size(); i++) {
            for (int j = 0; j < corpus.size(); j++) {
                compare(i, j);
            }
        }

        // The (span, padding) constructor. chr1 is 2000 long, so the last three cases run off each
        // end and the very last one cannot be placed at all.
        ctor("pad-inside", "chr1", 500, 600, 100);
        ctor("pad-zero", "chr1", 500, 600, 0);
        ctor("pad-off-front", "chr1", 10, 20, 100);
        ctor("pad-off-back", "chr1", 1950, 1990, 100);
        ctor("pad-whole-contig", "chr1", 500, 600, 5000);
        ctor("pad-unknown-contig", "chrX", 500, 600, 10);
        // A single base, which is the smallest legal active span.
        ctor("pad-one-base", "chr1", 500, 500, 10);

        // The (span, paddedSpan) constructor, where the padding is not symmetric and where the
        // containment precondition is broken.
        ctorPair("pair-asymmetric", "chr1", 500, 600, "chr1", 400, 1000);
        ctorPair("pair-equal", "chr1", 500, 600, "chr1", 500, 600);
        ctorPair("pair-not-contained", "chr1", 500, 600, "chr1", 550, 700);
        ctorPair("pair-other-contig", "chr1", 500, 600, "chr2", 400, 700);

        // add(), in the order the corpus declares. The corpus is deliberately NOT sorted by start
        // near its end, so this measures the out-of-order refusal on real records.
        addAllInCorpusOrder("add-corpus-order");
        // The same reads sorted, which is the order that is accepted.
        addSorted("add-sorted");
        // An unmapped read, which fails on SimpleInterval's own validation.
        addOne("add-unmapped", indexOf("flag_unmapped"));
        addOne("add-zero-start", indexOf("zero_start"));
        addOne("add-no-reference", indexOf("no_reference"));
        // A read on the other contig, which fails the overlap test rather than the interval one.
        addOne("add-other-contig", indexOf("chr2_mapped"));

        // Trimming. The region is built wide enough to hold most of chr1's reads, then trimmed by
        // both entry points.
        final String wide = "wide";
        region(wide, "chr1", 100, 1900, 100);
        trim(wide, "trim-inside", "chr1", 500, 1000, 50);
        trim(wide, "trim-javadoc-example", "chr1", 150, 225, 25);
        trim(wide, "trim-zero-padding", "chr1", 500, 1000, 0);
        trim(wide, "trim-padding-past-region", "chr1", 500, 1000, 5000);
        trim(wide, "trim-single-base", "chr1", 1000, 1000, 10);
        trim(wide, "trim-partly-outside", "chr1", 1800, 2500, 10);
        trim(wide, "trim-disjoint", "chr1", 1, 50, 0);
        trim(wide, "trim-other-contig", "chr2", 100, 200, 10);
        trimPair(wide, "trimpair-inside", "chr1", 500, 1000, "chr1", 400, 1100);
        trimPair(wide, "trimpair-padded-not-containing", "chr1", 500, 1000, "chr1", 600, 1100);
        trimPair(wide, "trimpair-padded-beyond-original", "chr1", 500, 1000, "chr1", 1, 2000);

        // The javadoc's own example, built exactly as it is written, so the row can be compared to
        // the sentence it contradicts.
        region("javadoc", "chr1", 100, 200, 50);
        trim("javadoc", "javadoc-trim", "chr1", 150, 225, 25);

        // A narrow region, to show that trim can drop every read.
        region("narrow", "chr1", 1000, 1010, 5);
        trim("narrow", "narrow-trim", "chr1", 1005, 1008, 0);
    }

    static int indexOf(final String name) {
        for (int i = 0; i < corpus.size(); i++) {
            if (corpus.get(i).getReadName().equals(name)) {
                return i;
            }
        }
        throw new IllegalStateException("no corpus record named " + name);
    }

    static GATKRead copy(final SAMRecord record) {
        return new SAMRecordToGATKReadAdapter(record.deepCopy());
    }

    static void compare(final int i, final int j) {
        final ReadCoordinateComparator comparator = new ReadCoordinateComparator(header);
        String result;
        try {
            result = Integer.toString(Integer.signum(
                    comparator.compare(copy(corpus.get(i)), copy(corpus.get(j)))));
        } catch (final Exception | AssertionError e) {
            result = "E:" + e.getClass().getName();
        }
        System.out.printf("cmp\t%d\t%d\t%s%n", i, j, result);
    }

    static String describe(final AssemblyRegion region) {
        return String.format("%s|%s|%b|%d", region.getSpan(), region.getPaddedSpan(),
                region.isActive(), region.size());
    }

    static void ctor(final String label, final String contig, final int start, final int end,
                     final int padding) {
        try {
            final AssemblyRegion region = new AssemblyRegion(
                    new SimpleInterval(contig, start, end), true, padding, header);
            System.out.printf("ctor\t%s\t%s%n", label, describe(region));
        } catch (final Exception | AssertionError e) {
            System.out.printf("ctor\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static void ctorPair(final String label, final String contig, final int start, final int end,
                         final String paddedContig, final int paddedStart, final int paddedEnd) {
        try {
            final AssemblyRegion region = new AssemblyRegion(
                    new SimpleInterval(contig, start, end),
                    new SimpleInterval(paddedContig, paddedStart, paddedEnd), false, header);
            System.out.printf("ctor\t%s\t%s%n", label, describe(region));
        } catch (final Exception | AssertionError e) {
            System.out.printf("ctor\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /** Add every corpus record in declaration order, reporting each outcome. */
    static void addAllInCorpusOrder(final String label) {
        final AssemblyRegion region = new AssemblyRegion(
                new SimpleInterval("chr1", 100, 1900), true, 100, header);
        for (int i = 0; i < corpus.size(); i++) {
            System.out.printf("add\t%s\t%d\t%s%n", label, i, outcome(region, corpus.get(i)));
        }
        System.out.printf("region\t%s\t%s%n", label, describe(region));
        printReads(label, region);
    }

    static void addSorted(final String label) {
        final AssemblyRegion region = new AssemblyRegion(
                new SimpleInterval("chr1", 100, 1900), true, 100, header);
        final List<SAMRecord> sorted = new ArrayList<>(corpus);
        sorted.sort((a, b) -> new ReadCoordinateComparator(header)
                .compare(new SAMRecordToGATKReadAdapter(a), new SAMRecordToGATKReadAdapter(b)));
        for (int i = 0; i < sorted.size(); i++) {
            System.out.printf("add\t%s\t%d\t%s%n", label, i, outcome(region, sorted.get(i)));
        }
        System.out.printf("region\t%s\t%s%n", label, describe(region));
        printReads(label, region);
    }

    static void addOne(final String label, final int index) {
        final AssemblyRegion region = new AssemblyRegion(
                new SimpleInterval("chr1", 100, 1900), true, 100, header);
        System.out.printf("add\t%s\t0\t%s%n", label, outcome(region, corpus.get(index)));
        System.out.printf("region\t%s\t%s%n", label, describe(region));
    }

    static String outcome(final AssemblyRegion region, final SAMRecord record) {
        try {
            region.add(copy(record));
            return "ok";
        } catch (final Exception | AssertionError e) {
            return "E:" + e.getClass().getName() + ":" + oneLine(e.getMessage());
        }
    }

    /** The named regions that trim() is then applied to, kept so each trim starts from the same one. */
    static final java.util.Map<String, SimpleInterval[]> BUILT = new java.util.LinkedHashMap<>();
    static final java.util.Map<String, Integer> PADDING = new java.util.LinkedHashMap<>();

    static void region(final String label, final String contig, final int start, final int end,
                       final int padding) {
        BUILT.put(label, new SimpleInterval[] {new SimpleInterval(contig, start, end)});
        PADDING.put(label, padding);
        final AssemblyRegion region = build(label);
        System.out.printf("region\t%s\t%s%n", label, describe(region));
        printReads(label, region);
    }

    /**
     * Rebuild the named region from scratch: trim() returns a new region, but the reads it clips
     * come from this one, so every trim case must start from an untouched copy.
     */
    static AssemblyRegion build(final String label) {
        final AssemblyRegion region = new AssemblyRegion(
                BUILT.get(label)[0], true, PADDING.get(label), header);
        final List<GATKRead> reads = new ArrayList<>();
        for (final SAMRecord record : corpus) {
            final GATKRead read = copy(record);
            if (read.isUnmapped() || read.getContig() == null || read.getStart() < 1) {
                continue;
            }
            if (!region.getPaddedSpan().overlaps(read)) {
                continue;
            }
            reads.add(read);
        }
        reads.sort(new ReadCoordinateComparator(header));
        region.addAll(reads);
        return region;
    }

    static void printReads(final String label, final AssemblyRegion region) {
        final List<GATKRead> reads = region.getReads();
        for (int i = 0; i < reads.size(); i++) {
            final GATKRead read = reads.get(i);
            System.out.printf("read\t%s\t%d\t%s|%d|%s|%s%n", label, i, read.getName(),
                    read.getStart(), read.getCigar().toString(),
                    new String(read.getBasesNoCopy()));
        }
    }

    static void trim(final String source, final String label, final String contig, final int start,
                     final int end, final int padding) {
        try {
            final AssemblyRegion trimmed = build(source)
                    .trim(new SimpleInterval(contig, start, end), padding);
            System.out.printf("trim\t%s\t%s%n", label, describe(trimmed));
            printReads(label, trimmed);
        } catch (final Exception | AssertionError e) {
            System.out.printf("trim\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static void trimPair(final String source, final String label, final String contig,
                         final int start, final int end, final String paddedContig,
                         final int paddedStart, final int paddedEnd) {
        try {
            final AssemblyRegion trimmed = build(source).trim(
                    new SimpleInterval(contig, start, end),
                    new SimpleInterval(paddedContig, paddedStart, paddedEnd));
            System.out.printf("trim\t%s\t%s%n", label, describe(trimmed));
            printReads(label, trimmed);
        } catch (final Exception | AssertionError e) {
            System.out.printf("trim\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
