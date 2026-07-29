/*
 * How AlignmentContextIteratorBuilder routes, and in what order it hands the samples over.
 *
 * Two facts here reach into every locus-based tool's output, and neither is visible from the
 * builder's setters:
 *
 *   - the samples are collected with `Collectors.toSet()`, which is a HashSet. LocusIteratorByState
 *     creates its per-sample managers in *that* iteration order, and concatenates their elements in
 *     the same order to build each pileup. So the order of elements in a multi-sample pileup is the
 *     bucket order of a Java HashSet over the sample names: deterministic, but neither sorted nor
 *     the header's order. A port using a sorted map, an insertion-ordered list or a different hash
 *     would agree on single-sample data and diverge as soon as a second sample appeared;
 *
 *   - `areIntervalsSpecified` is `intervals != null`, not `!intervals.isEmpty()`. So an *empty but
 *     non-null* interval list counts as specified, routes to IntervalOverlappingIterator, and that
 *     constructor's `Utils.nonEmpty` throws. An empty interval list is an error rather than an
 *     empty traversal, and only when emitEmptyLoci is off.
 *
 * The routing is probed across (emitEmptyLoci, intervals null / empty / present), and the sample
 * order across name sets chosen to collide, to resize the table past its 16 buckets, and to include
 * a name whose hash is negative.
 *
 * Output:
 *
 *     order\t<label>\t<names in HashSet iteration order>
 *     hash\t<name>\t<String.hashCode>
 *     route\t<label>\t<class of the returned iterator|E:class>
 *     ctx\t<label>\t<n>\t<contig>:<pos>\t<size>\t<read names>
 *     count\t<label>\t<number of contexts>
 *
 * Usage: AlignmentContextBuilderDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentContextIteratorBuilder;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class AlignmentContextBuilderDump {

    static final int CONTIG_LENGTH = 200;

    public static void main(final String[] args) {
        System.out.println("# AlignmentContextBuilderDump: routing and sample order");

        // The HashSet iteration order, over name sets chosen for their edges.
        order("two", "sampleA", "sampleB");
        order("two-reversed", "sampleB", "sampleA");
        order("three", "NA12878", "NA12891", "NA12892");
        order("digits", "1", "2", "3", "10", "11");
        // Thirteen names, which crosses the 0.75 load factor of a 16-bucket table and forces the
        // order-preserving resize.
        order("thirteen", "s01", "s02", "s03", "s04", "s05", "s06", "s07", "s08", "s09", "s10",
                "s11", "s12", "s13");
        // A name whose String.hashCode is negative, so the sign handling in the index is measured.
        order("negative", "zzzzzzzzzzzz", "sampleA");
        order("empty-name", "", "sampleA");

        // The routing, and what each route yields.
        final SAMFileHeader header = header();
        final List<GATKRead> reads = List.of(
                makeRead(header, "a1", "rg1", "10M", 101),
                makeRead(header, "b1", "rg2", "10M", 101),
                makeRead(header, "a2", "rg1", "10M", 120));

        route(header, reads, "noloci-nullintervals", false, null);
        route(header, reads, "noloci-emptyintervals", false, List.of());
        route(header, reads, "noloci-intervals", false,
                List.of(new SimpleInterval("chr1", 105, 108)));
        route(header, reads, "emptyloci-nullintervals", true, null);
        route(header, reads, "emptyloci-emptyintervals", true, List.of());
        route(header, reads, "emptyloci-intervals", true,
                List.of(new SimpleInterval("chr1", 105, 112)));
        // Two intervals with a gap, so the empty loci between them are the ones that appear.
        route(header, reads, "emptyloci-twointervals", true,
                List.of(new SimpleInterval("chr1", 105, 107),
                        new SimpleInterval("chr1", 115, 117)));
    }

    /** The order a HashSet built exactly as the builder builds it iterates in. */
    static void order(final String label, final String... names) {
        final Set<String> samples = Stream.of(names).collect(Collectors.toSet());
        final StringBuilder text = new StringBuilder();
        for (final String sample : samples) {
            if (text.length() > 0) {
                text.append('|');
            }
            text.append(sample.isEmpty() ? "<empty>" : sample);
        }
        System.out.printf("order\t%s\t%s%n", label, text);
        for (final String name : names) {
            System.out.printf("hash\t%s\t%d%n", name.isEmpty() ? "<empty>" : name, name.hashCode());
        }
    }

    static void route(final SAMFileHeader header, final List<GATKRead> reads, final String label,
                      final boolean emitEmptyLoci, final List<SimpleInterval> intervals) {
        final AlignmentContextIteratorBuilder builder = new AlignmentContextIteratorBuilder();
        builder.setEmitEmptyLoci(emitEmptyLoci);

        Iterator<AlignmentContext> iterator;
        try {
            iterator = builder.build(reads.iterator(), header, intervals,
                    header.getSequenceDictionary(), true);
        } catch (final Exception e) {
            System.out.printf("route\t%s\tE:%s%n", label, e.getClass().getName());
            return;
        }
        System.out.printf("route\t%s\t%s%n", label, iterator.getClass().getSimpleName());

        int index = 0;
        while (iterator.hasNext()) {
            final AlignmentContext context = iterator.next();
            final StringBuilder names = new StringBuilder();
            for (final PileupElement element : context.getBasePileup()) {
                if (names.length() > 0) {
                    names.append(',');
                }
                names.append(element.getRead().getName());
            }
            System.out.printf("ctx\t%s\t%d\t%s:%d\t%d\t%s%n",
                    label, index, context.getContig(), context.getPosition(), context.size(),
                    names.length() == 0 ? "-" : names.toString());
            index++;
            if (index > 60) {
                // A whole-reference traversal with emitEmptyLoci would print the contig; the cap
                // is recorded rather than silent, so a truncated row cannot read as an ending.
                System.out.printf("truncated\t%s\t%d%n", label, index);
                break;
            }
        }
        System.out.printf("count\t%s\t%d%n", label, index);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.addReadGroup(readGroup("rg1", "sampleA"));
        header.addReadGroup(readGroup("rg2", "sampleB"));
        return header;
    }

    static SAMReadGroupRecord readGroup(final String id, final String sample) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        group.setSample(sample);
        group.setPlatform("ILLUMINA");
        return group;
    }

    static GATKRead makeRead(final SAMFileHeader header, final String name, final String group,
                             final String cigar, final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigar(TextCigarCodec.decode(cigar));
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
        record.setAttribute("RG", group);
        return new SAMRecordToGATKReadAdapter(record);
    }
}
