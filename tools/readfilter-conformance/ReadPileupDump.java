/*
 * What one locus looks like to a tool, taken from the reference.
 *
 * AlignmentStateDump records where each read stops; PileupElementDump records what each stop looks
 * like; this records what the collection of them answers, which is what a caller actually consults.
 *
 * Four of those answers are decisions rather than aggregation:
 *
 *   - getBaseCounts counts `*` as an A. BaseUtils.baseIndexMap maps the wildcard character to
 *     Base.A.ordinal(), with a comment saying so, so a read carrying `*` inflates the A count
 *     instead of being skipped the way an N is;
 *   - it skips deletions by an explicit test and skips everything else that maps to -1 by the
 *     index check. Two exclusions, two reasons, and a port that merged them would still agree;
 *   - splitBySample throws when a read has no sample and unknownSampleName is null, naming the
 *     first such read;
 *   - fixPairOverlappingQualities truncates twice and differently. Agreeing bases sum into a byte
 *     and are then capped, so the cap has to test for a *negative* stored value as well as for one
 *     over 93; disagreeing bases multiply the winner by 0.8 in double arithmetic and cast back.
 *
 * The overlap fix is probed as a matrix of quality pairs rather than through a real fragment,
 * because FragmentCollection is a separate port and the arithmetic is what is measured here.
 *
 * Output:
 *
 *     pileup\t<label>\t<size>\t<bases>\t<quals>\t<A,C,G,T>\t<offsets>\t<pileup string>
 *     sorted\t<label>\t<read names in sortedIterator order>
 *     sample\t<label>\t<sample>=<size>...
 *     split\t<label>\t<ok|E:class>
 *     overlap\t<first base><second base>\t<q1>\t<q2>\t<new q1>,<new q2>
 *
 * Usage: ReadPileupDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ReadPileupDump {

    static final int CONTIG_LENGTH = 200;
    static final int LOCUS = 105;

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        // Two samples plus one read group that declares none, so getSampleName can return null for
        // a read that does have a read group.
        header.addReadGroup(readGroup("rg1", "sampleA"));
        header.addReadGroup(readGroup("rg2", "sampleB"));
        header.addReadGroup(readGroup("rg3", null));

        System.out.println("# ReadPileupDump: what one locus answers");

        // A pileup with every base, a deletion, an N and a wildcard.
        probe(header, "mixed", new String[][] {
            {"r1", "rg1", "10M", "101", "ACGTACGTAC"},
            {"r2", "rg1", "10M", "101", "CCCCCCCCCC"},
            {"r3", "rg2", "2M5D3M", "101", "ACGTA"},
            {"r4", "rg2", "10M", "101", "NNNNNNNNNN"},
            {"r5", "rg1", "10M", "101", "****:*****"},
        });

        // Reads starting at different positions, which is the only key sortedIterator uses.
        probe(header, "staggered", new String[][] {
            {"s3", "rg1", "10M", "100", "ACGTACGTAC"},
            {"s1", "rg1", "10M", "96", "ACGTACGTAC"},
            {"s2", "rg2", "10M", "100", "TTTTTTTTTT"},
            {"s4", "rg1", "10M", "98", "GGGGGGGGGG"},
        });

        // A read whose group declares no sample, which splitBySample refuses without a fallback.
        probe(header, "nosample", new String[][] {
            {"n1", "rg1", "10M", "101", "ACGTACGTAC"},
            {"n2", "rg3", "10M", "101", "TTTTTTTTTT"},
        });

        // The overlap fix, as arithmetic. The pairs cross the cap, the truncation and the tie.
        final int[][] pairs = {
            {30, 30}, {30, 31}, {31, 30}, {1, 1}, {0, 0}, {93, 93},
            {50, 60}, {60, 50}, {80, 80}, {40, 0}, {0, 40}, {2, 3},
        };
        for (final int[] pair : pairs) {
            probeOverlap(header, true, pair[0], pair[1]);
            probeOverlap(header, false, pair[0], pair[1]);
        }
    }

    static void probe(final SAMFileHeader header, final String label, final String[][] specs) {
        final List<GATKRead> reads = new ArrayList<>();
        for (final String[] spec : specs) {
            reads.add(makeRead(header, spec[0], spec[1], spec[2], Integer.parseInt(spec[3]),
                    spec[4]));
        }
        final ReadPileup pileup = new ReadPileup(new SimpleInterval("chr1", LOCUS, LOCUS), reads);

        final int[] counts = pileup.getBaseCounts();
        System.out.printf("pileup\t%s\t%d\t%s\t%s\t%d,%d,%d,%d\t%s\t%s%n",
                label,
                pileup.size(),
                new String(pileup.getBases()),
                join(pileup.getQuals()),
                counts[0], counts[1], counts[2], counts[3],
                pileup.getOffsets().toString().replace(" ", ""),
                pileup.getPileupString('A'));

        final StringBuilder sorted = new StringBuilder();
        pileup.sortedIterator().forEachRemaining(e -> {
            if (sorted.length() > 0) {
                sorted.append('|');
            }
            sorted.append(e.getRead().getName()).append('@').append(e.getRead().getStart());
        });
        System.out.printf("sorted\t%s\t%s%n", label, sorted);

        final StringBuilder samples = new StringBuilder();
        for (final String sample : new java.util.TreeSet<>(
                java.util.Objects.requireNonNullElse(namesOf(pileup, header), List.<String>of()))) {
            if (samples.length() > 0) {
                samples.append('|');
            }
            samples.append(sample).append('=')
                   .append(pileup.getPileupForSample("null".equals(sample) ? null : sample, header)
                           .size());
        }
        System.out.printf("sample\t%s\t%s%n", label, samples);

        String outcome;
        try {
            final Map<String, ReadPileup> split = pileup.splitBySample(header, null);
            outcome = "ok:" + split.size();
        } catch (final Exception e) {
            outcome = "E:" + e.getClass().getName();
        }
        System.out.printf("split\t%s\t%s%n", label, outcome);
    }

    /** The sample names present, with null rendered as the string "null" so it can be sorted. */
    static List<String> namesOf(final ReadPileup pileup, final SAMFileHeader header) {
        final List<String> names = new ArrayList<>();
        for (final String sample : pileup.getSamples(header)) {
            names.add(sample == null ? "null" : sample);
        }
        return names;
    }

    /** fixPairOverlappingQualities over one pair of qualities, agreeing or not. */
    static void probeOverlap(final SAMFileHeader header, final boolean sameBase,
                             final int firstQual, final int secondQual) {
        final GATKRead first = makeRead(header, "o1", "rg1", "10M", 101, "AAAAAAAAAA");
        final GATKRead second = makeRead(header, "o2", "rg1", "10M", 101,
                sameBase ? "AAAAAAAAAA" : "CCCCCCCCCC");
        final byte[] firstQuals = first.getBaseQualities();
        final byte[] secondQuals = second.getBaseQualities();
        firstQuals[4] = (byte) firstQual;
        secondQuals[4] = (byte) secondQual;
        first.setBaseQualities(firstQuals);
        second.setBaseQualities(secondQuals);

        final PileupElement firstElement = PileupElement.createPileupForReadAndOffset(first, 4);
        final PileupElement secondElement = PileupElement.createPileupForReadAndOffset(second, 4);
        // fixPairOverlappingQualities is package-private and @VisibleForTesting, and this dump
        // lives in the default package, so it is reached by reflection rather than by moving the
        // dump into htsjdk's package tree.
        try {
            final java.lang.reflect.Method method = ReadPileup.class.getDeclaredMethod(
                    "fixPairOverlappingQualities", PileupElement.class, PileupElement.class);
            method.setAccessible(true);
            method.invoke(null, firstElement, secondElement);
        } catch (final ReflectiveOperationException e) {
            throw new IllegalStateException("cannot reach fixPairOverlappingQualities", e);
        }

        System.out.printf("overlap\t%s\t%d\t%d\t%d,%d%n",
                sameBase ? "same" : "differ",
                firstQual, secondQual,
                first.getBaseQualities()[4] & 0xff,
                second.getBaseQualities()[4] & 0xff);
    }

    static String join(final byte[] values) {
        final StringBuilder text = new StringBuilder();
        for (final byte value : values) {
            if (text.length() > 0) {
                text.append(',');
            }
            text.append(value & 0xff);
        }
        return text.toString();
    }

    static SAMReadGroupRecord readGroup(final String id, final String sample) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        if (sample != null) {
            group.setSample(sample);
        }
        group.setPlatform("ILLUMINA");
        return group;
    }

    static GATKRead makeRead(final SAMFileHeader header, final String name, final String group,
                             final String cigar, final int start, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigar(TextCigarCodec.decode(cigar));
        record.setReadBases(bases.getBytes());
        final byte[] quals = new byte[bases.length()];
        for (int i = 0; i < quals.length; i++) {
            quals[i] = (byte) (20 + i);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", group);
        return new SAMRecordToGATKReadAdapter(record);
    }
}
