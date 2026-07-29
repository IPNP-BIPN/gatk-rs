/*
 * What the sample partition and the state bookkeeping do, taken from the reference.
 *
 * Between the per-read state machines and the pileup sits ReadStateManager, which decides which
 * reads enter, under which sample, and when they leave. Four of its behaviours change what a locus
 * contains, and none is visible from a tool's output:
 *
 *   - the sample map is a LinkedHashMap, and the reference says so in capitals: iteration is in the
 *     order the samples were given at construction, not sorted and not the header's order;
 *   - a read whose first stepForwardOnGenome() returns null is dropped in silence. A read that is
 *     all insertions and soft clips therefore never reaches any pileup, and upstream marks this a
 *     todo rather than an error;
 *   - collectPendingReads takes the left-most *genome position* among the states already in the
 *     system, which for a read part-way through a deletion is not where that read started, and
 *     admits only reads whose start equals it exactly on the same contig;
 *   - a read with no read group has the null sample, and a read whose sample was not declared at
 *     construction is a hard error rather than a new bucket.
 *
 * The classes are package-private, so the probe drives them through reflection rather than moving
 * this file into their package.
 *
 * Output:
 *
 *     step\t<label>\t<n>\t<total states>\t<sample>=<names at that sample>|...
 *     admit\t<label>\t<n>\t<reads admitted>
 *     error\t<label>\t<class>
 *
 * Usage: ReadStateDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class ReadStateDump {

    static final int CONTIG_LENGTH = 300;

    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                        new SAMSequenceRecord("chr2", CONTIG_LENGTH))));
        header.addReadGroup(readGroup("rg1", "sampleA"));
        header.addReadGroup(readGroup("rg2", "sampleB"));

        System.out.println("# ReadStateDump: the sample partition and the state advance");

        // Two samples, reads starting together and apart. The declared order is B then A, so the
        // LinkedHashMap ordering is measurable against the header's order, which is A then B.
        probe(header, "two-samples", List.of("sampleB", "sampleA"), new String[][] {
            {"a1", "rg1", "10M", "101"},
            {"b1", "rg2", "10M", "101"},
            {"a2", "rg1", "10M", "103"},
            {"b2", "rg2", "10M", "105"},
        }, 12);

        // A read that is all soft clips and insertions: its first step returns null and it never
        // enters, so the totals never count it.
        probe(header, "all-clipped", List.of("sampleA"), new String[][] {
            {"c1", "rg1", "10M", "101"},
            {"c2", "rg1", "5S5I", "101"},
            {"c3", "rg1", "10M", "101"},
        }, 12);

        // A read part-way through a deletion sets the boundary at its current genome position, not
        // at its start, so the read starting there is admitted mid-traversal.
        probe(header, "deletion-boundary", List.of("sampleA"), new String[][] {
            {"d1", "rg1", "2M10D8M", "101"},
            {"d2", "rg1", "10M", "105"},
            {"d3", "rg1", "10M", "113"},
        }, 20);

        // A read whose sample was never declared: a hard error rather than a new bucket.
        probe(header, "undeclared-sample", List.of("sampleA"), new String[][] {
            {"f1", "rg1", "10M", "101"},
            {"f2", "rg2", "10M", "101"},
        }, 5);

        // A read with no read group at all, declared as the null sample.
        probeWithNullSample(header);
    }

    static void probe(final SAMFileHeader header, final String label, final List<String> samples,
                      final String[][] specs, final int steps) throws Exception {
        final List<GATKRead> reads = new ArrayList<>();
        for (final String[] spec : specs) {
            reads.add(makeRead(header, spec[0], spec[1], spec[2], Integer.parseInt(spec[3])));
        }
        run(header, label, samples, reads, steps);
    }

    static void probeWithNullSample(final SAMFileHeader header) throws Exception {
        final List<GATKRead> reads = new ArrayList<>();
        reads.add(makeRead(header, "g1", "rg1", "10M", 101));
        reads.add(makeRead(header, "g2", null, "10M", 101));
        final List<String> samples = new ArrayList<>();
        samples.add("sampleA");
        samples.add(null);
        run(header, "null-sample", samples, reads, 12);
    }

    @SuppressWarnings("unchecked")
    static void run(final SAMFileHeader header, final String label, final List<String> samples,
                    final List<GATKRead> reads, final int steps) throws Exception {
        final Class<?> managerClass = Class.forName(
                "org.broadinstitute.hellbender.utils.locusiterator.ReadStateManager");
        final Class<?> infoClass = Class.forName(
                "org.broadinstitute.hellbender.utils.locusiterator.LIBSDownsamplingInfo");
        final Object noDownsampling = infoClass
                .getConstructor(boolean.class, int.class).newInstance(false, -1);

        final Constructor<?> constructor = managerClass.getDeclaredConstructor(
                Iterator.class, java.util.Collection.class, infoClass, SAMFileHeader.class);
        constructor.setAccessible(true);

        Object manager;
        try {
            manager = constructor.newInstance(reads.iterator(), samples, noDownsampling, header);
        } catch (final java.lang.reflect.InvocationTargetException e) {
            System.out.printf("error\t%s\t%s%n", label, e.getCause().getClass().getName());
            return;
        }

        final Method collect = managerClass.getDeclaredMethod("collectPendingReads");
        final Method update = managerClass.getDeclaredMethod("updateReadStates");
        final Method size = managerClass.getDeclaredMethod("size");
        final Method iterator = managerClass.getDeclaredMethod("iterator");
        collect.setAccessible(true);
        update.setAccessible(true);
        size.setAccessible(true);
        iterator.setAccessible(true);

        for (int step = 0; step < steps; step++) {
            try {
                collect.invoke(manager);
            } catch (final java.lang.reflect.InvocationTargetException e) {
                System.out.printf("error\t%s\t%s%n", label, e.getCause().getClass().getName());
                return;
            }

            final int total = (int) size.invoke(manager);
            final StringBuilder text = new StringBuilder();
            final Iterator<Map.Entry<String, Object>> entries =
                    (Iterator<Map.Entry<String, Object>>) iterator.invoke(manager);
            while (entries.hasNext()) {
                final Map.Entry<String, Object> entry = entries.next();
                if (text.length() > 0) {
                    text.append('|');
                }
                text.append(entry.getKey()).append('=').append(namesOf(entry.getValue()));
            }
            System.out.printf("step\t%s\t%d\t%d\t%s%n", label, step, total, text);

            update.invoke(manager);
        }
    }

    /** The read names currently held by one PerSampleReadStateManager, in its own order. */
    @SuppressWarnings("unchecked")
    static String namesOf(final Object perSample) throws Exception {
        final Method iterator = perSample.getClass().getDeclaredMethod("iterator");
        iterator.setAccessible(true);
        final Iterator<Object> states = (Iterator<Object>) iterator.invoke(perSample);
        final StringBuilder text = new StringBuilder();
        while (states.hasNext()) {
            final Object state = states.next();
            final Method getRead = state.getClass().getMethod("getRead");
            final Method getGenomePosition = state.getClass().getMethod("getGenomePosition");
            final GATKRead read = (GATKRead) getRead.invoke(state);
            if (text.length() > 0) {
                text.append(',');
            }
            text.append(read.getName()).append('@').append(getGenomePosition.invoke(state));
        }
        return text.length() == 0 ? "-" : text.toString();
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
        if (group != null) {
            record.setAttribute("RG", group);
        }
        return new SAMRecordToGATKReadAdapter(record);
    }
}
