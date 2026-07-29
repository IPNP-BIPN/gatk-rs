/*
 * Every pileup LocusIteratorByState yields, taken from the reference.
 *
 * This is what a LocusWalker iterates, so which elements land in which pileup is every
 * locus-based tool's output. Four decisions sit in the twenty lines that build each pileup:
 *
 *   - the two exclusions are not symmetric. A read whose current operator is N is skipped *before*
 *     the adaptor test; a read whose operator is D is skipped *inside* it. So a read sitting in an
 *     adaptor with a deletion is excluded once, by the adaptor, and a read with an N is excluded
 *     whatever the adaptor says. Reordering the two changes nothing on ordinary data and changes
 *     the pileup on exactly the reads carrying both;
 *   - the adaptor test is per base, so the same read contributes to some loci and not others;
 *   - the pileup is monolithic and concatenated in sample order, which the reference kept for the
 *     HaplotypeCaller's benefit and documents as such;
 *   - a locus with no surviving element yields no context at all, silently. Emitting empty loci is
 *     a different class's job.
 *
 * The four (includeDeletions, includeNs) combinations are all probed, because a tool picks them by
 * overriding two methods and the defaults differ: LocusWalker.includeDeletions() is true and
 * includeNs() is false.
 *
 * Output:
 *
 *     ctx\t<label>\t<n>\t<contig>:<pos>\t<size>\t<bases>\t<read names>
 *     count\t<label>\t<number of contexts>
 *
 * Usage: LocusIteratorDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;

public class LocusIteratorDump {

    static final int CONTIG_LENGTH = 400;

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.addReadGroup(readGroup("rg1", "sampleA"));
        header.addReadGroup(readGroup("rg2", "sampleB"));

        System.out.println("# LocusIteratorDump: every pileup LocusIteratorByState yields");

        // Plain coverage over two samples, staggered so the pileups grow and shrink.
        final String[][] plain = {
            {"a1", "rg1", "10M", "101", "0"},
            {"b1", "rg2", "10M", "104", "0"},
            {"a2", "rg1", "10M", "108", "0"},
        };
        // A deletion and a skip, so both exclusions have something to exclude.
        final String[][] indels = {
            {"d1", "rg1", "3M4D3M", "101", "0"},
            {"n1", "rg1", "3M4N3M", "101", "0"},
            {"m1", "rg2", "10M", "101", "0"},
        };
        // A read inside its adaptor: a short fragment on the reverse strand, so the adaptor
        // boundary falls inside the read and excludes its earlier bases only.
        final String[][] adaptor = {
            {"p1", "rg1", "10M", "101", "60"},
            {"q1", "rg1", "10M", "101", "0"},
        };

        for (final boolean deletions : new boolean[] {true, false}) {
            for (final boolean ns : new boolean[] {true, false}) {
                final String suffix = (deletions ? "D" : "d") + (ns ? "N" : "n");
                probe(header, "plain-" + suffix, plain, deletions, ns);
                probe(header, "indels-" + suffix, indels, deletions, ns);
                probe(header, "adaptor-" + suffix, adaptor, deletions, ns);
            }
        }
    }

    static void probe(final SAMFileHeader header, final String label, final String[][] specs,
                      final boolean includeDeletions, final boolean includeNs) {
        final List<GATKRead> reads = new ArrayList<>();
        for (final String[] spec : specs) {
            reads.add(makeRead(header, spec[0], spec[1], spec[2], Integer.parseInt(spec[3]),
                    Integer.parseInt(spec[4])));
        }
        // The sample order is the declaration order, which the pileup concatenation follows.
        final List<String> samples = List.of("sampleA", "sampleB");

        final LocusIteratorByState iterator = new LocusIteratorByState(
                reads.iterator(), LocusIteratorByState.NO_DOWNSAMPLING, samples, header,
                includeDeletions, includeNs);

        int index = 0;
        while (iterator.hasNext()) {
            final AlignmentContext context = iterator.next();
            final StringBuilder names = new StringBuilder();
            final StringBuilder bases = new StringBuilder();
            for (final PileupElement element : context.getBasePileup()) {
                if (names.length() > 0) {
                    names.append(',');
                }
                names.append(element.getRead().getName());
                bases.append((char) element.getBase());
            }
            System.out.printf("ctx\t%s\t%d\t%s:%d\t%d\t%s\t%s%n",
                    label, index, context.getContig(), context.getPosition(),
                    context.size(), bases, names);
            index++;
        }
        System.out.printf("count\t%s\t%d%n", label, index);
    }

    static SAMReadGroupRecord readGroup(final String id, final String sample) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        group.setSample(sample);
        group.setPlatform("ILLUMINA");
        return group;
    }

    /**
     * A read at chr1:start. A non-zero fragment length makes it a properly-paired reverse-strand
     * read with a mate ahead of it, which is what gives it a computable adaptor boundary.
     */
    static GATKRead makeRead(final SAMFileHeader header, final String name, final String group,
                             final String cigar, final int start, final int fragmentLength) {
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
        if (fragmentLength != 0) {
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
            record.setReadNegativeStrandFlag(true);
            record.setMateReferenceName("chr1");
            // The mate must start *inside* this read for the adaptor boundary to fall inside it:
            // for a reverse-strand read the boundary is mateStart - 1, so a mate before the read
            // puts the boundary before every locus and excludes nothing. The first version of this
            // fixture did exactly that and measured no exclusion at all.
            record.setMateAlignmentStart(start + 5);
            record.setMateNegativeStrandFlag(false);
            record.setInferredInsertSize(-fragmentLength);
        }
        return new SAMRecordToGATKReadAdapter(record);
    }
}
