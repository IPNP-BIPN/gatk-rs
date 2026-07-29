/*
 * Every PileupElement a read produces, taken from the reference.
 *
 * AlignmentStateDump records where a read stops. This records what the caller sees there, which is
 * a different set of answers: the indel-aware predicates that every pileup annotation and every
 * caller consults, and the base and quality a deletion reports instead of the read's own.
 *
 * The predicates are worth measuring together rather than one at a time, because two families of
 * them disagree on purpose:
 *
 *   - isBeforeInsertion and friends go through getAdjacentOperator, which looks at exactly the
 *     next cigar element whatever it is;
 *   - isBeforeDeletionStart goes through getNearestOnGenomeCigarElement, which skips everything
 *     that is not M, =, X or D.
 *
 * So on 3M2S3I3M the element before the clip is not "before an insertion", while a deletion behind
 * a clip is still seen. A port that used one walk for both would agree on every simple cigar and
 * diverge exactly here, which is why those cigars are in the list.
 *
 * The BI and BD tags are set on some reads and not others, because getBaseInsertionQual falls back
 * to a flat Q45 rather than to nothing: a pileup reports an insertion quality for a read that never
 * carried one, and the fallback array is as long as the read's *qualities*, not its bases.
 *
 * Output:
 *
 *     el\t<cigar>\t<n>\t<offset>|<base>|<qual>|<insQual>|<delQual>|<isDeletion>|<beforeDel>|<afterDel>|
 *          <beforeIns>|<afterIns>|<beforeClip>|<afterClip>|<atStart>|<atEnd>|<indelLength>|
 *          <insertedBases>|<prevOnGenome>|<nextOnGenome>|<between prev>|<between next>|<usable>
 *     count\t<cigar>\t<number of elements>
 *     forOffset\t<cigar>\t<offset>\t<ok:offset|E:class>
 *
 * Usage: PileupElementDump
 */

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.List;

public class PileupElementDump {

    static final int CONTIG_LENGTH = 200;

    /**
     * The cigars probed. The first block repeats AlignmentStateDump's, so the two suites disagree
     * only where the element disagrees with the stop; the second exists for the two navigation
     * families, where an insertion or deletion hides behind a soft clip.
     */
    static final String[] CIGARS = {
        "10M",
        "5M3D5M",
        "5M3I5M",
        "5M3D3I5M",
        "5M3I3D5M",
        "5M10N5M",
        "3S7M",
        "7M3S",
        "3S4M3S",
        "3H7M",
        "7M3H",
        "3H2S5M2S3H",
        "5M3I",
        "5M3S",
        "3I5M",
        "5M3P5M",
        "5=5X",
        "2M10D2M",
        // The two families of navigation, deliberately made to disagree.
        "3M2S3I3M",
        "3M3I2S3M",
        "3M2S3D3M",
        "3M2P3D3M",
        "3M2P3I3M",
        // An insertion and a deletion back to back, so getNextIndelCigarElement has to choose.
        "3M2I2D3M",
        "3M2D2I3M",
    };

    /** Cigars whose reads carry BI and BD tags, so the Q45 fallback is measured against a real one. */
    static final java.util.Set<String> TAGGED = java.util.Set.of("10M", "5M3D5M", "3S7M");

    public static void main(final String[] args) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));

        System.out.println("# PileupElementDump: what each stop looks like to a caller");
        for (final String cigar : CIGARS) {
            dump(header, cigar);
        }
    }

    static void dump(final SAMFileHeader header, final String cigarText) {
        final GATKRead read = makeRead(header, cigarText);
        final List<PileupElement> elements = new ArrayList<>();
        try {
            final AlignmentStateMachine machine = new AlignmentStateMachine(read);
            while (machine.stepForwardOnGenome() != null) {
                elements.add(machine.makePileupElement());
            }
        } catch (final Exception e) {
            // Recorded by the count row below: a cigar the machine refuses produces no elements.
        }

        for (int i = 0; i < elements.size(); i++) {
            final PileupElement element = elements.get(i);
            System.out.printf("el\t%s\t%d\t%d|%c|%d|%d|%d|%s|%s|%s|%s|%s|%s|%s|%s|%s|%d|%s|%s|%s|%s|%s|%s%n",
                    cigarText,
                    i,
                    element.getOffset(),
                    (char) element.getBase(),
                    element.getQual(),
                    element.getBaseInsertionQual(),
                    element.getBaseDeletionQual(),
                    element.isDeletion(),
                    element.isBeforeDeletionStart(),
                    element.isAfterDeletionEnd(),
                    element.isBeforeInsertion(),
                    element.isAfterInsertion(),
                    element.isBeforeSoftClip(),
                    element.isAfterSoftClip(),
                    element.isNextToSoftClip(),
                    element.atStartOfCurrentCigar(),
                    element.atEndOfCurrentCigar(),
                    element.getLengthOfImmediatelyFollowingIndel(),
                    element.getBasesOfImmediatelyFollowingInsertion() == null
                            ? "null" : element.getBasesOfImmediatelyFollowingInsertion(),
                    describe(element.getPreviousOnGenomeCigarElement()),
                    describe(element.getNextOnGenomeCigarElement()),
                    describeAll(element.getBetweenPrevPosition()),
                    describeAll(element.getBetweenNextPosition()),
                    PileupElement.isUsableBaseForAnnotation(element));
        }
        System.out.printf("count\t%s\t%d%n", cigarText, elements.size());

        // createPileupForReadAndOffset over every read offset, including the ones the alignment
        // never visits: inside a soft clip and inside an insertion, where the reference throws.
        for (int offset = 0; offset < read.getLength(); offset++) {
            String outcome;
            try {
                final PileupElement element =
                        PileupElement.createPileupForReadAndOffset(read, offset);
                outcome = "ok:" + element.getOffset() + ":" + element.getCurrentCigarOffset();
            } catch (final Exception e) {
                outcome = "E:" + e.getClass().getName();
            }
            System.out.printf("forOffset\t%s\t%d\t%s%n", cigarText, offset, outcome);
        }
    }

    static String describe(final CigarElement element) {
        return element == null ? "null" : element.getLength() + element.getOperator().toString();
    }

    static String describeAll(final List<CigarElement> elements) {
        if (elements.isEmpty()) {
            return "-";
        }
        final StringBuilder text = new StringBuilder();
        for (final CigarElement element : elements) {
            text.append(describe(element));
        }
        return text.toString();
    }

    /** A read at chr1:101, with BI and BD set on the cigars listed in TAGGED. */
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
            // Varied on purpose: a flat quality would hide an off-by-one in the offset.
            quals[i] = (byte) (20 + (i % 11));
        }
        record.setReadBases(bases);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        final GATKRead read = new SAMRecordToGATKReadAdapter(record);
        if (TAGGED.contains(cigarText)) {
            final byte[] insertion = new byte[length];
            final byte[] deletion = new byte[length];
            for (int i = 0; i < length; i++) {
                insertion[i] = (byte) (10 + i);
                deletion[i] = (byte) (30 + i);
            }
            ReadUtils.setInsertionBaseQualities(read, insertion);
            ReadUtils.setDeletionBaseQualities(read, deletion);
        }
        return read;
    }
}
