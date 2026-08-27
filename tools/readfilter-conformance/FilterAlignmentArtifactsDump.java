/*
 * FilterAlignmentArtifacts' two decidable rules, taken from the reference.
 *
 * Whether a variant is an artefact of alignment: the reads supporting it are reassembled into
 * unitigs, the unitigs are realigned, and if they map somewhere else just as well the variant is
 * filtered. Two of the three steps are decidable without an aligner, and those are what this
 * measures: which reads are taken to SUPPORT the variant, and which realignments are taken to be
 * the SAME joint alignment.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A READ SUPPORTS A SNP BY ITS BASES, compared from the offset the reference coordinate maps
 *     to, and a variant position that falls inside a DELETION never supports a SNP;
 *   - A READ THAT DOES NOT REACH THE POSITION SUPPORTS NOTHING;
 *   - AN INDEL IS MATCHED BY CIGAR OPERATOR RATHER THAN BY BASES, because indel representation is
 *     not unique: a deletion is supported by a D or an S, an insertion by an I or an S;
 *   - THE TOLERANCE LOOP STOPS ADVANCING once an element is within tolerance, so every element
 *     after the first match is also treated as being at the variant, whatever its length;
 *   - --indel-start-tolerance WIDENS THAT WINDOW;
 *   - ONE UNITIG MAKES EVERY ALIGNMENT ITS OWN JOINT ALIGNMENT, and no unitigs make none;
 *   - TWO UNITIGS JOIN ONLY ON THE SAME STRAND, and only within half the maximum fragment length
 *     on each side;
 *   - THE BEST-SCORING ALIGNMENT OF EACH UNITIG IS THE ONE KEPT for a joint alignment;
 *   - AND THE RESULT COMES OUT OF A HashSet OF LISTS whose elements define no equals, so its order
 *     is identity-hash order and is not reproducible: the report below sorts it.
 *
 * Output:
 *
 *     support\t<label>=<true|false>
 *     joint\t<label>=<one line per joint alignment: refId:refStart-refEnd/score/mismatches, sorted>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FilterAlignmentArtifactsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.tools.walkers.realignmentfilter.RealignmentEngine;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FilterAlignmentArtifactsDump {

    static final int CONTIG_LENGTH = 199980;

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                new SAMSequenceRecord("chr2", CONTIG_LENGTH))));
        return header;
    }

    static GATKRead read(final SAMFileHeader header, final String name, final int start,
                         final String cigar, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        return new SAMRecordToGATKReadAdapter(record);
    }

    static VariantContext variant(final int start, final String reference, final String alternate) {
        final Allele ref = Allele.create(reference, true);
        final Allele alt = Allele.create(alternate, false);
        return new VariantContextBuilder("dump", "chr1", start,
                start + reference.length() - 1, List.of(ref, alt)).make();
    }

    static void support(final String label, final GATKRead read, final VariantContext variant,
                        final int tolerance) {
        System.out.printf("support\t%s=%s%n", label,
                RealignmentEngine.supportsVariant(read, variant, tolerance));
    }

    /** One alignment, with only the fields the joint rule reads. */
    static BwaMemAlignment alignment(final int refId, final int refStart, final int refEnd,
                                     final int score, final int mismatches, final boolean reverse) {
        return new BwaMemAlignment(reverse ? 16 : 0, refId, refStart, refEnd, 0, refEnd - refStart,
                60, mismatches, score, 0, null, null, null, -1, -1, 0);
    }

    static String describe(final BwaMemAlignment alignment) {
        return alignment.getRefId() + ":" + alignment.getRefStart() + "-" + alignment.getRefEnd()
                + "/" + alignment.getAlignerScore() + "/" + alignment.getNMismatches()
                + ((alignment.getSamFlag() & 16) != 0 ? "/-" : "/+");
    }

    static void joint(final String label, final List<List<BwaMemAlignment>> unitigs,
                      final int maxFragmentLength) {
        // The reference returns a HashSet's iteration order over lists whose elements define no
        // equals, which is identity-hash order and not reproducible: sorted here so the golden is.
        final List<String> lines = new ArrayList<>();
        for (final List<BwaMemAlignment> group :
                RealignmentEngine.findJointAlignments(copy(unitigs), maxFragmentLength)) {
            final List<String> parts = new ArrayList<>();
            for (final BwaMemAlignment alignment : group) {
                parts.add(describe(alignment));
            }
            java.util.Collections.sort(parts);
            lines.add(String.join(",", parts));
        }
        java.util.Collections.sort(lines);
        System.out.printf("joint\t%s=%s%n", label, String.join(";", lines));
    }

    /** findJointAlignments SORTS its argument in place, so each call gets its own copy. */
    static List<List<BwaMemAlignment>> copy(final List<List<BwaMemAlignment>> unitigs) {
        final List<List<BwaMemAlignment>> out = new ArrayList<>();
        for (final List<BwaMemAlignment> unitig : unitigs) {
            out.add(new ArrayList<>(unitig));
        }
        return out;
    }

    public static void main(final String[] args) {
        System.out.println("# FilterAlignmentArtifactsDump: which reads support a variant, and "
                + "which realignments are the same one");

        final SAMFileHeader header = header();

        // A SNP at 1005, on a read whose base there is the alternate.
        final VariantContext snp = variant(1005, "A", "C");
        support("snp-matching", read(header, "r", 1000, "20M", "ACGTACGTACGTACGTACGT"), snp, 0);
        // The same read against a SNP whose alternate it does not carry.
        support("snp-not-matching", read(header, "r", 1000, "20M", "ACGTACGTACGTACGTACGT"),
                variant(1005, "A", "G"), 0);
        // A read that stops before the variant.
        support("snp-not-covered", read(header, "r", 1000, "3M", "ACG"), snp, 0);
        // A read whose alignment DELETES the variant position, which never supports a SNP.
        support("snp-in-deletion", read(header, "r", 1000, "3M5D12M", "ACGTACGTACGTACG"), snp, 0);

        // An indel is matched by CIGAR OPERATOR rather than by bases. The walk compares a running
        // sum of element lengths against the variant's offset in the read, which is 5 here, so at
        // a tolerance of 0 the sum has to land on 5 EXACTLY.
        final VariantContext deletion = variant(1005, "ACGTA", "A");
        final VariantContext insertion = variant(1005, "A", "ACGTA");
        support("del-exact-deletion", read(header, "r", 1000, "5M4D15M",
                "ACGTACGTACGTACGTACGT"), deletion, 0);
        // A soft clip stands in for a deletion as well as for an insertion.
        support("del-exact-clip", read(header, "r", 1000, "5M15S",
                "ACGTACGTACGTACGTACGT"), deletion, 0);
        support("ins-exact-clip", read(header, "r", 1000, "5M15S",
                "ACGTACGTACGTACGTACGT"), insertion, 0);
        // An I operator cannot be reached at tolerance 0 at all: the read index for the variant
        // position lands PAST the inserted bases, so the running sum is always short of it by the
        // insertion's own length. Only the soft clip above supports an insertion exactly.
        support("ins-exact-insertion", read(header, "r", 1000, "5M4I15M",
                "ACGTACGTACGTACGTACGTACGT"), insertion, 0);
        // With the tolerance wide enough for the first element to qualify, the sum freezes at 0
        // and the I is reached after all.
        support("ins-tolerant-insertion", read(header, "r", 1000, "5M4I15M",
                "ACGTACGTACGTACGTACGTACGT"), insertion, 100);
        // A deletion's operator does not support an insertion, or the other way round.
        support("ins-exact-deletion", read(header, "r", 1000, "5M4D15M",
                "ACGTACGTACGTACGTACGT"), insertion, 0);
        support("del-exact-insertion", read(header, "r", 1000, "5M4I15M",
                "ACGTACGTACGTACGTACGTACGT"), deletion, 0);
        // A sum that steps OVER the offset never lands on it: 6M takes it from 0 to 6.
        support("del-stepped-over", read(header, "r", 1000, "6M4D14M",
                "ACGTACGTACGTACGTACGT"), deletion, 0);

        // The tolerance, and the quirk in the loop it exposes: `readPosition` is advanced ONLY in
        // the else branch, so as soon as one element is within tolerance the sum freezes and every
        // element after it is treated as being at the variant too. At a tolerance of 5 the very
        // first element already qualifies, the sum never moves off 0, and the whole cigar is
        // scanned for a supporting operator however far away it really is.
        final GATKRead far = read(header, "r", 1000, "12M4D8M", "ACGTACGTACGTACGTACGT");
        support("tolerance-0", far, deletion, 0);
        support("tolerance-5", far, deletion, 5);
        support("tolerance-100", far, deletion, 100);
        // The same read against an insertion, which its D cannot support at any tolerance.
        support("tolerance-100-insertion", far, insertion, 100);

        // Joint alignments. No unitigs at all.
        joint("no-unitigs", List.of(), 1000);
        // One unitig: every alignment becomes its own joint alignment.
        joint("one-unitig", List.of(List.of(
                alignment(0, 1000, 1100, 100, 0, false),
                alignment(0, 5000, 5100, 90, 2, false),
                alignment(1, 1000, 1100, 80, 4, false))), 1000);
        // Two unitigs that overlap on the same strand.
        joint("two-same-strand", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false))), 1000);
        // The same pair on opposite strands, which never joins.
        joint("two-opposite-strands", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, true))), 1000);
        // Far apart, so the padding decides it: half the maximum fragment length on each side.
        joint("far-apart-narrow", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 3000, 3100, 95, 1, false))), 1000);
        joint("far-apart-wide", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 3000, 3100, 95, 1, false))), 100000);
        // Two candidates for the second unitig at one locus: the higher score is kept.
        joint("best-score-kept", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false),
                        alignment(0, 1060, 1160, 120, 3, false))), 1000);
        // Two separate loci, each of which joins: two joint alignments.
        joint("two-loci", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false),
                        alignment(0, 50000, 50100, 98, 1, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false),
                        alignment(0, 50050, 50150, 93, 2, false))), 1000);
        // A locus only one of the two unitigs reaches, which is dropped.
        joint("one-sided", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false),
                        alignment(0, 90000, 90100, 99, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false))), 1000);
        // On different contigs, which never overlap.
        joint("different-contigs", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(1, 1000, 1100, 95, 1, false))), 1000);
        // Three unitigs, which must ALL overlap.
        joint("three-unitigs", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false)),
                List.of(alignment(0, 1080, 1180, 90, 2, false))), 1000);
        joint("three-unitigs-one-missing", List.of(
                List.of(alignment(0, 1000, 1100, 100, 0, false)),
                List.of(alignment(0, 1050, 1150, 95, 1, false)),
                List.of(alignment(0, 90000, 90100, 90, 2, false))), 1000);
    }
}
