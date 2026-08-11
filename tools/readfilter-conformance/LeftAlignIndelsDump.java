/*
 * AlignmentUtils.leftAlignIndels and AlignmentUtils.normalizeAlleles, taken from the reference.
 *
 * The call LeftAlignIndels is built around, measured on its own because its interface is exactly
 * measurable: a cigar, a reference, a read and a start go in; a cigar and two counts of removed
 * deletion bases come out. Nothing about it needs a file.
 *
 * Four behaviours this is built to catch.
 *
 *   - THE CIGAR IS WALKED RIGHT TO LEFT, accumulating an indel's ranges and resolving them only at
 *     an alignment block or the start of the cigar. Two indels with too few matching bases between
 *     them therefore COLLIDE AND MERGE into one, which is a behaviour of the loop rather than of a
 *     named rule;
 *   - normalizeAlleles CAN SHIFT RIGHT. It trims shared bases off the end and then off the front,
 *     and the front trim does startShift--. So the returned start shift is NEGATIVE on
 *     non-parsimonious input, and leftAlignIndels reads that as newMatchOnLeftDueToTrimming. A
 *     port that typed the shift as unsigned would be wrong in a way no simple case shows. Both
 *     functions are dumped, the inner one on its own inputs;
 *   - A DELETION THAT REACHES THE START OF THE CIGAR IS DROPPED rather than emitted, and the
 *     reference bases it removed are reported instead: that count is what makes the tool move the
 *     read;
 *   - TWO VALIDATIONS REFUSE RATHER THAN CLAMP, and both are IllegalArgumentException out of a
 *     util rather than a UserException out of a tool.
 *
 * Output:
 *
 *     leftalign\t<label>\t<cigar>\t<ref>\t<read>\t<readStart>\t<result cigar>\t<leading>\t<trailing>
 *     leftalignerror\t<label>\t<class>:<message>
 *     normalize\t<label>\t<ref>\t<read>\t<refRange>\t<readRange>\t<maxShift>\t<trim>\t<startShift>\t<endShift>\t<refRange after>\t<readRange after>
 *
 * Usage: LeftAlignIndelsDump
 */

import htsjdk.samtools.TextCigarCodec;

import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;

import java.util.Arrays;
import java.util.List;

public class LeftAlignIndelsDump {

    public static void main(final String[] args) {
        System.out.println("# LeftAlignIndelsDump: AlignmentUtils.leftAlignIndels");

        // An indel already at the leftmost position, and one that cannot move at all: the cigar
        // comes back as it went in, and that is worth a row rather than an assumption.
        leftAlign("already-left", "2M2D4M", "AATTAACC", "AAAACC", 0);
        leftAlign("no-repeat", "3M1I3M", "AAACCC", "AAATCCC", 0);

        // A deletion inside a homopolymer: every position of the deleted base is the same
        // haplotype, and the convention is the leftmost.
        leftAlign("homopolymer-deletion", "4M1D3M", "AAAAAAAT", "AAAAAAT", 0);
        // The same for an insertion.
        leftAlign("homopolymer-insertion", "4M1I3M", "AAAAAAT", "AAAAAAAT", 0);
        // A dinucleotide repeat, where the shift is two bases at a time.
        leftAlign("dinucleotide-deletion", "6M2D4M", "ATATATATATGG", "ATATATATGG", 0);

        // Soft clips are not alignment blocks, so an indel may not be shifted into them. The
        // reference's own javadoc example, which is the case that says soft-clipped bases are
        // present in the read's byte array.
        leftAlign("javadoc-soft-clip", "2S2M2I", "GGAA", "TTAAAA", 2);
        leftAlign("soft-clip-both-ends", "2S4M1D3M2S", "AAAAAAAT", "TTAAAAAATTT", 0);
        // A hard clip consumes neither, and is not an alignment block either.
        leftAlign("hard-clip", "2H4M1D3M", "AAAAAAAT", "AAAAAAT", 0);

        // Two indels with too few matching bases between them to shift past each other.
        leftAlign("colliding-indels", "3M1D2M1D3M", "AAAAAAAAAT", "AAAAAAAT", 0);
        leftAlign("insertion-then-deletion", "3M2I2M2D3M", "AAAAAAAAT", "AAAAAAAAAT", 0);

        // A deletion that walks all the way to the start of the cigar is dropped, and the bases it
        // removed are reported instead. This is the row the tool's read-moving line exists for.
        leftAlign("deletion-to-the-start", "1M2D5M", "AAAAAAAA", "AAAAAA", 0);

        // A read that starts partway into the reference, which is what the tool passes as 0 and a
        // haplotype caller passes as something else.
        leftAlign("read-start-offset", "4M1D3M", "GGGGAAAAAAAT", "AAAAAAT", 4);

        // The two refusals.
        leftAlign("past-the-reference", "4M1D3M", "AAAA", "AAAAAAT", 0);
        leftAlign("cigar-misses-read-bases", "4M1D3M", "AAAAAAAAAAAA", "AAAAAAAAAAAA", 0);

        // normalizeAlleles on its own, because the shift it returns is the number leftAlignIndels
        // branches on and one of these is negative.
        //
        // The reference's own example: reference GAAT, read GAAAT, the insertion initially placed
        // before the T.
        normalize("javadoc-parsimonious", "GAAT", "GAAAT", 3, 3, 3, 4, 3, true);
        // The same alleles left un-parsimonious by one shared base at each end, which is the case
        // that trims and returns a NEGATIVE start shift.
        normalize("not-parsimonious", "GAAT", "GAAAT", 2, 4, 2, 5, 2, true);
        // With trimming off, the same input does not trim and only shifts.
        normalize("not-parsimonious-untrimmed", "GAAT", "GAAAT", 2, 4, 2, 5, 2, false);
        // maxShift is the wall: the alleles could move further and are not allowed to.
        normalize("max-shift-of-one", "GAAAAT", "GAAAAAT", 5, 5, 5, 6, 1, true);
        normalize("max-shift-of-zero", "GAAAAT", "GAAAAAT", 5, 5, 5, 6, 0, true);
        // Nothing shared to the left, so nothing moves.
        normalize("blocked-immediately", "GCAAT", "GCAAAT", 4, 4, 4, 5, 4, true);
        // THE NEGATIVE ONE. The alleles end differently, so the end trim stops at once, but they
        // start with the same base, so the front trim runs and does startShift--. Nothing can then
        // shift left, because the last bases still differ. The function returns -1.
        normalize("shifts-right", "CAG", "CAT", 1, 3, 1, 3, 1, true);
    }

    /** One call, with the cigar the reference gives back and the two counts beside it. */
    static void leftAlign(final String label, final String cigar, final String ref,
                          final String read, final int readStart) {
        try {
            final CigarBuilder.Result result = AlignmentUtils.leftAlignIndels(
                    TextCigarCodec.decode(cigar), ref.getBytes(), read.getBytes(), readStart);
            System.out.printf("leftalign\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d%n", label, cigar, ref,
                    read, readStart, result.getCigar().toString(),
                    result.getLeadingDeletionBasesRemoved(),
                    result.getTrailingDeletionBasesRemoved());
        } catch (final Exception | AssertionError e) {
            // The refusal is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("leftalignerror\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
        }
    }

    /**
     * One call of the inner function, with the ranges before and after.
     *
     * The ranges are a side effect: normalizeAlleles adjusts what it is given, and the caller reads
     * them afterwards rather than from the return value. Both are dumped for that reason.
     */
    static void normalize(final String label, final String ref, final String read,
                          final int refFrom, final int refTo, final int readFrom, final int readTo,
                          final int maxShift, final boolean trim) {
        final IndexRange refRange = new IndexRange(refFrom, refTo);
        final IndexRange readRange = new IndexRange(readFrom, readTo);
        final String before = String.format("[%d,%d)\t[%d,%d)", refFrom, refTo, readFrom, readTo);
        try {
            final var shifts = AlignmentUtils.normalizeAlleles(
                    Arrays.asList(ref.getBytes(), read.getBytes()),
                    List.of(refRange, readRange), maxShift, trim);
            System.out.printf("normalize\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t[%d,%d)\t[%d,%d)%n",
                    label, ref, read, before, maxShift, trim, shifts.getLeft(), shifts.getRight(),
                    refRange.getStart(), refRange.getEnd(), readRange.getStart(),
                    readRange.getEnd());
        } catch (final Exception | AssertionError e) {
            System.out.printf("normalize\t%s\t%s\t%s\t%s\t%d\t%s\t%s:%s%n", label, ref, read,
                    before, maxShift, trim, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
        }
    }
}
