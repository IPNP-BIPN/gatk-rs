/*
 * CigarUtils' clipping arithmetic and CigarBuilder's invariants, taken from the reference.
 *
 * Everything ReadClipper does passes through here: soft-clipping a tail, hard-clipping a region,
 * reverting soft clips. The cigar it produces is what the output BAM carries, so a rule missed
 * here is a wrong cigar on every clipped read.
 *
 * CigarBuilder is not a list: it merges consecutive identical operators, drops zero-length
 * elements, removes deletions that end up at either end, reorders a deletion that arrives after an
 * insertion, and refuses cigars whose sections are out of order. Each of those rewrites the cigar
 * a naive concatenation would have produced.
 *
 * Output:
 *
 *     clip\t<cigar>\t<start>\t<stopExclusive>\t<S|H>\t<result or E>
 *     shift\t<cigar>\t<numClipped>\t<refBasesClipped or E>
 *     revert\t<cigar>\t<result or E>
 *     build\t<elements, comma separated>\t<result or E>\t<leadingDelRemoved>\t<trailingDelRemoved>
 *
 * `E` is the reference throwing rather than answering: CigarBuilder validates, and a cigar that is
 * entirely soft-clipped or whose sections are out of order is a failure, not a value.
 *
 * Usage: CigarClipDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;

import java.util.function.Supplier;

public class CigarClipDump {

    /**
     * The cigars, chosen so each one reaches a different branch of the builder.
     *
     * `3M2D5M` puts a deletion where a left clip will strand it, `4M1I1D5M` puts a deletion
     * immediately after an insertion, and `2S2S6M` is two elements the builder merges into one.
     */
    static final String[] CIGARS = {
        "10M",
        "3S7M",
        "7M3S",
        "2S6M2S",
        "2H8M",
        "8M2H",
        "2H2S4M2S2H",
        "3M2D5M",
        "3M2I5M",
        "4M1I1D5M",
        "4M1D1I5M",
        "5M2N5M",
        "2S2S6M",
        "10S",
        "1M9S",
        "3M1D1I1D5M",
    };

    /** Element sequences fed to CigarBuilder directly, to reach its rewrites in isolation. */
    static final String[][] BUILDS = {
        {"3M", "4M"},                       // merged into 7M
        {"3M", "0M", "4M"},                 // zero-length dropped, then merged
        {"10S", "5D", "5M"},                // leading deletion after a clip: removed
        {"5M", "5D"},                       // trailing deletion: removed in make()
        {"5M", "5D", "10S"},                // trailing deletion before a right clip: removed
        {"5M", "2I", "3D", "5M"},           // deletion after insertion: reordered before it
        {"5M", "3D", "2I", "2D", "5M"},     // and merged into the preceding deletion
        {"10S", "2I", "5D", "5M"},          // deletion after a leading insertion: removed
        {"5M", "3D", "2I", "10S"},          // trailing deletion behind an insertion: removed
        {"10S"},                            // entirely soft-clipped: make() refuses
        {"5M", "3S", "5M"},                 // soft clip in the middle: section order refuses
        {"5H", "5M", "5H"},
    };

    public static void main(final String[] args) {
        System.out.println("# CigarClipDump: CigarUtils clipping and CigarBuilder, from the reference");

        for (final String text : CIGARS) {
            final Cigar cigar = TextCigarCodec.decode(text);
            final int readLength = cigar.getReadLength();

            // Left clips (start == 0) and right clips, which are the two shapes ReadClipper builds.
            for (final CigarOperator op : new CigarOperator[] {CigarOperator.S, CigarOperator.H}) {
                for (int stop = 1; stop <= readLength; stop++) {
                    emitClip(text, cigar, 0, stop, op);
                }
                for (int start = 0; start < readLength; start++) {
                    emitClip(text, cigar, start, readLength, op);
                }
            }

            for (int clipped = 0; clipped <= readLength; clipped++) {
                final int n = clipped;
                System.out.printf("shift\t%s\t%d\t%s%n",
                        text, n, call(() -> String.valueOf(CigarUtils.alignmentStartShift(cigar, n))));
            }

            System.out.printf("revert\t%s\t%s%n",
                    text, call(() -> CigarUtils.revertSoftClips(cigar).toString()));
        }

        for (final String[] elements : BUILDS) {
            final CigarBuilder builder = new CigarBuilder();
            final StringBuilder key = new StringBuilder();
            for (final String element : elements) {
                if (key.length() != 0) key.append(',');
                key.append(element);
            }
            final String result = call(() -> {
                for (final String element : elements) {
                    builder.add(one(element));
                }
                return builder.make().toString();
            });
            // The counters are read after make(), which is where the trailing removal happens.
            System.out.printf("build\t%s\t%s\t%s\t%s%n",
                    key, result,
                    call(() -> String.valueOf(builder.getLeadingDeletionBasesRemoved())),
                    call(() -> String.valueOf(builder.getTrailingDeletionBasesRemoved())));
        }
    }

    static void emitClip(final String text, final Cigar cigar, final int start, final int stop,
                         final CigarOperator op) {
        System.out.printf("clip\t%s\t%d\t%d\t%s\t%s%n", text, start, stop, op,
                call(() -> CigarUtils.clipCigar(cigar, start, stop, op).toString()));
    }

    /** One cigar element from its text, `3M` to (3, M). */
    static CigarElement one(final String text) {
        final int length = Integer.parseInt(text.substring(0, text.length() - 1));
        return new CigarElement(length, CigarOperator.characterToEnum(text.charAt(text.length() - 1)));
    }

    static String call(final Supplier<String> supplier) {
        try {
            return supplier.get();
        } catch (final Exception | AssertionError e) {
            return "E";
        }
    }
}
