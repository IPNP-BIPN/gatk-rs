/*
 * SAMPileupCodec, taken from the reference.
 *
 * The samtools mpileup format, which CheckPileup reads as its truth. One line is a locus: contig,
 * position, reference base, coverage, the read bases and their qualities.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE BASES COLUMN IS A LITTLE LANGUAGE, not a list of bases. `.` and `,` are the reference
 *     base, `*` is a deletion, `$` ends a read and consumes NO quality, `^` starts one and eats the
 *     NEXT CHARACTER as a mapping quality, and `+`/`-` introduce an indel whose length is the number
 *     that follows and whose bases are skipped entirely;
 *   - SO THE QUALITY STRING IS CONSUMED AT A DIFFERENT RATE FROM THE BASES STRING, one quality per
 *     emitted element and none for the markers, which is what makes the two counters `i` and `j`
 *     independent;
 *   - A COVERAGE OF ZERO RETURNS EARLY and never looks at the two columns after it, so a line whose
 *     bases and qualities are nonsense parses as long as its coverage says zero;
 *   - A COVERAGE THAT DISAGREES WITH THE NUMBER OF PARSED ELEMENTS IS A CodecLineParsingException,
 *     and its message has a typo the port must carry: "THe SAM pileup line". The field-count message
 *     carries one too, "this codes is only valid";
 *   - LEFTOVER QUALITIES ARE AN ERROR AND LEFTOVER BASES ARE NOT REACHABLE: the loop runs to the end
 *     of the bases, so only `j != qualities.length()` can fail that check;
 *   - AN N IS FOLDED TO N AND EVERY OTHER BASE GOES THROUGH simpleBaseToBaseIndex, so a lower-case
 *     base comes out UPPER CASE and anything outside ACGTN is refused by name;
 *   - THE FIELD COUNT MUST BE 4 TO 6, WHICH IS NOT THE FORMAT'S. A samtools mpileup line has six
 *     columns, and a SEVEN-column line, which the format allows, is REFUSED here. Worse, a FIVE
 *     column line passes the count check and then dies on `tokens[5]` with a raw
 *     ArrayIndexOutOfBoundsException rather than a CodecLineParsingException;
 *   - A DELETION PLACEHOLDER BECOMES THE LETTER `D`, which is not a base of the alphabet the rest of
 *     the codec accepts: `BaseUtils.Base.D.base` goes in without passing through parseBase;
 *   - AND canDecode IS BY EXTENSION ONLY, `.pileup` or `.mpileup`, case-insensitively, with one
 *     block-compressed extension stripped first, so `x.mpileup.gz` and `x.PILEUP` decode while
 *     `x.pileup.txt` and a bare `pileup` do not.
 *
 * Output:
 *
 *     decode\t<label>\t<the line>\t<chr:pos ref cov>\t<bases>\t<quals as numbers>
 *     error\t<label>\t<exception class>:<message>
 *     candecode\t<path>\t<true|false>
 *
 * Usage: SAMPileupCodecDump
 */

import org.broadinstitute.hellbender.utils.codecs.sampileup.SAMPileupCodec;
import org.broadinstitute.hellbender.utils.codecs.sampileup.SAMPileupFeature;

import java.util.Arrays;
import java.util.stream.Collectors;

public class SAMPileupCodecDump {

    public static void main(final String[] args) {
        System.out.println("# SAMPileupCodecDump: the samtools mpileup format, from the reference");

        final String[][] lines = {
                // Plain matches to the reference, in both directions.
                {"matches", "chr1\t10\tA\t4\t.,.,\tIIII"},
                // Explicit bases, upper and lower case, which fold to upper case.
                {"explicit", "chr1\t10\tA\t4\tACgt\tIIII"},
                // A deletion placeholder, which is a base of its own.
                {"deletion", "chr1\t10\tA\t3\t.*.\tIII"},
                // A read start, whose mapping quality character is eaten and emits nothing.
                {"read-start", "chr1\t10\tA\t2\t^I.,\tII"},
                // A read end, which consumes no quality at all.
                {"read-end", "chr1\t10\tA\t2\t.$,\tII"},
                // An insertion of two bases, skipped whole.
                {"insertion", "chr1\t10\tA\t2\t.+2AC,\tII"},
                // A deletion of ten bases, whose length is two digits.
                {"long-indel", "chr1\t10\tA\t2\t.-10ACGTACGTAC,\tII"},
                // An N in the bases, folded to N.
                {"n-base", "chr1\t10\tA\t2\t.N\tII"},
                // Coverage zero, which returns before the two columns after it are looked at.
                {"zero-coverage", "chr1\t10\tA\t0\t*\t*"},
                {"zero-coverage-nonsense", "chr1\t10\tA\t0\tZZZZ\t!!!!"},
                // A seventh column, which the format allows and this codec refuses.
                {"seven-columns", "chr1\t10\tA\t2\t.,\tII\t~~"},
                // The refusals.
                {"five-columns", "chr1\t10\tA\t2\t.,"},
                {"eight-columns", "chr1\t10\tA\t2\t.,\tII\t~~\textra"},
                {"coverage-mismatch", "chr1\t10\tA\t5\t.,\tII"},
                {"bad-position", "chr1\tten\tA\t2\t.,\tII"},
                {"bad-reference", "chr1\t10\tZ\t2\t.,\tII"},
                {"bad-base", "chr1\t10\tA\t2\t.Z\tII"},
                {"indel-without-length", "chr1\t10\tA\t2\t.+AC,\tII"},
                {"too-few-qualities", "chr1\t10\tA\t2\t.,\tI"},
                {"too-many-qualities", "chr1\t10\tA\t2\t.,\tIII"},
        };

        final SAMPileupCodec codec = new SAMPileupCodec();
        for (final String[] pair : lines) {
            try {
                final SAMPileupFeature feature = codec.decode(pair[1]);
                final String quals = Arrays.toString(feature.getBaseQuals())
                        .replace(" ", "");
                System.out.printf("decode\t%s\t%s\t%s:%d %c %d\t%s\t%s%n", pair[0],
                        ReferenceQueryDump.escape(pair[1]),
                        feature.getContig(), feature.getStart(), (char) feature.getRef(),
                        feature.size(),
                        feature.getBasesString(), quals);
            } catch (final Exception e) {
                // The message embeds the line, tabs and all, so it is escaped like any other
                // field: a raw tab here would split one row into several columns.
                System.out.printf("error\t%s\t%s:%s%n", pair[0], e.getClass().getName(),
                        ReferenceQueryDump.escape(e.getMessage()));
            }
        }

        for (final String path : new String[] {
                "x.pileup", "x.mpileup", "x.mpileup.gz", "x.pileup.bgz", "x.txt", "x.PILEUP",
                "x.pileup.txt", "pileup"}) {
            System.out.printf("candecode\t%s\t%s%n", path, codec.canDecode(path));
        }
    }
}
