/*
 * `SequenceDictionaryUtils.compareDictionaries` and `validateDictionaries`: whether two sequence
 * dictionaries may be used together, and what the refusal says when they may not.
 *
 * Measured because a CountVariants covering-array row disagreed on it and nothing else: four rows
 * pass a BAM to --input alongside a VCF to --variant, and the reference refuses the pair before
 * any record is read, where the port counted (gatk-rs#1038).
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE OUTCOME IS ONE OF EIGHT NAMED CASES, and three of them are only ever reached when the
 *     caller asked for the ordering to be checked;
 *   - A LENGTH OF ZERO IS EQUIVALENT TO ANY LENGTH, so a dictionary that declares a contig without
 *     one agrees with a dictionary that does;
 *   - REVERSING THE COMMON CONTIGS IS A SUPERSET without the ordering check and OUT_OF_ORDER with
 *     it, so the same pair is accepted or refused by an argument;
 *   - THE SAME CONTIGS AT DIFFERENT ABSOLUTE INDICES are DIFFERENT_INDICES, which is a case of its
 *     own and not OUT_OF_ORDER: the relative order can be right while the positions are not;
 *   - A COMMON SUBSET IS ACCEPTED unless the caller required a superset, and the refusal then
 *     names the contigs that are MISSING;
 *   - THE HUMAN ORDER CHECK NEEDS CHR1, CHR2 AND CHR10 by name AND by length, and it fires on
 *     either dictionary, naming whichever one is lexicographic;
 *   - AND AN EMPTY DICTIONARY HAS NO COMMON CONTIGS, which is the same refusal as two dictionaries
 *     that simply disagree.
 *
 * Output:
 *
 *     compare\t<case>\t<ordering checked>\t<the enum constant>
 *     validate\t<case>\t<superset required>\t<ordering checked>\t<ok, or class: message>
 *
 * Usage: SequenceDictionaryValidationDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;

import java.util.ArrayList;
import java.util.List;

public class SequenceDictionaryValidationDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A dictionary from `name:length` pairs, in the order they are given. */
    static SAMSequenceDictionary dict(final String... records) {
        final List<SAMSequenceRecord> sequences = new ArrayList<>();
        for (final String record : records) {
            final int colon = record.lastIndexOf(':');
            sequences.add(new SAMSequenceRecord(record.substring(0, colon),
                    Integer.parseInt(record.substring(colon + 1))));
        }
        return new SAMSequenceDictionary(sequences);
    }

    /** One pair, compared and validated under every combination of the two flags. */
    static void run(final String name, final SAMSequenceDictionary a, final SAMSequenceDictionary b) {
        for (final boolean ordering : new boolean[] {false, true}) {
            emit("compare", name + "\t" + ordering,
                    SequenceDictionaryUtils.compareDictionaries(a, b, ordering).name());
        }
        for (final boolean superset : new boolean[] {false, true}) {
            for (final boolean ordering : new boolean[] {false, true}) {
                String answer;
                try {
                    SequenceDictionaryUtils.validateDictionaries("reads", a, "features", b,
                            superset, ordering);
                    answer = "ok";
                } catch (final Exception e) {
                    answer = e.getClass().getName() + ": " + e.getMessage();
                }
                emit("validate", name + "\t" + superset + "\t" + ordering, answer);
            }
        }
    }

    /** hg19's chr1, chr2 and chr10, whose LENGTHS are what the human check recognises. */
    static final String CHR1 = "chr1:249250621";
    static final String CHR2 = "chr2:243199373";
    static final String CHR10 = "chr10:135534747";

    public static void main(final String[] args) {
        run("identical", dict("chr1:100", "chr2:200"), dict("chr1:100", "chr2:200"));
        // The first dictionary holds everything the second does, and more.
        run("superset", dict("chr1:100", "chr2:200", "chr3:300"), dict("chr1:100", "chr2:200"));
        // And the other way round, which is a common subset rather than a superset.
        run("common-subset", dict("chr1:100"), dict("chr1:100", "chr2:200"));
        // The row that started this: nothing in common at all.
        run("no-common-contigs", dict("chr1:100"), dict("chrA:100"));
        // A name they share and a length they do not.
        run("unequal-lengths", dict("chr1:100"), dict("chr1:200"));
        // A length of zero, which is "unknown" and equivalent to anything.
        run("zero-length", dict("chr1:0"), dict("chr1:200"));
        run("zero-length-both-sides", dict("chr1:0"), dict("chr1:0"));
        // The same contigs in the other order: a superset without the ordering check.
        run("reversed", dict("chr1:100", "chr2:200"), dict("chr2:200", "chr1:100"));
        // The same relative order at different absolute positions.
        run("different-indices", dict("chrX:10", "chr1:100", "chr2:200"), dict("chr1:100", "chr2:200"));
        // An empty dictionary, which shares nothing with anything.
        run("empty-first", dict(), dict("chr1:100"));
        run("empty-both", dict(), dict());
        // The human check: chr1, chr2 and chr10 by name AND by length, out of karyotypic order.
        run("lexicographic-human-first",
                dict(CHR1, CHR10, CHR2), dict(CHR1, CHR2, CHR10));
        run("lexicographic-human-second",
                dict(CHR1, CHR2, CHR10), dict(CHR1, CHR10, CHR2));
        // The same names in the same wrong order with the WRONG lengths, which the check does not
        // recognise as human at all.
        run("lexicographic-not-human",
                dict("chr1:100", "chr10:300", "chr2:200"), dict("chr1:100", "chr2:200", "chr10:300"));

        System.out.print(buf);
    }
}
