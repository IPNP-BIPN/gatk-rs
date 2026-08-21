/*
 * CompareReferences' table and analysis, taken from the reference.
 *
 * Several references compared by the MD5 of each sequence. The table is keyed by MD5 and not by
 * name, which is the whole idea, and the pairwise analysis that follows is a set of flags removed
 * and added in an order that matters.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE TABLE IS KEYED BY MD5, so two references whose sequences agree base for base share one
 *     row however differently they name them, and a name appearing twice under different bases is
 *     two rows;
 *   - THE COLUMN NAME IS THE FILE'S NAME, not its path, so two references in different directories
 *     with the same file name would collide;
 *   - A MISSING ENTRY IS `---`, and the row still carries the length of whichever reference had it;
 *   - EVERY PAIR STARTS AS EXACT_MATCH and loses it on the first disagreement, so the flags are a
 *     removal followed by additions rather than a classification;
 *   - DIFFER_IN_SEQUENCE_NAMES NEEDS BOTH SIDES PRESENT: one row whose two entries differ and are
 *     both non-empty, which is the same MD5 under two names;
 *   - DIFFER_IN_SEQUENCE IS COUNTED PER NAME, when a name is found in exactly two rows each of
 *     which has it on one side only;
 *   - SUPERSET AND SUBSET REPLACE DIFFER_IN_SEQUENCES_PRESENT, but only when the missing entries
 *     all point the same way AND no naming discrepancy was found;
 *   - THE MD5 MODE DECIDES WHAT IS READ: USE_DICT refuses a dictionary with no M5, and
 *     ALWAYS_RECALCULATE ignores the M5 that is there, warning when the two disagree;
 *   - THE ANALYSIS IS PRINTED TO STDOUT while the table goes to the output file, so a run's two
 *     halves land in two places;
 *   - AND THE STATUS LINES COME OUT IN EnumSet ORDER, which is the enum's declaration order and
 *     not the order the flags were added.
 *
 * The FULL_ALIGNMENT mode runs mummer and is not measured here; FIND_SNPS_ONLY is.
 *
 * Output:
 *
 *     fasta\t<label>=<the whole fasta, escaped>
 *     dict\t<label>=<the whole dictionary, escaped>
 *     table\t<label>=<the output table, escaped>
 *     stdout\t<label>=<what the run printed, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CompareReferencesDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.reference.CompareReferences;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CompareReferencesDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("compare-references-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CompareReferencesDump: references compared by sequence MD5");

        // The base reference: two sequences.
        final Path base = fasta(dir, "base",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n");
        // The same bases under different names, which share the base's rows by MD5.
        final Path renamed = fasta(dir, "renamed",
                ">1\nACGTACGTAC\nGTACGTACGT\n>2\nTTTTTTTTTT\n");
        // The same names over different bases for chr2, which is two rows for one name.
        final Path altered = fasta(dir, "altered",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTA\n");
        // The base's sequences and one more, so it is a superset.
        final Path extra = fasta(dir, "extra",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n>chr3\nGGGGGGGGGG\n");
        // The base without chr2, which is a subset.
        final Path fewer = fasta(dir, "fewer", ">chr1\nACGTACGTAC\nGTACGTACGT\n");

        // A dictionary with no M5 at all, for the USE_DICT refusal.
        final Path noMd5 = fastaWithoutMd5(dir, "no-md5",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n");
        // And one whose M5 is a lie, which ALWAYS_RECALCULATE overrides.
        final Path wrongMd5 = fastaWithWrongMd5(dir, "wrong-md5",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n");

        run(dir, "identical", base, List.of(base));
        run(dir, "renamed", base, List.of(renamed));
        run(dir, "altered", base, List.of(altered));
        run(dir, "superset", extra, List.of(base));
        run(dir, "subset", fewer, List.of(base));
        // Three at once, so the pairs are every combination in order.
        run(dir, "three", base, List.of(renamed, extra));
        // The table by sequence name, over references whose names map to one row each.
        run(dir, "by-name", base, List.of(renamed), "--display-sequences-by-name");
        // The MD5 modes.
        run(dir, "use-dict-missing", base, List.of(noMd5), "--md5-calculation-mode", "USE_DICT");
        run(dir, "always-recalculate", base, List.of(wrongMd5),
                "--md5-calculation-mode", "ALWAYS_RECALCULATE");
        run(dir, "use-dict-wrong", base, List.of(wrongMd5), "--md5-calculation-mode", "USE_DICT");
    }

    /** A fasta with its index and a dictionary carrying M5, as a user would have. */
    static Path fasta(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".fasta");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(file, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + file, "O=" + dir.resolve(label + ".dict")});
        System.out.printf("fasta\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        System.out.printf("dict\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(dir.resolve(label + ".dict")), dir)));
        return file;
    }

    /** The same, with every M5 field stripped from the dictionary. */
    static Path fastaWithoutMd5(final Path dir, final String label, final String text)
            throws Exception {
        final Path file = fasta(dir, label, text);
        final Path dict = dir.resolve(label + ".dict");
        final List<String> lines = new ArrayList<>();
        for (final String line : Files.readAllLines(dict, StandardCharsets.UTF_8)) {
            lines.add(line.replaceAll("\tM5:[0-9a-f]+", ""));
        }
        Files.writeString(dict, String.join("\n", lines) + "\n", StandardCharsets.UTF_8);
        System.out.printf("dict\t%s-stripped=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(dict), dir)));
        return file;
    }

    /** And one whose M5 is replaced by a digest of nothing at all. */
    static Path fastaWithWrongMd5(final Path dir, final String label, final String text)
            throws Exception {
        final Path file = fasta(dir, label, text);
        final Path dict = dir.resolve(label + ".dict");
        final List<String> lines = new ArrayList<>();
        for (final String line : Files.readAllLines(dict, StandardCharsets.UTF_8)) {
            lines.add(line.replaceAll("\tM5:[0-9a-f]+", "\tM5:d41d8cd98f00b204e9800998ecf8427e"));
        }
        Files.writeString(dict, String.join("\n", lines) + "\n", StandardCharsets.UTF_8);
        System.out.printf("dict\t%s-wrong=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(dict), dir)));
        return file;
    }

    static void run(final Path dir, final String label, final Path reference,
                    final List<Path> others, final String... extra) throws Exception {
        final Path out = dir.resolve("table-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", reference.toString(), "-O", out.toString()));
        for (final Path other : others) {
            argv.addAll(Arrays.asList("-refcomp", other.toString()));
        }
        argv.addAll(Arrays.asList(extra));

        final PrintStream realOut = System.out;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, StandardCharsets.UTF_8));
            new CompareReferences().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.setOut(realOut);
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        } finally {
            System.setOut(realOut);
        }
        System.out.printf("stdout\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(captured.toString(StandardCharsets.UTF_8), dir)));
        if (Files.exists(out)) {
            System.out.printf("table\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    /** The dump's own directory, whose absolute path reaches the dictionaries' UR fields. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
