/*
 * CompareIntervalLists, taken from the reference.
 *
 * Two interval files, each sorted and ALL-merged on its own, then walked against one another until
 * one of them runs out. The tool says one of three things: that they are equal, that they are not
 * and why, or nothing at all because it threw.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE COMPARISON IS NOT SYMMETRIC. `equateIntervals(master, test)` consumes master and pops
 *     test in step; a master that outlasts test pops an empty list and dies with a
 *     `NoSuchElementException`, while a test that outlasts master is reported politely;
 *   - A TEST INTERVAL WIDER THAN THE MASTER'S IS EQUAL, because only the master's remainder is
 *     pushed back: the test's overhang is never examined, so `chr1:10-20` against `chr1:1-100` is
 *     reported equal in that order and throws in the other;
 *   - THE REMAINDER IS PUSHED IN REVERSE ORDER, so a subtraction that leaves two pieces leaves
 *     them in coordinate order on the stack;
 *   - EACH FILE IS MERGED WITH THE `ALL` RULE BEFORE THE WALK, so abutting intervals in one file
 *     and the single interval they add up to in the other are equal;
 *   - EQUALITY IS PRINTED, INEQUALITY IS THROWN. The equal case writes `Intervals are equal` to
 *     standard output and returns zero; the unequal case is a `UserException` whose message
 *     carries the difference;
 *   - AND AN INTERVAL OFF THE END OF A CONTIG IS REFUSED BY THE PARSER, before either list is
 *     compared with anything.
 *
 * Output:
 *
 *     compare\t<label>=<the tool's own line, or the exception class and message>
 *
 * A file URI in a message has its directory masked: that is the working directory of the run and
 * not an answer of this tool's.
 *
 * Usage: CompareIntervalListsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.CompareIntervalLists;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;

public class CompareIntervalListsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("compare-interval-lists-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# CompareIntervalListsDump: whether two interval files agree");

        // The plain cases.
        run(dir, fasta, "equal", "chr1:1-100", "chr1:1-100");
        run(dir, fasta, "equal-two-contigs", "chr1:1-100\nchr2:1-50", "chr1:1-100\nchr2:1-50");
        // Abutting intervals in one file, merged by the ALL rule into what the other holds.
        run(dir, fasta, "abutting-merged", "chr1:1-50\nchr1:51-100", "chr1:1-100");
        // Overlapping intervals in one file, likewise merged.
        run(dir, fasta, "overlapping-merged", "chr1:1-60\nchr1:40-100", "chr1:1-100");
        // A test interval wider than the master's, which is equal in this order.
        run(dir, fasta, "test-wider", "chr1:10-20", "chr1:1-100");
        // The same pair the other way, where master outlasts test.
        run(dir, fasta, "master-wider", "chr1:1-100", "chr1:10-20");
        // An extra interval on the test side, reported rather than thrown.
        run(dir, fasta, "test-longer", "chr1:1-100", "chr1:1-100\nchr2:1-50");
        // An extra interval on the master side, which pops an empty list.
        run(dir, fasta, "master-longer", "chr1:1-100\nchr2:1-50", "chr1:1-100");
        // Disjoint intervals, which is the incompatible branch.
        run(dir, fasta, "disjoint", "chr1:1-100", "chr1:201-240");
        // Disjoint on different contigs, the same branch by another route.
        run(dir, fasta, "different-contigs", "chr1:1-100", "chr2:1-100");
        // A master interval that the test cuts in two, leaving a piece on each side.
        run(dir, fasta, "test-inside-master", "chr1:1-100", "chr1:40-60");
        // Empty on one side, which cannot even be given: an interval file with no lines.
        run(dir, fasta, "empty-master", "", "chr1:1-100");
        run(dir, fasta, "empty-test", "chr1:1-100", "");
        // Off the end of the contig, which the parser refuses.
        run(dir, fasta, "off-the-end", "chr1:1-1000", "chr1:1-100");
    }

    static void run(final Path dir, final Path fasta, final String label, final String first,
                    final String second) throws Exception {
        final Path one = write(dir, label + "-1.list", first);
        final Path two = write(dir, label + "-2.list", second);
        final PrintStream out = System.out;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        String answer;
        try {
            System.setOut(new PrintStream(captured));
            new CompareIntervalLists().instanceMain(new String[] {
                    "-R", fasta.toString(), "-L", one.toString(), "-L2", two.toString()});
            answer = captured.toString().trim();
        } catch (final Exception | AssertionError e) {
            answer = e.getClass().getName() + ":" + String.valueOf(e.getMessage());
        } finally {
            System.setOut(out);
        }
        // The malformed-file message names the file by URI, so the directory the run happened in
        // is masked out: it is the container's working directory and not this tool's answer.
        answer = answer.replaceAll("file:[^ ]*/([^/ ]+\\.list)", "file:<masked>/$1");
        System.out.printf("compare\t%s=%s%n", label, ReferenceQueryDump.escape(answer));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, (text.isEmpty() ? "" : text + "\n").getBytes());
        return path;
    }
}
