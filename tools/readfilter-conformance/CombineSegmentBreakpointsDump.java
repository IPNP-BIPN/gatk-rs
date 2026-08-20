/*
 * CombineSegmentBreakpoints, taken from the reference.
 *
 * Two segment files in, one out, cut at every breakpoint either of them carries and annotated from
 * both. The fourth and last tool on the annotated-interval collection.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE BREAKPOINTS OF BOTH FILES ARE POOLED PER CONTIG AND SORTED WITH STARTS BEFORE ENDS at
 *     equal positions, which is what keeps a single-base segment from vanishing;
 *   - A PIECE BETWEEN TWO SEGMENTS IS SHRUNK AT BOTH ENDS, and is dropped entirely when the two
 *     were adjacent, so the output covers only what one of the inputs covered;
 *   - A GAP IN ONE FILE IS STILL AN INTERVAL if the other file covers it, and its columns from the
 *     empty side come out as EMPTY STRINGS rather than being absent;
 *   - COLUMNS PRESENT IN BOTH FILES ARE SUFFIXED WITH THE LABELS, and the labels default to `1`
 *     and `2`; a column present in one file only keeps its name;
 *   - THE OUTPUT COLUMNS ARE SORTED ALPHABETICALLY, so a suffixed column can sort away from the
 *     column it came from;
 *   - `--columns-of-interest` FILTERS BOTH FILES, and naming a column neither file carries is a
 *     refusal that lists the missing ones;
 *   - THE FIRST RECORD OF EACH FILE DECIDES ITS COLUMNS, so a file of no rows takes the whole run
 *     down;
 *   - AND THE OUTPUT SAM HEADER IS THE MERGE OF THE TWO INPUTS' plus the reference's dictionary.
 *
 * Output:
 *
 *     combined\t<label>=<the whole output file, escaped>
 *
 * The M5 and UR fields of the sequence lines are masked, as in PreprocessIntervalsDump: both come
 * from the .dict the reference indexer wrote and neither is this tool's.
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CombineSegmentBreakpointsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.utils.CombineSegmentBreakpoints;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CombineSegmentBreakpointsDump {

    static final String FIRST =
            "CONTIG\tSTART\tEND\tCALL\tMEAN\n"
            + "chr1\t1\t100\t+\t0.5\n"
            + "chr1\t101\t200\t0\t0.0\n"
            + "chr2\t1\t100\t-\t-0.5\n";

    static final String SECOND =
            "CONTIG\tSTART\tEND\tCALL\tNAME\n"
            + "chr1\t50\t150\t+\tsecond-one\n"
            + "chr1\t180\t260\t-\tsecond-two\n"
            + "chr2\t200\t240\t+\tsecond-three\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("combine-segment-breakpoints-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# CombineSegmentBreakpointsDump: the pieces two segment files cut each other into");

        // `--columns-of-interest` is required, so every run names the columns it wants kept.
        run("default", dir, fasta, FIRST, SECOND, "--columns-of-interest", "CALL",
                "--columns-of-interest", "MEAN", "--columns-of-interest", "NAME");
        // The labels, which decide the suffix on every shared column.
        run("labels", dir, fasta, FIRST, SECOND, "--columns-of-interest", "CALL",
                "--columns-of-interest", "MEAN", "--columns-of-interest", "NAME",
                "--labels", "tumour", "--labels", "normal");
        // One column of interest, which is in both files.
        run("column-of-interest", dir, fasta, FIRST, SECOND, "--columns-of-interest", "CALL");
        // Two, one from each file.
        run("columns-from-both", dir, fasta, FIRST, SECOND, "--columns-of-interest", "MEAN",
                "--columns-of-interest", "NAME");
        // A column neither file carries.
        run("missing-column", dir, fasta, FIRST, SECOND, "--columns-of-interest", "ABSENT");
        // Identical files, whose breakpoints add nothing.
        run("identical", dir, fasta, FIRST, FIRST, "--columns-of-interest", "CALL",
                "--columns-of-interest", "MEAN");
        // Adjacent segments in one file against one segment spanning both, so the pieces line up
        // exactly and nothing is left between.
        run("adjacent", dir, fasta,
                "CONTIG\tSTART\tEND\tA\nchr1\t1\t100\ta1\nchr1\t101\t200\ta2\n",
                "CONTIG\tSTART\tEND\tB\nchr1\t1\t200\tb1\n",
                "--columns-of-interest", "A", "--columns-of-interest", "B");
        // Segments that do not touch at all, so the piece between them is dropped.
        run("disjoint", dir, fasta,
                "CONTIG\tSTART\tEND\tA\nchr1\t1\t50\ta1\n",
                "CONTIG\tSTART\tEND\tB\nchr1\t150\t200\tb1\n",
                "--columns-of-interest", "A", "--columns-of-interest", "B");
        // A single-base segment, which the start-before-end sort is there for.
        run("single-base", dir, fasta,
                "CONTIG\tSTART\tEND\tA\nchr1\t100\t100\ta1\n",
                "CONTIG\tSTART\tEND\tB\nchr1\t1\t200\tb1\n",
                "--columns-of-interest", "A", "--columns-of-interest", "B");
        // One file carrying a SAM header, the other not.
        run("one-sam-header", dir, fasta,
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:240\nCONTIG\tSTART\tEND\tA\nchr1\t1\t100\ta1\n",
                "CONTIG\tSTART\tEND\tB\nchr1\t1\t100\tb1\n",
                "--columns-of-interest", "A", "--columns-of-interest", "B");
        // A file of no rows, which takes the run down.
        run("empty-file", dir, fasta, FIRST, "CONTIG\tSTART\tEND\tB\n",
                "--columns-of-interest", "CALL", "--columns-of-interest", "MEAN");
        // Overlapping segments within one file, which the combining refuses.
        run("overlapping-input", dir, fasta,
                "CONTIG\tSTART\tEND\tA\nchr1\t1\t100\ta1\nchr1\t50\t150\ta2\n",
                "CONTIG\tSTART\tEND\tB\nchr1\t1\t100\tb1\n",
                "--columns-of-interest", "A", "--columns-of-interest", "B");
    }

    static void run(final String label, final Path dir, final Path fasta, final String first,
                    final String second, final String... extra) throws Exception {
        final Path one = dir.resolve(label + "-1.seg");
        final Path two = dir.resolve(label + "-2.seg");
        Files.write(one, first.getBytes());
        Files.write(two, second.getBytes());
        final Path out = dir.resolve("combined-" + label + ".seg");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "--segments", one.toString(),
                "--segments", two.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new CombineSegmentBreakpoints().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("combined\t%s=%s%n", label,
                ReferenceQueryDump.escape(PreprocessIntervalsDump.masked(
                        new String(Files.readAllBytes(out)))));
    }
}
