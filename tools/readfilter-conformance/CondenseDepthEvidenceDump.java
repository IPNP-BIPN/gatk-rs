/*
 * CondenseDepthEvidence's output, taken from the reference.
 *
 * Adjacent depth-evidence bins merged, which is nine lines of accumulator and three arguments that
 * do not mean what their names suggest.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE LENGTH TESTED IS THE ONE ALREADY ACCUMULATED, BEFORE THE NEXT BIN IS ADDED, so an
 *     interval can come out LONGER than --max-interval-size: with a max of 200 and bins of 100,
 *     the accumulator flushes at 200 and the output holds intervals of exactly 200, but with a max
 *     of 150 it flushes only once the accumulated length has reached 200;
 *   - A RECORD SHORTER THAN --min-interval-size IS DROPPED, not merged and not written, so a run
 *     can lose bins entirely and say nothing;
 *   - AND THE MINIMUM IS APPLIED AGAIN TO THE LAST ACCUMULATOR, so a trailing short interval is
 *     dropped too;
 *   - ADJACENCY IS `end + 1 == start` ON THE SAME CONTIG, so a one-base gap breaks the run and a
 *     contig change always does;
 *   - THE COUNTS ARE SUMMED ELEMENTWISE, one column per sample, and the merged record takes the
 *     FIRST bin's start and the LAST bin's end;
 *   - THE FILE IS ZERO-BASED HALF-OPEN ON DISK AND ONE-BASED INSIDE, the codec adding one to the
 *     start it reads and subtracting it again when it writes;
 *   - THE HEADER IS REWRITTEN FROM THE INPUT'S SAMPLE NAMES, so the output declares the same
 *     columns whatever the records held;
 *   - A MINIMUM ABOVE THE MAXIMUM IS A UserException raised before anything is read;
 *   - AND AN OUTPUT WHOSE EXTENSION IMPLIES ANOTHER FEATURE TYPE IS REFUSED BY NAME, the message
 *     naming the type the codec would have written.
 *
 * Output:
 *
 *     input\t<label>=<the whole .rd.txt, escaped>
 *     condensed\t<label>=<the whole output, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CondenseDepthEvidenceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.sv.CondenseDepthEvidence;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CondenseDepthEvidenceDump {

    static final String HEADER = "#Chr\tStart\tEnd\tsampleA\tsampleB\n";

    /** One bin, written zero-based and half-open as the file holds it. */
    static String bin(final String contig, final int start, final int end, final int a,
                      final int b) {
        return contig + "\t" + start + "\t" + end + "\t" + a + "\t" + b + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("condense-depth-evidence-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CondenseDepthEvidenceDump: adjacent depth bins merged");

        // Ten adjacent bins of a hundred bases, then a gap, then two more, then another contig.
        final StringBuilder records = new StringBuilder();
        for (int i = 0; i < 10; i++) {
            records.append(bin("chr1", i * 100, (i + 1) * 100, i + 1, 100 - i));
        }
        // A one-base gap, which breaks the run.
        records.append(bin("chr1", 1001, 1101, 11, 90));
        records.append(bin("chr1", 1101, 1201, 12, 89));
        // A contig change, which breaks it whatever the coordinates say.
        records.append(bin("chr2", 1201, 1301, 13, 88));
        records.append(bin("chr2", 1301, 1401, 14, 87));
        final String input = HEADER + records;

        run(dir, "defaults", input);
        // A maximum that divides the bins evenly, and one that does not.
        run(dir, "max-200", input, "--max-interval-size", "200");
        run(dir, "max-150", input, "--max-interval-size", "150");
        // A minimum that drops the short trailing runs.
        run(dir, "min-300", input, "--min-interval-size", "300");
        // Both at once, which is where the two interact.
        run(dir, "min-200-max-300", input,
                "--min-interval-size", "200", "--max-interval-size", "300");
        // A file of one bin, whose only interval is the last accumulator.
        run(dir, "single", HEADER + bin("chr1", 0, 100, 1, 2));
        // The same, dropped by the minimum.
        run(dir, "single-dropped", HEADER + bin("chr1", 0, 100, 1, 2),
                "--min-interval-size", "200");
        // A minimum above the maximum.
        run(dir, "min-above-max", input,
                "--min-interval-size", "500", "--max-interval-size", "100");
        // And an output whose extension implies another feature type.
        runNamed(dir, "wrong-output-type", input, "condensed.baf.txt");
    }

    static void run(final Path dir, final String label, final String input, final String... extra)
            throws Exception {
        runNamed(dir, label, input, "condensed-" + label + ".rd.txt", extra);
    }

    static void runNamed(final Path dir, final String label, final String input,
                         final String outputName, final String... extra) throws Exception {
        final Path in = dir.resolve(label + ".rd.txt");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", in.toString()});
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));

        final Path out = dir.resolve(outputName);
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--depth-evidence", in.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new CondenseDepthEvidence().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("condensed\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
