/*
 * CallCopyRatioSegments, taken from the reference.
 *
 * Pure arithmetic over a table: no reads, no reference, no walker. Copy-ratio segments in, called
 * segments out, with a length-weighted mean and standard deviation computed twice -- once to find
 * the outliers and once without them.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE NEUTRAL BOUNDS ARE INCLUSIVE AT BOTH ENDS, and they are compared against 2^log2 rather
 *     than against the log2 value itself, so a segment at exactly 0.9 or exactly 1.1 is neutral;
 *   - THE STATISTICS ARE COMPUTED TWICE. The first pass is over every copy-neutral segment; the
 *     outlier filter uses THAT mean and standard deviation, and the second pass -- over what
 *     survives -- is what the calling threshold is compared against;
 *   - THE OUTLIER FILTER IS INCLUSIVE TOO: a segment exactly two standard deviations out is kept;
 *   - THE STANDARD DEVIATION'S DENOMINATOR IS `((n - 1) / n) * totalLength`, which is neither the
 *     population nor the sample form, and with ONE copy-neutral segment it is zero -- the answer is
 *     an infinity or a NaN rather than an error;
 *   - A SEGMENT OUTSIDE THE BOUNDS IS CALLED BY DEVIATION FROM THE FILTERED MEAN, so a run whose
 *     copy-neutral set is empty calls everything NEUTRAL: the mean is NaN and every comparison
 *     against it is false;
 *   - THE LENGTH IS THE INTERVAL'S, `end - start + 1`, so a one-base segment still weighs one;
 *   - AND THE TOOL WRITES TWO FILES: the called segments and a legacy `.igv.seg` beside them, whose
 *     name is the output's plus a suffix and whose columns are IGV's rather than GATK's.
 *
 * Output:
 *
 *     called\t<label>\t<the called-segments file, escaped>
 *     legacy\t<label>\t<the .igv.seg file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CallCopyRatioSegmentsDump
 */

import org.broadinstitute.hellbender.tools.copynumber.CallCopyRatioSegments;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CallCopyRatioSegmentsDump {

    /** The table header every input shares: a SAM header, then the five columns. */
    static final String HEADER =
            "@HD\tVN:1.6\n"
            + "@SQ\tSN:chr1\tLN:1000\n"
            + "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
            + "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("callcopyratio-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CallCopyRatioSegmentsDump: copy-ratio segments called");

        // A spread of segments: three inside the neutral band, one just above it, one far above,
        // one far below. log2 of 1 is 0, of 0.9 is -0.152, of 1.1 is 0.1375.
        final String spread = HEADER
                + "chr1\t1\t100\t10\t0.0\n"
                + "chr1\t101\t200\t10\t-0.15200309344504997\n"   // exactly 0.9
                + "chr1\t201\t300\t10\t0.13750352374993502\n"    // exactly 1.1
                + "chr1\t301\t400\t10\t0.5\n"
                + "chr1\t401\t500\t10\t1.5\n"
                + "chr1\t501\t600\t10\t-1.5\n";
        run("spread", dir, write(dir, "spread.cr.seg", spread));

        // Every segment inside the band, so nothing is called anything but neutral.
        final String allNeutral = HEADER
                + "chr1\t1\t100\t10\t0.0\n"
                + "chr1\t101\t200\t10\t0.05\n"
                + "chr1\t201\t300\t10\t-0.05\n";
        run("all-neutral", dir, write(dir, "all-neutral.cr.seg", allNeutral));

        // Exactly one copy-neutral segment, so the standard deviation divides by zero.
        final String oneNeutral = HEADER
                + "chr1\t1\t100\t10\t0.0\n"
                + "chr1\t101\t200\t10\t1.5\n"
                + "chr1\t201\t300\t10\t-1.5\n";
        run("one-neutral", dir, write(dir, "one-neutral.cr.seg", oneNeutral));

        // No copy-neutral segment at all: the mean is a NaN and every comparison against it fails.
        final String noneNeutral = HEADER
                + "chr1\t1\t100\t10\t1.5\n"
                + "chr1\t101\t200\t10\t-1.5\n";
        run("none-neutral", dir, write(dir, "none-neutral.cr.seg", noneNeutral));

        // Segments of very different lengths, so the weighting is visible.
        final String weighted = HEADER
                + "chr1\t1\t10\t10\t0.0\n"
                + "chr1\t11\t910\t10\t0.09\n"
                + "chr1\t911\t920\t10\t-0.09\n"
                + "chr1\t921\t1000\t10\t0.6\n";
        run("weighted", dir, write(dir, "weighted.cr.seg", weighted));

        // The bounds moved, which changes which segments are neutral and therefore the statistics.
        run("wide-bounds", dir, write(dir, "wide.cr.seg", spread),
                "--neutral-segment-copy-ratio-lower-bound", "0.5",
                "--neutral-segment-copy-ratio-upper-bound", "2.0");
        // A calling threshold of zero, so every non-neutral segment is called one way or the other.
        run("zero-threshold", dir, write(dir, "zero.cr.seg", spread),
                "--calling-copy-ratio-z-score-threshold", "0.0");
        // A lower bound above the upper one, which is an argument error.
        run("inverted-bounds", dir, write(dir, "inverted.cr.seg", spread),
                "--neutral-segment-copy-ratio-lower-bound", "2.0",
                "--neutral-segment-copy-ratio-upper-bound", "1.0");
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        return path;
    }

    static void run(final String label, final Path dir, final Path input, final String... extra)
            throws Exception {
        final Path out = dir.resolve("called-" + label + ".called.seg");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new CallCopyRatioSegments().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("called\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
        final Path legacy = dir.resolve("called-" + label + ".called.igv.seg");
        System.out.printf("legacy\t%s\t%s%n", label, Files.exists(legacy)
                ? ReferenceQueryDump.escape(new String(Files.readAllBytes(legacy))) : "(none)");
    }
}
