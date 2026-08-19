/*
 * DenoiseReadCounts without a panel of normals, taken from the reference.
 *
 * With no panel and no GC annotations the tool does the one thing its documentation calls
 * "standardization": fractional coverage, divide by the sample median, log2, subtract the median
 * again. Four steps over one row of numbers, and every one of them has a rule.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE FIRST STEP DIVIDES BY THE SUM, so the row becomes fractional coverage and every later
 *     step is scale-free: doubling every count changes nothing in the output;
 *   - THE MEDIAN IS COMMONS-MATH'S, which for an even count interpolates between the two middle
 *     values rather than taking either;
 *   - THE LOG IS `Math.log(x) * INV_LOG_2` rather than a base-two logarithm, taken AFTER the
 *     division, so a value equal to the median is exactly zero and the output is centred on it;
 *   - THE MEDIAN IS SUBTRACTED TWICE IN EFFECT: once as a division before the log and once as a
 *     subtraction after it, and the second median is taken over the LOG values rather than over the
 *     original ones, so it is not the logarithm of the first;
 *   - A ZERO COUNT IS FLOORED RATHER THAN INFINITE. `safeLog2` answers `log2(1e-9)` for anything
 *     below `1e-9`, so an interval no read started in reads -29.897353 and not minus infinity, and
 *     the row's median is computed over that floor like any other number;
 *   - AND THE TWO OUTPUT FILES ARE IDENTICAL when there is no panel: the tool writes the
 *     standardized result as the denoised one as well, which is a fact about the files rather than
 *     about the arithmetic.
 *
 * Output:
 *
 *     standardized\t<label>\t<the standardized file, escaped>
 *     denoised\t<label>\t<the denoised file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: DenoiseReadCountsDump
 */

import org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DenoiseReadCountsDump {

    /** The read-counts table's header: a SAM header, then the three columns. */
    static final String HEADER =
            "@HD\tVN:1.6\n"
            + "@SQ\tSN:chr1\tLN:1000\n"
            + "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
            + "CONTIG\tSTART\tEND\tCOUNT\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("denoisereadcounts-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# DenoiseReadCountsDump: standardization without a panel of normals");

        // Five intervals with counts that are not all equal.
        final String plain = HEADER
                + "chr1\t1\t100\t10\n"
                + "chr1\t101\t200\t20\n"
                + "chr1\t201\t300\t30\n"
                + "chr1\t301\t400\t40\n"
                + "chr1\t401\t500\t50\n";
        run("plain", dir, write(dir, "plain.counts.tsv", plain));

        // The same counts doubled, which the fractional-coverage step makes identical.
        final String doubled = HEADER
                + "chr1\t1\t100\t20\n"
                + "chr1\t101\t200\t40\n"
                + "chr1\t201\t300\t60\n"
                + "chr1\t301\t400\t80\n"
                + "chr1\t401\t500\t100\n";
        run("doubled", dir, write(dir, "doubled.counts.tsv", doubled));

        // An even number of intervals, so the median interpolates.
        final String even = HEADER
                + "chr1\t1\t100\t10\n"
                + "chr1\t101\t200\t20\n"
                + "chr1\t201\t300\t30\n"
                + "chr1\t301\t400\t40\n";
        run("even", dir, write(dir, "even.counts.tsv", even));

        // A zero count, which becomes a negative infinity.
        final String withZero = HEADER
                + "chr1\t1\t100\t0\n"
                + "chr1\t101\t200\t20\n"
                + "chr1\t201\t300\t30\n";
        run("with-zero", dir, write(dir, "zero.counts.tsv", withZero));

        // Every count equal, so every standardized value is zero.
        final String flat = HEADER
                + "chr1\t1\t100\t7\n"
                + "chr1\t101\t200\t7\n"
                + "chr1\t201\t300\t7\n";
        run("flat", dir, write(dir, "flat.counts.tsv", flat));

        // One interval, where the median is the value itself.
        final String single = HEADER + "chr1\t1\t100\t42\n";
        run("single", dir, write(dir, "single.counts.tsv", single));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        return path;
    }

    static void run(final String label, final Path dir, final Path input, final String... extra)
            throws Exception {
        final Path standardized = dir.resolve("standardized-" + label + ".tsv");
        final Path denoised = dir.resolve("denoised-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(),
                "--standardized-copy-ratios", standardized.toString(),
                "--denoised-copy-ratios", denoised.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new DenoiseReadCounts().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("standardized\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(standardized))));
        System.out.printf("denoised\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(denoised))));
    }
}
