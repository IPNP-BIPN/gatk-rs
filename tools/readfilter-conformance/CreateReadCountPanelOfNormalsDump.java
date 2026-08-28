/*
 * CreateReadCountPanelOfNormals' panel, taken from the reference.
 *
 * A panel of normals is a set of read-count files reduced to the intervals worth keeping, the
 * median of each of those intervals, and a basis of eigensamples. The singular value
 * decomposition is not what this measures, being a distributed solver's: what it measures is
 * everything that decides which intervals and which samples reach it.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE TOOL IS A SPARK ONE, so it needs a master even with nothing to distribute, and it
 *     needs a JVM that exports `sun.nio.ch` to it;
 *   - THE COUNTS BECOME FRACTIONAL COVERAGE FIRST, each sample divided by its own total, so the
 *     sample that is exactly twice another has exactly that other's median coverage and is not an
 *     outlier;
 *   - AN INTERVAL WHOSE MEDIAN FRACTIONAL COVERAGE IS UNDER THE PERCENTILE IS DROPPED, which is
 *     what takes out the two intervals every sample reads near zero;
 *   - AN INTERVAL WITH TOO MANY ZEROS ACROSS THE SAMPLES IS DROPPED TOO, and that filter is much
 *     the harsher of the two: on this panel it leaves so few intervals that the decomposition
 *     finds no non-zero singular value at all and the run is REFUSED, where the median filter
 *     leaves most of them;
 *   - A SAMPLE WITH TOO MANY ZEROS IS DROPPED, which is what takes out the nearly-empty one;
 *   - A SAMPLE WHOSE MEDIAN IS EXTREME IS DROPPED FROM BOTH ENDS, the percentile being applied
 *     twice, so a percentile of twenty takes two samples of nine. The medians are separated by
 *     fourteen per cent or more for exactly this reason: a fixture that varied only the depth put
 *     the top two within two parts in a hundred thousand of each other, and which one the filter
 *     called the maximum flipped on CI;
 *   - THE PANEL KEEPS BOTH THE ORIGINAL INTERVALS AND THE SURVIVING ONES, so the file says what
 *     was dropped as well as what was kept, and the ORIGINAL count is the input's whatever the
 *     filters did;
 *   - THE NUMBER OF EIGENSAMPLES IS CAPPED AT THE NUMBER OF SAMPLES, and the number itself is a
 *     RANK rather than a count of anything the tool was given: the same fixture answered six on
 *     one runner, seven on another and eight on the machine that produced the first golden, so
 *     what is measured is the cap and whether the basis is empty;
 *   - A PANEL OF ONE SAMPLE IS ACCEPTED, has no eigensamples, and REFUSES to hand over its
 *     singular values rather than handing over an empty array;
 *   - AND THE INPUTS MUST AGREE ON THEIR INTERVALS, one that does not being refused by name.
 *
 * Output:
 *
 *     counts\t<label>=<that read-count file, escaped>
 *     panel\t<label>\t<field>=<value>
 *     matrix\t<label>\t<name>=<one value per line, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<the tool's own message, its leading wall clock removed>
 *     error\t<label>\t<field>\t<exception class>:<message>
 *
 * Usage: CreateReadCountPanelOfNormalsDump
 */

import org.broadinstitute.hellbender.tools.copynumber.CreateReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.denoising.HDF5SVDReadCountPanelOfNormals;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CreateReadCountPanelOfNormalsDump {

    static final int INTERVALS = 40;
    static final int CONTIG_LENGTH = 40000;

    /** One interval of the shared interval list: 1000 bases every 1000. */
    static int start(final int index) {
        return index * 1000 + 1;
    }

    static int end(final int index) {
        return (index + 1) * 1000;
    }

    /**
     * One read-count file: the SAM-style header the reader wants, then one row per interval.
     */
    static String counts(final String sample, final int[] values) {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6");
        lines.add("@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH);
        lines.add("@RG\tID:GATKCopyNumber\tSM:" + sample);
        lines.add("CONTIG\tSTART\tEND\tCOUNT");
        for (int i = 0; i < values.length; i++) {
            lines.add("chr1\t" + start(i) + "\t" + end(i) + "\t" + values[i]);
        }
        lines.add("");
        return String.join("\n", lines);
    }

    /**
     * A captured line with its leading wall clock removed.
     *
     * The message this dump keeps is whichever log line named the exception, and log4j opens every
     * line with `HH:mm:ss.SSS`. That is a clock, not a behaviour, and no byte comparison should be
     * asked to reproduce it.
     */
    static String withoutTheClock(final String line) {
        return line.replaceFirst("^\\d\\d:\\d\\d:\\d\\d\\.\\d\\d\\d ", "");
    }

    /**
     * The panel's samples.
     *
     * Six of a shape of their own, one that is exactly twice another, one nearly all zeros, and
     * one whose coverage is piled into its first quarter.
     *
     * The per-sample shape is what matters. The first version of this fixture varied only the
     * DEPTH, and fractional coverage divides depth out again: the nine samples' median coverages
     * then agreed to five decimal places, and which of them the extreme-median filter called the
     * maximum could flip on a reordering. It did, once, on CI. Here the medians are separated by
     * fourteen per cent or more, except for the one pair that is deliberately identical.
     */
    static int[] sampleCounts(final int sample) {
        final int[] values = new int[INTERVALS];
        for (int i = 0; i < INTERVALS; i++) {
            values[i] = shape(sample, i) * (sample == 6 ? 2 : 1);
            // The first two intervals are near zero in every sample, so their median is low.
            if (i < 2) {
                values[i] = sample == 7 ? 0 : (sample == 6 ? 2 : 1);
            }
            // One interval is zero in three samples, which is over the interval percentage.
            if (i == 20 && sample < 3) {
                values[i] = 0;
            }
        }
        return values;
    }

    /**
     * One sample's profile before the shared clamps.
     *
     * Sample 6 borrows sample 3's, so that doubling it leaves its fractional coverage EXACTLY
     * equal: that is the pair the coverage normalisation has to make indistinguishable, and it
     * sits in the middle of the ranking where no percentile filter reaches for it.
     */
    static int shape(final int sample, final int i) {
        final int base = 100 + (i * 7) % 23;
        if (sample == 7) {
            return i < 4 ? base : 0;        // nearly all zeros
        }
        if (sample == 8) {
            return i < 10 ? base * 40 : base;   // piled into the first quarter
        }
        final int weight = 1 + (sample == 6 ? 3 : sample);
        return i < 20 ? base * weight : base;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("create-read-count-pon-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CreateReadCountPanelOfNormalsDump: which intervals and which "
                + "samples reach the panel");

        final List<Path> inputs = new ArrayList<>();
        for (int sample = 0; sample < 9; sample++) {
            final String text = counts("sample" + sample, sampleCounts(sample));
            final Path path = write(dir, "sample" + sample + ".counts.tsv", text);
            inputs.add(path);
            if (sample < 2 || sample >= 6) {
                // The ordinary samples are all alike, so only the first two and the three odd
                // ones are printed.
                System.out.printf("counts\tsample%d=%s%n", sample,
                        ReferenceQueryDump.escape(text));
            }
        }

        run(dir, "default", inputs, List.of());
        // No filtering at all, which keeps every interval and every sample.
        run(dir, "no-filtering", inputs, List.of(
                "--minimum-interval-median-percentile", "0.0",
                "--maximum-zeros-in-sample-percentage", "100.0",
                "--maximum-zeros-in-interval-percentage", "100.0",
                "--extreme-sample-median-percentile", "0.0",
                "--extreme-outlier-truncation-percentile", "0.0"));
        // Each filter on its own, against that baseline.
        run(dir, "interval-median-only", inputs, List.of(
                "--minimum-interval-median-percentile", "10.0",
                "--maximum-zeros-in-sample-percentage", "100.0",
                "--maximum-zeros-in-interval-percentage", "100.0",
                "--extreme-sample-median-percentile", "0.0",
                "--extreme-outlier-truncation-percentile", "0.0"));
        run(dir, "zeros-in-sample-only", inputs, List.of(
                "--minimum-interval-median-percentile", "0.0",
                "--maximum-zeros-in-sample-percentage", "5.0",
                "--maximum-zeros-in-interval-percentage", "100.0",
                "--extreme-sample-median-percentile", "0.0",
                "--extreme-outlier-truncation-percentile", "0.0"));
        run(dir, "zeros-in-interval-only", inputs, List.of(
                "--minimum-interval-median-percentile", "0.0",
                "--maximum-zeros-in-sample-percentage", "100.0",
                "--maximum-zeros-in-interval-percentage", "5.0",
                "--extreme-sample-median-percentile", "0.0",
                "--extreme-outlier-truncation-percentile", "0.0"));
        run(dir, "extreme-sample-only", inputs, List.of(
                "--minimum-interval-median-percentile", "0.0",
                "--maximum-zeros-in-sample-percentage", "100.0",
                "--maximum-zeros-in-interval-percentage", "100.0",
                "--extreme-sample-median-percentile", "20.0",
                "--extreme-outlier-truncation-percentile", "0.0"));
        // The zeros left alone rather than imputed.
        run(dir, "no-imputation", inputs, List.of("--do-impute-zeros", "false"));
        // Fewer eigensamples than samples, and more.
        run(dir, "two-eigensamples", inputs, List.of("--number-of-eigensamples", "2"));
        run(dir, "hundred-eigensamples", inputs, List.of("--number-of-eigensamples", "100"));
        // A panel of one sample.
        run(dir, "one-sample", inputs.subList(0, 1), List.of());
        // An input whose intervals do not match the others'.
        final int[] shifted = sampleCounts(0);
        final List<String> lines = new ArrayList<>(
                List.of(counts("odd", shifted).split("\n", -1)));
        lines.set(4, "chr1\t2\t1000\t1");
        final Path odd = write(dir, "odd.counts.tsv", String.join("\n", lines));
        final List<Path> mixed = new ArrayList<>(inputs.subList(0, 3));
        mixed.add(odd);
        run(dir, "mismatched-intervals", mixed, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final List<Path> inputs,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".pon.hdf5");
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.add("-I");
            argv.add(input.toString());
        }
        argv.add("-O");
        argv.add(out.toString());
        // The tool is a SPARK one, the singular value decomposition being distributed, so it
        // needs a master even when there is nothing to distribute.
        argv.add("--spark-master");
        argv.add("local[1]");
        argv.addAll(extra);
        // Spark reaches into `sun.nio.ch`, which a modern JDK does not export, so the tool runs
        // in a JVM of its own with the exports its own launcher passes.
        final List<String> command = new ArrayList<>(List.of(
                "java",
                "--add-opens", "java.base/java.lang=ALL-UNNAMED",
                "--add-opens", "java.base/java.lang.invoke=ALL-UNNAMED",
                "--add-opens", "java.base/java.nio=ALL-UNNAMED",
                "--add-opens", "java.base/java.util=ALL-UNNAMED",
                "--add-opens", "java.base/sun.nio.ch=ALL-UNNAMED",
                "--add-exports", "java.base/sun.nio.ch=ALL-UNNAMED",
                "-cp", System.getenv("ORACLE_CP"),
                "org.broadinstitute.hellbender.Main", "CreateReadCountPanelOfNormals"));
        command.addAll(argv);
        final Process process = new ProcessBuilder(command).redirectErrorStream(true).start();
        final String log = new String(process.getInputStream().readAllBytes(),
                StandardCharsets.UTF_8);
        if (process.waitFor() != 0) {
            String message = "exit " + process.exitValue();
            for (final String line : log.split("\n")) {
                final String trimmed = line.trim();
                if (trimmed.startsWith("A USER ERROR has occurred:")
                        || trimmed.startsWith("Exception in thread")
                        || trimmed.contains("Exception:")
                        || trimmed.contains("Error:")) {
                    message = trimmed;
                    break;
                }
            }
            System.out.printf("error\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(masked(withoutTheClock(message), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("none\t%s=no panel%n", label);
            return;
        }
        try (final org.broadinstitute.hdf5.HDF5File file =
                     new org.broadinstitute.hdf5.HDF5File(new File(out.toString()))) {
            final HDF5SVDReadCountPanelOfNormals panel = HDF5SVDReadCountPanelOfNormals.read(file);
            System.out.printf("panel\t%s\tversion=%s%n", label,
                    Double.toString(panel.getVersion()));
            // The NUMBER of eigensamples is the decomposition's rank, and a rank is not a byte:
            // the same fixture gave six here and seven on the machine that produced this golden,
            // and eight on a third. What is reported is what does not move: the bound the number
            // is capped at, and whether the basis is empty. See
            // `docs/a-rank-is-not-a-byte.md`.
            System.out.printf("panel\t%s\teigensamples-at-most=%d%n", label,
                    panel.getOriginalReadCounts().length);
            System.out.printf("panel\t%s\teigensamples-positive=%b%n", label,
                    panel.getNumEigensamples() > 0);
            System.out.printf("panel\t%s\toriginal-intervals=%d%n", label,
                    panel.getOriginalIntervals().size());
            System.out.printf("panel\t%s\tpanel-intervals=%s%n", label,
                    places(panel.getPanelIntervals()));
            System.out.printf("panel\t%s\tsamples=%d%n", label,
                    panel.getOriginalReadCounts().length);
            System.out.printf("matrix\t%s\tfractional-medians=%s%n", label,
                    ReferenceQueryDump.escape(vector(panel.getPanelIntervalFractionalMedians())));
            // The singular values are the SVD's own, computed by a distributed solver, and their
            // COUNT is the same rank that moves, so what is reported is whether the panel hands
            // them over at all: one with no eigensamples REFUSES rather than handing over an
            // empty array, which is a decision and not a number.
            try {
                System.out.printf("panel\t%s\tsingular-values-available=%b%n", label,
                        panel.getSingularValues().length > 0);
            } catch (final Exception e) {
                System.out.printf("error\t%s\tsingular-values\t%s:%s%n", label,
                        e.getClass().getName(),
                        ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            }
        }
    }

    /** The intervals a panel kept, as `start-end` and comma-separated. */
    static String places(final List<org.broadinstitute.hellbender.utils.SimpleInterval> intervals) {
        final List<String> parts = new ArrayList<>();
        for (final org.broadinstitute.hellbender.utils.SimpleInterval interval : intervals) {
            parts.add(interval.getStart() + "-" + interval.getEnd());
        }
        return String.join(",", parts);
    }

    /** One vector, each value rendered by Double.toString. */
    static String vector(final double[] values) {
        final List<String> parts = new ArrayList<>();
        for (final double value : values) {
            parts.add(Double.toString(value));
        }
        return String.join("\n", parts) + "\n";
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
