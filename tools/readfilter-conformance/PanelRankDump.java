/*
 * The number of eigensamples CreateReadCountPanelOfNormals writes, and the singular values under it.
 *
 * `CreateReadCountPanelOfNormalsDump` carries the cap and not the count, because the count is the
 * rank Spark's solver decides, and that rank answered 8, 6 and 7 on one fixture
 * (docs/a-rank-is-not-a-byte.md). Measured since (docs/spark-parallel.md): what moved it is the ORDER
 * the driver sums the Gramian in. At `local[1]` the partitions arrive in their order, the LAPACK is
 * F2J's in a strict JVM, and twelve runs of one fixture gave one set of bits; at `local[4]` twelve
 * runs gave eleven sets of bits and two ranks. So this suite runs at `local[1]` and nothing else.
 *
 * Output, one row per run:
 *
 *     rank\t<label>\trequested=<r>\tsurviving=<n>\teigensamples=<e>
 *     values\t<label>\t<each singular value's bits, as sixteen hex digits, comma separated>
 *     error\t<label>\t<the first line that names what went wrong>
 *
 * Usage: PanelRankDump
 */

import org.broadinstitute.hellbender.tools.copynumber.denoising.HDF5SVDReadCountPanelOfNormals;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class PanelRankDump {

    /** The arguments that skip every filter, so a count matrix reaches the solver whole. */
    static final List<String> NO_FILTERS = List.of(
            "--minimum-interval-median-percentile", "0.0",
            "--maximum-zeros-in-sample-percentage", "100.0",
            "--maximum-zeros-in-interval-percentage", "100.0",
            "--extreme-sample-median-percentile", "0.0",
            "--extreme-outlier-truncation-percentile", "0.0");

    /** A count that looks random and needs no generator to rebuild: the port's test writes the same. */
    static long noise(final int sample, final int interval) {
        final long i = interval + 7L;
        return 200 + ((sample + 1L) * i * i * 2654435761L + sample * 97L) % 1009;
    }

    /** `samples` independent profiles. Their rank is full, so no value is near the threshold. */
    static int[][] independent(final int samples, final int intervals) {
        final int[][] counts = new int[samples][intervals];
        for (int s = 0; s < samples; s++) {
            for (int i = 0; i < intervals; i++) {
                counts[s][i] = (int) noise(s, i);
            }
        }
        return counts;
    }

    /**
     * `samples` copies of `families` profiles, each a whole multiple of its family's first: the
     * fractional coverage divides the multiple out, so the matrix has exactly that many distinct
     * rows and every value past them is rounding noise, which is what the threshold compares.
     */
    static int[][] families(final int samples, final int families, final int intervals) {
        final int[][] counts = new int[samples][intervals];
        for (int s = 0; s < samples; s++) {
            for (int i = 0; i < intervals; i++) {
                counts[s][i] = (int) noise(s % families, i) * (1 + s / families);
            }
        }
        return counts;
    }

    /** The other suite's nine samples: the fixture the rank moved on. */
    static int[][] nine() {
        final int[][] counts = new int[9][];
        for (int s = 0; s < 9; s++) {
            counts[s] = CreateReadCountPanelOfNormalsDump.sampleCounts(s);
        }
        return counts;
    }

    static String text(final String sample, final int[] values) {
        final StringBuilder out = new StringBuilder("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n");
        out.append("@RG\tID:GATKCopyNumber\tSM:").append(sample).append("\nCONTIG\tSTART\tEND\tCOUNT\n");
        for (int i = 0; i < values.length; i++) {
            out.append("chr1\t").append(i * 1000 + 1).append('\t').append((i + 1) * 1000).append('\t')
                    .append(values[i]).append('\n');
        }
        return out.toString();
    }

    static void run(final Path dir, final String label, final int[][] counts, final List<String> extra)
            throws Exception {
        // A hung Spark context must not hold the container: nothing here takes five minutes.
        final List<String> command = new ArrayList<>(List.of(
                "timeout", "300", "java",
                "--add-opens", "java.base/java.lang=ALL-UNNAMED",
                "--add-opens", "java.base/java.lang.invoke=ALL-UNNAMED",
                "--add-opens", "java.base/java.nio=ALL-UNNAMED",
                "--add-opens", "java.base/java.util=ALL-UNNAMED",
                "--add-opens", "java.base/sun.nio.ch=ALL-UNNAMED",
                "--add-exports", "java.base/sun.nio.ch=ALL-UNNAMED",
                "-cp", System.getenv("ORACLE_CP"),
                "org.broadinstitute.hellbender.Main", "CreateReadCountPanelOfNormals"));
        for (int s = 0; s < counts.length; s++) {
            final Path input = dir.resolve(label + "-" + s + ".counts.tsv");
            Files.writeString(input, text("sample" + s, counts[s]), StandardCharsets.UTF_8);
            command.add("-I");
            command.add(input.toString());
        }
        final Path out = dir.resolve(label + ".pon.hdf5");
        command.addAll(List.of("-O", out.toString(), "--spark-master", "local[1]"));
        command.addAll(extra);
        final Process process = new ProcessBuilder(command).redirectErrorStream(true).start();
        final String log = new String(process.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
        if (process.waitFor() != 0) {
            String message = "exit " + process.exitValue();
            for (final String line : log.split("\n")) {
                final String trimmed = line.trim();
                if (trimmed.startsWith("A USER ERROR has occurred:")
                        || trimmed.contains("No non-zero singular values")) {
                    message = "No non-zero singular values";
                    break;
                }
            }
            System.out.printf("error\t%s\t%s%n", label, message);
            return;
        }
        try (final org.broadinstitute.hdf5.HDF5File file =
                     new org.broadinstitute.hdf5.HDF5File(new File(out.toString()))) {
            final HDF5SVDReadCountPanelOfNormals panel = HDF5SVDReadCountPanelOfNormals.read(file);
            final int requested = extra.contains("--number-of-eigensamples")
                    ? Integer.parseInt(extra.get(extra.indexOf("--number-of-eigensamples") + 1)) : 20;
            final double[] values = panel.getSingularValues();
            System.out.printf("rank\t%s\trequested=%d\tsurviving=%d\teigensamples=%d%n", label, requested,
                    file.readStringArray("/panel/sample_filenames").length, panel.getNumEigensamples());
            final StringBuilder bits = new StringBuilder();
            for (final double value : values) {
                bits.append(bits.length() == 0 ? "" : ",")
                        .append(String.format("%016x", Double.doubleToRawLongBits(value)));
            }
            System.out.printf("values\t%s\t%s%n", label, bits);
        }
        Files.deleteIfExists(out);
    }

    static List<String> and(final List<String> first, final String... more) {
        final List<String> all = new ArrayList<>(first);
        all.addAll(List.of(more));
        return all;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("panel-rank-dump").toAbsolutePath();
        Files.createDirectories(dir);
        System.out.println("# PanelRankDump: the eigensamples a panel writes, at local[1]");

        // The fixture the rank moved on. The surviving samples are eight, and two of them equal.
        run(dir, "nine-defaults", nine(), List.of());
        run(dir, "nine-no-filters", nine(), NO_FILTERS);
        for (final int k : new int[]{4, 5, 6, 7, 8, 100}) {
            run(dir, "nine-k" + k, nine(), List.of("--number-of-eigensamples", Integer.toString(k)));
        }
        // Local ARPACK when k < n / 3, whose rank is the request while the values clear the threshold.
        run(dir, "nine-k2-lapack", nine(), List.of("--number-of-eigensamples", "2"));
        run(dir, "nine-k1-arpack", nine(), List.of("--number-of-eigensamples", "1"));

        // Rank-deficient panels: every value past the families' is noise.
        run(dir, "families-6-of-3", families(6, 3, 30), NO_FILTERS);
        run(dir, "families-12-of-3", families(12, 3, 30), NO_FILTERS);
        run(dir, "families-12-of-5", families(12, 5, 60), NO_FILTERS);
        run(dir, "families-25-of-4", families(25, 4, 40), NO_FILTERS);
        // More intervals than the hundred slices, so the partitions hold several rows and the sum
        // is grouped by them, which is a different association from row after row.
        run(dir, "families-8-of-3-wide", families(8, 3, 250), NO_FILTERS);
        run(dir, "independent-5-wide", independent(5, 250), NO_FILTERS);
        // A full-rank panel at the largest size this path takes: the request binds.
        run(dir, "independent-25", independent(25, 40), NO_FILTERS);
        // Past dbdsqr, which is divide and conquer, and past LAPACK, which is ARPACK: full rank, so
        // the count does not depend on the solver's bits.
        run(dir, "independent-30", independent(30, 50), and(NO_FILTERS, "--number-of-eigensamples", "30"));
        run(dir, "independent-70-arpack", independent(70, 90), NO_FILTERS);
    }
}
