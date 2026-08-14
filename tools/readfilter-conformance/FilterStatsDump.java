/*
 * FilterStats' summary file, taken from the reference.
 *
 * The `.filteringStats.tsv` that `FilterMutectCalls` writes beside its VCF: three metadata lines,
 * five columns, and one row per filter that accounted for anything. Six behaviours this is built to
 * catch.
 *
 *   - THE TWO ROUNDINGS ARE DIFFERENT. The columns go through `DataLine.set(name, value, 2)` and the
 *     metadata through `roundToNDecimalPlaces(x, 3)`, so the same number is written to two decimals
 *     in the table and three above it;
 *   - AND BOTH OF THEM TURN NaN INTO ZERO, the rounding being
 *     `Math.round((x + ulp(x)) * mult) / mult` and `Math.round` answering `0` for NaN. A run with no
 *     passing call divides by zero three times and writes `0.0` three times;
 *   - WHILE INFINITY BECOMES A VERY LARGE NUMBER rather than `Infinity`: `Math.round` saturates at
 *     `Long.MAX_VALUE`, which divided by the multiplier is `9.223372036854776E15`;
 *   - AND A NEGATIVE ROUNDS THE OTHER WAY, the rounding being `floor(x + 0.5)` after an ulp:
 *     `-1.005` comes out `-1.0` where `1.005` comes out `1.01`, and `-0.0005` comes out `0.0` with
 *     its sign gone;
 *   - THE METADATA IS WRITTEN BEFORE THE COLUMN LINE, the clustering pairs first and then
 *     `threshold`, `fdr` and `sensitivity` in that order;
 *   - THE ROWS ARE THE CALLER'S, in the caller's own order and with no filtering of any kind: the
 *     rule that drops a filter accounting for nothing lives in `FilteringOutputStats` rather than
 *     here;
 *   - AND A NAME THAT NEEDS QUOTING GETS IT, the writer being the same one every GATK table goes
 *     through.
 *
 * Output:
 *
 *     table\t<label>\t<the whole stats file, escaped>
 *
 * Usage: FilterStatsDump
 */

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterStats;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FilterStatsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filter-stats-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FilterStatsDump: three metadata lines, five columns");

        final List<Pair<String, String>> clustering =
                List.of(new ImmutablePair<>("clustering", "1"), new ImmutablePair<>("other", "two"));
        final List<Pair<String, String>> noClustering = new ArrayList<>();

        // An ordinary run: two filters, a threshold and three totals that divide cleanly enough.
        write(dir, "baseline",
                List.of(new FilterStats("weak_evidence", 3.0, 0.15, 2.0, 0.08),
                        new FilterStats("slippage", 1.0, 0.05, 4.0, 0.16)),
                clustering, 0.234, 20.0, 25.0, 3.0, 5.0);

        // Numbers that land on the rounding: the columns keep two decimals and the metadata three.
        write(dir, "rounding",
                List.of(new FilterStats("halves", 1.005, 2.675, 0.12345, 0.005)),
                noClustering, 0.123456, 3.0, 2.0, 1.0, 1.0);

        // Nothing passed: every metadata division is 0/0.
        write(dir, "no-calls", List.of(new FilterStats("weak_evidence", 0.0, 0.0, 0.0, 0.0)),
                noClustering, 0.5, 0.0, 0.0, 0.0, 0.0);

        // Nothing passed but something was expected to be false: a division by zero that is not 0/0.
        write(dir, "infinite-fdr", List.of(new FilterStats("weak_evidence", 1.0, 1.0, 0.0, 0.0)),
                noClustering, 0.5, 0.0, 1.0, 1.0, 0.0);

        // A negative threshold and a negative count, to see which way the rounding leans.
        write(dir, "negatives", List.of(new FilterStats("odd", -1.005, -0.125, -2.5, -0.0)),
                noClustering, -0.0005, 4.0, 2.0, 1.0, 1.0);

        // No rows at all.
        write(dir, "no-rows", new ArrayList<>(), clustering, 0.1, 10.0, 9.0, 1.0, 2.0);

        // A filter name the table writer has to quote, and one with a value-like separator.
        write(dir, "awkward-names",
                List.of(new FilterStats("has\ttab", 1.0, 0.1, 1.0, 0.1),
                        new FilterStats("has\"quote", 1.0, 0.1, 1.0, 0.1),
                        new FilterStats("has,comma", 1.0, 0.1, 1.0, 0.1)),
                noClustering, 0.2, 10.0, 9.0, 1.0, 2.0);
    }

    static void write(final Path dir, final String label, final List<FilterStats> stats,
                      final List<Pair<String, String>> clusteringMetadata, final double threshold,
                      final double totalCalls, final double expectedTruePositives,
                      final double expectedFalsePositives, final double expectedFalseNegatives)
            throws Exception {
        final Path file = dir.resolve(label + ".filteringStats.tsv");
        FilterStats.writeM2FilterSummary(stats, file, clusteringMetadata, threshold, totalCalls,
                expectedTruePositives, expectedFalsePositives, expectedFalseNegatives);
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
