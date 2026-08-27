/*
 * SeriesStats, taken from the reference.
 *
 * The accumulator the ground-truth tools keep their numbers in. It is a sorted map from value to
 * count rather than a list, and every statistic it reports is computed off that map, which is what
 * makes several of them not the statistic their name suggests.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A PERCENTILE IS AN OBSERVED VALUE, never an interpolation: the walk returns a BIN KEY, so the
 *     median of an even-sized set is one of its members rather than the mean of two;
 *   - THE PERCENTILE INDEX IS TRUNCATED, `(int)(count * p / 100)`, so the 50th of four values is
 *     the third smallest rather than the second;
 *   - A SINGLE VALUE SHORT-CIRCUITS to the LAST added rather than to the only bin, which is the
 *     same thing until it is not;
 *   - AN EMPTY SERIES REPORTS NaN for every statistic and zero for the counts;
 *   - THE STANDARD DEVIATION DIVIDES BY THE COUNT, not by the count less one, so it is the
 *     population deviation;
 *   - THE BINS ARE A TreeMap OVER Double, so -0.0 SORTS BEFORE 0.0 and the two are separate bins,
 *     while a NaN sorts after everything;
 *   - `getUniq` IS THE NUMBER OF BINS, so those two zeros count as two;
 *   - THE CSV IS WRITTEN AS INTEGERS only when every value arrived through `add(int)`, and one
 *     double is enough to switch the whole file to `%f`;
 *   - AND `add(int)` CALLS `add(double)` FIRST, so the integer count is incremented after the
 *     value is already binned.
 *
 * Output:
 *
 *     stat\t<label>=<count>,<uniq>,<min>,<max>,<mean>,<median>,<std>,<last>
 *     pct\t<label>=<the requested percentiles, comma separated>
 *     bins\t<label>=<key:count, comma separated, in map order>
 *     csv\t<label>=<the whole csv, escaped>
 *     digest\t<label>=<toDigest>
 *
 * Usage: SeriesStatsDump
 */

import org.broadinstitute.hellbender.tools.walkers.groundtruth.SeriesStats;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class SeriesStatsDump {

    static String d(final double value) {
        return String.format("%.10f", value);
    }

    static void report(final String label, final SeriesStats stats) {
        System.out.printf("stat\t%s=%d,%d,%s,%s,%s,%s,%s,%s%n", label,
                stats.getCount(), stats.getUniq(), d(stats.getMin()), d(stats.getMax()),
                d(stats.getMean()), d(stats.getMedian()), d(stats.getStd()), d(stats.getLast()));
        final List<String> percentiles = new ArrayList<>();
        for (final double p : new double[] {0, 10, 25, 50, 75, 90, 99, 100}) {
            percentiles.add(d(stats.getPercentile(p)));
        }
        System.out.printf("pct\t%s=%s%n", label, String.join(",", percentiles));
        final List<String> bins = new ArrayList<>();
        for (final Map.Entry<Double, AtomicInteger> entry : stats.getBins().entrySet()) {
            bins.add(entry.getKey() + ":" + entry.getValue().get());
        }
        System.out.printf("bins\t%s=%s%n", label, String.join(",", bins));
        System.out.printf("digest\t%s=%s%n", label, stats.toDigest());
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("series-stats-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SeriesStatsDump: the accumulator the ground-truth tools keep their "
                + "numbers in");

        // Nothing at all.
        report("empty", new SeriesStats());

        // One value, which short-circuits every percentile.
        final SeriesStats one = new SeriesStats();
        one.add(7);
        report("one", one);

        // An EVEN number of values whose middle two differ, which is where an interpolating median
        // and this one part company.
        final SeriesStats even = new SeriesStats();
        for (final int v : new int[] {1, 2, 8, 9}) {
            even.add(v);
        }
        report("even", even);

        // An odd number.
        final SeriesStats odd = new SeriesStats();
        for (final int v : new int[] {1, 2, 8, 9, 100}) {
            odd.add(v);
        }
        report("odd", odd);

        // Repeats, so a bin carries more than one and the walk has to step over it.
        final SeriesStats repeated = new SeriesStats();
        for (final int v : new int[] {5, 5, 5, 5, 5, 5, 5, 5, 5, 100}) {
            repeated.add(v);
        }
        report("repeated", repeated);

        // The two zeros, which a TreeMap over Double keeps apart.
        final SeriesStats zeros = new SeriesStats();
        zeros.add(0.0);
        zeros.add(-0.0);
        zeros.add(1.0);
        report("zeros", zeros);

        // A NaN, which sorts after everything.
        final SeriesStats withNan = new SeriesStats();
        withNan.add(1.0);
        withNan.add(Double.NaN);
        withNan.add(2.0);
        report("nan", withNan);

        // Doubles, whose deviation is not a whole number.
        final SeriesStats doubles = new SeriesStats();
        for (final double v : new double[] {1.5, 2.25, 2.25, 10.0}) {
            doubles.add(v);
        }
        report("doubles", doubles);

        // The CSV, integer keys and not.
        final Path intCsv = dir.resolve("int.csv");
        even.csvWrite(intCsv.toString());
        System.out.printf("csv\tint=%s%n",
                ReferenceQueryDump.escape(Files.readString(intCsv)));
        final Path doubleCsv = dir.resolve("double.csv");
        doubles.csvWrite(doubleCsv.toString());
        System.out.printf("csv\tdouble=%s%n",
                ReferenceQueryDump.escape(Files.readString(doubleCsv)));
        // One double among the integers, which switches the WHOLE file.
        final SeriesStats mixed = new SeriesStats();
        mixed.add(1);
        mixed.add(2);
        mixed.add(3.5);
        final Path mixedCsv = dir.resolve("mixed.csv");
        mixed.csvWrite(mixedCsv.toString());
        System.out.printf("csv\tmixed=%s%n",
                ReferenceQueryDump.escape(Files.readString(mixedCsv)));
        report("mixed", mixed);
        // And a series of integers whose VALUES are whole but which arrived as doubles.
        final SeriesStats wholeDoubles = new SeriesStats();
        wholeDoubles.add(1.0);
        wholeDoubles.add(2.0);
        final Path wholeCsv = dir.resolve("whole.csv");
        wholeDoubles.csvWrite(wholeCsv.toString());
        System.out.printf("csv\twhole=%s%n",
                ReferenceQueryDump.escape(Files.readString(wholeCsv)));
    }
}
