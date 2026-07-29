/*
 * What a FeatureDataSource's lookahead cache returns, taken from the reference.
 *
 * A tool asking "which variants are here" does not get what the file holds at that position: it
 * gets what the cache decided to keep. Two of the cache's behaviours are visible in the answers,
 * and neither is an optimisation a port may skip.
 *
 *   - `refillQueryCache` extends the query end by queryLookaheadBases and caches everything the
 *     reader returned, unfiltered. A later query inside that window is then answered from memory,
 *     and `getCachedFeaturesUpToStopPosition` tests only `start > stopPosition`. So a feature that
 *     ends before the query started is still returned, because nothing tests the end;
 *   - `trimToNewStartPosition` pops the features starting before the new start, keeps those still
 *     overlapping it, and pushes them back **in reverse**, preserving the file's order. A port
 *     that re-sorted would agree on every set and differ on every list.
 *
 * `cacheHit` is containment, not overlap, so a query one base past the cached end is a miss even
 * though nearly everything it needs is in memory. The lookahead is probed at 0 as well as at a
 * positive value, because IntervalWalker.initializeFeatures sets it to 0 deliberately.
 *
 * The features are a hand-written BED file rather than a VCF: the cache is what is being measured,
 * and a codec between the file and the cache would only add a second thing that can differ.
 *
 * Output:
 *
 *     query\t<lookahead>\t<n>\t<interval>\t<hit|miss>\t<returned names>\t<cache contents>\t<cached interval>
 *     trim\t<label>\t<ok|E>\t<cache contents after>
 *
 * Usage: FeatureCacheDump
 */

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class FeatureCacheDump {

    /**
     * A BED body chosen so the cache's edges are reachable: overlapping features, a long one that
     * spans several queries, and a gap.
     *
     * BED is half-open and 0-based, so `chr1 9 20` is chr1:10-20 once decoded.
     */
    static final String BED =
              "chr1\t9\t20\tf1\n"
            + "chr1\t14\t25\tf2\n"
            + "chr1\t19\t120\tf3\n"
            + "chr1\t49\t60\tf4\n"
            + "chr1\t54\t56\tf5\n"
            + "chr1\t99\t110\tf6\n"
            + "chr1\t149\t160\tf7\n"
            + "chr2\t9\t20\tg1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("featurecache-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);
        final Path bed = dir.resolve("features.bed");
        Files.write(bed, BED.getBytes());
        // Random access is a precondition of queryAndPrefetch: without an index the source throws
        // rather than falling back to a scan, and tells the user to run IndexFeatureFile.
        final Index index = IndexFactory.createLinearIndex(bed.toFile(), new BEDCodec());
        index.write(dir.resolve("features.bed.idx").toFile());

        System.out.println("# FeatureCacheDump: what the lookahead cache returns");

        // Zero lookahead is IntervalWalker's setting; 100 is closer to the default of 1,000 while
        // staying inside the fixture.
        for (final int lookahead : new int[] {0, 100}) {
            probeSequence(bed, lookahead);
        }
    }

    static void probeSequence(final Path bed, final int lookahead) throws Exception {
        try (final FeatureDataSource<BEDFeature> source =
                     new FeatureDataSource<>(bed.toString(), null, lookahead, null)) {

            final SimpleInterval[] queries = {
                // Increasing intervals, the pattern the cache is designed for.
                new SimpleInterval("chr1", 10, 12),
                new SimpleInterval("chr1", 11, 13),
                new SimpleInterval("chr1", 15, 16),
                // Inside the previous prefetch, so a hit when the lookahead is positive.
                new SimpleInterval("chr1", 20, 30),
                // Past the prefetch, so a miss either way.
                new SimpleInterval("chr1", 200, 210),
                // Backwards, which the cache is explicitly not designed for.
                new SimpleInterval("chr1", 50, 60),
                new SimpleInterval("chr1", 10, 12),
                // A point query inside a long feature that started well before it.
                new SimpleInterval("chr1", 100, 100),
                // A different contig, which can never be a hit.
                new SimpleInterval("chr2", 10, 12),
                // Off the end of the contig, which a Feature query is allowed to do.
                new SimpleInterval("chr1", 240, 260),
            };

            for (int i = 0; i < queries.length; i++) {
                final SimpleInterval query = queries[i];
                String outcome;
                String names = "";
                try {
                    final List<BEDFeature> features = source.queryAndPrefetch(query);
                    names = describe(features);
                    outcome = "ok";
                } catch (final Exception e) {
                    outcome = "E:" + e.getClass().getName();
                }
                System.out.printf("query\t%d\t%d\t%s\t%s\t%s%n",
                        lookahead, i, query.toString(), outcome, names);
            }
        }
    }

    static String describe(final List<BEDFeature> features) {
        if (features.isEmpty()) {
            return "-";
        }
        final StringBuilder text = new StringBuilder();
        for (final BEDFeature feature : features) {
            if (text.length() > 0) {
                text.append('|');
            }
            text.append(feature.getName())
                .append('@')
                .append(feature.getContig())
                .append(':')
                .append(feature.getStart())
                .append('-')
                .append(feature.getEnd());
        }
        return text.toString();
    }

}
