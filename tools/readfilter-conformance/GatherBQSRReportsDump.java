/*
 * GatherBQSRReports, taken from the reference.
 *
 * The gather that ends a scattered BQSR run: several recalibration reports in, one out, with every
 * observation and every error summed and the quantization recomputed from the total.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE READ GROUPS ARE THE UNION OF EVERY INPUT'S, gathered before anything is combined, so a
 *     shard that never saw a read group still carries its rows in the output;
 *   - THE COMBINE DOES NOT RECALCULATE THE QUANTIZATION, the tool does, from the summed tables and
 *     the FIRST report's quantizing level count;
 *   - THE EMPIRICAL QUALITY IS RECALCULATED, so a gathered row is not the sum of the shards'
 *     empirical qualities but the quality of the summed counts;
 *   - THE ARGUMENTS TABLE IS THE FIRST INPUT'S, so the gathered report names that shard's
 *     `--input` and its own known sites;
 *   - AN EMPTY SHARD IS SKIPPED BY THE COMBINE, and a gather of nothing but empty shards is a
 *     refusal saying there is no usable data;
 *   - THE SAME FILE TWICE IS READ TWICE, the inputs being a list, so the counts double;
 *   - AND THE OUTPUT IS A GATKReport OF FIVE TABLES, printed by the same writer BaseRecalibrator
 *     uses, so a gather of one shard is that shard requantized rather than a copy.
 *
 * Output:
 *
 *     shard\t<label>=<the whole recalibration table of one scatter, escaped>
 *     gathered\t<label>=<the whole gathered table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GatherBQSRReportsDump
 */

import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.tools.walkers.bqsr.GatherBQSRReports;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GatherBQSRReportsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("gather-bqsr-dump");
        BaseRecalibratorDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = BaseRecalibratorDump.writeReference(dir);
        final Path bam = dir.resolve("input.bam");
        BaseRecalibratorDump.buildFixture(bam.toFile());

        final Path sites = dir.resolve("sites.bed");
        Files.writeString(sites, "chr1\t9\t12\n", StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", sites.toString()});

        System.out.println("# GatherBQSRReportsDump: the gather that ends a scattered BQSR run");

        // Two shards of one run, split by interval, and one interval with no reads at all.
        final Path first = recalibrate(dir, bam, fasta, sites, "first", "chr1:1-20");
        final Path second = recalibrate(dir, bam, fasta, sites, "second", "chr1:21-95");
        final Path empty = recalibrate(dir, bam, fasta, sites, "empty", "chr1:60-95");

        gather("two-shards", dir, first, second);
        // The shards the other way round, which decides whose arguments table survives.
        gather("reversed", dir, second, first);
        // One shard, which is that shard requantized rather than copied.
        gather("one-shard", dir, first);
        // The same shard twice, which doubles every count.
        gather("same-twice", dir, first, first);
        // An empty shard beside a real one, which the combine skips.
        gather("with-empty", dir, first, empty);
        // Nothing but empty shards, which is the refusal.
        gather("all-empty", dir, empty, empty);
    }

    /** One scatter's recalibration table, and the row that carries it. */
    static Path recalibrate(final Path dir, final Path bam, final Path fasta, final Path sites,
                            final String label, final String interval) throws Exception {
        final Path output = dir.resolve("recal." + label + ".table");
        new BaseRecalibrator().instanceMain(new String[] {
                "-I", bam.toString(), "-R", fasta.toString(), "--known-sites", sites.toString(),
                "-L", interval, "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"});
        System.out.printf("shard\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(output)));
        return output;
    }

    static void gather(final String label, final Path dir, final Path... inputs) throws Exception {
        final Path output = dir.resolve("gathered." + label + ".table");
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.addAll(Arrays.asList("-I", input.toString()));
        }
        argv.addAll(Arrays.asList("-O", output.toString()));
        try {
            new GatherBQSRReports().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("gathered\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(output)));
    }
}
