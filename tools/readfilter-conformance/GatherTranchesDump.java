/*
 * GatherTranches, taken from the reference.
 *
 * The gather that ends a scattered VariantRecalibrator run: several VQSLOD tranche files in, one
 * truth-sensitivity tranche file out, with the counts summed per VQSLOD level and the requested
 * sensitivities matched to whichever merged level comes closest.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE SHARDS ARE POOLED INTO A TreeMap BY minVQSLod and walked in DESCENDING key order, so the
 *     merged list is ordered by VQSLOD however the files were given;
 *   - THE Ti/Tv OF A MERGED LEVEL IS A RATIO OF SUMS, not a mean of ratios: each shard's ratio is
 *     turned back into transitions and transversions, those are summed, and the merged ratio is
 *     their quotient;
 *   - A SHARD WITH NO KNOWN VARIANTS CONTRIBUTES A ZERO TO BOTH SUMS, so a merged level whose
 *     shards are all empty has a Ti/Tv of NaN, written as `NaN`;
 *   - THE REQUESTED SENSITIVITIES ARE SORTED IN PLACE and matched by a walk that adds the PREVIOUS
 *     tranche when the distance to the target grows, so the answer is the last level before the
 *     one that overshot;
 *   - THE FIRST MERGED TRANCHE IS NEVER A CANDIDATE on its own: it is consumed as the initial
 *     `currentTranche`, and only becomes an answer as somebody's `prevTranche`;
 *   - THE WALK STOPS AT THE FIRST TARGET IT CANNOT ADVANCE PAST, so asking for more sensitivities
 *     than there are levels writes fewer rows than asked;
 *   - A FILE OF A DIFFERENT VERSION IS REFUSED by the version line and not by its columns;
 *   - TWO SHARDS WHOSE VQSLOD LEVELS DIFFER ARE STILL POOLED, each level merged from whatever
 *     shards carry it, so a level only one file has is merged from one;
 *   - AND THE OUTPUT IS A TRUTH-SENSITIVITY FILE, whose header says VERSION 5 while the inputs
 *     must say version 6, whose first column is named `targetTruthSensitivity`, and whose rows
 *     carry the REQUESTED sensitivity rather than the achieved one.
 *
 * Output:
 *
 *     tranches\t<label>=<the whole gathered file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GatherTranchesDump
 */

import org.broadinstitute.hellbender.tools.walkers.vqsr.GatherTranches;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GatherTranchesDump {

    static final String HEADER =
            "# Variant quality score tranches file\n"
            + "# Version number 6\n"
            + "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,"
            + "accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";

    /** One row of a VQSLOD tranche file. */
    static String row(final double requested, final long known, final long novel,
                      final double knownTiTv, final double novelTiTv, final double minVQSLod,
                      final String model, final int accessible, final int called) {
        final double sensitivity = accessible == 0 ? 0.0 : (double) called / accessible;
        return String.format("%.4f,%d,%d,%.4f,%.4f,%.4f,VQSRTranche,%s,%d,%d,%.4f%n",
                requested, known, novel, knownTiTv, novelTiTv, minVQSLod, model, accessible,
                called, sensitivity);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("gather-tranches-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GatherTranchesDump: the gather that ends a scattered VQSR run");

        // Two shards carrying the same four VQSLOD levels, with different counts.
        final String firstShard = HEADER
                + row(4.0, 100, 20, 2.0, 1.5, 4.0, "SNP", 1000, 500)
                + row(2.0, 200, 50, 2.1, 1.6, 2.0, "SNP", 1000, 800)
                + row(0.0, 300, 90, 2.2, 1.7, 0.0, "SNP", 1000, 950)
                + row(-2.0, 400, 150, 2.3, 1.8, -2.0, "SNP", 1000, 990);
        final String secondShard = HEADER
                + row(4.0, 60, 10, 1.8, 1.4, 4.0, "SNP", 1000, 450)
                + row(2.0, 130, 30, 1.9, 1.5, 2.0, "SNP", 1000, 780)
                + row(0.0, 220, 70, 2.0, 1.6, 0.0, "SNP", 1000, 940)
                + row(-2.0, 330, 120, 2.1, 1.7, -2.0, "SNP", 1000, 985);
        // A shard carrying one level the others do not.
        final String extraLevel = HEADER
                + row(4.0, 10, 2, 2.0, 1.5, 4.0, "SNP", 1000, 400)
                + row(1.0, 90, 25, 2.0, 1.5, 1.0, "SNP", 1000, 700);
        // A shard whose known counts are all zero, which is the NaN ratio.
        final String noKnown = HEADER
                + row(4.0, 0, 20, 0.0, 1.5, 4.0, "SNP", 1000, 500)
                + row(2.0, 0, 50, 0.0, 1.6, 2.0, "SNP", 1000, 800);
        // A file of the wrong version.
        final String oldVersion = HEADER.replace("Version number 6", "Version number 5")
                + row(4.0, 100, 20, 2.0, 1.5, 4.0, "SNP", 1000, 500);

        final Path first = write(dir, "first.tranches", firstShard);
        final Path second = write(dir, "second.tranches", secondShard);
        final Path extra = write(dir, "extra.tranches", extraLevel);
        final Path zeroes = write(dir, "no-known.tranches", noKnown);
        final Path old = write(dir, "old-version.tranches", oldVersion);

        // The default tranche levels over two shards.
        run("two-shards", dir, "SNP", new Path[] {first, second});
        // The same two the other way round, which the TreeMap makes irrelevant.
        run("reversed", dir, "SNP", new Path[] {second, first});
        // One shard alone.
        run("one-shard", dir, "SNP", new Path[] {first});
        // A shard carrying a level the other does not.
        run("extra-level", dir, "SNP", new Path[] {first, extra});
        // The zero-known shard, whose merged Ti/Tv is NaN.
        run("no-known", dir, "SNP", new Path[] {zeroes});
        // Fewer requested sensitivities than there are levels.
        run("one-level", dir, "SNP", new Path[] {first, second}, "--truth-sensitivity-tranche", "99.0");
        // More requested than there are levels.
        run("many-levels", dir, "SNP", new Path[] {first, second},
                "--truth-sensitivity-tranche", "100.0", "--truth-sensitivity-tranche", "99.9",
                "--truth-sensitivity-tranche", "99.5", "--truth-sensitivity-tranche", "99.0",
                "--truth-sensitivity-tranche", "98.0", "--truth-sensitivity-tranche", "95.0",
                "--truth-sensitivity-tranche", "90.0");
        // The levels given out of order, which the tool sorts.
        run("unsorted-levels", dir, "SNP", new Path[] {first, second},
                "--truth-sensitivity-tranche", "90.0", "--truth-sensitivity-tranche", "99.9");
        // The other two modes, which only change the model column.
        run("indel-mode", dir, "INDEL", new Path[] {first, second});
        run("both-mode", dir, "BOTH", new Path[] {first, second});
        // The version refusal.
        run("old-version", dir, "SNP", new Path[] {old});
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final String label, final Path dir, final String mode, final Path[] inputs,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("gathered-" + label + ".tranches");
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.addAll(Arrays.asList("-I", input.toString()));
        }
        argv.addAll(Arrays.asList("--mode", mode, "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new GatherTranches().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("tranches\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
