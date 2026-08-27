/*
 * ModelSegments' segmentation, taken from the reference.
 *
 * The tool genotypes the allelic counts, segments the denoised copy ratios and the het sites
 * together, and then fits a Markov chain to each segment. The chain is not measured: the
 * multi-sample mode does the genotyping and the segmentation and stops, which is what every run
 * here uses.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE MULTI-SAMPLE MODE WRITES ONE FILE, a Picard interval list, and none of the model's own;
 *   - THE ALLELE-FRACTION SIGNAL DOMINATES THE COPY-RATIO ONE AT THE DEFAULT PENALTY: a fixture
 *     whose copy ratios step at 40000 and 60000 and whose allele fractions step at 30000 and
 *     80000 is cut at 30000 and 80000 alone, and the copy-ratio steps appear only when the
 *     penalty is lowered;
 *   - WITHOUT THE ALLELIC COUNTS THE COPY RATIOS DECIDE, and the same fixture is then cut at
 *     40000;
 *   - --number-of-changepoints-penalty-factor DECIDES HOW MANY BREAKS SURVIVE: a thousand leaves
 *     one segment and a hundredth leaves five, which is every step both signals hold;
 *   - --maximum-number-of-segments-per-chromosome CAPS THEM whatever the penalty says, and the
 *     cap keeps the LAST break rather than the first;
 *   - THE WINDOW SIZES ARE A SET: naming one twice gives the same segmentation as naming it once;
 *   - --minimum-total-allele-count-case DROPS A SITE BEFORE GENOTYPING, and a floor above every
 *     site's total leaves the run with no het at all, so it segments exactly as the copy-ratio
 *     run does;
 *   - --genotyping-homozygous-log-ratio-threshold AT -30 CHANGES NOTHING on a fixture whose het
 *     sites are unambiguous;
 *   - THE INTERVAL LIST CARRIES THE SEQUENCE DICTIONARY of the inputs;
 *   - TWO SAMPLES WHOSE COPY-RATIO INTERVALS DISAGREE ARE REFUSED, by a message that says they
 *     must be identical across all case samples;
 *   - AND THE OUTPUT DIRECTORY IS CREATED IF IT IS NOT THERE.
 *
 * Output:
 *
 *     counts\t<name>=<that input file, escaped>
 *     out\t<label>\t<file name>=<that output file, escaped>
 *     files\t<label>=<the output file names, comma-separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ModelSegmentsDump
 */

import org.broadinstitute.hellbender.tools.copynumber.ModelSegments;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

public class ModelSegmentsDump {

    static final int INTERVALS = 120;
    static final int CONTIG_LENGTH = 400000;

    static int start(final int index) {
        return index * 1000 + 1;
    }

    static int end(final int index) {
        return (index + 1) * 1000;
    }

    /** A denoised copy-ratio file, whose values step at the indices given. */
    static String copyRatios(final String sample, final double baseline, final int[] steps,
                             final double stepSize) {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6");
        lines.add("@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH);
        lines.add("@RG\tID:GATKCopyNumber\tSM:" + sample);
        lines.add("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO");
        for (int i = 0; i < INTERVALS; i++) {
            double value = baseline;
            for (final int step : steps) {
                if (i >= step) {
                    value += stepSize;
                }
            }
            // A small deterministic wobble, so the segmenter sees variance rather than a
            // perfectly flat signal.
            value += ((i * 37) % 7 - 3) * 0.002;
            lines.add(String.format("chr1\t%d\t%d\t%.6f", start(i), end(i), value));
        }
        lines.add("");
        return String.join("\n", lines);
    }

    /** An allelic-count file, whose alternate fraction steps at the indices given. */
    static String allelicCounts(final String sample, final int[] steps, final int total) {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6");
        lines.add("@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH);
        lines.add("@RG\tID:GATKCopyNumber\tSM:" + sample);
        lines.add("CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE");
        for (int i = 0; i < INTERVALS; i++) {
            double fraction = 0.5;
            for (final int step : steps) {
                if (i >= step) {
                    fraction -= 0.15;
                }
            }
            final int alt = (int) Math.round(total * fraction);
            lines.add(String.format("chr1\t%d\t%d\t%d\tA\tC", start(i) + 500, total - alt, alt));
        }
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("model-segments-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ModelSegmentsDump: the genotyping and the segmentation, which is "
                + "all the multi-sample mode does");

        // One copy-ratio step at 40 and one allele-fraction step at 80, so the two signals break
        // the sequence in different places.
        final String ratiosA = copyRatios("sampleA", 0.0, new int[] {40}, 0.6);
        final String countsA = allelicCounts("sampleA", new int[] {80}, 60);
        final Path ratiosAPath = write(dir, "a.cr.tsv", ratiosA);
        final Path countsAPath = write(dir, "a.ac.tsv", countsA);
        System.out.printf("counts\tratios-a=%s%n", ReferenceQueryDump.escape(ratiosA));
        System.out.printf("counts\tcounts-a=%s%n", ReferenceQueryDump.escape(countsA));

        // A second sample whose steps sit elsewhere.
        final String ratiosB = copyRatios("sampleB", 0.1, new int[] {60}, 0.5);
        final String countsB = allelicCounts("sampleB", new int[] {30}, 60);
        final Path ratiosBPath = write(dir, "b.cr.tsv", ratiosB);
        final Path countsBPath = write(dir, "b.ac.tsv", countsB);

        // Two samples, segmented together.
        run(dir, "two-samples", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString()));
        // The same pair with the copy-ratio signal alone.
        run(dir, "copy-ratios-only", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString()));
        // The penalty raised until one segment survives, and lowered.
        run(dir, "penalty-high", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--number-of-changepoints-penalty-factor", "1000"));
        run(dir, "penalty-low", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--number-of-changepoints-penalty-factor", "0.01"));
        // The cap, under a penalty that would otherwise give more.
        run(dir, "capped-segments", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--number-of-changepoints-penalty-factor", "0.01",
                "--maximum-number-of-segments-per-chromosome", "2"));
        // The window sizes as a set: naming one twice changes nothing.
        run(dir, "one-window", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--window-size", "16"));
        run(dir, "one-window-twice", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--window-size", "16", "--window-size", "16"));
        // The genotyping thresholds.
        run(dir, "het-threshold-low", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--genotyping-homozygous-log-ratio-threshold", "-30.0"));
        run(dir, "minimum-count-high", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", ratiosBPath.toString(),
                "--allelic-counts", countsAPath.toString(),
                "--allelic-counts", countsBPath.toString(),
                "--minimum-total-allele-count-case", "1000"));

        // Two samples whose intervals disagree.
        final String shorter = copyRatios("sampleC", 0.0, new int[] {40}, 0.6)
                .replace("chr1\t119001\t120000", "chr1\t119002\t120000");
        final Path shorterPath = write(dir, "c.cr.tsv", shorter);
        run(dir, "mismatched-intervals", List.of(
                "--denoised-copy-ratios", ratiosAPath.toString(),
                "--denoised-copy-ratios", shorterPath.toString()));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label);
        final List<String> argv = new ArrayList<>(List.of(
                "-O", out.toString(),
                "--output-prefix", label));
        argv.addAll(extra);
        try {
            new ModelSegments().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("error\t%s\tno output directory%n", label);
            return;
        }
        final TreeSet<String> names = new TreeSet<>();
        try (final var stream = Files.list(out)) {
            stream.forEach(path -> names.add(path.getFileName().toString()));
        }
        System.out.printf("files\t%s=%s%n", label, String.join(",", names));
        for (final String name : names) {
            System.out.printf("out\t%s\t%s=%s%n", label, name,
                    ReferenceQueryDump.escape(masked(Files.readString(out.resolve(name)), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
