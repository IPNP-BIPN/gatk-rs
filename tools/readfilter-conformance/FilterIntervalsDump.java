/*
 * FilterIntervals, taken from the reference.
 *
 * The intervals a copy-number panel is allowed to use. Three filters over a mask, and the order and
 * the inclusivity of each is what decides the answer.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE FILTERS SHARE ONE MASK AND RUN IN ORDER. An interval failed by the GC filter is not
 *     considered by the count filters at all, and -- more importantly -- IS NOT IN THE POPULATION
 *     the percentiles are computed over;
 *   - THE ANNOTATION BOUNDS ARE INCLUSIVE, so an interval at exactly the minimum GC content of 0.1
 *     or exactly the maximum of 0.9 survives;
 *   - THE LOW-COUNT FILTER IS STRICTLY GREATER. An interval fails when the number of samples below
 *     the count threshold is `> percentage * numSamples / 100`, so at the default 50 per cent a
 *     single sample out of two is not enough;
 *   - THE PERCENTILE FILTER TAKES ITS PERCENTILES OVER THE SURVIVORS ONLY, per sample, and its
 *     bounds are inclusive as well;
 *   - A PERCENTILE OF ZERO IS SHORT-CIRCUITED TO A THRESHOLD OF ZERO rather than evaluated, which
 *     is what makes `--extreme-count-filter-minimum-percentile 0` mean "no lower bound" instead of
 *     "the smallest value";
 *   - WITH NEITHER ANNOTATIONS NOR COUNTS THE TOOL REFUSES, and with both it intersects them;
 *   - A CONTIG LEFT WITH A SINGLE INTERVAL LOSES IT. After every other filter the tool counts the
 *     survivors per contig and removes any contig's only one, so a run that filters down to one
 *     interval ends with NONE and fails -- which is why the `solitary` case is a refusal;
 *   - AND THE OUTPUT IS A PICARD INTERVAL LIST, `@HD`/`@SQ` and then five tab-separated columns,
 *     not the TSV every other copy-number tool writes.
 *
 * Output:
 *
 *     list\t<label>\t<the whole .interval_list, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FilterIntervalsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.FilterIntervals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FilterIntervalsDump {

    /** The SAM header both input tables carry. */
    static final String HEADER =
            "@HD\tVN:1.6\n"
            + "@SQ\tSN:chr1\tLN:1000\n"
            + "@RG\tID:GATKCopyNumber\tSM:SAMPLE\n";

    /**
     * Five intervals whose GC content walks across the default bounds: below, exactly at the
     * minimum, in the middle, exactly at the maximum, and above.
     */
    static final String ANNOTATED = HEADER
            + "CONTIG\tSTART\tEND\tGC_CONTENT\n"
            + "chr1\t1\t100\t0.050000\n"
            + "chr1\t101\t200\t0.100000\n"
            + "chr1\t201\t300\t0.500000\n"
            + "chr1\t301\t400\t0.900000\n"
            + "chr1\t401\t500\t0.950000\n";

    /** The same five intervals, counted in two samples. */
    static final String COUNTS_ONE = HEADER
            + "CONTIG\tSTART\tEND\tCOUNT\n"
            + "chr1\t1\t100\t5\n"
            + "chr1\t101\t200\t50\n"
            + "chr1\t201\t300\t100\n"
            + "chr1\t301\t400\t150\n"
            + "chr1\t401\t500\t5000\n";

    static final String COUNTS_TWO = HEADER.replace("SM:SAMPLE", "SM:SAMPLE2")
            + "CONTIG\tSTART\tEND\tCOUNT\n"
            + "chr1\t1\t100\t5\n"
            + "chr1\t101\t200\t60\n"
            + "chr1\t201\t300\t110\n"
            + "chr1\t301\t400\t160\n"
            + "chr1\t401\t500\t6000\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filterintervals-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path annotated = write(dir, "annotated.tsv", ANNOTATED);
        final Path countsOne = write(dir, "counts1.tsv", COUNTS_ONE);
        final Path countsTwo = write(dir, "counts2.tsv", COUNTS_TWO);

        System.out.println("# FilterIntervalsDump: the intervals a panel may use");
        // Annotations alone, at the default bounds: the two extremes go, the two on the bounds stay.
        run("annotations", dir, "--annotated-intervals", annotated.toString());
        // Bounds tightened so that ONE interval survives the GC filter. It is then removed as
        // well, because a contig left with a single interval has that interval filtered out, and
        // the run fails with nothing left at all.
        run("solitary", dir, "--annotated-intervals", annotated.toString(),
                "--minimum-gc-content", "0.2", "--maximum-gc-content", "0.8");
        // Bounds that leave two, which survive.
        run("tight-gc", dir, "--annotated-intervals", annotated.toString(),
                "--minimum-gc-content", "0.2", "--maximum-gc-content", "0.92");
        // Counts alone from one sample, at the default low-count threshold of 10.
        run("counts-one-sample", dir, "-I", countsOne.toString());
        // Two samples, where the low-count filter needs MORE than half of them.
        run("counts-two-samples", dir, "-I", countsOne.toString(), "-I", countsTwo.toString());
        // A low-count threshold that catches the middle intervals as well.
        run("high-threshold", dir, "-I", countsOne.toString(),
                "--low-count-filter-count-threshold", "120");
        // The percentile bounds opened, so the extreme interval survives.
        run("wide-percentiles", dir, "-I", countsOne.toString(),
                "--extreme-count-filter-minimum-percentile", "0",
                "--extreme-count-filter-maximum-percentile", "100");
        // Both inputs, which are intersected before anything is filtered.
        run("both", dir, "--annotated-intervals", annotated.toString(), "-I", countsOne.toString());
        // Neither input, which is the tool's own refusal.
        run("neither", dir);
        // And no intervals at all, which is Barclay's.
        run("no-intervals", dir, "--annotated-intervals", annotated.toString());
        // A `-L` naming the whole contig, which equals none of the table's intervals: the
        // intersection is by list equality and leaves nothing at all.
        run("whole-contig", dir, "--annotated-intervals", annotated.toString());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        return path;
    }

    static void run(final String label, final Path dir, final String... extra) throws Exception {
        final Path out = dir.resolve("filtered-" + label + ".interval_list");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-O", out.toString(),
                "--interval-merging-rule", "OVERLAPPING_ONLY"));
        // `-L` is required, and the intersection is `ListUtils.intersection` -- a LIST
        // intersection by equality rather than a genomic one. So `-L chr1` asks for `chr1:1-1000`,
        // which equals none of the table's intervals and leaves nothing: the intervals have to be
        // named exactly. `whole-contig` is the case that records that.
        if (!label.equals("no-intervals")) {
            if (label.equals("whole-contig")) {
                argv.add("-L");
                argv.add("chr1");
            } else {
                for (int start = 1; start <= 401; start += 100) {
                    argv.add("-L");
                    argv.add("chr1:" + start + "-" + (start + 99));
                }
            }
        }
        argv.addAll(Arrays.asList(extra));
        try {
            new FilterIntervals().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("list\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
