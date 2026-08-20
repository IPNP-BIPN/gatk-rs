/*
 * SplitIntervals, taken from the reference.
 *
 * The scatter modes are already pinned by `interval-list-scatter`; this is the tool around them:
 * which intervals it starts from, how many files it writes, and what each is called.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE FILE NAME'S WIDTH COMES FROM A LOGARITHM OF `scatterCount - 1`. At a scatter count of
 *     one that is `log10(0)`, which is negative infinity: the cast floors it to
 *     `Integer.MIN_VALUE`, the `+ 1` is still enormously negative, and the `max` with
 *     `--interval-file-num-digits` is what saves the name. So the default four digits are used
 *     for every scatter count up to ten thousand and only a larger one widens the name;
 *   - NO INTERVALS AT ALL MEANS THE WHOLE REFERENCE, filtered by `--min-contig-size`, which is a
 *     filter on the CONTIG and not on the intervals: with `-L` given, the size is ignored;
 *   - `--dont-mix-contigs` SPLITS AFTER THE SCATTER, so the requested count is a lower bound and
 *     the files can outnumber it;
 *   - THAT SPLIT ORDERS THE SUBLISTS BY THE CONTIG'S DICTIONARY INDEX, not by the order the
 *     intervals arrived in;
 *   - THE MERGING RULE IS THE COMMON DEFAULT, which is ALL, so adjacent input intervals are merged
 *     before anything is scattered unless OVERLAPPING_ONLY is asked for;
 *   - THE PREFIX AND EXTENSION ARE CONCATENATED RAW, so a prefix with no separator and an
 *     extension with no dot both land exactly as given;
 *   - AND A SCATTER COUNT OF ZERO IS THE TOOL'S OWN REFUSAL while a negative number of digits is
 *     the parser's.
 *
 * Output:
 *
 *     files\t<label>=<name>,<name>,...
 *     list\t<label>,<name>=<the whole interval list, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * The M5 and UR fields of the sequence lines are masked, as in PreprocessIntervalsDump: both come
 * from the .dict the reference indexer wrote and neither is this tool's.
 *
 * Usage: SplitIntervalsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.SplitIntervals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SplitIntervalsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("split-intervals-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# SplitIntervalsDump: the shards a scatter-gather run is given");

        // The whole reference, scattered every way.
        run("whole-genome-1", dir, fasta, "--scatter-count", "1");
        run("whole-genome-3", dir, fasta, "--scatter-count", "3");
        run("whole-genome-5", dir, fasta, "--scatter-count", "5");
        // A contig filter, which only applies when no intervals were given.
        run("min-contig-size", dir, fasta, "--scatter-count", "2", "--min-contig-size", "241");
        run("min-contig-size-with-intervals", dir, fasta, "--scatter-count", "2",
                "--min-contig-size", "241", "-L", "chr1:1-100");
        // Explicit intervals under each mode, at a count the modes disagree about.
        for (final String mode : new String[] {"INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
                "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW", "INTERVAL_COUNT",
                "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER"}) {
            run("mode-" + mode, dir, fasta, "--scatter-count", "3", "--subdivision-mode", mode,
                    "-L", "chr1:1-100", "-L", "chr1:141-240", "-L", "chr2:1-50",
                    "--interval-merging-rule", "OVERLAPPING_ONLY");
        }
        // Adjacent intervals under the common default merging rule, which merges them, and under
        // OVERLAPPING_ONLY, which does not.
        run("adjacent-merged", dir, fasta, "--scatter-count", "2", "-L", "chr1:1-50",
                "-L", "chr1:51-100");
        run("adjacent-kept", dir, fasta, "--scatter-count", "2", "-L", "chr1:1-50",
                "-L", "chr1:51-100", "--interval-merging-rule", "OVERLAPPING_ONLY");
        // The contig split, which happens after the scatter and can produce more files than asked.
        run("dont-mix-contigs", dir, fasta, "--scatter-count", "2", "--dont-mix-contigs",
                "-L", "chr1:1-100", "-L", "chr2:1-100", "-L", "chr2:201-240",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        // The same intervals given second contig first, to show the sublists ordered by the
        // dictionary rather than by arrival.
        run("dont-mix-contigs-reversed", dir, fasta, "--scatter-count", "1", "--dont-mix-contigs",
                "-L", "chr2:1-100", "-L", "chr1:1-100",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        // The naming: a prefix, an extension without a dot, and every width the formula can reach.
        run("named", dir, fasta, "--scatter-count", "2", "--interval-file-prefix", "shard-",
                "--extension", ".list", "-L", "chr1:1-100");
        run("one-digit", dir, fasta, "--scatter-count", "2", "--interval-file-num-digits", "1",
                "-L", "chr1:1-100");
        run("one-digit-twelve", dir, fasta, "--scatter-count", "12",
                "--interval-file-num-digits", "1", "-L", "chr1:1-240");
        run("one-digit-one", dir, fasta, "--scatter-count", "1", "--interval-file-num-digits", "1",
                "-L", "chr1:1-100");
        run("eight-digits", dir, fasta, "--scatter-count", "2", "--interval-file-num-digits", "8",
                "-L", "chr1:1-100");
        // A scatter count larger than the intervals can honour, which comes back short.
        run("more-shards-than-intervals", dir, fasta, "--scatter-count", "10",
                "--subdivision-mode", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
                "-L", "chr1:1-100", "-L", "chr1:201-240",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        // The two refusals.
        run("zero-scatter-count", dir, fasta, "--scatter-count", "0", "-L", "chr1:1-100");
        run("zero-digits", dir, fasta, "--scatter-count", "2", "--interval-file-num-digits", "0",
                "-L", "chr1:1-100");
    }

    static void run(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label);
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new SplitIntervals().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        final List<Path> written = new ArrayList<>();
        try (final java.util.stream.Stream<Path> stream = Files.list(out)) {
            stream.sorted().forEach(written::add);
        }
        final List<String> names = new ArrayList<>();
        for (final Path path : written) {
            names.add(path.getFileName().toString());
        }
        System.out.printf("files\t%s=%s%n", label, String.join(",", names));
        for (final Path path : written) {
            System.out.printf("list\t%s,%s=%s%n", label, path.getFileName(),
                    ReferenceQueryDump.escape(PreprocessIntervalsDump.masked(
                            new String(Files.readAllBytes(path)))));
        }
    }
}
