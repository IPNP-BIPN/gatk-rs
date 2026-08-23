/*
 * PrintReadCounts's output, taken from the reference.
 *
 * A depth-evidence file (.rd.txt), or a counts file (.counts.tsv), rewritten as one counts file per
 * sample for the CNV callers. It is a `FeatureWalker`, and the two inputs it accepts arrive with two
 * different kinds of header, which is where nearly everything measurable lives.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - ONE OUTPUT FILE PER SAMPLE, named by RAW CONCATENATION: `outputPrefix + sampleName +
 *     ".counts.tsv"`, so a prefix that does not end in a separator glues onto the sample name;
 *   - AN .rd.txt HEADER CARRIES SAMPLE NAMES BUT A NULL DICTIONARY, so the run needs
 *     --sequence-dictionary and refuses without one, with a message that misspells the argument:
 *     "No dictionary available.  Supply one with --sequence-dictonary.";
 *   - THE OUTPUT IS ONE-BASED WHERE THE INPUT IS ZERO-BASED: the codec adds one to the start on
 *     read and the tool writes `getStart()` rather than re-encoding, so a bin written `0 100`
 *     comes back `1 100`;
 *   - THE OUTPUT HEADER IS BUILT, NOT COPIED: a fresh SAMFileHeader over the dictionary plus one
 *     read group whose ID is always GATKCopyNumber and whose SM is the sample, then the column line
 *     CONTIG START END COUNT;
 *   - SO A .counts.tsv INPUT LOSES EVERYTHING ELSE ITS HEADER CARRIED, its @PG, its @CO and its
 *     read group's own ID among them, while its records pass through unchanged and the dictionary
 *     used is its own rather than the one --sequence-dictionary names;
 *   - A .counts.tsv NAMES ITS SAMPLE THROUGH ITS READ GROUPS, and the refusal happens while the
 *     feature reader parses the header, so it reaches the caller WRAPPED TWICE, as a GATKException
 *     "Error initializing feature reader for path ..." over a TribbleException.MalformedFeatureFile
 *     "Unable to parse header with error: ..., for input source: ...";
 *   - AND A READ GROUP WITH NO SM DOES NOT PRODUCE readSampleName's "does not contain a sample
 *     name": the distinct list holds one null, which passes the emptiness test, and the refusal
 *     comes later from Utils.nonEmpty as "The string is null: string must not be null or empty",
 *     so that message is unreachable from a header that has read groups at all;
 *   - --output-file-list WRITES <sample>\t<filename> WITH THE SAME CONCATENATED NAME, and is closed
 *     before a single record is written;
 *   - A DUPLICATED COLUMN NAME IN AN .rd.txt HEADER IS NOT REFUSED: two writers open the same path,
 *     the list names it twice, and one file survives carrying the second writer's counts;
 *   - A RECORD WITH FEWER COUNTS THAN THE HEADER HAS SAMPLES IS AN ArrayIndexOutOfBoundsException
 *     rather than a UserException, and one with more counts silently drops the extra;
 *   - WHAT THAT CRASH LEAVES BEHIND IS HALF A HEADER: SAMTextHeaderCodec flushes as it encodes, so
 *     the @HD, @SQ and @RG lines are on disk, while the tool's own column line is still in the
 *     BufferedWriter and dies with it;
 *   - -L SUBSETS THE RECORDS BUT NOT THE FILES, every sample still getting one;
 *   - AND ANY OTHER FEATURE TYPE IS REFUSED BY FeatureWalker BEFORE THE TOOL RUNS, with
 *     "contains features of the wrong type".
 *
 * Output:
 *
 *     input\t<label>\t<name>=<the whole file, escaped>
 *     out\t<label>\t<file name>=<the whole file, escaped>
 *     list\t<label>=<the whole output file list, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PrintReadCountsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.sv.PrintReadCounts;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class PrintReadCountsDump {

    /** One depth bin, written zero-based and half-open as an .rd.txt holds it. */
    static String bin(final String contig, final int start, final int end, final int... counts) {
        final StringBuilder line = new StringBuilder(contig + "\t" + start + "\t" + end);
        for (final int count : counts) {
            line.append('\t').append(count);
        }
        return line.append('\n').toString();
    }

    static String rdHeader(final String... samples) {
        return "#Chr\tStart\tEnd\t" + String.join("\t", samples) + "\n";
    }

    /** One count, written one-based and closed as a .counts.tsv holds it. */
    static String count(final String contig, final int start, final int end, final int count) {
        return contig + "\t" + start + "\t" + end + "\t" + count + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("print-read-counts-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintReadCountsDump: depth evidence rewritten as one counts file per sample");

        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "counts", List.of("chr1", "chr2"));
        final String samHeader = Files.readString(dict);

        // Two samples over two bins, which is the shape the tool is for.
        final String twoSamples = rdHeader("alpha", "beta")
                + bin("chr1", 0, 100, 11, 21)
                + bin("chr1", 100, 200, 12, 22)
                + bin("chr2", 0, 100, 13, 23);
        run(dir, "rd-two-samples", "rd.txt", twoSamples, dict, true);

        // The same input with no dictionary at all, which an .rd.txt header cannot supply.
        run(dir, "rd-no-dictionary", "rd.txt", twoSamples, null, false);

        // A prefix that does not end in a separator, which is concatenated raw.
        runPrefixed(dir, "rd-glued-prefix", "rd.txt", twoSamples, dict, "sample-", false);

        // One sample, to show the smallest output whole.
        run(dir, "rd-one-sample", "rd.txt",
                rdHeader("solo") + bin("chr1", 0, 100, 7), dict, false);

        // The same column name twice: two writers, one path.
        run(dir, "rd-duplicate-sample", "rd.txt",
                rdHeader("twin", "twin") + bin("chr1", 0, 100, 11, 22), dict, true);

        // A record with fewer counts than the header has samples, and one with more.
        run(dir, "rd-short-record", "rd.txt",
                rdHeader("alpha", "beta") + bin("chr1", 0, 100, 11), dict, false);
        run(dir, "rd-long-record", "rd.txt",
                rdHeader("alpha") + bin("chr1", 0, 100, 11, 99), dict, false);

        // Intervals, which drop records but never files.
        runIntervals(dir, "rd-intervals", "rd.txt", twoSamples, dict, "chr2");

        // A counts file, whose header names its own sample and its own dictionary, and whose
        // records are already one-based.
        final String countsIn = samHeader
                + "@RG\tID:not-the-cnv-id\tSM:gamma\n"
                + "@PG\tID:something\tPN:something\n"
                + "@CO\ta comment the output will not carry\n"
                + "CONTIG\tSTART\tEND\tCOUNT\n"
                + count("chr1", 1, 100, 31)
                + count("chr1", 101, 200, 32);
        run(dir, "counts-round-trip", "counts.tsv", countsIn, dict, true);

        // Two read groups naming two samples, and a read group naming none.
        run(dir, "counts-two-samples", "counts.tsv",
                samHeader + "@RG\tID:one\tSM:gamma\n@RG\tID:two\tSM:delta\n"
                        + "CONTIG\tSTART\tEND\tCOUNT\n" + count("chr1", 1, 100, 31),
                dict, false);
        run(dir, "counts-no-sample", "counts.tsv",
                samHeader + "@RG\tID:one\n" + "CONTIG\tSTART\tEND\tCOUNT\n" + count("chr1", 1, 100, 31),
                dict, false);

        // A counts file whose header does not come first.
        run(dir, "counts-header-late", "counts.tsv",
                "CONTIG\tSTART\tEND\tCOUNT\n" + samHeader + count("chr1", 1, 100, 31), dict, false);

        // Another SV feature type entirely, which never reaches the tool.
        run(dir, "wrong-feature-type", "baf.txt",
                "chr1\t100\t0.50\talpha\n", dict, false);
    }

    static void run(final Path dir, final String label, final String extension, final String text,
                    final Path dictionary, final boolean withList) throws Exception {
        runPrefixed(dir, label, extension, text, dictionary, "", withList);
    }

    static void runPrefixed(final Path dir, final String label, final String extension,
                            final String text, final Path dictionary, final String namePrefix,
                            final boolean withList) throws Exception {
        execute(dir, label, extension, text, dictionary, namePrefix, withList, new String[0]);
    }

    static void runIntervals(final Path dir, final String label, final String extension,
                             final String text, final Path dictionary, final String interval)
            throws Exception {
        execute(dir, label, extension, text, dictionary, "", false, new String[] {"-L", interval});
    }

    static void execute(final Path dir, final String label, final String extension,
                        final String text, final Path dictionary, final String namePrefix,
                        final boolean withList, final String[] extra) throws Exception {
        final Path work = dir.resolve(label);
        Files.createDirectories(work);
        final Path input = work.resolve("input." + extension);
        Files.writeString(input, text, StandardCharsets.UTF_8);
        try {
            new IndexFeatureFile().instanceMain(new String[] {"-I", input.toString()});
        } catch (final Exception e) {
            System.out.printf("index\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
        }
        System.out.printf("input\t%s\tinput.%s=%s%n", label, extension,
                ReferenceQueryDump.escape(masked(text, dir)));

        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--input-counts", input.toString(),
                "--output-prefix", work.toString() + "/" + namePrefix));
        if (dictionary != null) {
            argv.addAll(Arrays.asList("--sequence-dictionary", dictionary.toString()));
        }
        final Path list = work.resolve("outputs.list");
        if (withList) {
            argv.addAll(Arrays.asList("--output-file-list", list.toString()));
        }
        argv.addAll(Arrays.asList(extra));

        Exception failure = null;
        try {
            new PrintReadCounts().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception e) {
            failure = e;
        }
        if (failure != null) {
            System.out.printf("error\t%s\t%s:%s%n", label, failure.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(failure.getMessage()), dir)));
            for (Throwable cause = failure.getCause(); cause != null; cause = cause.getCause()) {
                System.out.printf("cause\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            }
        }
        try (final Stream<Path> written = Files.list(work)) {
            final List<Path> counts = written
                    .filter(p -> p.getFileName().toString().endsWith(".counts.tsv"))
                    .filter(p -> !p.equals(input))
                    .sorted(Comparator.comparing(p -> p.getFileName().toString()))
                    .toList();
            for (final Path path : counts) {
                System.out.printf("out\t%s\t%s=%s%n", label, path.getFileName(),
                        ReferenceQueryDump.escape(masked(Files.readString(path), dir)));
            }
        }
        if (Files.exists(list)) {
            System.out.printf("list\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(list), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
