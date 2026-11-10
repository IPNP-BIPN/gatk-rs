/*
 * ConvertCountsToDepthFile's output, taken from the reference (new in GATK 4.7.0.0).
 *
 * A read-counts file (.counts.tsv) rewritten as a one-sample depth-evidence file (.rd.txt). It is a
 * `FeatureWalker` over `SimpleCount`, so the counts file arrives with its whole SAM header, and the
 * output goes through `FeatureOutputCodecFinder` and `FeatureOutputStream`, which choose the format
 * and the compression from the output's NAME.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE COORDINATES CHANGE BASE: a .counts.tsv is one-based and closed, and `DepthEvidenceCodec`
 *     writes `getStart() - 1`, so a count `1 100` comes out `0 100`;
 *   - THE HEADER IS ONE COLUMN LINE NAMING ONE SAMPLE, `--sample-name` (or `-n`) when given and
 *     the counts file's read-group SM otherwise, whatever else that header carried;
 *   - THE TOOL'S OWN "Sample name is not provided" IS UNREACHABLE: the counts codec builds a
 *     `SimpleSampleLocatableMetadata`, whose `Utils.nonEmpty` refuses a header without a sample
 *     while the feature reader is still being opened, so the refusal arrives WRAPPED TWICE
 *     (GATKException over TribbleException.MalformedFeatureFile) and `--sample-name` cannot
 *     rescue it;
 *   - A HEADER WITH NO READ GROUP IS REFUSED THE SAME WAY, by `MetadataUtils.readSampleName`;
 *   - A BLOCK-COMPRESSED NAME (.gz) IS BGZF AT --compression-level, WITH A TABIX INDEX BESIDE IT
 *     built over the counts header's dictionary, and a plain name ignores the level altogether;
 *   - --compression-level IS RANGE-CHECKED BY THE PARSER, [0, 9];
 *   - AN OUTPUT NAME NO CODEC CLAIMS IS THE FINDER'S REFUSAL, and one claimed by another feature
 *     type is the tool's own BadInput, naming the output by its URI string;
 *   - -L SUBSETS THE RECORDS, a count overlapping two intervals being written once, and NEEDS
 *     --sequence-dictionary: the counts header's dictionary is not one the engine consults;
 *   - AND ANY OTHER FEATURE TYPE IS REFUSED BY FeatureWalker BEFORE THE TOOL RUNS, with
 *     "contains features of the wrong type", as is an input that does not exist.
 *
 * Output:
 *
 *     input\t<label>=<the whole input file, escaped>
 *     out\t<label>=<the whole text output, escaped>
 *     gz\t<label>=<the block-compressed output, base64>
 *     depth\t<label>=<the same, decompressed and escaped>
 *     tbi\t<label>=<its tabix index, base64>
 *     error\t<label>\t<exception class>:<message>
 *     cause\t<label>\t<exception class>:<message>, one per cause, outermost first
 *
 * Usage: ConvertCountsToDepthFileDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.sv.ConvertCountsToDepthFile;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.zip.GZIPInputStream;

public class ConvertCountsToDepthFileDump {

    static final String COLUMNS = "CONTIG\tSTART\tEND\tCOUNT\n";

    /** One count, written one-based and closed as a .counts.tsv holds it. */
    static String count(final String contig, final int start, final int end, final int count) {
        return contig + "\t" + start + "\t" + end + "\t" + count + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("convert-counts-to-depth-file-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConvertCountsToDepthFileDump: a counts file rewritten as depth evidence");

        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "counts", List.of("chr1", "chr2"));
        final String samHeader = Files.readString(dict);

        final String records = count("chr1", 1, 100, 31)
                + count("chr1", 101, 200, 32)
                + count("chr1", 201, 300, 0)
                + count("chr2", 1, 100, 33);
        // A read group whose own ID the output does not carry, and a @PG and @CO it drops too.
        final String counts = samHeader
                + "@RG\tID:not-the-cnv-id\tSM:gamma\n"
                + "@PG\tID:something\tPN:something\n"
                + "@CO\ta comment the output will not carry\n"
                + COLUMNS + records;

        run(dir, "plain", counts, "out.rd.txt");
        // The sample named on the command line wins over the header's, by either name.
        run(dir, "sample-name", counts, "out.rd.txt", "--sample-name", "delta");
        run(dir, "sample-name-short", counts, "out.rd.txt", "-n", "epsilon");
        // A level a text output never reads.
        run(dir, "plain-level-9", counts, "out.rd.txt", "--compression-level", "9");
        // A header and no records.
        run(dir, "empty", samHeader + "@RG\tID:one\tSM:gamma\n" + COLUMNS, "out.rd.txt");

        // Block compressed at the default level, at the two ends of the range, and empty.
        run(dir, "gz", counts, "out.rd.txt.gz");
        run(dir, "gz-level-1", counts, "out.rd.txt.gz", "--compression-level", "1");
        run(dir, "gz-level-9", counts, "out.rd.txt.gz", "--compression-level", "9");
        run(dir, "gz-empty", samHeader + "@RG\tID:one\tSM:gamma\n" + COLUMNS, "out.rd.txt.gz");
        // Out of the declared range.
        run(dir, "level-10", counts, "out.rd.txt", "--compression-level", "10");

        // Intervals: the second count overlaps two of them, and is written once. The counts
        // header's dictionary is not the engine's, so without --sequence-dictionary -L is refused.
        run(dir, "intervals-no-dictionary", counts, "out.rd.txt", "-L", "chr2");
        run(dir, "intervals", counts, "out.rd.txt", "--sequence-dictionary", dict.toString(),
                "-L", "chr1:150-160", "-L", "chr1:180-250", "-L", "chr2");

        // A read group with no SM, with and without --sample-name, and no read group at all.
        final String noSample = samHeader + "@RG\tID:one\n" + COLUMNS + records;
        run(dir, "no-sample", noSample, "out.rd.txt");
        run(dir, "no-sample-named", noSample, "out.rd.txt", "--sample-name", "delta");
        run(dir, "no-read-groups", samHeader + COLUMNS + records, "out.rd.txt",
                "--sample-name", "delta");

        // An output named for another feature type, and one no codec claims.
        run(dir, "wrong-output-type", counts, "out.baf.txt");
        run(dir, "no-output-codec", counts, "out.txt");

        // An input of another feature type, and one that does not exist.
        runInput(dir, "wrong-feature-type", "input.rd.txt", "#Chr\tStart\tEnd\tgamma\nchr1\t0\t100\t31\n",
                "out.rd.txt");
        final Path work = dir.resolve("missing-input");
        Files.createDirectories(work);
        execute(dir, "missing-input", work, work.resolve("absent.counts.tsv"), "out.rd.txt");
    }

    static void run(final Path dir, final String label, final String text, final String output,
                    final String... extra) throws Exception {
        runInput(dir, label, "input.counts.tsv", text, output, extra);
    }

    static void runInput(final Path dir, final String label, final String inputName,
                         final String text, final String output, final String... extra)
            throws Exception {
        final Path work = dir.resolve(label);
        Files.createDirectories(work);
        final Path input = work.resolve(inputName);
        Files.writeString(input, text, StandardCharsets.UTF_8);
        try {
            new IndexFeatureFile().instanceMain(new String[] {"-I", input.toString()});
        } catch (final Exception e) {
            System.out.printf("index\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
        }
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(masked(text, dir)));
        execute(dir, label, work, input, output, extra);
    }

    static void execute(final Path dir, final String label, final Path work, final Path input,
                        final String output, final String... extra) throws Exception {
        final Path out = work.resolve(output);
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--counts-filename", input.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new ConvertCountsToDepthFile().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                System.out.printf("cause\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            }
        }
        if (!Files.exists(out)) {
            return;
        }
        if (output.endsWith(".gz")) {
            System.out.printf("gz\t%s=%s%n", label, base64(out));
            try (InputStream in = new GZIPInputStream(Files.newInputStream(out))) {
                System.out.printf("depth\t%s=%s%n", label, ReferenceQueryDump.escape(
                        new String(in.readAllBytes(), StandardCharsets.UTF_8)));
            }
            final Path tbi = Path.of(out + ".tbi");
            if (Files.exists(tbi)) {
                System.out.printf("tbi\t%s=%s%n", label, base64(tbi));
            }
        } else {
            System.out.printf("out\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
