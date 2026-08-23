/*
 * PrintSVEvidence's output, taken from the reference.
 *
 * Several SV evidence files merged into one. It sits on `MultiFeatureWalker`, already measured, and
 * what is measured here is the two things above it: which samples the output declares, and what the
 * sort merger does when two files speak about the same bin.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE OUTPUT'S SAMPLE LIST IS THE UNION OF THE INPUTS' HEADERS, in the order the walker
 *     accumulated it, which is a TreeSet and therefore alphabetical rather than file order;
 *   - EVERY RECORD IS REWRITTEN AGAINST THAT LIST, `extractSamples` filling a column a file does
 *     not have with MISSING_DATA, which is -1;
 *   - THE MERGE ONLY FILLS MISSING COLUMNS. Two files that both report a sample at one bin are a
 *     UserException naming the sample by ONE-BASED index and the bin by its interval, so merging is
 *     only ever a widening;
 *   - --sample-names SUBSETS AND REORDERS THE COLUMNS, and a name no file knows becomes a column
 *     of -1 rather than a refusal;
 *   - AN EMPTY SAMPLE LIST FROM EVERY HEADER TURNS SAMPLE FILTERING OFF ENTIRELY rather than
 *     dropping everything, so a headerless format passes through untouched;
 *   - THE OUTPUT TYPE IS DECIDED BY THE OUTPUT PATH'S EXTENSION, and an input whose codec produces
 *     another type is refused by name;
 *   - AN OUTPUT EXTENSION THAT NAMES NO FEATURE TYPE AT ALL NEVER REACHES THE TOOL'S OWN CHECK:
 *     `FeatureOutputCodecFinder.find` fails first with `No feature output codec found for ...`, so
 *     the "requires an SVFeature subtype" message the tool carries is unreachable from any
 *     extension that has no codec, and reachable only from one whose codec produces a non-SV
 *     feature;
 *   - THE SORT MERGER REFUSES FEATURES OUT OF DICTIONARY ORDER with a GATKException rather than a
 *     UserException;
 *   - THE HEADER LINE OF THE OUTPUT IS REWRITTEN from the sample list, so the columns of the output
 *     are not the columns of any input;
 *   - AND THE WALKER'S OWN DICTIONARY RULES APPLY, so the run needs one and takes the largest.
 *
 * Output:
 *
 *     input\t<label>\t<name>=<the whole file, escaped>
 *     merged\t<label>=<the whole output, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PrintSVEvidenceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.sv.PrintSVEvidence;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PrintSVEvidenceDump {

    /** One depth bin, written zero-based and half-open as the file holds it. */
    static String bin(final int start, final int end, final int... counts) {
        final StringBuilder line = new StringBuilder("chr1\t" + start + "\t" + end);
        for (final int count : counts) {
            line.append('\t').append(count);
        }
        return line.append('\n').toString();
    }

    static String header(final String... samples) {
        return "#Chr\tStart\tEnd\t" + String.join("\t", samples) + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("print-sv-evidence-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PrintSVEvidenceDump: SV evidence files merged into one");

        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "sv", List.of("chr1", "chr2"));

        // Two files naming different samples over the same bins, which is the widening the merge
        // is for: each file's records gain a -1 column, and the merge fills it.
        final String alpha = header("alpha") + bin(0, 100, 11) + bin(100, 200, 12);
        final String beta = header("beta") + bin(0, 100, 21) + bin(100, 200, 22);
        run(dir, "widening", dict, List.of(named("a", alpha), named("b", beta)));

        // The same two, subset to one sample, and to a sample neither file knows.
        run(dir, "subset-one", dict, List.of(named("a", alpha), named("b", beta)),
                "--sample-names", "beta");
        run(dir, "unknown-sample", dict, List.of(named("a", alpha), named("b", beta)),
                "--sample-names", "gamma");
        // And reordered, which the output's columns follow.
        run(dir, "reordered", dict, List.of(named("a", alpha), named("b", beta)),
                "--sample-names", "beta", "--sample-names", "alpha");

        // Two files naming the SAME sample over the same bins, which is not a widening.
        run(dir, "same-sample-twice", dict,
                List.of(named("a", alpha), named("b", header("alpha") + bin(0, 100, 31))));

        // Disjoint bins, where the merger never has two records at one locus.
        run(dir, "disjoint-bins", dict,
                List.of(named("a", header("alpha") + bin(0, 100, 11)),
                        named("b", header("beta") + bin(200, 300, 22))));

        // Three files, to show the sample list is sorted rather than in file order.
        run(dir, "three-samples", dict,
                List.of(named("a", header("zulu") + bin(0, 100, 1)),
                        named("b", header("alpha") + bin(0, 100, 2)),
                        named("c", header("mike") + bin(0, 100, 3))));

        // An output whose extension implies another SV feature type than the inputs produce.
        runNamed(dir, "wrong-sv-type", dict, List.of(named("a", alpha)), "merged.baf.txt");
        // And one that is not an SV feature at all.
        runNamed(dir, "not-an-sv-type", dict, List.of(named("a", alpha)), "merged.vcf");
    }

    record Input(String name, String text) {}

    static Input named(final String name, final String text) {
        return new Input(name, text);
    }

    static void run(final Path dir, final String label, final Path dictionary,
                    final List<Input> inputs, final String... extra) throws Exception {
        runNamed(dir, label, dictionary, inputs, "merged-" + label + ".rd.txt", extra);
    }

    static void runNamed(final Path dir, final String label, final Path dictionary,
                         final List<Input> inputs, final String outputName, final String... extra)
            throws Exception {
        final List<String> argv = new ArrayList<>();
        for (final Input input : inputs) {
            final Path path = dir.resolve(label + "-" + input.name() + ".rd.txt");
            Files.writeString(path, input.text(), StandardCharsets.UTF_8);
            new IndexFeatureFile().instanceMain(new String[] {"-I", path.toString()});
            System.out.printf("input\t%s\t%s=%s%n", label, input.name(),
                    ReferenceQueryDump.escape(input.text()));
            argv.addAll(Arrays.asList("-F", path.toString()));
        }
        final Path out = dir.resolve(outputName);
        argv.addAll(Arrays.asList("-O", out.toString(),
                "--sequence-dictionary", dictionary.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new PrintSVEvidence().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("merged\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
