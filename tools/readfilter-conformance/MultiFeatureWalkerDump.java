/*
 * ExampleMultiFeatureWalker's output, taken from the reference.
 *
 * The tool is forty-eight lines and prints every feature it is handed; all of its substance is
 * `MultiFeatureWalker`, the engine class under `SiteDepthtoBAF`, `CondenseDepthEvidence`,
 * `PrintSVEvidence` and the rest of the SV family. So what is measured here is the merge.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE MERGE IS A PriorityQueue, WHICH IS A BINARY HEAP AND NOT STABLE, so two features that
 *     compare equal come out in heap order rather than in the order their files were named, and
 *     the order is decided by the sequence of insertions rather than by anything the user wrote;
 *   - THE INPUTS ARE ITERATED AS A Set, `FeatureManager.getAllInputs()`, so even the insertion
 *     order is not the command line's;
 *   - THE COMPARISON IS CONTIG INDEX, THEN START, THEN END, and the contig index comes from the
 *     dictionary rather than from the name, so the dictionary's order is the file order;
 *   - A CONTIG ABSENT FROM THE DICTIONARY HAS INDEX -1, WHICH SORTS BEFORE EVERY NAMED CONTIG, and
 *     the consequence is not a missing-contig message: a file holding chr1 then chr2, under a
 *     dictionary naming chr1 and chr3, is reported as `inputs are not sorted at chr2:101`. The
 *     file is sorted. The dictionary is incomplete. The diagnostic says neither;
 *   - AND THAT PATH IS ONLY REACHABLE WHEN THE INPUT ALSO HOLDS A CONTIG THE DICTIONARY NAMES:
 *     with no overlap at all the run is refused earlier, by the engine's dictionary comparison,
 *     as `IncompatibleSequenceDictionaries ... No overlapping contigs found`;
 *   - THE SORT CHECK IS PER INPUT AND FIRES LATE: `next()` compares the replacement drawn from the
 *     same input against the entry it just returned, so an unsorted file is caught only when its
 *     next record is pulled, and the message names the locus of the NEW feature;
 *   - AND THAT CHECK IS UNREACHABLE THROUGH AN INDEXED FILE, because an unsorted file cannot be
 *     indexed at all: `IndexFeatureFile` refuses first, with a `CouldNotIndexFile` wrapping
 *     tribble's own complaint, which names the two starts and reads the other way round;
 *   - THE MOST COMPREHENSIVE DICTIONARY WINS, the larger of any two, and the smaller must be a
 *     subset in the same relative order;
 *   - A CONTIG IN THE SMALLER DICTIONARY THAT IS ABSENT FROM THE LARGER IS A UserException naming
 *     both sources;
 *   - CONTIGS IN A DIFFERENT ORDER ARE A UserException naming the contig, the one it should have
 *     followed, and both sources;
 *   - AND NO DICTIONARY AT ALL IS REFUSED, because a `.rd.txt` header carries sample names and a
 *     null dictionary, so the file cannot supply one.
 *
 * Output:
 *
 *     input\t<label>\t<name>=<the whole .rd.txt, escaped>
 *     dict\t<label>\t<name>=<the whole .dict, escaped>
 *     features\t<label>=<everything the tool printed, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MultiFeatureWalkerDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.examples.ExampleMultiFeatureWalker;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MultiFeatureWalkerDump {

    static final String HEADER = "#Chr\tStart\tEnd\tsampleA\n";

    /** One depth bin, written zero-based and half-open as the file holds it. */
    static String bin(final String contig, final int start, final int end, final int count) {
        return contig + "\t" + start + "\t" + end + "\t" + count + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("multi-feature-walker-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MultiFeatureWalkerDump: several locus-sorted feature files merged");

        // Three dictionaries: the full one, one holding fewer contigs, and one holding the same
        // contigs in another order. The last two exist to drive `betterDictionary`.
        final Path full = writeDictionary(dir, "full", List.of("chr1", "chr2", "chr3"));
        final Path subset = writeDictionary(dir, "subset", List.of("chr1", "chr3"));
        final Path shuffled = writeDictionary(dir, "shuffled", List.of("chr2", "chr1", "chr3"));
        final Path other = writeDictionary(dir, "other", List.of("chr1", "chrX"));

        // Two files whose loci alternate, which is the merge doing what it says.
        run(dir, "interleaved", full,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1) + bin("chr1", 300, 400, 3)),
                        named("b", HEADER + bin("chr1", 200, 300, 2) + bin("chr1", 400, 500, 4))));

        // Two files whose records are the SAME interval, so the queue decides and nothing else can.
        run(dir, "tie-two", full,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1)),
                        named("b", HEADER + bin("chr1", 100, 200, 2))));

        // Three of them, where a heap's order is least like the order they were named in.
        run(dir, "tie-three", full,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1)),
                        named("b", HEADER + bin("chr1", 100, 200, 2)),
                        named("c", HEADER + bin("chr1", 100, 200, 3))));

        // The same start, different ends, which is the third key of the comparison.
        run(dir, "same-start", full,
                List.of(named("a", HEADER + bin("chr1", 100, 500, 1)),
                        named("b", HEADER + bin("chr1", 100, 200, 2))));

        // Two contigs, so the dictionary's order rather than the name's decides.
        run(dir, "two-contigs", full,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1) + bin("chr3", 100, 200, 3)),
                        named("b", HEADER + bin("chr2", 100, 200, 2))));

        // A contig the dictionary does not name, whose index is -1.
        run(dir, "unknown-contig", subset,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1) + bin("chr3", 100, 200, 3)),
                        named("b", HEADER + bin("chr2", 100, 200, 2))));

        // A contig the dictionary does not name, in a file that also names one it does, which is
        // the only way past the compatibility check and into the index of -1.
        run(dir, "partly-unknown-contig", subset,
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1) + bin("chr2", 100, 200, 2))));

        // One file that goes backwards. It cannot be indexed at all, so the walker's own check is
        // reachable only from an unindexed input, and both refusals are worth having.
        run(dir, "unsorted-index", full,
                List.of(named("a", HEADER + bin("chr1", 300, 400, 3) + bin("chr1", 100, 200, 1))));
        run(dir, "unsorted", full,
                List.of(unindexed("a", HEADER + bin("chr1", 300, 400, 3) + bin("chr1", 100, 200, 1)),
                        unindexed("b", HEADER + bin("chr1", 200, 300, 2))));

        // No dictionary at all: the .rd.txt header carries sample names and a null dictionary.
        runRaw(dir, "no-dictionary",
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1))), List.of());

        // A reference whose dictionary disagrees with the master one, both ways round.
        runRaw(dir, "dict-contig-absent",
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1))),
                List.of("--sequence-dictionary", other.toString(),
                        "--reference", dir.resolve("full.fasta").toString()));
        runRaw(dir, "dict-out-of-order",
                List.of(named("a", HEADER + bin("chr1", 100, 200, 1))),
                List.of("--sequence-dictionary", shuffled.toString(),
                        "--reference", dir.resolve("full.fasta").toString()));
    }

    record Input(String name, String text, boolean indexed) {}

    static Input named(final String name, final String text) {
        return new Input(name, text, true);
    }

    /** The same, left unindexed, which is the only way the walker's own sort check is reachable. */
    static Input unindexed(final String name, final String text) {
        return new Input(name, text, false);
    }

    /** A fasta of the named contigs, and the .dict CreateSequenceDictionary writes beside it. */
    static Path writeDictionary(final Path dir, final String label, final List<String> contigs)
            throws Exception {
        final Path fasta = dir.resolve(label + ".fasta");
        final StringBuilder bases = new StringBuilder();
        for (final String contig : contigs) {
            bases.append(">").append(contig).append("\n");
            for (int i = 0; i < 20; i++) {
                bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
            }
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve(label + ".dict");
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dict});
        System.out.printf("dict\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(dict), dir)));
        return dict;
    }

    static void run(final Path dir, final String label, final Path dictionary,
                    final List<Input> inputs) throws Exception {
        runRaw(dir, label, inputs, List.of("--sequence-dictionary", dictionary.toString()));
    }

    static void runRaw(final Path dir, final String label, final List<Input> inputs,
                       final List<String> extra) throws Exception {
        final List<String> argv = new ArrayList<>();
        for (final Input input : inputs) {
            final Path path = dir.resolve(label + "-" + input.name() + ".rd.txt");
            Files.writeString(path, input.text(), StandardCharsets.UTF_8);
            System.out.printf("input\t%s\t%s=%s%n", label, input.name(),
                    ReferenceQueryDump.escape(input.text()));
            if (input.indexed()) {
                try {
                    new IndexFeatureFile().instanceMain(new String[] {"-I", path.toString()});
                } catch (final Exception | AssertionError e) {
                    // An unsorted file cannot be indexed at all, which is a refusal of its own and
                    // the reason the walker's check needs an unindexed input to be reached.
                    System.out.printf("index-error\t%s\t%s\t%s:%s%n", label, input.name(),
                            e.getClass().getName(),
                            ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
                    return;
                }
            }
            argv.addAll(Arrays.asList("-F", path.toString()));
        }
        argv.addAll(extra);

        // The tool prints each feature to System.out, which is where this dump writes too, so the
        // stream is swapped for the duration of the run and put back before anything is reported.
        final PrintStream saved = System.out;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        String error = null;
        try {
            System.setOut(new PrintStream(captured, true, StandardCharsets.UTF_8));
            try {
                new ExampleMultiFeatureWalker().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception | AssertionError e) {
                error = e.getClass().getName() + ":" + masked(String.valueOf(e.getMessage()), dir);
            }
        } finally {
            System.setOut(saved);
        }
        if (error != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(error));
            return;
        }
        System.out.printf("features\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(captured.toString(StandardCharsets.UTF_8), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
