/*
 * TagGermlineEvents, taken from the reference.
 *
 * A tumour segment file and a matched normal one, and the question asked of each tumour segment:
 * does the normal carry the same event? The answer is written as one more column.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE NORMAL IS MERGED BY CALL FIRST, run by run, and the merged region keeps ONLY the call
 *     annotation: everything else on those rows is dropped before the comparison;
 *   - ONLY NON-NEUTRAL MERGED NORMAL REGIONS ARE CONSIDERED, so a normal that is neutral
 *     everywhere tags nothing however well it lines up;
 *   - A NORMAL REGION TAGS ITS TUMOUR SEGMENTS WHEN BOTH BREAKPOINTS ARE SEEN WITHIN THE PADDING,
 *     OR when a merged tumour run reciprocally overlaps it past the threshold. The two are an OR,
 *     so a shifted pair of breakpoints can still tag through the overlap;
 *   - THE BREAKPOINT SEARCH IS OVER THE UNMERGED TUMOUR SEGMENTS while the reciprocal overlap is
 *     over the MERGED ones, so the same pair of files can answer differently for the two;
 *   - THE PER-SEGMENT FILTER IS A DIFFERENT TEST AGAIN: a tumour segment is tagged when one of its
 *     own breakpoints is within the padding OR the normal region strictly contains it, AND the
 *     intersection is STRICTLY MORE than its own length times the threshold;
 *   - THE INTERSECTION TEST IS STRICT WHILE THE RECIPROCAL ONE IS NOT, `>` against `>=`;
 *   - THE DEFAULT TAG IS `0`, which is what `CalledCopyRatioSegment.Call.NEUTRAL` renders as, so
 *     an untagged segment carries a value rather than an empty column;
 *   - AN EMPTY CALL ANYWHERE IS A REFUSAL, on either side, and so is a padding below zero or a
 *     threshold outside the unit interval;
 *   - AND THE OUTPUT IS THE TUMOUR FILE'S OWN ANNOTATIONS PLUS THE TAG, sorted alphabetically, so
 *     the tag column lands where its name sorts rather than at the end.
 *
 * Output:
 *
 *     tagged\t<label>=<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: TagGermlineEventsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.utils.TagGermlineEvents;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TagGermlineEventsDump {

    /** A tumour file of four segments with two annotations besides the call. */
    static final String TUMOUR =
            "CONTIG\tSTART\tEND\tCALL\tMEAN_LOG2_COPY_RATIO\tNUM_POINTS\n"
            + "chr1\t1\t100\t+\t0.7\t10\n"
            + "chr1\t101\t200\t+\t0.7\t10\n"
            + "chr1\t201\t300\t0\t0.0\t10\n"
            + "chr2\t1\t100\t-\t-0.8\t10\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("tag-germline-events-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# TagGermlineEventsDump: which tumour segments the normal explains");

        // A normal whose amplified run covers both tumour segments exactly.
        run("exact-match", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t1\t100\t+\n"
                + "chr1\t101\t200\t+\n"
                + "chr1\t201\t300\t0\n"
                + "chr2\t1\t100\t0\n");
        // The same, with the normal's breakpoints shifted past the default padding. It still
        // tags, through the reciprocal overlap rather than through the breakpoints, which is the
        // OR in the open.
        run("shifted", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t20\t180\t+\n"
                + "chr1\t201\t300\t0\n"
                + "chr2\t1\t100\t0\n");
        // The same shift with a padding wide enough to see the breakpoints.
        runWith("shifted-padded", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t20\t180\t+\n"
                + "chr1\t201\t300\t0\n"
                + "chr2\t1\t100\t0\n",
                "--endpoint-padding", "25");
        // A normal that is neutral everywhere, which tags nothing.
        run("all-neutral", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t1\t300\t0\n"
                + "chr2\t1\t100\t0\n");
        // A normal deletion matching the tumour's deletion on the second contig.
        run("second-contig", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t1\t300\t0\n"
                + "chr2\t1\t100\t-\n");
        // A normal region strictly containing one tumour segment, which reaches the third branch
        // of the per-segment filter.
        run("normal-contains", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t1\t200\t+\n"
                + "chr1\t201\t300\t0\n"
                + "chr2\t1\t100\t0\n");
        // A reciprocal overlap without matching breakpoints, at a threshold low enough to see it.
        runWith("reciprocal", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t40\t240\t+\n"
                + "chr1\t241\t300\t0\n"
                + "chr2\t1\t100\t0\n",
                "--reciprocal-threshold", "0.5");
        // The same pair at the default threshold.
        run("reciprocal-default", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t40\t240\t+\n"
                + "chr1\t241\t300\t0\n"
                + "chr2\t1\t100\t0\n");
        // A threshold of zero, which makes every overlap reciprocal.
        runWith("zero-threshold", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\n"
                + "chr1\t40\t240\t+\n"
                + "chr1\t241\t300\t0\n"
                + "chr2\t1\t100\t0\n",
                "--reciprocal-threshold", "0.0");
        // A different call column name.
        runWith("other-call-column", dir, fasta,
                "CONTIG\tSTART\tEND\tcall_state\n"
                + "chr1\t1\t100\t+\n"
                + "chr1\t101\t200\t0\n",
                "CONTIG\tSTART\tEND\tcall_state\n"
                + "chr1\t1\t100\t+\n"
                + "chr1\t101\t200\t0\n",
                "--input-call-header", "call_state");
        // The refusals: an empty call on each side, a negative padding, a threshold above one.
        run("empty-tumour-call", dir, fasta,
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t\n",
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n");
        run("empty-normal-call", dir, fasta,
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n",
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t\n");
        runWith("negative-padding", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n", "--endpoint-padding", "-1");
        runWith("threshold-above-one", dir, fasta, TUMOUR,
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n", "--reciprocal-threshold", "1.5");
        // Overlapping tumour segments, which the tagger refuses outright.
        run("overlapping-tumour", dir, fasta,
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\nchr1\t50\t150\t+\n",
                "CONTIG\tSTART\tEND\tCALL\nchr1\t1\t100\t+\n");
    }

    static void run(final String label, final Path dir, final Path fasta, final String tumour,
                    final String normal) throws Exception {
        runWith(label, dir, fasta, tumour, normal);
    }

    static void runWith(final String label, final Path dir, final Path fasta, final String tumour,
                        final String normal, final String... extra) throws Exception {
        final Path tumourFile = dir.resolve(label + "-tumour.seg");
        final Path normalFile = dir.resolve(label + "-normal.seg");
        Files.write(tumourFile, tumour.getBytes());
        Files.write(normalFile, normal.getBytes());
        final Path out = dir.resolve("tagged-" + label + ".seg");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "--segments", tumourFile.toString(),
                "--called-matched-normal-seg-file", normalFile.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new TagGermlineEvents().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("tagged\t%s=%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
