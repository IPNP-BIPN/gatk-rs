/*
 * PreprocessIntervals, taken from the reference.
 *
 * The bins every copy-number tool counts over: the input intervals padded, de-overlapped, split
 * into bins and stripped of the bins that are all N.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE OVERLAP IS RESOLVED AT THE MIDPOINT OF THE ORIGINAL INTERVALS, not of the padded ones:
 *     `(originalThisEnd + originalNextStart) / 2`, an integer division that rounds toward zero and
 *     therefore leans left;
 *   - THE PASS IS SEQUENTIAL AND IN PLACE. Interval i+1 is rewritten and then compared against i+2
 *     in the next step, so three intervals padded into one another resolve left to right and the
 *     original of i+1 is still what the midpoint of the second pair is taken from;
 *   - THE MIDPOINT CANNOT FALL OUTSIDE EITHER INTERVAL, however large the padding is, because the
 *     inputs arrive sorted and merged: the midpoint of two originals is at least the first's end
 *     and at most the second's start, so the `SimpleInterval` the resolution builds never throws.
 *     The `heavy-padding` row is two intervals two bases apart padded by sixty, which is as close
 *     as that can be driven;
 *   - PADDING CLAMPS AT ONE AND AT THE CONTIG LENGTH, both from the sequence dictionary rather
 *     than from the reference itself;
 *   - THE BINS ARE LAID FROM THE INTERVAL'S START, so the last bin of an interval is short rather
 *     than the first, and an interval shorter than one bin is one bin of its own length;
 *   - A BIN LENGTH OF ZERO MEANS NO BINNING AT ALL, and the padded intervals come out as they are;
 *   - THE N FILTER IS `allMatch`, so a bin of one non-N base survives and an EMPTY bin cannot
 *     occur; the decode is case-insensitive, so a bin of lower-case n is dropped too;
 *   - AND NO INTERVALS AT ALL MEANS THE WHOLE REFERENCE, one interval per contig, which is not the
 *     same as passing every contig by name: the padding then has nothing to clamp.
 *
 * Output:
 *
 *     list\t<label>\t<the whole interval list, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PreprocessIntervalsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.PreprocessIntervals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PreprocessIntervalsDump {

    /**
     * Two contigs of 240 bases, four lines of sixty each. chr1 is upper-case ACGT to 60,
     * lower-case to 120, a run of sixty Ns to 180 and upper-case again to 240; chr2 is all N but
     * for an `AC` at 150 and 151, with its second line written in lower case.
     */
    public static final String FASTA = ">chr1\n"
            + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            + "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n"
            + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            + ">chr2\n"
            + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            + "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn\n"
            + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNACNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"
            + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("preprocess-intervals-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# PreprocessIntervalsDump: the bins a copy-number run counts over");

        // No intervals at all: the whole reference, one interval per contig.
        run("whole-genome", dir, fasta, "--bin-length", "50", "--padding", "0");
        // The same, with the padding that has nothing to clamp.
        run("whole-genome-padded", dir, fasta, "--bin-length", "50");
        // A bin length that does not divide the interval, so the last bin is short.
        run("uneven-bins", dir, fasta, "--bin-length", "30", "--padding", "0", "-L", "chr1:1-100");
        // An interval shorter than one bin.
        run("short-interval", dir, fasta, "--bin-length", "50", "--padding", "0", "-L", "chr1:10-20");
        // No binning at all.
        run("no-bins", dir, fasta, "--bin-length", "0", "--padding", "0", "-L", "chr1:10-20",
                "-L", "chr1:150-170");
        // Padding clamped at both ends of a contig.
        run("clamped", dir, fasta, "--bin-length", "0", "--padding", "20", "-L", "chr1:5-10",
                "-L", "chr1:230-235");
        // Two intervals whose padding overlaps: the midpoint of the ORIGINALS decides, and the
        // division rounds down.
        run("overlap", dir, fasta, "--bin-length", "0", "--padding", "20", "-L", "chr1:10-20",
                "-L", "chr1:41-50");
        // The same pair one base apart, so the midpoint lands on an odd sum and rounds left.
        run("overlap-odd", dir, fasta, "--bin-length", "0", "--padding", "20", "-L", "chr1:10-20",
                "-L", "chr1:42-50");
        // Three intervals padded into one another, resolved left to right.
        run("overlap-three", dir, fasta, "--bin-length", "0", "--padding", "20", "-L", "chr1:10-20",
                "-L", "chr1:41-50", "-L", "chr1:71-80");
        // The same three, binned, so the bins of a de-overlapped interval are visible.
        run("overlap-three-binned", dir, fasta, "--bin-length", "10", "--padding", "20",
                "-L", "chr1:10-20", "-L", "chr1:41-50", "-L", "chr1:71-80");
        // The N run: bins entirely inside it are dropped, the bins around it are kept.
        run("n-run", dir, fasta, "--bin-length", "20", "--padding", "0", "-L", "chr1:101-200");
        // A contig that is all N but for two bases: only the bins holding them survive.
        run("almost-all-n", dir, fasta, "--bin-length", "10", "--padding", "0", "-L", "chr2");
        // The lower-case n stretch, which decodes to N and is dropped like any other, leaving a
        // list of no intervals at all.
        run("lower-case-n", dir, fasta, "--bin-length", "20", "--padding", "0", "-L", "chr2:61-120");
        // A bin straddling the edge of the N run, kept for the one base outside it.
        run("straddling-bin", dir, fasta, "--bin-length", "0", "--padding", "0",
                "-L", "chr1:120-140");
        // Two intervals two bases apart padded by sixty, which is as far as the resolution can be
        // pushed: the midpoint still lands between them.
        run("heavy-padding", dir, fasta, "--bin-length", "0", "--padding", "60",
                "-L", "chr1:100-101", "-L", "chr1:102-103");
        // Intervals on two contigs, which do not overlap whatever the padding.
        run("two-contigs", dir, fasta, "--bin-length", "0", "--padding", "20",
                "-L", "chr1:230-235", "-L", "chr2:1-5");
        // The three argument refusals: a merging rule other than OVERLAPPING_ONLY, the common
        // interval padding, and the common interval exclusion padding.
        runRaw("merging-rule", dir, fasta, "--bin-length", "0", "-L", "chr1:10-20",
                "--interval-merging-rule", "ALL");
        runRaw("interval-padding", dir, fasta, "--bin-length", "0", "-L", "chr1:10-20",
                "--interval-merging-rule", "OVERLAPPING_ONLY", "--interval-padding", "10");
        runRaw("interval-exclusion-padding", dir, fasta, "--bin-length", "0", "-L", "chr1:10-20",
                "--interval-merging-rule", "OVERLAPPING_ONLY", "--interval-exclusion-padding", "10");
        // A negative bin length and a negative padding, which the argument parser bounds.
        runRaw("negative-bin-length", dir, fasta, "--bin-length", "-1", "-L", "chr1:10-20",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        runRaw("negative-padding", dir, fasta, "--bin-length", "0", "--padding", "-1",
                "-L", "chr1:10-20", "--interval-merging-rule", "OVERLAPPING_ONLY");
    }

    /**
     * The sequence lines with their checksum and their URI blanked.
     *
     * Both come from the `.dict` the reference indexer wrote and neither is this tool's: the M5 is
     * a digest of the contig and the UR is wherever the fasta happened to sit. Masking them is what
     * keeps the golden about the bins.
     */
    public static String masked(final String text) {
        return text.replaceAll("\tM5:[0-9a-f]+", "\tM5:<masked>")
                .replaceAll("\tUR:file:[^\t\n]+", "\tUR:<masked>");
    }

    /**
     * A run with the merging rule the tool insists on, which every run that is not about the
     * argument checks has to pass: the common default is ALL and the tool refuses it.
     */
    static void run(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList(extra));
        if (Arrays.asList(extra).contains("-L")) {
            argv.addAll(Arrays.asList("--interval-merging-rule", "OVERLAPPING_ONLY"));
        }
        runRaw(label, dir, fasta, argv.toArray(new String[0]));
    }

    static void runRaw(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final Path out = dir.resolve("bins-" + label + ".interval_list");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new PreprocessIntervals().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("list\t%s\t%s%n", label,
                ReferenceQueryDump.escape(masked(new String(Files.readAllBytes(out)))));
    }
}
