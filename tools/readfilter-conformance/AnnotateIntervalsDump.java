/*
 * AnnotateIntervals, taken from the reference.
 *
 * GC content per interval, and two optional annotations read from BED tracks. Everything here is
 * arithmetic over the reference and over feature overlaps, and two of the three have an off-by-one
 * or a missing-value rule that a port would otherwise invent.
 *
 * Six behaviours this is built to catch.
 *
 *   - GC CONTENT IS OVER ACGT ALONE, and the denominator is `gc + at` rather than the interval's
 *     length. The reader has already turned every IUPAC code into `N`, so an interval holding
 *     ambiguity codes has a SMALLER denominator than its length, and one holding nothing else is a
 *     NaN rather than a zero;
 *   - THE BED TRACK IS ZERO-BASED AND THE INTERVAL IS ONE-BASED, and the reference bridges that by
 *     passing `feature.getEnd() - 1` into an overlap computed against the interval's one-based
 *     coordinates. So a BED line covering exactly one interval contributes its length minus one;
 *   - A MISSING OR NaN SCORE IS ONE, not zero and not a skip, so a BED file with no score column
 *     annotates every overlap at full weight;
 *   - THE ANNOTATION IS A LENGTH-WEIGHTED AVERAGE over the INTERVAL's length, so a track covering
 *     half an interval at score 1 gives about a half;
 *   - AN OVERLAPPING TRACK IS REFUSED. The tool merges the track's features with
 *     `OVERLAPPING_ONLY` and compares the counts, so two features that touch are fine and two that
 *     overlap are a `UserException.BadInput`;
 *   - THE OUTPUT'S COLUMNS FOLLOW THE ANNOTATORS THAT RAN, so a run without a track has three
 *     columns and one with both has five, in the order the annotators were added;
 *   - AND THE ENGINE'S DEFAULT MERGING RULE IS REFUSED. This tool requires
 *     `--interval-merging-rule OVERLAPPING_ONLY` and throws before reading anything without it,
 *     which is a validation no other tool in this port makes.
 *
 * Output:
 *
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AnnotateIntervalsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AnnotateIntervalsDump {

    /** A mappability track: scored, unscored, and one that only partly covers an interval. */
    static final String MAPPABILITY =
            "chr1\t0\t12\tfull\t1.0\n"
            + "chr1\t12\t18\thalf\t0.5\n"
            + "chr1\t24\t30\tunscored\n";

    /** Two features that overlap, which the tool refuses. */
    static final String OVERLAPPING =
            "chr1\t0\t12\tone\t1.0\n"
            + "chr1\t6\t18\ttwo\t1.0\n";

    /** Two features that touch but do not overlap, which it accepts. */
    static final String TOUCHING =
            "chr1\t0\t12\tone\t1.0\n"
            + "chr1\t12\t18\ttwo\t1.0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("annotateintervals-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path mappability = write(dir, "mappability.bed", MAPPABILITY);
        final Path overlapping = write(dir, "overlapping.bed", OVERLAPPING);
        final Path touching = write(dir, "touching.bed", TOUCHING);

        System.out.println("# AnnotateIntervalsDump: GC content and BED-weighted annotations");

        // The first twelve bases are plain `ACGT`; the next twelve are soft-masked; then the IUPAC
        // stretch, which the reader has already turned into Ns.
        run("gc-only", dir, fasta, "-L", "chr1:1-12", "-L", "chr1:13-24", "-L", "chr1:25-36");
        // An interval inside the N run, whose GC content is a NaN.
        run("all-n", dir, fasta, "-L", "chr1:29-32");
        // One base, which is either 0 or 1 and never anything between.
        run("single-base", dir, fasta, "-L", "chr1:1-1", "-L", "chr1:2-2");
        // The mappability track, whose three features cover the first interval fully, the second
        // partly, and the third with no score at all.
        run("mappability", dir, fasta, "-L", "chr1:1-12", "-L", "chr1:13-24", "-L", "chr1:25-30",
                "--mappability-track", mappability.toString());
        // Both tracks at once, which decides the column order.
        run("both-tracks", dir, fasta, "-L", "chr1:1-12",
                "--mappability-track", mappability.toString(),
                "--segmental-duplication-track", touching.toString());
        // A track whose features overlap.
        run("overlapping-track", dir, fasta, "-L", "chr1:1-12",
                "--mappability-track", overlapping.toString());
        // A track whose features touch, which is allowed.
        run("touching-track", dir, fasta, "-L", "chr1:1-12",
                "--mappability-track", touching.toString());
        // And the engine's default merging rule, which this tool will not run under.
        run("default-merging", dir, fasta, "-L", "chr1:1-12");
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        // A track is QUERIED by interval, so it needs an index -- except on the overlap check,
        // which iterates the whole file and therefore fires before the index is ever wanted.
        if (!name.startsWith("overlapping")) {
            new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                    .instanceMain(new String[] {"-I", path.toString()});
        }
        return path;
    }

    static void run(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final Path out = dir.resolve("annotated-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-O", out.toString()));
        // Every case but `default-merging` passes it: this tool REFUSES the engine's default
        // merging rule, which is `ALL`, and says so before it reads anything.
        if (!label.equals("default-merging")) {
            argv.add("--interval-merging-rule");
            argv.add("OVERLAPPING_ONLY");
        }
        argv.addAll(Arrays.asList(extra));
        try {
            new AnnotateIntervals().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
