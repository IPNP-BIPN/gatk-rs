/*
 * MergeAnnotatedRegionsByAnnotation, taken from the reference.
 *
 * The sibling of MergeAnnotatedRegions: neighbouring regions merge when the named annotations
 * agree and the gap between them is within a distance, rather than when they overlap.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE DISTANCE OF ABUTTING REGIONS IS ONE, NOT ZERO, whatever the method's own comment says:
 *     it answers zero for overlapping regions only, `IntervalUtils.overlaps` not counting an
 *     abuttal, and otherwise the gap between the nearer endpoints. So `1-100` and `101-200` are a
 *     distance of one apart and a maximum of zero merges nothing at all;
 *   - THE COMPARISON IS AGAINST THE MERGED REGION, so a chain walks: each merge widens the region
 *     that the next candidate is measured against, and a run of short hops merges however far it
 *     travels in total;
 *   - THE ANNOTATIONS ARE STILL RECONCILED THE SAME WAY, so the columns NOT named on the command
 *     line are split, deduplicated, sorted and rejoined with the separator while the named ones
 *     match by construction;
 *   - AN ANNOTATION NAMED THAT THE FILE DOES NOT CARRY IS A REFUSAL, and the message names the
 *     missing ones;
 *   - THE OUTPUT HEADER IS BUILT FROM THE FIRST MERGED REGION'S ANNOTATIONS, so a file of no
 *     regions dies rather than writing an empty result;
 *   - THE THREE LOCATABLE COLUMN NAMES ARE ARGUMENTS, so the output can carry any spelling and the
 *     `@CO` lines follow it;
 *   - AND NOTHING HERE WRITES A SAM HEADER OF THE INPUT'S: this tool goes through the writer
 *     directly rather than through the collection, so the sequence lines of the input are gone.
 *
 * Output:
 *
 *     merged\t<label>=<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MergeAnnotatedRegionsByAnnotationDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.utils.MergeAnnotatedRegionsByAnnotation;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MergeAnnotatedRegionsByAnnotationDump {

    /** Six regions with two annotations, two of them far apart and two abutting. */
    static final String SEGMENTS =
            "CONTIG\tSTART\tEND\tCALL\tNAME\n"
            + "chr1\t1\t100\t+\ta\n"
            + "chr1\t101\t200\t+\tb\n"
            + "chr1\t1201\t1300\t+\tc\n"
            + "chr1\t1301\t1400\t-\td\n"
            + "chr2\t1\t100\t+\te\n"
            + "chr2\t50000\t50100\t+\tf\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("merge-by-annotation-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# MergeAnnotatedRegionsByAnnotationDump: the segments an annotation merge leaves behind");

        // The default distance, which is a million bases: everything of the same call merges.
        run("default-distance", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL");
        // A distance of zero, which merges nothing: even the abutting pair is one apart.
        run("zero-distance", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL",
                "--max-merge-distance", "0");
        // A distance that reaches the thousand-base gap but not the fifty-thousand one.
        run("middle-distance", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL",
                "--max-merge-distance", "1001");
        // Matching on the other annotation, which never repeats, so nothing merges.
        run("match-on-name", dir, fasta, SEGMENTS, "--annotations-to-match", "NAME");
        // Matching on both, which is the same answer by another route.
        run("match-on-both", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL",
                "--annotations-to-match", "NAME");
        // The three output column names.
        run("renamed-columns", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL",
                "--output-contig-column", "chrom", "--output-start-column", "chromStart",
                "--output-end-column", "chromEnd");
        // An input carrying a SAM header, which this tool does not carry through.
        run("sam-header", dir, fasta,
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:240\nCONTIG\tSTART\tEND\tCALL\tNAME\n"
                + "chr1\t1\t100\t+\ta\n",
                "--annotations-to-match", "CALL");
        // An annotation the file does not carry.
        run("missing-annotation", dir, fasta, SEGMENTS, "--annotations-to-match", "ABSENT");
        // A negative distance, which the merge refuses.
        run("negative-distance", dir, fasta, SEGMENTS, "--annotations-to-match", "CALL",
                "--max-merge-distance", "-1");
    }

    static void run(final String label, final Path dir, final Path fasta, final String input,
                    final String... extra) throws Exception {
        final Path in = dir.resolve(label + ".seg");
        Files.write(in, input.getBytes());
        final Path out = dir.resolve("merged-" + label + ".seg");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "--segments", in.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new MergeAnnotatedRegionsByAnnotation().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("merged\t%s=%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
