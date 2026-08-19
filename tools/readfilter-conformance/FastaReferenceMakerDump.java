/*
 * FastaReferenceMaker, taken from the reference.
 *
 * The fourth member of the reference archetype, and the first that WRITES a reference rather than
 * reading one. It is a ReferenceWalker whose apply appends one base, and whose output is a FASTA
 * with its .fai and .dict beside it.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE OUTPUT SEQUENCES ARE NUMBERED, not named after the contig they came from: the name is a
 *     counter starting at one, and the contig and coordinates go in the DESCRIPTION;
 *   - A GAP STARTS A NEW SEQUENCE, and the test is withinDistanceOf(interval, 1), so two abutting
 *     intervals are one output sequence and two intervals a base apart are two;
 *   - CROSSING A CONTIG IS A GAP TOO, so a run with no -L over a two-contig reference writes two
 *     sequences and not one;
 *   - THE DESCRIPTION'S START IS THE SEQUENCE'S FIRST BASE and its end is the last position seen,
 *     so it is the span actually written rather than the interval asked for;
 *   - THE BASES ARE THE CACHING READER'S, so lower case comes back upper-cased and IUPAC codes come
 *     back as N: a FASTA written from this tool is not the FASTA it read;
 *   - AND --line-width IS PER RUN, recorded in every index row.
 *
 * Output:
 *
 *     fasta\t<label>\t<the FASTA text, escaped>
 *     fai\t<label>\t<the .fai text, escaped>
 *     dict\t<label>\t<the .dict text, escaped>
 *     error\t<label>\t<exception class>
 *
 * Usage: FastaReferenceMakerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.fasta.FastaReferenceMaker;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FastaReferenceMakerDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("fastareferencemaker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# FastaReferenceMakerDump: a reference written out of a reference");

        // No -L: the whole reference, which is two contigs and therefore two output sequences.
        run("all", dir, fasta);
        run("one-interval", dir, fasta, "-L", "chr1:1-12");
        // Two intervals a base apart, which is a gap and two sequences.
        run("gap", dir, fasta, "-L", "chr1:1-5", "-L", "chr1:7-12");
        // Two abutting intervals, which the merging rule joins before the walker sees them.
        run("abutting", dir, fasta, "-L", "chr1:1-5", "-L", "chr1:6-12");
        // One interval per contig: a gap, and the description names each contig.
        run("two-contigs", dir, fasta, "-L", "chr1:1-6", "-L", "chr2:1-6");
        // The soft-masked and IUPAC stretches, which come back upper-cased and as N.
        run("masked", dir, fasta, "-L", "chr1:13-24");
        run("iupac", dir, fasta, "-L", "chr1:25-36");
        // A narrow line width, so the wrapping is visible.
        run("narrow-lines", dir, fasta, "-L", "chr1:1-12", "--line-width", "5");
        // A line width of zero, which the writer refuses.
        run("zero-width", dir, fasta, "-L", "chr1:1-12", "--line-width", "0");
        // A single base, which is a sequence of one.
        run("one-base", dir, fasta, "-L", "chr1:7-7");
    }

    static void run(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".fasta");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new FastaReferenceMaker().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s%n", label, e.getClass().getName());
            return;
        }
        System.out.printf("fasta\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
        final Path fai = dir.resolve("out-" + label + ".fasta.fai");
        System.out.printf("fai\t%s\t%s%n", label, Files.exists(fai)
                ? ReferenceQueryDump.escape(new String(Files.readAllBytes(fai))) : "(none)");
        final Path dict = dir.resolve("out-" + label + ".dict");
        System.out.printf("dict\t%s\t%s%n", label, Files.exists(dict)
                ? ReferenceQueryDump.escape(new String(Files.readAllBytes(dict))) : "(none)");
    }
}
