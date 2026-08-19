/*
 * ShiftFasta, taken from the reference.
 *
 * The fifth member of the reference archetype, and the first that needs nothing new: it reads the
 * reference this port already reads and writes through the FASTA writer it already has. What it
 * produces is four files, and three of them are formats nothing else in the port writes.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE DEFAULT OFFSET IS HALF THE CONTIG, integer division, computed per contig rather than once
 *     for the reference;
 *   - THE SHIFTED SEQUENCE IS THE TWO HALVES SWAPPED, tail first: startSequence then appendBases of
 *     everything from the offset, then appendBases of everything before it, which is two calls and
 *     therefore a line break that does not fall at the join;
 *   - A CONTIG WHOSE OFFSET IS 0 OR THE WHOLE LENGTH IS SKIPPED ENTIRELY, not copied unshifted, so
 *     it is absent from the output FASTA and from the chain file;
 *   - THE CHAIN FILE IS TWO RECORDS PER CONTIG with a chain id that counts across contigs, and the
 *     score field of the first is the shift-back offset while the second's is the offset minus one;
 *   - THE INTERVAL FILES ARE OFF BY ONE FROM EACH OTHER BY THE CONTIG'S PARITY: the shifted end
 *     gains `contigLength % 2` and the unshifted one does not, and both starts are `shiftOffset/2`,
 *     which is zero for a contig shifted by one;
 *   - AND AN OFFSET LIST OF THE WRONG LENGTH IS A UserException.BadInput, checked before any contig
 *     is written.
 *
 * Output:
 *
 *     fasta\t<label>\t<the shifted FASTA, escaped>
 *     fai\t<label>\t<its .fai, escaped>
 *     dict\t<label>\t<its .dict, escaped>
 *     chain\t<label>\t<the shift-back chain file, escaped>
 *     intervals\t<label>\t<the .intervals file, escaped>
 *     shifted\t<label>\t<the .shifted.intervals file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ShiftFastaDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.fasta.ShiftFasta;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ShiftFastaDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("shiftfasta-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# ShiftFastaDump: a reference shifted, and the three files beside it");

        // No offsets: each contig is shifted by half its own length.
        run("halves", dir, fasta);
        // An explicit offset per contig.
        run("explicit", dir, fasta, "--shift-offset-list", "5", "--shift-offset-list", "7");
        // A zero offset, which skips that contig rather than copying it.
        run("skip-first", dir, fasta, "--shift-offset-list", "0", "--shift-offset-list", "7");
        // An offset equal to the contig's length, which is skipped for the same reason.
        run("skip-whole", dir, fasta, "--shift-offset-list", "43", "--shift-offset-list", "7");
        // A shift of one, whose interval start is zero.
        run("shift-one", dir, fasta, "--shift-offset-list", "1", "--shift-offset-list", "1");
        // A narrow line width, so the two appends and the wrapping are both visible.
        run("narrow-lines", dir, fasta, "--line-width", "7");
        // An offset list that does not match the number of contigs.
        run("wrong-length", dir, fasta, "--shift-offset-list", "5");
    }

    static void run(final String label, final Path dir, final Path fasta, final String... extra)
            throws Exception {
        final Path out = dir.resolve("shifted-" + label + ".fasta");
        final Path chain = dir.resolve("shifted-" + label + ".chain");
        final String intervals = dir.resolve("shifted-" + label).toString();
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString(),
                "--shift-back-output", chain.toString(),
                "--interval-file-name", intervals));
        argv.addAll(Arrays.asList(extra));
        try {
            new ShiftFasta().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        emit("fasta", label, out);
        emit("fai", label, dir.resolve("shifted-" + label + ".fasta.fai"));
        emit("dict", label, dir.resolve("shifted-" + label + ".dict"));
        emit("chain", label, chain);
        emit("intervals", label, Path.of(intervals + ".intervals"));
        emit("shifted", label, Path.of(intervals + ".shifted.intervals"));
    }

    static void emit(final String kind, final String label, final Path path) throws Exception {
        System.out.printf("%s\t%s\t%s%n", kind, label, Files.exists(path)
                ? ReferenceQueryDump.escape(new String(Files.readAllBytes(path))) : "(none)");
    }
}
