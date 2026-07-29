/*
 * ReferenceContext's answers, taken from the reference.
 *
 * A walker never queries the reference directly: it is handed a ReferenceContext, which carries an
 * interval (where the walker is) and a window (what it can see). Everything that compares a read to
 * the reference reads its bases through that window, so a window off by one base is an off-by-one
 * in every annotation, invisible in the tool's own output.
 *
 * Three of its behaviours are decisions rather than arithmetic, and the cases below are chosen to
 * expose each:
 *
 *   - near a contig edge the window is silently smaller than requested, and numWindowLeadingBases
 *     reports what was obtained rather than what was asked for;
 *   - getBases(leading, trailing) expands the *window*, not the interval, so the expansions
 *     compose: setWindow(10,10) followed by getBases(5,5) returns fifteen bases each side;
 *   - getKmerAround returns null, not a shorter kmer, when the contig edge prevents expansion.
 *
 * The fixture is ReferenceQueryDump's FASTA, shared rather than copied so the two suites describe
 * the same bytes.
 *
 * Output:
 *
 *     fasta\t<escaped FASTA text>
 *     fai\t<escaped .fai text>
 *     case\t<label>\t<what was built and asked>
 *     result\t<label>\t<answer, or E where the reference threw>
 *
 * Usage: ReferenceContextDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.Supplier;

public class ReferenceContextDump {

    static ReferenceDataSource source;

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("refcontext");
        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        final Path fai = dir.resolve("ref.fasta.fai");
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# ReferenceContextDump: what a walker sees around its position");
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(fai))));

        try (final ReferenceDataSource dataSource = new ReferenceFileSource(fasta)) {
            source = dataSource;

            // chr1 is 43 bases over four lines, chr2 is 24 over two. The interesting positions are
            // the first base, the last base, and the line boundaries in between.
            for (final int[] span : new int[][] {{1, 1}, {5, 5}, {12, 13}, {20, 25}, {43, 43}}) {
                for (final int[] window : new int[][] {{0, 0}, {3, 3}, {10, 10}, {50, 50},
                        {0, 5}, {5, 0}}) {
                    context("chr1", span[0], span[1], window[0], window[1]);
                }
            }
            for (final int[] span : new int[][] {{1, 1}, {24, 24}}) {
                context("chr2", span[0], span[1], 5, 5);
            }

            // A negative window is a GATKException, not a clamp.
            emit("negativeWindow", "new ReferenceContext(chr1:5-5, -1, 0)", () -> {
                new ReferenceContext(source, new SimpleInterval("chr1", 5, 5), -1, 0);
                return "no throw";
            });
            // A contig the reference does not have is a UserException from trimToContigLength.
            emit("unknownContig", "new ReferenceContext(chr9:5-5, 3, 3)", () -> {
                new ReferenceContext(source, new SimpleInterval("chr9", 5, 5), 3, 3);
                return "no throw";
            });
            // The window must contain the interval.
            emit("windowInsideInterval", "window chr1:6-7 inside interval chr1:5-10", () -> {
                new ReferenceContext(source, new SimpleInterval("chr1", 5, 10),
                        new SimpleInterval("chr1", 6, 7));
                return "no throw";
            });
            // No interval at all: every query answers empty rather than failing.
            emit("noInterval", "new ReferenceContext()", () -> {
                final ReferenceContext context = new ReferenceContext();
                return String.format("bases=%s lead=%d trail=%d backing=%b",
                        new String(context.getBases()),
                        context.numWindowLeadingBases(),
                        context.numWindowTrailingBases(),
                        context.hasBackingDataSource());
            });

            // The copy constructor carries the window *sizes*, so a context built at a contig edge
            // propagates the cropped size rather than the requested one.
            emit("copyFromEdge", "ReferenceContext(chr1:1-1 window 10,10) then interval chr1:20-20",
                    () -> {
                        final ReferenceContext edge = new ReferenceContext(
                                source, new SimpleInterval("chr1", 1, 1), 10, 10);
                        final ReferenceContext moved = new ReferenceContext(
                                edge, new SimpleInterval("chr1", 20, 20));
                        return String.format("edgeLead=%d edgeTrail=%d movedWindow=%s bases=%s",
                                edge.numWindowLeadingBases(), edge.numWindowTrailingBases(),
                                moved.getWindow(), new String(moved.getBases()));
                    });
        }
    }

    /**
     * One context, and every accessor on it as its own row.
     *
     * Each accessor is guarded separately on purpose: a composite row would collapse to `E` the
     * moment one call throws, and which call throws is exactly what the golden is for.
     */
    static void context(final String contig, final int start, final int end,
                        final int leading, final int trailing) {
        final String label = String.format("%s:%d-%d+%d,%d", contig, start, end, leading, trailing);
        System.out.printf("case\t%s\t%s%n", label, String.format(
                "interval %s:%d-%d, setWindow(%d, %d)", contig, start, end, leading, trailing));

        final Supplier<ReferenceContext> build = () -> new ReferenceContext(
                source, new SimpleInterval(contig, start, end), leading, trailing);

        result(label + "/window", () -> String.valueOf(build.get().getWindow()));
        result(label + "/lead", () -> String.valueOf(build.get().numWindowLeadingBases()));
        result(label + "/trail", () -> String.valueOf(build.get().numWindowTrailingBases()));
        result(label + "/bases", () -> new String(build.get().getBases()));
        result(label + "/forward", () -> new String(build.get().getForwardBases()));
        result(label + "/base", () -> String.valueOf((char) build.get().getBase()));
        result(label + "/expand5", () -> new String(build.get().getBases(5, 5)));
        result(label + "/expand0", () -> new String(build.get().getBases(0, 0)));
        result(label + "/kmer3", () -> String.valueOf(build.get().getKmerAround(start, 3)));
        result(label + "/kmer20", () -> String.valueOf(build.get().getKmerAround(start, 20)));
        // A window given explicitly rather than as two counts.
        result(label + "/explicit", () -> new String(new ReferenceContext(
                source, new SimpleInterval(contig, start, end),
                new SimpleInterval(contig, start - leading, end + trailing)).getBases()));
    }

    static void result(final String label, final Supplier<String> op) {
        String result;
        try {
            result = op.get();
        } catch (final Exception | AssertionError e) {
            result = "E";
        }
        System.out.printf("result\t%s\t%s%n", label, result);
    }

    static void emit(final String label, final String what, final Supplier<String> op) {
        String result;
        try {
            result = op.get();
        } catch (final Exception | AssertionError e) {
            result = "E";
        }
        System.out.printf("case\t%s\t%s%n", label, what);
        System.out.printf("result\t%s\t%s%n", label, result);
    }
}
