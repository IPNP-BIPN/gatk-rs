/*
 * What a ReferenceWalker hands to apply(), and what CountBasesInReference makes of it, taken from
 * the reference.
 *
 * The first member of the reference-utility archetype. The walker is thirty lines and everything
 * interesting about it is in what it stands on: which intervals a tool with NO required interval
 * argument traverses, how they are split, and which bases come back.
 *
 * Six behaviours this is built to catch.
 *
 *   - AN ABSENT -L IS THE WHOLE REFERENCE, one apply() per base of every contig in the sequence
 *     dictionary. IntervalWalker cannot reach that branch, because its requiresIntervals() is true
 *     and Barclay rejects the run before any interval is parsed. This walker's is false, so the
 *     branch is reachable and measured here for the first time;
 *   - -XL ALONE IS THEREFORE LEGAL, and subtracts from that whole-reference list;
 *   - THE INTERVALS ARE SORTED INTO DICTIONARY ORDER, not left in the order they were given, so
 *     -L chr2 -L chr1 traverses chr1 first;
 *   - EVERY LOCUS IS ONE BASE, IntervalLocusIterator at size one, so an interval of eleven bases is
 *     eleven apply() calls and not one;
 *   - THE BASE IS THE CACHING READER'S, not the file's: lower case comes back upper-cased and every
 *     IUPAC ambiguity code comes back as N. That is what makes CountBasesInReference's answer a
 *     count of five symbols rather than of fifteen;
 *   - AND getReferenceWindow IS A METHOD, not an argument. A walker that widens it sees more bases
 *     per locus while still being called once per base, and the widened window is clipped at the
 *     contig ends rather than running off them.
 *
 * The FASTA is ReferenceQueryDump's, which was built awkward on purpose: mixed case, a line width
 * that does not divide the length, an N run, and IUPAC codes.
 *
 * Output:
 *
 *     apply\t<label>\t<index>\t<contig>:<start>-<end>|<bases>
 *     count\t<label>\t<number of apply calls>
 *     summary\t<label>\t<ok|E:class>
 *     counts\t<label>\t<what CountBasesInReference wrote, escaped>
 *
 * Usage: ReferenceWalkerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceWalker;
import org.broadinstitute.hellbender.tools.walkers.fasta.CountBasesInReference;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReferenceWalkerDump {

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String> APPLIED = new ArrayList<>();

    @CommandLineProgramProperties(
            summary = "Records what a ReferenceWalker hands to apply()",
            oneLineSummary = "ReferenceWalker traversal probe",
            programGroup = ReferenceProgramGroup.class)
    public static class ProbeWalker extends ReferenceWalker {
        @Override
        public void apply(final ReferenceContext reference, final ReadsContext reads,
                          final FeatureContext features) {
            final SimpleInterval window = reference.getWindow();
            APPLIED.add(String.format("%s:%d-%d|%s", window.getContig(), window.getStart(),
                    window.getEnd(), new String(reference.getBases())));
        }
    }

    /** The same walker with a widened window, which is a method rather than an argument. */
    @CommandLineProgramProperties(
            summary = "Records what a ReferenceWalker with a widened window hands to apply()",
            oneLineSummary = "ReferenceWalker window probe",
            programGroup = ReferenceProgramGroup.class)
    public static final class WideWalker extends ProbeWalker {
        @Override
        protected SimpleInterval getReferenceWindow(final SimpleInterval locus) {
            return new SimpleInterval(locus.getContig(),
                    Math.max(1, locus.getStart() - 2), locus.getEnd() + 2);
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("referencewalker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        // GATK refuses a reference with no sequence dictionary, so the harness makes one the same
        // way a user would.
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# ReferenceWalkerDump: what a ReferenceWalker hands to apply()");

        // No -L at all, which is the branch IntervalWalker cannot reach.
        traverse("all", fasta);
        traverse("chr1-window", fasta, "-L", "chr1:5-15");
        traverse("chr2", fasta, "-L", "chr2");
        // Given out of dictionary order, and traversed in it.
        traverse("out-of-order", fasta, "-L", "chr2:3-5", "-L", "chr1:1-3");
        // Excluded without any -L: legal here, and a whole-reference list minus the exclusion.
        traverse("excluded-only", fasta, "-XL", "chr1:5-40");
        // Padding, which widens the interval before it is split into loci.
        traverse("padded", fasta, "-L", "chr1:10-12", "--interval-padding", "2");
        // Two intervals that touch, which the default merging rule joins into one.
        traverse("abutting", fasta, "-L", "chr1:1-3", "-L", "chr1:4-6");
        // An interval that runs past the contig, which is a user error rather than a clip.
        traverse("past-the-end", fasta, "-L", "chr2:20-40");

        // The widened window: still one apply per base, but each sees five bases, clipped at the
        // contig ends.
        traverseWith(WideWalker::new, "wide-window", fasta, "-L", "chr1:1-4");

        // And the tool itself.
        countBases("all", fasta);
        countBases("chr1-window", fasta, "-L", "chr1:5-15");
        // The IUPAC run alone, which is where the reader's N substitution shows up in the counts.
        countBases("iupac", fasta, "-L", "chr1:25-36");
    }

    static void traverse(final String label, final Path fasta, final String... extra) {
        traverseWith(ProbeWalker::new, label, fasta, extra);
    }

    static void traverseWith(final java.util.function.Supplier<ProbeWalker> factory,
                             final String label, final Path fasta, final String... extra) {
        APPLIED.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-R", fasta.toString()));
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            factory.get().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E:" + e.getClass().getName();
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i));
        }
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
        System.out.printf("summary\t%s\t%s%n", label, summary);
    }

    static void countBases(final String label, final Path fasta, final String... extra)
            throws Exception {
        final Path out = fasta.getParent().resolve("counts-" + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        String text;
        try {
            new CountBasesInReference().instanceMain(argv.toArray(new String[0]));
            text = new String(Files.readAllBytes(out));
        } catch (final Exception | AssertionError e) {
            text = "E:" + e.getClass().getName();
        }
        System.out.printf("counts\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
    }
}
