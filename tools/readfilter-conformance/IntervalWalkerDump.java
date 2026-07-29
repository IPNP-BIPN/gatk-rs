/*
 * What an IntervalWalker is handed, taken from the reference.
 *
 * IntervalWalker.traverse() is five lines: it iterates userIntervals and calls apply on each. The
 * behaviour worth measuring is therefore not the loop but the list, which IntervalArgumentCollection
 * builds out of -L, -XL, --interval-padding, --interval-exclusion-padding, --interval-set-rule and
 * --interval-merging-rule. Every case below is a combination of those, and the row records the
 * intervals that actually reached apply.
 *
 * The cases are chosen for the decisions that are invisible from the argument names:
 *
 *   - padding is clamped to the contig, not extended past it, and an interval padded entirely off
 *     the contig disappears rather than throwing;
 *   - padding merges with ALL regardless of --interval-merging-rule, because getIntervalsWithFlanks
 *     hard-codes it, so two padded arguments that come to abut are joined even under
 *     OVERLAPPING_ONLY;
 *   - INTERSECTION is a running fold whose accumulator starts empty, and mergeListsBySetOperator
 *     short-circuits on an empty side, so a single -L under INTERSECTION is the identity rather
 *     than an empty result;
 *   - an empty intersection is an error, not an empty traversal, and so is a -XL that removes
 *     everything -L asked for;
 *   - -XL with no -L is subtraction from the whole reference;
 *   - -XL inside an interval splits it in two.
 *
 * Each case also records what apply saw beyond the interval: how many reads the ReadsContext
 * returned and what window the ReferenceContext arrived with, so the walker's three arguments are
 * measured together rather than the interval alone.
 *
 * Output:
 *
 *     apply\t<label>\t<n>\t<interval>|<reads>|<reference window>|<bases>
 *     summary\t<label>\t<ok|E>
 *     count\t<label>\t<number of apply calls>
 *
 * The fixture is ReadWalkerDump's, reused rather than rebuilt: the two suites then disagree about
 * a traversal only when the traversals differ.
 *
 * Usage: IntervalWalkerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.IntervalWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class IntervalWalkerDump {

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String> APPLIED = new ArrayList<>();

    @CommandLineProgramProperties(
            summary = "Records what an IntervalWalker hands to apply()",
            oneLineSummary = "IntervalWalker traversal probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class ProbeWalker extends IntervalWalker {
        @Override
        public void apply(final SimpleInterval interval, final ReadsContext reads,
                          final ReferenceContext reference, final FeatureContext features) {
            int readCount = 0;
            for (final GATKRead ignored : reads) {
                readCount++;
            }
            APPLIED.add(String.format("%s|%d|%s|%s",
                    interval.toString(),
                    readCount,
                    reference.getWindow() == null ? "null" : reference.getWindow().toString(),
                    new String(reference.getBases())));
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("intervalwalker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# IntervalWalkerDump: what an IntervalWalker hands to apply()");
        // The fixture travels in this golden too, rather than being borrowed from the ReadWalker
        // one. It is the same fixture, built by the same method, but a suite that reads its input
        // out of another suite's golden is a suite that passes because the other one changed.
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(dir.resolve("ref.fasta.fai")))));
        System.out.printf("bam\t%s%n", Base64.getEncoder().encodeToString(
                Files.readAllBytes(bam)));
        System.out.printf("bai\t%s%n", Base64.getEncoder().encodeToString(
                Files.readAllBytes(dir.resolve("reads.bai"))));

        // One whole contig, and one interval inside it.
        traverse("chr1", bam, fasta, "-L", "chr1");
        traverse("chr1:100-160", bam, fasta, "-L", "chr1:100-160");
        // Dictionary order decides, not the order the arguments were given in.
        traverse("chr2-then-chr1", bam, fasta, "-L", "chr2", "-L", "chr1");

        // Merging: abutting, overlapping, and separated by exactly one base.
        traverse("abutting-all", bam, fasta, "-L", "chr1:1-100", "-L", "chr1:101-200");
        traverse("abutting-overlapping-only", bam, fasta,
                "-L", "chr1:1-100", "-L", "chr1:101-200",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        traverse("overlapping-overlapping-only", bam, fasta,
                "-L", "chr1:1-100", "-L", "chr1:50-150",
                "--interval-merging-rule", "OVERLAPPING_ONLY");
        traverse("one-base-gap", bam, fasta, "-L", "chr1:1-10", "-L", "chr1:12-20");

        // Padding: ordinary, clamped at the start of the contig, clamped at its end.
        traverse("padded", bam, fasta, "-L", "chr1:50-60", "--interval-padding", "20");
        traverse("padded-clamped-start", bam, fasta, "-L", "chr1:1-5", "--interval-padding", "20");
        traverse("padded-clamped-end", bam, fasta, "-L", "chr1:195-200",
                "--interval-padding", "20");
        // Padding merges with ALL inside getIntervalsWithFlanks whatever the merging rule says, so
        // these two arguments come out as one interval even here.
        traverse("padded-overlapping-only", bam, fasta,
                "-L", "chr1:1-50", "-L", "chr1:60-100", "--interval-padding", "5",
                "--interval-merging-rule", "OVERLAPPING_ONLY");

        // The set rule, including the fold-over-one-argument identity and the empty intersection.
        traverse("intersection", bam, fasta,
                "-L", "chr1:1-100", "-L", "chr1:50-150", "--interval-set-rule", "INTERSECTION");
        traverse("intersection-single", bam, fasta,
                "-L", "chr1:1-100", "--interval-set-rule", "INTERSECTION");
        traverse("intersection-empty", bam, fasta,
                "-L", "chr1", "-L", "chr2", "--interval-set-rule", "INTERSECTION");

        // Exclusion: from the whole reference, from one interval's middle, from either end, and
        // taking everything away.
        traverse("exclude-only", bam, fasta, "-XL", "chr1");
        traverse("exclude-middle", bam, fasta, "-L", "chr1", "-XL", "chr1:50-100");
        traverse("exclude-head", bam, fasta, "-L", "chr1", "-XL", "chr1:1-50");
        traverse("exclude-tail", bam, fasta, "-L", "chr1", "-XL", "chr1:150-200");
        traverse("exclude-everything", bam, fasta, "-L", "chr1", "-XL", "chr1");
        traverse("exclude-other-contig", bam, fasta, "-L", "chr1", "-L", "chr2", "-XL", "chr2");
        traverse("exclude-padded", bam, fasta, "-L", "chr1:1-100", "-XL", "chr1:40-60",
                "--interval-exclusion-padding", "10");

        // -L unmapped on a walker that requires intervals: the request is separated out of the
        // interval list, so what is left may be nothing at all.
        traverse("unmapped", bam, fasta, "-L", "unmapped");
        traverse("unmapped-and-chr1", bam, fasta, "-L", "unmapped", "-L", "chr1:1-20");

        // No reference: the window is the interval and only the bases are empty.
        traverse("chr1-noref", bam, null, "-L", "chr1:100-160");
    }

    static void traverse(final String label, final Path bam, final Path fasta,
                         final String... extra) {
        APPLIED.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-I", bam.toString()));
        if (fasta != null) {
            argv.add("-R");
            argv.add(fasta.toString());
        }
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            new ProbeWalker().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            // The class is recorded, the message is not: the message carries the run's paths and
            // would make the golden unreproducible, while the class distinguishes a rejected
            // argument from a failed traversal, which is a difference a port can get wrong.
            summary = "E:" + e.getClass().getName();
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i));
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
    }
}
