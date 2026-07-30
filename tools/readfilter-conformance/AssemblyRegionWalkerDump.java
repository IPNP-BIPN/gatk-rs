/*
 * What an AssemblyRegionWalker hands to apply(), taken from the reference.
 *
 * The last of the five walker archetypes, and the one HaplotypeCaller and Mutect2 are built on. The
 * probe is a real AssemblyRegionWalker run through the real command line rather than a
 * reconstruction of its parts, so the defaults it carries are measured rather than transcribed:
 * min 50, max 300, padding 100, threshold 0.002, propagation 50, all from an
 * AssemblyRegionArgumentCollection the walker instantiates without overriding anything.
 *
 * Four behaviours are what this suite is for.
 *
 *   - one read shard per CONTIG, not one per interval. Two -L arguments on the same contig share a
 *     shard and can therefore share a region; two on different contigs cannot;
 *   - apply receives contexts over the region's PADDED span, not its active span, so a tool reads
 *     reference bases outside the territory it is allowed to call in;
 *   - --force-active rewrites isActive AFTER the regions have been cut, so it changes the flag on
 *     every region without changing a single boundary. A port that folded it into the evaluator
 *     would produce one region covering everything instead;
 *   - the default read filters are WellformedReadFilter and MappedReadFilter, the same pair as a
 *     LocusWalker's and not a ReadWalker's single filter.
 *
 * The fixture is ReadWalkerDump's, so a divergence between this suite and the locus or read walker
 * suites is a divergence between the traversals rather than between their inputs.
 *
 * Output:
 *
 *     apply\t<label>\t<n>\t<span>|<paddedSpan>|<isActive>|<nReads>|<refBases>
 *     aread\t<label>\t<n>\t<read names, comma-separated>
 *     summary\t<label>\t<ok|E:class>
 *     count\t<label>\t<number of apply calls>
 *
 * Usage: AssemblyRegionWalkerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.AssemblyRegionEvaluator;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.AssemblyRegionWalker;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.function.Supplier;

public class AssemblyRegionWalkerDump {

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String[]> APPLIED = new ArrayList<>();

    @CommandLineProgramProperties(
            summary = "Records what an AssemblyRegionWalker hands to apply()",
            oneLineSummary = "AssemblyRegionWalker traversal probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static class ProbeWalker extends AssemblyRegionWalker {

        @Override
        public AssemblyRegionEvaluator assemblyRegionEvaluator() {
            // A probe, not a caller: activity is a declared function of the pileup depth, so the
            // regions are decided by the traversal rather than by anything statistical.
            return (pileup, reference, features) -> new ActivityProfileState(
                    new SimpleInterval(pileup.getContig(),
                            (int) pileup.getPosition(), (int) pileup.getPosition()),
                    pileup.getBasePileup().size() >= 1 ? 1.0 : 0.0);
        }

        @Override
        public boolean shouldTrackPileupsForAssemblyRegions() {
            return false;
        }

        @Override
        public void apply(final AssemblyRegion region, final ReferenceContext reference,
                          final FeatureContext features) {
            final StringJoiner names = new StringJoiner(",");
            for (final GATKRead read : region.getReads()) {
                names.add(read.getName());
            }
            final byte[] bases = reference.getBases();
            APPLIED.add(new String[] {
                String.format("%s|%s|%b|%d|%d", region.getSpan(), region.getPaddedSpan(),
                        region.isActive(), region.size(), bases.length),
                names.toString(),
            });
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("assemblyregionwalker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# AssemblyRegionWalkerDump: what an AssemblyRegionWalker hands to apply()");

        // The defaults, over the whole reference and over each contig.
        traverse("all", bam, fasta);
        traverse("chr1", bam, fasta, "-L", "chr1");
        traverse("chr2", bam, fasta, "-L", "chr2");
        // Two intervals on one contig: one shard, and the padding can join them.
        traverse("two-on-one-contig", bam, fasta, "-L", "chr1:10-40", "-L", "chr1:150-190");
        // Two intervals on two contigs: two shards, so nothing can join them.
        traverse("two-contigs", bam, fasta, "-L", "chr1:10-40", "-L", "chr2:10-40");
        // A narrow interval, where the padding reaches well outside it.
        traverse("narrow", bam, fasta, "-L", "chr1:100-110");

        // The region-size arguments, at their defaults and away from them.
        traverse("small-regions", bam, fasta, "-L", "chr1",
                "--min-assembly-region-size", "5", "--max-assembly-region-size", "20");
        traverse("min-above-max", bam, fasta, "-L", "chr1",
                "--min-assembly-region-size", "300", "--max-assembly-region-size", "50");
        traverse("zero-padding", bam, fasta, "-L", "chr1", "--assembly-region-padding", "0");
        traverse("large-padding", bam, fasta, "-L", "chr1:100-110",
                "--assembly-region-padding", "500");

        // force-active, which rewrites the flag after the regions have been cut.
        traverse("force-active", bam, fasta, "-L", "chr1", "--force-active", "true");
        // A threshold above every probability the probe emits, so nothing is active.
        traverse("threshold-above-all", bam, fasta, "-L", "chr1",
                "--active-probability-threshold", "2.0");
        // The two together, which is the only way to see force-active do anything: at the default
        // threshold every region is already active, so the flag it rewrites was already true.
        traverse("force-active-above-threshold", bam, fasta, "-L", "chr1",
                "--active-probability-threshold", "2.0", "--force-active", "true");
        // The same pair with small regions, so the flag is rewritten on several regions at once
        // and the boundaries can be compared against the run that did not force anything.
        traverse("force-active-small-regions", bam, fasta, "-L", "chr1",
                "--active-probability-threshold", "2.0", "--force-active", "true",
                "--min-assembly-region-size", "5", "--max-assembly-region-size", "20");
        traverse("threshold-above-all-small-regions", bam, fasta, "-L", "chr1",
                "--active-probability-threshold", "2.0",
                "--min-assembly-region-size", "5", "--max-assembly-region-size", "20");
        // A propagation distance of zero, which changes when a region may be popped at all.
        traverse("zero-propagation", bam, fasta, "-L", "chr1",
                "--max-prob-propagation-distance", "0");

        // The downsampler, which the walker creates only when the argument is above zero.
        traverse("max-starts-1", bam, fasta, "-L", "chr1", "--max-reads-per-alignment-start", "1");
        traverse("max-starts-0", bam, fasta, "-L", "chr1", "--max-reads-per-alignment-start", "0");
        traverse("max-starts-negative", bam, fasta, "-L", "chr1",
                "--max-reads-per-alignment-start", "-1");

        // The default filters, confirmed by disabling them: the reads that reappear are the ones
        // WellformedReadFilter and MappedReadFilter were removing.
        traverse("no-filters", bam, fasta, "-L", "chr1", "--disable-tool-default-read-filters");
    }

    static void traverse(final String label, final Path bam, final Path fasta,
                         final String... extra) {
        traverseWith(ProbeWalker::new, label, bam, fasta, extra);
    }

    static void traverseWith(final Supplier<ProbeWalker> factory, final String label,
                             final Path bam, final Path fasta, final String... extra) {
        APPLIED.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-I", bam.toString()));
        if (fasta != null) {
            argv.add("-R");
            argv.add(fasta.toString());
        }
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            factory.get().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            // The message travels too: "no-filters" throws an IllegalStateException whose text is
            // the only thing that says which layer refused, and a class name alone would leave that
            // unattributed.
            summary = "E:" + e.getClass().getName() + ":" + oneLine(e.getMessage());
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i)[0]);
            System.out.printf("aread\t%s\t%d\t%s%n", label, i, APPLIED.get(i)[1]);
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
