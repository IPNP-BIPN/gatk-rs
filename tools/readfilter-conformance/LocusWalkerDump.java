/*
 * What a LocusWalker hands to apply(), taken from the reference.
 *
 * The traversal the largest GATK archetype inherits. The probe is a real LocusWalker run through
 * the real command line rather than a reconstruction of its parts, so the defaults it carries are
 * measured rather than transcribed. Four of them are not a ReadWalker's:
 *
 *   - two default read filters, not one: WellformedReadFilter *and* MappedReadFilter. An unmapped
 *     read carrying its mate's position reaches a ReadWalker and never reaches a LocusWalker;
 *   - includeDeletions is true and includeNs is false;
 *   - emitEmptyLoci is false, so an uncovered position inside -L produces no apply call at all;
 *   - --max-depth-per-sample defaults to 0, meaning no downsampling, and a negative value is a bad
 *     argument rather than a synonym for unlimited.
 *
 * The fixture is ReadWalkerDump's, so a divergence between the two suites is a divergence between
 * the traversals rather than between their inputs. It carries an unmapped read with a mate
 * position, which is exactly the read the second default filter removes.
 *
 * Output:
 *
 *     apply\t<label>\t<n>\t<contig>:<pos>|<depth>|<bases>|<reference base>
 *     summary\t<label>\t<ok|E:class>
 *     count\t<label>\t<number of apply calls>
 *
 * Usage: LocusWalkerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.LocusWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LocusWalkerDump {

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String> APPLIED = new ArrayList<>();

    @CommandLineProgramProperties(
            summary = "Records what a LocusWalker hands to apply()",
            oneLineSummary = "LocusWalker traversal probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static class ProbeWalker extends LocusWalker {
        @Override
        public void apply(final AlignmentContext context, final ReferenceContext reference,
                          final FeatureContext features) {
            final StringBuilder bases = new StringBuilder();
            for (final PileupElement element : context.getBasePileup()) {
                bases.append((char) element.getBase());
            }
            final byte[] referenceBases = reference.getBases();
            APPLIED.add(String.format("%s:%d|%d|%s|%s",
                    context.getContig(),
                    context.getPosition(),
                    context.size(),
                    bases.length() == 0 ? "-" : bases.toString(),
                    referenceBases.length == 0 ? "-" : new String(referenceBases)));
        }
    }

    /** The same walker with emitEmptyLoci overridden, which is a method rather than an argument. */
    @CommandLineProgramProperties(
            summary = "Records what a LocusWalker hands to apply(), emitting empty loci",
            oneLineSummary = "LocusWalker empty-loci probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class EmptyLociWalker extends ProbeWalker {
        @Override
        public boolean emitEmptyLoci() {
            return true;
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("locuswalker-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# LocusWalkerDump: what a LocusWalker hands to apply()");

        traverse("all", bam, fasta);
        traverse("chr1:100-130", bam, fasta, "-L", "chr1:100-130");
        // An interval whose middle is uncovered, so the default (no empty loci) skips positions.
        traverse("gap", bam, fasta, "-L", "chr1:20-70");
        // The same interval with emitEmptyLoci overridden, which fills the gap with zero-depth
        // contexts. It is a method a tool overrides, not a command-line argument.
        traverseWith(EmptyLociWalker::new, "gap-emptyloci", bam, fasta, "-L", "chr1:20-70");
        // Deletions excluded, which is the flag a tool overrides rather than a command-line one,
        // so this run exists to confirm the *default* is to include them.
        traverse("chr2", bam, fasta, "-L", "chr2");
        // No reference: the ReferenceContext arrives with an empty base rather than a null one.
        traverse("all-noref", bam, null);
        // A negative depth cap, which is a bad argument rather than unlimited.
        traverse("negative-depth", bam, fasta, "--max-depth-per-sample", "-1");
        // A depth cap that would downsample. Recorded because the port refuses it, and a refusal
        // has to be measured against what the reference actually does.
        traverse("depth-1", bam, fasta, "--max-depth-per-sample", "1");
    }

    static void traverse(final String label, final Path bam, final Path fasta,
                         final String... extra) {
        traverseWith(ProbeWalker::new, label, bam, fasta, extra);
    }

    static void traverseWith(final java.util.function.Supplier<ProbeWalker> factory,
                             final String label, final Path bam, final Path fasta,
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
            factory.get().instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E:" + e.getClass().getName();
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i));
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
    }
}
