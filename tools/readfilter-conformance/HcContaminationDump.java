/*
 * Contamination downsampling in HaplotypeCaller and RampedHaplotypeCaller, end to end (G3.2q, #1351).
 *
 * --contamination-fraction-to-filter and --contamination-fraction-per-sample-file: at each event,
 * after the reads are narrowed to those overlapping it, AlleleLikelihoods.contaminationDownsampling
 * removes (int) (informative reads * fraction) of each sample's reads, chosen by
 * AlleleBiasedDownsamplingUtils so the remaining allele counts look most like a het or a hom, the
 * reads of an allele drawn by MathUtils.sampleIndicesWithoutReplacement from
 * Utils.getRandomDataGenerator(). The annotations then read a fresh marginalization of the event's
 * reads unless --use-filtered-reads-for-annotations. The per-sample file is read by
 * loadContaminationFile, whose refusals are measured too.
 *
 * Nothing resets that generator between runs in one JVM, so each run here reseeds it where a
 * fresh JVM starts. Utils.getRandomGenerator(), which the callers also draw from, is left to run
 * on across the runs, as in every other dump.
 *
 * Read set: HcEndToEndDump's (samples s1 and s2). Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     file\t<name>=<a contamination file, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message, the dump's directory as <dir>>
 *
 * Usage: HcContaminationDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;

public class HcContaminationDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-contamination-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcContaminationDump: contamination downsampling in the callers");

        // HcEndToEndDump's contigs and reads, from the same seed.
        final Random random = new Random(20261010L);
        final String[] references = new String[HcEndToEndDump.CONTIGS.length];
        for (int c = 0; c < HcEndToEndDump.CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < HcEndToEndDump.LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\te2e:%s=%s%n", HcEndToEndDump.CONTIGS[c], references[c]);
        }
        final Path e2eDir = dir.resolve("e2e");
        Files.createDirectories(e2eDir);
        final Path fasta = HcEndToEndDump.writeReference(e2eDir, references);
        final Path bam = HcLinkedGraphDump.bam(e2eDir, "e2e", HcEndToEndDump.sam(references, random));

        final Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> hc = HaplotypeCaller::new;
        final Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> ramped = RampedHaplotypeCaller::new;
        final String perSample = file(dir, "per-sample", "s1\t0.2\ns2\t0.0\nother\t0.5\n");
        final String spaced = file(dir, "spaced", "  s2\t0.3  \n\n");
        final String threeColumns = file(dir, "three-columns", "s1\t0.1\tx\n");
        final String emptyField = file(dir, "empty-field", "\t0.1\n");
        final String unparsable = file(dir, "unparsable", "s1\tlots\n");
        final String duplicate = file(dir, "duplicate", "s1\t0.1\ns1\t0.2\n");
        final String duplicateOther = file(dir, "duplicate-other", "x\t0.1\nx\t0.2\ns2\t0.15\n");
        final String outOfRange = file(dir, "out-of-range", "s1\t1.5\n");
        final String notANumber = file(dir, "nan", "s1\tNaN\n");
        final String fraction = "--contamination-fraction-to-filter";
        final String perSampleFile = "--contamination-fraction-per-sample-file";

        run(dir, "hc-0.05", hc, bam, fasta, List.of(fraction, "0.05"));
        run(dir, "hc-0.1", hc, bam, fasta, List.of(fraction, "0.1"));
        run(dir, "hc-0.3", hc, bam, fasta, List.of(fraction, "0.3"));
        run(dir, "hc-all", hc, bam, fasta, List.of(fraction, "1.0"));
        run(dir, "hc-zero", hc, bam, fasta, List.of(fraction, "0.0"));
        run(dir, "hc-filtered-annotations", hc, bam, fasta,
                List.of(fraction, "0.2", "--use-filtered-reads-for-annotations"));
        run(dir, "hc-gvcf", hc, bam, fasta,
                List.of(fraction, "0.2", "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        run(dir, "hc-per-sample", hc, bam, fasta, List.of(perSampleFile, perSample));
        run(dir, "hc-per-sample-default", hc, bam, fasta, List.of(perSampleFile, spaced, fraction, "0.1"));
        run(dir, "hc-duplicate-other", hc, bam, fasta, List.of(perSampleFile, duplicateOther));
        run(dir, "hc-nan", hc, bam, fasta, List.of(perSampleFile, notANumber));
        run(dir, "hc-three-columns", hc, bam, fasta, List.of(perSampleFile, threeColumns));
        run(dir, "hc-empty-field", hc, bam, fasta, List.of(perSampleFile, emptyField));
        run(dir, "hc-unparsable", hc, bam, fasta, List.of(perSampleFile, unparsable));
        run(dir, "hc-duplicate", hc, bam, fasta, List.of(perSampleFile, duplicate));
        run(dir, "hc-out-of-range", hc, bam, fasta, List.of(perSampleFile, outOfRange));
        run(dir, "hc-missing", hc, bam, fasta, List.of(perSampleFile, dir.resolve("absent.txt").toString()));
        run(dir, "ramped-0.2", ramped, bam, fasta, List.of(fraction, "0.2"));
        run(dir, "ramped-gvcf", ramped, bam, fasta,
                List.of(fraction, "0.2", "-ERC", "GVCF", "--sample-name", "s2", "-L", "chr2"));
    }

    /** One run, from where a fresh JVM's random data generator starts. */
    static void run(final Path dir, final String label, final Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> tool,
                    final Path bam, final Path fasta, final List<String> extra) throws Exception {
        Utils.getRandomDataGenerator().reSeed(47382911L);  // Utils.GATK_RANDOM_SEED, private
        HcAllAnnotationsDump.run(dir, label, tool, bam, fasta, new ArrayList<>(extra));
    }

    /** A contamination file, printed. */
    static String file(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name + ".txt");
        Files.writeString(path, text);
        System.out.printf("file\t%s=%s%n", name, ReferenceQueryDump.escape(text));
        return path.toString();
    }
}
