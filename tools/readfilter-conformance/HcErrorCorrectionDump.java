/*
 * Read error correction in the callers, end to end (G3.2n, #1348).
 *
 * --error-correction-log-odds turns on PileupReadErrorCorrector: a copy of each region's reads is
 * piled up, and where Mutect2's log likelihood ratio of the bases other than the plurality is below
 * the threshold they are corrected to it for the assembly, except near indels and soft clips and at
 * read ends with three edits in fifteen bases. --error-correct-reads builds the k-mer corrector,
 * which is never given a k-mer and corrects nothing, but finalizes the region with a minimum tail
 * quality of 6 and checks --kmer-length-for-read-error-correction.
 *
 * Read set: HcEndToEndDump's, then the same reads with three per cent of their bases replaced by
 * another at a random quality, from a fixed seed. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcErrorCorrectionDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcErrorCorrectionDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-error-correction-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcErrorCorrectionDump: read error correction in the callers");

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
        final Path e2eFasta = HcEndToEndDump.writeReference(e2eDir, references);
        final String e2eSam = HcEndToEndDump.sam(references, random);
        final Path e2eBam = HcLinkedGraphDump.bam(e2eDir, "e2e", e2eSam);
        final Path noisyDir = dir.resolve("noisy");
        Files.createDirectories(noisyDir);
        final Path noisyFasta = HcEndToEndDump.writeReference(noisyDir, references);
        final Path noisyBam = HcLinkedGraphDump.bam(noisyDir, "noisy", noisy(e2eSam, new Random(20261014L)));

        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> hc = HaplotypeCaller::new;
        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> ramped = RampedHaplotypeCaller::new;
        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> m2 = Mutect2::new;
        final String odds = "--error-correction-log-odds";
        HcAllAnnotationsDump.run(dir, "hc-noisy", hc, noisyBam, noisyFasta, List.of());
        HcAllAnnotationsDump.run(dir, "hc-kmer-corrector", hc, noisyBam, noisyFasta, List.of("--error-correct-reads"));
        HcAllAnnotationsDump.run(dir, "hc-kmer-corrector-tuned", hc, noisyBam, noisyFasta,
                List.of("--error-correct-reads", "--kmer-length-for-read-error-correction", "10",
                        "--min-observations-for-kmer-to-be-solid", "2"));
        HcAllAnnotationsDump.run(dir, "hc-kmer-corrector-zero", hc, noisyBam, noisyFasta,
                List.of("--error-correct-reads", "--kmer-length-for-read-error-correction", "0"));
        HcAllAnnotationsDump.run(dir, "hc-kmer-corrector-zero-pileup", hc, noisyBam, noisyFasta,
                List.of("--error-correct-reads", "--kmer-length-for-read-error-correction", "0", odds, "0.0"));
        for (final String threshold : List.of("-5.0", "0.0", "2.0", "10.0", "1000.0")) {
            HcAllAnnotationsDump.run(dir, "hc-pileup-corrector-" + threshold, hc, noisyBam, noisyFasta, List.of(odds, threshold));
        }
        HcAllAnnotationsDump.run(dir, "hc-pileup-corrector-clean", hc, e2eBam, e2eFasta, List.of(odds, "1000.0"));
        HcAllAnnotationsDump.run(dir, "hc-pileup-corrector-gvcf", hc, noisyBam, noisyFasta,
                List.of(odds, "10.0", "--error-correct-reads", "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        HcAllAnnotationsDump.run(dir, "hc-pileup-corrector-linked", hc, noisyBam, noisyFasta,
                List.of(odds, "10.0", "--linked-de-bruijn-graph"));
        HcAllAnnotationsDump.run(dir, "ramped-pileup-corrector", ramped, noisyBam, noisyFasta, List.of(odds, "10.0"));
        HcAllAnnotationsDump.run(dir, "ramped-kmer-corrector", ramped, noisyBam, noisyFasta, List.of("--error-correct-reads"));
        HcAllAnnotationsDump.run(dir, "mutect2-noisy", m2, noisyBam, noisyFasta, List.of());
        HcAllAnnotationsDump.run(dir, "mutect2-pileup-corrector", m2, noisyBam, noisyFasta, List.of(odds, "10.0"));
        HcAllAnnotationsDump.run(dir, "mutect2-pileup-corrector-tn", m2, noisyBam, noisyFasta,
                List.of(odds, "2.0", "--normal-sample", "s2"));
        HcAllAnnotationsDump.run(dir, "mutect2-kmer-corrector", m2, noisyBam, noisyFasta, List.of("--error-correct-reads"));
        HcAllAnnotationsDump.run(dir, "mutect2-kmer-corrector-zero", m2, noisyBam, noisyFasta,
                List.of("--error-correct-reads", "--kmer-length-for-read-error-correction", "0"));
    }

    /**
     * The SAM text with three per cent of each read's bases replaced by another base, at a quality
     * drawn from 5 to 34. Header lines are kept as they are.
     */
    static String noisy(final String sam, final Random random) {
        final StringBuilder out = new StringBuilder();
        for (final String line : sam.split("\n")) {
            if (line.startsWith("@")) {
                out.append(line).append('\n');
                continue;
            }
            final String[] fields = line.split("\t", -1);
            final char[] bases = fields[9].toCharArray();
            final char[] quals = fields[10].toCharArray();
            for (int i = 0; i < bases.length; i++) {
                if (random.nextInt(100) < 3) {
                    final String others = "ACGT".replace(String.valueOf(bases[i]), "");
                    bases[i] = others.charAt(random.nextInt(others.length()));
                    if (quals.length == bases.length) {
                        quals[i] = (char) (33 + 5 + random.nextInt(30));
                    }
                }
            }
            fields[9] = new String(bases);
            fields[10] = new String(quals);
            out.append(String.join("\t", fields)).append('\n');
        }
        return out.toString();
    }
}
