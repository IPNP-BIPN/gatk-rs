/*
 * The callers' --bam-output, end to end (G3.2m, #1347).
 *
 * The writer adds, for each region called, its haplotypes as reads (HC_<n>_<index>, qualities of
 * 33, the HC tag holding the haplotype's hash, the artificial read group, CR the callable region)
 * and then its reads as scored: realigned to their best haplotype, the soft-clip tags dropped, AR
 * the location of their best haplotype and XA the alleles they were informative for. Not
 * presorted, so htsjdk sorts by coordinate as it closes. Each writer type, GVCF mode,
 * --disable-optimizations, the PDHMM, RampedHaplotypeCaller, and Mutect2.
 *
 * Read set: HcEndToEndDump's. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     bam\t<label>=<the bamout as SAM text, header included, escaped>
 *     bai\t<label>=<whether an index was written beside it>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcBamoutDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcBamoutDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-bamout-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcBamoutDump: the callers' --bam-output");

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
        final Path e2eBam = HcLinkedGraphDump.bam(e2eDir, "e2e", HcEndToEndDump.sam(references, random));

        bamout(dir, "hc-bamout", HaplotypeCaller::new, e2eBam, e2eFasta, List.of());
        bamout(dir, "hc-bamout-all", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--bam-writer-type", "ALL_POSSIBLE_HAPLOTYPES"));
        bamout(dir, "hc-bamout-no-haplotypes", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--bam-writer-type", "NO_HAPLOTYPES", "--create-output-bam-index", "false"));
        bamout(dir, "hc-bamout-no-reads", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--bam-writer-type", "CALLED_HAPLOTYPES_NO_READS"));
        bamout(dir, "hc-bamout-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        bamout(dir, "hc-bamout-disable-optimizations", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--disable-optimizations", "-L", "chr2"));
        bamout(dir, "hc-bamout-pdhmm", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--use-pdhmm", "--bam-writer-type", "ALL_POSSIBLE_HAPLOTYPES", "-L", "chr1:1100-1350"));
        bamout(dir, "ramped-bamout", RampedHaplotypeCaller::new, e2eBam, e2eFasta, List.of());
        bamout(dir, "mutect2-bamout", Mutect2::new, e2eBam, e2eFasta, List.of());
        bamout(dir, "mutect2-bamout-all-tn", Mutect2::new, e2eBam, e2eFasta,
                List.of("--bam-writer-type", "ALL_POSSIBLE_HAPLOTYPES", "--normal-sample", "s2"));
    }

    /** A run with --bam-output, then the bamout as SAM text and whether it was indexed. */
    static void bamout(final Path dir, final String label,
                       final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> tool,
                       final Path bam, final Path fasta, final List<String> extra) throws Exception {
        final Path out = dir.resolve("bamout-" + label + ".bam");
        HcAllAnnotationsDump.run(dir, label, tool, bam, fasta,
                HcPileupDump.concat(extra, List.of("--bam-output", out.toString())));
        if (!Files.exists(out)) {
            return;
        }
        final Path sam = dir.resolve("bamout-" + label + ".sam");
        try (htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault()
                     .validationStringency(htsjdk.samtools.ValidationStringency.SILENT).open(out);
             htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory()
                     .makeSAMWriter(reader.getFileHeader(), true, sam.toFile())) {
            for (final htsjdk.samtools.SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("bam\t%s=%s%n", label, ReferenceQueryDump.escape(Files.readString(sam)));
        System.out.printf("bai\t%s=%s%n", label, Files.exists(dir.resolve("bamout-" + label + ".bai")));
    }
}
