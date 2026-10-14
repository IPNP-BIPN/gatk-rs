/*
 * The callers' debug outputs and switches, end to end (G3.2m, #1347).
 *
 * --haplotype-debug-histogram-output counts, as HaplotypeCaller's engine shuts down, the haplotypes
 * each assembly recovered and the k-mer size of each graph that assembled variation; Mutect2 keeps
 * them and never prints them. --assembly-region-out writes each region the walker hands over, an end marker at its start and
 * then the region with its size, 1 when active and -1 when not, in IGV's TAB format. --keep-rg
 * filters every other read group's reads before genotyping, --just-determine-active-regions stops
 * HaplotypeCaller before anything (the ramped engine only once inactive regions have their
 * reference model) and --dont-genotype once the trimmed assembly is known, which GVCF mode refuses.
 *
 * Read set: HcEndToEndDump's. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     file\t<label>=<the debug file the run wrote, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcDebugOutputsDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcDebugOutputsDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-debug-outputs-dump").toAbsolutePath();
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

        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> hc = HaplotypeCaller::new;
        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> ramped = RampedHaplotypeCaller::new;
        final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> m2 = Mutect2::new;
        withFile(dir, "hc-region-out", hc, e2eBam, e2eFasta, "--assembly-region-out", List.of());
        withFile(dir, "hc-region-out-gvcf", hc, e2eBam, e2eFasta, "--assembly-region-out",
                List.of("-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"));
        withFile(dir, "ramped-region-out", ramped, e2eBam, e2eFasta, "--assembly-region-out", List.of("-L", "chr2"));
        withFile(dir, "mutect2-region-out", m2, e2eBam, e2eFasta, "--assembly-region-out", List.of());
        final String histogram = "--haplotype-debug-histogram-output";
        withFile(dir, "hc-histogram", hc, e2eBam, e2eFasta, histogram, List.of());
        withFile(dir, "hc-histogram-linked", hc, e2eBam, e2eFasta, histogram, List.of("--linked-de-bruijn-graph"));
        withFile(dir, "hc-histogram-empty", hc, e2eBam, e2eFasta, histogram, List.of("-L", "chr2:10-20"));
        withFile(dir, "ramped-histogram", ramped, e2eBam, e2eFasta, histogram, List.of());
        withFile(dir, "mutect2-histogram", m2, e2eBam, e2eFasta, histogram, List.of());
        HcAllAnnotationsDump.run(dir, "hc-keep-rg", hc, e2eBam, e2eFasta, List.of("--keep-rg", "rg1"));
        HcAllAnnotationsDump.run(dir, "ramped-keep-rg", ramped, e2eBam, e2eFasta, List.of("--keep-rg", "rg2"));
        HcAllAnnotationsDump.run(dir, "hc-just-active", hc, e2eBam, e2eFasta, List.of("--just-determine-active-regions"));
        HcAllAnnotationsDump.run(dir, "hc-just-active-gvcf", hc, e2eBam, e2eFasta,
                List.of("--just-determine-active-regions", "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"));
        HcAllAnnotationsDump.run(dir, "ramped-just-active-gvcf", ramped, e2eBam, e2eFasta,
                List.of("--just-determine-active-regions", "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"));
        HcAllAnnotationsDump.run(dir, "hc-dont-genotype", hc, e2eBam, e2eFasta, List.of("--dont-genotype"));
        HcAllAnnotationsDump.run(dir, "ramped-dont-genotype", ramped, e2eBam, e2eFasta, List.of("--dont-genotype"));
        HcAllAnnotationsDump.run(dir, "hc-dont-genotype-gvcf", hc, e2eBam, e2eFasta,
                List.of("--dont-genotype", "-ERC", "GVCF"));
    }

    /** A run writing one debug file, then the file. */
    static void withFile(final Path dir, final String label,
                         final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> tool,
                         final Path bam, final Path fasta, final String argument, final List<String> extra) throws Exception {
        final Path file = dir.resolve("debug-" + label + ".txt");
        HcAllAnnotationsDump.run(dir, label, tool, bam, fasta, HcPileupDump.concat(extra, List.of(argument, file.toString())));
        if (Files.exists(file)) {
            System.out.printf("file\t%s=%s%n", label, ReferenceQueryDump.escape(Files.readString(file)));
        }
    }
}
