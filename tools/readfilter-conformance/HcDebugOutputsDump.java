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
 *     dot\t<label>\t<name>=<a dot file the run wrote in the working directory, escaped>
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
        final String variantsOut = "--debug-assembly-variants-out";
        withFile(dir, "hc-assembly-variants", hc, e2eBam, e2eFasta, variantsOut, List.of());
        withFile(dir, "hc-assembly-variants-chr2", hc, e2eBam, e2eFasta, variantsOut, List.of("-L", "chr2"));
        withFile(dir, "ramped-assembly-variants", ramped, e2eBam, e2eFasta, variantsOut, List.of("-L", "chr2"));
        withFile(dir, "mutect2-assembly-variants", m2, e2eBam, e2eFasta, variantsOut, List.of());
        withFile(dir, "mutect2-assembly-variants-chr2", m2, e2eBam, e2eFasta, variantsOut, List.of("-L", "chr2"));
        final String alm = "--alm-path";
        withFile(dir, "hc-alm", hc, e2eBam, e2eFasta, alm, List.of("--alm-interval", "chr2", "-L", "chr2"));
        withFile(dir, "hc-alm-window", hc, e2eBam, e2eFasta, alm, List.of("--alm-interval", "chr1:1000-1400"));
        withFile(dir, "hc-alm-no-interval", hc, e2eBam, e2eFasta, alm, List.of());
        withFile(dir, "hc-alm-pdhmm", hc, e2eBam, e2eFasta, alm,
                List.of("--alm-interval", "chr1:1100-1350", "--use-pdhmm", "-L", "chr1:1150-1300"));
        withFile(dir, "ramped-alm", ramped, e2eBam, e2eFasta, alm, List.of("--alm-interval", "chr2", "-L", "chr2"));
        withFile(dir, "mutect2-alm", m2, e2eBam, e2eFasta, alm, List.of("--alm-interval", "chr2", "-L", "chr2"));
        final String hmmResults = "--pair-hmm-results-file";
        withFile(dir, "hc-pairhmm-results", hc, e2eBam, e2eFasta, hmmResults, List.of("-L", "chr2:380-420"));
        withFile(dir, "ramped-pairhmm-results", ramped, e2eBam, e2eFasta, hmmResults, List.of("-L", "chr2:380-420"));
        withFile(dir, "mutect2-pairhmm-results", m2, e2eBam, e2eFasta, hmmResults, List.of("-L", "chr2:380-420"));
        final String state = "--debug-assembly-region-state";
        withFile(dir, "hc-region-state", hc, e2eBam, e2eFasta, state, List.of("-L", "chr2"));
        withFile(dir, "hc-region-state-filter", hc, e2eBam, e2eFasta, state,
                List.of("-L", "chr1:1100-1300", "--flow-filter-alleles"));
        withFile(dir, "ramped-region-state", ramped, e2eBam, e2eFasta, state, List.of("-L", "chr2"));
        final String genotyper = "--debug-genotyper-output";
        withFile(dir, "hc-genotyper", hc, e2eBam, e2eFasta, genotyper, List.of("-L", "chr2"));
        withFile(dir, "hc-genotyper-cis", hc, e2eBam, e2eFasta, genotyper, List.of("-L", "chr1:1150-1300"));
        withFile(dir, "hc-genotyper-dynamic", hc, e2eBam, e2eFasta, genotyper,
                List.of("-L", "chr1:900-1050", "--enable-dynamic-read-disqualification-for-genotyping"));
        withFile(dir, "hc-genotyper-legacy", hc, e2eBam, e2eFasta, genotyper,
                List.of("-L", "chr1:700-800", "--enable-legacy-assembly-region-trimming"));
        withFile(dir, "hc-genotyper-dragen", hc, e2eBam, e2eFasta, genotyper, List.of("-L", "chr1:900-1050", "--dragen-mode"));
        withFile(dir, "hc-genotyper-pdhmm", hc, e2eBam, e2eFasta, genotyper, List.of("-L", "chr1:1150-1300", "--use-pdhmm"));
        withFile(dir, "hc-genotyper-dragen-snp", hc, e2eBam, e2eFasta, genotyper,
                List.of("-L", "chr1:1150-1300", "--dragen-mode"));
        withFile(dir, "hc-genotyper-disqualify", hc, e2eBam, e2eFasta, genotyper,
                List.of("-L", "chr1:900-1050", "--expected-mismatch-rate-for-read-disqualification", "0.0"));
        withFile(dir, "hc-genotyper-disqualify-dynamic", hc, e2eBam, e2eFasta, genotyper,
                List.of("-L", "chr1:900-1050", "--expected-mismatch-rate-for-read-disqualification", "0.0",
                        "--enable-dynamic-read-disqualification-for-genotyping", "--dynamic-read-disqualification-threshold", "0.02"));
        withFile(dir, "ramped-genotyper", ramped, e2eBam, e2eFasta, genotyper, List.of("-L", "chr2"));
        withDots(dir, "hc-graph", hc, e2eBam, e2eFasta, List.of());
        withDots(dir, "hc-graph-debug", hc, e2eBam, e2eFasta,
                List.of("-L", "chr1:900-1300", "--debug-graph-transformations"));
        withDots(dir, "hc-graph-adaptive", hc, e2eBam, e2eFasta,
                List.of("-L", "chr1:900-1300", "--adaptive-pruning", "--debug-graph-transformations"));
        withDots(dir, "hc-graph-linked", hc, e2eBam, e2eFasta,
                List.of("-L", "chr1:900-1300", "--linked-de-bruijn-graph", "--assembly-region-padding", "30",
                        "--debug-graph-transformations"));
        withDots(dir, "hc-graph-linked-long-name", hc, e2eBam, e2eFasta,
                List.of("-L", "chr2", "--linked-de-bruijn-graph", "--debug-graph-transformations"));
        withDots(dir, "mutect2-graph", m2, e2eBam, e2eFasta,
                List.of("-L", "chr1:1150-1300", "--debug-graph-transformations"));
        withDots(dir, "ramped-graph", ramped, e2eBam, e2eFasta,
                List.of("-L", "chr2", "--debug-graph-transformations"));
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
    /**
     * --graph-output, and with --debug-graph-transformations every step's graph, which the
     * assembler writes in the working directory: each printed by name, in name order, then
     * deleted.
     */
    static void withDots(final Path dir, final String label,
                         final java.util.function.Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> tool,
                         final Path bam, final Path fasta, final List<String> extra) throws Exception {
        withFile(dir, label, tool, bam, fasta, "--graph-output", extra);
        final Path cwd = Path.of("").toAbsolutePath();
        final List<Path> dots;
        try (java.util.stream.Stream<Path> listing = Files.list(cwd)) {
            dots = listing.filter(p -> p.getFileName().toString().endsWith(".dot")).sorted().toList();
        }
        for (final Path dot : dots) {
            System.out.printf("dot\t%s\t%s=%s%n", label, dot.getFileName(), ReferenceQueryDump.escape(Files.readString(dot)));
            Files.delete(dot);
        }
    }

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
