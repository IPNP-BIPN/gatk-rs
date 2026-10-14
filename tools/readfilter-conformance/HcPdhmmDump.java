/*
 * The callers with the partially determined PairHMM, end to end (G3.2g, #1338).
 *
 * --use-pdhmm replaces the assembled haplotypes, once trimmed, by partially determined ones:
 * every assembled event (less the bad pileup events, with the good ones) is grouped with the
 * events it cannot share a haplotype with, overlapping or realigning into another, and each locus
 * in turn is determined, at each of its alleles and the reference, while every other event is left
 * undetermined in the haplotype's alternate bytes. LoglessPDPairHMM scores them, a base flagged as
 * a possible SNP matching either base and a possible deletion skipped or walked; reads are not
 * realigned; a haplotype supports alleles only at its determined position. Also: the overlap
 * optimisation, determined haplotypes in place of PD ones, the GGA fallback, the results file,
 * and --dragen-378-concordance-mode, which sets all of it with pileup detection.
 *
 * The read sets are HcPileupDump's. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     file\t<label>=<the PDHMM results file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcPdhmmDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcPdhmmDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-pdhmm-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcPdhmmDump: the callers with the partially determined PairHMM");

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

        // HcDragenEndToEndDump's contig and reads.
        final Random strRandom = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(strRandom.nextInt(4)));
        }
        HcDragenEndToEndDump.splice(ref, 401, "AAAAAAAAAA");
        HcDragenEndToEndDump.splice(ref, 701, "CACACACACACACA");
        HcDragenEndToEndDump.splice(ref, 1001, "TTTTTTTT");
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String strReference = ref.toString();
        System.out.printf("fasta\tstr:%s=%s%n", HcDragenEndToEndDump.CONTIG, strReference);
        final Path strDir = dir.resolve("str");
        Files.createDirectories(strDir);
        final Path strFasta = HcDragenEndToEndDump.writeReference(strDir, strReference);
        final Path strBam = HcLinkedGraphDump.bam(strDir, "str", HcDragenEndToEndDump.sam(strReference));

        // The paired set.
        final Random pairedRandom = new Random(20261014L);
        final StringBuilder paired = new StringBuilder();
        for (int i = 0; i < HcPileupDump.PAIRED_LENGTH; i++) {
            paired.append("ACGT".charAt(pairedRandom.nextInt(4)));
        }
        final String pairedReference = paired.toString();
        System.out.printf("fasta\tpaired:chr1=%s%n", pairedReference);
        final Path pairedDir = dir.resolve("paired");
        Files.createDirectories(pairedDir);
        final Path pairedFasta = HcPileupDump.writeReference(pairedDir, pairedReference);
        final Path pairedBam = HcLinkedGraphDump.bam(pairedDir, "paired", HcPileupDump.pairedSam(pairedReference, pairedRandom));

        final String on = "--pileup-detection";
        final String indels = "--pileup-detection-enable-indel-pileup-calling";
        final String pdhmm = "--use-pdhmm";
        final List<String> noAssembly = List.of("--kmer-size", "101", "--dont-increase-kmer-sizes-for-cycles");
        final String dragen = "--dragen-378-concordance-mode";

        HcAllAnnotationsDump.run(dir, "hc-pdhmm", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(pdhmm));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-pileup", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(pdhmm, on, indels));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-overlap", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(pdhmm, on, indels, "--use-pdhmm-overlap-optimization"));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-determined", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(pdhmm, on, indels, "--make-determined-haps-from-pd-code"));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-no-assembly", HaplotypeCaller::new, e2eBam, e2eFasta,
                HcPileupDump.concat(List.of(pdhmm, on, indels), noAssembly));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-no-fallback", HaplotypeCaller::new, e2eBam, e2eFasta,
                HcPileupDump.concat(List.of(pdhmm, on, indels, "--fallback-gga-if-pdhmm-fails", "false"), noAssembly));
        HcAllAnnotationsDump.run(dir, "hc-pdhmm-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(pdhmm, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"));
        HcAllAnnotationsDump.run(dir, "hc-dragen-378", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(dragen));
        HcAllAnnotationsDump.run(dir, "hc-dragen-378-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(dragen, "-ERC", "GVCF", "--sample-name", "s2", "-L", "chr1:1150-1300"));
        HcAllAnnotationsDump.run(dir, "hc-dragen-378-with-3412", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(dragen, "--dragen-mode"));
        HcAllAnnotationsDump.run(dir, "ramped-pdhmm", RampedHaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(pdhmm, on, indels));
        HcAllAnnotationsDump.run(dir, "mutect2-pdhmm", Mutect2::new, e2eBam, e2eFasta, List.of(pdhmm, on, indels));
        withResults(dir, "hc-pdhmm-results", e2eBam, e2eFasta, List.of(pdhmm, on, indels, "-L", "chr2:360-440"));
        withResults(dir, "hc-pdhmm-overlap-results", e2eBam, e2eFasta,
                List.of(pdhmm, on, indels, "--use-pdhmm-overlap-optimization", "-L", "chr1:1190-1240"));
        HcAllAnnotationsDump.run(dir, "str-hc-pdhmm", HaplotypeCaller::new, strBam, strFasta, List.of(pdhmm, on, indels));
        HcAllAnnotationsDump.run(dir, "str-hc-dragen-378", HaplotypeCaller::new, strBam, strFasta, List.of(dragen));
        HcAllAnnotationsDump.run(dir, "str-hc-pdhmm-determined", HaplotypeCaller::new, strBam, strFasta,
                List.of(pdhmm, on, indels, "--make-determined-haps-from-pd-code"));
        HcAllAnnotationsDump.run(dir, "paired-hc-pdhmm", HaplotypeCaller::new, pairedBam, pairedFasta,
                List.of(pdhmm, on, indels));
        HcAllAnnotationsDump.run(dir, "paired-hc-dragen-378", HaplotypeCaller::new, pairedBam, pairedFasta,
                List.of(dragen));
        HcAllAnnotationsDump.run(dir, "paired-hc-dragen-378-no-assembly", HaplotypeCaller::new, pairedBam, pairedFasta,
                HcPileupDump.concat(List.of(dragen), noAssembly));
    }

    /** A HaplotypeCaller run, then its PDHMM results file. */
    static void withResults(final Path dir, final String label, final Path bam, final Path fasta,
                            final List<String> extra) throws Exception {
        final Path results = dir.resolve("results-" + label + ".txt");
        HcAllAnnotationsDump.run(dir, label, HaplotypeCaller::new, bam, fasta,
                HcPileupDump.concat(extra, List.of("--pdhmm-results-file", results.toString())));
        System.out.printf("file\t%s=%s%n", label, ReferenceQueryDump.escape(Files.readString(results)));
    }

}
