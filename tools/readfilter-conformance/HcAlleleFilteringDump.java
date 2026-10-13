/*
 * --flow-filter-alleles in the callers, end to end (G3.2f, #1335).
 *
 * AlleleFiltering drops haplotypes before genotyping: alleles no haplotype carries together are
 * linked when three bases apart or less, or twenty apart and equal up to one hmer's length once
 * inserted in the reference; within each linked group the weakest allele by quality against its
 * inverse (the het or hom-alt PL over the hom-inverse one, the smallest over samples), or the most
 * strand-biased, has its haplotypes disabled, until none qualifies. A lone allele goes only on a
 * quality of one or less, unless --flow-filter-lone-alleles. The haplotypes kept come out in the
 * order of a HashSet over their bases, whether or not any went, and every allele weighed where a
 * candidate was found marks its start SUSP_NOISY_ADJACENT_TP_VARIANT.
 *
 * Over HcEndToEndDump's reads: the defaults, a quality threshold of 300 and of 2000, a strand odds
 * ratio threshold of 0.5, lone alleles removed, RampedHaplotypeCaller, and GVCF mode on one
 * sample. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcAlleleFilteringDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Random;

public class HcAlleleFilteringDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-allele-filtering-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcAlleleFilteringDump: --flow-filter-alleles in the callers");

        final Random random = new Random(20261010L);
        final String[] references = new String[HcEndToEndDump.CONTIGS.length];
        for (int c = 0; c < HcEndToEndDump.CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < HcEndToEndDump.LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\t%s=%s%n", HcEndToEndDump.CONTIGS[c], references[c]);
        }
        final Path fasta = HcEndToEndDump.writeReference(dir, references);
        final Path bam = HcLinkedGraphDump.bam(dir, "input", HcEndToEndDump.sam(references, random));

        final String filter = "--flow-filter-alleles";
        HcAllAnnotationsDump.run(dir, "hc-filter", HaplotypeCaller::new, bam, fasta, List.of(filter));
        HcAllAnnotationsDump.run(dir, "hc-filter-qual-300", HaplotypeCaller::new, bam, fasta,
                List.of(filter, "--flow-filter-alleles-qual-threshold", "300"));
        HcAllAnnotationsDump.run(dir, "hc-filter-qual-2000", HaplotypeCaller::new, bam, fasta,
                List.of(filter, "--flow-filter-alleles-qual-threshold", "2000"));
        HcAllAnnotationsDump.run(dir, "hc-filter-sor", HaplotypeCaller::new, bam, fasta,
                List.of(filter, "--flow-filter-alleles-sor-threshold", "0.5"));
        HcAllAnnotationsDump.run(dir, "hc-filter-lone", HaplotypeCaller::new, bam, fasta,
                List.of(filter, "--flow-filter-lone-alleles", "--flow-filter-alleles-qual-threshold", "300"));
        HcAllAnnotationsDump.run(dir, "ramped-filter", RampedHaplotypeCaller::new, bam, fasta,
                List.of(filter, "--flow-filter-alleles-qual-threshold", "300"));
        HcAllAnnotationsDump.run(dir, "hc-filter-gvcf", HaplotypeCaller::new, bam, fasta,
                List.of(filter, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
    }
}
