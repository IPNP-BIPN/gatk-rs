/*
 * The callers force-calling, --alleles, end to end (G3.2l, #1346).
 *
 * A locus a given record overlaps is active whatever its pileup. A region's given events are the
 * records over its padded span, filtered ones dropped unless --force-call-filtered-alleles, split
 * one per alternate (right trimmed when split); HaplotypeCaller keeps those starting in the region,
 * the trimmer sees them, and each one the assembly lacks is inserted into the five best
 * haplotypes. An alternate consistent with a given event is genotyped and emitted whatever its
 * confidence. RampedHaplotypeCaller only genotypes them; Mutect2 inserts them and keeps them past
 * its tumour and normal odds.
 *
 * Read sets: HcEndToEndDump's and HcPileupDump's paired one, each with a VCF of given alleles made
 * here from the reference: alleles the reads carry, alleles they do not (a SNP, a deletion, an
 * insertion), a multiallelic record whose MNP alternate shares a suffix, a filtered record, and
 * one in a region with no reads. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     vcf\t<set>=<the given alleles' VCF, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcAllelesDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcAllelesDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-alleles-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcAllelesDump: the callers force-calling");

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

        final Path e2eAlleles = allelesVcf(e2eDir, "e2e", new String[]{"chr1", "chr2"}, references, new String[][]{
                {"chr1", "300", "0", "C", "."},        // carried by the reads
                {"chr1", "410", "0", "T", "."},        // carried by none
                {"chr1", "455", "3", "", "."},         // a deletion of three bases nobody carries
                {"chr1", "610", "0", "+GA", "."},      // an insertion nobody carries
                {"chr1", "700", "2", "MNP", "."},      // multiallelic: a SNP and an MNP sharing a suffix
                {"chr1", "850", "0", "G", "LowQual"},  // filtered
                {"chr1", "1215", "0", "A", "."},       // between the two cis SNPs
                {"chr2", "120", "0", "C", "."},
                {"chr2", "795", "0", "T", "PASS"}});
        final Path pairedAlleles = allelesVcf(pairedDir, "paired", new String[]{"chr1"}, new String[]{pairedReference}, new String[][]{
                {"chr1", "400", "0", "A", "."},
                {"chr1", "430", "0", "G", "."},
                {"chr1", "600", "0", "+TG", "."},
                {"chr1", "800", "0", "C", "."}});

        final String alleles = "--alleles";
        final String e2eGiven = e2eAlleles.toString();
        final String pairedGiven = pairedAlleles.toString();
        HcAllAnnotationsDump.run(dir, "hc-alleles", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(alleles, e2eGiven));
        HcAllAnnotationsDump.run(dir, "hc-alleles-filtered", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "--force-call-filtered-alleles"));
        HcAllAnnotationsDump.run(dir, "hc-alleles-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1:250-750"));
        HcAllAnnotationsDump.run(dir, "hc-alleles-no-assembly", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "--kmer-size", "101", "--dont-increase-kmer-sizes-for-cycles"));
        HcAllAnnotationsDump.run(dir, "hc-alleles-pileup-pdhmm", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "--pileup-detection", "--use-pdhmm"));
        HcAllAnnotationsDump.run(dir, "hc-alleles-intervals", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "-L", "chr1:400-470", "-L", "chr2:100-140"));
        HcAllAnnotationsDump.run(dir, "ramped-alleles", RampedHaplotypeCaller::new, e2eBam, e2eFasta, List.of(alleles, e2eGiven));
        HcAllAnnotationsDump.run(dir, "mutect2-alleles", Mutect2::new, e2eBam, e2eFasta, List.of(alleles, e2eGiven));
        HcAllAnnotationsDump.run(dir, "mutect2-alleles-filtered-tn", Mutect2::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "--force-call-filtered-alleles", "--normal-sample", "s2"));
        HcAllAnnotationsDump.run(dir, "mutect2-alleles-gvcf", Mutect2::new, e2eBam, e2eFasta,
                List.of(alleles, e2eGiven, "-ERC", "GVCF", "-L", "chr2"));
        HcAllAnnotationsDump.run(dir, "paired-hc-alleles", HaplotypeCaller::new, pairedBam, pairedFasta, List.of(alleles, pairedGiven));
        HcAllAnnotationsDump.run(dir, "paired-mutect2-alleles", Mutect2::new, pairedBam, pairedFasta, List.of(alleles, pairedGiven));
    }

    /**
     * A VCF of given alleles, from {contig, position, deleted length, alt, filter}: alt "+XY"
     * inserts after the reference base, "" with a length deletes, "MNP" makes a multiallelic
     * record of a SNP and a two-base MNP sharing its second base.
     */
    static Path allelesVcf(final Path dir, final String set, final String[] contigs, final String[] references,
                           final String[][] records) throws Exception {
        final StringBuilder vcf = new StringBuilder("##fileformat=VCFv4.2\n");
        vcf.append("##FILTER=<ID=LowQual,Description=\"Low quality\">\n");
        for (int c = 0; c < contigs.length; c++) {
            vcf.append("##contig=<ID=").append(contigs[c]).append(",length=").append(references[c].length()).append(">\n");
        }
        vcf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final String[] r : records) {
            final String ref = references[java.util.Arrays.asList(contigs).indexOf(r[0])];
            final int pos = Integer.parseInt(r[1]);
            final int deleted = Integer.parseInt(r[2]);
            final String base = ref.substring(pos - 1, pos);
            final String refAllele;
            final String alt;
            if (r[3].equals("MNP")) {
                refAllele = ref.substring(pos - 1, pos + 1);
                final char other = ref.charAt(pos - 1) == 'A' ? 'G' : 'A';
                alt = other + refAllele.substring(1) + "," + (ref.charAt(pos - 1) == 'C' ? 'T' : 'C') + (ref.charAt(pos) == 'G' ? "T" : "G");
            } else if (r[3].startsWith("+")) {
                refAllele = base;
                alt = base + r[3].substring(1);
            } else if (deleted > 0) {
                refAllele = ref.substring(pos - 1, pos + deleted);
                alt = base;
            } else {
                refAllele = base;
                alt = r[3].equals(base) ? (base.equals("A") ? "C" : "A") : r[3];
            }
            vcf.append(r[0]).append('\t').append(pos).append("\t.\t").append(refAllele).append('\t').append(alt)
                    .append("\t.\t").append(r[4]).append("\t.\n");
        }
        final Path path = dir.resolve("given.vcf");
        Files.writeString(path, vcf.toString());
        final htsjdk.tribble.index.Index index = htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                path, new htsjdk.variant.vcf.VCFCodec());
        index.write(dir.resolve("given.vcf.idx"));
        System.out.printf("vcf\t%s=%s%n", set, ReferenceQueryDump.escape(vcf.toString()));
        return path;
    }

}
