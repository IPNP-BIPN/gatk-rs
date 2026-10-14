/*
 * --ploidy-regions in HaplotypeCaller and RampedHaplotypeCaller, end to end (G3.2p, #1350).
 *
 * Each BED feature is a SimpleCount whose name is the ploidy; every interval must lie on the
 * reads' dictionary. A locus takes the largest ploidy of the intervals overlapping it for its
 * activity engine (at least two), and a region the largest overlapping its span for its genotyping
 * engine and its reference-confidence model; anything else keeps --sample-ploidy.
 *
 * Read set: HcEndToEndDump's. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     bed\t<name>=<a ploidy BED, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcPloidyDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcPloidyDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-ploidy-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcPloidyDump: --ploidy-regions in the callers");

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
        final String haploid = bed(dir, "haploid", "chr1\t250\t800\t1\n");
        final String triploid = bed(dir, "triploid", "chr1\t900\t1300\t3\nchr2\t0\t500\t4\n");
        final String overlapping = bed(dir, "overlapping", "chr1\t250\t600\t1\nchr1\t500\t800\t3\n");
        final String zero = bed(dir, "zero", "chr2\t0\t500\t0\n");
        final String partial = bed(dir, "partial", "chr1\t990\t1000\t4\n");
        final String offContig = bed(dir, "off-contig", "chr9\t0\t100\t2\n");
        final String pastEnd = bed(dir, "past-end", "chr1\t0\t999999\t2\n");
        final String notANumber = bed(dir, "not-a-number", "chr1\t0\t100\tdiploid\n");
        final String negative = bed(dir, "negative", "chr1\t0\t100\t-1\n");
        final String regions = "--ploidy-regions";
        HcAllAnnotationsDump.run(dir, "hc-haploid", hc, e2eBam, e2eFasta, List.of(regions, haploid));
        HcAllAnnotationsDump.run(dir, "hc-triploid", hc, e2eBam, e2eFasta, List.of(regions, triploid));
        HcAllAnnotationsDump.run(dir, "hc-overlapping", hc, e2eBam, e2eFasta, List.of(regions, overlapping));
        HcAllAnnotationsDump.run(dir, "hc-zero", hc, e2eBam, e2eFasta, List.of(regions, zero, "-L", "chr2"));
        HcAllAnnotationsDump.run(dir, "hc-partial", hc, e2eBam, e2eFasta, List.of(regions, partial, "-L", "chr1:900-1300"));
        HcAllAnnotationsDump.run(dir, "hc-haploid-gvcf", hc, e2eBam, e2eFasta,
                List.of(regions, haploid, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1:200-900"));
        HcAllAnnotationsDump.run(dir, "hc-haploid-bp", hc, e2eBam, e2eFasta,
                List.of(regions, haploid, "-ERC", "BP_RESOLUTION", "--sample-name", "s1", "-L", "chr1:200-900"));
        HcAllAnnotationsDump.run(dir, "hc-triploid-bp", hc, e2eBam, e2eFasta,
                List.of(regions, triploid, "-ERC", "BP_RESOLUTION", "--sample-name", "s2", "-L", "chr1:950-1050"));
        HcAllAnnotationsDump.run(dir, "hc-triploid-sample-ploidy", hc, e2eBam, e2eFasta,
                List.of(regions, triploid, "--sample-ploidy", "1"));
        HcAllAnnotationsDump.run(dir, "hc-off-contig", hc, e2eBam, e2eFasta, List.of(regions, offContig));
        HcAllAnnotationsDump.run(dir, "hc-past-end", hc, e2eBam, e2eFasta, List.of(regions, pastEnd));
        HcAllAnnotationsDump.run(dir, "hc-not-a-number", hc, e2eBam, e2eFasta, List.of(regions, notANumber));
        HcAllAnnotationsDump.run(dir, "hc-negative", hc, e2eBam, e2eFasta, List.of(regions, negative));
        HcAllAnnotationsDump.run(dir, "ramped-triploid", ramped, e2eBam, e2eFasta, List.of(regions, triploid));
        HcAllAnnotationsDump.run(dir, "ramped-haploid-gvcf", ramped, e2eBam, e2eFasta,
                List.of(regions, haploid, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1:200-900"));
    }

    /** A ploidy BED, printed. */
    static String bed(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name + ".bed");
        Files.writeString(path, text);
        System.out.printf("bed\t%s=%s%n", name, ReferenceQueryDump.escape(text));
        return path.toString();
    }
}
