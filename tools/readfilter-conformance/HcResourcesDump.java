/*
 * HaplotypeCaller's variant resources, end to end (G3.2o, #1349).
 *
 * --dbsnp fills the ID column of a call whose event a dbSNP record shares and sets DB;
 * --comparison sets a flag named by its tag (or path); a comparison that is the dbSNP file is
 * renamed DB and one tagged DB is refused. --population-callset turns the reference-confidence
 * model's posteriors on: every block and every call it emits gains PP, its call and GQ from the
 * posteriors, and a call PG and recomputed AC, AF and AN, the priors coming from the callset's
 * records that start there, or from --num-reference-samples-if-no-call.
 *
 * Read set: HcEndToEndDump's. The resources are written from a first plain run's calls: a dbSNP
 * file with every other call under an rsID and the rest with another alternate, a comparison with
 * every third call, and a population callset with counts at every call and at a few reference
 * positions. Each is indexed. Runs in ONE JVM, in this order.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     vcf\t<name>=<a resource VCF, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcResourcesDump
 */
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcResourcesDump {


    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-resources-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcResourcesDump: the variant resources of HaplotypeCaller");

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
        HcAllAnnotationsDump.run(dir, "hc-plain", hc, e2eBam, e2eFasta, List.of());
        final List<String[]> calls = new ArrayList<>();
        for (final String line : Files.readAllLines(dir.resolve("out-hc-plain.vcf"))) {
            if (!line.startsWith("#")) {
                calls.add(line.split("\t"));
            }
        }
        final String header = "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                + "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele count\">\n"
                + "##contig=<ID=chr1,length=" + references[0].length() + ">\n"
                + "##contig=<ID=chr2,length=" + references[1].length() + ">\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        final StringBuilder dbsnp = new StringBuilder(header);
        final StringBuilder comp = new StringBuilder(header);
        final List<String[]> populationRecords = new ArrayList<>();
        for (int n = 0; n < calls.size(); n++) {
            final String[] c = calls.get(n);
            final String alt = c[4].split(",")[0];
            final String otherAlt = "ACGT".replace(c[3].substring(0, 1), "").substring(0, 1);
            dbsnp.append(String.join("\t", c[0], c[1], "rs" + (1000 + n), c[3],
                    n % 2 == 0 ? alt : otherAlt, ".", ".", ".")).append('\n');
            if (n % 3 == 0) {
                comp.append(String.join("\t", c[0], c[1], ".", c[3], alt, ".", ".", ".")).append('\n');
            }
            populationRecords.add(new String[]{c[0], c[1], ".", c[3], alt, ".", ".",
                    "AC=" + (n % 4 + 1) + ";AN=40;MLEAC=" + (n % 5 + 1)});
        }
        for (final int position : new int[]{300, 310, 320}) {
            final String ref = references[0].substring(position - 1, position);
            final String alt = "ACGT".replace(ref, "").substring(0, 1);
            populationRecords.add(new String[]{"chr1", String.valueOf(position), ".", ref, alt, ".", ".",
                    "AC=10;AN=20;MLEAC=12"});
        }
        populationRecords.sort(java.util.Comparator.<String[], String>comparing(r -> r[0])
                .thenComparingInt(r -> Integer.parseInt(r[1])));
        final StringBuilder population = new StringBuilder(header);
        for (final String[] r : populationRecords) {
            population.append(String.join("\t", r)).append('\n');
        }
        final String db = resource(dir, "dbsnp", dbsnp.toString()).toString();
        final String co = resource(dir, "comp", comp.toString()).toString();
        final String po = resource(dir, "population", population.toString()).toString();

        HcAllAnnotationsDump.run(dir, "hc-dbsnp", hc, e2eBam, e2eFasta, List.of("--dbsnp", db));
        HcAllAnnotationsDump.run(dir, "hc-dbsnp-comp", hc, e2eBam, e2eFasta,
                List.of("--dbsnp", db, "--comparison", co, "--comparison:truth", db));
        HcAllAnnotationsDump.run(dir, "hc-comp-is-dbsnp", hc, e2eBam, e2eFasta,
                List.of("--dbsnp", db, "--comparison", db));
        HcAllAnnotationsDump.run(dir, "hc-comp-tagged-db", hc, e2eBam, e2eFasta,
                List.of("--comparison:DB", co));
        HcAllAnnotationsDump.run(dir, "hc-dbsnp-gvcf", hc, e2eBam, e2eFasta,
                List.of("--dbsnp", db, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        HcAllAnnotationsDump.run(dir, "hc-population-vcf", hc, e2eBam, e2eFasta, List.of("--population-callset", po));
        HcAllAnnotationsDump.run(dir, "hc-population-gvcf", hc, e2eBam, e2eFasta,
                List.of("--population-callset", po, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        HcAllAnnotationsDump.run(dir, "hc-population-gvcf-numref", hc, e2eBam, e2eFasta,
                List.of("--population-callset", po, "--num-reference-samples-if-no-call", "5",
                        "-ERC", "GVCF", "--sample-name", "s2", "-L", "chr1"));
        HcAllAnnotationsDump.run(dir, "hc-population-bp", hc, e2eBam, e2eFasta,
                List.of("--population-callset", po, "-ERC", "BP_RESOLUTION", "--sample-name", "s1", "-L", "chr1:250-450"));
        HcAllAnnotationsDump.run(dir, "ramped-resources", ramped, e2eBam, e2eFasta,
                List.of("--dbsnp", db, "--comparison:truth", co, "--population-callset", po,
                        "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
    }

    /** A resource VCF, printed and indexed beside itself. */
    static Path resource(final Path dir, final String name, final String text) throws Exception {
        final Path vcf = dir.resolve(name + ".vcf");
        Files.writeString(vcf, text);
        System.out.printf("vcf\t%s=%s%n", name, ReferenceQueryDump.escape(text));
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[]{"-I", vcf.toString()});
        return vcf;
    }
}
