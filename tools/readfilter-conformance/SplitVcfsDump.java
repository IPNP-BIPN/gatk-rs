/*
 * SplitVcfs, taken from the reference.
 *
 * One VCF in, two out: the indels in one file and the SNPs in the other, and everything that is
 * neither in neither.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE INDEL TEST RUNS FIRST, so a record that is both by htsjdk's reckoning lands in the indel
 *     file and never reaches the SNP one;
 *   - A MIXED RECORD IS IN NEITHER FILE. `isIndel()` and `isSNP()` are both false for a record
 *     carrying a SNP and an indel alternate, and the tool silently counts it out;
 *   - A MONOMORPHIC RECORD IS ALSO IN NEITHER, its type being NO_VARIATION;
 *   - A SYMBOLIC ALTERNATE IS IN NEITHER, while a SPANNING-DELETION `*` ALTERNATE BESIDE A SNP
 *     LEAVES THE RECORD A SNP and it goes to the SNP file: `*` is one base long, so the type test
 *     never sees a length difference;
 *   - AN MNP IS IN NEITHER, being neither a SNP nor an indel;
 *   - STRICT IS ON BY DEFAULT, so a record that is neither is a refusal naming the type it found
 *     unless the user turns it off; the refusal comes after the earlier records are written;
 *   - BOTH FILES GET THE INPUT'S OWN HEADER, samples and all, however few records they end up
 *     holding;
 *   - AND CREATE_INDEX IS ON BY DEFAULT, so an input whose header declares no contigs is a refusal
 *     before either file is opened.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     snps\t<label>=<the whole SNP output, escaped>
 *     indels\t<label>=<the whole indel output, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SplitVcfsDump
 */

import picard.vcf.SplitVcfs;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SplitVcfsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
            + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
            + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            + "##contig=<ID=chr1,length=240>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n";

    /** One record of each type the tool can meet. */
    static final String EVERY_TYPE = HEADER
            // A plain SNP.
            + "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"
            // A biallelic SNP with two alternates, still a SNP.
            + "chr1\t20\t.\tA\tC,G\t50\tPASS\tAC=1,1\tGT\t1/2\n"
            // An insertion and a deletion.
            + "chr1\t30\t.\tA\tACGT\t50\tPASS\tAC=1\tGT\t0/1\n"
            + "chr1\t40\t.\tACGT\tA\t50\tPASS\tAC=1\tGT\t0/1\n"
            // A SNP and an indel in one record, which is MIXED.
            + "chr1\t50\t.\tA\tC,ACGT\t50\tPASS\tAC=1,1\tGT\t1/2\n"
            // An MNP.
            + "chr1\t60\t.\tAC\tGT\t50\tPASS\tAC=1\tGT\t0/1\n"
            // Monomorphic: no alternate at all.
            + "chr1\t70\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\n"
            // A symbolic alternate.
            + "chr1\t80\t.\tA\t<DEL>\t50\tPASS\tEND=90\tGT\t0/1\n"
            // A spanning deletion alternate beside a SNP.
            + "chr1\t100\t.\tA\tC,*\t50\tPASS\tAC=1,1\tGT\t1/2\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("split-vcfs-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SplitVcfsDump: the SNPs in one file and the indels in the other");

        // STRICT defaults to TRUE, so the file of every type refuses unless it is turned off.
        run("every-type", dir, EVERY_TYPE, "STRICT=false");
        run("strict", dir, EVERY_TYPE);
        // Only SNPs, so the indel file is a header and nothing else.
        run("snps-only", dir, HEADER
                + "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"
                + "chr1\t20\t.\tG\tT\t50\tPASS\tAC=1\tGT\t0/1\n");
        // Only indels.
        run("indels-only", dir, HEADER
                + "chr1\t30\t.\tA\tACGT\t50\tPASS\tAC=1\tGT\t0/1\n");
        // No records at all: two headers.
        run("no-records", dir, HEADER);
        // Strict over a file that is all SNPs, which does not refuse.
        run("strict-all-snps", dir, HEADER
                + "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n");
        // A header with no contigs, which the on-by-default index refuses.
        run("no-contigs", dir,
                HEADER.replace("##contig=<ID=chr1,length=240>\n", "")
                        + "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n");
        // The same with the index off.
        run("no-contigs-no-index", dir,
                HEADER.replace("##contig=<ID=chr1,length=240>\n", "")
                        + "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n",
                "CREATE_INDEX=false");
    }

    static void run(final String label, final Path dir, final String input, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        final Path snps = dir.resolve("snps-" + label + ".vcf");
        final Path indels = dir.resolve("indels-" + label + ".vcf");
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "I=" + in, "SNP_OUTPUT=" + snps, "INDEL_OUTPUT=" + indels));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new SplitVcfs().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("snps\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(snps)));
        System.out.printf("indels\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(indels)));
    }
}
