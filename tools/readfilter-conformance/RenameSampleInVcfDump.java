/*
 * RenameSampleInVcf, taken from the reference.
 *
 * A single-sample VCF written back out with its one sample column renamed. Small, and every one of
 * its corners is in what the WRITER does rather than in what the tool does.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE HEADER IS REBUILT FROM ITS METADATA IN INPUT ORDER and one sample name, so a line the
 *     input carried survives and its ORDER is the writer's, not the input's;
 *   - THE SAMPLE NAME IS NOT VALIDATED. A name with a space, a tab-free oddity or an empty string
 *     goes straight into the column header;
 *   - A MULTI-SAMPLE INPUT IS REFUSED, but a SITES-ONLY input is NOT: a VCF with no sample column
 *     at all passes the size test and then the rename gives it one, with every record carrying a
 *     missing genotype;
 *   - OLD_SAMPLE_NAME IS CHECKED AGAINST THE FIRST SAMPLE only, and its refusal names the sample
 *     that was there;
 *   - THE RECORDS ARE WRITTEN BACK THROUGH THE PARSER, so a record's INFO field order, its float
 *     formatting and its missing-value spellings are the writer's rather than the input's;
 *   - A GENOTYPE'S OWN FIELDS SURVIVE THE RENAME, the sample being renamed in the header and the
 *     genotype carried through by position;
 *   - THE FILE'S OTHER HEADER LINES ARE NOT TOUCHED, contigs and INFO declarations included, so an
 *     input declaring nothing still writes records that use undeclared keys;
 *   - AND THE QUAL COLUMN IS REFORMATTED: `50.00` comes back `50` and `10.5` comes back `10.50`,
 *     because the writer prints an integral quality as an integer and everything else to two
 *     places. The output carries no provenance line of its own either: nothing says it was
 *     renamed.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     renamed\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: RenameSampleInVcfDump
 */

import picard.vcf.RenameSampleInVcf;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class RenameSampleInVcfDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
            + "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">\n"
            + "##contig=<ID=chr1,length=240>\n";

    static final String RECORDS =
            "chr1\t10\trs1\tA\tC\t50.00\tPASS\tAF=0.5;DB\tGT:AD:DP:GQ\t0/1:10,5:15:99\n"
            + "chr1\t20\t.\tT\tG,C\t.\tLowQual\tAF=0.25,0.125\tGT:AD:DP\t1/2:1,2,3:6\n"
            + "chr1\t30\t.\tG\t.\t10.5\t.\t.\tGT\t./.\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("rename-sample-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# RenameSampleInVcfDump: a single-sample VCF with its column renamed");

        final String single = HEADER + column("NA12878") + RECORDS;
        run("plain", dir, single, "NEW_SAMPLE_NAME=renamed");
        // The old name asserted, and asserted wrongly.
        run("old-name-right", dir, single, "NEW_SAMPLE_NAME=renamed", "OLD_SAMPLE_NAME=NA12878");
        run("old-name-wrong", dir, single, "NEW_SAMPLE_NAME=renamed", "OLD_SAMPLE_NAME=OTHER");
        // Names the tool does not validate.
        run("name-with-space", dir, single, "NEW_SAMPLE_NAME=two words");
        run("name-that-is-a-number", dir, single, "NEW_SAMPLE_NAME=12345");
        run("same-name", dir, single, "NEW_SAMPLE_NAME=NA12878");
        // A sites-only VCF, which has no sample to rename and gets one anyway.
        run("sites-only", dir,
                HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        + "chr1\t10\trs1\tA\tC\t50.00\tPASS\tAF=0.5;DB\n",
                "NEW_SAMPLE_NAME=renamed");
        // Two samples, which is the refusal.
        run("two-samples", dir,
                HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n"
                        + "chr1\t10\trs1\tA\tC\t50.00\tPASS\tAF=0.5\tGT\t0/1\t1/1\n",
                "NEW_SAMPLE_NAME=renamed");
        // A file declaring no INFO lines at all, whose records still use them.
        run("undeclared-info", dir,
                "##fileformat=VCFv4.2\n"
                        + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                        + "##contig=<ID=chr1,length=240>\n"
                        + column("NA12878")
                        + "chr1\t10\t.\tA\tC\t50.00\t.\tXX=1\tGT\t0/1\n",
                "NEW_SAMPLE_NAME=renamed");
        // No records at all.
        run("no-records", dir, HEADER + column("NA12878"), "NEW_SAMPLE_NAME=renamed");
    }

    static String column(final String sample) {
        return "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n";
    }

    static void run(final String label, final Path dir, final String input, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        final Path out = dir.resolve("renamed-" + label + ".vcf");
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "I=" + in, "O=" + out, "CREATE_INDEX=false"));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new RenameSampleInVcf().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("renamed\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
