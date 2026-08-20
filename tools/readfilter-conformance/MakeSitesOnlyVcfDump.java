/*
 * MakeSitesOnlyVcf, taken from the reference.
 *
 * A VCF with its genotype columns dropped, or kept for a named subset. The tool is a header rebuilt
 * with the samples the user asked for and a per-record subset, and everything else is the writer.
 *
 * Eight behaviours this is built to catch.
 *
 *   - SAMPLE IS A TreeSet, so the output's sample columns are in ALPHABETICAL order however the
 *     user typed them and whatever order the input had;
 *   - NO SAMPLE AT ALL IS THE DEFAULT and leaves the file with eight columns: no FORMAT column and
 *     no sample columns, not an empty FORMAT;
 *   - A NAME THE INPUT DOES NOT CARRY BECOMES A COLUMN OF ITS OWN. The header is built from the
 *     names the user asked for, not from the names the file has, so the output declares a sample
 *     that never existed and every record writes `./.` under it;
 *   - THE ANNOTATIONS ARE KEPT AS THEY WERE. The INFO fields of a subset record are the whole
 *     file's, not recomputed, so AC and AN can disagree with the genotypes that are left;
 *   - THE ALLELES ARE RESET FROM THE ORIGINAL RECORD, `builder.alleles(ctx.getAlleles())`, so an
 *     alternate no remaining genotype calls still appears in the ALT column;
 *   - THE GENOTYPES ARE NO LONGER LAZY once subset, so the FORMAT column is recomputed from the
 *     genotypes that remain and can be SHORTER than the input's;
 *   - A FILE WITH NO SAMPLES AT ALL passes through unchanged;
 *   - AND CREATE_INDEX IS ON BY DEFAULT for this tool, so a file whose header declares no contigs
 *     is a refusal rather than an unindexed output.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     sites\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MakeSitesOnlyVcfDump
 */

import picard.vcf.MakeSitesOnlyVcf;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MakeSitesOnlyVcfDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
            + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
            + "##contig=<ID=chr1,length=240>\n";

    static final String THREE_SAMPLES = HEADER
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tzeta\talpha\tmiddle\n"
            + "chr1\t10\trs1\tA\tC\t50\tPASS\tAC=2;AN=6\tGT:AD:DP\t0/1:10,5:15\t0/0:20,0:20\t1/1:0,9:9\n"
            + "chr1\t20\t.\tT\tG,C\t.\tLowQual\tAC=1,1;AN=6\tGT:AD\t0/1:5,5,0\t0/2:4,0,4\t0/0:9,0,0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("make-sites-only-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MakeSitesOnlyVcfDump: a VCF with its genotypes dropped");

        // The default: no samples at all.
        run("sites-only", dir, THREE_SAMPLES);
        // One sample, and the same one asked for twice.
        run("one-sample", dir, THREE_SAMPLES, "S=alpha");
        run("one-sample-twice", dir, THREE_SAMPLES, "S=alpha", "S=alpha");
        // Two samples given in the order the file does not have them, which the TreeSet sorts.
        run("two-samples-unsorted", dir, THREE_SAMPLES, "S=zeta", "S=alpha");
        // Every sample, which is the file with its columns reordered.
        run("all-samples", dir, THREE_SAMPLES, "S=zeta", "S=alpha", "S=middle");
        // A sample the file does not carry, alone and beside one it does.
        run("absent-sample", dir, THREE_SAMPLES, "S=absent");
        run("absent-and-present", dir, THREE_SAMPLES, "S=absent", "S=alpha");
        // A file that already has no samples.
        run("already-sites-only", dir,
                HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        + "chr1\t10\trs1\tA\tC\t50\tPASS\tAC=2;AN=6\n");
        // A header declaring no contigs, which the on-by-default index refuses.
        run("no-contigs", dir,
                HEADER.replace("##contig=<ID=chr1,length=240>\n", "")
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\talpha\n"
                        + "chr1\t10\trs1\tA\tC\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n");
        // The same file with the index turned off, which writes.
        run("no-contigs-no-index", dir,
                HEADER.replace("##contig=<ID=chr1,length=240>\n", "")
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\talpha\n"
                        + "chr1\t10\trs1\tA\tC\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n",
                "CREATE_INDEX=false");
    }

    static void run(final String label, final Path dir, final String input, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        final Path out = dir.resolve("sites-" + label + ".vcf");
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final List<String> argv = new ArrayList<>(Arrays.asList("I=" + in, "O=" + out));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new MakeSitesOnlyVcf().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("sites\t%s=%s%n", label, ReferenceQueryDump.escape(Files.readString(out)));
    }
}
