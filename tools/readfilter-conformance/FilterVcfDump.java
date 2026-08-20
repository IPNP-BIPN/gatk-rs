/*
 * FilterVcf, taken from the reference.
 *
 * Three site filters and two genotype filters applied to every record, written back out with the
 * FILTER column and every genotype's FT set.
 *
 * Nine behaviours this is built to catch.
 *
 *   - EVERY GENOTYPE IS GIVEN AN FT AND A PASSING ONE NEVER REACHES THE FILE. The iterator sets
 *     `FT=PASS` on every genotype it did not filter, and the writer drops it, so the FORMAT column
 *     gains `FT` only on the records where something was filtered;
 *   - THE SITE FILTERS ARE COLLECTED INTO A SET and the writer sorts it, so two filters come out
 *     alphabetically rather than in the order the filters ran;
 *   - `LowQD` IS NOT APPLIED WHEN QD IS ABSENT. The default is -1 and the test is `qd >= 0 && qd <
 *     minimum`, so a record with no QD passes however low the threshold;
 *   - `StrandBias` IS APPLIED WHEN FS IS ABSENT ONLY IF THE THRESHOLD IS NEGATIVE, the default
 *     being 0 and the test `fs > max`;
 *   - THE ALLELE BALANCE FILTER GROUPS BY THE GENOTYPE'S ALLELE LIST, so two samples with the same
 *     het call share one tally and a third with a different call has its own;
 *   - IT IGNORES A HET GENOTYPE WITH NO AD, and answers nothing at all when the record has no het
 *     genotype;
 *   - THE GENOTYPE FILTERS READ getGQ() AND getDP(), which are -1 when absent, so a genotype
 *     carrying neither is filtered by both at any non-negative threshold;
 *   - `AllGtsFiltered` IS A SITE FILTER SET BY THE GENOTYPES: a record every one of whose
 *     genotypes was filtered is itself filtered, and that filter replaces the PASS the site checks
 *     would otherwise have written;
 *   - AND THE HEADER THE WRITER GETS IS THE READER'S OWN OBJECT, mutated in place, so the four
 *     lines the tool adds are there even though `writeHeader` was handed the input's header.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     filtered\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FilterVcfDump
 */

import picard.vcf.filter.FilterVcf;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FilterVcfDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher strand\">\n"
            + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">\n"
            + "##contig=<ID=chr1,length=240>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\tthree\n";

    /** One record per behaviour, all three samples present. */
    static final String RECORDS =
            // Everything passes.
            "chr1\t10\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:AD:DP:GQ\t0/1:10,10:20:99\t"
            + "0/0:20,0:20:99\t1/1:0,20:20:99\n"
            // A het whose allele balance is far off, shared by two samples.
            + "chr1\t20\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:AD:DP:GQ\t0/1:19,1:20:99\t"
            + "0/1:18,2:20:99\t0/0:20,0:20:99\n"
            // High FS and low QD together, which is two site filters at once.
            + "chr1\t30\t.\tA\tC\t50\t.\tFS=99.0;QD=0.5\tGT:AD:DP:GQ\t0/1:10,10:20:99\t"
            + "0/0:20,0:20:99\t0/0:20,0:20:99\n"
            // No QD and no FS at all.
            + "chr1\t40\t.\tA\tC\t50\t.\t.\tGT:AD:DP:GQ\t0/1:10,10:20:99\t0/0:20,0:20:99\t"
            + "0/0:20,0:20:99\n"
            // One genotype with neither GQ nor DP, one with a low GQ, one with a low DP.
            + "chr1\t50\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:AD\t0/1:10,10\t"
            + "0/0:20,0\t1/1:0,20\n"
            + "chr1\t60\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:AD:DP:GQ\t0/1:10,10:20:5\t"
            + "0/0:20,0:2:99\t1/1:0,20:20:99\n"
            // A het with no AD, which the allele balance filter skips.
            + "chr1\t70\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:DP:GQ\t0/1:20:99\t0/0:20:99\t"
            + "0/0:20:99\n"
            // No het genotype at all.
            + "chr1\t80\t.\tA\tC\t50\t.\tFS=1.0;QD=20.0\tGT:AD:DP:GQ\t0/0:20,0:20:99\t"
            + "1/1:0,20:20:99\t0/0:20,0:20:99\n"
            // A record that already carries a filter, which is replaced rather than added to.
            + "chr1\t90\t.\tA\tC\t50\tmine\tFS=1.0;QD=20.0\tGT:AD:DP:GQ\t0/1:10,10:20:99\t"
            + "0/0:20,0:20:99\t0/0:20,0:20:99\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filter-vcf-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FilterVcfDump: a VCF with its sites and genotypes filtered");

        final String input = HEADER + RECORDS;
        run("defaults", dir, input);
        // Each threshold moved in turn.
        run("min-ab", dir, input, "MIN_AB=0.4");
        run("max-fs", dir, input, "MAX_FS=0.5");
        run("min-qd", dir, input, "MIN_QD=25");
        run("min-gq", dir, input, "MIN_GQ=50");
        run("min-dp", dir, input, "MIN_DP=10");
        // A negative FS threshold, which filters the record that has no FS at all.
        run("negative-fs", dir, input, "MAX_FS=-1");
        // Every threshold at once.
        run("everything", dir, input, "MIN_AB=0.45", "MAX_FS=0.5", "MIN_QD=25", "MIN_GQ=50",
                "MIN_DP=10");
        // A file with no records.
        run("no-records", dir, HEADER);
        // A header with no contigs, which a .vcf output refuses.
        run("no-contigs", dir,
                HEADER.replace("##contig=<ID=chr1,length=240>\n", "") + RECORDS);
    }

    static void run(final String label, final Path dir, final String input, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final Path out = dir.resolve("filtered-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "I=" + in, "O=" + out, "CREATE_INDEX=false"));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new FilterVcf().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("filtered\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
