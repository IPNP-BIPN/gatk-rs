/*
 * CalculateAverageCombinedAnnotations' output, taken from the reference.
 *
 * Annotations that GenomicsDB summed across samples divided by the number of het and hom-var
 * samples. Nine lines of arithmetic, and most of the behaviour is in what it does NOT do.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE DIVISOR IS THE SECOND AND THIRD FIELDS OF `RAW_GENOTYPE_COUNT`, het plus hom-var, read
 *     as DOUBLES from a string list, so a count of `0,2,1` divides by 3;
 *   - A RECORD WHOSE DIVISOR IS ZERO IS WRITTEN THROUGH UNTOUCHED, gaining no AVERAGE_ field at
 *     all rather than one holding zero or a missing value;
 *   - A RECORD MISSING `RAW_GENOTYPE_COUNT` ENTIRELY IS A UserException that names the site, and
 *     it comes AFTER the header and any earlier records have been written;
 *   - AN ANNOTATION THE RECORD DOES NOT CARRY IS SKIPPED for that record, so two records can come
 *     out with different sets of AVERAGE_ fields;
 *   - THE ORIGINAL ANNOTATION IS KEPT, the average being added beside it, and THE AVERAGE IS
 *     WRITTEN BY THE ENCODER'S OWN FORMAT: `30.0 / 3` comes out `3.00`, two decimals, while the
 *     source annotation keeps the single decimal it was written with;
 *   - THE HEADER GAINS ONE INFO LINE PER REQUESTED ANNOTATION whether any record carries it or
 *     not, Number=1 and Type=Float, with a description quoting the source annotation twice;
 *   - AND THE LINE IS ADDED EVEN FOR AN ANNOTATION THE INPUT NEVER DECLARED, so the output can
 *     declare `AVERAGE_XX` while `XX` itself is undeclared;
 *   - AN EMPTY --summed-annotation-to-divide IS REFUSED before the header is written, but the
 *     argument is mandatory, so the parser refuses first and the tool's own message is
 *     unreachable.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     averaged\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CalculateAverageCombinedAnnotationsDump
 */

import org.broadinstitute.hellbender.tools.CalculateAverageCombinedAnnotations;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CalculateAverageCombinedAnnotationsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Allele specific QD\">\n"
            + "##INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=\"Counts of genotypes\">\n"
            + "##INFO=<ID=SUMMED,Number=1,Type=Float,Description=\"A summed annotation\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n";

    static String record(final int position, final String info) {
        return "chr1\t" + position + "\t.\tA\tC\t50\t.\t" + info + "\tGT\t0/1\t0/0\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("calculate-average-annotations-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CalculateAverageCombinedAnnotationsDump: summed annotations averaged");

        final String records =
                // Three het-or-homvar samples: 0 hom-ref, 2 het, 1 hom-var.
                record(100, "RAW_GT_COUNT=0,2,1;SUMMED=30.0;AS_QD=9.0")
                // A divisor of zero, which writes the record through untouched.
                + record(200, "RAW_GT_COUNT=5,0,0;SUMMED=30.0")
                // A record carrying only one of the two requested annotations.
                + record(300, "RAW_GT_COUNT=0,1,0;SUMMED=7.0")
                // One carrying neither, whose divisor is still fine.
                + record(400, "RAW_GT_COUNT=0,3,0")
                // A divisor whose fields are not whole numbers, read as doubles.
                + record(500, "RAW_GT_COUNT=0,1,1;SUMMED=5.0;AS_QD=4.0");
        final String input = HEADER + records;
        run(dir, "two-annotations", input, List.of("SUMMED", "AS_QD"));
        // One annotation only, so AS_QD keeps its value and gains no average.
        run(dir, "one-annotation", input, List.of("SUMMED"));
        // An annotation the file never declares, whose AVERAGE_ line is added anyway.
        run(dir, "undeclared", input, List.of("NOT_THERE"));
        // A file with no records at all, which is a header rewrite.
        run(dir, "no-records", HEADER, List.of("SUMMED"));
        // A record with no RAW_GT_COUNT, which is refused after earlier records were written.
        run(dir, "missing-counts",
                HEADER + record(100, "RAW_GT_COUNT=0,2,1;SUMMED=30.0") + record(200, "SUMMED=1.0"),
                List.of("SUMMED"));
    }

    static void run(final Path dir, final String label, final String input,
                    final List<String> annotations) throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", in.toString()});
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));

        final Path out = dir.resolve("averaged-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-V", in.toString(), "-O", out.toString()));
        for (final String annotation : annotations) {
            argv.addAll(Arrays.asList("--summed-annotation-to-divide", annotation));
        }
        try {
            new CalculateAverageCombinedAnnotations().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            // The writer is open before the failing record is reached, so what was written stays.
            if (Files.exists(out)) {
                System.out.printf("averaged\t%s-partial=%s%n", label,
                        ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
            }
            return;
        }
        System.out.printf("averaged\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replaceAll("##GATKCommandLine=<[^\n]*>", "##GATKCommandLine=<MASKED>")
                .replace(dir.toString(), "<dir>");
    }
}
