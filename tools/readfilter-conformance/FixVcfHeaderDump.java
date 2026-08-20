/*
 * FixVcfHeader, taken from the reference.
 *
 * A VCF whose header does not declare everything its records use, written back out with the
 * missing declarations invented, or with a header taken wholesale from another file.
 *
 * Nine behaviours this is built to catch.
 *
 *   - AN INVENTED LINE IS ALWAYS `Number=.` AND `Type=String`, whatever the value looks like, and
 *     its description is a fixed sentence naming the tool;
 *   - THE STANDARD FORMAT LINES ARE ADDED WHATEVER THE FILE USES. `addStandardFormatLines` puts
 *     GT, AD, DP, GQ, PL and the rest into the output header even for a file that carries none of
 *     them, so the output declares more than the input ever needed;
 *   - A FILTER ON A RECORD IS A DECLARATION TO INVENT, and `PASS` is not: the filters of a passing
 *     record are empty, so nothing is added for it;
 *   - CHECK_FIRST_N_RECORDS STOPS THE SEARCH EARLY, so a key that appears only in a later record
 *     is not declared and the write then fails on it;
 *   - THE HEADER FILE REPLACES THE HEADER WHOLESALE. Nothing of the input's header survives except
 *     its samples, and the tool warns about lines it is adding without ever comparing the two;
 *   - ENFORCE_SAME_SAMPLES IS ON BY DEFAULT and compares the SORTED sample lists PAIRWISE, naming
 *     the index of the first that differs;
 *   - WITH IT OFF, THE INPUT'S SAMPLES ARE KEPT and the header file's are discarded, so a header
 *     file with no samples at all still writes a file with the input's columns;
 *   - A HEADER FILE AND CHECK_FIRST_N_RECORDS TOGETHER ARE A COMMAND-LINE REFUSAL;
 *   - AND THE WRITER IS BUILT WITH ALLOW_MISSING_FIELDS_IN_HEADER UNSET, so anything the fixing
 *     missed is a refusal at the record that uses it rather than a silently written file.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     header\t<label>=<the replacement header file, escaped>
 *     fixed\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FixVcfHeaderDump
 */

import picard.vcf.FixVcfHeader;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FixVcfHeaderDump {

    static final String COLUMNS =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n";

    /** A header declaring nothing but the contig and GT. */
    static final String BARE =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##contig=<ID=chr1,length=240>\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("fix-vcf-header-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FixVcfHeaderDump: a VCF header rebuilt from the records under it");

        // Undeclared INFO keys, an undeclared FORMAT key and an undeclared filter.
        final String undeclared = BARE + COLUMNS
                + "chr1\t10\t.\tA\tC\t50\tPASS\tXX=1;YY=a,b\tGT:ZZ\t0/1:7\n"
                + "chr1\t20\t.\tA\tG\t50\tmy_filter\tXX=2\tGT\t0/1\n";
        run("undeclared", dir, undeclared, null);
        // Only the first record examined, so the second's filter is never declared.
        run("first-one-record", dir, undeclared, null, "N=1");
        // A file whose header already declares everything.
        run("nothing-missing", dir,
                BARE.replace("##contig", "##INFO=<ID=XX,Number=1,Type=Integer,Description=\"x\">\n##contig")
                        + COLUMNS + "chr1\t10\t.\tA\tC\t50\tPASS\tXX=1\tGT\t0/1\n", null);
        // A file with no records at all, which still gains the standard FORMAT lines.
        run("no-records", dir, BARE + COLUMNS, null);

        // A replacement header that declares everything, with the same sample.
        final String replacement = "##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=ZZ,Number=1,Type=Integer,Description=\"z\">\n"
                + "##FILTER=<ID=my_filter,Description=\"mine\">\n"
                + "##INFO=<ID=XX,Number=1,Type=Integer,Description=\"x\">\n"
                + "##INFO=<ID=YY,Number=.,Type=String,Description=\"y\">\n"
                + "##contig=<ID=chr1,length=240>\n"
                + COLUMNS;
        run("replacement-header", dir, undeclared, replacement);
        // The same header carrying a different sample name.
        run("different-sample", dir, undeclared,
                replacement.replace("\tNA1\n", "\tOTHER\n"));
        // The same, with the sample check turned off, which keeps the input's samples.
        run("different-sample-unenforced", dir, undeclared,
                replacement.replace("\tNA1\n", "\tOTHER\n"), "ENFORCE_SAME_SAMPLES=false");
        // A header file with no samples at all, unenforced.
        run("sites-only-header", dir, undeclared,
                replacement.replace("\tFORMAT\tNA1\n", "\n"), "ENFORCE_SAME_SAMPLES=false");
        // A replacement header that does not declare one of the keys, which the writer refuses.
        run("incomplete-header", dir, undeclared,
                replacement.replace("##INFO=<ID=YY,Number=.,Type=String,Description=\"y\">\n", ""));
        // Both a header file and a record limit, which is the command-line refusal.
        run("header-and-limit", dir, undeclared, replacement, "N=1");
    }

    static void run(final String label, final Path dir, final String input, final String header,
                    final String... extra) throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final List<String> argv = new ArrayList<>(Arrays.asList("I=" + in));
        if (header != null) {
            final Path headerPath = dir.resolve(label + "-header.vcf");
            Files.writeString(headerPath, header, StandardCharsets.UTF_8);
            System.out.printf("header\t%s=%s%n", label, ReferenceQueryDump.escape(header));
            argv.add("H=" + headerPath);
        }
        final Path out = dir.resolve("fixed-" + label + ".vcf");
        argv.addAll(Arrays.asList("O=" + out, "CREATE_INDEX=false"));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new FixVcfHeader().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("fixed\t%s=%s%n", label, ReferenceQueryDump.escape(Files.readString(out)));
    }
}
