/*
 * VariantsToTable's output, taken from the reference.
 *
 * A VCF in, a tab-separated table out, and almost every column is a decision rather than a copy.
 *
 * Ten behaviours this is built to catch, four of which are not what the arguments suggest.
 *
 *   - NO FIELDS AT ALL MEANS EVERY FIELD IN THE HEADER, not an error: the mandatory columns except
 *     INFO, then every INFO line, then every FORMAT line with GT forced to the front;
 *   - AND ASKING FOR NO GENOTYPE FIELD EMPTIES THE SAMPLE LIST, so `-F CHROM` alone produces one
 *     column and not one per sample;
 *   - A MISSING FIELD IS THE STRING `NA` and not an empty cell, unless --error-if-missing-data
 *     turns it into a UserException naming the record;
 *   - A FIELD ENDING IN `*` IS A WILDCARD over the INFO keys, joined by commas in SORTED order,
 *     which is not the order the record wrote them in;
 *   - QUAL IS `Double.toString` OF THE PHRED SCORE, so an integer quality comes out `50.0`;
 *   - FILTER IS `PASS` FOR AN UNFILTERED RECORD, even one whose column said `.`;
 *   - --split-multi-allelic PRODUCES ONE ROW PER ALTERNATE, and a value is spread across the rows
 *     ONLY when it is a list of exactly the right length: anything else is repeated whole, so a
 *     Number=R field lands identically in every row while a Number=A field is split;
 *   - AND -ASF SUBSETS AN R-TYPE FIELD BY DROPPING ITS FIRST ENTRY, which is the reference's, so
 *     the same field read as -F and as -ASF gives different columns;
 *   - -ASGF AD IS SPECIAL-CASED to `<ref depth>,<alt depth>` per row, which no other field does;
 *   - AND --moltenize REPLACES THE TABLE with four columns and one line per value, numbering the
 *     records from 1 while the header line it printed says `RecordID`.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     line\t<label>\t<one line of the output table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VariantsToTableDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VariantsToTableDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AS_SB,Number=R,Type=Integer,Description=\"Per allele strand bias\">\n"
                    + "##INFO=<ID=AS_QD,Number=A,Type=Float,Description=\"Per alternate quality by depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FILTER=<ID=LowQD,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variantstotable-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VariantsToTableDump: a VCF as a table, from the reference");

        final Path input = writeVcf(dir, "records",
                // A biallelic SNP with an rsID and an integer quality.
                "chr1\t100\trs100\tA\tC\t50\t.\tDP=30;AC=1;AS_SB=10,5;AS_QD=20.0"
                        + "\tGT:GQ:AD\t0/1:60:15,15\t0/0:60:30,0",
                // A multi-allelic record, for the splitting.
                "chr1\t200\t.\tA\tC,G\t60\tPASS\tDP=30;AC=1,2;AS_SB=10,5,3;AS_QD=20.0,30.0"
                        + "\tGT:GQ:AD\t0/1:60:10,10,10\t1/2:60:0,15,15",
                // A filtered record, which is dropped unless --show-filtered.
                "chr1\t300\t.\tA\tC\t70\tLowQD\tDP=30;AC=1;AS_SB=10,5;AS_QD=20.0"
                        + "\tGT:GQ:AD\t0/1:60:15,15\t0/0:60:30,0",
                // An indel with no AC at all, for the missing-field column.
                "chr1\t400\t.\tACGT\tA\t80.5\t.\tDP=30;AS_SB=10,5;AS_QD=20.0"
                        + "\tGT:GQ:AD:FT\t0/1:60:15,15:LowQD\t0/0:60:30,0:PASS");

        // No fields at all, which takes everything the header declares.
        run(dir, "no-fields", input);

        // The mandatory columns and the derived ones.
        run(dir, "mandatory", input, "-F", "CHROM", "-F", "POS", "-F", "ID", "-F", "REF",
                "-F", "ALT", "-F", "QUAL", "-F", "FILTER");
        run(dir, "derived", input, "-F", "TYPE", "-F", "EVENTLENGTH", "-F", "TRANSITION",
                "-F", "HET", "-F", "HOM-REF", "-F", "HOM-VAR", "-F", "NO-CALL",
                "-F", "NSAMPLES", "-F", "NCALLED", "-F", "MULTI-ALLELIC");

        // An INFO field, one that is missing from a record, and the wildcard.
        run(dir, "info", input, "-F", "POS", "-F", "DP", "-F", "AC");
        run(dir, "missing", input, "-F", "POS", "-F", "ZZ");
        run(dir, "missing-is-an-error", input, "-F", "POS", "-F", "ZZ",
                "--error-if-missing-data", "true");
        run(dir, "wildcard", input, "-F", "POS", "-F", "AS_*");

        // The genotype fields.
        run(dir, "genotype", input, "-F", "POS", "-GF", "GT", "-GF", "GQ");
        run(dir, "genotype-filter", input, "-F", "POS", "-GF", "FT");
        run(dir, "genotype-depths", input, "-F", "POS", "-GF", "AD");

        // The filtered record.
        run(dir, "show-filtered", input, "-F", "POS", "-F", "FILTER", "--show-filtered", "true");

        // The splitting, with and without the allele-specific arguments.
        run(dir, "split", input, "-F", "POS", "-F", "ALT", "-F", "AC", "-F", "AS_SB",
                "--split-multi-allelic", "true");
        run(dir, "split-allele-specific", input, "-F", "POS", "-F", "ALT", "-ASF", "AS_SB",
                "-ASF", "AS_QD", "--split-multi-allelic", "true");
        run(dir, "allele-specific-unsplit", input, "-F", "POS", "-ASF", "AS_SB");
        run(dir, "split-genotype-depths", input, "-F", "POS", "-ASGF", "AD",
                "--split-multi-allelic", "true");

        // And the shape that is not a table.
        run(dir, "moltenize", input, "-F", "POS", "-F", "DP", "-GF", "GT", "--moltenize", "true");
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path input,
                    final String... arguments) {
        final Path output = dir.resolve(label + "-out.table");
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(),
                "-O", output.toString()));
        all.addAll(List.of(arguments));
        try {
            new VariantsToTable().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, output);
    }

    static void print(final String label, final Path output) {
        final List<String> lines;
        try {
            lines = Files.readAllLines(output, StandardCharsets.UTF_8);
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
            return;
        }
        for (final String line : lines) {
            System.out.printf("line\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
