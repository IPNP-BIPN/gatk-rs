/*
 * Which records SelectVariants keeps, taken from the reference.
 *
 * The tool filters in three places and they are not interchangeable: a first round before anything
 * is subset, two genotype-count gates in the middle, and `--exclude-non-variants` after the subset,
 * with the JEXL expressions on either side of that depending on one flag.
 *
 * Ten behaviours this is built to catch, three of which are defects rather than surprises.
 *
 *   - THE FILTERED-GENOTYPE FRACTION IS AN INTEGER DIVISION. `numFilteredSamples / samples.size()`
 *     is int over int, assigned to a double, so the fraction is 0 unless every sample is filtered
 *     and 1 when they all are. --max-fraction-filtered-genotypes therefore does nothing at any
 *     value above 0, while --max-nocall-fraction, one line below, casts and works;
 *   - AND THE SAME GATE IS ONLY CONSULTED WHEN AN ARGUMENT ASKS FOR IT, so a record with filtered
 *     genotypes survives a run that never mentions them;
 *   - --select-type-to-exclude BEATS --select-type-to-include, since the exclusions are removed
 *     from the inclusion set after it is built, and giving neither selects every type;
 *   - --max-indel-size AND --min-indel-size ARE BOTH ABOUT ABSOLUTE LENGTH and both reject the
 *     record rather than the allele: a record with one 1bp and one 20bp indel is rejected by
 *     --max-indel-size 10 even though one of its alternates would pass;
 *   - --exclude-non-variants RUNS AFTER THE SUBSET, so it drops a record whose alternates are only
 *     carried by samples that were not selected;
 *   - AND A SPANNING DELETION ALONE COUNTS AS NON-VARIANT, so `A -> *` is dropped by the same flag
 *     that keeps `A -> C`;
 *   - THE JEXL EXPRESSIONS ARE OR-ED, not and-ed, and --invert-select inverts EACH of them before
 *     the or, which is not the complement of the whole;
 *   - --select-genotype IS TRUE WHEN ANY GENOTYPE MATCHES, and its result is or-ed with the INFO
 *     expressions rather than and-ed;
 *   - --apply-jexl-filters-first MOVES THE EXPRESSIONS BEFORE THE SUBSET, so an expression over AC
 *     sees the record's own annotation rather than the recomputed one, and the same command line
 *     keeps different records;
 *   - AND --exclude-filtered LOOKS AT THE FILTER COLUMN ONLY, which SelectVariants never writes.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     kept\t<label>\t<comma-joined positions of the records written>
 *     vcfline\t<label>\t<one record line of the output VCF, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SelectVariantsFiltersDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SelectVariantsFiltersDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                    + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FILTER=<ID=LowQD,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("selectvariants-filters-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SelectVariantsFiltersDump: which records survive, from the reference");

        final Path input = writeVcf(dir, "records",
                // 100: a plain SNP, unfiltered, with an rsID.
                "chr1\t100\trs100\tA\tC\t50\t.\tDP=50;AC=2;AF=0.250;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t0/1:60\t0/0:60\t0/0:60",
                // 200: a SNP whose FILTER column already says something.
                "chr1\t200\trs200\tA\tC\t50\tLowQD\tDP=50;AC=1;AF=0.125;AN=8;QD=1.0"
                        + "\tGT:GQ\t0/1:60\t0/0:60\t0/0:60\t0/0:60",
                // 300: a one-base insertion.
                "chr1\t300\t.\tA\tAG\t50\t.\tDP=50;AC=1;AF=0.125;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t0/0:60\t0/0:60\t0/0:60",
                // 400: a twenty-base deletion, for the indel size gates.
                "chr1\t400\t.\tACGTACGTACGTACGTACGTA\tA\t50\t.\tDP=50;AC=1;AF=0.125;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t0/0:60\t0/0:60\t0/0:60",
                // 500: multi-allelic, one SNP and one insertion, so its type is MIXED.
                "chr1\t500\t.\tA\tC,AG\t50\t.\tDP=50;AC=1,1;AF=0.125,0.125;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t0/2:60\t0/0:60\t0/0:60",
                // 600: only s2 and s3 carry the alternate, so subsetting to s0 and s1 makes it
                // non-variant.
                "chr1\t600\t.\tA\tC\t50\t.\tDP=50;AC=2;AF=0.250;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\t0/1:60",
                // 700: a spanning deletion alone.
                "chr1\t700\t.\tA\t*\t50\t.\tDP=50;AC=1;AF=0.125;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t0/0:60\t0/0:60\t0/0:60",
                // 800: two of the four genotypes are filtered, one is a no-call.
                "chr1\t800\t.\tA\tC\t50\t.\tDP=50;AC=1;AF=0.125;AN=8;QD=20.0"
                        + "\tGT:GQ:FT\t0/1:60:PASS\t0/0:60:LowQD\t0/0:60:LowQD\t./.:60:PASS",
                // 900: three no-calls out of four.
                "chr1\t900\t.\tA\tC\t50\t.\tDP=50;AC=1;AF=0.125;AN=8;QD=20.0"
                        + "\tGT:GQ\t0/1:60\t./.:60\t./.:60\t./.:60");

        // Nothing at all, for the baseline.
        run(dir, "no-filter", input);

        // The types, each way round.
        run(dir, "type-snp", input, "--select-type-to-include", "SNP");
        run(dir, "type-indel", input, "--select-type-to-include", "INDEL");
        run(dir, "type-exclude-snp", input, "--select-type-to-exclude", "SNP");
        run(dir, "type-include-and-exclude", input, "--select-type-to-include", "SNP",
                "--select-type-to-include", "INDEL", "--select-type-to-exclude", "SNP");

        // The allele-count restriction.
        run(dir, "biallelic", input, "--restrict-alleles-to", "BIALLELIC");
        run(dir, "multiallelic", input, "--restrict-alleles-to", "MULTIALLELIC");

        // The indel sizes, which reject the record rather than the allele.
        run(dir, "max-indel-size", input, "--max-indel-size", "10");
        run(dir, "min-indel-size", input, "--min-indel-size", "5");

        // The rsIDs.
        run(dir, "keep-ids", input, "--keep-ids", "rs100");
        run(dir, "exclude-ids", input, "--exclude-ids", "rs100");

        // The FILTER column.
        run(dir, "exclude-filtered", input, "--exclude-filtered", "true");

        // The two flags that look at the genotypes, and the fraction that does not work.
        run(dir, "max-filtered-genotypes", input, "--max-filtered-genotypes", "1");
        run(dir, "max-fraction-filtered-genotypes", input,
                "--max-fraction-filtered-genotypes", "0.1");
        run(dir, "max-nocall-number", input, "--max-nocall-number", "1");
        run(dir, "max-nocall-fraction", input, "--max-nocall-fraction", "0.1");

        // Non-variant after the subset, and the spanning deletion.
        run(dir, "exclude-non-variants", input, "--exclude-non-variants", "true");
        run(dir, "exclude-non-variants-subset", input, "-sn", "s0", "-sn", "s1",
                "--exclude-non-variants", "true");

        // The JEXL expressions, or-ed, inverted, and moved before the subset.
        run(dir, "select-one", input, "-select", "QD > 10.0");
        run(dir, "select-two", input, "-select", "QD > 10.0", "-select", "AC > 1");
        run(dir, "select-inverted", input, "-select", "QD > 10.0", "--invert-select", "true");
        run(dir, "select-two-inverted", input, "-select", "QD > 10.0", "-select", "AC > 1",
                "--invert-select", "true");
        run(dir, "select-genotype", input, "-select-genotype", "GQ > 55");
        // The same expression over AC, before and after the subset recomputes it.
        run(dir, "select-ac-after-subset", input, "-sn", "s0", "-select", "AC > 1");
        run(dir, "select-ac-before-subset", input, "-sn", "s0", "-select", "AC > 1",
                "--apply-jexl-filters-first", "true");
        // And an expression the engine cannot compile.
        run(dir, "select-unparseable", input, "-select", "QD >");
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
        final Path output = dir.resolve(label + "-out.vcf");
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(),
                "-O", output.toString()));
        all.addAll(List.of(arguments));
        try {
            new SelectVariants().instanceMain(all.toArray(new String[0]));
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
        final List<String> kept = new ArrayList<>();
        for (final String line : lines) {
            if (line.startsWith("#")) {
                continue;
            }
            kept.add(line.split("\t")[1]);
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
        System.out.printf("kept\t%s\t%s%n", label, String.join(",", kept));
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
