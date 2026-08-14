/*
 * FilterVariantTranches, taken from the reference.
 *
 * A two-pass walker: the first pass collects the input's own score for every record that overlaps a
 * resource, the second turns those scores into cutoffs and filters against them. Eight behaviours
 * this is built to catch.
 *
 *   - THE CUTOFF IS AN INDEX INTO THE RESOURCE SCORES, TRUNCATED. The scores are sorted DESCENDING
 *     and the index is `(int)((t / 100.0) * (size - 1))`, so five scores and a tranche of 50 give
 *     index 2 rather than the median by any other definition, and a higher tranche means a lower
 *     cutoff;
 *   - THE SCORE STORED IS THE INPUT RECORD'S, NOT THE RESOURCE'S. The first pass reads
 *     `variant.getAttribute(infoKey)` after matching against the resource, so the cutoffs are
 *     quantiles of the INPUT's scores at resource sites;
 *   - AND IT IS STORED ONCE: the loop `return`s on its first match, so a record overlapping two
 *     resources contributes one score;
 *   - MEMBERSHIP IS DECIDED BY THE FIRST CUTOFF ALONE, `score <= cutoffs.get(0)`, which is the
 *     cutoff of the SMALLEST tranche because the list was sorted; and the comparison is `<=`, so a
 *     score exactly on the cutoff is filtered;
 *   - THE TRANCHE LIST IS DEDUPLICATED AND SORTED, and refused outright unless every value is in
 *     [0, 100): 100 itself is a `CommandLineException`;
 *   - THE FILTER NAME IS BUILT WITH `%.2f`, `<info key>_<type>_Tranche_<t1>_<t2>`, and the last
 *     tranche always runs to `100.00`;
 *   - A RECORD WITH NO SCORE IS STILL WRITTEN AND STILL PASSES, `passFilters()` being applied to
 *     anything the tool did not filter itself;
 *   - THE FOUR REFUSALS OF `afterFirstPass` each say something different: no scored variant at all,
 *     no resource overlap at all, SNPs with no SNP resource, indels with no indel resource; and an
 *     info key the input header does not declare is refused before any of them;
 *   - AND `--invalidate-previous-filters` clears both the header's FILTER lines and each record's.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     vcfline\t<run>\t<one record line of the output vcf, escaped>
 *     filter\t<run>\t<one ##FILTER line of the output vcf, escaped>
 *     error\t<run>\t<exception class>:<message>
 *
 * Usage: FilterVariantTranchesDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.vqsr.FilterVariantTranches;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FilterVariantTranchesDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"the score to filter on\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filter-variant-tranches-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FilterVariantTranchesDump: five scores, two tranches, three bands");

        // Five SNPs and three indels, all of them at resource sites, plus one record with no score
        // and one already filtered.
        final Path variants = writeVcf(dir, "variants", HEADER,
                "chr1\t100\t.\tA\tC\t50\t.\tSCORE=5.0\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\t.\tSCORE=4.0\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\t.\tSCORE=3.0\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\t.\tSCORE=2.0\tGT\t0/1",
                "chr1\t500\t.\tA\tC\t50\tweak\tSCORE=1.0\tGT\t0/1",
                "chr1\t600\t.\tACC\tA\t50\t.\tSCORE=9.0\tGT\t0/1",
                "chr1\t700\t.\tACC\tA\t50\t.\tSCORE=8.0\tGT\t0/1",
                "chr1\t800\t.\tACC\tA\t50\t.\tSCORE=7.0\tGT\t0/1",
                // No score at all: counted nowhere and filtered by nothing.
                "chr1\t900\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path resource = writeVcf(dir, "resource", HEADER,
                "chr1\t100\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t400\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t500\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t600\t.\tACC\tA\t50\t.\t.\tGT\t0/1",
                "chr1\t700\t.\tACC\tA\t50\t.\t.\tGT\t0/1",
                "chr1\t800\t.\tACC\tA\t50\t.\t.\tGT\t0/1");

        // The same input with nothing scored.
        final Path unscored = writeVcf(dir, "unscored", HEADER,
                "chr1\t100\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        // A resource that overlaps nothing.
        final Path elsewhere = writeVcf(dir, "resource-elsewhere", HEADER,
                "chr1\t950\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        // A resource with indels alone, against an input holding SNPs too.
        final Path indelsOnly = writeVcf(dir, "resource-indels-only", HEADER,
                "chr1\t600\t.\tACC\tA\t50\t.\t.\tGT\t0/1",
                "chr1\t700\t.\tACC\tA\t50\t.\t.\tGT\t0/1",
                "chr1\t800\t.\tACC\tA\t50\t.\t.\tGT\t0/1");

        run(dir, "two-tranches", variants, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--snp-tranche", "99.0",
                        "--indel-tranche", "50.0"));
        // The same two tranches, given out of order and one of them twice.
        run(dir, "unsorted-and-repeated", variants, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "99.0", "--snp-tranche", "50.0",
                        "--snp-tranche", "99.0", "--indel-tranche", "50.0"));
        run(dir, "invalidate-previous-filters", variants, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--snp-tranche", "99.0",
                        "--indel-tranche", "50.0", "--invalidate-previous-filters"));
        // One tranche only: every filtered record lands in the same band.
        run(dir, "one-tranche", variants, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--indel-tranche", "50.0"));

        run(dir, "tranche-of-a-hundred", variants, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "100.0", "--indel-tranche", "50.0"));
        run(dir, "info-key-not-in-header", variants, resource,
                List.of("--info-key", "MISSING", "--snp-tranche", "50.0", "--indel-tranche", "50.0"));
        run(dir, "nothing-scored", unscored, resource,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--indel-tranche", "50.0"));
        run(dir, "no-overlap", variants, elsewhere,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--indel-tranche", "50.0"));
        run(dir, "snps-without-snp-resources", variants, indelsOnly,
                List.of("--info-key", "SCORE", "--snp-tranche", "50.0", "--indel-tranche", "50.0"));
    }

    static Path writeVcf(final Path dir, final String label, final String header,
                         final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path variants, final Path resource,
                    final List<String> extra) {
        final Path output = dir.resolve(label + ".out.vcf");
        final List<String> all = new ArrayList<>(List.of(
                "-V", variants.toString(),
                "--resource", resource.toString(),
                "-O", output.toString()));
        all.addAll(extra);
        try {
            new FilterVariantTranches().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                System.out.printf("cause\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
            }
            return;
        }
        try {
            for (final String line : Files.readAllLines(output, StandardCharsets.UTF_8)) {
                if (line.startsWith("##FILTER=")) {
                    System.out.printf("filter\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                } else if (!line.startsWith("#")) {
                    System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                }
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
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
