/*
 * CountFalsePositives' table, taken from the reference.
 *
 * The tool counts unfiltered records into two buckets and divides each by the target territory,
 * and every one of those three words hides a decision.
 *
 * Seven behaviours this is built to catch.
 *
 *   - EVERYTHING THAT IS NOT AN INDEL IS A SNP. `apply` is `if (variant.isIndel()) indel++ else
 *     snp++`, and `isIndel()` is true for the INDEL type alone, so an MNP, a symbolic allele, a
 *     mixed record and a record with no alternate at all are all counted in the `snp` column;
 *   - THE ID IS THE FILE NAME, not a sample: `drivingVariantFile.getBaseName()`, which is the name
 *     with one extension removed, so a `.vcf.gz` input keeps its `.vcf`;
 *   - THE TERRITORY IS THE MERGED INTERVALS' and not the ones typed, so two overlapping -L
 *     arguments contribute their union once;
 *   - AND IT IS COUNTED IN BASES, `SimpleInterval.size()` being end - start + 1;
 *   - THE RATES ARE PER MEGABASE, `count / territory * 1e6`, computed in that order, and written
 *     through `DataLine.set(double)`, whose dead rounding branch means an integral rate keeps its
 *     `.0`;
 *   - -L IS REQUIRED, `requiresIntervals` being true, and its absence is refused before any record
 *     is read;
 *   - AND AN OUTPUT THAT CANNOT BE OPENED IS REFUSED AFTER THE WHOLE TRAVERSAL, in a plain
 *     UserException whose message is "Encountered an IO exception while opening" and the path.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CountFalsePositivesDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.CountFalsePositives;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CountFalsePositivesDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the block\">\n"
                    + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak_evidence,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "##contig=<ID=chr2,length=900>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("countfalsepositives-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CountFalsePositivesDump: two buckets, and what falls in the second");

        // One record of every shape the two branches can see: a SNP, an insertion, a deletion, an
        // MNP, a symbolic deletion, a mixed record, a record with no alternate, and a filtered SNP
        // and a filtered indel that neither branch may count.
        final Path mixed = writeVcf(dir, "mixed",
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tACC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tACC\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t400\t.\tAC\tGT\t50\tPASS\t.\tGT\t0/1",
                "chr1\t500\t.\tA\t<DEL>\t50\tPASS\tEND=520\tGT\t0/1",
                "chr1\t600\t.\tA\tC,ACC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t700\t.\tA\t.\t50\tPASS\t.\tGT\t0/0",
                "chr1\t800\t.\tA\tC\t50\tweak_evidence\t.\tGT\t0/1",
                "chr1\t900\t.\tA\tACC\t50\tweak_evidence\t.\tGT\t0/1");

        // Nothing that passes, so both counts are zero and both rates with them.
        final Path allFiltered = writeVcf(dir, "all-filtered",
                "chr1\t100\t.\tA\tC\t50\tweak_evidence\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tACC\t50\tweak_evidence\t.\tGT\t0/1");

        // One PASS SNP, over a territory small enough that the rate is not integral.
        final Path oneSnp = writeVcf(dir, "one-snp",
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");

        // The same records under a name with two extensions, which is what the id is taken from.
        final Path twoExtensions = writeVcf(dir, "two-extensions.vcf",
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");

        // The whole of chr1, which is the baseline for every count.
        run(dir, "whole-contig", mixed, "whole-contig.table", "-L", "chr1");
        run(dir, "all-filtered", allFiltered, "all-filtered.table", "-L", "chr1");

        // The territory is the merged intervals: 1-200 and 150-400 are one interval of 400 bases,
        // and 1-100 with 900-1000 are two of 101 each.
        run(dir, "overlapping-intervals", oneSnp, "overlapping.table",
                "-L", "chr1:1-200", "-L", "chr1:150-400");
        run(dir, "disjoint-intervals", oneSnp, "disjoint.table",
                "-L", "chr1:1-100", "-L", "chr1:900-1000");

        // A territory of three bases, whose rate is a repeating decimal, and one of exactly a
        // megabase's worth of denominators for an integral one.
        run(dir, "small-territory", oneSnp, "small.table", "-L", "chr1:98-100");
        run(dir, "integral-rate", oneSnp, "integral.table", "-L", "chr1:1-100");

        // Selection also decides the counts, not just the denominator.
        run(dir, "interval-selects-some", mixed, "some.table", "-L", "chr1:1-350");

        // The id, from a file name carrying two extensions.
        run(dir, "two-extensions", twoExtensions, "two-extensions.table", "-L", "chr1");

        // The two refusals.
        run(dir, "no-intervals", mixed, "no-intervals.table");
        run(dir, "output-is-a-directory", mixed, ".", "-L", "chr1");
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

    /** One run, reported as the table it wrote or as the exception it threw. */
    static void run(final Path dir, final String label, final Path input, final String output,
                    final String... arguments) throws Exception {
        final Path file = dir.resolve(output);
        final List<String> all = new ArrayList<>(List.of(
                "-V", input.toString(), "-O", file.toString()));
        all.addAll(List.of(arguments));
        try {
            new CountFalsePositives().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
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
