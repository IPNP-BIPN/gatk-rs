/*
 * VariantFiltration's output, taken from the reference.
 *
 * A variant transform whose whole job is the FILTER column: JEXL expressions at the site and at the
 * genotype, a clustered-SNP test that looks at neighbours, and a mask read from a second file.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE FILTER SET IS A LinkedHashSet SEEDED FROM THE RECORD'S OWN FILTERS, so an existing
 *     filter survives a run that adds another. THE ORDER IN THE FILE IS NOT THAT ORDER: the writer
 *     sorts, so a record carrying `OldFilter` and gaining `LowQD` comes out `LowQD;OldFilter`, and
 *     two expressions given as LowQD then LowDP come out `LowDP;LowQD`. The insertion order the
 *     tool builds is invisible;
 *   - AN EMPTY SET BECOMES PASS AND NOT NOTHING. The comment in the source says it: "making this
 *     empty set effectively converts the VC to PASS, whereas an unfiltered VC has null filters", so
 *     a record that arrived with no FILTER at all leaves with `PASS`;
 *   - EXCEPT UNDER --invalidate-previous-filters, where an empty set leaves the record UNFILTERED,
 *     a dot rather than PASS: the same emptiness, two different columns;
 *   - THE CLUSTER TEST IGNORES NON-SNPs entirely, both as the candidate and as neighbours, and a
 *     window below 1 disables it;
 *   - A GENOTYPE FILTER WRITES FT AND LEAVES THE CALL, unless --set-filtered-genotype-to-no-call,
 *     which replaces the call itself. AND THE FT COLUMN IS PER RECORD, NOT PER SAMPLE: it appears
 *     only where at least one genotype of that record was filtered, and the unfiltered samples of
 *     that record then carry `PASS` while a record with nothing filtered has no FT field at all;
 *   - AND AN EXISTING FT IS SPLIT AND CARRIED, so a genotype that arrived filtered keeps its old
 *     names beside the new ones;
 *   - --invert-filter-expression FLIPS THE MATCH, not the result: an expression that matched now
 *     does not, so the FILTER column carries the complement rather than the negation of the file;
 *   - A MISSING FIELD IS NOT A FALSE unless --missing-values-evaluate-as-failing, which is what
 *     decides whether a record lacking the annotation is filtered or passed;
 *   - AND THE MASK IS A SECOND FILE, whose absence of overlap is what filters when
 *     --filter-not-in-mask is given, so the same mask means opposite things depending on one flag.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VariantFiltrationDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VariantFiltrationDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FILTER=<ID=OldFilter,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variantfiltration-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VariantFiltrationDump: the FILTER column, from the reference");

        // Records covering: a low and a high annotation, a record with no FILTER, one that already
        // carries a filter, one that already says PASS, an indel, and a genotype already filtered.
        final Path input = writeVcf(dir, "records",
                "chr1\t100\t.\tA\tC\t50\t.\tDP=5;QD=1.0\tGT:GQ\t0/1:20\t0/0:60",
                "chr1\t200\t.\tA\tC\t50\tOldFilter\tDP=50;QD=20.0\tGT:GQ\t0/1:5\t1/1:60",
                "chr1\t300\t.\tA\tC\t50\tPASS\tDP=50;QD=20.0\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t400\t.\tACG\tA\t50\t.\tDP=5;QD=1.0\tGT:GQ\t0/1:20\t0/0:60",
                "chr1\t500\t.\tA\tC\t50\t.\tDP=50\tGT:GQ:FT\t0/1:20:OldGtFilter\t0/0:60:PASS");

        // Three SNPs within 20 bases, then one far away, for the cluster test, with an indel in the
        // middle to show it is ignored.
        final Path clustered = writeVcf(dir, "clustered",
                "chr1\t1000\t.\tA\tC\t50\t.\tDP=50\tGT\t0/1\t0/0",
                "chr1\t1005\t.\tACG\tA\t50\t.\tDP=50\tGT\t0/1\t0/0",
                "chr1\t1010\t.\tA\tC\t50\t.\tDP=50\tGT\t0/1\t0/0",
                "chr1\t1015\t.\tA\tC\t50\t.\tDP=50\tGT\t0/1\t0/0",
                "chr1\t9000\t.\tA\tC\t50\t.\tDP=50\tGT\t0/1\t0/0");

        // A mask file, overlapping one record of the first input.
        final Path mask = writeVcf(dir, "mask",
                "chr1\t200\t.\tA\tG\t50\t.\tDP=1\tGT\t0/1\t0/0");

        // The site expression on its own, then inverted.
        run(dir, "site-filter", input, "--filter-name", "LowQD", "--filter-expression", "QD < 2.0");
        run(dir, "site-filter-inverted", input, "--filter-name", "LowQD",
                "--filter-expression", "QD < 2.0", "--invert-filter-expression", "true");
        // Two expressions, so the order of the FILTER column can be seen.
        run(dir, "two-filters", input,
                "--filter-name", "LowQD", "--filter-expression", "QD < 2.0",
                "--filter-name", "LowDP", "--filter-expression", "DP < 10");
        // A record missing the annotation entirely, with and without the flag.
        run(dir, "missing-values", input, "--filter-name", "LowQD",
                "--filter-expression", "QD < 2.0", "--missing-values-evaluate-as-failing", "true");

        // Genotype filters, with and without the no-call replacement.
        run(dir, "genotype-filter", input, "--genotype-filter-name", "LowGQ",
                "--genotype-filter-expression", "GQ < 30");
        run(dir, "genotype-filter-nocall", input, "--genotype-filter-name", "LowGQ",
                "--genotype-filter-expression", "GQ < 30", "--set-filtered-genotype-to-no-call", "true");

        // The cluster test at two window sizes and one that disables it.
        run(dir, "cluster", clustered, "--cluster-size", "3", "--cluster-window-size", "20");
        run(dir, "cluster-narrow", clustered, "--cluster-size", "3", "--cluster-window-size", "5");
        run(dir, "cluster-disabled", clustered, "--cluster-size", "3", "--cluster-window-size", "0");

        // The mask, both ways round.
        run(dir, "mask", input, "--mask", mask.toString(), "--mask-name", "InMask");
        run(dir, "mask-inverted", input, "--mask", mask.toString(), "--mask-name", "NotInMask",
                "--filter-not-in-mask", "true");

        // And the two ways an empty filter set can be written.
        run(dir, "no-filters", input, "--filter-name", "Never", "--filter-expression", "QD > 1000.0");
        run(dir, "no-filters-invalidated", input, "--filter-name", "Never",
                "--filter-expression", "QD > 1000.0", "--invalidate-previous-filters", "true");
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
            new VariantFiltration().instanceMain(all.toArray(new String[0]));
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
            if (line.startsWith("##GATKCommandLine")) {
                System.out.printf("commandline\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(line.replaceAll("Date=\"[^\"]*\"", "Date=\"MASKED\"")));
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
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
