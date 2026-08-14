/*
 * The order and the shape of what SelectVariants writes, taken from the reference.
 *
 * Records do not go straight to the writer. They go into a PriorityQueue ordered by start, which is
 * drained only as far as the record currently being read, because trimming can move a record to the
 * RIGHT and a file written in input order would then be unsorted.
 *
 * Eight behaviours this is built to catch, three of which are only visible in the order.
 *
 *   - THE OUTPUT IS RESORTED BY START. A record whose alleles share a prefix moves right when it is
 *     trimmed, and a record that came after it in the file can then come before it in the output;
 *   - THE QUEUE IS DRAINED WITH `<=`, so a pending record at the same start as the record being
 *     read is written FIRST, before that record is even filtered;
 *   - AND A CONTIG CHANGE DRAINS IT ENTIRELY, whatever the starts say, which is the only thing
 *     keeping the queue from comparing positions across contigs;
 *   - TRIMMING ONLY HAPPENS ON THE SUBSETTING PATH: a run that selects every sample leaves the
 *     alleles untouched, so the same file comes out in the same order, and adding one `-sn` both
 *     shortens the alleles and reorders the file;
 *   - --set-filtered-gt-to-nocall REPLACES THE CALL AND KEEPS THE FT, so the genotype comes out
 *     `./.` still carrying the filter name that made it so;
 *   - AND IT RECOMPUTES AC, AN AND AF ITSELF, but only when it actually replaced something:
 *     `setFilteredGenotypeToNocall` calls calculateChromosomeCounts with removeStaleValues, so a
 *     whole-cohort run that recomputes nothing else still rewrites the counts of a record whose
 *     genotypes it no-called, and ADDS an AF the input never had;
 *   - --drop-info-annotation REMOVES A FIELD THE TOOL ITSELF WOULD HAVE WRITTEN, so dropping AC on
 *     a subset run removes the recomputed one rather than the original;
 *   - AND --drop-genotype-annotation REBUILDS EVERY GENOTYPE THROUGH `noAttributes`, which drops
 *     the named key and keeps the rest, GT, GQ, DP, AD and PL included since those are not
 *     extended attributes.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     order\t<label>\t<comma-joined positions of the records written, in file order>
 *     vcfline\t<label>\t<one record line of the output VCF, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SelectVariantsOutputDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SelectVariantsOutputDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                    + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FORMAT=<ID=XX,Number=1,Type=Integer,Description=\"An annotation to drop\">\n"
                    + "##FILTER=<ID=LowGQ,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "##contig=<ID=chr2,length=90000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("selectvariants-output-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SelectVariantsOutputDump: the order and shape of the output, from the reference");

        final Path input = writeVcf(dir, "records",
                // 100: four bases of shared prefix, so trimming moves it to 104.
                "chr1\t100\t.\tACGTA\tACGTAC\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7",
                // 101: a plain SNP, one base after the record above and four before where it lands.
                "chr1\t101\t.\tA\tC\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7",
                // 104: the same start the first record trims to, to show which of the two is
                // written first.
                "chr1\t104\t.\tA\tG\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7",
                // 200: a record with two filtered genotypes, for the no-call replacement.
                "chr1\t200\t.\tA\tC\t50\t.\tDP=30;AC=2;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:FT:XX\t0/1:60:10:LowGQ:7\t0/1:60:10:PASS:7\t0/0:60:10:LowGQ:7",
                // A second contig, which drains the queue whatever the positions say.
                "chr2\t50\t.\tACGTA\tACGTAC\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7",
                "chr2\t51\t.\tA\tC\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7");

        // The whole cohort, which never subsets and so never trims.
        run(dir, "all-samples", input);
        // One sample dropped, which is what turns the trimming on.
        run(dir, "subset", input, "-sn", "s0", "-sn", "s1");
        // The same subset with the trimming turned off again.
        run(dir, "subset-preserve", input, "-sn", "s0", "-sn", "s1",
                "--preserve-alleles", "true");

        // The filtered genotypes, replaced and not.
        run(dir, "set-filtered-to-nocall", input, "--set-filtered-gt-to-nocall", "true");
        run(dir, "set-filtered-to-nocall-subset", input, "-sn", "s0", "-sn", "s1",
                "--set-filtered-gt-to-nocall", "true");

        // The annotations, dropped from the INFO and from the genotypes.
        run(dir, "drop-info", input, "--drop-info-annotation", "QD");
        run(dir, "drop-info-recomputed", input, "-sn", "s0", "-sn", "s1",
                "--drop-info-annotation", "AC");
        run(dir, "drop-genotype", input, "--drop-genotype-annotation", "XX");
        run(dir, "drop-both", input, "--drop-info-annotation", "QD",
                "--drop-genotype-annotation", "XX");
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
        final List<String> order = new ArrayList<>();
        for (final String line : lines) {
            if (line.startsWith("#")) {
                continue;
            }
            final String[] field = line.split("\t");
            order.add(field[0] + ":" + field[1]);
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
        System.out.printf("order\t%s\t%s%n", label, String.join(",", order));
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
