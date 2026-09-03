/*
 * The header SelectVariants writes, taken from the reference.
 *
 * The five existing SelectVariants suites all skip the `#` lines, and for a good reason: the
 * default header carries a `##GATKCommandLine` line with the wall-clock Date of the run, so a
 * golden of the whole header would flake by construction. That is exactly why the header is the
 * part of this tool nothing measures, and it is not a small part: `createVCFHeaderLineList` merges,
 * adds, REPLACES and removes lines, and the replacements are the ones a port gets wrong.
 *
 * Six behaviours this is built to catch.
 *
 *   - AC, AF, AN AND DP ARE REPLACED, not merged. Whatever the input declared for those four INFO
 *     ids is removed and htsjdk's own standard line is added, so an input whose `##INFO=<ID=AC>`
 *     says `Number=A,Type=Integer,Description="Allele count"` comes out with the standard
 *     description instead. AF is added even when the input never had it;
 *   - --keep-original-ac ADDS THREE LINES (AC_Orig, AF_Orig, AN_Orig) and --keep-original-dp adds
 *     one (DP_Orig), whether or not any record ends up carrying them;
 *   - --drop-info-annotation AND --drop-genotype-annotation REMOVE THE HEADER LINE TOO, and they
 *     run AFTER the four replacements, so dropping AC removes the standard line that had just been
 *     put there;
 *   - THE SAMPLE COLUMNS ARE THE SELECTED SET, sorted, which is what `new VCFHeader(lines, samples)`
 *     does with the TreeSet `createSampleNameInclusionList` returns;
 *   - THE CONTIG LINES ARE REWRITTEN FROM THE DICTIONARY even with no reference, because
 *     `updateHeaderContigLines` runs whenever the input header has one;
 *   - AND `--add-output-vcf-command-line false` REMOVES THE `##source` LINE AS WELL, because both
 *     come from `getDefaultToolVCFHeaderLines` and it returns an empty set.
 *
 * The last case runs with the command line ON and prints it with the Date and the Version elided,
 * so that where the line goes and what it holds are both measured without the golden moving.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     header\t<label>\t<one header line of the output VCF, escaped>
 *     samples\t<label>\t<the sample columns, comma-joined>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SelectVariantsHeaderDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SelectVariantsHeaderDump {

    /*
     * The input's own AC, AN and DP lines are DELIBERATELY not the standard ones: a port that
     * merges them through instead of replacing them produces this file's descriptions rather than
     * htsjdk's, and no other case can tell the difference. AF is absent for the same reason, since
     * the tool adds it whether or not the input had one.
     */
    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Not the standard depth\">\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Not the standard count\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Not the standard number\">\n"
                    + "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Quality by depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=XX,Number=1,Type=Integer,Description=\"An annotation to drop\">\n"
                    + "##FILTER=<ID=LowGQ,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "##contig=<ID=chr2,length=90000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("selectvariants-header-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SelectVariantsHeaderDump: the header the tool writes, from the reference");

        final Path input = writeVcf(dir, "records",
                "chr1\t100\t.\tA\tC\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7",
                "chr2\t50\t.\tA\tG\t50\t.\tDP=30;AC=1;AN=6;QD=20.0"
                        + "\tGT:GQ:DP:XX\t0/1:60:10:7\t0/0:60:10:7\t0/0:60:10:7");

        // The four replacements, with nothing else asked for.
        run(dir, "plain", input);
        // The three lines and the one line the keep-original arguments add.
        run(dir, "keep-original-ac", input, "--keep-original-ac", "true");
        run(dir, "keep-original-dp", input, "--keep-original-dp", "true");
        run(dir, "keep-original-both", input, "--keep-original-ac", "true",
                "--keep-original-dp", "true");
        // The drops, which run after the replacements: `AC` removes the line just added.
        run(dir, "drop-info-qd", input, "--drop-info-annotation", "QD");
        run(dir, "drop-info-ac", input, "--drop-info-annotation", "AC");
        run(dir, "drop-genotype", input, "--drop-genotype-annotation", "XX");
        // A drop of something the header never declared, which removes nothing and refuses nothing.
        run(dir, "drop-absent", input, "--drop-info-annotation", "NOPE");
        // The sample columns, which are the selected set rather than the file's.
        run(dir, "subset", input, "-sn", "s2", "-sn", "s0");
        // Sites only, where the sample columns go entirely.
        run(dir, "sites-only", input, "--sites-only-vcf-output", "true");
        // The command line, with its Date and Version elided so the golden holds still.
        runWithCommandLine(dir, "command-line", input);
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
        final List<String> all = new ArrayList<>(List.of("--add-output-vcf-command-line", "false"));
        all.addAll(List.of(arguments));
        execute(dir, label, input, all);
    }

    static void runWithCommandLine(final Path dir, final String label, final Path input) {
        execute(dir, label, input, List.of());
    }

    static void execute(final Path dir, final String label, final Path input,
                        final List<String> arguments) {
        final Path output = dir.resolve(label + "-out.vcf");
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(),
                "-O", output.toString()));
        all.addAll(arguments);
        try {
            new SelectVariants().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, output);
    }

    /** The `##GATKCommandLine` line with the two fields that move replaced by a fixed word. */
    static String elide(final String line) {
        return line.replaceAll("Version=\"[^\"]*\"", "Version=\"ELIDED\"")
                .replaceAll("Date=\"[^\"]*\"", "Date=\"ELIDED\"");
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
            if (!line.startsWith("#")) {
                break;
            }
            if (line.startsWith("#CHROM")) {
                final String[] field = line.split("\t");
                final List<String> samples = new ArrayList<>();
                for (int index = 9; index < field.length; index++) {
                    samples.add(field[index]);
                }
                System.out.printf("samples\t%s\t%s%n", label, String.join(",", samples));
                continue;
            }
            System.out.printf("header\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(elide(line)));
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
