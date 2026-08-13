/*
 * SelectVariants' sample selection, taken from the reference.
 *
 * Which genotype columns survive `-sn`, `-se`, `-xl-sn` and `-xl-se`, and which combinations refuse.
 * The list is built once in `createSampleNameInclusionList` before any record is read, and the
 * decisions it makes are not the ones the argument names suggest.
 *
 * Eight behaviours this is built to catch, four of which are worse than they look.
 *
 *   - THE EXPRESSIONS ARE UNANCHORED `find()` AND NOT `matches()`: Utils.filterCollectionByExpressions
 *     compiles each expression and calls `find()`, so `-se s1` selects `xs10` as well as `s1`. A
 *     sample list written as names rather than patterns silently reaches further than it reads;
 *   - AN EXPRESSION THAT MATCHES NOTHING SELECTS EVERYTHING, because the empty accumulated set is
 *     the same emptiness as "no sample was requested": `-se zzz` gets every sample in the file,
 *     not none;
 *   - --allow-nonoverlapping-command-line-samples DOES THE SAME, when every name given is missing:
 *     the missing names are removed from the accumulated set, the set is then empty, and the empty
 *     set is read as "all". Asking for one sample that does not exist and permitting it outputs
 *     THE WHOLE COHORT;
 *   - AND THE MISSING-NAME CHECK LOOKS AT `-sn` ONLY, never at the expressions, so an expression
 *     matching nothing is never reported while a name matching nothing is a refusal;
 *   - THE OUTPUT ORDER IS THE SORTED ORDER, NOT THE COMMAND LINE'S: both the header's sample set
 *     and the accumulated set are TreeSets, so `-sn tumor -sn s0` writes s0 before tumor;
 *   - EXCLUSION TAKES PRECEDENCE OVER INCLUSION and empties the set, which is then a
 *     UserException rather than an empty output, but only when something was explicitly included:
 *     excluding every sample with no include given hits the same refusal by the second half of the
 *     `noSamplesSpecified` conjunction;
 *   - THE REFUSAL FOR MISSING NAMES IS A FOUR-PARAGRAPH BadInput, joined with `%n%n`, listing the
 *     names comma-separated in the order they were given rather than sorted;
 *   - AND AN UNCOMPILABLE EXPRESSION IS A PatternSyntaxException, thrown from the pattern compiler
 *     and not caught into a UserException, so the tool reports the regex engine's own message.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     expressions\t<label>\t<comma-joined result of Utils.filterCollectionByExpressions>
 *     samples\t<label>\t<comma-joined sample columns of the output header>
 *     vcfline\t<label>\t<one record line of the output VCF, escaped>
 *     error\t<label>\t<exception class>:<message, escaped>
 *
 * Usage: SampleSelectionDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;
import org.broadinstitute.hellbender.utils.Utils;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class SampleSelectionDump {

    /*
     * Six samples, in an order that is not their sorted order, and named so that one is a substring
     * of another: `s1` occurs inside `xs10`, which is what makes the unanchored match visible.
     */
    static final List<String> SAMPLES = List.of("tumor", "s1", "NA12891", "xs10", "s0", "NA12878");

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    + String.join("\t", SAMPLES) + "\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sampleselection-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SampleSelectionDump: which genotype columns survive, from the reference");

        // The expression matcher on its own, so the unanchored match is visible without a tool.
        expressions("substring", "s1");
        expressions("anchored", "^s1$");
        expressions("prefix", "^NA");
        expressions("nothing", "zzz");
        expressions("exact-name", "NA12878");
        expressions("two", "^s", "tumor");
        expressions("everything", ".");

        final Path input = writeVcf(dir, "records",
                "chr1\t100\t.\tA\tC\t50\t.\tDP=50\tGT:GQ\t0/1:20\t0/0:60\t0/1:30\t1/1:40\t0/0:50\t0/1:60",
                "chr1\t200\t.\tA\tC\t50\t.\tDP=50\tGT:GQ\t0/0:20\t0/0:60\t0/0:30\t0/0:40\t0/0:50\t0/0:60");

        // No selection at all, which is the untouched record.
        run(dir, "all-samples", input);
        // One name, then two given in an order that is not the sorted one.
        run(dir, "one-name", input, "-sn", "s1");
        run(dir, "two-names-reversed", input, "-sn", "tumor", "-sn", "s0");
        // The expressions, unanchored and anchored.
        run(dir, "expression-substring", input, "-se", "s1");
        run(dir, "expression-anchored", input, "-se", "^s1$");
        run(dir, "expression-prefix", input, "-se", "^NA");
        // An expression matching nothing, which is not the same as selecting nothing.
        run(dir, "expression-matches-nothing", input, "-se", "zzz");
        // A name and an expression together.
        run(dir, "name-and-expression", input, "-sn", "tumor", "-se", "^NA");

        // A name that is not in the header, refused and then permitted.
        run(dir, "missing-name", input, "-sn", "ghost");
        run(dir, "missing-name-allowed", input, "-sn", "ghost",
                "--allow-nonoverlapping-command-line-samples", "true");
        run(dir, "missing-and-present-allowed", input, "-sn", "ghost", "-sn", "s1",
                "--allow-nonoverlapping-command-line-samples", "true");
        // Two missing names, to see the order and the separator of the message.
        run(dir, "two-missing-names", input, "-sn", "zeta", "-sn", "alpha");

        // Exclusion on its own, and exclusion beating inclusion.
        run(dir, "exclude-one", input, "-xl-sn", "s0");
        run(dir, "exclude-expression", input, "-xl-se", "^s");
        run(dir, "exclude-what-was-included", input, "-sn", "s0", "-xl-sn", "s0");
        run(dir, "exclude-some-of-what-was-included", input, "-sn", "s0", "-sn", "s1",
                "-xl-sn", "s0");
        // Excluding every sample, with nothing included.
        final List<String> excludeAll = new ArrayList<>();
        for (final String sample : SAMPLES) {
            excludeAll.add("-xl-sn");
            excludeAll.add(sample);
        }
        run(dir, "exclude-everything", input, excludeAll.toArray(new String[0]));
        // An exclusion naming a sample that is not there, which is not checked at all.
        run(dir, "exclude-missing-name", input, "-xl-sn", "ghost");

        // And an expression the regex engine refuses.
        run(dir, "uncompilable-expression", input, "-se", "[");
    }

    static void expressions(final String label, final String... patterns) {
        final Set<String> filters = new LinkedHashSet<>(List.of(patterns));
        try {
            final Set<String> matched = Utils.filterCollectionByExpressions(SAMPLES, filters, false);
            System.out.printf("expressions\t%s\t%s%n", label, String.join(",", matched));
        } catch (final Exception e) {
            System.out.printf("error\texpressions-%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
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
        for (final String line : lines) {
            if (line.startsWith("#CHROM")) {
                final String[] field = line.split("\t", -1);
                final List<String> samples = new ArrayList<>();
                for (int i = 9; i < field.length; i++) {
                    samples.add(field[i]);
                }
                System.out.printf("samples\t%s\t%s%n", label, String.join(",", samples));
                continue;
            }
            if (line.startsWith("#")) {
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
