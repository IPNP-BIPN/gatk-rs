/*
 * A GenomicsDB workspace read back: GenotypeGVCFs, GnarlyGenotyper and SelectVariants over
 * `gendb://`, through Main.
 *
 * The reader is native. `GenomicsDBFeatureReader` asks GenomicsDB's `GenomicsDBBCFGenerator` for
 * VCF text, which is `VariantQueryProcessor::scan_and_operate` cutting the cells into gVCF
 * intervals and `BroadCombinedGVCFOperator` merging the calls of each one, printed by htslib. What
 * that text does that a CombineGVCFs output does not is the reason for the cases:
 *
 *   - A SPANNING DELETION'S ALLELE COMES FIRST at the site under it (`*,T,<NON_REF>`), its INFO
 *     dropped and its likelihoods reduced to `*` and `<NON_REF>`;
 *   - IDS ARE THE UNION OF THE CALLS that overlap a position, in libstdc++'s hash order;
 *   - THE ALLELE-SPECIFIC RAW ANNOTATIONS ARE COMBINED, element-wise or as histograms;
 *   - FORMAT DP IS LAST, and a sample with no call has `.` in every field, GT included;
 *   - A QUERY CLIPS THE BLOCKS TO ITS INTERVAL, so `-L` is not a filter over the whole stream;
 *   - AND THE CALLERS FORCE --call-genotypes, which SelectVariants takes from its own argument.
 *
 * The inputs are the shared corpus's `cg_a`, `cg_b` and `cg_c` and their reference, imported into
 * a workspace of the run's own. The import's command line, which the workspace's header carries
 * into every output, records its `Date`, which is masked.
 *
 * Output:
 *
 *     input\t<name>\t<base64 of a fixture file>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     out\t<label>\t<the output VCF, Date and cwd masked, escaped>
 *     error\t<label>\t<exception class>:<message, escaped, cwd masked>
 *
 * Usage: GenomicsDBQueryDump (with the shared fixtures at /work/fixtures)
 */

import org.broadinstitute.hellbender.Main;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class GenomicsDBQueryDump {

    static final Path FIXTURES = Path.of("/work/fixtures");
    static final Path DIR = Path.of("genomicsdb-query-dump");
    static final String WS = "gendb://" + DIR.resolve("ws");
    static final String REF = DIR.resolve("cgv_ref.fasta").toString();

    public static void main(final String[] args) throws Exception {
        Files.createDirectories(DIR.resolve("out"));
        for (final String name : List.of("cg_a.g.vcf", "cg_b.g.vcf", "cg_c.g.vcf", "cgv_ref.fasta",
                "cgv_ref.fasta.fai", "cgv_ref.dict")) {
            final byte[] bytes = Files.readAllBytes(FIXTURES.resolve(name));
            Files.write(DIR.resolve(name), bytes);
            System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(bytes));
        }
        for (final String name : List.of("cg_a", "cg_b", "cg_c")) {
            new Main().instanceMain(new String[] {"IndexFeatureFile", "-I", DIR.resolve(name + ".g.vcf").toString()});
        }
        run("import", null, "GenomicsDBImport", "-V", DIR.resolve("cg_a.g.vcf").toString(),
                "-V", DIR.resolve("cg_b.g.vcf").toString(), "-V", DIR.resolve("cg_c.g.vcf").toString(),
                "-L", "chr1", "-L", "chr2", "--genomicsdb-workspace-path", DIR.resolve("ws").toString());

        select("select");
        select("select-intervals", "-L", "chr1:120-215");
        select("select-two-contigs", "-L", "chr1:95-105", "-L", "chr2:15-25");
        select("select-call-genotypes", "--call-genotypes");
        select("select-shared-posixfs", "--genomicsdb-shared-posixfs-optimizations");
        select("select-sample", "--sample-name", "sB", "--exclude-non-variants");
        select("select-max-alternates", "--genomicsdb-max-alternate-alleles", "7");
        select("select-too-few-alternates", "--genomicsdb-max-alternate-alleles", "6");
        genotype("genotype");
        genotype("genotype-intervals", "-L", "chr1:120-215");
        genotype("genotype-allele-specific", "-G", "StandardAnnotation", "-G", "AS_StandardAnnotation");
        genotype("genotype-all-sites", "--include-non-variant-sites", "-L", "chr1:95-145");
        genotype("genotype-merged-intervals", "-L", "chr1:95-105", "-L", "chr1:205-215", "--merge-input-intervals");
        genotype("genotype-too-few-alternates", "--genomicsdb-max-alternate-alleles", "5");
        gnarly("gnarly");
        gnarly("gnarly-intervals", "-L", "chr1:95-145");
        gnarly("gnarly-keep-all-sites", "--keep-all-sites");

        run("no-reference", "x", "SelectVariants", "-V", WS, "-O", out("no-reference"),
                "--add-output-vcf-command-line", "false");
        run("no-workspace", "x", "SelectVariants", "-R", REF, "-V", "gendb://" + DIR.resolve("none"),
                "-O", out("no-workspace"), "--add-output-vcf-command-line", "false");

        // The importer leaves threads behind that instanceMain, unlike main, does not exit past.
        System.out.flush();
        System.exit(0);
    }

    static String out(final String label) {
        return DIR.resolve("out").resolve(label + ".vcf").toString();
    }

    static void select(final String label, final String... rest) throws Exception {
        caller(label, "SelectVariants", rest);
    }

    static void genotype(final String label, final String... rest) throws Exception {
        caller(label, "GenotypeGVCFs", rest);
    }

    static void gnarly(final String label, final String... rest) throws Exception {
        caller(label, "GnarlyGenotyper", rest);
    }

    static void caller(final String label, final String tool, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(List.of(tool, "-R", REF, "-V", WS, "-O", out(label),
                "--add-output-vcf-command-line", "false"));
        argv.addAll(Arrays.asList(rest));
        run(label, out(label), argv.toArray(new String[0]));
    }

    static void run(final String label, final String output, final String... argv) throws Exception {
        System.out.printf("argv\t%s\t%s%n", label, ReferenceQueryDump.escape(String.join("\t", argv)));
        String failure = null;
        try {
            new Main().instanceMain(argv);
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            failure = e.getClass().getName() + ":" + String.valueOf(cause.getMessage());
        }
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(masked(failure)));
            return;
        }
        if (output != null && Files.exists(Path.of(output))) {
            final String text = Files.readString(Path.of(output), StandardCharsets.UTF_8);
            System.out.printf("out\t%s\t%s%n", label, ReferenceQueryDump.escape(masked(text)));
        }
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath() + "/", "<cwd>/")
                .replaceAll("Date=\"[^\"]*\"", "Date=\"<masked>\"");
    }
}
