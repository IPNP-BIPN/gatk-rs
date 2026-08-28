/*
 * AlleleFrequencyQC's verdict, taken from the reference.
 *
 * The tool is VariantEval with every knob preset: one module, two stratifiers, and a logarithmic
 * allele-frequency scale. What it adds is a chi-squared statistic over the AF bins and the
 * p-value that statistic gives, which is the whole of the tool's own arithmetic.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE COMPARISON SET IS GIVEN TWICE, once as `--comp` and once as a second eval track, which
 *     is what puts two rows in each bin: with one eval track every bin holds a single row and the
 *     statistic is nought whatever the data says;
 *   - THE STATISTIC IS NOT PEARSON'S: the expected count in the denominator is replaced by a
 *     CONSTANT variance, squared, so the same difference costs the same wherever it happens;
 *   - THE VARIANCE IS SQUARED AND DIVIDES THE WHOLE SUM, so raising --allowed-variance from a
 *     hundredth to a tenth divides the statistic by a HUNDRED and not by ten;
 *   - THE BIN LADDER IS FIXED AND NOT THE DATA'S: the logarithmic scale emits SIXTY-ONE bins
 *     whatever the file holds, each with one row per eval track, so the degrees of freedom are
 *     sixty in every run here and a bin no variant reached contributes a term of nought;
 *   - WHICH MAKES THE `fewer than two entries` GUARD UNREACHABLE ON THIS PATH: every bin has
 *     exactly the two rows, so the guard is carried in the port and never fires against the
 *     reference;
 *   - A COMPARISON SITE THE CALL SET HAS NOTHING AT STILL CONTRIBUTES, its bin holding the
 *     comparison frequency against a nought, which moves the statistic by that square alone;
 *   - THE P-VALUE IS THE UPPER TAIL, so a perfect match is one and a large statistic is nought;
 *   - THE ROWS ARE FILTERED TO `called` BEFORE THE GROUPING, so the filtered variants' own rows
 *     never reach the statistic;
 *   - THE SAMPLE NAME COMES FROM A `##sampleAlias` HEADER LINE and not from the genotype columns,
 *     so a VCF without one fails on a null rather than falling back;
 *   - --debug-file KEEPS THE VARIANT EVAL REPORT the statistic was read out of, which is otherwise
 *     written to a temporary file and deleted;
 *   - THE METRICS ARE WRITTEN BEFORE THE PLOT IS ATTEMPTED, so a run whose R script fails has
 *     already answered, which is why this dump reads the file rather than the exit status;
 *   - AND A FILE WITH ONE VARIANT IN IT IS NOT A DEGENERATE CASE: the ladder is the same
 *     sixty-one bins, the two tracks agree in the one bin they reach, and the statistic is nought
 *     against a p-value of one.
 *
 * Output:
 *
 *     vcf\t<label>=<that input file, escaped>
 *     metrics\t<label>\t<the metrics table without its comments, escaped>
 *     report\t<label>\t<the debug file's own table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AlleleFrequencyQCDump
 */

import org.broadinstitute.hellbender.tools.walkers.varianteval.AlleleFrequencyQC;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AlleleFrequencyQCDump {

    /** The eval set carries genotypes and the alias; the comparison set is sites only. */
    static String evalVcf(final List<String> sites, final boolean withAlias) {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + VariantEvalDump.CONTIG_LENGTH + ">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">"));
        if (withAlias) {
            lines.add("##sampleAlias=NA12878");
        }
        lines.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1");
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    static String compVcf(final List<String> sites) {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + VariantEvalDump.CONTIG_LENGTH + ">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    /** One eval site, whose genotype is the allele fraction the evaluator sums. */
    static String called(final int position, final String genotype) {
        return "chr1\t" + position + "\t.\tA\tG\t100.00\tPASS\t.\tGT\t" + genotype;
    }

    static String filtered(final int position, final String genotype) {
        return "chr1\t" + position + "\t.\tA\tG\t100.00\tLowQual\t.\tGT\t" + genotype;
    }

    /** One comparison site, whose AF is the bin the pair is stratified into. */
    static String site(final int position, final String frequency) {
        return "chr1\t" + position + "\t.\tA\tG\t100.00\tPASS\tAF=" + frequency;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("allele-frequency-qc-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AlleleFrequencyQCDump: the chi-squared statistic over the allele "
                + "frequency bins, and the p-value it gives");

        final Path fasta = VariantEvalDump.writeReference(dir);

        // Four allele frequencies an order of magnitude apart, so the logarithmic scale puts each
        // in a bin of its own, and the eval genotypes that match them.
        final List<String> frequencies = List.of("0.005", "0.05", "0.25", "0.5");
        final List<String> compSites = new ArrayList<>();
        final List<String> matchingSites = new ArrayList<>();
        final List<String> divergingSites = new ArrayList<>();
        for (int i = 0; i < frequencies.size(); i++) {
            final int position = 1000 * (i + 1);
            compSites.add(site(position, frequencies.get(i)));
            matchingSites.add(called(position, "0/1"));
            // The last bin's call disagrees: a homozygous variant where the others are het.
            divergingSites.add(called(position, i == frequencies.size() - 1 ? "1/1" : "0/1"));
        }
        final Path comp = VariantEvalDump.writeIndexed(dir, "comp.vcf", compVcf(compSites), "comp");
        final Path matching = VariantEvalDump.writeIndexed(dir, "matching.vcf",
                evalVcf(matchingSites, true), "matching");
        final Path diverging = VariantEvalDump.writeIndexed(dir, "diverging.vcf",
                evalVcf(divergingSites, true), "diverging");

        run(dir, "het-calls", fasta, matching, comp, List.of());
        run(dir, "one-bin-homozygous", fasta, diverging, comp, List.of());
        // The variance, which is squared and divides the whole sum.
        run(dir, "variance-tenth", fasta, diverging, comp,
                List.of("-allowed-variance", "0.1"));
        run(dir, "variance-hundredth", fasta, diverging, comp,
                List.of("-allowed-variance", "0.001"));
        // The threshold, which decides whether the tool complains and nothing else.
        run(dir, "threshold-above-the-pvalue", fasta, matching, comp,
                List.of("-pvalue-threshold", "0.99"));

        // A filtered variant, whose rows never reach the statistic.
        final List<String> withFiltered = new ArrayList<>(divergingSites);
        withFiltered.add(filtered(5000, "1/1"));
        final Path filteredEval = VariantEvalDump.writeIndexed(dir, "filtered.vcf",
                evalVcf(withFiltered, true), "filtered");
        run(dir, "a-filtered-variant", fasta, filteredEval, comp, List.of());

        // A comparison site the call set has nothing at, whose bin has one entry and not two.
        final List<String> extraComp = new ArrayList<>(compSites);
        extraComp.add(site(9000, "0.001"));
        final Path lonely = VariantEvalDump.writeIndexed(dir, "lonely.vcf", compVcf(extraComp),
                "comp-with-an-extra-bin");
        run(dir, "a-bin-with-one-entry", fasta, diverging, lonely, List.of());

        // A single variant, whose bin is the only one either track reaches. The ladder is still
        // sixty-one bins wide, so the statistic is nought rather than undefined.
        final Path oneComp = VariantEvalDump.writeIndexed(dir, "one-comp.vcf",
                compVcf(List.of(site(1000, "0.5"))), "comp-one-bin");
        final Path oneEval = VariantEvalDump.writeIndexed(dir, "one-eval.vcf",
                evalVcf(List.of(called(1000, "0/1")), true), "eval-one-bin");
        run(dir, "one-variant", fasta, oneEval, oneComp, List.of());

        // A VCF with no alias header at all.
        final Path noAlias = VariantEvalDump.writeIndexed(dir, "no-alias.vcf",
                evalVcf(matchingSites, false), "no-alias");
        run(dir, "no-sample-alias", fasta, noAlias, comp, List.of());
    }

    /** One run, with the metrics it wrote and the report it was read out of. */
    static void run(final Path dir, final String label, final Path fasta, final Path eval,
                    final Path comp, final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".txt");
        final Path debug = dir.resolve("debug-" + label + ".txt");
        final List<String> argv = new ArrayList<>(List.of(
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--eval", eval.toString(),
                "--comp", comp.toString(),
                // The comparison set is given AGAIN as a second, tagged eval track. That is what
                // puts two rows in each allele-frequency bin, one from the call set's genotypes
                // and one from the comparison set's own AF field, and the statistic is the
                // difference between them: with one eval track every bin holds a single row and
                // the sum is nought whatever the data says.
                "-eval:thousand_genomes", comp.toString(),
                "-L", comp.toString(),
                "-debug-file", debug.toString()));
        argv.addAll(extra);
        try {
            new AlleleFrequencyQC().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            // The plot is the LAST thing the tool does, and it wants R packages no oracle image
            // here carries. The metrics and the report are both written before it, so a run that
            // got that far is measured and the plot's own fate is left out: the golden then says
            // the same thing whether or not the image has R, which is what keeps it a measurement
            // of the tool rather than of the container.
            if (!Files.exists(out)) {
                System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
                return;
            }
        }
        System.out.printf("metrics\t%s\t%s%n", label,
                ReferenceQueryDump.escape(withoutComments(Files.readString(out,
                        StandardCharsets.UTF_8), dir)));
        System.out.printf("report\t%s\t%s%n", label,
                ReferenceQueryDump.escape(withoutComments(Files.readString(debug,
                        StandardCharsets.UTF_8), dir)));
    }

    /** The comment lines carry the command line and the run's own clock. */
    static String withoutComments(final String text, final Path dir) {
        final List<String> kept = new ArrayList<>();
        for (final String line : text.split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                kept.add(masked(line, dir));
            }
        }
        return String.join("\n", kept);
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
