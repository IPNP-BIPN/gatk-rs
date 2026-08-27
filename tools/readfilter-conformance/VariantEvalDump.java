/*
 * VariantEval's report, taken from the reference.
 *
 * How a call set is counted against a comparison one. The tool writes a GATKReport of one table per
 * evaluation module, each stratified by whatever stratifiers were asked for, and what a row counts
 * depends on the stratification as much as on the data.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE REPORT IS ONE TABLE PER MODULE, and --list names the modules without reading a record;
 *   - THE STANDARD STRATIFIERS APPLY UNLESS TURNED OFF, so every table has an `all` row and a row
 *     per novelty even when none were asked for;
 *   - NOVELTY IS DECIDED BY THE dbSNP TRACK AND NOT BY --comp: the same file given as `--comp`
 *     leaves every site novel and given as `--dbsnp` splits them, which is the whole of the
 *     difference between the two arguments here;
 *   - A STRATIFIER MULTIPLIES THE ROWS rather than adding a column;
 *   - `CountVariants` COUNTS BY TYPE and `TiTvVariantEvaluator` BY SUBSTITUTION, so the same file
 *     is summarised twice under different questions;
 *   - AN INDEL'S LENGTH IS SIGNED, insertions positive and deletions negative;
 *   - A MULTIALLELIC SITE IS ITS OWN CATEGORY rather than being counted once per alternate;
 *   - --select-expression ADDS A NAMED SUBSET as another stratum;
 *   - AND AN UNKNOWN MODULE NAME IS REFUSED, naming what it could not find.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     out\t<label>=<the whole report, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VariantEvalDump
 */

import org.broadinstitute.hellbender.tools.walkers.varianteval.VariantEval;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VariantEvalDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
    }

    static String site(final int position, final String reference, final String alternate,
                       final String genotype) {
        return "chr1\t" + position + "\t.\t" + reference + "\t" + alternate
                + "\t100.00\tPASS\t.\tGT\t" + genotype;
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variant-eval-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VariantEvalDump: how a call set is counted against a comparison one");

        final Path fasta = writeReference(dir);

        // The eval set: transitions, transversions, an insertion, a deletion, a multiallelic.
        final String eval = vcf(List.of(
                site(1000, "A", "G", "0/1"),          // transition
                site(2000, "C", "T", "0/1"),          // transition
                site(3000, "A", "C", "0/1"),          // transversion
                site(4000, "G", "T", "1/1"),          // transversion
                site(5000, "A", "ACGT", "0/1"),       // insertion of three
                site(6000, "ACGT", "A", "0/1"),       // deletion of three
                site(7000, "A", "C,G", "1/2")));      // multiallelic
        final Path evalPath = writeIndexed(dir, "eval.vcf", eval, "eval");

        // The comparison set: three of the seven positions, so the rest are novel.
        final String comp = vcf(List.of(
                site(1000, "A", "G", "0/1"),
                site(3000, "A", "C", "0/1"),
                site(5000, "A", "ACGT", "0/1")));
        final Path compPath = writeIndexed(dir, "comp.vcf", comp, "comp");

        run(dir, "no-comp", fasta, List.of("--eval", evalPath.toString()));
        // The SAME comparison file given as dbSNP, which is what the novelty stratifier reads.
        run(dir, "dbsnp", fasta, List.of("--eval", evalPath.toString(),
                "--dbsnp", compPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "CountVariants"));
        run(dir, "with-comp", fasta,
                List.of("--eval", evalPath.toString(), "--comp", compPath.toString()));
        // One module at a time, with the standard stratifiers turned off.
        run(dir, "count-variants", fasta, List.of("--eval", evalPath.toString(),
                "--comp", compPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "CountVariants"));
        run(dir, "titv", fasta, List.of("--eval", evalPath.toString(),
                "--comp", compPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "TiTvVariantEvaluator"));
        run(dir, "indel-length", fasta, List.of("--eval", evalPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "IndelLengthHistogram"));
        run(dir, "multiallelic", fasta, List.of("--eval", evalPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "MultiallelicSummary"));
        // The standard stratifiers off, which removes the novelty rows.
        run(dir, "no-standard-stratifiers", fasta, List.of("--eval", evalPath.toString(),
                "--comp", compPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "CountVariants",
                "--do-not-use-all-standard-stratifications", "true"));
        // A stratifier that multiplies the rows.
        run(dir, "stratify-by-type", fasta, List.of("--eval", evalPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "CountVariants",
                "-ST", "VariantType"));
        // A named subset.
        run(dir, "select-expression", fasta, List.of("--eval", evalPath.toString(),
                "--do-not-use-all-standard-modules", "true", "-EV", "CountVariants",
                "-select", "QUAL > 50", "-select-name", "highqual"));
        // An unknown module, and an unknown stratifier.
        run(dir, "unknown-module", fasta, List.of("--eval", evalPath.toString(),
                "-EV", "NoSuchEvaluator"));
        run(dir, "unknown-stratifier", fasta, List.of("--eval", evalPath.toString(),
                "-ST", "NoSuchStratifier"));
        // --list is LAST because it ENDS THE PROCESS: it prints the module names and exits the
        // JVM. The marker below is written BEFORE the call and the one after it never is, which is
        // how the golden records that.
        System.out.println("none\tlist=about to run, which ends the process");
        run(dir, "list", fasta, List.of("--eval", evalPath.toString(), "--list", "true"));
        System.out.println("none\tafter-list=this line is never reached");
    }

    static Path writeIndexed(final Path dir, final String name, final String text,
                             final String label) throws Exception {
        System.out.printf("vcf\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        final Path path = write(dir, name, text);
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".txt");
        final List<String> argv = new ArrayList<>(List.of(
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new VariantEval().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            // Some runs return without writing anything and without throwing: --list prints the
            // modules and stops, and an unrecognised module name does the same. Reported rather
            // than passed over in silence.
            System.out.printf("none\t%s=no output file%n", label);
            return;
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(List.of(
                new htsjdk.samtools.SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
