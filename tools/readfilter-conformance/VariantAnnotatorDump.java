/*
 * VariantAnnotator's copied annotations, taken from the reference.
 *
 * How one VCF is annotated from another. A resource file is tagged with a name, an expression names
 * one of its fields, and the value lands on any record at the same position. What is measurable is
 * which records are annotated, under what key, and what happens when the two files disagree about
 * the alleles.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE RESOURCE IS NAMED BY ITS TAG and the annotation's key is `<tag>.<field>`, so the same
 *     file under two tags is two annotations;
 *   - AN EXPRESSION CAN NAME AN INFO FIELD, `ID`, `ALT` OR `FILTER`, not just an INFO one;
 *   - A RECORD WITH NO RESOURCE AT ITS POSITION IS LEFT ALONE rather than annotated with nothing;
 *   - A PER-ALLELE FIELD CANNOT CROSS TO A DIFFERENT ALTERNATE and is withheld by default, while
 *     a SCALAR one and the ID, ALT and FILTER fields cross freely;
 *   - --resource-allele-concordance IS WHAT WITHHOLDS THE SCALAR ONE too;
 *   - --comparison ADDS ITS OWN ANNOTATION under the tag alone;
 *   - AN EXPRESSION NAMING A FIELD THE RESOURCE DOES NOT HAVE ADDS NOTHING, silently;
 *   - AN EXPRESSION WHOSE TAG NAMES NO RESOURCE IS REFUSED, which the unknown FIELD is not;
 *   - AND THE OUTPUT KEEPS EVERY ANNOTATION THE INPUT ALREADY CARRIED.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VariantAnnotatorDump
 */

import org.broadinstitute.hellbender.tools.walkers.annotator.VariantAnnotator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class VariantAnnotatorDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FILTER=<ID=LOW,Description=\"Low\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "##INFO=<ID=NOTE,Number=1,Type=String,Description=\"A note\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
    }

    static String site(final int position, final String id, final String reference,
                       final String alternate, final String filter, final String info) {
        return "chr1\t" + position + "\t" + id + "\t" + reference + "\t" + alternate
                + "\t100.00\t" + filter + "\t" + info + "\tGT\t0/1";
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("variant-annotator-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VariantAnnotatorDump: how one VCF is annotated from another");

        final Path fasta = writeReference(dir);

        // The input: four sites, one of which the resource does not mention.
        final String input = vcf(List.of(
                site(1000, ".", "A", "C", "PASS", "NOTE=kept"),
                site(2000, ".", "G", "T", "PASS", "."),
                site(3000, ".", "C", "A", "PASS", "."),
                site(4000, ".", "T", "G", "PASS", ".")));
        final Path inputPath = write(dir, "input.vcf", input);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(input));

        // The resource: the same first three positions, one of them with a DIFFERENT alternate.
        // Each record carries a PER-ALLELE field (AC, Number=A) and a SCALAR one (NOTE,
        // Number=1), so a value that cannot be mapped to a different alternate is told apart from
        // one that can.
        final String resource = vcf(List.of(
                site(1000, "rs1", "A", "C", "PASS", "AC=5;AF=0.25;NOTE=first"),
                site(2000, "rs2", "G", "A", "LOW", "AC=7;AF=0.35;NOTE=second"),
                site(3000, ".", "C", "A", "PASS", "AC=9;AF=0.45;NOTE=third")));
        final Path resourcePath = write(dir, "resource.vcf", resource);
        // The resource is QUERIED by interval, so it needs an index beside it.
        htsjdk.tribble.index.IndexFactory.createLinearIndex(resourcePath.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(resourcePath.toFile());
        System.out.printf("vcf\tresource=%s%n", ReferenceQueryDump.escape(resource));

        run(dir, "one-expression", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.AC"));
        // The scalar field, whose value does not depend on which alternate it belongs to.
        run(dir, "scalar-expression", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.NOTE"));
        run(dir, "scalar-concordance", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.NOTE",
                        "--resource-allele-concordance", "true"));
        // Two expressions from one resource, and the same file under a second tag.
        run(dir, "two-expressions", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.AC", "-E", "res.AF"));
        run(dir, "two-tags", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(),
                        "--resource:other", resourcePath.toString(),
                        "-E", "res.AC", "-E", "other.AC"));
        // The three non-INFO fields an expression may name.
        run(dir, "id-alt-filter", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(),
                        "-E", "res.ID", "-E", "res.ALT", "-E", "res.FILTER"));
        // Allele concordance, which withholds the value where the alternates differ.
        run(dir, "allele-concordance", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.AC",
                        "--resource-allele-concordance", "true"));
        // A field the resource does not carry, and a tag that names no resource.
        run(dir, "unknown-field", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "res.MISSING"));
        run(dir, "unknown-tag", inputPath, fasta,
                List.of("--resource:res", resourcePath.toString(), "-E", "nothere.AC"));
        // A comparison file, which annotates under the tag alone.
        run(dir, "comparison", inputPath, fasta,
                List.of("--comp:cmp", resourcePath.toString()));
        // No resource at all, which copies the input through.
        run(dir, "no-resource", inputPath, fasta, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new VariantAnnotator().instanceMain(argv.toArray(new String[0]));
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
            return;
        }
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(out).split("\n", -1)) {
            if (!line.startsWith("##") && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(body.toString(), dir)));
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
