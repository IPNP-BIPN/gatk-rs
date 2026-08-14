/*
 * AnnotateVcfWithExpectedAlleleFraction's output VCF, taken from the reference.
 *
 * A VariantWalker that writes its input back out with one Float INFO field, the dot product of each
 * sample's genotype weight with its mixing fraction.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE DEFAULT TOOL HEADER LINES NEVER REACH THE FILE:
 *
 *         final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
 *         headerLines.addAll(getDefaultToolVCFHeaderLines());
 *         vcfWriter.writeHeader(vcfHeader);
 *
 *     the set is added to AFTER the header is built from it, so `##source=` and `##GATKCommandLine`
 *     are written by the sibling tool AnnotateVcfWithBamDepth and not by this one. The golden is
 *     what says whether the constructor copies the set or keeps a reference to it;
 *   - THE WEIGHT IS 1.0, 0.5 OR ZERO and nothing else: `isHomVar` then `isHet`, so a no-call, a
 *     half-call and a hom ref all weigh zero, and a 1/2 call weighs 0.5 like any other het;
 *   - THE MIXING FRACTIONS ARE READ IN THE HEADER'S SAMPLE ORDER, `getSampleNamesInOrder()`, which
 *     is SORTED rather than the order the columns appear in: a table read into the wrong order
 *     would still produce a number, just the wrong one;
 *   - A SAMPLE MISSING FROM THE TABLE IS A NULL UNBOXING, `mapToDouble(map::get)` on a null, so the
 *     refusal is a NullPointerException rather than anything the tool words itself;
 *   - A SAMPLE LISTED TWICE IN THE TABLE IS AN IllegalStateException out of `Collectors.toMap`;
 *   - THE FRACTIONS ARE NOT REQUIRED TO SUM TO ONE, and nothing checks it;
 *   - AND AF_EXP IS WRITTEN AS A FLOAT, so the value goes through htsjdk's own formatter and 0.15
 *     is not necessarily what the arithmetic produced.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     fractions\t<label>\t<the whole mixing fractions table, escaped>
 *     vcfline\t<label>\t<one line of the output vcf, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AnnotateVcfWithExpectedAlleleFractionDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithExpectedAlleleFraction;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AnnotateVcfWithExpectedAlleleFractionDump {

    /** The columns are declared in the header's order, which is not the sorted one. */
    static final String SAMPLES = "zebra\talpha\tmike";

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"sorts before AF_EXP\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=200>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + SAMPLES + "\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("annotatevcfwithexpectedallelefraction-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AnnotateVcfWithExpectedAlleleFractionDump: a dot product, and a header built too early");

        // One record per genotype shape, and one that carries an INFO field already.
        final Path variants = writeVcf(dir, "variants",
                // zebra het, the others hom ref: 0.5 * zebra's fraction.
                "chr1\t20\t.\tA\tC\t50\tPASS\tAA=1\tGT\t0/1\t0/0\t0/0",
                // zebra hom var, alpha het: 1.0 * zebra + 0.5 * alpha.
                "chr1\t40\t.\tA\tC\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0",
                // Every sample hom var, so the answer is the sum of the fractions.
                "chr1\t60\t.\tA\tC\t50\tPASS\t.\tGT\t1/1\t1/1\t1/1",
                // A no-call, a half-call and a hom ref, all of which weigh nothing.
                "chr1\t80\t.\tA\tC\t50\tPASS\t.\tGT\t./.\t./1\t0/0",
                // A multi-allelic 1/2, which is a het like any other.
                "chr1\t100\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2\t0/0\t0/0");

        // Fractions that do not sum to one, and whose order is not the header's.
        final Path fractions = writeFractions(dir, "fractions",
                "alpha\t0.2", "mike\t0.1", "zebra\t0.3");
        // The same three summing to one, to show nothing checks it either way.
        final Path normalized = writeFractions(dir, "normalized",
                "zebra\t0.5", "alpha\t0.3", "mike\t0.2");
        // A table missing one of the VCF's samples.
        final Path missing = writeFractions(dir, "missing-sample",
                "zebra\t0.5", "alpha\t0.5");
        // A table naming one sample twice.
        final Path duplicated = writeFractions(dir, "duplicate-sample",
                "zebra\t0.5", "zebra\t0.4", "alpha\t0.05", "mike\t0.05");

        run(dir, "annotated", variants, fractions, "annotated.vcf");
        run(dir, "normalized", variants, normalized, "normalized.vcf");
        run(dir, "missing-sample", variants, missing, "missing-sample.vcf");
        run(dir, "duplicate-sample", variants, duplicated, "duplicate-sample.vcf");
        run(dir, "output-is-a-directory", variants, fractions, ".");
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

    /** The table `MixingFraction.readMixingFractions` reads, which is the writer's own format. */
    static Path writeFractions(final Path dir, final String label, final String... rows)
            throws Exception {
        final StringBuilder text = new StringBuilder("SAMPLE\tMIXING_FRACTION\n");
        for (final String row : rows) {
            text.append(row).append("\n");
        }
        final Path file = dir.resolve(label + ".table");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        System.out.printf("fractions\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path input, final Path fractions,
                    final String output, final String... arguments) throws Exception {
        final Path file = dir.resolve(output);
        final List<String> all = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", file.toString(),
                "--mixing-fractions", fractions.toString()));
        all.addAll(List.of(arguments));
        try {
            new AnnotateVcfWithExpectedAlleleFraction().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, file);
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
