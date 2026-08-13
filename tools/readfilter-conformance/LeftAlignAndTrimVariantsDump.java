/*
 * LeftAlignAndTrimVariants' output, taken from the reference.
 *
 * The tool the last three bricks were built for. Its own logic is short: choose a window, call
 * leftAlignAndTrim, and remember the record it just wrote.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE WINDOW IS BOUNDED BY THE PREVIOUS RECORD, `min(maxLeadingBases, distanceToLastVariant -
 *     1)`, measured from the previous record's END to this one's START;
 *   - AND THE PREVIOUS RECORD IS THE ONE AS WRITTEN, NOT AS READ. lastVariant is the ALIGNED
 *     record, which has moved LEFT, so aligning a variant RELAXES the bound on the next one: two
 *     indels a base apart in the input are eight apart by the time the second is measured, and the
 *     second moves after all;
 *   - A NEW CONTIG IS Integer.MAX_VALUE AWAY, so the bound only applies within a contig;
 *   - AN INDEL LONGER THAN --max-indel-length IS WRITTEN UNCHANGED, and it still becomes
 *     lastVariant, so it bounds the window of the record after it though it was never aligned;
 *   - AND THE LENGTH TESTED IS THE LARGEST ABSOLUTE indel length of the record, `getIndelLengths`
 *     mapped through abs and maxed, so a deletion of 300 and an insertion of 1 in one record are
 *     tested as 300;
 *   - --dont-trim-alleles PASSES `!dontTrimAlleles` DOWN, so the record moves without shrinking;
 *   - --split-multi-allelics SPLITS FIRST AND ALIGNS EACH PIECE, so one input record can leave
 *     several output records at several positions, and each piece bounds the next through the same
 *     lastVariant;
 *   - AND THE HEADER'S CONTIG LINES COME FROM THE REFERENCE, not from the input, so a contig the
 *     input never mentions is in the output header.
 *
 * Output:
 *
 *     reference\t<contig>\t<the bases the tool aligned against>
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: LeftAlignAndTrimVariantsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.LeftAlignAndTrimVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class LeftAlignAndTrimVariantsDump {

    /** The same shape as the leftAlignAndTrim suite: runs an indel can walk through. */
    static final String CHR1 =
            "GGGGGGGGGG"                       // 1..10
            + "AAAAAAAA"                       // 11..18
            + "GGGGGGGGGGGG"                   // 19..30
            + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" // 31..60
            + "GGGGGGGGGG"                     // 61..70
            + "CACACACACACACACACACA"           // 71..90
            + "GGGGGGGGGG";                    // 91..100

    /** A second contig the inputs never mention, so the header shows where its line came from. */
    static final String CHR2 = "TTTTTTTTTTAAAAAAAAAATTTTTTTTTT";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("leftalignandtrimvariants-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# LeftAlignAndTrimVariantsDump: the tool, from the reference");
        System.out.printf("reference\tchr1\t%s%n", CHR1);
        System.out.printf("reference\tchr2\t%s%n", CHR2);

        final Path fasta = writeReference(dir);

        // One deletion in each run, far apart, so nothing bounds anything.
        final Path apart = writeVcf(dir, "apart",
                "chr1\t17\t.\tAA\tA\t.\t.\t.",
                "chr1\t59\t.\tTT\tT\t.\t.\t.");
        // Two indels one base apart, so the second gets a window of zero.
        final Path adjacent = writeVcf(dir, "adjacent",
                "chr1\t17\t.\tAA\tA\t.\t.\t.",
                "chr1\t19\t.\tGG\tG\t.\t.\t.");
        // The same two with one base more between them.
        final Path nearly = writeVcf(dir, "nearly-adjacent",
                "chr1\t17\t.\tAA\tA\t.\t.\t.",
                "chr1\t20\t.\tGG\tG\t.\t.\t.");
        // A record on each contig, so the bound does not cross.
        final Path contigs = writeVcf(dir, "contigs",
                "chr1\t59\t.\tTT\tT\t.\t.\t.",
                "chr2\t19\t.\tAA\tA\t.\t.\t.");
        // A long deletion that CAN move, so skipping it is visible, and a short indel after it.
        final Path longIndel = writeVcf(dir, "long-indel",
                "chr1\t34\t.\t" + CHR1.substring(33, 59) + "\tT\t.\t.\t.",
                "chr1\t69\t.\tGG\tG\t.\t.\t.");
        // A multiallelic record whose two alternates walk to different places.
        final Path multi = writeVcf(dir, "multiallelic",
                "chr1\t17\t.\tAA\tA,AAA\t.\t.\t.",
                "chr1\t59\t.\tTT\tT\t.\t.\t.");
        // Alleles sharing a trailing base, so trimming is visible.
        final Path untrimmed = writeVcf(dir, "untrimmed",
                "chr1\t18\t.\tAGG\tAG\t.\t.\t.");

        run(dir, "apart", apart, fasta);
        run(dir, "adjacent", adjacent, fasta);
        run(dir, "nearly-adjacent", nearly, fasta);
        run(dir, "contigs", contigs, fasta);
        run(dir, "long-indel", longIndel, fasta);
        run(dir, "long-indel-allowed", longIndel, fasta, "--max-indel-length", "500");
        // Short enough to skip the long deletion, which is then written untouched AND still
        // bounds the window of the record after it.
        run(dir, "long-indel-skipped", longIndel, fasta, "--max-indel-length", "5");
        run(dir, "multiallelic", multi, fasta);
        run(dir, "multiallelic-split", multi, fasta, "--split-multi-allelics", "true");
        run(dir, "untrimmed", untrimmed, fasta);
        run(dir, "untrimmed-no-trim", untrimmed, fasta, "--dont-trim-alleles", "true");
        // A window narrow enough to stop the walk, which the record does not say.
        run(dir, "narrow-window", apart, fasta, "--max-leading-bases", "2");
        run(dir, "zero-window", apart, fasta, "--max-leading-bases", "0");
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final String... arguments) {
        final Path output = dir.resolve(label + "-out.vcf");
        final List<String> all = new ArrayList<>(List.of(
                "-V", input.toString(), "-R", fasta.toString(), "-O", output.toString(),
                "--suppress-reference-path", "true"));
        all.addAll(List.of(arguments));
        try {
            new LeftAlignAndTrimVariants().instanceMain(all.toArray(new String[0]));
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

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n>chr2\n" + CHR2 + "\n",
                StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict,
                "@HD\tVN:1.6\tSO:unsorted\n"
                        + "@SQ\tSN:chr1\tLN:" + CHR1.length() + "\tM5:0\tUR:file:" + fasta + "\n"
                        + "@SQ\tSN:chr2\tLN:" + CHR2.length() + "\tM5:0\tUR:file:" + fasta + "\n",
                StandardCharsets.UTF_8);
        return fasta;
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
