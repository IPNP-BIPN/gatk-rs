/*
 * CombineGVCFs' merged bands, taken from the reference.
 *
 * How several single-sample GVCFs become one multi-sample one. Every sample's reference blocks are
 * cut wherever any other sample has an edge, so the output's records are the union of every
 * input's boundaries rather than any one input's.
 *
 * Nine behaviours this is built to catch.
 *
 *   - EVERY SAMPLE'S BOUNDARIES CUT EVERY OTHER SAMPLE'S BLOCKS, so a block that no input broke is
 *     still broken in the output;
 *   - A VARIANT IN ONE SAMPLE CUTS THE OTHERS' BLOCKS at its own position, and the samples that
 *     have no variant there are written as reference at that base;
 *   - THE MERGED BLOCK'S QUALITY IS EACH SAMPLE'S OWN, so the output carries one genotype per
 *     sample rather than a summary;
 *   - --convert-to-base-pair-resolution BREAKS EVERY BLOCK INTO SINGLE BASES;
 *   - --break-bands-at-multiples-of BREAKS THEM AT A GRID instead, whatever the data says;
 *   - THE TWO ARE NOT MUTUALLY EXCLUSIVE: given together, base-pair resolution wins and the grid
 *     is ignored, which the run below shows by producing byte for byte what base-pair resolution
 *     produced on its own;
 *   - <NON_REF> SURVIVES ON EVERY RECORD, and the alternates of a variant are the union of the
 *     samples that carried one;
 *   - A SAMPLE THAT ENDS EARLY KEEPS ITS COLUMN AND LOSES ITS FIELDS: it is written `./.` with
 *     nothing after it, rather than being padded with reference or dropped;
 *   - NO GENOTYPE IS CALLED: every sample is `./.` on every record, and --call-genotypes is what
 *     changes that;
 *   - AND THE SAME SAMPLE NAME IN TWO INPUTS IS REFUSED.
 *
 * Output:
 *
 *     vcf\t<label>=<that gvcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CombineGVCFsDump
 */

import org.broadinstitute.hellbender.tools.walkers.CombineGVCFs;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CombineGVCFsDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header(final String sample) {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##ALT=<ID=NON_REF,Description=\"Any other allele\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth\">",
                "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample));
    }

    static String block(final int start, final int end, final int quality) {
        return "chr1\t" + start + "\t.\tA\t<NON_REF>\t.\t.\tEND=" + end
                + "\tGT:DP:GQ:MIN_DP:PL\t0/0:20:" + quality + ":20:0," + quality + ","
                + (quality * 10);
    }

    static String variant(final int position, final String alternate, final int quality) {
        return "chr1\t" + position + "\t.\tA\t" + alternate + ",<NON_REF>\t50.00\t.\t."
                + "\tGT:AD:DP:GQ:PL\t0/1:8,4,0:12:" + quality + ":" + (quality * 10) + ",0,"
                + (quality * 10) + "," + (quality * 10) + ",0," + (quality * 10);
    }

    static Path writeGvcf(final Path dir, final String sample, final List<String> records)
            throws Exception {
        final List<String> lines = header(sample);
        lines.addAll(records);
        lines.add("");
        final String text = String.join("\n", lines);
        System.out.printf("vcf\t%s=%s%n", sample, ReferenceQueryDump.escape(text));
        final Path path = write(dir, sample + ".g.vcf", text);
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("combine-gvcfs-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CombineGVCFsDump: how several single-sample GVCFs become one");

        final Path fasta = writeReference(dir);

        // s1: one long block, one variant, one more block.
        final Path s1 = writeGvcf(dir, "s1", List.of(
                block(1000, 1199, 40),
                variant(1200, "C", 60),
                block(1201, 1400, 40)));
        // s2: blocks whose edges fall INSIDE s1's, and a variant elsewhere.
        final Path s2 = writeGvcf(dir, "s2", List.of(
                block(1000, 1099, 30),
                block(1100, 1299, 50),
                variant(1300, "G", 70),
                block(1301, 1400, 30)));
        // s3: a sample that stops early, so it is absent from the later records.
        final Path s3 = writeGvcf(dir, "s3", List.of(
                block(1000, 1150, 20)));

        run(dir, "three-samples", List.of(s1, s2, s3), fasta, List.of());
        run(dir, "two-samples", List.of(s1, s2), fasta, List.of());
        // Single bases everywhere.
        run(dir, "base-pair-resolution", List.of(s1, s2, s3), fasta,
                List.of("--convert-to-base-pair-resolution", "true"));
        // A grid, whatever the data says.
        run(dir, "break-bands-100", List.of(s1, s2, s3), fasta,
                List.of("--break-bands-at-multiples-of", "100"));
        run(dir, "break-bands-50", List.of(s1, s2, s3), fasta,
                List.of("--break-bands-at-multiples-of", "50"));
        // The two together.
        run(dir, "both-band-arguments", List.of(s1, s2, s3), fasta,
                List.of("--convert-to-base-pair-resolution", "true",
                        "--break-bands-at-multiples-of", "100"));
        // Genotypes called rather than left as no-calls.
        run(dir, "call-genotypes", List.of(s1, s2, s3), fasta,
                List.of("--call-genotypes", "true"));
        // One input on its own, which is the control.
        run(dir, "one-sample", List.of(s1), fasta, List.of());
        // The same sample name twice.
        run(dir, "duplicate-sample", List.of(s1, s1), fasta, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final List<Path> inputs, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".g.vcf");
        final List<String> argv = new ArrayList<>();
        for (final Path input : inputs) {
            argv.add("-V");
            argv.add(input.toString());
        }
        argv.addAll(List.of("-O", out.toString(), "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new CombineGVCFs().instanceMain(argv.toArray(new String[0]));
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
