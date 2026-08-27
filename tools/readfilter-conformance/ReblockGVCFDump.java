/*
 * ReblockGVCF's bands, taken from the reference.
 *
 * How a GVCF's reference blocks are coarsened and its weak variants demoted. The tool rewrites both
 * the blocks and the variants, and which of them survive is decided by a handful of thresholds
 * rather than by the data alone.
 *
 * Ten behaviours this is built to catch.
 *
 *   - ADJACENT REFERENCE BLOCKS IN THE SAME GQ BAND ARE MERGED, and the merged block carries the
 *     LOWEST genotype quality of the blocks that made it;
 *   - --gqb DEFINES THE BANDS' EDGES, so a different set of bounds merges a different set of
 *     blocks;
 *   - --rgq-threshold DEMOTES A VARIANT TO A REFERENCE BLOCK when its reference likelihood is
 *     below the threshold, rather than dropping it;
 *   - --drop-low-quals AND --rgq-threshold TOUCH DIFFERENT RECORDS: the first removes a GQ0
 *     reference BLOCK and leaves a weak variant alone, the second demotes the weak VARIANT and
 *     leaves the block;
 *   - A DEMOTED VARIANT BECOMES A ONE-BASE BLOCK OF ITS OWN and does NOT merge with the blocks
 *     either side of it, though they are in the same band;
 *   - --keep-all-alts CHANGES NOTHING for a variant that is already biallelic with <NON_REF>, so
 *     it is measured as the control that says the default was not trimming anything here;
 *   - THE TWO ANNOTATION ARGUMENTS ARE NOT SYMMETRIC: --format-annotations-to-remove takes a
 *     FORMAT key and --annotations-to-keep an INFO one, and asking the latter for a FORMAT key is
 *     REFUSED rather than ignored;
 *   - --floor-blocks WRITES THE BAND'S LOWER BOUND rather than the observed quality, and drops
 *     MIN_DP and PL from the block's FORMAT with it: the block keeps GT, DP and GQ alone;
 *   - AND THE <NON_REF> ALLELE SURVIVES EVERYTHING, because a GVCF without it is not a GVCF.
 *
 * Output:
 *
 *     vcf\tinput=<the input gvcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ReblockGVCFDump
 */

import org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ReblockGVCFDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##ALT=<ID=NON_REF,Description=\"Any other allele\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">",
                "##INFO=<ID=EXTRA,Number=1,Type=String,Description=\"An unrecognised annotation\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth\">",
                "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
                "##FORMAT=<ID=SPARE,Number=1,Type=Integer,Description=\"A spare format field\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
    }

    /** A reference block, which is a record whose only alternate is <NON_REF>. */
    static String block(final int start, final int end, final int quality, final int depth) {
        return "chr1\t" + start + "\t.\tA\t<NON_REF>\t.\t.\tEND=" + end
                + "\tGT:DP:GQ:MIN_DP:PL:SPARE\t0/0:" + depth + ":" + quality + ":" + depth + ":0,"
                + quality + "," + (quality * 10) + ":7";
    }

    /** A variant, whose PL[0] is what the reference-quality threshold reads. */
    static String variant(final int position, final String alternate, final int refLikelihood,
                          final int quality) {
        return "chr1\t" + position + "\t.\tA\t" + alternate + ",<NON_REF>\t50.00\t.\tEXTRA=note"
                + "\tGT:AD:DP:GQ:PL:SPARE\t0/1:8,4,0:12:" + quality + ":" + refLikelihood
                + ",0," + (quality * 10) + "," + refLikelihood + ",0," + refLikelihood + ":7";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("reblock-gvcf-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ReblockGVCFDump: how a GVCF's reference blocks are coarsened");

        final Path fasta = writeReference(dir);

        final List<String> sites = new ArrayList<>(header());
        // Three adjacent blocks: the first two in one band, the third far above it.
        sites.add(block(1000, 1099, 25, 20));
        sites.add(block(1100, 1199, 35, 22));
        sites.add(block(1200, 1299, 80, 30));
        // A confident variant, which survives everything.
        sites.add(variant(1300, "C", 900, 90));
        // Blocks either side of a WEAK variant, so a demoted one can merge with them.
        sites.add(block(1301, 1399, 25, 20));
        sites.add(variant(1400, "G", 5, 3));
        sites.add(block(1401, 1499, 25, 20));
        // A block at quality zero, which --drop-low-quals removes.
        sites.add(block(1500, 1599, 0, 2));
        sites.add("");
        final String input = String.join("\n", sites);
        final Path inputPath = write(dir, "input.g.vcf", input);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(input));

        run(dir, "default", inputPath, fasta, List.of());
        // A different set of band edges.
        run(dir, "one-band", inputPath, fasta, List.of("--gvcf-gq-bands", "60"));
        run(dir, "many-bands", inputPath, fasta,
                List.of("--gvcf-gq-bands", "10", "--gvcf-gq-bands", "20", "--gvcf-gq-bands", "30", "--gvcf-gq-bands", "40", "--gvcf-gq-bands", "50",
                        "--gvcf-gq-bands", "60"));
        // The reference-quality threshold, which demotes rather than drops.
        run(dir, "rgq-threshold", inputPath, fasta, List.of("--rgq-threshold", "10"));
        // And the argument that drops instead.
        run(dir, "drop-low-quals", inputPath, fasta, List.of("--drop-low-quals", "true"));
        run(dir, "drop-and-threshold", inputPath, fasta,
                List.of("--drop-low-quals", "true", "--rgq-threshold", "10"));
        // Every alternate and the full likelihood array.
        run(dir, "keep-all-alts", inputPath, fasta, List.of("--keep-all-alts", "true"));
        // The band's lower bound instead of the observed quality.
        run(dir, "floor-blocks", inputPath, fasta, List.of("--floor-blocks", "true"));
        // The two annotation lists, which are not the same list.
        run(dir, "keep-annotation", inputPath, fasta,
                List.of("--annotations-to-keep", "EXTRA"));
        run(dir, "remove-annotation", inputPath, fasta,
                List.of("--format-annotations-to-remove", "SPARE"));
        // A FORMAT key asked for as an INFO one, which is not an error and does nothing.
        run(dir, "keep-format-key", inputPath, fasta,
                List.of("--annotations-to-keep", "SPARE"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".g.vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new ReblockGVCF().instanceMain(argv.toArray(new String[0]));
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
