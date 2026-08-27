/*
 * GenotypeGVCFs' calls, taken from the reference.
 *
 * How a GVCF becomes a genotyped VCF. Every reference block is dropped, every variant is
 * re-genotyped from its likelihoods against a calling threshold, and <NON_REF> is removed from the
 * alleles that survive.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A REFERENCE BLOCK IS NEVER WRITTEN: the output is a VCF and not a GVCF;
 *   - <NON_REF> IS REMOVED FROM THE ALTERNATES of every site that is written;
 *   - A SITE WHOSE BEST GENOTYPE IS REFERENCE IS DROPPED, and --include-non-variant-sites writes it
 *     anyway, which is the only way to see what the calling threshold removed;
 *   - THE CALLING THRESHOLD DOES NOT DECIDE EMISSION: what is written is decided by the CALL, so
 *     a site whose called genotype is reference is dropped at a threshold of 2 as surely as at 50,
 *     and only --include-non-variant-sites brings it back. The two threshold runs are here as the
 *     controls that say so;
 *   - AN ALTERNATE NO SAMPLE CARRIES IS DROPPED and the likelihoods are re-indexed around it;
 *   - THE GENOTYPES ARE CALLED FROM THE LIKELIHOODS, not copied, so a sample's GT can differ from
 *     the one the GVCF carried;
 *   - THE SITE ANNOTATIONS ARE COMPUTED FROM THE CALLED GENOTYPES, so AC, AN and AF describe the
 *     output;
 *   - --keep-combined-raw-annotations RETAINS THE RAW FIELDS beside the finalised ones;
 *   - AND --sample-ploidy CHANGES NOTHING when the likelihood arrays are diploid, which is the
 *     control that says the ploidy is read off the data rather than off the argument.
 *
 * Output:
 *
 *     vcf\tinput=<the gvcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GenotypeGVCFsDump
 */

import org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class GenotypeGVCFsDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header() {
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
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"));
    }

    static String block(final int start, final int end) {
        return "chr1\t" + start + "\t.\tA\t<NON_REF>\t.\t.\tEND=" + end
                + "\tGT:DP:GQ:MIN_DP:PL\t0/0:20:40:20:0,40,400";
    }

    /** A variant, whose PL array decides what it is genotyped as. */
    static String variant(final int position, final String alternates, final String genotype,
                          final String likelihoods) {
        final int alleles = alternates.split(",").length + 2;
        final StringBuilder depths = new StringBuilder("8");
        for (int i = 1; i < alleles; i++) {
            depths.append(",").append(i == 1 ? 4 : 0);
        }
        return "chr1\t" + position + "\t.\tA\t" + alternates + ",<NON_REF>\t.\t.\t."
                + "\tGT:AD:DP:GQ:PL\t" + genotype + ":" + depths + ":12:50:" + likelihoods;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("genotype-gvcfs-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GenotypeGVCFsDump: how a GVCF becomes a genotyped VCF");

        final Path fasta = writeReference(dir);

        final List<String> sites = new ArrayList<>(header());
        // A reference block, which never reaches the output.
        // Five bases only: --include-non-variant-sites expands a block into one record per BASE,
        // so a longer one would bury every other run's output in the report.
        sites.add(block(1000, 1004));
        // A confident heterozygote: the middle likelihood is the best.
        sites.add(variant(1100, "C", "0/1", "900,0,900,900,900,900"));
        // A confident homozygous variant.
        sites.add(variant(1200, "C", "1/1", "900,900,0,900,900,900"));
        // A site whose best likelihood is the REFERENCE, which the threshold drops.
        sites.add(variant(1300, "C", "0/0", "0,60,600,60,600,600"));
        // A marginal site, three quality points above a reference call.
        sites.add(variant(1400, "C", "0/1", "3,0,900,3,3,3"));
        // Two alternates, only the first of which the sample carries.
        sites.add(variant(1500, "C,G", "0/1",
                "900,0,900,900,900,900,900,900,900,900"));
        sites.add("");
        final String input = String.join("\n", sites);
        final Path inputPath = write(dir, "input.g.vcf", input);
        htsjdk.tribble.index.IndexFactory.createLinearIndex(inputPath.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(inputPath.toFile());
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(input));

        run(dir, "default", inputPath, fasta, List.of());
        run(dir, "all-sites", inputPath, fasta,
                List.of("--include-non-variant-sites", "true"));
        // The calling threshold, above and below the marginal site.
        run(dir, "call-threshold-2", inputPath, fasta,
                List.of("--standard-min-confidence-threshold-for-calling", "2"));
        run(dir, "call-threshold-50", inputPath, fasta,
                List.of("--standard-min-confidence-threshold-for-calling", "50"));
        // The raw annotations kept beside the finalised ones.
        run(dir, "keep-combined", inputPath, fasta,
                List.of("--keep-combined-raw-annotations", "true"));
        // A ploidy the genotyper is told to assume, which does not match the likelihood arrays.
        run(dir, "ploidy-one", inputPath, fasta, List.of("--sample-ploidy", "1"));
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
            new GenotypeGVCFs().instanceMain(argv.toArray(new String[0]));
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
            System.out.printf("none\t%s=no output file%n", label);
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
