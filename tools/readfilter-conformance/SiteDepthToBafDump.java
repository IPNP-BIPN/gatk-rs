/*
 * SiteDepthtoBAF's output, taken from the reference.
 *
 * Per-sample allele depths at a set of sites turned into B-allele fractions, with three filters and
 * an adjustment that moves every value at a locus at once.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE REF AND ALT ALLELES COME FROM THE SITES VCF, not from the depths, so which of the four
 *     counts is "ref" and which is "alt" is decided by another file entirely;
 *   - A SITE WHOSE TOTAL DEPTH IS BELOW --min-total-depth PRODUCES NOTHING, and the total is the
 *     sum of ALL FOUR counts rather than of the two the alleles name;
 *   - THE HET TEST IS A CHI-SQUARED FIT WITH ONE DEGREE OF FREEDOM against an expectation of half
 *     the total, computed from the ref and alt counts alone, and a sample whose fit probability is
 *     below --min-het-probability is dropped;
 *   - A LOCUS WITH EXACTLY ONE SURVIVING SAMPLE IS WRITTEN AS 0.5, whatever its measured fraction:
 *     the value is replaced, not adjusted;
 *   - A LOCUS WITH TWO OR MORE IS ADJUSTED BY ITS OWN MEDIAN, every value moved by `0.5 - median`,
 *     so the numbers written are not the fractions measured;
 *   - AND THAT WHOLE LOCUS IS DROPPED WHEN THE STANDARD DEVIATION EXCEEDS --max-std, which is the
 *     POPULATION deviation of the surviving values;
 *   - THE MEDIAN OF AN EVEN COUNT IS THE MEAN OF THE TWO MIDDLE VALUES;
 *   - THE VALUE IS WRITTEN THROUGH `DecimalFormat("#.00")`, so 0.5 comes out `.50` with no leading
 *     zero and a half-way value rounds the way that formatter rounds;
 *   - THE SITES FILE IS WALKED IN LOCKSTEP and must not run out first, nor disagree about the
 *     locus, and either is a UserException naming both positions;
 *   - AND THE SITE DEPTH FILE IS ZERO-BASED ON DISK AND ONE-BASED INSIDE, like every other SV
 *     evidence format.
 *
 * Output:
 *
 *     depths\t<label>=<the whole .sd.txt, escaped>
 *     sites\t<label>=<the whole sites vcf, escaped>
 *     baf\t<label>=<the whole .baf.txt, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SiteDepthToBafDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.sv.SiteDepthtoBAF;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SiteDepthToBafDump {

    static final String VCF_HEADER =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=10000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    /** One site-depth row: contig, zero-based start, sample, then the four counts A C G T. */
    static String depth(final int start, final String sample, final int a, final int c,
                        final int g, final int t) {
        return "chr1\t" + start + "\t" + sample + "\t" + a + "\t" + c + "\t" + g + "\t" + t + "\n";
    }

    static String site(final int position, final String ref, final String alt) {
        return "chr1\t" + position + "\t.\t" + ref + "\t" + alt + "\t.\t.\t.\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("site-depth-to-baf-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SiteDepthToBafDump: allele depths turned into B-allele fractions");

        // The evidence files carry no dictionary of their own, so the tool is given one: a fasta
        // of the single contig the fixtures use, and the .dict CreateSequenceDictionary writes.
        final Path fasta = dir.resolve("ref.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < 1000; i++) {
            bases.append("ACGTACGTAC");
            if (i % 6 == 5) {
                bases.append("\n");
            }
        }
        bases.append("\n");
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        // Locus 100: three samples, a clean het each, so the locus is adjusted by its median.
        // Locus 200: one sample only, which is written as .50 whatever it measured.
        // Locus 300: a sample below the depth floor, and one that is far from het.
        // Locus 400: three samples whose fractions are far apart, so the locus is dropped.
        final String depths =
                depth(99, "s1", 10, 0, 10, 0)
                + depth(99, "s2", 12, 0, 8, 0)
                + depth(99, "s3", 8, 0, 12, 0)
                + depth(199, "s1", 11, 0, 9, 0)
                + depth(299, "s1", 3, 0, 2, 0)
                + depth(299, "s2", 19, 0, 1, 0)
                + depth(399, "s1", 10, 0, 10, 0)
                + depth(399, "s2", 16, 0, 4, 0)
                + depth(399, "s3", 4, 0, 16, 0)
                // Locus 500: three samples whose median is 0.6, so every value moves by -0.1.
                + depth(499, "s1", 6, 0, 14, 0)
                + depth(499, "s2", 8, 0, 12, 0)
                + depth(499, "s3", 10, 0, 10, 0)
                // Locus 600: three whose spread is wide enough to fail --max-std at a low het
                // floor, so the whole locus is dropped rather than adjusted.
                + depth(599, "s1", 5, 0, 15, 0)
                + depth(599, "s2", 10, 0, 10, 0)
                + depth(599, "s3", 15, 0, 5, 0);
        // The sites: ref A, alt G, so the counts that matter are the first and the third.
        final String sites = VCF_HEADER
                + site(100, "A", "G") + site(200, "A", "G") + site(300, "A", "G")
                + site(400, "A", "G") + site(500, "A", "G") + site(600, "A", "G");

        run(dir, "defaults", depths, sites);
        // The chi-squared floor moved, which changes which samples survive.
        run(dir, "het-0.9", depths, sites, "--min-het-probability", "0.9");
        run(dir, "het-0.05", depths, sites, "--min-het-probability", "0.05");
        // The depth floor moved, which drops the shallow sample at locus 300.
        run(dir, "depth-20", depths, sites, "--min-total-depth", "20");
        // The deviation limit loosened, which lets locus 400 through.
        run(dir, "std-0.5", depths, sites, "--max-std", "0.5");
        // A floor low enough to keep three samples at loci 500 and 600, where the median
        // adjustment and the deviation limit are both visible.
        run(dir, "het-0.01", depths, sites, "--min-het-probability", "0.01");
        run(dir, "het-0.01-std-0.5", depths, sites,
                "--min-het-probability", "0.01", "--max-std", "0.5");
        // A site whose ref is not a single base, and one where the sites file runs out.
        run(dir, "bad-ref", depths, VCF_HEADER + site(100, "N", "G") + site(200, "A", "G")
                + site(300, "A", "G") + site(400, "A", "G")
                + site(500, "A", "G") + site(600, "A", "G"));
        run(dir, "short-sites", depths, VCF_HEADER + site(100, "A", "G"));
        // And one where the sites file names a locus the depths do not.
        run(dir, "wrong-locus", depths, VCF_HEADER + site(150, "A", "G") + site(200, "A", "G")
                + site(300, "A", "G") + site(400, "A", "G")
                + site(500, "A", "G") + site(600, "A", "G"));
    }

    static void run(final Path dir, final String label, final String depths, final String sites,
                    final String... extra) throws Exception {
        final Path depthFile = dir.resolve(label + ".sd.txt");
        Files.writeString(depthFile, depths, StandardCharsets.UTF_8);
        final Path sitesFile = dir.resolve(label + "-sites.vcf");
        Files.writeString(sitesFile, sites, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", sitesFile.toString()});
        System.out.printf("depths\t%s=%s%n", label, ReferenceQueryDump.escape(depths));
        System.out.printf("sites\t%s=%s%n", label, ReferenceQueryDump.escape(sites));

        final Path out = dir.resolve("baf-" + label + ".baf.txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--site-depth", depthFile.toString(),
                "--baf-sites-vcf", sitesFile.toString(),
                "-O", out.toString(),
                "--sequence-dictionary", dir.resolve("ref.dict").toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new SiteDepthtoBAF().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("baf\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
