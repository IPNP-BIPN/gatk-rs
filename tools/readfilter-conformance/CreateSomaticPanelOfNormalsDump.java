/*
 * CreateSomaticPanelOfNormals' panel, taken from the reference.
 *
 * Which sites recur across normal samples often enough to be called artefacts, and with what
 * allele-fraction shape. A site is kept when enough samples carry an alternate that is more likely
 * to be an artefact than to be germline, and the shape is a beta fitted to those samples' counts.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A SITE WITH NO ALTERNATE ALLELE IS SKIPPED, and so is one whose only alternate is the
 *     spanning deletion;
 *   - THE GERMLINE TEST IS SKIPPED ENTIRELY FOR A MULTIALLELIC SITE: every genotype counts,
 *     whatever its counts say, which is a documented TODO rather than a rule;
 *   - A GENOTYPE WITH NO ALTERNATE READ NEVER COUNTS, whatever the germline frequency;
 *   - WITH NO GERMLINE RESOURCE EVERY ALTERNATE COUNTS, because a frequency of zero is below the
 *     negligible threshold and the germline probability is returned as zero;
 *   - A HIGH GERMLINE FREQUENCY REMOVES A HET-LOOKING GENOTYPE and leaves a low-fraction one, which
 *     is the whole point of the resource;
 *   - --min-sample-count IS COMPARED AGAINST THE SURVIVORS, not against the samples;
 *   - FRACTION IS OVER ALL SAMPLES IN THE HEADER, not over the survivors, so it falls when a
 *     sample is added that carries nothing;
 *   - THE BETA IS FITTED BY A BRENT SEARCH over a scale, whose base shape is the empirical mean and
 *     whose answer therefore carries the optimiser's tolerances;
 *   - AND A GENOTYPE WITHOUT AD CONTRIBUTES NOTHING TO THE FIT even when it counted for the site.
 *
 * Output:
 *
 *     vcf\t<label>=<the whole input vcf, escaped>
 *     germline\tmain=<the germline resource vcf, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CreateSomaticPanelOfNormalsDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.CreateSomaticPanelOfNormals;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CreateSomaticPanelOfNormalsDump {

    static final int CONTIG_LENGTH = 199980;

    static List<String> header(final List<String> samples) {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + String.join("\t", samples)));
    }

    /** One site. `genotypes` are `GT:AD` strings, or `.` for a genotype with no AD at all. */
    static String site(final int position, final String reference, final String alternate,
                       final String... genotypes) {
        final List<String> columns = new ArrayList<>();
        for (final String genotype : genotypes) {
            columns.add(genotype);
        }
        return "chr1\t" + position + "\t.\t" + reference + "\t" + alternate + "\t.\t.\t.\tGT:AD\t"
                + String.join("\t", columns);
    }

    static String buildVcf(final List<String> samples, final List<String> sites) {
        final List<String> lines = header(samples);
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    /** Three normals, one site per behaviour. */
    static String buildMain() {
        final List<String> sites = new ArrayList<>();
        // No alternate at all.
        sites.add(site(1000, "A", ".", "0/0:20,0", "0/0:20,0", "0/0:20,0"));
        // Only the spanning deletion.
        sites.add(site(2000, "A", "*", "0/1:18,2", "0/1:18,2", "0/0:20,0"));
        // Two samples carrying a clear low-fraction alternate, one carrying none.
        sites.add(site(3000, "A", "C", "0/1:18,2", "0/1:17,3", "0/0:20,0"));
        // Only ONE sample carrying it, which the default minimum of two refuses.
        sites.add(site(4000, "A", "C", "0/1:18,2", "0/0:20,0", "0/0:20,0"));
        // All three carrying it, so the fraction is one.
        sites.add(site(5000, "A", "C", "0/1:18,2", "0/1:17,3", "0/1:16,4"));
        // Two samples at a half fraction, which looks germline where the resource says so.
        sites.add(site(6000, "A", "C", "0/1:10,10", "0/1:11,9", "0/0:20,0"));
        // Multiallelic, where the germline test is skipped for every genotype.
        sites.add(site(7000, "A", "C,G", "0/1:10,10,0", "0/1:11,9,0", "0/0:20,0,0"));
        // One sample with no AD at all beside one with counts.
        sites.add(site(8000, "A", "C", "0/1:18,2", "0/1", "0/1:16,4"));
        // Very deep counts, which sharpen the fitted beta.
        sites.add(site(9000, "A", "C", "0/1:1800,200", "0/1:1700,300", "0/0:2000,0"));
        return buildVcf(List.of("n1", "n2", "n3"), sites);
    }

    /** The same sites with a fourth sample that carries nothing, to move FRACTION. */
    static String buildFourSamples() {
        final List<String> sites = new ArrayList<>();
        sites.add(site(3000, "A", "C", "0/1:18,2", "0/1:17,3", "0/0:20,0", "0/0:20,0"));
        sites.add(site(5000, "A", "C", "0/1:18,2", "0/1:17,3", "0/1:16,4", "0/0:20,0"));
        return buildVcf(List.of("n1", "n2", "n3", "n4"), sites);
    }

    /** A germline resource giving 6000 a high frequency and 3000 none. */
    static String buildGermline() {
        return String.join("\n",
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t6000\t.\tA\tC\t.\t.\tAF=0.4",
                "chr1\t7000\t.\tA\tC\t.\t.\tAF=0.4",
                "chr1\t9000\t.\tA\tC\t.\t.\tAF=0.4",
                "");
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("panel-of-normals-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CreateSomaticPanelOfNormalsDump: which sites recur across normals "
                + "often enough to be artefacts");

        final String main = buildMain();
        final String four = buildFourSamples();
        final String germline = buildGermline();
        final Path mainPath = write(dir, "normals.vcf", main);
        final Path fourPath = write(dir, "four.vcf", four);
        final Path germlinePath = write(dir, "germline.vcf", germline);
        // The germline resource is QUERIED by interval, so it needs an index beside it: without
        // one the run is refused before a record is read.
        htsjdk.tribble.index.Index index = htsjdk.tribble.index.IndexFactory.createLinearIndex(
                germlinePath.toFile(), new htsjdk.variant.vcf.VCFCodec());
        index.writeBasedOnFeatureFile(germlinePath.toFile());
        System.out.printf("vcf\tmain=%s%n", ReferenceQueryDump.escape(main));
        System.out.printf("vcf\tfour=%s%n", ReferenceQueryDump.escape(four));
        System.out.printf("germline\tmain=%s%n", ReferenceQueryDump.escape(germline));

        run(dir, "default", mainPath, List.of());
        // With the germline resource, which removes the half-fraction site.
        run(dir, "germline", mainPath, List.of("--germline-resource", germlinePath.toString()));
        // A germline probability threshold high enough to keep everything again.
        run(dir, "germline-permissive", mainPath, List.of(
                "--germline-resource", germlinePath.toString(),
                "--max-germline-probability", "1.0"));
        // And one low enough to drop what the default keeps.
        run(dir, "germline-strict", mainPath, List.of(
                "--germline-resource", germlinePath.toString(),
                "--max-germline-probability", "0.0001"));
        // A minimum of one sample, which keeps the singleton site.
        run(dir, "min-one", mainPath, List.of("--min-sample-count", "1"));
        // A minimum of three, which drops the two-sample sites.
        run(dir, "min-three", mainPath, List.of("--min-sample-count", "3"));
        // A fourth sample carrying nothing, which moves FRACTION without moving anything else.
        run(dir, "four-samples", fourPath, List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(),
                "-O", out.toString()));
        argv.addAll(extra);
        try {
            new CreateSomaticPanelOfNormals().instanceMain(argv.toArray(new String[0]));
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

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
