/*
 * CalculateGenotypePosteriors' output, taken from the reference.
 *
 * Genotype likelihoods turned into posteriors under a Dirichlet-multinomial prior built from allele
 * counts. Where those counts come from is the whole of the tool: a supporting panel, the input
 * samples themselves, an assumed number of hom-ref samples, or nothing at all, and each choice is
 * a different prior.
 *
 * Ten behaviours this is built to catch.
 *
 *   - A SITE MISSING FROM EVERY PANEL FALLS BACK TO FLAT PRIORS unless something supplies counts:
 *     useFlatPriors is true when the resources are empty AND the input samples are not used AND
 *     --num-reference-samples-if-no-call is zero, so the PP comes out equal to the PL;
 *   - THE INPUT SAMPLES ARE ONLY USED FOR A MISSING SITE IF THERE ARE TEN OF THEM, or if reference
 *     samples were assumed: `vc1.getNSamples() >= 10 || numRefSamplesFromMissingResources != 0`;
 *   - THE PANEL IS READ THROUGH MLEAC FIRST AND AC SECOND, which --default-to-allele-count flips;
 *   - --num-reference-samples-if-no-call IS DOUBLED, because the count is of diploid samples and
 *     what enters the prior is chromosomes;
 *   - THE SNP PRIOR IS CHOSEN BY ALLELE LENGTH, not by the site's type: an allele the same length
 *     as the reference takes the SNP pseudocount and every other allele the indel one, so one site
 *     can mix the two;
 *   - THE NON-REF ALLELE TAKES THE LARGER OF THE TWO PRIORS plus every count belonging to an
 *     allele the input did not carry;
 *   - GENOTYPES ARE RECALLED FROM THE POSTERIORS, so a sample's GT can change and its GQ is
 *     recomputed from the PP rather than the PL;
 *   - A HOM-REF BLOCK IS LEFT ALONE: a record whose only alternate is <NON_REF> keeps its AC, AF
 *     and AN and gains no GENOTYPE_PRIOR;
 *   - A MALFORMED MLEAC HEADER IS A REFUSAL rather than a fallback, checked for count type A and
 *     type Integer separately;
 *   - AND A VCF WITH NO GENOTYPES AT ALL IS REFUSED IN onTraversalStart.
 *
 * Output:
 *
 *     vcf\t<name>=<the whole input, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     header\t<label>=<the output's own header lines, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CalculateGenotypePosteriorsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.CalculateGenotypePosteriors;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CalculateGenotypePosteriorsDump {

    /**
     * The measured input.
     *
     * Three samples, so the ten-sample rule for missing sites is NOT met and the fallback to flat
     * priors is reachable. A biallelic SNP, a multiallelic SNP, an indel, a mixed-length site, a
     * hom-ref block ending in NON_REF, and a site the panel does not carry.
     */
    static final String INPUT = String.join("\n",
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=100000>",
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
            "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele count\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
            "##ALT=<ID=NON_REF,Description=\"Non-reference\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
            // A biallelic SNP the panel carries.
            "chr1\t100\tsnp\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300"
                    + "\t0/1:40:20:40,0,400\t1/1:20:20:300,20,0",
            // Multiallelic, so the prior is a triangle over three alleles.
            "chr1\t200\tmulti\tA\tG,T\t50\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300,30,300,300"
                    + "\t0/1:40:20:40,0,400,40,400,400\t1/2:20:20:300,200,100,200,0,100",
            // An indel, which takes the indel pseudocount.
            "chr1\t300\tindel\tAT\tA\t50\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300"
                    + "\t0/1:40:20:40,0,400\t0/0:20:20:0,20,200",
            // Reference-length and non-reference-length alternates at one site.
            "chr1\t400\tmixed\tAT\tGT,A\t50\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300,30,300,300"
                    + "\t0/1:40:20:40,0,400,40,400,400\t0/2:20:20:100,200,300,0,200,300",
            // A hom-ref block, whose only alternate is NON_REF.
            "chr1\t500\thomref\tA\t<NON_REF>\t.\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300"
                    + "\t0/0:40:20:0,40,400\t0/0:20:20:0,20,200",
            // A site no panel carries.
            "chr1\t600\tabsent\tC\tT\t50\t.\t.\tGT:GQ:DP:PL\t0/0:30:20:0,30,300"
                    + "\t0/1:40:20:40,0,400\t1/1:20:20:300,20,0",
            "");

    /** The supporting panel: AC and MLEAC that disagree, so which one is read is observable. */
    static final String PANEL = String.join("\n",
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=100000>",
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
            "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele count\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tG\t.\t.\tAC=200;AN=2000;MLEAC=20",
            "chr1\t200\t.\tA\tG,T\t.\t.\tAC=100,300;AN=2000;MLEAC=10,30",
            "chr1\t300\t.\tAT\tA\t.\t.\tAC=400;AN=2000;MLEAC=40",
            "chr1\t400\t.\tAT\tGT,A\t.\t.\tAC=100,200;AN=2000;MLEAC=10,20",
            "chr1\t500\t.\tA\tG\t.\t.\tAC=50;AN=2000;MLEAC=5",
            "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("calculate-genotype-posteriors-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CalculateGenotypePosteriorsDump: likelihoods turned into posteriors");

        final Path input = write(dir, "input.vcf", INPUT);
        final Path panel = write(dir, "panel.vcf", PANEL);
        // A supporting panel is QUERIED BY INTERVAL, so it has to be indexed: without an index the
        // run is refused rather than falling back to a linear read. Measured once, below.
        new IndexFeatureFile().instanceMain(new String[] {"-I", panel.toString()});
        final Path unindexed = write(dir, "unindexed.vcf", PANEL);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(INPUT));
        System.out.printf("vcf\tpanel=%s%n", ReferenceQueryDump.escape(PANEL));

        // No panel at all, which is where the flat-prior fallback lives.
        run(dir, "no-panel", input, List.of());
        // With the panel, which is the ordinary path.
        run(dir, "panel", input, List.of("--supporting", panel.toString()));
        // The panel plus assumed reference samples, which also turns the input samples on.
        run(dir, "ref-samples", input, List.of("--num-reference-samples-if-no-call", "1000"));
        run(dir, "panel-and-ref", input, List.of("--supporting", panel.toString(),
                "--num-reference-samples-if-no-call", "1000"));
        // AC instead of MLEAC, which the panel makes observable.
        run(dir, "default-to-ac", input, List.of("--supporting", panel.toString(),
                "--default-to-allele-count", "true"));
        // External counts only.
        run(dir, "ignore-input", input, List.of("--supporting", panel.toString(),
                "--ignore-input-samples", "true"));
        // And the discovered counts turned off for missing sites.
        run(dir, "discovered-off", input, List.of("--supporting", panel.toString(),
                "--num-reference-samples-if-no-call", "1000",
                "--discovered-allele-count-priors-off", "true"));
        // Flat priors for indels, which reaches the mixed-length site too.
        run(dir, "flat-indels", input, List.of("--supporting", panel.toString(),
                "--use-flat-priors-for-indels", "true"));
        // The pseudocounts themselves.
        run(dir, "priors", input, List.of("--supporting", panel.toString(),
                "--global-prior-snp", "0.01", "--global-prior-indel", "0.0001"));
        // Nothing applied at all.
        run(dir, "skip-population", input, List.of("--supporting", panel.toString(),
                "--skip-population-priors", "true"));
        // The same panel without an index.
        run(dir, "unindexed-panel", input, List.of("--supporting", unindexed.toString()));

        // A panel carrying AC but no MLEAC at all, which is what makes the preference observable
        // from the outside.
        final Path acOnly = write(dir, "ac-only.vcf", PANEL.replace(";MLEAC=20", "")
                .replace(";MLEAC=10,30", "").replace(";MLEAC=40", "")
                .replace(";MLEAC=10,20", "").replace(";MLEAC=5", "")
                .replace("##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele count\">\n", ""));
        new IndexFeatureFile().instanceMain(new String[] {"-I", acOnly.toString()});
        run(dir, "panel-ac-only", input, List.of("--supporting", acOnly.toString()));

        // A VCF with no genotype columns.
        final Path sites = write(dir, "sites.vcf", String.join("\n",
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=100000>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t100\t.\tA\tG\t50\t.\t.",
                ""));
        run(dir, "no-genotypes", sites, List.of());

        // Two malformed MLEAC headers, each checked separately and in order.
        run(dir, "mleac-count", write(dir, "mleac-count.vcf",
                INPUT.replace("##INFO=<ID=MLEAC,Number=A,Type=Integer",
                        "##INFO=<ID=MLEAC,Number=1,Type=Integer")), List.of());
        run(dir, "mleac-type", write(dir, "mleac-type.vcf",
                INPUT.replace("##INFO=<ID=MLEAC,Number=A,Type=Integer",
                        "##INFO=<ID=MLEAC,Number=A,Type=Float")), List.of());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final List<String> extra)
            throws Exception {
        // NOT `label + ".vcf"`: the run labelled `panel` would then write to `panel.vcf`, which is
        // the supporting panel itself. createVCFWriter opens the output in onTraversalStart, so
        // the tool truncates its own input before reading it, and every later run reads the
        // truncated file. The first version of this dump did exactly that.
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-V", input.toString(), "-O", out.toString()));
        argv.addAll(extra);
        try {
            new CalculateGenotypePosteriors().instanceMain(argv.toArray(new String[0]));
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
        final String text = Files.readString(out);
        final StringBuilder body = new StringBuilder();
        final StringBuilder header = new StringBuilder();
        for (final String line : text.split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            if (line.startsWith("##")) {
                // The command line the tool stamps carries paths and a version; the header lines
                // that matter are the ones it ADDS, which are the FORMAT and INFO ones.
                if (line.startsWith("##FORMAT") || line.startsWith("##INFO")) {
                    header.append(line).append("\n");
                }
            } else {
                body.append(line).append("\n");
            }
        }
        System.out.printf("header\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(header.toString(), dir)));
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(body.toString(), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
