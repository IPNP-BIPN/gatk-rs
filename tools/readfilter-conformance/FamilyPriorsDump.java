/*
 * CalculateGenotypePosteriors' family priors, taken from the reference.
 *
 * The other half of the tool: a trio's three genotypes recomputed together, weighting every
 * combination by how many Mendelian violations it implies. The population half is measured by the
 * calculate-genotype-posteriors suite; this one is about the pedigree.
 *
 * Nine behaviours this is built to catch.
 *
 *   - FAMILY PRIORS ONLY APPLY TO BIALLELIC SITES: apply tests `variant.isBiallelic()` before
 *     calling the family engine at all, so a triallelic site keeps its likelihoods whatever the
 *     pedigree says;
 *   - THE NON-VIOLATION COEFFICIENT IS NOT ONE: it is `1 - 10*deNovoPrior - deNovoPrior^2`, so
 *     every consistent combination is scaled DOWN by ten times the de novo prior;
 *   - A TRIO NEEDS ALL THREE MEMBERS IN THE VCF: setTrios keeps a family only when exactly three
 *     of its members are present AND one of them has both parents among them;
 *   - NO PEDIGREE MEANS NO FAMILY PRIORS, and the tool says so and carries on rather than
 *     refusing;
 *   - AN UNCALLED PARENT BECOMES A UNIFORM THIRD, and the trio is still processed as a
 *     parent/child pair, but AN UNCALLED CHILD STOPS THE WHOLE TRIO;
 *   - THE JOINT TAGS ARE COMPUTED ONLY WHEN ALL THREE ARE CALLED, and are -1 otherwise;
 *   - THE JOINT LIKELIHOOD IS TAKEN AT THE POSTERIOR'S ARGMAX, not at its own, so JL reports the
 *     likelihood of the configuration the prior chose;
 *   - THE POSTERIORS OVERWRITE PP, which the population half then READS BACK rather than the PL,
 *     so the two halves compose in one direction only;
 *   - AND --de-novo-prior MOVES BOTH SIDES OF THE WEIGHTING, the violation term and the
 *     non-violation one.
 *
 * Output:
 *
 *     vcf\t<name>=<the whole input, escaped>
 *     ped\t<name>=<the whole pedigree, escaped>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FamilyPriorsDump
 */

import org.broadinstitute.hellbender.tools.walkers.variantutils.CalculateGenotypePosteriors;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FamilyPriorsDump {

    /**
     * A trio and an unrelated sample.
     *
     * The sites are chosen so that every branch of the family engine is reached: a consistent
     * inheritance, a Mendelian violation, an uncalled father, an uncalled child, a triallelic site
     * the engine never sees, and a site where the child's likelihoods are flat.
     */
    static final String INPUT = String.join("\n",
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=100000>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmom\tdad\tkid\tother",
            // Consistent: both parents het, child het.
            "chr1\t100\tconsistent\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/1:40:20:40,0,400"
                    + "\t0/1:40:20:40,0,400\t0/1:40:20:40,0,400\t0/0:30:20:0,30,300",
            // A violation: both parents hom-ref, child hom-var.
            "chr1\t200\tviolation\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/0:40:20:0,40,400"
                    + "\t0/0:40:20:0,40,400\t1/1:40:20:400,40,0\t0/0:30:20:0,30,300",
            // A weak violation, where the child's own likelihoods barely prefer the violating call.
            "chr1\t250\tweak\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/0:40:20:0,40,400"
                    + "\t0/0:40:20:0,40,400\t1/1:3:20:5,3,0\t0/0:30:20:0,30,300",
            // The father uncalled, which makes it a parent/child pair.
            "chr1\t300\tno-father\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/1:40:20:40,0,400"
                    + "\t./.:.:.:.\t0/1:40:20:40,0,400\t0/0:30:20:0,30,300",
            // The child uncalled, which stops the trio.
            "chr1\t400\tno-child\tA\tG\t50\t.\t.\tGT:GQ:DP:PL\t0/1:40:20:40,0,400"
                    + "\t0/1:40:20:40,0,400\t./.:.:.:.\t0/0:30:20:0,30,300",
            // Triallelic, which the family engine never sees.
            "chr1\t500\ttriallelic\tA\tG,T\t50\t.\t.\tGT:GQ:DP:PL"
                    + "\t0/1:40:20:40,0,400,40,400,400\t0/1:40:20:40,0,400,40,400,400"
                    + "\t1/2:40:20:400,200,100,200,0,100\t0/0:30:20:0,30,300,30,300,300",
            "");

    /** mom and dad are kid's parents; `other` is its own family. */
    static final String PED = String.join("\n",
            "fam1\tkid\tdad\tmom\t1\t2",
            "fam1\tdad\t0\t0\t1\t1",
            "fam1\tmom\t0\t0\t2\t1",
            "fam2\tother\t0\t0\t1\t1",
            "");

    /** The same pedigree with the mother missing, so no family has three members in the VCF. */
    static final String PED_NO_TRIO = String.join("\n",
            "fam1\tkid\tdad\t0\t1\t2",
            "fam1\tdad\t0\t0\t1\t1",
            "fam2\tother\t0\t0\t1\t1",
            "");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("family-priors-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FamilyPriorsDump: a trio's genotypes recomputed together");

        final Path input = write(dir, "input.vcf", INPUT);
        final Path ped = write(dir, "family.ped", PED);
        final Path pedNoTrio = write(dir, "no-trio.ped", PED_NO_TRIO);
        System.out.printf("vcf\tinput=%s%n", ReferenceQueryDump.escape(INPUT));
        System.out.printf("ped\tfamily=%s%n", ReferenceQueryDump.escape(PED));
        System.out.printf("ped\tno-trio=%s%n", ReferenceQueryDump.escape(PED_NO_TRIO));

        // Family priors alone, so nothing the population half does is in the way.
        run(dir, "family-only", input, List.of("--pedigree", ped.toString(),
                "--skip-population-priors", "true"));
        // The de novo prior, which moves both sides of the weighting.
        run(dir, "denovo-high", input, List.of("--pedigree", ped.toString(),
                "--skip-population-priors", "true", "--de-novo-prior", "0.001"));
        run(dir, "denovo-low", input, List.of("--pedigree", ped.toString(),
                "--skip-population-priors", "true", "--de-novo-prior", "1e-9"));
        // Both halves, in the order the tool applies them: family first, then population, which
        // reads the PP the family half just wrote.
        run(dir, "both-halves", input, List.of("--pedigree", ped.toString()));
        // Family priors skipped by argument, and by having no pedigree at all.
        run(dir, "skip-family", input, List.of("--pedigree", ped.toString(),
                "--skip-family-priors", "true", "--skip-population-priors", "true"));
        run(dir, "no-pedigree", input, List.of("--skip-population-priors", "true"));
        // A pedigree with no complete trio in the VCF, which is the same as none.
        run(dir, "no-trio", input, List.of("--pedigree", pedNoTrio.toString(),
                "--skip-population-priors", "true"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final List<String> extra)
            throws Exception {
        // Never `label + ".vcf"`: an output whose name collides with an input would be truncated
        // before it is read. The population-prior dump learned that the hard way.
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
