/*
 * SelectVariants' --discordance and --concordance, taken from the reference.
 *
 * Both compare the record being read against whatever a second VCF has at the same position, and
 * both change meaning entirely depending on whether any sample was named.
 *
 * Nine behaviours this is built to catch, three of which are not comparisons of genotypes at all.
 *
 *   - WITHOUT -sn, --discordance IS "THE OTHER FILE HAS NOTHING HERE" and nothing else: the
 *     genotypes are never looked at, so a record whose calls disagree with the other file's is
 *     CONCORDANT as long as some record is there;
 *   - AND WITHOUT -sn, --concordance IS "THE OTHER FILE HAS SOMETHING HERE", equally blind;
 *   - WITH -sn THEY COMPARE ALLELE SETS AND NOT GENOTYPES. `haveSameGenotypes` asks each list to
 *     contain the other, so 0/1 and 1/0 are the same call, and so are 1/1 and a haploid 1: the
 *     multiplicity is invisible;
 *   - A FILTERED GENOTYPE NEVER MATCHES ANYTHING, INCLUDING ANOTHER FILTERED ONE. The first two
 *     clauses of `haveSameGenotypes` are `g1.isCalled() && g2.isFiltered()` and its mirror, and
 *     `isCalled()` is about ALLELES rather than filters: a genotype with alleles is called whether
 *     or not it carries an FT. So two identically filtered genotypes take the first clause and are
 *     declared different, and the third clause, `both filtered && excludeFiltered`, IS DEAD CODE;
 *   - AND --exclude-filtered CHANGES NOTHING HERE FOR THE SAME REASON. `sampleHasVariant` reads
 *     `!isHomRef && (isCalled || (isFiltered && !excludeFiltered))`, whose second half can only
 *     matter for a genotype that is not called, and a filtered genotype still is. Measured: the
 *     flag moves no record in either direction;
 *   - DISCORDANCE ONLY CONSIDERS SAMPLES THAT CARRY A VARIANT: a hom-ref genotype is skipped
 *     entirely, so a record where the selected sample is 0/0 is never discordant no matter what
 *     the other file says;
 *   - CONCORDANCE REQUIRES EVERY SELECTED SAMPLE TO MATCH, so it is not the negation of
 *     discordance: a record can be neither, and a record can be both;
 *   - AND A DIFFERENT ALTERNATE IS A DIFFERENT ALLELE SET, so the same `0/1` on both sides is a
 *     mismatch when the two records do not call the same alternate.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     kept\t<label>\t<comma-joined positions of the records written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SelectVariantsConcordanceDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SelectVariantsConcordanceDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FILTER=<ID=LowGQ,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("selectvariants-concordance-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SelectVariantsConcordanceDump: comparing two files, from the reference");

        final Path input = writeVcf(dir, "records",
                // 100: the other file has the same calls.
                "chr1\t100\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                // 200: the other file has the same alleles written the other way round.
                "chr1\t200\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                // 300: the other file calls a different genotype.
                "chr1\t300\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                // 400: the other file has nothing here.
                "chr1\t400\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                // 500: the selected sample is hom-ref, which discordance skips.
                "chr1\t500\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/0:60\t0/1:60",
                // 600: this file's genotype is filtered, the other file's is not.
                "chr1\t600\t.\tA\tC\t50\t.\tDP=30\tGT:GQ:FT\t0/1:60:LowGQ\t0/1:60:PASS",
                // 700: both files' genotypes are filtered.
                "chr1\t700\t.\tA\tC\t50\t.\tDP=30\tGT:GQ:FT\t0/1:60:LowGQ\t0/1:60:PASS",
                // 800: the other file calls `0/1` here too, but of a different alternate.
                "chr1\t800\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                // 900: a haploid call against a diploid one made of the same allele.
                "chr1\t900\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t1/1:60\t0/1:60");

        final Path comp = writeVcf(dir, "comparison",
                "chr1\t100\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t200\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t1/0:60\t1/0:60",
                "chr1\t300\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t1/1:60\t1/1:60",
                "chr1\t500\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t1/1:60\t1/1:60",
                "chr1\t600\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t700\t.\tA\tC\t50\t.\tDP=30\tGT:GQ:FT\t0/1:60:LowGQ\t0/1:60:LowGQ",
                "chr1\t800\t.\tA\tG\t50\t.\tDP=30\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t900\t.\tA\tC\t50\t.\tDP=30\tGT:GQ\t1:60\t1:60");

        // Neither flag, for the baseline.
        run(dir, "no-comparison", input);

        // Without a sample, both flags are about presence alone.
        run(dir, "discordance", input, "--discordance", comp.toString());
        run(dir, "concordance", input, "--concordance", comp.toString());

        // With a sample, both compare allele sets.
        run(dir, "discordance-one-sample", input, "--discordance", comp.toString(), "-sn", "s0");
        run(dir, "concordance-one-sample", input, "--concordance", comp.toString(), "-sn", "s0");
        run(dir, "discordance-both-samples", input, "--discordance", comp.toString(),
                "-sn", "s0", "-sn", "s1");
        run(dir, "concordance-both-samples", input, "--concordance", comp.toString(),
                "-sn", "s0", "-sn", "s1");

        // And the flag that changes what "the same genotype" means.
        run(dir, "discordance-exclude-filtered", input, "--discordance", comp.toString(),
                "-sn", "s0", "--exclude-filtered", "true");
        run(dir, "concordance-exclude-filtered", input, "--concordance", comp.toString(),
                "-sn", "s0", "--exclude-filtered", "true");
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

    static void run(final Path dir, final String label, final Path input,
                    final String... arguments) {
        final Path output = dir.resolve(label + "-out.vcf");
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(),
                "-O", output.toString()));
        all.addAll(List.of(arguments));
        try {
            new SelectVariants().instanceMain(all.toArray(new String[0]));
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
        final List<String> kept = new ArrayList<>();
        for (final String line : lines) {
            if (line.startsWith("#")) {
                continue;
            }
            kept.add(line.split("\t")[1]);
        }
        System.out.printf("kept\t%s\t%s%n", label, String.join(",", kept));
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
