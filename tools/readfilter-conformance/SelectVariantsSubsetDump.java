/*
 * What SelectVariants does to a record once it knows which samples to keep, taken from the
 * reference.
 *
 * `subsetGenotypesBySampleNames` and `addAnnotations` between them decide the whole of the output
 * record: which columns survive, which alternate alleles survive, which annotations are recomputed
 * and which are carried, and whether the alleles are trimmed afterwards.
 *
 * Nine behaviours this is built to catch, four of which are not what the argument names suggest.
 *
 *   - THE INFO DP IS RECOMPUTED FROM THE GENOTYPES AND OVERWRITES THE RECORD'S OWN, and it sums
 *     only the UNFILTERED genotypes of the SELECTED samples. A record whose DP said 200 comes out
 *     saying whatever the kept columns add up to, and a genotype carrying FT contributes nothing;
 *   - AND IT IS WRITTEN ONLY IF SOME KEPT GENOTYPE HAD A DP. Where none does, the record keeps
 *     the DP it arrived with, so the same argument either replaces the field or leaves it,
 *     depending on the FORMAT column;
 *   - AC, AF AND AN ARE RECOMPUTED by calculateChromosomeCounts, from the kept genotypes. AN is
 *     the called chromosome count, so a no-call reduces it; AC is a list when there is more than
 *     one alternate and a scalar when there is one, which changes the INFO field's shape;
 *   - AF IS A DOUBLE AND THE WRITER FORMATS IT, so 0.5 comes out `0.500` and 0.0 comes out `0.00`:
 *     two different numbers of decimals in the same column of the same file;
 *   - MLEAC AND MLEAF ARE STRIPPED whenever the record is rewritten, because they describe a
 *     calling that no longer applies. AC and AF are not stripped, they are replaced;
 *   - --remove-unused-alternates DROPS AN ALTERNATE WITH NO CALLS and rewrites PL and AD through
 *     AlleleSubsettingUtils, so the genotype fields shrink with the allele list;
 *   - THE RESULT IS THEN TRIMMED unless --preserve-alleles, so removing an alternate can also
 *     shorten the reference allele and MOVE NOTHING while changing every allele's spelling;
 *   - --keep-original-ac WRITES AC_Orig, AF_Orig AND AN_Orig FROM THE INPUT'S OWN ATTRIBUTES, not
 *     from a recount, and reorders the per-allele ones when the allele list changed. A record
 *     without AC in its INFO gets no AC_Orig at all;
 *   - AND A RECORD IS RETURNED UNTOUCHED when nothing was selected and no alternate was removed,
 *     which is what keeps a whole-cohort run from gaining annotations it did not have.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one record line of the output VCF, escaped>
 *     samples\t<label>\t<comma-joined sample columns of the output header>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SelectVariantsSubsetDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class SelectVariantsSubsetDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                    + "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Max likelihood allele count\">\n"
                    + "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Max likelihood allele frequency\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n"
                    + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">\n"
                    + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
                    + "##FILTER=<ID=LowGQ,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("selectvariants-subset-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SelectVariantsSubsetDump: subsetting a record, from the reference");

        final Path input = writeVcf(dir, "records",
                // Biallelic, and the only carriers of the alternate are s2 and s3: subsetting to
                // s0 and s1 leaves an alternate nobody calls.
                "chr1\t100\t.\tA\tC\t50\t.\tDP=200;AC=3;AF=0.375;AN=8;MLEAC=3;MLEAF=0.375"
                        + "\tGT:GQ:DP:AD:PL\t0/0:60:10:10,0:0,30,300\t0/0:60:20:20,0:0,60,600"
                        + "\t0/1:60:30:15,15:300,0,300\t1/1:60:40:0,40:800,80,0",
                // Multi-allelic with a 1/2 call, so the allele list and every PL depend on which
                // samples are kept.
                "chr1\t200\t.\tA\tC,G\t50\t.\tDP=200;AC=2,2;AF=0.250,0.250;AN=8;MLEAC=2,2;MLEAF=0.250,0.250"
                        + "\tGT:GQ:DP:AD:PL\t0/1:50:10:5,5,0:100,0,200,150,250,300"
                        + "\t0/2:50:20:10,0,10:100,150,300,0,250,200\t1/2:50:30:0,15,15:300,200,100,150,0,250"
                        + "\t0/0:50:40:40,0,0:0,90,900,90,900,900",
                // A filtered genotype and a no-call, for the DP sum and the called-chromosome
                // count.
                "chr1\t300\t.\tA\tC\t50\t.\tDP=200;AC=1;AF=0.167;AN=6;MLEAC=1;MLEAF=0.167"
                        + "\tGT:GQ:DP:AD:PL:FT\t0/1:60:10:5,5:100,0,200:PASS\t0/0:60:20:20,0:0,60,600:LowGQ"
                        + "\t./.:.:30:.:.:PASS\t0/0:60:40:40,0:0,90,900:PASS",
                // An indel whose alternate disappears with s2 and s3, so removing it trims the
                // remaining alleles.
                "chr1\t400\t.\tACGT\tACG,A\t50\t.\tDP=200;AN=8"
                        + "\tGT:GQ:DP\t0/0:60:10\t0/0:60:20\t0/1:60:30\t0/2:60:40",
                // No DP in any genotype, so the record's own DP survives.
                "chr1\t500\t.\tA\tC\t50\t.\tDP=200;AN=8\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\t1/1:60");

        // The whole cohort, which is the untouched record.
        run(dir, "all-samples", input);
        // Two samples, then the same two with each of the flags that changes what is written.
        run(dir, "subset-two", input, "-sn", "s0", "-sn", "s1");
        run(dir, "subset-two-remove-unused", input, "-sn", "s0", "-sn", "s1",
                "--remove-unused-alternates", "true");
        run(dir, "subset-two-preserve", input, "-sn", "s0", "-sn", "s1",
                "--remove-unused-alternates", "true", "--preserve-alleles", "true");
        run(dir, "subset-two-keep-original", input, "-sn", "s0", "-sn", "s1",
                "--keep-original-ac", "true", "--keep-original-dp", "true");
        run(dir, "subset-two-remove-unused-keep-original", input, "-sn", "s0", "-sn", "s1",
                "--remove-unused-alternates", "true", "--keep-original-ac", "true");
        // One sample, and the pair that keeps every alternate alive.
        run(dir, "subset-one", input, "-sn", "s2");
        run(dir, "subset-carriers", input, "-sn", "s2", "-sn", "s3");
        // Every sample, with the flag that subsets alleles anyway.
        run(dir, "all-remove-unused", input, "--remove-unused-alternates", "true");
        // An exclusion, which selects everything else.
        run(dir, "exclude-one", input, "-xl-sn", "s3");
        // And the flag that drops a record whose alternates nobody calls.
        run(dir, "subset-two-exclude-non-variants", input, "-sn", "s0", "-sn", "s1",
                "--exclude-non-variants", "true");
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
        for (final String line : lines) {
            if (line.startsWith("#CHROM")) {
                final String[] field = line.split("\t", -1);
                final List<String> samples = new ArrayList<>();
                for (int i = 9; i < field.length; i++) {
                    samples.add(field[i]);
                }
                System.out.printf("samples\t%s\t%s%n", label, String.join(",", samples));
                continue;
            }
            if (line.startsWith("#")) {
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
