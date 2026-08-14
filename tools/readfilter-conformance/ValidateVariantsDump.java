/*
 * ValidateVariants' refusals, taken from the reference.
 *
 * The tool writes nothing: its whole output is whether it threw and what it said, which makes the
 * message the only thing there is to be identical about.
 *
 * Nine behaviours this is built to catch, three of which are about which check runs at all.
 *
 *   - THE DEFAULT IS "ALL", AND "ALL" WITHOUT A REFERENCE IS TWO CHECKS RATHER THAN FOUR: the
 *     alternate alleles and the chromosome counts. The reference-base check needs a reference and
 *     the ID check needs a dbSNP file, so a plain run silently tests less than the name says;
 *   - AND ASKING FOR REF WITHOUT A REFERENCE IS A MissingReference REFUSAL, thrown before any
 *     record is read;
 *   - --validation-type-to-exclude REMOVES ONE CHECK FROM THAT SET, so excluding CHR_COUNTS leaves
 *     the allele check alone and vice versa;
 *   - THE ALLELE CHECK IS ABOUT GENOTYPES, not about the ALT column: an alternate no sample calls
 *     is the failure, and the message names it;
 *   - THE COUNT CHECK COMPARES AC AND AN AGAINST THE GENOTYPES and reports the two numbers;
 *   - --validate-GVCF ADDS THREE CHECKS OF ITS OWN, not two: the <NON_REF> allele in every record,
 *     the records being ordered, and THE FILE COVERING THE WHOLE REFERENCE. The last one counts
 *     every uncovered locus and names the first gap, so a two-record GVCF over a 1900-base
 *     reference is refused for the 1898 loci it does not describe, whatever its records say;
 *   - AND THE RECORD-LEVEL CHECKS FIRE FIRST, since they run per record and the coverage one runs
 *     at the end of the traversal: a file that is both incomplete and missing <NON_REF> is refused
 *     for the allele;
 *   - --warn-on-errors TURNS EVERY REFUSAL INTO A LOG LINE, so the tool exits 0 on a file it would
 *     otherwise refuse and the caller learns nothing from the exit code;
 *   - --do-not-validate-filtered-records SKIPS A FILTERED RECORD ENTIRELY, so a broken record hides
 *     behind a FILTER column;
 *   - AND A RECORD WITH NO GENOTYPES AT ALL PASSES BOTH CHECKS, since there is nothing to disagree
 *     with.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     ok\t<label>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ValidateVariantsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import picard.sam.CreateSequenceDictionary;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ValidateVariantsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                    + "##FILTER=<ID=LowQD,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "##contig=<ID=chr2,length=900>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("validatevariants-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ValidateVariantsDump: what the validator refuses, from the reference");

        // A reference of two contigs of A, which is what the REF check reads and what every
        // validation type other than ALLELES and CHR_COUNTS needs to exist at all.
        final Path reference = writeReference(dir);

        // A file that passes both checks a plain run makes.
        final Path good = writeVcf(dir, "good",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=2;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t200\t.\tA\tC,G\t50\t.\tAC=1,1;AN=4;DP=30\tGT:GQ\t0/1:60\t0/2:60");

        // An alternate no genotype calls.
        final Path unusedAlternate = writeVcf(dir, "unused-alternate",
                "chr1\t100\t.\tA\tC,G\t50\t.\tAC=2,0;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60");

        // An AC that disagrees with the genotypes.
        final Path badCounts = writeVcf(dir, "bad-counts",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=4;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60");

        // The same, hidden behind a FILTER column.
        final Path filteredBadCounts = writeVcf(dir, "filtered-bad-counts",
                "chr1\t100\t.\tA\tC\t50\tLowQD\tAC=4;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60");

        // A record with no genotype columns at all.
        final Path sitesOnly = writeSitesOnly(dir, "sites-only",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=4;AN=4;DP=30");

        // A GVCF-shaped file, and one missing the <NON_REF> allele.
        final Path gvcf = writeVcf(dir, "gvcf",
                "chr1\t100\t.\tA\tC,<NON_REF>\t50\t.\tAC=2,0;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60",
                "chr1\t200\t.\tA\t<NON_REF>\t50\t.\tAN=4;DP=30\tGT:GQ\t0/0:60\t0/0:60");
        final Path notGvcf = writeVcf(dir, "not-gvcf",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=2;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60");

        // Records out of order within a contig, and the same start on two contigs.
        final Path unordered = writeUnsorted(dir, "unordered",
                "chr1\t300\t.\tA\t<NON_REF>\t50\t.\tAN=4;DP=30\tGT:GQ\t0/0:60\t0/0:60",
                "chr1\t100\t.\tA\t<NON_REF>\t50\t.\tAN=4;DP=30\tGT:GQ\t0/0:60\t0/0:60");
        final Path acrossContigs = writeVcf(dir, "across-contigs",
                "chr1\t300\t.\tA\t<NON_REF>\t50\t.\tAN=4;DP=30\tGT:GQ\t0/0:60\t0/0:60",
                "chr2\t100\t.\tA\t<NON_REF>\t50\t.\tAN=4;DP=30\tGT:GQ\t0/0:60\t0/0:60");

        // A record whose reference base is not the reference's.
        final Path wrongReferenceBase = writeVcf(dir, "wrong-reference-base",
                "chr1\t100\t.\tT\tC\t50\t.\tAC=2;AN=4;DP=30\tGT:GQ\t0/1:60\t0/1:60");

        // The baseline, and each check on its own.
        run(dir, "good", good);
        run(dir, "unused-alternate", unusedAlternate);
        run(dir, "unused-alternate-excluded", unusedAlternate,
                "--validation-type-to-exclude", "ALLELES");
        run(dir, "bad-counts", badCounts);
        run(dir, "bad-counts-excluded", badCounts, "--validation-type-to-exclude", "CHR_COUNTS");
        run(dir, "bad-counts-warn-only", badCounts, "--warn-on-errors", "true");
        run(dir, "filtered-bad-counts", filteredBadCounts);
        run(dir, "filtered-bad-counts-skipped", filteredBadCounts,
                "--do-not-validate-filtered-records", "true");
        run(dir, "sites-only", sitesOnly);

        // The reference-base check, which needs the reference every run now has.
        run(dir, "wrong-reference-base", wrongReferenceBase);
        run(dir, "wrong-reference-base-excluded", wrongReferenceBase,
                "--validation-type-to-exclude", "REF");

        // The GVCF checks.
        run(dir, "gvcf", gvcf, "--validate-GVCF", "true");
        run(dir, "not-gvcf", notGvcf, "--validate-GVCF", "true");
        run(dir, "unordered", unordered, "--validate-GVCF", "true");
        run(dir, "across-contigs", acrossContigs, "--validate-GVCF", "true");
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        return write(dir, label, true, records);
    }

    /** The same file, left unsorted, so the indexer is not asked to look at it. */
    static Path writeUnsorted(final Path dir, final String label, final String... records)
            throws Exception {
        return write(dir, label, false, records);
    }

    static Path write(final Path dir, final String label, final boolean index,
                      final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        if (index) {
            new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        }
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static Path writeSitesOnly(final Path dir, final String label, final String... records)
            throws Exception {
        final String header = HEADER.replace("\tFORMAT\ts0\ts1", "")
                .replace("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", "")
                .replace("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n", "");
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static Path reference;

    /** Two contigs of `A`, with the `.fai` and the `.dict` the tool asks for. */
    static Path writeReference(final Path dir) throws Exception {
        final StringBuilder text = new StringBuilder(">chr1\n");
        for (int line = 0; line < 10; line++) {
            text.append("A".repeat(100)).append("\n");
        }
        text.append(">chr2\n");
        for (int line = 0; line < 9; line++) {
            text.append("A".repeat(100)).append("\n");
        }
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        new CreateSequenceDictionary().instanceMain(new String[] {
            "-R", fasta.toString(), "-O", dir.resolve("reference.dict").toString()});
        reference = fasta;
        return fasta;
    }

    static void run(final Path dir, final String label, final Path input,
                    final String... arguments) {
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(),
                "-R", reference.toString()));
        all.addAll(List.of(arguments));
        try {
            new ValidateVariants().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("ok\t%s%n", label);
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
