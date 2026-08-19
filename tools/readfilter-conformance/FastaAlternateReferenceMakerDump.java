/*
 * FastaAlternateReferenceMaker, taken from the reference.
 *
 * The sixth member of the reference archetype: FastaReferenceMaker with a VCF applied to it. Every
 * base of the traversal is either copied, replaced, dropped or expanded, and which one depends on
 * the variant at that locus and on state carried from the loci before it.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A DELETION KEEPS ITS FIRST BASE AND DROPS THE NEXT n, through a counter that survives from
 *     one apply to the next, so the bases removed are the ones AFTER the record's position;
 *   - AN INSERTION EMITS THE WHOLE ALTERNATE ALLELE, its anchor base included, so one locus can
 *     contribute several bases and the output is longer than the interval;
 *   - A FILTERED RECORD IS SKIPPED, and so is a record whose start is not this locus, which is what
 *     keeps a deletion's interior from being reconsidered;
 *   - THE FIRST CONCRETE ALTERNATE IS THE ONE USED, so a record whose first alternate is the
 *     spanning deletion `*` uses the second, and one with no concrete alternate emits nothing;
 *   - THE SNP MASK WRITES N, and --snp-mask-priority decides whether it beats a called SNP at the
 *     same site or loses to it;
 *   - --use-iupac-sample REPLACES A HET WITH ITS AMBIGUITY CODE, taken from that sample's genotype
 *     rather than from the alternate allele;
 *   - AND THE ARGUMENT CHECKS FIRE BEFORE THE TRAVERSAL: a priority flag with no mask is a
 *     CommandLineException and a sample the VCF does not have is a UserException.BadInput.
 *
 * Output:
 *
 *     fasta\t<label>\t<the FASTA text, escaped>
 *     fai\t<label>\t<its .fai, escaped>
 *     error\t<label>\t<exception class>
 *
 * Usage: FastaAlternateReferenceMakerDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.fasta.FastaAlternateReferenceMaker;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FastaAlternateReferenceMakerDump {

    /**
     * The variants, against chr1 of ReferenceQueryDump.FASTA:
     * `ACGTACGTACGTacgtNNNNacgtACGTRYKMSWBDHVNACGT`.
     *
     *  - 2 C>T, a plain SNP;
     *  - 5 A>AGG, a simple insertion;
     *  - 8 TAC>T, a simple deletion of the two bases AFTER position 8;
     *  - 15 N>A, filtered, and therefore not applied;
     *  - 20 N>*,C, whose first alternate is the spanning deletion;
     *  - 30 N>G, a het for sample NA1 and a hom var for NA2. The reference base there is an IUPAC
 *    code, which the reader answers as N, and a VCF may not carry one as a reference allele.
     */
    static final String VARIANTS =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##contig=<ID=chr1,length=43>\n"
            + "##contig=<ID=chr2,length=24>\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n"
            + "chr1\t2\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t1/1\n"
            + "chr1\t5\t.\tA\tAGG\t50\tPASS\t.\tGT\t0/1\t1/1\n"
            + "chr1\t8\t.\tTAC\tT\t50\tPASS\t.\tGT\t0/1\t1/1\n"
            + "chr1\t15\t.\tN\tA\t50\tLowQual\t.\tGT\t0/1\t1/1\n"
            + "chr1\t20\t.\tN\t*,C\t50\tPASS\t.\tGT\t0/1\t1/1\n"
            + "chr1\t30\t.\tN\tG\t50\tPASS\t.\tGT\t0/1\t1/1\n";

    /** The mask: one SNP at position 3, and one at 30 where the called SNP also sits. */
    static final String MASK =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=43>\n"
            + "##contig=<ID=chr2,length=24>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t3\t.\tG\tA\t50\tPASS\t.\n"
            + "chr1\t30\t.\tN\tT\t50\tPASS\t.\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("fastaalternate-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReferenceQueryDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path variants = dir.resolve("variants.vcf");
        Files.write(variants, VARIANTS.getBytes());
        final Path mask = dir.resolve("mask.vcf");
        Files.write(mask, MASK.getBytes());
        // A Feature file is queried by interval, so both need an index before the walker can ask.
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", variants.toString()});
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", mask.toString()});

        System.out.println("# FastaAlternateReferenceMakerDump: a reference with a VCF applied");

        // The whole of chr1, so every record is reached.
        run("plain", dir, fasta, variants, "-L", "chr1");
        // The mask, which loses to a called SNP at 30 unless it is given priority.
        run("masked", dir, fasta, variants, "-L", "chr1", "--snp-mask", mask.toString());
        run("mask-priority", dir, fasta, variants, "-L", "chr1",
                "--snp-mask", mask.toString(), "--snp-mask-priority", "true");
        // A het rendered as its IUPAC code, and a hom var rendered as its allele.
        run("iupac-het", dir, fasta, variants, "-L", "chr1", "--use-iupac-sample", "NA1");
        run("iupac-homvar", dir, fasta, variants, "-L", "chr1", "--use-iupac-sample", "NA2");
        // A window that starts inside the deletion, so the counter has nothing to carry.
        run("after-deletion", dir, fasta, variants, "-L", "chr1:9-20");
        // The two argument checks.
        run("priority-without-mask", dir, fasta, variants, "-L", "chr1",
                "--snp-mask-priority", "true");
        run("unknown-sample", dir, fasta, variants, "-L", "chr1", "--use-iupac-sample", "NOBODY");
    }

    static void run(final String label, final Path dir, final Path fasta, final Path variants,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".fasta");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-O", out.toString(), "-V", variants.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new FastaAlternateReferenceMaker().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("fasta\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
        final Path fai = dir.resolve("out-" + label + ".fasta.fai");
        System.out.printf("fai\t%s\t%s%n", label, Files.exists(fai)
                ? ReferenceQueryDump.escape(new String(Files.readAllBytes(fai))) : "(none)");
    }
}
