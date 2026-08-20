/*
 * ValidateBasicSomaticShortMutations, taken from the reference.
 *
 * A VariantWalker that asks, for each discovery call, whether a separate tumour-normal pair of bams
 * confirms it. The counting and the arithmetic under the tool are already pinned by
 * `allele-pileup-counter` and `somatic-validation-power`; this is the walker itself, the three
 * files it writes, and the places where it decides something the layers below do not.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE VALIDATION TOTAL IS REF PLUS ALT, NOT THE PILEUP DEPTH. Both counts come from the
 *     counter's map, so a read that is neither the reference nor the first alternate is in neither,
 *     and the depth the power is computed at is smaller than the pileup;
 *   - THE REFERENCE COUNT IS READ WITHOUT A DEFAULT while the alternate count has one, which is
 *     only safe because the counter always holds the reference key;
 *   - A NULL RESULT IS ADDED TO NOTHING AND THEN DEREFERENCED. The tool guards the list append
 *     with a null check and reads `getNumAltSupportingReadsInNormal()` off the same reference on the
 *     next line, so a genotype that is validatable but whose result is null throws;
 *   - THE NORMAL ARTIFACT TEST IS STRICTLY GREATER, so a count exactly at
 *     --max-validation-normal-count still validates;
 *   - AN ARTIFACT IS POWERED WHATEVER ITS POWER IS, `normalArtifact || power > minPower`, so an
 *     artifact with no power at all is still counted a false positive;
 *   - THE SUMMARY SPLITS ON `VariantContext.isSNP()`, which is the record's own type and not the
 *     genotype's, so a multiallelic SNP record counts as a SNP while the validation looked at one
 *     alternate;
 *   - A GENOTYPE THAT CANNOT BE VALIDATED IS WRITTEN TO THE ANNOTATED VCF AND NOWHERE ELSE, with
 *     JUDGMENT=SKIPPED and no POWER or VAL_AD;
 *   - A MISSING VALIDATION CONTROL SAMPLE SKIPS THE RECORD SILENTLY, writing nothing at all, not
 *     even a skipped judgment;
 *   - AND THE FILTER STRING IS THE RECORD'S FILTERS SORTED AND JOINED WITH SEMICOLONS, then
 *     concatenated with the genotype's own filters with no separator between the two.
 *
 * The symbolic-reference branch at the top of apply() is unreachable from a VCF: a record's
 * reference allele is bases by construction. It is left unmeasured for the same reason the
 * counter's own symbolic branch is.
 *
 * Output:
 *
 *     table\t<label>\t<the validation tsv, escaped>
 *     summary\t<label>\t<the concordance summary tsv, escaped>
 *     vcfline\t<label>\t<one line of the annotated vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ValidateBasicSomaticShortMutationsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.ValidateBasicSomaticShortMutations;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ValidateBasicSomaticShortMutationsDump {

    static final int CONTIG_LENGTH = 200;
    /** The discovery sample's name in the VCF, and the two sample names in the bam. */
    static final String DISCOVERY = "discovery";
    static final String CASE = "valcase";
    static final String CONTROL = "valcontrol";

    /** The reference is `ACGT` repeating, so a position's base is fixed by its offset. */
    static char referenceBase(final int position) {
        return "ACGT".charAt((position - 1) % 4);
    }

    static final String VCF_HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FILTER=<ID=strand_bias,Description=\"strand\">\n"
            + "##FILTER=<ID=weak_evidence,Description=\"weak\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
            + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##contig=<ID=chr2,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + DISCOVERY + "\n";

    /** The seven records of the main run, one per behaviour the walker decides. */
    static final String VARIANTS = VCF_HEADER
            // Plenty of validation alternate reads and a clean normal: validated.
            + "chr1\t20\t.\tT\tA\t50\tPASS\t.\tGT:AD\t0/1:40,10\n"
            // No validation alternate reads at all: unvalidated, and powered enough to be a false
            // positive.
            + "chr1\t24\t.\tT\tG\t50\tPASS\t.\tGT:AD\t0/1:40,10\n"
            // Three alternate reads in the validation normal, which is over the default maximum of
            // one: an artifact, and therefore powered whatever its power.
            + "chr1\t30\t.\tC\tA\t50\tPASS\t.\tGT:AD\t0/1:40,10\n"
            // An insertion.
            + "chr1\t34\t.\tC\tCGGG\t50\tPASS\t.\tGT:AD\t0/1:40,10\n"
            // A deletion.
            + "chr1\t40\t.\tTA\tT\t50\tPASS\t.\tGT:AD\t0/1:40,10\n"
            // Multiallelic: the AD has three entries, so the genotype cannot be validated at all.
            + "chr1\t44\t.\tT\tA,C\t50\tPASS\t.\tGT:AD\t0/1:30,10,5\n"
            // No AD: not validatable either, by the other half of the same test.
            + "chr1\t48\t.\tT\tA\t50\tPASS\t.\tGT\t0/1\n"
            // Filtered on the record and on the genotype, with exactly one alternate read in the
            // validation normal, which is not over the maximum.
            + "chr1\t56\t.\tT\tA\t50\tweak_evidence;strand_bias\t.\tGT:AD:FT\t0/1:40,10:base_qual\n";

    /** The one record whose result is null while its genotype is validatable. */
    static final String ZERO_AD = VCF_HEADER + "chr1\t20\t.\tT\tA\t50\tPASS\t.\tGT:AD\t0/1:0,0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("validate-basic-somatic-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("validation.bam");
        buildFixture(bam.toFile());

        final Path variants = write(dir, "variants.vcf", VARIANTS);
        final Path zeroAd = write(dir, "zero-ad.vcf", ZERO_AD);

        System.out.println("# ValidateBasicSomaticShortMutationsDump: the validation table a tool writes");

        run("default", dir, fasta, bam, variants);
        // Every unvalidated record is now powered, so the false-positive counts move.
        run("min-power-zero", dir, fasta, bam, variants, "--min-power", "0.0");
        // The artifact threshold raised past the three reads in the normal.
        run("normal-count-three", dir, fasta, bam, variants, "--max-validation-normal-count", "3");
        // The quality cutoff dropped to zero, which lets the low-quality alternate reads count.
        run("cutoff-zero", dir, fasta, bam, variants, "--min-base-quality-cutoff", "0");
        // A cutoff above every base in the fixture: every pileup is empty of usable reads.
        run("cutoff-fifty", dir, fasta, bam, variants, "--min-base-quality-cutoff", "50");
        // A control sample name the bam does not carry: every record is skipped in silence.
        runAs("missing-control", dir, fasta, bam, variants, CASE, "absent");
        // A case sample name the bam does not carry, which is a null pileup rather than a skip.
        runAs("missing-case", dir, fasta, bam, variants, "absent", CONTROL);
        // The validatable genotype whose result is null.
        run("zero-ad", dir, fasta, bam, zeroAd);
        // A window holding no records at all, which still writes a table and a summary.
        run("no-records", dir, fasta, bam, variants, "-L", "chr1:150-160");
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", path.toString()});
        return path;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final Path variants, final String... extra) throws Exception {
        runAs(label, dir, fasta, bam, variants, CASE, CONTROL, extra);
    }

    static void runAs(final String label, final Path dir, final Path fasta, final Path bam,
                      final Path variants, final String caseName, final String controlName,
                      final String... extra) throws Exception {
        final Path table = dir.resolve("table-" + label + ".tsv");
        final Path summary = dir.resolve("summary-" + label + ".tsv");
        final Path annotated = dir.resolve("annotated-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-I", bam.toString(),
                "-V", variants.toString(),
                "-O", table.toString(),
                "--" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_DISCOVERY_VCF_LONG_NAME, DISCOVERY,
                "--" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CASE, caseName,
                "--" + ValidateBasicSomaticShortMutations.SAMPLE_NAME_VALIDATION_CONTROL, controlName,
                "--" + ValidateBasicSomaticShortMutations.ANNOTATED_VCF_LONG_NAME, annotated.toString(),
                "--summary", summary.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new ValidateBasicSomaticShortMutations().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(table))));
        System.out.printf("summary\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(summary))));
        for (final String line : Files.readAllLines(annotated)) {
            // The command line carries a date and the tool's own version, neither of which a port
            // reproduces; the other header lines are compared as they are.
            if (line.startsWith("##GATKCommandLine")) {
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
    }

    /**
     * The validation pair: two samples over eight loci, plus the three reads the default filters
     * drop.
     */
    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : new String[] {"chr1", "chr2"}) {
            dictionary.addSequence(new SAMSequenceRecord(contig, CONTIG_LENGTH));
        }
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : new String[][] {{"rgcase", CASE}, {"rgcontrol", CONTROL}}) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample(group[1]);
            record.setPlatform("ILLUMINA");
            header.addReadGroup(record);
        }

        final List<SAMRecord> records = new ArrayList<>();
        // Position 20: validated, with four low-quality alternate reads the default cutoff drops.
        substitutions(records, header, CASE, 20, 'T', 22, 'A', 8, 30);
        substitutions(records, header, CASE, 20, 'T', 0, 'A', 4, 5);
        substitutions(records, header, CONTROL, 20, 'T', 20, 'A', 0, 30);
        // Position 24: no alternate reads in the validation tumour at all.
        substitutions(records, header, CASE, 24, 'T', 30, 'G', 0, 30);
        substitutions(records, header, CONTROL, 24, 'T', 20, 'G', 0, 30);
        // Position 30: three alternate reads in the validation normal, which is an artifact.
        substitutions(records, header, CASE, 30, 'C', 20, 'A', 10, 30);
        substitutions(records, header, CONTROL, 30, 'C', 20, 'A', 3, 30);
        // Position 34: an insertion of GGG.
        substitutions(records, header, CASE, 34, 'C', 20, 'C', 0, 30);
        insertions(records, header, CASE, 34, "GGG", 6);
        substitutions(records, header, CONTROL, 34, 'C', 20, 'C', 0, 30);
        // Position 40: a deletion of one base.
        substitutions(records, header, CASE, 40, 'T', 20, 'T', 0, 30);
        deletions(records, header, CASE, 40, 1, 6);
        substitutions(records, header, CONTROL, 40, 'T', 20, 'T', 0, 30);
        // Position 44: the multiallelic record's locus.
        substitutions(records, header, CASE, 44, 'T', 20, 'A', 5, 30);
        substitutions(records, header, CONTROL, 44, 'T', 20, 'A', 0, 30);
        // Position 48: the record with no AD.
        substitutions(records, header, CASE, 48, 'T', 10, 'A', 0, 30);
        substitutions(records, header, CONTROL, 48, 'T', 10, 'A', 0, 30);
        // Position 56: one alternate read in the normal, which is exactly the maximum and
        // therefore not an artifact.
        substitutions(records, header, CASE, 56, 'T', 20, 'A', 8, 30);
        substitutions(records, header, CONTROL, 56, 'T', 20, 'A', 1, 30);

        // The three reads the tool's own filters drop, all alternate at position 20: a duplicate,
        // a vendor-failed read and one with a mapping quality of zero.
        final SAMRecord duplicate = substitution(header, CASE, "dup", 20, 'A', 30);
        duplicate.setDuplicateReadFlag(true);
        records.add(duplicate);
        final SAMRecord vendor = substitution(header, CASE, "vendor", 20, 'A', 30);
        vendor.setReadFailsVendorQualityCheckFlag(true);
        records.add(vendor);
        final SAMRecord unmapped = substitution(header, CASE, "mapq0", 20, 'A', 30);
        unmapped.setMappingQuality(0);
        records.add(unmapped);

        records.sort((left, right) -> Integer.compare(left.getAlignmentStart(), right.getAlignmentStart()));
        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** The ten reference bases a 10M read starting five before the locus carries. */
    static String window(final int locus) {
        final StringBuilder bases = new StringBuilder();
        for (int position = locus - 5; position <= locus + 4; position++) {
            bases.append(referenceBase(position));
        }
        return bases.toString();
    }

    static void substitutions(final List<SAMRecord> records, final SAMFileHeader header,
                              final String sample, final int locus, final char referenceBase,
                              final int referenceCount, final char alternateBase,
                              final int alternateCount, final int quality) {
        for (int i = 0; i < referenceCount; i++) {
            records.add(substitution(header, sample, "r" + i, locus, referenceBase, quality));
        }
        for (int i = 0; i < alternateCount; i++) {
            records.add(substitution(header, sample, "a" + quality + "-" + i, locus, alternateBase,
                    quality));
        }
    }

    static SAMRecord substitution(final SAMFileHeader header, final String sample, final String name,
                                  final int locus, final char base, final int quality) {
        final StringBuilder bases = new StringBuilder(window(locus));
        bases.setCharAt(5, base);
        return record(header, sample, name, locus, "10M", bases.toString(), quality);
    }

    static void insertions(final List<SAMRecord> records, final SAMFileHeader header,
                           final String sample, final int locus, final String inserted,
                           final int count) {
        final String reference = window(locus);
        final String bases = reference.substring(0, 6) + inserted + reference.substring(6);
        for (int i = 0; i < count; i++) {
            records.add(record(header, sample, "i" + i, locus, "6M" + inserted.length() + "I4M",
                    bases, 30));
        }
    }

    static void deletions(final List<SAMRecord> records, final SAMFileHeader header,
                          final String sample, final int locus, final int length, final int count) {
        final StringBuilder bases = new StringBuilder();
        for (int position = locus - 5; position <= locus; position++) {
            bases.append(referenceBase(position));
        }
        for (int position = locus + 1 + length; position <= locus + 4 + length; position++) {
            bases.append(referenceBase(position));
        }
        for (int i = 0; i < count; i++) {
            records.add(record(header, sample, "d" + i, locus, "6M" + length + "D4M",
                    bases.toString(), 30));
        }
    }

    static SAMRecord record(final SAMFileHeader header, final String sample, final String name,
                            final int locus, final String cigar, final String bases,
                            final int quality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(sample + "-" + locus + "-" + name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(locus - 5);
        record.setCigar(TextCigarCodec.decode(cigar));
        record.setReadBases(bases.getBytes());
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) quality);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", sample.equals(CASE) ? "rgcase" : "rgcontrol");
        return record;
    }
}
