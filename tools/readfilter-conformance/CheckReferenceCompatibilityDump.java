/*
 * CheckReferenceCompatibility's table, taken from the reference.
 *
 * A BAM or a VCF checked against several references. It leans on the same ReferenceSequenceTable as
 * CompareReferences, which is already ported, and adds a verdict per reference: the interesting
 * part is that the verdict is reached TWO different ways depending on one property of the input.
 *
 * Ten behaviours this is built to catch.
 *
 *   - WHETHER THE INPUT'S DICTIONARY CARRIES MD5s DECIDES THE WHOLE ALGORITHM: with them, the MD5
 *     table is built and the pair analysis is read; without, the dictionaries are compared by name
 *     and length alone;
 *   - AND A VCF CAN NEVER TAKE THE MD5 PATH. `VCFContigHeaderLine.getSAMSequenceRecord` copies the
 *     ID, the length and `assembly` and DROPS `M5`, so a header carrying M5 for every contig still
 *     produces a dictionary with none: every VCF run here lands on the name-and-length branch,
 *     including the one whose M5 is a lie;
 *   - AND `dictionaryHasMD5s` NEEDS EVERY SEQUENCE TO HAVE ONE, so one missing M5 sends the whole
 *     run down the other path;
 *   - THE TWO PATHS PRODUCE DIFFERENT SUMMARIES FOR THE SAME VERDICT, the name-and-length one
 *     saying outright that mismatching bases cannot be ruled out;
 *   - COMPATIBLE_SUBSET NEEDS THE STATUS SET TO BE EXACTLY {SUBSET}: any other flag beside it and
 *     the verdict is NOT_COMPATIBLE, whose summary quotes the EnumSet;
 *   - AND THE TWO PATHS DISAGREE ON WHAT A SUBSET IS: the MD5 path reads ReferencePair's SUBSET,
 *     the other reads SequenceDictionaryUtils' SUPERSET, which is the same relation named from the
 *     other side;
 *   - THE MISSING SEQUENCES ARE LISTED FROM THE REFERENCE'S SIDE, as a Java list's toString;
 *   - THE OUTPUT CARRIES A COMMENT LINE naming the input file before the header;
 *   - A BAM AND A VCF TOGETHER IS A BadInput, and so is neither;
 *   - AND A REFERENCE IS COMPARED AGAINST THE INPUT ONLY, never against another reference, so the
 *     table has one row per reference whatever they are to each other.
 *
 * Output:
 *
 *     fasta\t<label>=<the whole fasta, escaped>
 *     dict\t<label>=<the whole dictionary, escaped>
 *     input\t<label>=<the whole vcf, escaped>
 *     table\t<label>=<the output table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CheckReferenceCompatibilityDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.reference.CheckReferenceCompatibility;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CheckReferenceCompatibilityDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("check-reference-compatibility-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CheckReferenceCompatibilityDump: an input checked against references");

        final Path base = fasta(dir, "base",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n");
        final Path altered = fasta(dir, "altered",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTA\n");
        final Path extra = fasta(dir, "extra",
                ">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTTTTTTT\n>chr3\nGGGGGGGGGG\n");

        // A VCF whose header carries the base reference's contigs WITH their MD5s, and one that
        // carries them without.
        final String withMd5 = vcf(dir, "with-md5", true);
        final String withoutMd5 = vcf(dir, "without-md5", false);
        // And one whose contigs match by name and length but whose M5 is a lie.
        final String lying = vcf(dir, "lying", true).replaceAll(
                "M5=[0-9a-f]+", "M5=d41d8cd98f00b204e9800998ecf8427e");
        write(dir, "lying.vcf", lying);
        System.out.printf("input\tlying=%s%n", ReferenceQueryDump.escape(lying));

        run(dir, "md5-exact", "with-md5.vcf", List.of(base));
        run(dir, "md5-altered", "with-md5.vcf", List.of(altered));
        run(dir, "md5-subset", "with-md5.vcf", List.of(extra));
        // Every reference at once, which is one row each and no comparison between them.
        run(dir, "md5-all", "with-md5.vcf", List.of(base, altered, extra));
        // The same three without MD5s, which is the other algorithm entirely.
        run(dir, "no-md5-exact", "without-md5.vcf", List.of(base));
        run(dir, "no-md5-altered", "without-md5.vcf", List.of(altered));
        run(dir, "no-md5-subset", "without-md5.vcf", List.of(extra));
        // A lying M5, which the MD5 path believes.
        run(dir, "md5-lying", "lying.vcf", List.of(base));

        // A BAM, which is the only input whose dictionary can carry MD5s at all.
        final Path bamWithMd5 = dir.resolve("with-md5.bam");
        buildBam(bamWithMd5, true);
        final Path bamWithoutMd5 = dir.resolve("without-md5.bam");
        buildBam(bamWithoutMd5, false);
        runReads(dir, "bam-md5-exact", bamWithMd5, List.of(base));
        runReads(dir, "bam-md5-altered", bamWithMd5, List.of(altered));
        runReads(dir, "bam-md5-subset", bamWithMd5, List.of(extra));
        runReads(dir, "bam-md5-all", bamWithMd5, List.of(base, altered, extra));
        runReads(dir, "bam-no-md5", bamWithoutMd5, List.of(base));
        // A BAM and a VCF together, which the tool refuses outright.
        runBoth(dir, "both-inputs", bamWithMd5, "with-md5.vcf", List.of(base));
        // And no input at all.
        runWithoutInput(dir, "no-input", List.of(base));
    }

    /** A vcf whose contig lines come from the base reference, with or without their M5. */
    static String vcf(final Path dir, final String label, final boolean withMd5) throws Exception {
        final String chr1Md5 = withMd5 ? ",M5=a965a71aa3690f605935c54d320905ab" : "";
        final String chr2Md5 = withMd5 ? ",M5=820c922e12cfa860eb181d6269c77e63" : "";
        final String text =
                "##fileformat=VCFv4.2\n"
                + "##contig=<ID=chr1,length=20" + chr1Md5 + ">\n"
                + "##contig=<ID=chr2,length=10" + chr2Md5 + ">\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "chr1\t5\t.\tA\tC\t.\t.\t.\n";
        write(dir, label + ".vcf", text);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        return text;
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path file = dir.resolve(name);
        Files.writeString(file, text, StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        return file;
    }

    static Path fasta(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".fasta");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(file, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + file, "O=" + dir.resolve(label + ".dict")});
        System.out.printf("fasta\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        System.out.printf("dict\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(dir.resolve(label + ".dict")), dir)));
        return file;
    }

    static void run(final Path dir, final String label, final String vcf, final List<Path> references)
            throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList("-V", dir.resolve(vcf).toString()));
        finish(dir, label, argv, references);
    }

    /** A BAM whose header carries the base reference's contigs, with or without their M5. */
    static void buildBam(final Path file, final boolean withMd5) {
        final SAMSequenceRecord chr1 = new SAMSequenceRecord("chr1", 20);
        final SAMSequenceRecord chr2 = new SAMSequenceRecord("chr2", 10);
        if (withMd5) {
            chr1.setMd5("a965a71aa3690f605935c54d320905ab");
            chr2.setMd5("820c922e12cfa860eb181d6269c77e63");
        }
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(chr1, chr2)));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        try (SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, file.toFile())) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName("r1");
            record.setReferenceName("chr1");
            record.setAlignmentStart(1);
            record.setCigarString("5M");
            record.setReadBases("ACGTA".getBytes(StandardCharsets.UTF_8));
            record.setBaseQualities(new byte[] {30, 30, 30, 30, 30});
            record.setMappingQuality(60);
            writer.addAlignment(record);
        }
    }

    static void runReads(final Path dir, final String label, final Path bam,
                         final List<Path> references) throws Exception {
        finish(dir, label, new ArrayList<>(Arrays.asList("-I", bam.toString())), references);
    }

    static void runBoth(final Path dir, final String label, final Path bam, final String vcf,
                        final List<Path> references) throws Exception {
        finish(dir, label, new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-V", dir.resolve(vcf).toString())), references);
    }

    static void runWithoutInput(final Path dir, final String label, final List<Path> references)
            throws Exception {
        finish(dir, label, new ArrayList<>(), references);
    }

    static void finish(final Path dir, final String label, final List<String> argv,
                       final List<Path> references) throws Exception {
        final Path out = dir.resolve("table-" + label + ".tsv");
        argv.addAll(Arrays.asList("-O", out.toString()));
        for (final Path reference : references) {
            argv.addAll(Arrays.asList("-refcomp", reference.toString()));
        }
        try {
            new CheckReferenceCompatibility().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("table\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
