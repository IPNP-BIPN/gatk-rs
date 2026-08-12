/*
 * BaseRecalibrator's output, taken from the reference.
 *
 * The tool that writes the table ApplyBQSR reads, and the one that closes the BQSR cycle. Its output
 * is a GATKReport of five tables, so this dump carries that whole text and every run is compared
 * character for character.
 *
 * Seven behaviours this is built to catch.
 *
 *   - ITS DEFAULT READ FILTERS ARE SEVEN, NOT ONE. `getStandardBQSRReadFilterList` is six BQSR
 *     filters plus WellformedReadFilter: mapping quality not zero, mapping quality available,
 *     mapped, not secondary, not duplicate, passes vendor quality. That is a fifth pattern of
 *     `getDefaultReadFilters` across the ported tools, and it is what decides which reads are
 *     counted at all;
 *   - --known-sites IS REQUIRED and takes any Feature file, so a BED and a VCF naming the same
 *     interval must produce the same table;
 *   - THE TWO ADDITIONAL COVARIATE TABLES SHARE ONE REPORT TABLE and one row counter. The context
 *     and cycle tables are written into a single RecalTable2, and `rowIndex` is not reset between
 *     them, which the reference's own XXX comment calls knowledge about the ordering of tables;
 *   - THE ROW KEYS ARE INTEGERS AND THE SORT IS SORT_BY_COLUMN, so the rows come out ordered by
 *     their VALUES and not by the order they were written;
 *   - THE EmpiricalQuality COLUMN COMPUTES AND CACHES, because `getEmpiricalQuality()` is called
 *     while the table is being written and the datum keeps the answer;
 *   - THE ARGUMENTS TABLE CARRIES THE COVARIATE CLASS NAMES, which is what RecalibrationReport
 *     checks when it reads the table back;
 *   - AND THE QUANTIZATION TABLE IS COMPUTED FROM THE FINAL TABLES, after finalizeData, so it
 *     depends on every read the traversal kept.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     reference\t<the FASTA text>
 *     sites\t<label>\t<the known sites file's text>
 *     table\t<label>\t<the whole recalibration table, escaped>
 *     error\t<label>\t<exception>\t<message>
 *
 * Usage: BaseRecalibratorDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BaseRecalibratorDump {

    /** Long enough that a ten-base read has room, and varied enough that mismatches are real. */
    static final String REFERENCE =
            "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGCTCTAGCATCGATCGATCGATTAGCTAGCTAGCTAACCGGTTACGT";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I is recorded inside the report's argument
        // table, so an absolute temporary path would make the output unstable.
        final Path dir = Path.of("baserecalibrator-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# BaseRecalibratorDump: BaseRecalibrator's output, from the reference");
        System.out.printf("reference\t%s%n", REFERENCE);

        final Path fasta = writeReference(dir);
        final Path bam = dir.resolve("input.bam");
        buildFixture(bam.toFile());
        System.out.printf("fixture\tinput\t%s%n", RecordTransformDump.base64(bam));
        System.out.printf("fixtureindex\tinput\t%s%n",
                RecordTransformDump.base64(dir.resolve("input.bai")));

        // A BED and a VCF naming the same interval, so the two must agree.
        final Path bed = dir.resolve("sites.bed");
        Files.writeString(bed, "chr1\t9\t12\n", StandardCharsets.UTF_8);
        System.out.printf("sites\tbed\t%s%n",
                ReferenceQueryDump.escape(Files.readString(bed)));
        final Path vcf = dir.resolve("sites.vcf");
        Files.writeString(vcf,
                "##fileformat=VCFv4.2\n"
                        + "##contig=<ID=chr1,length=" + REFERENCE.length() + ">\n"
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        + "chr1\t10\t.\tA\tC\t.\t.\t.\n"
                        + "chr1\t11\t.\tC\tG\t.\t.\t.\n"
                        + "chr1\t12\t.\tG\tT\t.\t.\t.\n",
                StandardCharsets.UTF_8);
        System.out.printf("sites\tvcf\t%s%n",
                ReferenceQueryDump.escape(Files.readString(vcf)));

        // A known-sites file must support random access, so both are indexed the way the tool's own
        // message tells the user to.
        for (final Path sites : new Path[] {bed, vcf}) {
            new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                    .instanceMain(new String[] {"-I", sites.toString()});
        }

        run(dir, bam, fasta, bed, "bed-sites", new String[] {});
        run(dir, bam, fasta, vcf, "vcf-sites", new String[] {});
        // The indel tables, which the hidden argument turns on.
        run(dir, bam, fasta, bed, "indel-tables",
                new String[] {"--compute-indel-bqsr-tables", "true"});
        // BAQ, which is off by default.
        run(dir, bam, fasta, bed, "baq-enabled", new String[] {"--enable-baq", "true"});
        // A different quantization level count, which changes only the Quantized table.
        run(dir, bam, fasta, bed, "quantizing-4", new String[] {"--quantizing-levels", "4"});
        // A different preserve threshold, which changes which bases are counted.
        run(dir, bam, fasta, bed, "preserve-20",
                new String[] {"--preserve-qscores-less-than", "20"});
        // A larger context, which widens the context table's keys.
        run(dir, bam, fasta, bed, "context-3",
                new String[] {"--mismatches-context-size", "3"});
    }

    /**
     * Reads shaped so the seven default filters have something to drop, and so the counting loop
     * sees mismatches, an insertion and a deletion.
     */
    static void buildFixture(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Kept: plain matches and mismatches.
            writer.addAlignment(read(header, "match", 1, "10M", "ACGTACGTAC", 0, 60));
            writer.addAlignment(read(header, "mismatch", 1, "10M", "ACGTTCGTAC", 0, 60));
            writer.addAlignment(read(header, "deletion", 5, "4M2D6M", "ACGTACGTAC", 0, 60));
            // Over the known sites, so the skip array has something to remove.
            writer.addAlignment(read(header, "over-sites", 8, "10M", "ACGTACGTAC", 0, 60));
            writer.addAlignment(read(header, "insertion", 9, "4M2I4M", "ACGTACGTAC", 16, 60));
            // Dropped, one per filter, so the default list is visible in what the table does not
            // hold.
            writer.addAlignment(read(header, "mapq-zero", 20, "10M", "ACGTACGTAC", 0, 0));
            writer.addAlignment(read(header, "mapq-unavailable", 24, "10M", "ACGTACGTAC", 0, 255));
            writer.addAlignment(read(header, "secondary", 28, "10M", "ACGTACGTAC", 0x100, 60));
            writer.addAlignment(read(header, "duplicate", 32, "10M", "ACGTACGTAC", 0x400, 60));
            writer.addAlignment(read(header, "vendor-fail", 36, "10M", "ACGTACGTAC", 0x200, 60));
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", REFERENCE.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        group.setPlatform("ILLUMINA");
        group.setPlatformUnit("unit-rg1");
        header.addReadGroup(group);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases, final int flags,
                          final int mapq) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        for (int i = 0; i < quals.length; i++) {
            // A gradient across the usable threshold, so the skip array has both answers.
            quals[i] = (byte) (2 + i * 4);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(mapq);
        record.setAttribute("RG", "rg1");
        return record;
    }

    /** One run of the tool, with the whole recalibration table it wrote. */
    static void run(final Path dir, final Path bam, final Path fasta, final Path sites,
                    final String label, final String[] extra) throws Exception {
        final Path output = dir.resolve("recal." + label + ".table");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-R", fasta.toString(),
                "--known-sites", sites.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new BaseRecalibrator().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(output)));
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + REFERENCE + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + REFERENCE.length()
                + "\n", StandardCharsets.UTF_8);
        return fasta;
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.delete(entry);
            }
        }
    }
}
