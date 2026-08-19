/*
 * CollectAllelicCounts, taken from the reference.
 *
 * The counts ModelSegments and CalculateContamination's copy-number cousins read. A locus walker
 * with four extra filters, a base quality threshold applied INSIDE the collector rather than by a
 * filter, and a table whose header is a whole SAM header.
 *
 * Seven behaviours this is built to catch.
 *
 *   - EVERY LOCUS IN THE INTERVAL IS A ROW. `emitEmptyLoci` is true, so a position with no reads at
 *     all is a row of zeroes rather than a gap, and the table's length is the interval's;
 *   - THE ALTERNATE COUNT IS TOTAL MINUS REFERENCE, not the count of the alternate base. The
 *     reference's own comment says so. A locus with three different non-reference bases has all
 *     three in the alternate count;
 *   - THE ALTERNATE BASE IS THE MOST COMMON NON-REFERENCE ONE, chosen by a sort that is NOT stable
 *     against ties in any documented way, and it is `N` when the alternate count is zero;
 *   - THE TOTAL IS OVER ACGT ONLY, so an `N` in the pileup is in neither count and a deletion is
 *     filtered before the counting;
 *   - THE BASE QUALITY THRESHOLD IS THE COLLECTOR'S, not a read filter's, and its default is 20 --
 *     so a base at quality 19 is dropped while its read is kept, and `--minimum-base-quality 0`
 *     brings it back;
 *   - A REFERENCE BASE THAT IS NOT ACGT SKIPS THE LOCUS ENTIRELY, so an `N` in the reference
 *     produces no row at all even though `emitEmptyLoci` is true;
 *   - AND THE TABLE'S HEADER IS A SAM HEADER, `@HD`, `@SQ` and `@RG` lines and all, followed by the
 *     column names. The sample name comes from the read group.
 *
 * Output:
 *
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CollectAllelicCountsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.CollectAllelicCounts;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectAllelicCountsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("collectallelic-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("allelic.bam");
        buildFixture(bam.toFile());

        System.out.println("# CollectAllelicCountsDump: reference and alternate counts per locus");

        // Ten loci, of which the last four have no reads at all.
        run("default", dir, fasta, bam, "-L", "chr1:10-19");
        // The base quality threshold, which is the collector's rather than a filter's.
        run("base-quality-zero", dir, fasta, bam, "-L", "chr1:10-19",
                "--minimum-base-quality", "0");
        run("base-quality-thirty", dir, fasta, bam, "-L", "chr1:10-19",
                "--minimum-base-quality", "30");
        // The N run in the reference, where the locus is skipped rather than emitted empty.
        run("reference-n", dir, fasta, bam, "-L", "chr1:123-130");
        // A window with no reads at all, which is all zeroes and an N alternate.
        run("no-reads", dir, fasta, bam, "-L", "chr1:60-64");
        // The mapping quality filter this tool adds, at its default of 30.
        run("low-mapping-quality", dir, fasta, bam, "-L", "chr1:30-33");
    }

    /**
     * Reads over chr1:10-19 and chr1:30-33, whose reference bases repeat `ACGT`.
     *
     *  - `r001` and `r002` carry the reference bases;
     *  - `r003` carries `A` at every position, so it is the alternate wherever the reference is not
     *    `A`;
     *  - `r004` carries `C` at every position, giving a second non-reference base;
     *  - `r005` carries `N` at every position, which the total does not count;
     *  - `r006` has one base at quality 19, just under the default threshold;
     *  - and `r007` covers 30-33 at mapping quality 20, below this tool's threshold of 30.
     */
    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("NA1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(read(header, "r001", 10, "ACGTAC", 60, 30, -1));
        records.add(read(header, "r002", 10, "ACGTAC", 60, 30, -1));
        records.add(read(header, "r003", 10, "AAAAAA", 60, 30, -1));
        records.add(read(header, "r004", 10, "CCCCCC", 60, 30, -1));
        records.add(read(header, "r005", 10, "NNNNNN", 60, 30, -1));
        // One base at quality 19, at the first position.
        records.add(read(header, "r006", 10, "ACGTAC", 60, 30, 0));
        records.add(read(header, "r007", 30, "ACGT", 20, 30, -1));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String bases, final int mappingQuality, final int baseQuality,
                          final int lowQualityIndex) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(bases.length() + "M");
        record.setMappingQuality(mappingQuality);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) baseQuality);
        if (lowQualityIndex >= 0) {
            qualities[lowQualityIndex] = 19;
        }
        record.setBaseQualities(qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("counts-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-I", bam.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new CollectAllelicCounts().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
