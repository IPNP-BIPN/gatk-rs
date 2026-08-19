/*
 * CollectReadCounts, taken from the reference.
 *
 * Read counts per interval, which is what the copy-number panel and the denoiser consume. A read
 * walker whose entire body is a lookup, and the lookup is the interesting part.
 *
 * Six behaviours this is built to catch.
 *
 *   - A READ IS COUNTED BY ITS START AND BY NOTHING ELSE. The overlap is computed against
 *     `SimpleInterval(contig, read.getStart(), read.getStart())`, so a read spanning three
 *     intervals is counted in ONE, and a read that starts before an interval and covers all of it
 *     is counted in NEITHER;
 *   - THE OVERLAP DETECTOR IS REBUILT PER CONTIG and holds only that contig's intervals, so a read
 *     can only ever match an interval on its own contig;
 *   - EVERY REQUESTED INTERVAL IS A ROW, in the order the interval argument produced, so an
 *     interval no read starts in is a row of zero rather than a gap;
 *   - THE FILTERS ARE THE WALKER'S PLUS FOUR, one of them a mapping quality of 30, so a read at 20
 *     is not counted anywhere;
 *   - THE DEFAULT OUTPUT IS HDF5 AND THE TSV IS OPT-IN. The TSV carries a SAM header, the sample as
 *     a read group, and three columns;
 *   - AND THE TOOL REQUIRES INTERVALS. Without `-L` it refuses before reading anything.
 *
 * Output:
 *
 *     table\t<label>\t<the whole TSV, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CollectReadCountsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectReadCountsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("collectreadcounts-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("counts.bam");
        buildFixture(bam.toFile());

        System.out.println("# CollectReadCountsDump: read counts per interval");

        // Four intervals of ten bases each. The reads start at 10, 12, 25, 40, 55 and 70.
        run("default", dir, fasta, bam,
                "-L", "chr1:10-19", "-L", "chr1:20-29", "-L", "chr1:30-39", "-L", "chr1:40-49");
        // An interval that a long read covers entirely but starts before, which counts nowhere.
        run("covered-not-started", dir, fasta, bam, "-L", "chr1:56-59");
        // An interval on the other contig, where no read is.
        run("other-contig", dir, fasta, bam, "-L", "chr2:1-100");
        // Intervals given out of order, which the argument layer sorts before the tool sees them.
        run("out-of-order", dir, fasta, bam, "-L", "chr1:40-49", "-L", "chr1:10-19");
        // The mapping quality filter, which this tool sets to 30.
        run("low-mapping-quality", dir, fasta, bam, "-L", "chr1:70-79");
        // No intervals at all.
        run("no-intervals", dir, fasta, bam);
    }

    /**
     * Six reads on chr1 and none on chr2.
     *
     *  - `r001` starts at 10 and `r002` at 12, both inside the first interval;
     *  - `r003` starts at 25, inside the second;
     *  - `r004` starts at 40, inside the fourth, and spans into the next interval, which it is not
     *    counted in;
     *  - `r005` starts at 55 and is forty bases long, so it COVERS 56-59 entirely without starting
     *    in it;
     *  - and `r006` starts at 70 at mapping quality 20, below this tool's threshold of 30.
     */
    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("SAMPLE");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(read(header, "r001", 10, 10, 60));
        records.add(read(header, "r002", 12, 10, 60));
        records.add(read(header, "r003", 25, 10, 60));
        records.add(read(header, "r004", 40, 20, 60));
        records.add(read(header, "r005", 55, 40, 60));
        records.add(read(header, "r006", 70, 10, 20));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final int length, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(length + "M");
        record.setMappingQuality(mappingQuality);
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < length; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        record.setReadBases(bases.toString().getBytes());
        final byte[] qualities = new byte[length];
        java.util.Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("counts-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-I", bam.toString(), "-O", out.toString(),
                "--format", "TSV",
                "--interval-merging-rule", "OVERLAPPING_ONLY"));
        argv.addAll(Arrays.asList(extra));
        try {
            new CollectReadCounts().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
