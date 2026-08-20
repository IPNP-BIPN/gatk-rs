/*
 * CallableLoci, taken from the reference.
 *
 * Every locus of the requested intervals is given one of six states, the runs of equal state are
 * written as BED, and the six counts are written as a summary.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE RUN TEST DOES NOT COMPARE CONTIGS. A run continues when the next locus starts at the
 *     previous end plus one and carries the same state, and nothing there asks whether the contig
 *     changed: two intervals on different contigs whose coordinates happen to run on come out as
 *     ONE BED line, under the first contig's name;
 *   - A DELETION COUNTS TOWARD THE QC DEPTH WHATEVER ITS BASE QUALITY, the test being
 *     `qual >= minBaseQuality || isDeletion()`, so a deleted position can be CALLABLE on reads
 *     whose bases would not have counted;
 *   - THE POOR-MAPPING-QUALITY TEST COMES BEFORE THE DEPTH TESTS and is a ratio at or above the
 *     fraction, so a locus can be POOR_MAPPING_QUALITY while its QC depth would have called it
 *     LOW_COVERAGE;
 *   - THE LOW-MAPQ COUNT IS `<=` AND THE PASSING COUNT IS `>=`, so a read at exactly the low
 *     threshold counts as low and a read at exactly the minimum counts as passing: with the
 *     defaults, both at once is impossible, but a run that sets them equal has reads counted in
 *     both;
 *   - EXCESSIVE_COVERAGE IS TESTED ON THE RAW DEPTH, not the QC depth, and only when a maximum was
 *     given;
 *   - AN N IN THE REFERENCE IS REF_N WHATEVER THE PILEUP IS, tested before any depth is counted;
 *   - THE BED IS ZERO-BASED ON ITS START AND ONE-BASED ON ITS END, and the summary is
 *     `%30s %d`, so both files are fixed-width in a way a port has to reproduce exactly;
 *   - AND THE TOOL REFUSES MORE THAN ONE SAMPLE, naming the samples it found.
 *
 * Output:
 *
 *     bed\t<label>=<the whole BED file, escaped>
 *     summary\t<label>=<the whole summary file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CallableLociDump
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
import org.broadinstitute.hellbender.tools.walkers.coverage.CallableLoci;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CallableLociDump {

    static final int CONTIG_LENGTH = 240;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("callable-loci-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        buildFixture(bam.toFile(), "sample");
        final Path twoSamples = dir.resolve("two-samples.bam");
        buildTwoSampleFixture(twoSamples.toFile());

        System.out.println("# CallableLociDump: the state every locus is given");

        // The plain run over the covered stretch, at the defaults.
        run("default", dir, fasta, bam, "-L", "chr1:1-60");
        // The N run, which is REF_N whatever covers it.
        run("n-run", dir, fasta, bam, "-L", "chr1:118-125");
        // The uncovered tail.
        run("no-coverage", dir, fasta, bam, "-L", "chr1:200-210");
        // Every locus one base at a time, which is the other output format.
        run("state-per-base", dir, fasta, bam, "-L", "chr1:1-40", "--format", "STATE_PER_BASE");
        // Two intervals on different contigs whose coordinates run on, which the run test merges.
        // Both stretches are REF_N and the second starts where the first ended plus one, so the
        // two contigs come out as ONE line under chr1's name.
        run("contig-run-on", dir, fasta, bam, "-L", "chr1:171-180", "-L", "chr2:181-190");
        // The same pair with different states, which the run test separates as it should.
        run("contig-run-on-differing", dir, fasta, bam, "-L", "chr1:200-210", "-L", "chr2:211-220");
        // A maximum depth, which is tested on the raw depth.
        run("max-depth", dir, fasta, bam, "-L", "chr1:1-60", "--max-depth", "3");
        // A minimum depth low enough that the low-quality reads still make it callable.
        run("min-depth-one", dir, fasta, bam, "-L", "chr1:1-60", "--min-depth", "1");
        // The mapping-quality thresholds set equal, so a read at the threshold counts in both.
        run("thresholds-equal", dir, fasta, bam, "-L", "chr1:1-60", "--max-low-mapq", "20",
                "--min-mapping-quality", "20");
        // The low-mapping-quality path, with a fraction low enough to reach it.
        run("poor-mapping-quality", dir, fasta, bam, "-L", "chr1:1-60",
                "--min-depth-for-low-mapq", "1", "--max-fraction-of-reads-with-low-mapq", "0.01");
        // A base-quality cutoff above every base, where only the deletions still count.
        run("high-base-quality", dir, fasta, bam, "-L", "chr1:1-60", "--min-base-quality", "60",
                "--min-depth", "1");
        // Two samples, which the tool refuses.
        run("two-samples", dir, fasta, twoSamples, "-L", "chr1:1-10");
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path bed = dir.resolve("out-" + label + ".bed");
        final Path summary = dir.resolve("summary-" + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-I", bam.toString(),
                "-O", bed.toString(),
                "--summary", summary.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new CallableLoci().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("bed\t%s=%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(bed))));
        System.out.printf("summary\t%s=%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(summary))));
    }

    /**
     * One sample over chr1:1-60 and the N run, with mapping and base qualities laid out so that
     * every state is reachable by moving one threshold.
     */
    static void buildFixture(final File bam, final String sample) {
        final SAMFileHeader header = header(new String[][] {{"rg1", sample}});
        final List<SAMRecord> records = new ArrayList<>();
        // Five reads of mapping quality 60 over 1-20, which is callable at the defaults.
        for (int i = 0; i < 5; i++) {
            records.add(read(header, "rg1", "high" + i, 1, "20M", (byte) 30, 60));
        }
        // Three reads of mapping quality 60 over 21-40, which is below the default minimum depth.
        for (int i = 0; i < 3; i++) {
            records.add(read(header, "rg1", "thin" + i, 21, "20M", (byte) 30, 60));
        }
        // Five reads over 41-60 whose mapping quality is 1, which is the default low threshold.
        for (int i = 0; i < 5; i++) {
            records.add(read(header, "rg1", "low" + i, 41, "20M", (byte) 30, 1));
        }
        // Two reads over 41-60 at mapping quality 20, so the low fraction is five sevenths.
        for (int i = 0; i < 2; i++) {
            records.add(read(header, "rg1", "mid" + i, 41, "20M", (byte) 30, 20));
        }
        // A read carrying a deletion over 5-8, which counts toward the QC depth whatever the
        // base quality cutoff is.
        records.add(read(header, "rg1", "del", 1, "4M4D12M", (byte) 5, 60));
        // The N run at 121-180, covered so that REF_N is not merely an absence of reads.
        for (int i = 0; i < 4; i++) {
            records.add(read(header, "rg1", "ns" + i, 115, "20M", (byte) 30, 60));
        }
        write(bam, header, records);
    }

    /** The same reads under two read groups of different samples. */
    static void buildTwoSampleFixture(final File bam) {
        final SAMFileHeader header = header(new String[][] {{"rg1", "first"}, {"rg2", "second"}});
        final List<SAMRecord> records = new ArrayList<>();
        records.add(read(header, "rg1", "a", 1, "20M", (byte) 30, 60));
        records.add(read(header, "rg2", "b", 1, "20M", (byte) 30, 60));
        write(bam, header, records);
    }

    static SAMFileHeader header(final String[][] groups) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : new String[] {"chr1", "chr2"}) {
            dictionary.addSequence(new SAMSequenceRecord(contig, CONTIG_LENGTH));
        }
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : groups) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample(group[1]);
            record.setPlatform("ILLUMINA");
            header.addReadGroup(record);
        }
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String group, final String name,
                          final int start, final String cigar, final byte quality,
                          final int mappingQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigar(TextCigarCodec.decode(cigar));
        int length = 0;
        for (final htsjdk.samtools.CigarElement element : record.getCigar().getCigarElements()) {
            if (element.getOperator().consumesReadBases()) {
                length += element.getLength();
            }
        }
        final byte[] bases = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] quals = new byte[length];
        Arrays.fill(quals, quality);
        record.setBaseQualities(quals);
        record.setMappingQuality(mappingQuality);
        record.setAttribute("RG", group);
        return record;
    }

    static void write(final File bam, final SAMFileHeader header, final List<SAMRecord> records) {
        records.sort((left, right) -> Integer.compare(left.getAlignmentStart(), right.getAlignmentStart()));
        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }
}
