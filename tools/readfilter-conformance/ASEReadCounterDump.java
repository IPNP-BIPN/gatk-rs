/*
 * ASEReadCounter, taken from the reference.
 *
 * A table of reference and alternate counts at heterozygous sites, and the interesting part is not
 * the counting but the ORDER of the tests that decide which reads are counted at all: five buckets
 * and four `continue`s, so a read that would fail two of them is only ever charged to the first.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE COUNTING CASCADE IS ORDERED. Improper pair, then low mapping quality, then low base
 *     quality, then other-base; each `continue`s, so a low-quality base on an improperly paired read
 *     counts as an improper pair and nothing else, and `rawDepth` is the only counter that sees
 *     every read;
 *   - OVERLAPPING MATES ARE ONE READ BY DEFAULT. `COUNT_FRAGMENTS_REQUIRE_SAME_BASE` keeps one
 *     element per read name and DISCARDS THE PAIR ENTIRELY when the two disagree, while
 *     `COUNT_FRAGMENTS` keeps the better-quality one either way and `COUNT_READS` keeps both;
 *   - A DISCARDED PAIR STAYS DISCARDED. The filter remembers the names it deleted, so a third
 *     element with that name later in the pileup does not resurrect it;
 *   - THE SITE FILTERS ARE THREE DIFFERENT ANSWERS: a non-biallelic site is skipped with a warning,
 *     a site with no het genotype is skipped with another, and TWO records at one position is a
 *     UserException that stops the run;
 *   - A SITE BELOW --min-depth-of-non-filtered-base PRODUCES NO LINE AT ALL, not a line of zeroes;
 *   - THE SEPARATOR IS THE ONLY DIFFERENCE BETWEEN TABLE AND RTABLE, and CSV is a comma;
 *   - AND THE HEADER IS PRINTED BEFORE THE TRAVERSAL, so a run that produces no rows still writes
 *     one line.
 *
 * Output:
 *
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ASEReadCounterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.ASEReadCounter;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ASEReadCounterDump {

    /**
     * Het sites against chr1 of ReadWalkerDump.FASTA, whose bases repeat `ACGT`.
     *
     * The reference base at each of these positions is `T`; the reads carry `G` there except for
     * `s001`, which carries the reference base, so both counters have something to count.
     *
     *  - 12 T>G het, where the overlapping pair agrees;
     *  - 16 T>G het, where the overlapping pair disagrees;
     *  - 20 T>G het, where one read's base quality is 20;
     *  - 24 T>G, hom var rather than het;
     *  - 28 T>A,C, triallelic.
     */
    static final String VARIANTS =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n"
            + "chr1\t12\trs1\tT\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t16\trs2\tT\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t20\trs3\tT\tG\t50\tPASS\t.\tGT\t0/1\n"
            + "chr1\t24\trs4\tT\tG\t50\tPASS\t.\tGT\t1/1\n"
            + "chr1\t28\trs5\tT\tA,C\t50\tPASS\t.\tGT\t0/1\n";

    /** A second file with a record at a position the first also carries. */
    static final String DUPLICATE_SITE =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n"
            + "chr1\t12\trs9\tG\tA\t50\tPASS\t.\tGT\t0/1\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("asereadcounter-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("ase.bam");
        buildFixture(bam.toFile());

        final Path variants = write(dir, "variants.vcf", VARIANTS);
        final Path duplicate = write(dir, "duplicate.vcf", DUPLICATE_SITE);

        System.out.println("# ASEReadCounterDump: reference and alternate counts at het sites");

        run("default", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40");
        run("count-reads", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--count-overlap-reads-handling", "COUNT_READS");
        run("count-fragments", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--count-overlap-reads-handling", "COUNT_FRAGMENTS");
        // The quality thresholds, which move reads from the counted buckets into their own.
        run("min-mapq", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--min-mapping-quality", "40");
        run("min-baseq", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--min-base-quality", "25");
        // A depth threshold that removes lines rather than zeroing them.
        run("min-depth", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--min-depth-of-non-filtered-base", "3");
        // The two other formats.
        run("csv", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--output-format", "CSV");
        run("table", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:1-40",
                "--output-format", "TABLE");
        // A window with no sites at all, which still writes the header.
        run("no-sites", dir, fasta, bam, "-V", variants.toString(), "-L", "chr1:100-110");
        // Two records at one position, which is a refusal rather than a skip.
        run("two-records", dir, fasta, bam, "-V", variants.toString(),
                "-V", duplicate.toString(), "-L", "chr1:12-12");
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", path.toString()});
        return path;
    }

    /**
     * Reads over chr1:10-29, built so that every branch of the cascade is taken by something.
     *
     * The pair `p001` overlaps at 12 with the same base and at 16 with different ones. `s001`
     * carries the REFERENCE base at 12 and is mapped at quality 30, below the
     * `min-mapping-quality` case's threshold, so that case moves it out of `refCount`. `s002` carries a
     * base quality of 20 at position 20, below the `min-base-quality` case's. `s003` is paired and
     * not properly paired, which is the first bucket of the cascade.
     */
    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("NA1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        // The overlapping pair: both cover 10-29, agreeing at 12 (`G`) and disagreeing at 16,
        // where the second mate carries `A` against the first's `G`.
        records.add(pair(header, "p001", 10, "ACGTACGTACGTACGTACGT", true, 60, 30));
        records.add(pair(header, "p001", 10, "ACGTACATACGTACGTACGT", false, 60, 30));
        // A single read at mapping quality 30.
        records.add(single(header, "s001", 10, "ACTTACGTACGTACGTACGT", 30, 30));
        // A single read whose base at position 20 is at quality 20.
        final SAMRecord lowBase = single(header, "s002", 10, "ACGTACGTACGTACGTACGT", 60, 30);
        final byte[] qualities = lowBase.getBaseQualities();
        qualities[10] = 20;
        lowBase.setBaseQualities(qualities);
        records.add(lowBase);
        // A paired read that is not properly paired.
        final SAMRecord improper = single(header, "s003", 10, "ACGTACGTACGTACGTACGT", 60, 30);
        improper.setReadPairedFlag(true);
        improper.setFirstOfPairFlag(true);
        improper.setProperPairFlag(false);
        improper.setMateReferenceName("chr1");
        improper.setMateAlignmentStart(10);
        records.add(improper);

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord single(final SAMFileHeader header, final String name, final int start,
                            final String bases, final int mappingQuality, final int baseQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(bases.length() + "M");
        record.setMappingQuality(mappingQuality);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static SAMRecord pair(final SAMFileHeader header, final String name, final int start,
                          final String bases, final boolean first, final int mappingQuality,
                          final int baseQuality) {
        final SAMRecord record = single(header, name, start, bases, mappingQuality, baseQuality);
        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);
        record.setFirstOfPairFlag(first);
        record.setSecondOfPairFlag(!first);
        record.setMateReferenceName("chr1");
        record.setMateAlignmentStart(start);
        record.setInferredInsertSize(first ? 20 : -20);
        return record;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("table-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-I", bam.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new ASEReadCounter().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
