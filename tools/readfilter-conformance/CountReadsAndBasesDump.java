/*
 * CountReads and CountBases, taken from the reference.
 *
 * The two smallest tools in GATK, and the second and third members of the reference-and-coverage
 * archetype: their whole implementation is an `apply` that increments a counter. What is worth
 * measuring is therefore not the arithmetic but everything the engine does before `apply` is
 * reached, which is the same thing for both and is what the marginal cost of an archetype member
 * actually is.
 *
 * Five behaviours this is built to catch.
 *
 *   - THE DEFAULT FILTERS ARE NOT NOTHING. `WellformedReadFilter` runs on every read walker, so the
 *     fixture's read with no read group and its read with an N operator never reach apply and are
 *     not counted. A port that counted the file's records would answer a larger number;
 *   - COUNTBASES COUNTS getLength(), THE SEQUENCE, NOT THE SPAN. The fixture holds a read with a
 *     ten-base deletion, whose reference span is twenty and whose length is ten, and a read with an
 *     empty cigar, whose length is still its bases;
 *   - -L RESTRICTS BOTH, and an unmapped read parked at its mate's coordinate is inside the
 *     interval while an unplaced read is outside every one of them;
 *   - AN ADDED --read-filter IS ADDITIONAL, not a replacement, so NotDuplicateReadFilter removes
 *     the duplicate from a count the default filters had kept;
 *   - AND THE TWO TOOLS AGREE ON WHICH READS THEY SEE, which is what makes them one archetype: the
 *     same eight reads reach apply in both, and only the increment differs.
 *
 * Output:
 *
 *     reads\t<label>\t<what CountReads wrote, escaped>
 *     bases\t<label>\t<what CountBases wrote, escaped>
 *     error\t<label>\t<exception class>
 *
 * Usage: CountReadsAndBasesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.CountBases;
import org.broadinstitute.hellbender.tools.CountReads;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CountReadsAndBasesDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("countreads-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# CountReadsAndBasesDump: the two smallest tools in GATK");

        both("all", dir, bam, null);
        both("all-withref", dir, bam, fasta);
        both("chr1", dir, bam, null, "-L", "chr1");
        both("chr2", dir, bam, null, "-L", "chr2");
        // A window holding the deletion read, whose span is twice its length.
        both("deletion", dir, bam, null, "-L", "chr1:140-160");
        // The duplicate removed by a filter the default set does not include.
        both("no-duplicates", dir, bam, null, "--read-filter", "NotDuplicateReadFilter");
        // Every filter switched off, which is what counts the records the file actually holds.
        both("no-filters", dir, bam, null, "--disable-tool-default-read-filters");
        // An interval on a contig the dictionary does not have.
        both("unknown-contig", dir, bam, null, "-L", "chrZ");

        // A second fixture whose reads are not all the same length, which is what makes
        // CountBases's getLength() visible as something other than ten times CountReads.
        final Path varied = dir.resolve("varied.bam");
        buildVariedFixture(varied.toFile());
        both("varied", dir, varied, null);
        both("varied-nofilters", dir, varied, null, "--disable-tool-default-read-filters");
    }

    /**
     * Reads of fifteen, five and zero bases, plus one whose span is longer than its sequence.
     *
     * The zero-base read is the interesting one: `getLength()` is the SEQUENCE length, so it
     * contributes nothing to CountBases while still being one read to CountReads.
     */
    static void buildVariedFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", ReadWalkerDump.CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(varied(header, "v001", 10, "15M", "ACGTACGTACGTACG"));
        records.add(varied(header, "v002", 40, "5M", "ACGTA"));
        // Five bases either side of a ten-base deletion: a span of twenty, a length of ten.
        records.add(varied(header, "v003", 60, "5M10D5M", "ACGTACGTAC"));
        // No sequence at all, which is a read of length zero rather than a refusal.
        records.add(varied(header, "v004", 100, "*", ""));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord varied(final SAMFileHeader header, final String name, final int start,
                            final String cigar, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadBases(bases.isEmpty()
                ? SAMRecord.NULL_SEQUENCE : bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(bases.isEmpty() ? SAMRecord.NULL_QUALS : qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void both(final String label, final Path dir, final Path bam, final Path fasta,
                     final String... extra) throws Exception {
        run("reads", label, dir, bam, fasta, extra);
        run("bases", label, dir, bam, fasta, extra);
    }

    static void run(final String kind, final String label, final Path dir, final Path bam,
                    final Path fasta, final String... extra) throws Exception {
        final Path out = dir.resolve(kind + "-" + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-O", out.toString()));
        if (fasta != null) {
            argv.add("-R");
            argv.add(fasta.toString());
        }
        argv.addAll(Arrays.asList(extra));

        try {
            if (kind.equals("reads")) {
                new CountReads().instanceMain(argv.toArray(new String[0]));
            } else {
                new CountBases().instanceMain(argv.toArray(new String[0]));
            }
            System.out.printf("%s\t%s\t%s%n", kind, label,
                    ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s-%s\t%s%n", kind, label, e.getClass().getName());
        }
    }
}
