/*
 * ReadAnonymizer's output, taken from the reference.
 *
 * The ninth whole tool of the record-transform archetype, and the first whose transform CHANGES THE
 * LENGTH OF THE READ. Every base that disagrees with the reference is replaced by the reference base,
 * and the cigar is rewritten to say so.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A DELETION ADDS BASES TO THE READ. `case X: case D:` puts the reference bases in, so a `4M2D6M`
 *     read comes out twelve bases long from ten. No other record-transform tool ported so far makes
 *     the read longer;
 *   - AN INSERTION REMOVES THEM, and its cigar contribution is `currentNewCigarOp` with a count of
 *     ZERO, so it merges into whatever came before rather than ending that element;
 *   - EVERY M AND X AND D BECOMES ONE OPERATOR, `=` by default and `M` under --use-simple-cigar, so
 *     consecutive elements of different kinds collapse into one;
 *   - AND THE LAST ELEMENT IS ADDED UNCONDITIONALLY, so a read whose cigar ends in an insertion
 *     emits whatever the accumulator held, which can be a zero-length element;
 *   - THE QUALITY OF A REPLACED BASE BECOMES --ref-base-quality AND THE OTHERS ARE KEPT, so a run of
 *     matching bases keeps its gradient and a mismatch does not;
 *   - EVERY ATTRIBUTE IS CLEARED EXCEPT THE READ GROUP, so a read's tags do not survive;
 *   - AND ITS WRITER IS NOT PRESORTED. `createSAMWriter(output, false)` says the reads are NOT in
 *     order, where every other record-transform tool ported so far passes true. An index is still
 *     written, because `--create-output-bam-index` is a separate argument whose default is true, so
 *     the flag decides the sorting and not the index. The header still says `SO:coordinate`.
 *
 * Its default read filters are a sixth pattern: seven filters and NOT WellformedReadFilter, namely
 * valid alignment start, valid alignment end, read length equals cigar length, sequence is stored,
 * matching bases and quals, mapped, and alignment agrees with the header.
 *
 * Output, one row per (label, kind):
 *
 *     reference\t<the reference bases>
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     header\t<label>\t<the output header, escaped>
 *     commandline\t<label>\t<the @PG command line>
 *     output\t<label>\t<the output BAM, base64>
 *     index\t<label>\t<the index, base64, or absent>
 *     error\t<label>\t<exception>\t<message>
 *
 * Usage: ReadAnonymizerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.walkers.ReadAnonymizer;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ReadAnonymizerDump {

    static final String REFERENCE =
            "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGCTCTAGCATCGATCGATCGATTAGCTAGCTAGCTAACCGGTTACGT";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("readanonymizer-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ReadAnonymizerDump: ReadAnonymizer's output, from the reference");
        System.out.printf("reference\t%s%n", REFERENCE);

        final Path fasta = writeReference(dir);
        final Path bam = dir.resolve("input.bam");
        buildFixture(bam.toFile());
        System.out.printf("fixture\tinput\t%s%n", RecordTransformDump.base64(bam));
        System.out.printf("fixtureindex\tinput\t%s%n",
                RecordTransformDump.base64(dir.resolve("input.bai")));

        run(dir, bam, fasta, "plain", new String[] {});
        run(dir, bam, fasta, "simple-cigar", new String[] {"--use-simple-cigar", "true"});
        run(dir, bam, fasta, "ref-qual-0", new String[] {"--ref-base-quality", "0"});
        run(dir, bam, fasta, "ref-qual-60", new String[] {"--ref-base-quality", "60"});
        // Above the declared maximum, which the argument parser refuses.
        run(dir, bam, fasta, "ref-qual-61", new String[] {"--ref-base-quality", "61"});
        run(dir, bam, fasta, "no-index",
                new String[] {"--create-output-bam-index", "false"});
        run(dir, bam, fasta, "no-program-record",
                new String[] {"--add-output-sam-program-record", "false"});
    }

    /**
     * Reads whose cigars make the transform change the read's length in both directions, plus one
     * per default filter so the list is visible in what the output does not hold.
     */
    static void buildFixture(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Every base matches: nothing is replaced and every quality is kept.
            writer.addAlignment(read(header, "exact", 1, "10M", "ACGTACGTAC", 0, 60));
            // One mismatch, so one quality becomes the reference quality and the rest do not.
            writer.addAlignment(read(header, "mismatch", 1, "10M", "ACGTTCGTAC", 0, 60));
            // A deletion, which ADDS two bases to the read.
            writer.addAlignment(read(header, "deletion", 3, "4M2D6M", "ACGTACGTAC", 0, 60));
            // An insertion, which removes two and contributes a zero-length cigar element.
            writer.addAlignment(read(header, "insertion", 5, "4M2I4M", "ACGTACGTAC", 0, 60));
            // A soft clip, whose bases and qualities are kept as they are.
            writer.addAlignment(read(header, "soft-clipped", 7, "3S7M", "ACGTACGTAC", 0, 60));
            // An X, which the transform treats exactly as a D would be treated.
            writer.addAlignment(read(header, "explicit-x", 9, "4M2X4M", "ACGTACGTAC", 0, 60));
            // Already `=`, which is kept rather than rewritten.
            writer.addAlignment(read(header, "explicit-eq", 11, "10=", "GGGGCCCCAA", 0, 60));
            // A cigar ending in an insertion, which is where the unconditional last element shows.
            writer.addAlignment(read(header, "trailing-insertion", 13, "8M2I", "ACGTACGTAC", 0, 60));
            // A cigar starting with one, where the accumulator is still empty.
            writer.addAlignment(read(header, "leading-insertion", 15, "2I8M", "ACGTACGTAC", 0, 60));
            // Tags that the transform clears, all but the read group.
            final SAMRecord tagged = read(header, "tagged", 17, "10M", "ACGTACGTAC", 0, 60);
            tagged.setAttribute("NM", 3);
            tagged.setAttribute("OQ", "!!!!!!!!!!");
            writer.addAlignment(tagged);
            // Dropped by the filters: unmapped, and a cigar that does not match the read's length.
            writer.addAlignment(read(header, "unmapped", 19, "10M", "ACGTACGTAC", 4, 60));
            writer.addAlignment(read(header, "length-mismatch", 21, "5M", "ACGTACGTAC", 0, 60));
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
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases, final int flags, final int mapq) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        for (int i = 0; i < quals.length; i++) {
            // A gradient, so a kept quality is distinguishable from a replaced one.
            quals[i] = (byte) (10 + i * 3);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(mapq);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void run(final Path dir, final Path input, final Path fasta, final String label,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve("ReadAnonymizer." + label + ".bam");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-R", fasta.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new ReadAnonymizer().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }

        String commandLine = "";
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(output.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMProgramRecord record : header.getProgramRecords()) {
                if (record.getCommandLine() != null) {
                    commandLine = record.getCommandLine();
                }
            }
            System.out.printf("header\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(header.getSAMString()));
        }
        System.out.printf("commandline\t%s\t%s%n", label, commandLine);
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
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
