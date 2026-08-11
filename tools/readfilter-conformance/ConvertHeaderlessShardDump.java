/*
 * ConvertHeaderlessHadoopBamShardToBam, taken from the reference.
 *
 * The thirteenth whole tool of the record-transform archetype, and the FIRST THAT IS NOT A GATKTool
 * AT ALL: it extends CommandLineProgram and implements doWork(). No traversal, no read filter, no
 * @PG record, no engine. Three files in, one file out, and the middle one is copied rather than
 * read.
 *
 * Five behaviours this is built to catch.
 *
 *   - THE SHARD IS COPIED BYTE FOR BYTE. `FileUtils.copyFile(bamShard, outStream)` between a header
 *     block and a terminator, so the output's data blocks are the shard's own bytes: not
 *     decompressed, not re-deflated, not re-blocked. A port that read the records and wrote them
 *     back would produce a valid BAM with different bytes, which is the failure this dump exists to
 *     make visible;
 *   - THE HEADER IS ENCODED WITH keepExistingVersionNumber = true, which is the OTHER branch from
 *     the one PrintReadsHeader takes and the one htsjdk-rs#164 says the ordinary BAM writer does not
 *     take. This is the block-copy reheader path, and here `true` is correct;
 *   - THE HEADER BLOCK CARRIES NO TERMINATOR. `writeBAMHeaderToStream` calls
 *     `blockCompressedOutputStream.flush()` rather than closing it, so the empty gzip block appears
 *     exactly once, at the very end, after the copied shard;
 *   - THE SEQUENCE DICTIONARY IS WRITTEN TWICE, once inside the header text and once as the binary
 *     block after it, and the second is written from the same header rather than parsed back from
 *     the first;
 *   - AND NOTHING IS APPENDED. No @PG, because there is no GATKTool to append one, so the output's
 *     header is the donor's header unchanged apart from the version rule above.
 *
 * The shard the dump feeds in is the data blocks of a real BAM with its header block removed, which
 * is what a Spark tool with --sharded-output produces.
 *
 * Output:
 *
 *     deflater\t<class>
 *     donor\t<base64 bam whose header is used>
 *     shard\t<base64 headerless shard>
 *     shardlength\t<bytes>
 *     output\t<label>\t<base64 bam>
 *     outputheader\t<label>\t<escaped SAM header read back>
 *     reads\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<cigar>
 *     copiedverbatim\t<label>\t<true|false>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: ConvertHeaderlessShardDump
 */

import htsjdk.samtools.BAMRecordCodec;
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
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.zip.DeflaterFactory;

import org.broadinstitute.hellbender.tools.ConvertHeaderlessHadoopBamShardToBam;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class ConvertHeaderlessShardDump {

    public static void main(final String[] args) throws Exception {
        // The factory is static and the first writer wins. This dump calls no Picard entry point,
        // so nothing should replace it; the pin makes that a fact rather than a hope.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose, as in the other dumps of this archetype.
        final Path dir = Path.of("shard-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConvertHeaderlessShardDump: ConvertHeaderlessHadoopBamShardToBam");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // A complete BAM, which is both the header donor and the source of the shard.
        final Path donor = dir.resolve("donor.bam");
        buildDonor(donor.toFile());
        System.out.printf("donor\t%s%n", base64(donor));

        // The shard: the donor's records encoded through BAMRecordCodec into a BGZF stream with no
        // header block and no terminator, which is what a Spark tool writes with
        // --sharded-output true. Slicing a finished BAM would not do: htsjdk packs the header and
        // small record data into ONE block, so there is no boundary between them to cut at.
        final Path shard = dir.resolve("shard.bam");
        writeShard(donor, shard);
        final byte[] shardBytes = Files.readAllBytes(shard);
        System.out.printf("shard\t%s%n", base64(shard));
        System.out.printf("shardlength\t%d%n", shardBytes.length);

        convert(dir, shard, donor, "plain");
        // A shard of nothing at all: the output is a header block and a terminator.
        final Path empty = dir.resolve("empty_shard.bam");
        Files.write(empty, new byte[0]);
        convert(dir, empty, donor, "emptyshard");
        // A donor that is not a BAM at all.
        convert(dir, shard, shard, "badheader");

        // Whether the shard's bytes survive verbatim inside each output, which is the claim.
        verbatim(dir, "plain", shardBytes);
        verbatim(dir, "emptyshard", new byte[0]);

        reads(dir, "plain");
    }

    /**
     * A headerless shard: the donor's records, BGZF-compressed, with no header block in front and
     * no terminator behind.
     *
     * `flush()` rather than `close()`, for the same reason the tool's own header writer flushes:
     * closing would append the empty gzip block that marks the end of a BAM, and a shard is by
     * definition not the end of one.
     */
    static void writeShard(final Path donor, final Path shard) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(donor.toFile())) {
            final BlockCompressedOutputStream out =
                    new BlockCompressedOutputStream(shard.toFile());
            final BAMRecordCodec codec = new BAMRecordCodec(reader.getFileHeader());
            codec.setOutputStream(out);
            for (final SAMRecord record : reader) {
                codec.encode(record);
            }
            out.flush();
        }
    }

    /** The total length of the BGZF block starting at `offset`, from its BSIZE field. */
    static int blockLength(final byte[] data, final int offset) {
        // The BSIZE subfield sits at bytes 16-17 of a BGZF block and holds total length minus one.
        final int bsize = (data[offset + 16] & 0xFF) | ((data[offset + 17] & 0xFF) << 8);
        return bsize + 1;
    }

    static void buildDonor(final File file) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000), new SAMSequenceRecord("chr2", 2000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        existing.setCommandLine("upstream --in a.bam");
        header.addProgramRecord(existing);
        header.addComment("a donor comment");

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, file)) {
            for (int i = 0; i < 4; i++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("r" + i);
                record.setReferenceName("chr1");
                record.setAlignmentStart(100 + i * 100);
                record.setCigarString("10M");
                record.setReadBases("ACGTACGTAC".getBytes());
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    static void convert(final Path dir, final Path shard, final Path header, final String label)
            throws Exception {
        final Path output = dir.resolve("Converted." + label + ".bam");
        // --use-jdk-deflater is on CommandLineProgram, not on GATKTool, so it reaches this tool
        // too, and it has to be passed: CommandLineProgram installs IntelDeflaterFactory into the
        // static BlockCompressedOutputStream default when it is absent, whatever this dump pinned
        // before the run. Measured: without it the header block is 177 bytes and the JDK deflater
        // produces 173 at every level, which is #160's mechanism seen from the other side.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--bam-shard", shard.toString(),
                "--bam-with-header", header.toString(),
                "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));

        try {
            new ConvertHeaderlessHadoopBamShardToBam().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
            return;
        }

        // What the factory was when the tool actually wrote, rather than what the dump set.
        System.out.printf("deflaterinrun\t%s\t%s%n", label,
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        System.out.printf("output\t%s\t%s%n", label, base64(output));
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            System.out.printf("outputheader\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(reader.getFileHeader().getSAMString()));
        }
    }

    /**
     * Whether the shard's bytes appear verbatim in the output, and where.
     *
     * The claim is not "the reads round-trip" but "the bytes are copied", and only a byte search
     * can tell the two apart.
     */
    static void verbatim(final Path dir, final String label, final byte[] shard) throws Exception {
        final Path output = dir.resolve("Converted." + label + ".bam");
        if (!Files.exists(output)) {
            return;
        }
        final byte[] bytes = Files.readAllBytes(output);
        final int terminator = BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length;
        final int headerBlock = blockLength(bytes, 0);
        final byte[] middle = Arrays.copyOfRange(bytes, headerBlock, bytes.length - terminator);
        System.out.printf("copiedverbatim\t%s\t%b%n", label, Arrays.equals(middle, shard));
        System.out.printf("layout\t%s\t%d\t%d\t%d%n", label, headerBlock, middle.length,
                terminator);
    }

    static void reads(final Path dir, final String label) throws Exception {
        final Path output = dir.resolve("Converted." + label + ".bam");
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            for (final SAMRecord record : reader) {
                System.out.printf("reads\t%s\t%s\t%d\t%s\t%d\t%s%n", label, record.getReadName(),
                        record.getFlags(), record.getReferenceName(), record.getAlignmentStart(),
                        record.getCigarString());
            }
        }
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }
}
