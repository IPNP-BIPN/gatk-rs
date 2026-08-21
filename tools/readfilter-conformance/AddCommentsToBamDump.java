/*
 * AddCommentsToBam's output, taken from the reference.
 *
 * A BAM whose header gains @CO lines and whose records are copied block for block. The whole tool
 * is fifteen lines, and every one of them is a decision this dump measures.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE COMMENTS ARE APPENDED IN ORDER AND AFTER ANY THE FILE ALREADY HAD, `addComment` pushing
 *     onto the header's own list;
 *   - `@CO` IS ADDED BY htsjdk AND NOT BY THE TOOL, AND IT IS NOT DOUBLED: `addComment` prefixes
 *     only a comment that does not already start with `@CO\t`, so `@CO\talready prefixed` comes
 *     out exactly once. Measured: the guess was that it doubled;
 *   - A COMMENT HOLDING A TAB SURVIVES WHOLE, so a comment can forge extra header fields;
 *   - A COMMENT HOLDING A NEWLINE NEVER REACHES THE TOOL'S OWN CHECK: the parser refuses it first,
 *     as an IllegalArgumentException naming the character, so the PicardException the tool carries
 *     for that case is unreachable from the command line;
 *   - THE SAM REFUSAL IS ON THE PATH'S SUFFIX, not on the file's contents, so a BAM named `.sam` is
 *     refused by the tool, while a sam named `.bam` gets past that check and is refused LATER, by
 *     the block copy, for having no valid GZIP block at its end;
 *   - THE RECORDS ARE COPIED AS COMPRESSED BLOCKS, so the output's record bytes are the input's
 *     bytes and only the header block is rewritten;
 *   - THE HEADER IS REWRITTEN WHOLE, so its sequence dictionary, read groups and program records
 *     come back in htsjdk's own order rather than the input's;
 *   - NO COMMENT AT ALL IS STILL A REWRITE, which is not a copy: the header block is re-encoded;
 *   - AND CREATE_MD5_FILE AND CREATE_INDEX WRITE BESIDE THE OUTPUT without changing its bytes.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<label>\t<the input BAM, base64>
 *     output\t<label>\t<the rewritten BAM, base64>
 *     sam\t<label>=<the rewritten BAM as text, escaped>
 *     md5\t<label>\t<the .md5's contents>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AddCommentsToBamDump
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
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import picard.sam.AddCommentsToBam;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AddCommentsToBamDump {

    public static void main(final String[] args) throws Exception {
        // Static, and it decides every output byte.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("add-comments-to-bam-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AddCommentsToBamDump: a BAM header gaining @CO lines");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // A BAM whose header already carries a comment, a read group and a program record.
        final Path bam = dir.resolve("reads.bam");
        buildBam(bam);
        System.out.printf("fixture\treads\t%s%n", RecordTransformDump.base64(bam));
        // The same bytes under a `.sam` name, which the tool refuses on the suffix alone.
        final Path named = dir.resolve("misnamed.sam");
        Files.copy(bam, named);
        // And a real sam file under a `.bam` name, which gets past the suffix check.
        final Path samAsBam = dir.resolve("really-sam.bam");
        buildSam(samAsBam);
        System.out.printf("fixture\treally-sam\t%s%n", RecordTransformDump.base64(samAsBam));

        run(dir, "one-comment", bam, List.of("a comment"));
        run(dir, "two-comments", bam, List.of("first", "second"));
        // No comment at all, which is still a rewrite of the header block.
        run(dir, "no-comment", bam, List.of());
        // A comment that already carries the prefix, which is NOT added again.
        run(dir, "already-prefixed", bam, List.of("@CO\talready prefixed"));
        // A comment holding a tab, which reaches the file whole.
        run(dir, "with-tab", bam, List.of("key\tvalue"));
        // A comment holding a newline, which the parser refuses before the tool sees it.
        run(dir, "with-newline", bam, List.of("first line\nsecond line"));
        // The digest and the index, neither of which changes the output's bytes.
        run(dir, "md5", bam, List.of("a comment"), "CREATE_MD5_FILE=true");
        run(dir, "indexed", bam, List.of("a comment"), "CREATE_INDEX=true");
        // The suffix refusal, and the file that is a sam but is not named one.
        run(dir, "named-sam", named, List.of("a comment"));
        run(dir, "sam-named-bam", samAsBam, List.of("a comment"));
    }

    static void run(final Path dir, final String label, final Path input,
                    final List<String> comments, final String... extra) throws Exception {
        final Path out = dir.resolve("commented-" + label + ".bam");
        final List<String> argv = new ArrayList<>(Arrays.asList("I=" + input, "O=" + out));
        for (final String comment : comments) {
            argv.add("C=" + comment);
        }
        argv.add("USE_JDK_DEFLATER=true");
        argv.add("USE_JDK_INFLATER=true");
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new AddCommentsToBam().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s\t%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(out));
        System.out.printf("sam\t%s=%s%n", label, ReferenceQueryDump.escape(asText(out)));
        final Path digest = dir.resolve("commented-" + label + ".bam.md5");
        if (Files.exists(digest)) {
            System.out.printf("md5\t%s\t%s%n", label, Files.readString(digest));
        }
    }

    static String asText(final Path bam) {
        final StringBuilder text = new StringBuilder();
        try (SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(new File(bam.toString()))) {
            text.append(reader.getFileHeader().getSAMString());
            for (final SAMRecord record : reader) {
                text.append(record.getSAMString());
            }
        } catch (final Exception e) {
            text.append("error: ").append(e);
        }
        return text.toString();
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000),
                new SAMSequenceRecord("chr2", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        group.setLibrary("lib1");
        header.addReadGroup(group);
        final SAMProgramRecord program = new SAMProgramRecord("upstream");
        program.setProgramVersion("1.0");
        header.addProgramRecord(program);
        // A comment the file already had, which the tool's own comments follow.
        header.addComment("an existing comment");
        return header;
    }

    static void buildBam(final Path file) {
        final SAMFileHeader header = header();
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, file.toFile())) {
            for (final SAMRecord record : records(header)) {
                writer.addAlignment(record);
            }
        }
    }

    static void buildSam(final Path file) {
        final SAMFileHeader header = header();
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeSAMWriter(header, true, file.toFile())) {
            for (final SAMRecord record : records(header)) {
                writer.addAlignment(record);
            }
        }
    }

    static List<SAMRecord> records(final SAMFileHeader header) {
        final List<SAMRecord> records = new ArrayList<>();
        for (final int start : new int[] {100, 200, 300}) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName("r" + start);
            record.setReferenceName("chr1");
            record.setAlignmentStart(start);
            record.setCigarString("10M");
            final byte[] bases = new byte[10];
            Arrays.fill(bases, (byte) 'A');
            record.setReadBases(bases);
            final byte[] qualities = new byte[10];
            Arrays.fill(qualities, (byte) 30);
            record.setBaseQualities(qualities);
            record.setMappingQuality(60);
            record.setAttribute("RG", "rg1");
            records.add(record);
        }
        return records;
    }

    /** The dump's own directory, whose absolute path reaches the refusals. */
    static String masked(final String text, final Path dir) {
        return text.replace(dir.toAbsolutePath().toString(), "<dir>");
    }
}
