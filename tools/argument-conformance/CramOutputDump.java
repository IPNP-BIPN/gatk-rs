/*
 * What a read walker writes for an output named `.cram`, `.sam` or anything else: the bytes of
 * every file it leaves, or the refusal it ends with.
 *
 * `GATKTool.createSAMWriter` hands the output's NAME to htsjdk's `SAMFileWriterFactory.makeWriter`,
 * which picks the format from the extension: a `CRAMFileWriter` for `.cram` (htsjdk 5.0.0's CRAM
 * 3.1 under the NORMAL profile, indexed by `CRAMBAIIndexer` into `<output>.bai`), a
 * `SAMTextWriter` for `.sam` (never indexed), and a `BAMFileWriter` for everything else (indexed
 * by `<output without .bam>.bai`). Before it, a `.cram` output without a reference is refused, and
 * `--output-cram-version` is read for every output, a BAM's included.
 *
 * Every path is relative to the working directory, because the command line is recorded in the
 * output's `@PG` record. A CRAM's file definition also carries the output's absolute URI, cut to
 * twenty bytes, which no relative path fixes; the port's test puts the golden's twenty bytes back
 * before it compares.
 *
 * Output:
 *
 *     input\t<relative path>\t<base64>
 *     file\t<case>\t<file name>\t<base64>
 *     error\t<case>\t<class>:<message>
 *
 * Usage: CramOutputDump (in a directory holding the corpus's `fixtures/`)
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class CramOutputDump {

    static final Path FIXTURES = Paths.get("fixtures");
    static final Path OUT = Paths.get("cram-output");
    static final String[] INPUTS = {
        "reads.bam", "reads.bai", "reads.bam.bai",
        "reference.fasta", "reference.fasta.fai", "reference.dict",
    };

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramOutputDump: what a read walker writes, chosen by the output's name");
        for (final String name : INPUTS) {
            final Path path = FIXTURES.resolve(name);
            if (Files.exists(path)) {
                System.out.printf("input\t%s\t%s%n", path, base64(path));
            }
        }

        final String reads = FIXTURES.resolve("reads.bam").toString();
        final String reference = FIXTURES.resolve("reference.fasta").toString();
        run("cram", "out.cram", "-I", reads, "-R", reference);
        run("cram-3.0", "out.cram", "-I", reads, "-R", reference, "--output-cram-version", "3.0");
        run("cram-3.1", "out.cram", "-I", reads, "-R", reference, "--output-cram-version", "3.1");
        run("cram-md5", "out.cram", "-I", reads, "-R", reference, "--create-output-bam-md5");
        run("cram-no-index", "out.cram", "-I", reads, "-R", reference,
                "--create-output-bam-index", "false");
        run("cram-interval", "out.cram", "-I", reads, "-R", reference, "-L", "chr1:100-400");
        run("cram-no-reference", "out.cram", "-I", reads);
        run("sam", "out.sam", "-I", reads);
        run("sam-md5", "out.sam", "-I", reads, "--create-output-bam-md5");
        run("other-extension", "out.tmp", "-I", reads);
        run("bad-cram-version", "out.bam", "-I", reads, "--output-cram-version", "3.2");

        // Last, because writing a BAM outside GATK loads htsjdk's `Defaults` before `Main` has put
        // GATK's compression level 2 in the system properties, and every later BAM would be level 5.
        // Beside the fixtures rather than among them, which the container may mount read-only.
        final Path queryname = Paths.get("inputs", "queryname.bam");
        Files.createDirectories(queryname.getParent());
        writeQueryname(FIXTURES.resolve("reads.bam"), queryname);
        System.out.printf("input\t%s\t%s%n", queryname, base64(queryname));
        run("cram-queryname", "out.cram", "-I", queryname.toString(), "-R", reference);
    }

    static void run(final String name, final String output, final String... extra) throws Exception {
        final Path dir = OUT.resolve(name);
        Files.createDirectories(dir);
        final List<String> argv = new ArrayList<>(List.of("PrintReads"));
        argv.addAll(Arrays.asList(extra));
        argv.addAll(List.of("-O", dir.resolve(output).toString()));

        // The tool's own logging carries a clock; what is measured is the files and the refusal.
        final PrintStream original = System.out;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            new org.broadinstitute.hellbender.Main().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception e) {
            System.setOut(original);
            emit("error", name, e.getClass().getName() + ":" + e.getMessage());
        } finally {
            System.setOut(original);
        }
        try (final var entries = Files.list(dir)) {
            for (final Path file : entries.sorted().toList()) {
                System.out.printf("file\t%s\t%s\t%s%n", name, file.getFileName(), base64(file));
            }
        }
    }

    /** The corpus's reads sorted by name, a header that is not coordinate-sorted over the same reference. */
    static void writeQueryname(final Path input, final Path output) throws Exception {
        try (final htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault()
                .open(input.toFile())) {
            final htsjdk.samtools.SAMFileHeader header = reader.getFileHeader().clone();
            header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.queryname);
            try (final htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory()
                    .makeBAMWriter(header, false, output.toFile())) {
                for (final htsjdk.samtools.SAMRecord record : reader) {
                    writer.addAlignment(record);
                }
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }

    static void emit(final String kind, final String name, final String payload) {
        System.out.printf("%s\t%s\t%s%n", kind, name,
                payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }
}
