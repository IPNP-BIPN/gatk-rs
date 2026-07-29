/*
 * PrintReads' output bytes, taken from the reference.
 *
 * This is the first whole tool measured here, and the claim it carries is different from every
 * suite before it: not "the right reads in the right order" but "the right bytes". The output BAM
 * travels in the golden in full, base64, index included, so the port is compared against the file
 * the reference wrote rather than against a description of it.
 *
 * Two things the tool does to the header are worth naming, because a port that skipped either
 * would produce a valid BAM that is not this BAM:
 *
 *   - a @PG record is appended, with ID = the tool name, VN = the GATK version, CL = the whole
 *     command line and PN = the tool name. Its ID collides deliberately: a second PrintReads over
 *     an output of the first gets ID `PrintReads.1`, because createProgramGroupID appends
 *     consecutive integers until the ID is free;
 *   - getHeaderForSAMWriter mutates the *reads* header in place rather than copying it.
 *
 * The command line lands in the golden as its own row, because it is an input to the writer and
 * not something a port can invent: it carries the temporary paths of the run that produced it.
 *
 * Which deflater wrote these bytes is not a detail, and it is not what the flag suggests.
 * `--use-jdk-deflater true` does **not** restore the JDK deflater: GATK only ever *installs* the
 * Intel one, in `if (!useJdkDeflater) setDefaultDeflaterFactory(new IntelDeflaterFactory())`, and
 * that setter is static and global. Picard's CommandLineProgram has the same shape, so a single
 * earlier CreateSequenceDictionary call in the same JVM leaves the Intel deflater installed for
 * everything that follows, and the flag is then a no-op. The first version of this dump did
 * exactly that: its output was GKL-compressed while claiming to be JDK-compressed, 708 deflate
 * bytes where zlib gives 698 at level 5 and 706 at level 4, matching no zlib setting at all.
 *
 * So the factory is installed explicitly here, and the golden records which one produced it.
 * htsjdk-rs reproduces the JDK deflater exactly; GKL-exact deflate is separate work, and until it
 * exists a byte claim over BGZF has to name the deflater it is a claim about.
 *
 * Output:
 *
 *     deflater\t<the DeflaterFactory class that produced every byte below>
 *     bam\t<base64 input BAM>       bai\t<base64 input index>
 *     commandline\t<label>\t<the CL string the tool recorded>
 *     header\t<label>\t<output SAM header text, \n escaped>
 *     output\t<label>\t<base64 of the whole output BAM>
 *     index\t<label>\t<base64 of the output .bai>
 *
 * Usage: PrintReadsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.PrintReads;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class PrintReadsDump {

    public static void main(final String[] args) throws Exception {
        // Before anything else, and before the fixture is written: the factory is static, and
        // whoever touches it first wins for the life of the JVM.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Files.createTempDirectory("printreads");
        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# PrintReadsDump: the bytes PrintReads writes");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        System.out.printf("bam\t%s%n", base64(bam));
        System.out.printf("bai\t%s%n", base64(dir.resolve("reads.bai")));

        run(dir, bam, "all", new String[] {});
        run(dir, bam, "chr1", new String[] {"-L", "chr1"});
        run(dir, bam, "chr1:100-160", new String[] {"-L", "chr1:100-160"});
        run(dir, bam, "nofilter",
                new String[] {"--disable-tool-default-read-filters", "true"});
        run(dir, bam, "nodup", new String[] {"--read-filter", "NotDuplicateReadFilter"});
        // No index requested: the output is one file rather than two, and the BAM's own bytes are
        // unchanged by the absence of its index.
        run(dir, bam, "noindex", new String[] {"--create-output-bam-index", "false"});
    }

    static void run(final Path dir, final Path bam, final String label, final String[] extra)
            throws Exception {
        final Path output = dir.resolve("out." + label.replace(':', '_') + ".bam");
        // --use-jdk-deflater is the knob that decides which bytes come out. GATK's default is the
        // Intel GKL deflater, whose output htsjdk-rs does not yet reproduce; the JDK deflater's
        // it does. Naming it here is what makes the byte claim a claim about something.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        new PrintReads().instanceMain(argv.toArray(new String[0]));

        // The @PG the tool appended, read back from what it wrote.
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
        System.out.printf("output\t%s\t%s%n", label, base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? base64(index) : "absent");
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }
}
