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
 * The Intel deflater is disabled for this dump (-Dsamjdk.try_use_intel_deflater=false), which the
 * suite declares. htsjdk-rs reproduces the JDK deflater exactly; GKL-exact deflate is a separate
 * piece of work, and until it exists a byte claim over BGZF output has to name which deflater it
 * is a claim about.
 *
 * Output:
 *
 *     bam\t<base64 input BAM>       fai\t<escaped>
 *     bai\t<base64 input index>     fasta\t<escaped>
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
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.PrintReads;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class PrintReadsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("printreads");
        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        System.out.println("# PrintReadsDump: the bytes PrintReads writes");
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(dir.resolve("ref.fasta.fai")))));
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
