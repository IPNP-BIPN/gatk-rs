/*
 * Mutect2 with --mitochondria-mode, end to end (G3.2e, #1334).
 *
 * The mode is argument presets, set where the command line did not set them: a default allele
 * frequency of 4e-3, emission and initial LODs of 0, every dangling branch recovered, a pruning
 * LOD threshold of log(10^-4), and the OriginalAlignment annotation. Over HcDragenEndToEndDump's
 * reads: the mode on the first sample, the mode with an explicit --tumor-lod-to-emit it must not
 * override, the mode in GVCF, and the mode over both samples with the second as the normal.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text of the one-sample BAM, escaped>
 *     sam\tboth=<the SAM text of the two-sample BAM, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: Mutect2MitochondriaDump
 */

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;

public class Mutect2MitochondriaDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-mitochondria-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2MitochondriaDump: Mutect2 in mitochondria mode");

        // HcDragenEndToEndDump's contig, built the same way from the same seed.
        final Random random = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        HcDragenEndToEndDump.splice(ref, 401, "AAAAAAAAAA");
        HcDragenEndToEndDump.splice(ref, 701, "CACACACACACACA");
        HcDragenEndToEndDump.splice(ref, 1001, "TTTTTTTT");
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String reference = ref.toString();
        System.out.printf("fasta\t%s=%s%n", HcDragenEndToEndDump.CONTIG, reference);
        final Path fasta = HcDragenEndToEndDump.writeReference(dir, reference);

        final String both = HcDragenEndToEndDump.sam(reference);
        final String sam = java.util.Arrays.stream(both.split("\n"))
                .filter(line -> !line.contains("rg2"))
                .collect(java.util.stream.Collectors.joining("\n")) + "\n";
        final Path twoSamples = dir.resolve("both.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(
                     writeSam(dir, "both.sam", both));
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, twoSamples)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        final Path samPath = dir.resolve("input.sam");
        Files.writeString(samPath, sam, StandardCharsets.UTF_8);
        final Path bam = dir.resolve("input.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(samPath);
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam));
        System.out.printf("sam\tboth=%s%n", ReferenceQueryDump.escape(both));

        run(dir, "mitochondria", Mutect2::new, bam, fasta, List.of("--mitochondria-mode"));
        run(dir, "mitochondria-explicit-emission", Mutect2::new, bam, fasta,
                List.of("--mitochondria-mode", "--tumor-lod-to-emit", "3.0"));
        run(dir, "mitochondria-gvcf", Mutect2::new, bam, fasta,
                List.of("--mitochondria-mode", "-ERC", "GVCF", "-L", "chr1:540-600"));
        run(dir, "mitochondria-tumor-normal", Mutect2::new, twoSamples, fasta,
                List.of("--mitochondria-mode", "--normal-sample", "s2"));
    }

    static Path writeSam(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Supplier<CommandLineProgram> tool, final Path bam,
                    final Path fasta, final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            tool.get().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(HcEndToEndDump.masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        final StringBuilder body = new StringBuilder();
        boolean inBody = false;
        for (final String line : Files.readAllLines(out)) {
            inBody |= line.startsWith("#CHROM");
            if ((inBody || (line.startsWith("##") && !line.startsWith("##GATKCommandLine"))) && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }
}
