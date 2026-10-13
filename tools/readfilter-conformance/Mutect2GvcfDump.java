/*
 * Mutect2 in somatic reference-confidence mode, end to end (G3.2d, #1333).
 *
 * With --emit-ref-confidence the emission and activity thresholds drop to minus infinity, so every
 * region is active and every alternate is emitted, <NON_REF> among them. A base with no call is a
 * record whose likelihood is the somatic log odds of its alternate reads (no indel model, so no GQ
 * or PL, and a TLOD in FORMAT), and SomaticGVCFWriter merges them into blocks by TLOD band rather
 * than by GQ: the bands are decimal LODs, binned to the precision their smallest gap needs, and a
 * block carries the smallest TLOD it saw.
 *
 * Over the reads of HcDragenEndToEndDump's first sample alone, since the mode takes one sample:
 * GVCF over the contig, GVCF and BP_RESOLUTION over a window, GVCF with three bands a quarter apart
 * (two decimal places), GVCF with --enable-all-annotations, and the two-sample BAM, refused.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text of the one-sample BAM, escaped>
 *     sam\tboth=<the SAM text of the two-sample BAM, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: Mutect2GvcfDump
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

public class Mutect2GvcfDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-gvcf-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2GvcfDump: Mutect2 in somatic reference-confidence mode");

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

        run(dir, "gvcf", Mutect2::new, bam, fasta, List.of("-ERC", "GVCF"));
        run(dir, "gvcf-window", Mutect2::new, bam, fasta, List.of("-ERC", "GVCF", "-L", "chr1:280-420"));
        run(dir, "bp-resolution", Mutect2::new, bam, fasta, List.of("-ERC", "BP_RESOLUTION", "-L", "chr1:290-310"));
        run(dir, "gvcf-quarter-bands", Mutect2::new, bam, fasta,
                List.of("-ERC", "GVCF", "-L", "chr1:1-600", "-LODB", "-1.5", "-LODB", "-1.25", "-LODB", "-1.0"));
        run(dir, "gvcf-all-annotations", Mutect2::new, bam, fasta,
                List.of("-ERC", "GVCF", "-L", "chr1:280-420", "--enable-all-annotations"));
        run(dir, "two-samples", Mutect2::new, twoSamples, fasta, List.of("-ERC", "GVCF", "-L", "chr1:280-320"));
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
