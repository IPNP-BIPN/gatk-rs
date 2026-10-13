/*
 * HaplotypeCaller, RampedHaplotypeCaller and Mutect2 with --enable-all-annotations, end to end
 * (G3.2c, #1332).
 *
 * Every annotation GATK's plugin discovers, computed in the context each caller gives it: the
 * allele-specific raw data reference-confidence mode writes, AllelePseudoDepth's DD and DF, the
 * flow annotations under the default flow order (CycleSkipStatus, which needs an actual one,
 * writing nothing), the pedigree annotations with no pedigree (nothing), OriginalAlignment's OCM
 * where Mutect2's TLOD lets it run, and the jumbo INFO annotations, AssemblyComplexity and
 * HaplotypeFilteringAnnotation, over the caller's read- or fragment-by-haplotype matrix.
 *
 * What a reader would not guess: HaplotypeCaller refuses allele-specific annotations in VCF mode;
 * Mutect2 writes DD, DF and BQHIST in the Allele.hashCode order its HashMap marginalises into; and
 * GQ_MEAN reads the genotype's real GQ, which only the VCF text caps at 99.
 *
 * The reads and reference are HcDragenEndToEndDump's.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcAllAnnotationsDump
 */

import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
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

public class HcAllAnnotationsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-all-annotations-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcAllAnnotationsDump: the callers with --enable-all-annotations");

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

        final String sam = HcDragenEndToEndDump.sam(reference);
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

        run(dir, "hc-vcf", HaplotypeCaller::new, bam, fasta, List.of("--enable-all-annotations"));
        run(dir, "hc-gvcf", HaplotypeCaller::new, bam, fasta,
                List.of("--enable-all-annotations", "-ERC", "GVCF", "--sample-name", "s1"));
        run(dir, "hc-bp-resolution", HaplotypeCaller::new, bam, fasta,
                List.of("--enable-all-annotations", "-ERC", "BP_RESOLUTION", "--sample-name", "s2", "-L", "chr1:550-570"));
        run(dir, "hc-gvcf-dragen-mode", HaplotypeCaller::new, bam, fasta,
                List.of("--enable-all-annotations", "-ERC", "GVCF", "--sample-name", "s2", "--dragen-mode"));
        run(dir, "ramped-gvcf", RampedHaplotypeCaller::new, bam, fasta,
                List.of("--enable-all-annotations", "-ERC", "GVCF", "--sample-name", "s1"));
        run(dir, "mutect2", Mutect2::new, bam, fasta, List.of("--enable-all-annotations"));
        run(dir, "mutect2-tumor-normal", Mutect2::new, bam, fasta,
                List.of("--enable-all-annotations", "--normal-sample", "s2"));
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
                    ReferenceQueryDump.escape(HcEndToEndDump.masked(String.valueOf(cause.getMessage()), dir)
                            // PedReader names its reader by Object.toString, an identity hash.
                            .replaceAll("java\\.io\\.InputStreamReader@[0-9a-f]+", "java.io.InputStreamReader@<hash>")));
            return;
        }
        final StringBuilder body = new StringBuilder();
        boolean inBody = false;
        for (final String line : Files.readAllLines(out)) {
            inBody |= line.startsWith("#CHROM");
            if ((inBody || line.startsWith("##INFO") || line.startsWith("##FORMAT")) && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }
}
