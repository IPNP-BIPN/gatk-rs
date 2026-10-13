/*
 * HaplotypeCaller and RampedHaplotypeCaller with --disable-optimizations, end to end (G4.2, #1331).
 *
 * The flag lifts three shortcuts of callRegion: a region the trimmer finds no variation in, an
 * assembly with no variation, and a genotyping region with no reads left would each answer the
 * reference model at once. Lifted, the first is not genotyped either: it asks the trimmer's result
 * for a variant region it does not have, and the run fails with "There is no variation present."
 * In VCF mode no region reaches that point on these reads, since only regions with calls are
 * active enough to be kept; in GVCF mode every region is, and the run fails.
 *
 * One more run measures what is not this flag's but surfaced beside it in the arrays: the
 * writer HaplotypeCaller builds ignores --lenient.
 *
 * The reads and reference are HcDragenEndToEndDump's: one contig with three short tandem repeats
 * and two samples.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output VCF's ##INFO lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcDisableOptimizationsDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;

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

public class HcDisableOptimizationsDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-disable-optimizations-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcDisableOptimizationsDump: HaplotypeCaller with --disable-optimizations");

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

        run(dir, "vcf", false, bam, fasta, List.of("--disable-optimizations"));
        run(dir, "vcf-dragen-mode", false, bam, fasta, List.of("--disable-optimizations", "--dragen-mode"));
        run(dir, "gvcf", false, bam, fasta, List.of("--disable-optimizations", "-ERC", "GVCF", "--sample-name", "s1"));
        run(dir, "no-active-region", false, bam, fasta, List.of("--disable-optimizations", "-L", "chr1:100-250"));
        // HaplotypeCaller's writer ignores --lenient: the reference-model records carry AD, which
        // the header lacks once the default annotations are off, and the run is refused anyway.
        run(dir, "lenient-ignored", false, bam, fasta, List.of("-ERC", "BP_RESOLUTION", "--sample-name", "s1",
                "--disable-tool-default-annotations", "--lenient", "-L", "chr1:290-310"));
        run(dir, "ramped-vcf", true, bam, fasta, List.of("--disable-optimizations"));
        run(dir, "ramped-gvcf", true, bam, fasta, List.of("--disable-optimizations", "-ERC", "GVCF", "--sample-name", "s1"));
    }

    static void run(final Path dir, final String label, final boolean ramped, final Path bam, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            if (ramped) {
                new RampedHaplotypeCaller().instanceMain(argv.toArray(new String[0]));
            } else {
                new HaplotypeCaller().instanceMain(argv.toArray(new String[0]));
            }
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
            if ((inBody || line.startsWith("##INFO")) && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }
}
