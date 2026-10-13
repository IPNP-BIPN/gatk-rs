/*
 * Mutect2 end to end, taken from the reference.
 *
 * HcEndToEndDump's two samples and two contigs, read as a tumor (s1) and, when -normal names it,
 * a normal (s2), called by the tool the way a user runs it. Everything the somatic caller adds
 * over HaplotypeCaller is in one output:
 *
 *   - THE ACTIVE REGIONS COME FROM THE TUMOR'S LOG ODDS, not from a genotyping engine: per locus,
 *     the likeliest alternate (a base, another substitution, or an indel) and the odds of a
 *     somatic allele fraction against sequencing error;
 *   - THE CALLS ARE SOMATIC: no diploid genotype, a tumor allele fraction (AF) from the somatic
 *     likelihoods model, TLOD, and with a normal NLOD and NALOD;
 *   - A NORMAL SILENCES A SITE where it carries the tumor's alternate in more than a third of its
 *     reads with enough quality;
 *   - THE DEFAULT ANNOTATIONS ARE MUTECT2'S, with the allele-specific strand table, the median
 *     base and mapping qualities and fragment lengths, the event count and the read position;
 *   - AND THE STATS FILE carries the callable sites FilterMutectCalls needs.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output VCF from its #CHROM line on, escaped>
 *     stats\t<label>=<the .stats file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: Mutect2EndToEndDump
 */

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

public class Mutect2EndToEndDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-end-to-end-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2EndToEndDump: Mutect2 end to end");

        final Random random = new Random(20261010L);
        final String[] references = new String[HcEndToEndDump.CONTIGS.length];
        for (int c = 0; c < HcEndToEndDump.CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < HcEndToEndDump.LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\t%s=%s%n", HcEndToEndDump.CONTIGS[c], references[c]);
        }
        final Path fasta = HcEndToEndDump.writeReference(dir, references);

        final String sam = HcEndToEndDump.sam(references, random);
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

        run(dir, "tumor-only", bam, fasta, List.of());
        run(dir, "tumor-normal", bam, fasta, List.of("-normal", "s2"));
        run(dir, "tumor-normal-interval", bam, fasta,
                List.of("-normal", "s2", "-L", "chr1:200-800", "-L", "chr2"));
        run(dir, "genotype-germline", bam, fasta,
                List.of("-normal", "s2", "--genotype-germline-sites", "true"));
        run(dir, "initial-tumor-lod-10", bam, fasta,
                List.of("-normal", "s2", "--initial-tumor-lod", "10", "--tumor-lod-to-emit", "10"));
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new Mutect2().instanceMain(argv.toArray(new String[0]));
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
            if (inBody && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
        final Path stats = dir.resolve("out-" + label + ".vcf.stats");
        System.out.printf("stats\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.exists(stats) ? Files.readString(stats) : "<none>"));
    }
}
