/*
 * RampedHaplotypeCaller's off ramps end to end, taken from the reference.
 *
 * HcEndToEndDump's reads and reference, called by RampedHaplotypeCaller the way a user runs it,
 * with each off ramp in turn. An off ramp halts every region at one step, writes what the caller
 * holds there into a zip, and calls nothing, so the VCF keeps only its header.
 *
 * Five behaviours this is built to catch.
 *
 *   - PRE_ASSEMBLER_OFF WRITES ONLY THE REGIONS, before anything is assembled: each region's
 *     contig, start and end in info.json, and no other entry;
 *   - POST_ASSEMBLER_OFF WRITES THE ASSEMBLED HAPLOTYPES as a CSV and as a BAM of artificial
 *     haplotype reads, with the trimmer's flanks and the padded and genotyping spans;
 *   - THE HAPLOTYPE BAM CARRIES THE SOURCE HEADER plus the ArtificialHaplotypeRG read group and
 *     a HaplotypeBAMWriter program record, is coordinate sorted, and is indexed as it is written;
 *   - PRE_FILTER_OFF WRITES THE READ LIKELIHOODS against every haplotype, per sample, before the
 *     filter that drops poorly modelled reads;
 *   - AND info.json LISTS ITS KEYS IN HashMap ORDER, a lone member inline, and is written last.
 *
 * The post-assembler run writes its BAM through the JDK deflater, so the bytes do not depend on
 * whether the GKL deflater is available on the machine that replays them. It runs FIRST: a run
 * without --use-jdk-deflater installs GKL as BlockCompressedOutputStream's static default and
 * nothing resets it, so in one JVM every later run would write through GKL whatever it asks for.
 *
 * Runs, in this order and in ONE JVM.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output VCF from its #CHROM line on, escaped>
 *     entry\t<label>=<the zip's entry names, comma separated, in zip order>
 *     content\t<label>/<entry>=<a text entry, escaped>
 *     bytes\t<label>/<entry>=<a binary entry (the BAM and its index), as lower-case hex>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: RampsEndToEndDump
 */

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
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class RampsEndToEndDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("ramps-end-to-end-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# RampsEndToEndDump: RampedHaplotypeCaller's off ramps end to end");

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

        run(dir, "post-assembler-off", bam, fasta, "POST_ASSEMBLER_OFF",
                List.of("--use-jdk-deflater", "true"));
        run(dir, "pre-assembler-off", bam, fasta, "PRE_ASSEMBLER_OFF",
                List.of("-L", "chr1:200-600", "-L", "chr2"));
        run(dir, "pre-filter-off", bam, fasta, "PRE_FILTER_OFF",
                List.of("-L", "chr1:1150-1350"));
        run(dir, "no-ramp-file", bam, fasta, "PRE_FILTER_OFF", null);
    }

    /** One run with an off ramp of `type`; a null `extra` leaves the ramp file out. */
    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final String type, final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".vcf");
        final Path zip = dir.resolve(label + ".zip");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--off-ramp-type", type));
        if (extra != null) {
            argv.addAll(List.of("--off-ramp-file", zip.toString()));
            argv.addAll(extra);
        }
        try {
            new RampedHaplotypeCaller().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(
                            HcEndToEndDump.masked(String.valueOf(cause.getMessage()), dir)));
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

        try (ZipFile file = new ZipFile(zip.toFile())) {
            final List<String> names = new ArrayList<>();
            for (final ZipEntry entry : Collections.list(file.entries())) {
                names.add(entry.getName());
            }
            System.out.printf("entry\t%s=%s%n", label, String.join(",", names));
            for (final ZipEntry entry : Collections.list(file.entries())) {
                final byte[] bytes = file.getInputStream(entry).readAllBytes();
                if (entry.getName().endsWith(".bam") || entry.getName().endsWith(".bai")) {
                    final StringBuilder hex = new StringBuilder();
                    for (final byte b : bytes) {
                        hex.append(String.format("%02x", b & 0xff));
                    }
                    System.out.printf("bytes\t%s/%s=%s%n", label, entry.getName(), hex);
                } else {
                    System.out.printf("content\t%s/%s=%s%n", label, entry.getName(),
                            ReferenceQueryDump.escape(
                                    new String(bytes, StandardCharsets.UTF_8)));
                }
            }
        }
    }
}
