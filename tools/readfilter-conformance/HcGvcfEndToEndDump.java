/*
 * HaplotypeCaller in reference-confidence mode, end to end, taken from the reference.
 *
 * The reads of HcEndToEndDump's first sample alone (GVCF mode refuses more than one), over the
 * same two contigs, called by the tool with --emit-ref-confidence. What GVCF mode changes, all of
 * it in one output:
 *
 *   - EVERY BASE OF THE TRAVERSAL IS ACCOUNTED FOR: an inactive region, an active one that
 *     assembles to the reference, and the flanks the trimmer cut off a variant region are each
 *     given a reference-confidence record per base, from the ref-vs-any likelihoods of its pileup
 *     and the reads with no plausible indel, the less confident of the two;
 *   - A <NON_REF> ALLELE IS ADDED TO EVERY CALL, its likelihood the median of each read's
 *     non-best concrete alleles, and the calling threshold is zero, so a site is written whatever
 *     its quality;
 *   - THE ANNOTATIONS CHANGE: StrandBiasBySample is added, ChromosomeCounts, FisherStrand,
 *     StrandOddsRatio and QualByDepth are removed, and the reducible ones are written raw;
 *   - PHYSICAL PHASING IS ON, so the two SNPs in cis carry PGT, PID and PS;
 *   - AND GVCF COMBINES THE REFERENCE RECORDS INTO BLOCKS by GQ band, where BP_RESOLUTION writes
 *     every one of them.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output from its #CHROM line on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcGvcfEndToEndDump
 */

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
import java.util.stream.Collectors;

public class HcGvcfEndToEndDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-gvcf-end-to-end-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcGvcfEndToEndDump: HaplotypeCaller in reference-confidence mode");

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

        // The first sample's reads only.
        final String sam = HcEndToEndDump.sam(references, random).lines()
                .filter(line -> !line.contains("rg2"))
                .map(line -> line + "\n")
                .collect(Collectors.joining());
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

        HcEndToEndDump.run(dir, "gvcf", bam, fasta, List.of("-ERC", "GVCF"));
        HcEndToEndDump.run(dir, "bp-resolution", bam, fasta,
                List.of("-ERC", "BP_RESOLUTION", "-L", "chr1:1150-1350"));
        HcEndToEndDump.run(dir, "gvcf-interval", bam, fasta,
                List.of("-ERC", "GVCF", "-L", "chr1:200-600", "-L", "chr2"));
        HcEndToEndDump.run(dir, "gvcf-no-phasing", bam, fasta,
                List.of("-ERC", "GVCF", "--do-not-run-physical-phasing", "true"));
        HcEndToEndDump.run(dir, "gvcf-bands", bam, fasta,
                List.of("-ERC", "GVCF", "-GQB", "20", "-GQB", "60"));
    }
}
