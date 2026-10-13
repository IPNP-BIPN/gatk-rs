/*
 * The callers with the linked de Bruijn graph, end to end (G3.2h, #1337).
 *
 * --linked-de-bruijn-graph swaps the assembler's sequence graph for JunctionTreeLinkedDeBruijnGraph:
 * every k-mer size in turn (the configured ones, then six larger) is built as a read threading
 * graph with no non-unique k-mers and cycles allowed, its dangling ends recovered along the stored
 * reference path, and each read threaded through it again to record, at every vertex leading into
 * a merge, which way it went at each later fork. JunctionTreeKBestHaplotypeFinder then walks from
 * the reference source, following those junction trees where they have three threads of evidence
 * and the edge weights where not, and ends a path at the reference sink once no open tree vouches
 * for going on; with artificial haplotype recovery (the default), a fork edge no path took is
 * grafted onto the best path through its vertex. The first k-mer size that yields haplotypes is
 * kept, in the order found, after the reference. --disable-artificial-haplotype-recovery is
 * refused in this mode, the only one that would read it: the reference's test is the wrong way
 * round, and the refusal names an argument that does not exist.
 *
 * Paths of exactly equal score leave the finder's queue in an order the reference takes from
 * identity hashes; the calls do not depend on it here.
 *
 * Two read sets: HcEndToEndDump's (two samples, two SNPs thirty bases apart in cis for one and in
 * trans for the other, indels, soft clips) and HcDragenEndToEndDump's (three short tandem repeats,
 * which give the graph cycles). Runs in ONE JVM, in this order, as HcEndToEndDump's are.
 *
 * Output:
 *
 *     fasta\t<set>:<contig>=<the contig's bases>
 *     sam\t<set>=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcLinkedGraphDump
 */

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
import java.util.List;
import java.util.Random;

public class HcLinkedGraphDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-linked-graph-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcLinkedGraphDump: the callers with the linked de Bruijn graph");

        // HcEndToEndDump's contigs and reads, from the same seed.
        final Random random = new Random(20261010L);
        final String[] references = new String[HcEndToEndDump.CONTIGS.length];
        for (int c = 0; c < HcEndToEndDump.CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < HcEndToEndDump.LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\te2e:%s=%s%n", HcEndToEndDump.CONTIGS[c], references[c]);
        }
        final Path e2eDir = dir.resolve("e2e");
        Files.createDirectories(e2eDir);
        final Path e2eFasta = HcEndToEndDump.writeReference(e2eDir, references);
        final Path e2eBam = bam(e2eDir, "e2e", HcEndToEndDump.sam(references, random));

        // HcDragenEndToEndDump's contig and reads.
        final Random strRandom = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(strRandom.nextInt(4)));
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
        final String strReference = ref.toString();
        System.out.printf("fasta\tstr:%s=%s%n", HcDragenEndToEndDump.CONTIG, strReference);
        final Path strDir = dir.resolve("str");
        Files.createDirectories(strDir);
        final Path strFasta = HcDragenEndToEndDump.writeReference(strDir, strReference);
        final Path strBam = bam(strDir, "str", HcDragenEndToEndDump.sam(strReference));

        final String linked = "--linked-de-bruijn-graph";
        HcAllAnnotationsDump.run(dir, "hc-linked", HaplotypeCaller::new, e2eBam, e2eFasta, List.of(linked));
        HcAllAnnotationsDump.run(dir, "hc-linked-no-recovery", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(linked, "--disable-artificial-haplotype-recovery"));
        HcAllAnnotationsDump.run(dir, "hc-recovery-without-linked", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of("--disable-artificial-haplotype-recovery"));
        HcAllAnnotationsDump.run(dir, "hc-linked-kmer-15", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(linked, "--kmer-size", "15", "--dont-increase-kmer-sizes-for-cycles"));
        HcAllAnnotationsDump.run(dir, "hc-linked-recover-all", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(linked, "--recover-all-dangling-branches"));
        HcAllAnnotationsDump.run(dir, "hc-linked-gvcf", HaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(linked, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"));
        HcAllAnnotationsDump.run(dir, "ramped-linked", RampedHaplotypeCaller::new, e2eBam, e2eFasta,
                List.of(linked));
        HcAllAnnotationsDump.run(dir, "mutect2-linked", Mutect2::new, e2eBam, e2eFasta, List.of(linked));
        HcAllAnnotationsDump.run(dir, "mutect2-linked-tumor-normal", Mutect2::new, e2eBam, e2eFasta,
                List.of(linked, "--normal-sample", "s2"));
        HcAllAnnotationsDump.run(dir, "str-hc-linked", HaplotypeCaller::new, strBam, strFasta, List.of(linked));
        HcAllAnnotationsDump.run(dir, "str-hc-linked-no-recovery", HaplotypeCaller::new, strBam, strFasta,
                List.of(linked, "--disable-artificial-haplotype-recovery"));
        HcAllAnnotationsDump.run(dir, "str-mutect2-linked", Mutect2::new, strBam, strFasta, List.of(linked));
    }

    static Path bam(final Path dir, final String set, final String sam) throws Exception {
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
        System.out.printf("sam\t%s=%s%n", set, ReferenceQueryDump.escape(sam));
        return bam;
    }
}
