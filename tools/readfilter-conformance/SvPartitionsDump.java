/*
 * FindBreakpointEvidenceSpark and StructuralVariationDiscoveryPipelineSpark through Main on BAMs
 * that Disq cuts into several partitions, and under dictionaries that make KmerCleaner's table large.
 *
 * SvEvidenceDump's fixtures and SvEvidenceDump's rows, with a small `--bam-partition-size`: a local
 * file is cut at that size or at 32 MiB, whichever is less, the last split taking a tenth over.
 * reads.bam is some hundreds of KiB, so these make a handful to some thousands of partitions, most of
 * them empty at the small end.
 *
 * Output and Usage: as SvEvidenceDump, whose rows these are not, so that its golden stands.
 */

import htsjdk.samtools.BAMSBIIndexer;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Base64;
import java.util.List;

public class SvPartitionsDump {

    public static void main(final String[] args) throws Exception {
        PrintReadsDump.emptyDirectory(SvEvidenceDump.DIR);
        Files.createDirectories(SvEvidenceDump.DIR);
        System.out.println("# SvPartitionsDump: the last two Spark tools through Main, on several partitions");
        SvEvidenceDump.fixtures();
        fixtures();
        for (final String name : new String[] {"sv.fasta", "sv.fasta.fai", "sv.dict", "sv.fasta.img", "reads.bam",
                "reads.bai", "shallow.bam", "shallow.bai", "other.bam", "other.bai", "not-sorted.bam", "kill.txt",
                "kill31.txt", "cnv.vcf", "external.bed", "external-long.bed", "gaps.txt", "cross.txt", "sbi.bam",
                "sbi.bam.sbi", "wide.bam", "wide.bai", "wrapped.bam", "wrapped.bai"}) {
            final byte[] bytes = Files.readAllBytes(SvEvidenceDump.DIR.resolve(name));
            System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(
                    name.endsWith(".img") ? SvEvidenceDump.withoutAddresses(bytes) : bytes));
        }
        rows();
        System.out.flush();
        // Spark's threads would keep the container up.
        System.exit(0);
    }

    static String at(final String name) {
        return SvEvidenceDump.at(name);
    }

    static void fixtures() throws Exception {
        // Features across the whole of the contigs, long ones reaching over partitions.
        Files.writeString(SvEvidenceDump.DIR.resolve("external-long.bed"), "chr1\t500\t600\nchr1\t1000\t12000\n"
                + "chr1\t13000\t13100\nchr2\t500\t8000\nchr2\t8500\t8600\n");
        Files.copy(SvEvidenceDump.DIR.resolve("reads.bam"), SvEvidenceDump.DIR.resolve("sbi.bam"));
        BAMSBIIndexer.createIndex(SvEvidenceDump.DIR.resolve("sbi.bam"), 4096L);
        // The reads under a dictionary that makes KmerCleaner's table over 2^27 entries.
        reheadered("wide.bam", 70_000_000);
        // And one whose table size narrows to a negative int.
        reheadered("wrapped.bam", 1_200_000_000);
    }

    static void reheadered(final String name, final int firstContigLength) throws Exception {
        final List<SAMRecord> records = new ArrayList<>();
        final SAMFileHeader header;
        try (SamReader reader = SamReaderFactory.makeDefault().open(SvEvidenceDump.DIR.resolve("reads.bam"))) {
            header = reader.getFileHeader().clone();
            header.getSequence(0).setSequenceLength(firstContigLength);
            for (final SAMRecord record : reader) {
                record.setHeader(header);
                records.add(record);
            }
        }
        SvEvidenceDump.write(SvEvidenceDump.DIR.resolve(name), header, records);
    }

    static void find(final String label, final String input, final String... rest) throws Exception {
        SvEvidenceDump.find(label, input, rest);
    }

    static void pipeline(final String label, final String input, final String... rest) throws Exception {
        SvEvidenceDump.pipeline(label, input, rest);
    }

    /**
     * The BAMs the pipeline writes, through the JDK's deflater so that a host without GKL's level 2 can
     * replay them. These come first, as in SvEvidenceDump: a JVM's deflater stays Intel's once a run has set it.
     */
    static void jdkDeflaterRows() throws Exception {
        pipeline("split-pipeline", "reads.bam", "--use-jdk-deflater", "--bam-partition-size", "65536");
        pipeline("split-pipeline-small", "reads.bam", "--use-jdk-deflater", "--bam-partition-size", "7000",
                "--exp-interpret");
        pipeline("split-pipeline-threshold", "reads.bam", "--use-jdk-deflater", "--bam-partition-size", "20000",
                "--imprecise-variant-evidence-threshold", "1");
        pipeline("split-pipeline-other", "other.bam", "--use-jdk-deflater", "--bam-partition-size", "65536",
                "--exp-interpret");
        pipeline("split-pipeline-external", "reads.bam", "--use-jdk-deflater", "--bam-partition-size", "65536",
                "--external-evidence", at("external-long.bed"));
        pipeline("split-pipeline-cnv", "reads.bam", "--use-jdk-deflater", "--bam-partition-size", "65536",
                "--cnv-calls", at("cnv.vcf"));
    }

    static void rows() throws Exception {
        jdkDeflaterRows();
        for (final String size : new String[] {"131072", "65536", "20000", "7000", "1500"}) {
            find("split-" + size, "reads.bam", "--bam-partition-size", size);
        }
        // The slop: a tenth over the size is one split, and more than that is two.
        final long bytes = Files.size(SvEvidenceDump.DIR.resolve("reads.bam"));
        find("split-slop-under", "reads.bam", "--bam-partition-size", String.valueOf(bytes * 100 / 105));
        find("split-slop-over", "reads.bam", "--bam-partition-size", String.valueOf(bytes * 100 / 115));
        // Every 64 KiB block its own split or none, and splits with no block starting in them.
        find("split-300", "reads.bam", "--bam-partition-size", "300");
        // `(int) splitSize`: negative and zero leave the default, and the rest is what is left.
        find("split-2g", "reads.bam", "--bam-partition-size", "2147483648");
        find("split-4g", "reads.bam", "--bam-partition-size", "4294967296");
        find("split-4g-and-64k", "reads.bam", "--bam-partition-size", "4295032832");
        find("split-over-block", "reads.bam", "--bam-partition-size", "1073741824");
        find("split-nio", "reads.bam", "--bam-partition-size", "65536", "--use-nio");
        find("split-nio-small", "reads.bam", "--bam-partition-size", "7000", "--use-nio");
        find("split-nio-default", "reads.bam", "--use-nio");
        find("split-sbi", "sbi.bam", "--bam-partition-size", "65536");
        find("split-sbi-small", "sbi.bam", "--bam-partition-size", "7000");
        find("split-intervals", "reads.bam", "--bam-partition-size", "65536", "-L", "chr1:2000-12000");
        find("split-intervals-two", "reads.bam", "--bam-partition-size", "20000", "-L", "chr1:1-9000", "-L", "chr2");
        find("split-shallow", "shallow.bam", "--bam-partition-size", "20000");
        find("split-other", "other.bam", "--bam-partition-size", "65536");
        find("split-k31", "reads.bam", "--bam-partition-size", "65536", "--k-size", "31",
                "--kmers-to-ignore", at("kill31.txt"));
        find("split-cleaner", "reads.bam", "--bam-partition-size", "65536", "--cleaner-max-intervals", "1",
                "--cleaner-min-kmer-count", "6", "--cleaner-max-copy-number", "2");
        find("split-min-kmers", "reads.bam", "--bam-partition-size", "65536", "--min-kmers-per-interval", "1000000");
        find("split-thin", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "0.05",
                "--min-coherent-evidence-coverage-ratio", "0.02");
        find("split-coherent", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "5.0",
                "--min-coherent-evidence-coverage-ratio", "0.1");
        find("split-weight", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "5.0",
                "--min-coherent-evidence-coverage-ratio", "5.0");
        find("split-external", "reads.bam", "--bam-partition-size", "65536", "--external-evidence", at("external.bed"));
        find("split-external-long", "reads.bam", "--bam-partition-size", "20000",
                "--external-evidence", at("external-long.bed"));
        find("split-external-empty-partition", "reads.bam", "--bam-partition-size", "300",
                "--external-evidence", at("external.bed"));
        find("split-gaps", "reads.bam", "--bam-partition-size", "65536", "--exclusion-intervals", at("gaps.txt"),
                "--exclusion-interval-padding", "200");
        find("split-cross", "reads.bam", "--bam-partition-size", "65536", "--cross-contigs-to-ignore", at("cross.txt"));
        find("split-depth", "reads.bam", "--bam-partition-size", "65536", "--high-depth-coverage-peak-factor", "1",
                "--high-depth-coverage-factor", "1");
        find("split-interval-only", "reads.bam", "--bam-partition-size", "65536", "--interval-only-assembly");
        find("split-queryname", "reads.bam", "--bam-partition-size", "65536",
                "--assembled-contigs-output-order", "queryname");
        find("split-not-sorted", "not-sorted.bam", "--bam-partition-size", "7000");
        // KmerCleaner's table, from twice the dictionary over the partitions, narrowed to an int.
        find("table-wide", "wide.bam");
        find("table-wide-split", "wide.bam", "--bam-partition-size", "65536");
        find("table-wide-k31", "wide.bam", "--k-size", "31", "--kmers-to-ignore", at("kill31.txt"));
        find("table-wrapped", "wrapped.bam");
        find("table-wrapped-split", "wrapped.bam", "--bam-partition-size", "65536");
    }
}
