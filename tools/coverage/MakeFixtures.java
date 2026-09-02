/*
 * The corpus a covering array is run over, written inside the pinned container.
 *
 * The fixtures are NOT committed: they are rebuilt from this program on every run, which keeps
 * them deterministic and keeps binary files out of the tree. Three of them are the three shapes
 * IndexFeatureFile's index kinds fall into: a plain VCF (linear index), a BED (linear index over a
 * different codec) and a block-compressed VCF (tabix). The fourth is a small coordinate-sorted BAM
 * with its index, which is what a read walker needs to run at all.
 *
 * Usage: MakeFixtures <directory>
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class MakeFixtures {

    static String vcf() {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        for (int position = 100; position <= 5000; position += 700) {
            text.append("chr1\t").append(position).append("\trs").append(position)
                    .append("\tA\tC\t100\tPASS\t.\tGT\t0/1\n");
        }
        return text.toString();
    }

    static String bed() {
        final StringBuilder text = new StringBuilder();
        for (int start = 100; start <= 5000; start += 700) {
            text.append("chr1\t").append(start).append('\t').append(start + 50)
                    .append("\tregion").append(start).append('\n');
        }
        return text.toString();
    }

    /** A small coordinate-sorted BAM: eight reads on one contig, one of them a duplicate. */
    static void bam(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < 8; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setFlags(index == 7 ? 0x400 : 0);
                record.setReferenceName("chr1");
                record.setAlignmentStart(100 + index * 700);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        // The deflater is pinned exactly as the oracle contract pins it for goldens: a fixture
        // that is not byte-reproducible would make a coverage measurement unrepeatable.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final Path dir = Paths.get(args[0]);
        Files.createDirectories(dir);
        Files.writeString(dir.resolve("reads.vcf"), vcf(), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("regions.bed"), bed(), StandardCharsets.UTF_8);
        try (final OutputStream out =
                     new BlockCompressedOutputStream(dir.resolve("reads.vcf.gz").toFile())) {
            out.write(vcf().getBytes(StandardCharsets.UTF_8));
        }
        bam(dir.resolve("reads.bam"));
        // The same VCF with a Tribble index beside it. A feature walker refuses `-L` against an
        // input with no random access, so an array whose only VCF were unindexed would compare two
        // refusals on every interval row and never reach a traversal.
        final Path indexed = dir.resolve("indexed.vcf");
        Files.writeString(indexed, vcf(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        indexed, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("indexed.vcf.idx"));
        // A reference, with the .fai and .dict beside it that GATK requires: the writer names all
        // three, so the naming is htsjdk's rather than this harness's.
        try (final htsjdk.samtools.reference.FastaReferenceWriter reference =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(dir.resolve("reference.fasta"))
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            final StringBuilder bases = new StringBuilder();
            for (int i = 0; i < 100000; i++) {
                bases.append("ACGT".charAt(i % 4));
            }
            reference.startSequence("chr1").appendBases(bases.toString());
        }
        // A second reference on a contig the corpus does not carry, so `--reference` has a value
        // that agrees with the reads and a value that does not.
        try (final htsjdk.samtools.reference.FastaReferenceWriter other =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(dir.resolve("other.fasta"))
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            final StringBuilder bases = new StringBuilder();
            for (int i = 0; i < 1000; i++) {
                bases.append("ACGT".charAt(i % 4));
            }
            other.startSequence("chrOther").appendBases(bases.toString());
        }

        // Two sequence dictionaries for `--sequence-dictionary`: one that agrees with the corpus's
        // own contig and one that shares nothing with it, so the argument has a row that is
        // accepted and a row that is refused.
        Files.writeString(dir.resolve("matching.dict"),
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("other.dict"),
                "@HD\tVN:1.6\n@SQ\tSN:chrOther\tLN:1000\n", StandardCharsets.UTF_8);
        System.out.println("wrote " + dir);
    }
}
