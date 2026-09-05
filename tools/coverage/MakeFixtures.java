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
import java.util.ArrayList;
import java.util.List;

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

    /**
     * A POPULATION VCF: biallelic SNPs carrying the `AF` info field `GetPileupSummaries` reads.
     *
     * The corpus's own `reads.vcf` declares no `AF` at all, and that tool refuses such a file
     * before its first locus, so an array built on it would measure one refusal on every row. This
     * one carries the field in its header and on every record, and the frequencies straddle the
     * tool's default window (0.01 to 0.2, both bounds STRICT): 0.005 is below it, 0.2 is exactly
     * the upper bound and therefore excluded, and the rest are inside.
     *
     * The records sit on the reads of `reads.bam`, one per read, because a site the reads do not
     * cover produces no pileup and therefore no row. The last is a triallelic site, which the tool
     * skips whatever its frequency is: without it nothing in the array distinguishes the
     * biallelic-SNP test from the frequency test.
     */
    static String populationVcf() {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        final double[] frequencies = {0.005, 0.02, 0.05, 0.1, 0.15, 0.2, 0.03};
        int index = 0;
        for (int position = 100; position <= 4300; position += 700) {
            text.append("chr1\t").append(position).append("\trs").append(position)
                    .append("\tA\tC\t100\tPASS\tAF=").append(frequencies[index]).append('\n');
            index++;
        }
        text.append("chr1\t5000\trs5000\tA\tC,G\t100\tPASS\tAF=0.05,0.03\n");
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

    /**
     * A second coordinate-sorted BAM, so that `--input` has two values rather than one.
     *
     * An argument with a single fixture value is held at it and no row can notice whether it
     * matters, which is the difference between an argument that is covered and one that is only
     * present. The reads differ from `reads.bam` in the three ways the corpus needs: a different
     * count, positions that fall on the other side of both interval fixtures, and two records that
     * the default read filters disagree about. `1D9M` is a well-formed cigar that
     * `GoodCigarReadFilter` refuses for its leading deletion, and the unmapped record is what
     * `MappedReadFilter` is there to remove; without them every filter fixture would be inert.
     */
    static void bamTwo(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg2");
        group.setSample("sample2");
        group.setLibrary("lib2");
        group.setPlatformUnit("unit2");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        final int[] starts = {200, 900, 1600, 50500, 51200};
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < starts.length; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:2:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setFlags(0);
                record.setReferenceName("chr1");
                record.setAlignmentStart(starts[index]);
                // One record carries the leading deletion, and it consumes nine read bases.
                final boolean clipped = index == 2;
                record.setCigarString(clipped ? "1D9M" : "10M");
                record.setMappingQuality(60);
                record.setReadString(clipped ? "ACGTACGTA" : "ACGTACGTAC");
                record.setBaseQualityString(clipped ? "IIIIIIIII" : "IIIIIIIIII");
                record.setAttribute("RG", "rg2");
                writer.addAlignment(record);
            }
            final SAMRecord unmapped = new SAMRecord(header);
            unmapped.setReadName("HWI:2:FC:1:1:9:9");
            unmapped.setReadUnmappedFlag(true);
            unmapped.setReadString("ACGTACGTAC");
            unmapped.setBaseQualityString("IIIIIIIIII");
            unmapped.setAttribute("RG", "rg2");
            writer.addAlignment(unmapped);
        }
    }

    /**
     * A BAM of PAIRS, which is what a tool asking for mate information needs.
     *
     * `PrintDistantMates` reads every record's mate, and the two BAMs above are unpaired: the
     * REFERENCE itself answers `Cannot get mate information for an unpaired read` on eight of that
     * tool's twenty-one rows, so a corpus without pairs cannot measure it at all. Three pairs, and
     * they differ in the one way the tool selects on: the first two mates sit beside each other,
     * the second pair straddles most of the contig, and the third is a pair whose mate is
     * unmapped.
     */
    static void pairs(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg3");
        group.setSample("sample3");
        group.setLibrary("lib3");
        group.setPlatformUnit("unit3");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        // start, mate start: near, far, and far again. The writer is PRESORTED, so the six
        // records are built first and written in coordinate order rather than pair by pair.
        final int[][] pairs = {{100, 300}, {1000, 60000}, {2000, 90000}};
        final java.util.List<SAMRecord> records = new java.util.ArrayList<>();
        for (int pair = 0; pair < pairs.length; pair++) {
            for (int end = 0; end < 2; end++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("PAIR:" + (pair + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(pairs[pair][end]);
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(pairs[pair][1 - end]);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                record.setAttribute("RG", "rg3");
                record.setReadPairedFlag(true);
                record.setProperPairFlag(pair == 0);
                record.setMateUnmappedFlag(false);
                record.setFirstOfPairFlag(end == 0);
                record.setSecondOfPairFlag(end == 1);
                records.add(record);
            }
        }
        records.sort(java.util.Comparator.comparingInt(SAMRecord::getAlignmentStart));
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            records.forEach(writer::addAlignment);
        }
    }

    /**
     * The samtools mpileup files `CheckPileup` compares GATK's own pileup against.
     *
     * The truth is not written by hand: it is the REFERENCE's own `Pileup` output, converted line
     * by line into samtools' six columns. The two tools share their five default read filters and
     * `reads.bam` carries no pairs, so nothing between them changes a base or a quality, and a
     * file built any other way would be a guess about the traversal rather than a record of it.
     *
     * Two files, because one is a value and not a covered argument: `truth.pileup` agrees with the
     * traversal at every locus and `wrong.pileup` disagrees at the first, so `--pileup` has a row
     * that validates and a row that is refused.
     *
     * Each is indexed the way the tool's own message tells the user to. A `FeatureInput` is queried
     * by interval, so without an index every run dies before it reads a locus.
     */
    static void pileups(final Path dir) throws Exception {
        final Path raw = dir.resolve("pileup.txt");
        new org.broadinstitute.hellbender.tools.walkers.qc.Pileup().instanceMain(new String[] {
                "--input", dir.resolve("reads.bam").toString(),
                "--reference", dir.resolve("reference.fasta").toString(),
                "--output", raw.toString(),
        });
        final List<String> lines = new ArrayList<>();
        for (final String line : Files.readAllLines(raw, StandardCharsets.UTF_8)) {
            if (line.isBlank()) {
                continue;
            }
            // `contig position referenceBase bases quals`, which is the same five fields samtools
            // writes with the DEPTH inserted before the bases.
            final String[] fields = line.split(" ");
            lines.add(String.join("\t", fields[0], fields[1], fields[2],
                    String.valueOf(fields[3].length()), fields[3], fields[4]));
        }
        Files.write(dir.resolve("truth.pileup"), lines, StandardCharsets.UTF_8);
        // One locus disagreeing, which is a `Bases not equal` and not a size or a location: the
        // three comparisons are ordered, and the one the array should reach is the deepest.
        final List<String> wrong = new ArrayList<>(lines);
        final String[] first = wrong.get(0).split("\t");
        first[4] = first[4].replace('A', 'T').replace('C', 'G');
        wrong.set(0, String.join("\t", first));
        Files.write(dir.resolve("wrong.pileup"), wrong, StandardCharsets.UTF_8);
        Files.delete(raw);
        for (final String label : new String[] {"truth", "wrong"}) {
            new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                    .instanceMain(new String[] {"-I", dir.resolve(label + ".pileup").toString()});
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
        bamTwo(dir.resolve("reads2.bam"));
        pairs(dir.resolve("pairs.bam"));
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
        // A recalibration table, produced by the REFERENCE rather than written here: ApplyBQSR
        // reads one and nothing else in the corpus can make a valid one. BaseRecalibrator is a
        // GATK tool like any other, so it is run the way the array runs a tool.
        final int status = new org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator()
                .instanceMain(new String[] {
                        "--input", dir.resolve("reads.bam").toString(),
                        "--reference", dir.resolve("reference.fasta").toString(),
                        // The INDEXED copy: known sites are queried by interval, so an
                        // unindexed VCF is refused before a read is looked at.
                        "--known-sites", dir.resolve("indexed.vcf").toString(),
                        "--output", dir.resolve("recal.table").toString(),
                }) == null ? 1 : 0;
        System.out.println("recalibrator status " + status);
        // A second table over the second BAM, for the same reason `reads2.bam` exists: with one
        // value `--bqsr-recal-file` is held at it and the argument is present rather than covered.
        final int otherStatus =
                new org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator()
                        .instanceMain(new String[] {
                                "--input", dir.resolve("reads2.bam").toString(),
                                "--reference", dir.resolve("reference.fasta").toString(),
                                "--known-sites", dir.resolve("indexed.vcf").toString(),
                                "--output", dir.resolve("recal2.table").toString(),
                        }) == null ? 1 : 0;
        System.out.println("second recalibrator status " + otherStatus);
        pileups(dir);
        // The population VCF, indexed: a `FeatureInput` is queried by interval, so an unindexed
        // one is refused before the traversal starts.
        final Path population = dir.resolve("population.vcf");
        Files.writeString(population, populationVcf(), StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", population.toString()});
        System.out.println("wrote " + dir);
    }
}
