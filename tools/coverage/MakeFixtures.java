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
     * A coordinate-sorted BAM whose reads carry `N` in their cigars, which is what
     * `SplitNCigarReads` splits on, and one whose cigar has none.
     *
     * A read with k `N` elements becomes k+1 reads, so a file of plain `10M` reads comes out of
     * the tool unchanged and measures the traversal rather than the split. The `N` here is ten
     * bases of reference between two matched sections, which is a splice a reference of repeats
     * still supports.
     */
    static void spliced(final Path bam) {
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
        final String[] cigars = {"4M10N6M", "3M5N3M5N4M", "10M"};
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < cigars.length; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(5 + index * 700);
                record.setCigarString(cigars[index]);
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /**
     * A coordinate-sorted BAM that reads like bisulfite sequencing over the corpus's reference.
     *
     * The reference is `ACGT` repeated, so a C sits at every fourth base and a G beside it. The
     * tool counts, at a reference C, the FORWARD reads that kept the C against those that read T,
     * and at a reference G the REVERSE reads that kept the G against those that read A. Both
     * strands are here and both conversions, so the array has records to compare rather than an
     * empty VCF: four forward reads of which two are converted, and four reverse reads of which
     * two are.
     */
    static void methylation(final Path bam) {
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
                final boolean reverse = index >= 4;
                final boolean converted = index % 4 >= 2;
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(5);
                record.setCigarString("12M");
                record.setMappingQuality(60);
                // The reference from position five is ACGTACGTACGT. A converted forward read reads
                // T where the reference has C, and a converted reverse read reads A where it has G.
                final String bases;
                if (!converted) {
                    bases = "ACGTACGTACGT";
                } else if (reverse) {
                    bases = "ACATACATACAT";
                } else {
                    bases = "ATGTATGTATGT";
                }
                record.setReadString(bases);
                record.setBaseQualityString("IIIIIIIIIIII");
                record.setReadNegativeStrandFlag(reverse);
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /**
     * A coordinate-sorted BAM with TWO read groups, differing in sample and in library.
     *
     * `SplitReads` writes one file per key, so a file with a single read group is one file
     * whichever splitter is asked for, and an array over it would compare one output to itself.
     * Two groups make `--split-sample`, `--split-read-group` and `--split-library-name` each
     * produce two files, and their keys differ from one another.
     */
    static void twoGroups(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : new String[][] {
                {"rg1", "sample1", "lib1"}, {"rg2", "sample2", "lib2"}}) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample(group[1]);
            record.setLibrary(group[2]);
            record.setPlatformUnit("unit1");
            record.setPlatform("ILLUMINA");
            header.addReadGroup(record);
        }
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < 6; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(100 + index * 700);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("II##IIII##");
                record.setAttribute("RG", index % 2 == 0 ? "rg1" : "rg2");
                writer.addAlignment(record);
            }
        }
    }

    /**
     * A coordinate-sorted BAM whose reads carry INDELS inside the reference's repeat, which is
     * what `LeftAlignIndels` needs to move anything.
     *
     * The reference this corpus carries is `ACGT` repeated, so an indel of a whole four-base unit
     * can walk left through the repeat: a deletion that reaches the front of the read's window is
     * DROPPED and the read moves right by the bases it removed, and an insertion that reaches it
     * is kept where it is. Both branches are here, next to the two kinds of read that never reach
     * the call at all: one whose cigar is a single element, and an unmapped one.
     *
     * A file of plain `10M` reads would leave every row of the array with the input unchanged,
     * which is an array that measures the traversal and not the tool.
     */
    static void indels(final Path bam) {
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
        // Position 5 is an `A`, so a read starting there is in phase with the repeat and an indel
        // of one whole unit leaves the alignment just as good four bases to the left.
        final String[][] reads = {
                {"4M4D6M", "ACGTACGTAC"},
                {"4M4I6M", "ACGTACGTACGTAC"},
                {"10M", "ACGTACGTAC"},
        };
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < reads.length; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(5 + index * 700);
                record.setCigarString(reads[index][0]);
                record.setMappingQuality(60);
                record.setReadString(reads[index][1]);
                record.setBaseQualityString("I".repeat(reads[index][1].length()));
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
            final SAMRecord unmapped = new SAMRecord(header);
            unmapped.setReadName("HWI:1:FC:1:1:9:9");
            unmapped.setReadUnmappedFlag(true);
            unmapped.setReadString("ACGTACGTAC");
            unmapped.setBaseQualityString("IIIIIIIIII");
            unmapped.setAttribute("RG", "rg1");
            writer.addAlignment(unmapped);
        }
    }

    /**
     * A coordinate-sorted BAM whose every read carries `OQ`, which is what
     * `RevertBaseQualityScores` needs to do anything at all.
     *
     * That tool ABORTS on the first read without the tag rather than skipping it, so a corpus
     * holding only `reads.bam` measures the refusal and nothing else. Here every read has one, and
     * the original qualities DIFFER from the current ones -- `2` against `I`, which is quality two
     * against forty -- so a row that reverts is a different file rather than the same one.
     */
    static void bamWithOriginalQualities(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg4");
        group.setSample("sample4");
        group.setLibrary("lib4");
        group.setPlatformUnit("unit4");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (int index = 0; index < 6; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("OQ:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setFlags(0);
                record.setReferenceName("chr1");
                record.setAlignmentStart(150 + index * 800);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                // The original qualities the tool restores, and they are not the current ones.
                record.setAttribute("OQ", "##########");
                record.setAttribute("RG", "rg4");
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
        // The tabix index of that block-compressed VCF, which is the only `.tbi` in the corpus and
        // the only thing `DumpTabixIndex` can be given that it does not refuse. It is written by
        // the reference's own `IndexFeatureFile` rather than here, so what the array reads is an
        // index GATK produced.
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", dir.resolve("reads.vcf.gz").toString()});
        // A file NAMED like a tabix index and not gzipped at all. `DumpTabixIndex` checks the name
        // before it opens anything, so a `.vcf` handed to it is refused for its name and never
        // reaches the gzip layer; this one gets past the name and fails inside `java.util.zip`,
        // which is the second of the two refusals in the dump-tabix-index golden.
        Files.writeString(dir.resolve("plain.tbi"), vcf(), StandardCharsets.UTF_8);
        bam(dir.resolve("reads.bam"));
        bamTwo(dir.resolve("reads2.bam"));
        bamWithOriginalQualities(dir.resolve("reads_oq.bam"));
        indels(dir.resolve("indels.bam"));
        twoGroups(dir.resolve("groups.bam"));
        spliced(dir.resolve("spliced.bam"));
        methylation(dir.resolve("methyl.bam"));
        // The `-XF` file `ClipReads` reads: a FASTA of sequences to clip, which is a different
        // argument from `-X` and takes its names from the records rather than numbering them.
        Files.writeString(dir.resolve("clip.fasta"),
                ">adapterOne\nACGTACGT\n>adapterTwo\nTTTTGGGG\n", StandardCharsets.UTF_8);
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
        // The pileup summaries `CalculateContamination` reads, produced by the REFERENCE's own
        // `GetPileupSummaries` over the corpus. The chain is the point: one tool's output is the
        // other's input, so the second tool is measured on a table the first really writes rather
        // than on one this harness invented. A second table over `reads2.bam` gives
        // `--matched-normal` a value that is not the input.
        for (final String[] pair : new String[][] {
                {"reads.bam", "summaries.table"}, {"reads2.bam", "summaries2.table"}}) {
            new org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries()
                    .instanceMain(new String[] {
                            "--input", dir.resolve(pair[0]).toString(),
                            "--variant", population.toString(),
                            "--intervals", "chr1:1-60000",
                            "--output", dir.resolve(pair[1]).toString(),
                    });
        }
        System.out.println("wrote " + dir);
    }
}
