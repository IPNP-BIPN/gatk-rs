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
import java.util.Comparator;
import java.util.List;

public class MakeFixtures {

    /**
     * A VCF of INDELS at measured distances, which is what `RemoveNearbyIndels` needs.
     *
     * Every other VCF in this corpus holds SNPs alone, and a run over one of those emits every
     * record whatever the spacing is: the array would measure the traversal and not the tool. The
     * records here are, in order: an isolated indel; a PAIR ten bases apart, which any spacing above
     * ten removes; a SNP between them, which survives its neighbours being dropped; a RUN of three
     * indels, which the buffer loses whole because it measures the next one against an indel it has
     * already thrown away; and a last indel a hundred bases past the run, which a spacing of 200
     * takes with it and a spacing of 10 does not.
     */
    static String indelVcf() {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        final int[][] records = {
                // position, kind: 0 is a deletion, 1 an insertion, 2 a snp
                {1000, 0}, {2000, 1}, {2010, 0}, {2015, 2}, {3000, 1}, {3005, 0}, {3010, 1},
                {3110, 0},
        };
        for (final int[] record : records) {
            final String reference = record[1] == 0 ? "ACGT" : "A";
            final String alternate = record[1] == 0 ? "A" : (record[1] == 1 ? "ACGT" : "C");
            text.append("chr1\t").append(record[0]).append("\trs").append(record[0])
                    .append('\t').append(reference).append('\t').append(alternate)
                    .append("\t100\tPASS\t.\tGT\t0/1\n");
        }
        return text.toString();
    }

    /**
     * A TWO-sample VCF whose records are singleton hets, alternating which sample carries them.
     *
     * `CalculateMixingFractions` fills one bucket per sample and then divides each bucket's alt
     * fraction by the SUM of every sample's, so a one-sample file can only ever answer `1.0` or
     * `NaN`. With two samples the table has two rows that add up, and their ORDER is a
     * `HashMap`'s iteration order rather than the header's, which is the property the tool's own
     * golden pins and which no command line could reach while the corpus had one sample.
     */
    static String duoVcf() {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n");
        // The sites sit one base INTO each read rather than at its start. `reads.bam`'s reads are
        // `ACGTACGTAC` beginning at 100, 800, and so on, so the base at the read's own start is the
        // reference `A` and the base one further in is `C`, which is this file's alternate. A site
        // at the start counts a total and no alt, every fraction is `0/0`, and the table is NaN
        // whatever else the row says: the first version of this fixture did exactly that.
        //
        // The hets are dealt five to `sample1` and three to `sample2`, so the two mixing fractions
        // are different numbers rather than one number twice.
        final boolean[] toFirstSample = {true, true, false, true, false, true, false, true};
        int index = 0;
        for (int position = 101; position <= 5001; position += 700) {
            text.append("chr1\t").append(position).append("\trs").append(position)
                    .append("\tA\tC\t100\tPASS\t.\tGT\t")
                    .append(toFirstSample[index] ? "0/1\t0/0" : "0/0\t0/1").append('\n');
            index++;
        }
        return text.toString();
    }

    /**
     * A copy-ratio segment file, which is what `CallCopyRatioSegments` reads.
     *
     * Written by hand rather than produced by `ModelSegments`, which would need read counts and
     * allelic counts this corpus does not carry. The format is htsjdk's SAM header, a column line
     * and the rows, and the reference validates all three when it reads the file: a header without
     * the `@RG` line's sample name, or a column line that is not this one, is refused there rather
     * than accepted quietly.
     *
     * The segments straddle the copy-neutral window on purpose: two sit inside it, one is far
     * below and one far above, so the calls are `0`, `-` and `+` rather than one letter repeated,
     * and the length-weighted statistics have something to weigh.
     */
    static String copyRatioSegments() {
        final StringBuilder text = new StringBuilder();
        text.append("@HD\tVN:1.6\n");
        text.append("@SQ\tSN:chr1\tLN:100000\n");
        text.append("@RG\tID:GATKCopyNumber\tSM:sample1\n");
        text.append("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n");
        final int[][] segments = {{1, 10000, 100}, {10001, 20000, 120}, {20001, 30000, 90},
                {30001, 40000, 80}, {40001, 50000, 110}};
        final double[] means = {0.01, -0.02, -1.5, 1.2, 0.03};
        for (int index = 0; index < segments.length; index++) {
            text.append("chr1\t").append(segments[index][0]).append('\t')
                    .append(segments[index][1]).append('\t').append(segments[index][2])
                    .append('\t').append(means[index]).append('\n');
        }
        return text.toString();
    }

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

    /**
     * `population.vcf` with every allele frequency moved, which is what gives
     * `EvaluateInfoFieldConcordance` a difference to average.
     *
     * A file compared with itself produces a delta of zero at every true positive, so the mean and
     * the standard deviation are zero whatever the arithmetic does: the first version of that
     * tool's array had exactly that, and measured the walk rather than the numbers. The shift is
     * not uniform, because a constant offset would make the standard deviation zero as well.
     */
    static String shiftedPopulationVcf() {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
        text.append("##contig=<ID=chr1,length=100000>\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        final double[] frequencies = {0.01, 0.03, 0.04, 0.2, 0.1, 0.25, 0.02};
        int index = 0;
        for (int position = 100; position <= 4300; position += 700) {
            text.append("chr1\t").append(position).append("\trs").append(position)
                    .append("\tA\tC\t100\tPASS\tAF=").append(frequencies[index]).append('\n');
            index++;
        }
        text.append("chr1\t5000\trs5000\tA\tC,G\t100\tPASS\tAF=0.07,0.01\n");
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
     * The same BAM, with an `M5` on its one `@SQ` line taken from `reference.fasta`'s own
     * dictionary.
     *
     * `CheckReferenceCompatibility` takes one of two paths depending on a single property of its
     * input: with an MD5 on EVERY sequence it compares bases through `CompareReferences`' table,
     * and without one it compares names and lengths alone and says so in every summary. No BAM in
     * this corpus carries an M5, so without this one the first path is unreachable from a command
     * line, and adding it to `reads.bam` would change a header that several goldens print.
     */
    static void md5Bam(final Path bam, final Path reference) {
        final SAMFileHeader header = new SAMFileHeader();
        // The reference's own `.dict`, M5 included, which is what makes this BAM's dictionary
        // agree with `reference.fasta` base for base rather than by name alone.
        header.setSequenceDictionary(htsjdk.samtools.reference.ReferenceSequenceFileFactory
                .getReferenceSequenceFile(reference).getSequenceDictionary());
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
                // Every row of the array excludes an interval near the start of the contig, so
                // reads placed there are filtered out and the tool writes a header and no record.
                // 1005 is in phase with the ACGT repeat and outside both excluded ranges.
                record.setAlignmentStart(1005);
                record.setCigarString("12M");
                record.setMappingQuality(60);
                // The reference from position 1005 is ACGTACGTACGT. A converted forward read reads
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
     * A coordinate-sorted BAM whose reads PILE UP: eight of them at one locus, four at the next.
     *
     * `--max-depth-per-sample` thins a pileup deeper than its target, and every other file in this
     * corpus has a depth of one, so the argument decided nothing and could not be measured. Eight
     * at one position is deep enough for the reservoir to draw and the leveller to level at any
     * target the array uses.
     *
     * The bases differ read by read, which is what makes a thinned pileup VISIBLE: eight identical
     * reads would print the same column whichever four survived.
     */
    static void deep(final Path bam) {
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
            final String[] bases = {
                    "ACGTACGTAC", "ACGTACGTAG", "ACGTACGTAT", "ACGTACGTAA",
                    "CCGTACGTAC", "GCGTACGTAC", "TCGTACGTAC", "ACGTACGTCC",
            };
            for (int index = 0; index < bases.length; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:1:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                // Eight at 1005 and four more at 1105, so a run sees one deep locus and one that is
                // deep for half as long: the leveller's plan depends on the stacks it is given.
                record.setAlignmentStart(index < 8 ? 1005 : 1105);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString(bases[index]);
                record.setBaseQualityString("IIIIIIIIII");
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
            for (int index = 0; index < 4; index++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("HWI:1:FC:1:2:" + (index + 1) + ":" + (index + 1));
                record.setReferenceName("chr1");
                record.setAlignmentStart(1105);
                record.setCigarString("10M");
                record.setMappingQuality(60);
                record.setReadString(bases[index]);
                record.setBaseQualityString("IIIIIIIIII");
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
    /**
     * One coordinate-sorted BAM carrying TWO samples over the same locus, which is the shape
     * `GetNormalArtifactData` reads: the split is by sample name and not by file.
     *
     * Six reads per sample, forty bases each, all starting at chr1:101 so that every read covers
     * every locus of the window. Forty rather than ten because Mutect2's chain refuses a read
     * shorter than thirty, and one start rather than six because a locus the normal does not reach
     * has no alternate and produces nothing.
     *
     * Two sites carry an alternate, and they are different on purpose:
     *
     *   chr1:101  one normal read and two tumour reads carry `C`, so the tumour p-value is tiny,
     *             the keep probability is all but one and the locus becomes a row;
     *   chr1:121  one normal read carries `G` and no tumour read does, so the p-value is one, the
     *             keep probability falls to its floor of 0.05, and whether the locus survives is
     *             the seeded draw's answer rather than the counts'.
     *
     * Every other locus matches the reference in both samples, which is what leaves the table
     * short enough to read.
     */
    /** One record of `truth.vcf` or `calls.vcf`. */
    static String concordanceRecord(final int position, final String reference, final String alternate,
                                    final String filter) {
        return String.format("chr1\t%d\t.\t%s\t%s\t50\t%s\t.\tGT\t0/1%n", position, reference,
                alternate, filter).replace(System.lineSeparator(), "\n");
    }

    /** The corpus reference's bases at a one-based position: `ACGT`, repeated. */
    static String referenceBases(final int position, final int length) {
        final StringBuilder bases = new StringBuilder();
        for (int index = 0; index < length; index++) {
            bases.append("ACGT".charAt((position + index - 1) % 4));
        }
        return bases.toString();
    }

    /** One heterozygous record of `shiftable.vcf`. */
    static String shiftableRecord(final int position, final String reference, final String alternate) {
        return String.format("chr1\t%d\t.\t%s\t%s\t100\tPASS\t.\tGT\t0/1%n", position, reference,
                alternate).replace(System.lineSeparator(), "\n");
    }

    static void tumorAndNormal(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String[] group : new String[][] {
                {"rgn", "normal", "libn"}, {"rgt", "tumor", "libt"}}) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample(group[1]);
            record.setLibrary(group[2]);
            record.setPlatformUnit("unit1");
            record.setPlatform("ILLUMINA");
            header.addReadGroup(record);
        }
        // The reference repeats `ACGT`, and position 101 is an `A`, so a read of `ACGT` ten times
        // over matches it base for base.
        final String matching = "ACGT".repeat(10);
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            for (final String[] read : new String[][] {
                    {"n0", "rgn", matching},
                    {"n1", "rgn", matching},
                    {"n2", "rgn", matching},
                    {"n3", "rgn", matching},
                    {"n4", "rgn", "C" + matching.substring(1)},
                    {"n5", "rgn", matching.substring(0, 20) + "G" + matching.substring(21)},
                    {"t0", "rgt", matching},
                    {"t1", "rgt", matching},
                    {"t2", "rgt", matching},
                    {"t3", "rgt", matching},
                    {"t4", "rgt", "C" + matching.substring(1)},
                    {"t5", "rgt", "C" + matching.substring(1)}}) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(read[0]);
                record.setReferenceName("chr1");
                record.setAlignmentStart(101);
                record.setCigarString("40M");
                record.setMappingQuality(60);
                record.setReadString(read[2]);
                record.setBaseQualityString("I".repeat(40));
                record.setAttribute("RG", read[1]);
                writer.addAlignment(record);
            }
        }
    }

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
    /**
     * A QUERY-NAME sorted BAM of pairs, which is what the two tools that are no walkers require.
     *
     * `PostProcessReadsForRSEM` refuses anything else outright, and `TransferReadTags` refuses an
     * aligned file whose header does not say `SO:queryname`, so every other BAM in this corpus
     * reaches one line of either tool and stops. The four groups here are the four answers
     * `passesRSEMFilter` gives:
     *
     *   - `PAIR:1` is a proper pair of single-`M` reads, which passes, and it carries a pair of
     *     SECONDARY alignments whose mate positions point at each other, which is the one shape
     *     `groupSecondaryReads` keeps;
     *   - `PAIR:2`'s second read is unmapped, which is the `notBothMapped` count;
     *   - `PAIR:3`'s first read has an insertion in its cigar, which is the `unsupportedCigar`
     *     count: RSEM takes one `M` element and nothing else;
     *   - and `PAIR:4` is a first-of-pair with no second, which warns and is dropped.
     *
     * A group holding only a SECOND-of-pair is deliberately absent: it dereferences null in the
     * reference and would end every row of both arrays at the same line.
     *
     * There is no index. A queryname-sorted BAM cannot have one, which is also why the writer is
     * not asked for it here.
     */
    static void queryNameSorted(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, bam.toFile())) {
            // PAIR:1, the group that passes: two primaries and two secondaries.
            writer.addAlignment(mate(header, "PAIR:1", 100, 300, "10M", true, false, false));
            writer.addAlignment(mate(header, "PAIR:1", 500, 700, "10M", true, true, false));
            writer.addAlignment(mate(header, "PAIR:1", 300, 100, "10M", false, false, false));
            writer.addAlignment(mate(header, "PAIR:1", 700, 500, "10M", false, true, false));
            // PAIR:2, whose second read is unmapped.
            writer.addAlignment(mate(header, "PAIR:2", 1000, 1200, "10M", true, false, false));
            writer.addAlignment(mate(header, "PAIR:2", 1200, 1000, "10M", false, false, true));
            // PAIR:3, whose first read carries an insertion.
            writer.addAlignment(mate(header, "PAIR:3", 2000, 2200, "5M1I4M", true, false, false));
            writer.addAlignment(mate(header, "PAIR:3", 2200, 2000, "10M", false, false, false));
            // PAIR:4, a first of pair with no second.
            writer.addAlignment(mate(header, "PAIR:4", 3000, 3200, "10M", true, false, false));
        }
    }

    /** One record of {@link #queryNameSorted}, with the flags the group it belongs to needs. */
    static SAMRecord mate(final SAMFileHeader header, final String name, final int start,
                          final int mateStart, final String cigar, final boolean first,
                          final boolean secondary, final boolean unmapped) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadString("ACGTACGTAC");
        record.setBaseQualityString("IIIIIIIIII");
        record.setAttribute("RG", "rg1");
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(first);
        record.setSecondOfPairFlag(!first);
        if (unmapped) {
            record.setReadUnmappedFlag(true);
            record.setMateReferenceName("chr1");
            record.setMateAlignmentStart(mateStart);
            return record;
        }
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setMateReferenceName("chr1");
        record.setMateAlignmentStart(mateStart);
        record.setProperPairFlag(!secondary);
        record.setSecondaryAlignment(secondary);
        return record;
    }

    /**
     * The unmapped, query-name sorted file `TransferReadTags` copies tags FROM.
     *
     * One record per query name of {@link #queryNameSorted}, each carrying `RX` and none carrying
     * `MI`: the tool asks for the tags named on its command line and refuses a read whose value is
     * absent, so the pair of tag names is a row that answers and a row that refuses.
     *
     * The names are a SUPERSET of the aligned file's on purpose. The traversal plays this file
     * forward until it catches up with the aligned read, so a name here that the aligned file does
     * not carry is skipped, while the reverse is the tool's `IllegalStateException`.
     */
    /**
     * Reads grouped by molecule as fgbio leaves them, for `DownsampleByDuplicateSet`: `MI` is
     * `<molecule>/<strand>`, the file is ordered by molecule and its header says unsorted, so the
     * writer keeps the order it is given. Molecule 1 has two reads, one per strand; 2 has four; 3
     * has three, which the walker rejects as unpaired; 4 has two reads on one strand, which a
     * per-strand minimum of one rejects; 5..9 have two each, so the seeded draws have something to
     * keep and drop; and 10, the last, has three, which escapes every rule because the last set is
     * never asked. With `backwards` the molecules 4 and 5 swap places, which the walker refuses.
     */
    static void moleculeGroups(final Path bam, final boolean backwards) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        final String[][] molecules = {
                {"1", "A", "B"}, {"2", "A", "A", "B", "B"}, {"3", "A", "B", "A"}, {"4", "A", "A"},
                {"5", "A", "B"}, {"6", "A", "B"}, {"7", "B", "A"}, {"8", "A", "B"}, {"9", "A", "B"},
                {"10", "A", "B", "B"}};
        if (backwards) {
            final String[] fourth = molecules[3];
            molecules[3] = molecules[4];
            molecules[4] = fourth;
        }
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, bam.toFile())) {
            int serial = 0;
            for (final String[] molecule : molecules) {
                for (int i = 1; i < molecule.length; i++) {
                    final SAMRecord record = new SAMRecord(header);
                    record.setReadName("read" + (serial++));
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(1000 + 10 * Integer.parseInt(molecule[0]));
                    record.setMappingQuality(60);
                    record.setCigarString("10M");
                    record.setReadString("ACGTACGTAC");
                    record.setBaseQualityString("IIIIIIIIII");
                    record.setAttribute("RG", "rg1");
                    record.setAttribute("MI", molecule[0] + "/" + molecule[i]);
                    writer.addAlignment(record);
                }
            }
        }
    }

    /**
     * Reads over four sites for `ASEReadCounter`, on `reference.fasta`, whose base at position p is
     * "ACGT"[(p - 1) % 4]. Every read is ten bases from four before its site, the reference's own
     * bases except the one AT the site. At 1001 (ref A, alt G): two reads of A and two of G, one A
     * at mapping quality 10, one G at base quality 2, a C that is neither allele, a proper pair
     * whose mates agree on A and one whose mates disagree (A and G), and an improper pair. 1102
     * (ref C, alt T) carries three reads; 1203 is multiallelic and 1304 homozygous in the VCF.
     */
    static void aseBam(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        // name, site, base at the site, mapping quality, base quality at the site, pairing
        // (0 unpaired, 1 proper first, 2 proper second, 3 improper first)
        final Object[][] reads = {
                {"a1", 1001, 'A', 60, 30, 0}, {"a2", 1001, 'A', 60, 30, 0},
                {"g1", 1001, 'G', 60, 30, 0}, {"g2", 1001, 'G', 60, 30, 0},
                {"lowmq", 1001, 'A', 10, 30, 0}, {"lowbq", 1001, 'G', 60, 2, 0},
                {"other", 1001, 'C', 60, 30, 0},
                {"pairsame", 1001, 'A', 60, 30, 1}, {"pairsame", 1001, 'A', 60, 25, 2},
                {"pairdiff", 1001, 'A', 60, 30, 1}, {"pairdiff", 1001, 'G', 60, 25, 2},
                {"improper", 1001, 'G', 60, 30, 3},
                {"c1", 1102, 'C', 60, 30, 0}, {"t1", 1102, 'T', 60, 30, 0},
                {"t2", 1102, 'T', 60, 30, 0},
                {"m1", 1203, 'A', 60, 30, 0}, {"h1", 1304, 'C', 60, 30, 0}};
        final List<SAMRecord> records = new ArrayList<>();
        for (final Object[] read : reads) {
            final int site = (Integer) read[1];
            final int start = site - 4;
            final StringBuilder bases = new StringBuilder();
            final StringBuilder quals = new StringBuilder();
            for (int p = start; p < start + 10; p++) {
                bases.append(p == site ? (char) (Character) read[2] : "ACGT".charAt((p - 1) % 4));
                quals.append(p == site ? (char) (33 + (Integer) read[4]) : 'I');
            }
            final SAMRecord record = new SAMRecord(header);
            record.setReadName((String) read[0]);
            record.setReferenceName("chr1");
            record.setAlignmentStart(start);
            record.setCigarString("10M");
            record.setMappingQuality((Integer) read[3]);
            record.setReadString(bases.toString());
            record.setBaseQualityString(quals.toString());
            record.setAttribute("RG", "rg1");
            final int pairing = (Integer) read[5];
            if (pairing != 0) {
                record.setReadPairedFlag(true);
                record.setProperPairFlag(pairing != 3);
                record.setFirstOfPairFlag(pairing != 2);
                record.setSecondOfPairFlag(pairing == 2);
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(start);
                record.setMateNegativeStrandFlag(pairing == 1);
                record.setReadNegativeStrandFlag(pairing == 2);
            }
            records.add(record);
        }
        records.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true,
                             bam.toFile())) {
            records.forEach(writer::addAlignment);
        }
    }

    static void umi(final Path bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, bam.toFile())) {
            final String[] names = {"PAIR:1", "PAIR:2", "PAIR:3", "PAIR:4", "PAIR:5"};
            for (final String name : names) {
                for (final boolean first : new boolean[] {true, false}) {
                    final SAMRecord record = new SAMRecord(header);
                    record.setReadName(name);
                    record.setReadUnmappedFlag(true);
                    record.setMateUnmappedFlag(true);
                    record.setReadPairedFlag(true);
                    record.setFirstOfPairFlag(first);
                    record.setSecondOfPairFlag(!first);
                    record.setReadString("ACGTACGTAC");
                    record.setBaseQualityString("IIIIIIIIII");
                    record.setAttribute("RG", "rg1");
                    record.setAttribute("RX", "ACG-TGC");
                    writer.addAlignment(record);
                }
            }
        }
    }

    /**
     * `reads.bam`'s eight reads, at the same positions and under the same NAMES, with a different
     * quality array.
     *
     * `CompareBaseQualities` walks two files in lockstep and refuses the pair as soon as two reads
     * disagree by name, so a comparison needs two files that hold the same reads. Every other pair
     * in this corpus differs by name on the first record, which is one refusal on every row and no
     * comparison at all. The qualities here are `I` where reads.bam has `I` on six bases and `#` on
     * four, so the matrix has off-diagonal entries: the tool returns 1 rather than 0, and
     * `--throw-on-diff` turns that into a refusal.
     */
    static void requalified(final Path bam) {
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
                record.setBaseQualityString("IIIIII####");
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

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
        // The known sites `BaseRecalibrator` reads, as a BED naming the same loci the population
        // VCF does. Both formats are registered for that argument and the reference's own two runs
        // over them produce the same table, so the pair is what makes the argument's two values
        // comparable rather than merely different.
        final StringBuilder bed = new StringBuilder();
        for (int position = 100; position <= 4300; position += 700) {
            // A BED is half-open and zero-based, so the same one-based locus starts one lower.
            bed.append("chr1\t").append(position - 1).append('\t').append(position).append('\n');
        }
        Files.writeString(dir.resolve("known.bed"), bed.toString(), StandardCharsets.UTF_8);
        // Indexed, because `--known-sites` is QUERIED by interval: an unindexed file is refused
        // with `must support random access to enable queries by interval`, and a corpus that
        // carried one would compare two refusals rather than two tables.
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", dir.resolve("known.bed").toString()});
        // The annotation `GtfToBed` reads, which is the gtf-to-bed golden's own: the
        // reference's Gencode codec refuses anything less than a full one -- a hand-written
        // file of gene and transcript lines is `Decoded feature is not valid: null`, because
        // every transcript needs its exon line and its type, name and havana attributes.
        Files.writeString(dir.resolve("annotation.gtf"),
                "chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id \"GENE_B.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\ttranscript\t50\t250\t.\t+\t.\tgene_id \"GENE_B.1\"; transcript_id \"TX_B1.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TX_B1.1\"; level 2; tag \"basic\"; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\texon\t50\t250\t.\t+\t.\tgene_id \"GENE_B.1\"; transcript_id \"TX_B1.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TX_B1.1\"; exon_number 1; exon_id \"TX_B1.1.1\"; level 2;\n"
                        + "chr1\tHAVANA\ttranscript\t120\t180\t.\t+\t.\tgene_id \"GENE_B.1\"; transcript_id \"TX_B2.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TX_B2.1\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\texon\t120\t180\t.\t+\t.\tgene_id \"GENE_B.1\"; transcript_id \"TX_B2.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TX_B2.1\"; exon_number 1; exon_id \"TX_B2.1.1\"; level 2;\n"
                        + "chr1\tHAVANA\tgene\t300\t400\t.\t+\t.\tgene_id \"GENE_A.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\ttranscript\t300\t400\t.\t+\t.\tgene_id \"GENE_A.1\"; transcript_id \"TX_A1.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TX_A1.1\"; level 2; tag \"basic\"; tag \"basic\"; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id \"GENE_A.1\"; transcript_id \"TX_A1.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TX_A1.1\"; exon_number 1; exon_id \"TX_A1.1.1\"; level 2;\n"
                        + "chr1\tHAVANA\tgene\t300\t400\t.\t+\t.\tgene_id \"GENE_C.1\"; gene_type \"protein_coding\"; gene_name \"gamma\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\ttranscript\t300\t500\t.\t+\t.\tgene_id \"GENE_C.1\"; transcript_id \"TX_C1.1\"; gene_type \"protein_coding\"; gene_name \"gamma\"; transcript_type \"protein_coding\"; transcript_name \"TX_C1.1\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr1\tHAVANA\texon\t300\t500\t.\t+\t.\tgene_id \"GENE_C.1\"; transcript_id \"TX_C1.1\"; gene_type \"protein_coding\"; gene_name \"gamma\"; transcript_type \"protein_coding\"; transcript_name \"TX_C1.1\"; exon_number 1; exon_id \"TX_C1.1.1\"; level 2;\n"
                        + "chr2\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id \"GENE_D.1\"; gene_type \"protein_coding\"; gene_name \"delta\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr2\tHAVANA\ttranscript\t10\t20\t.\t+\t.\tgene_id \"GENE_D.1\"; transcript_id \"TX_D1.1\"; gene_type \"protein_coding\"; gene_name \"delta\"; transcript_type \"protein_coding\"; transcript_name \"TX_D1.1\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n"
                        + "chr2\tHAVANA\texon\t10\t20\t.\t+\t.\tgene_id \"GENE_D.1\"; transcript_id \"TX_D1.1\"; gene_type \"protein_coding\"; gene_name \"delta\"; transcript_type \"protein_coding\"; transcript_name \"TX_D1.1\"; exon_number 1; exon_id \"TX_D1.1.1\"; level 2;\n",
                StandardCharsets.UTF_8);
        // The dictionary that annotation is sorted by, which is the golden's: both contigs,
        // so the corpus's own matching.dict is the value that refuses the chr2 gene.
        Files.writeString(dir.resolve("gtf.dict"),
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1040\n@SQ\tSN:chr2\tLN:1040\n",
                StandardCharsets.UTF_8);
        bam(dir.resolve("reads.bam"));
        bamTwo(dir.resolve("reads2.bam"));
        bamWithOriginalQualities(dir.resolve("reads_oq.bam"));
        indels(dir.resolve("indels.bam"));
        twoGroups(dir.resolve("groups.bam"));
        tumorAndNormal(dir.resolve("tumor_normal.bam"));
        deep(dir.resolve("deep.bam"));
        spliced(dir.resolve("spliced.bam"));
        methylation(dir.resolve("methyl.bam"));
        // The `-XF` file `ClipReads` reads: a FASTA of sequences to clip, which is a different
        // argument from `-X` and takes its names from the records rather than numbering them.
        Files.writeString(dir.resolve("clip.fasta"),
                ">adapterOne\nACGTACGT\n>adapterTwo\nTTTTGGGG\n", StandardCharsets.UTF_8);
        pairs(dir.resolve("pairs.bam"));
        queryNameSorted(dir.resolve("qname.bam"));
        umi(dir.resolve("umi.bam"));
        aseBam(dir.resolve("ase.bam"));
        moleculeGroups(dir.resolve("molecules.bam"), false);
        moleculeGroups(dir.resolve("molecules_backwards.bam"), true);
        requalified(dir.resolve("requal.bam"));
        // The indel VCF, INDEXED: a variant walker refuses `-L` over an input with no random
        // access, so an unindexed one would answer a refusal on every interval row.
        final Path indels = dir.resolve("indels.vcf");
        Files.writeString(indels, indelVcf(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        indels, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("indels.vcf.idx"));
        // The same VCF with a Tribble index beside it. A feature walker refuses `-L` against an
        // input with no random access, so an array whose only VCF were unindexed would compare two
        // refusals on every interval row and never reach a traversal.
        // The two-sample VCF, indexed for the same reason: a variant walker refuses `-L` over an
        // input with no random access.
        Files.writeString(dir.resolve("segments.cr.seg"), copyRatioSegments(),
                StandardCharsets.UTF_8);
        final Path duo = dir.resolve("duo.vcf");
        Files.writeString(duo, duoVcf(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        duo, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("duo.vcf.idx"));
        // The mixing fractions `AnnotateVcfWithExpectedAlleleFraction` reads, produced by the
        // REFERENCE's own `CalculateMixingFractions` over the two-sample VCF and the corpus's
        // reads. The chain is the point: one tool's output is the other's input, so the second is
        // measured on a table the first really writes.
        new org.broadinstitute.hellbender.tools.walkers.validation.CalculateMixingFractions()
                .instanceMain(new String[] {
                        "--variant", dir.resolve("duo.vcf").toString(),
                        "--input", dir.resolve("reads.bam").toString(),
                        "--intervals", "chr1:1-6000",
                        "--output", dir.resolve("mixing.table").toString(),
                });
        // A GVCF: two reference BLOCKS with an `END` and a variant between them, all carrying
        // `<NON_REF>`. `ValidateVariants --validate-GVCF` needs one, and it needs the blocks to
        // stop short of the contig: the coverage check counts every locus no record covers, so a
        // file over chr1:1-1000 and an interval of chr1:1-6000 leave a gap the message names.
        // The reference bases are the corpus's own repeat, so the REF check passes on every row.
        final Path blocks = dir.resolve("blocks.g.vcf");
        Files.writeString(blocks,
                "##fileformat=VCFv4.2\n"
                        + "##contig=<ID=chr1,length=100000>\n"
                        + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele\">\n"
                        + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n"
                        + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
                        + "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=500\tGT\t0/0\n"
                        + "chr1\t501\t.\tA\tC,<NON_REF>\t50\t.\t.\tGT\t0/1\n"
                        + "chr1\t502\t.\tC\t<NON_REF>\t.\t.\tEND=1000\tGT\t0/0\n",
                StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        blocks, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("blocks.g.vcf.idx"));
        // Indels that can MOVE. The corpus reference is `ACGT` repeated, so a deletion or an
        // insertion of one whole repeat unit is equivalent at every offset of the repeat, and
        // `LeftAlignAndTrimVariants` walks it left as far as its window allows. The corpus's own
        // indels.vcf cannot show that: its alleles do not match the reference, so nothing moves and
        // an array over it measures the traversal rather than the alignment.
        //
        // Five records, and each is a different branch: a deletion that walks, an insertion that
        // walks, a second deletion close behind the first so that the distance to the record
        // already written is what bounds it, a deletion longer than the default
        // `--max-indel-length` which is written untouched and still bounds the next, and a SNV,
        // which the alignment returns before it reads a base.
        final StringBuilder shiftable = new StringBuilder();
        shiftable.append("##fileformat=VCFv4.2\n")
                .append("##contig=<ID=chr1,length=100000>\n")
                .append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                .append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        shiftable.append(shiftableRecord(2005, referenceBases(2005, 5), referenceBases(2005, 1)));
        shiftable.append(shiftableRecord(2020, referenceBases(2020, 1),
                referenceBases(2020, 1) + referenceBases(2021, 4)));
        shiftable.append(shiftableRecord(2024, referenceBases(2024, 5), referenceBases(2024, 1)));
        shiftable.append(shiftableRecord(5000, referenceBases(5000, 301), referenceBases(5000, 1)));
        shiftable.append(shiftableRecord(6000, referenceBases(6000, 1),
                referenceBases(6000, 1).equals("A") ? "C" : "A"));
        final Path shifts = dir.resolve("shiftable.vcf");
        Files.writeString(shifts, shiftable.toString(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        shifts, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("shiftable.vcf.idx"));
        // The pair `Concordance` walks: a truth callset and an evaluation of it, arranged so that
        // every one of the five concordance states happens once.
        //
        //   chr1:1001  called and agreeing                       true positive
        //   chr1:2001  called with another alternate             false positive AND false negative
        //   chr1:3001  called at a truth locus and FILTERED      filtered false negative
        //   chr1:4001  in truth and not called at all            false negative
        //   chr1:5001  called nowhere near truth and FILTERED    filtered true negative
        //   chr1:6001  called nowhere near truth, unfiltered     false positive
        //
        // The filtered true negative carries TWO filters, so neither of them is unique to it: the
        // filter-analysis table counts uniqueness per RECORD, not per filter.
        final String truthBody =
                concordanceRecord(1001, "A", "C", ".")
                        + concordanceRecord(2001, "A", "ACGT", ".")
                        + concordanceRecord(3001, "A", "C", ".")
                        + concordanceRecord(4001, "A", "C", ".");
        final String evalBody =
                concordanceRecord(1001, "A", "C", "PASS")
                        + concordanceRecord(2001, "A", "AG", "PASS")
                        + concordanceRecord(3001, "A", "C", "LOW_QUAL")
                        + concordanceRecord(5001, "A", "C", "ARTIFACT;LOW_QUAL")
                        + concordanceRecord(6001, "A", "C", "PASS");
        final String vcfHeader = "##fileformat=VCFv4.2\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##FILTER=<ID=LOW_QUAL,Description=\"Low quality\">\n"
                + "##FILTER=<ID=ARTIFACT,Description=\"Artifact\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n";
        for (final String[] pair : new String[][] {{"truth.vcf", truthBody}, {"calls.vcf", evalBody}}) {
            final Path path = dir.resolve(pair[0]);
            Files.writeString(path, vcfHeader + pair[1], StandardCharsets.UTF_8);
            htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                            path, new htsjdk.variant.vcf.VCFCodec(),
                            htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                    .write(dir.resolve(pair[0] + ".idx"));
        }
        // The discovery callset `ValidateBasicSomaticShortMutations` validates, written against
        // `tumor_normal.bam`: the calls are at the two sites that BAM carries an alternate at, and
        // the genotype carries the AD the validator needs. A call with no AD is SKIPPED, which is a
        // judgment of its own and the only one a callset without depths can produce.
        Files.writeString(dir.resolve("somatic.vcf"),
                "##fileformat=VCFv4.2\n"
                        + "##contig=<ID=chr1,length=100000>\n"
                        + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                        + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n"
                        + "chr1\t101\t.\tA\tC\t50\tPASS\t.\tGT:AD\t0/1:4,2\n"
                        + "chr1\t121\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:5,1\n"
                        + "chr1\t141\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n",
                StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(
                        dir.resolve("somatic.vcf"), new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(dir.resolve("somatic.vcf.idx"));
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
        // The same bases as `other.fasta` under a different contig name. `CompareReferences` keys
        // its table by the sequence's MD5 and not by its name, so this reference and that one land
        // on ONE row carrying two names, which is the DIFFER_IN_SEQUENCE_NAMES answer. Without it
        // every accepted row of that tool's array produced the same table: the pair
        // reference/other has nothing in common, and a reference compared with itself is refused.
        try (final htsjdk.samtools.reference.FastaReferenceWriter renamed =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(dir.resolve("renamed.fasta"))
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            final StringBuilder bases = new StringBuilder();
            for (int i = 0; i < 1000; i++) {
                bases.append("ACGT".charAt(i % 4));
            }
            renamed.startSequence("chrRenamed").appendBases(bases.toString());
        }
        // The reads again, under `reference.fasta`'s own dictionary: the one input in this corpus
        // whose `@SQ` line carries an M5, which is the only way a command line reaches
        // `CheckReferenceCompatibility`'s MD5 path.
        md5Bam(dir.resolve("md5header.bam"), dir.resolve("reference.fasta"));

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
        // The same population VCF with its frequencies moved, indexed like the rest: a feature
        // input is queried by interval, so an unindexed one is refused before the traversal.
        final Path shifted = dir.resolve("shifted.vcf");
        Files.writeString(shifted, shiftedPopulationVcf(), StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", shifted.toString()});

        // The bins as an interval LIST, which is the only way one `-L` value can name more than
        // one of them. `FilterIntervals` intersects the requested intervals with its inputs' by
        // list equality and then removes a contig's only survivor, so a window naming a single bin
        // always ends with nothing: the file names four, and a second file names two.
        Files.writeString(dir.resolve("bins.interval_list"),
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"
                        + "chr1\t1\t1000\t+\t.\n"
                        + "chr1\t2001\t3000\t+\t.\n"
                        + "chr1\t4001\t5000\t+\t.\n"
                        + "chr1\t6001\t7000\t+\t.\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("bins2.interval_list"),
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"
                        + "chr1\t1\t1000\t+\t.\n"
                        + "chr1\t2001\t3000\t+\t.\n",
                StandardCharsets.UTF_8);

        // The two files `FilterIntervals` reads, produced by the REFERENCE's own tools: the
        // annotated intervals `AnnotateIntervals` writes and the counts `CollectReadCounts` writes.
        // Both need the copy-number interval rule, which is why they carry it here: those tools
        // refuse anything but OVERLAPPING_ONLY.
        new org.broadinstitute.hellbender.tools.copynumber.AnnotateIntervals()
                .instanceMain(new String[] {
                        "--reference", dir.resolve("reference.fasta").toString(),
                        // FOUR windows rather than one: `FilterIntervals` removes a contig's only
                        // surviving interval, so a table of one row filters to none and the run is
                        // then refused for having nothing left.
                        "--intervals", "chr1:1-1000",
                        "--intervals", "chr1:2001-3000",
                        "--intervals", "chr1:4001-5000",
                        "--intervals", "chr1:6001-7000",
                        "--interval-merging-rule", "OVERLAPPING_ONLY",
                        "--output", dir.resolve("annotated.tsv").toString(),
                });
        new org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts()
                .instanceMain(new String[] {
                        "--input", dir.resolve("reads.bam").toString(),
                        "--reference", dir.resolve("reference.fasta").toString(),
                        "--intervals", "chr1:1-1000",
                        "--intervals", "chr1:2001-3000",
                        "--intervals", "chr1:4001-5000",
                        "--intervals", "chr1:6001-7000",
                        "--interval-merging-rule", "OVERLAPPING_ONLY",
                        "--format", "TSV",
                        "--output", dir.resolve("counts.tsv").toString(),
                });
        // The same table under the name the SV codec recognises. `SimpleCountCodec.canDecode` tests
        // for the extension `.counts.tsv`, and a file called exactly `counts.tsv` does not have it:
        // `PrintReadCounts` refuses it for having no suitable codec, which is a name away from the
        // file it was written to read.
        Files.copy(dir.resolve("counts.tsv"), dir.resolve("sv.counts.tsv"),
                java.nio.file.StandardCopyOption.REPLACE_EXISTING);

        // Two depth-evidence files for `CondenseDepthEvidence`, written in the codec's own layout:
        // a header of column names and zero-based half-open bins. The first is a run of ten
        // adjacent hundred-base bins, a one-base gap, two more, and a contig change, so every
        // maximum and minimum in the array cuts it somewhere different. The second has three
        // samples, fifty-base bins, and a count above 2^31, which `Integer.parseUnsignedInt` reads
        // and the merge then sums as a wrapped int.
        final StringBuilder depth = new StringBuilder("#Chr\tStart\tEnd\tsA\tsB\n");
        for (int i = 0; i < 10; i++) {
            depth.append(String.format("chr1\t%d\t%d\t%d\t%d%n", i * 100, (i + 1) * 100, i + 1, 100 - i));
        }
        depth.append("chr1\t1001\t1101\t11\t90\n");
        depth.append("chr1\t1101\t1201\t12\t89\n");
        depth.append("chr2\t1201\t1301\t13\t88\n");
        Files.writeString(dir.resolve("depth.rd.txt"), depth.toString(), StandardCharsets.UTF_8);
        final StringBuilder depth2 = new StringBuilder("#Chr\tStart\tEnd\tzulu\talpha\tmike\n");
        for (int i = 0; i < 16; i++) {
            depth2.append(String.format("chr1\t%d\t%d\t%d\t%d\t%s%n",
                    i * 50, (i + 1) * 50, i, 2 * i, i == 3 ? "3000000000" : Integer.toString(7 * i)));
        }
        depth2.append("chr2\t0\t50\t1\t2\t3\n");
        Files.writeString(dir.resolve("depth2.rd.txt"), depth2.toString(), StandardCharsets.UTF_8);
        // For `PrintSVEvidence`: a third file naming one more sample at three of depth2's bins, so
        // a row that merges the two widens those bins rather than interleaving them; the pair is
        // named by a `.list`, which Barclay expands for a collection argument. A second `.list`
        // names samples for `--sample-names`, one of them twice and one that no file carries. And
        // a dictionary naming both contigs the depth files use, since the corpus's own name one.
        Files.writeString(dir.resolve("depth3.rd.txt"),
                "#Chr\tStart\tEnd\tbravo\n"
                        + "chr1\t0\t50\t101\n"
                        + "chr1\t100\t150\t102\n"
                        + "chr2\t0\t50\t103\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("evidence.list"),
                // The paths the ROWS see: the corpus is written here and read under /work/fixtures.
                "/work/fixtures/depth2.rd.txt\n/work/fixtures/depth3.rd.txt\n",
                StandardCharsets.UTF_8);
        // The same pair with depth2 named twice, which `FeatureManager` keeps as one input.
        Files.writeString(dir.resolve("evidence_twice.list"),
                "/work/fixtures/depth2.rd.txt\n/work/fixtures/depth3.rd.txt\n"
                        + "/work/fixtures/depth2.rd.txt\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("samples.list"), "zulu\nbravo\nnobody\nzulu\n",
                StandardCharsets.UTF_8);
        // For `SiteDepthtoBAF`: allele depths for two samples at four chr1 sites and one chr2 site,
        // zero-based on disk, and two sites VCFs over the same loci, the second with every ref and
        // alt swapped, so the same depths give each sample the other fraction. Each VCF also holds
        // an indel between two sites, which `BAFSiteIterator` skips, and declares the two contigs
        // `sv.dict` does, since the tool asserts the VCF's dictionary is the walk's. The depths are
        // chosen so each threshold in the array keeps a different set: one site fails the
        // chi-squared test, one has a total under 30, and one has samples far enough apart that a
        // tight --max-std drops the whole locus.
        final String sdA = "chr1\t99\ts1\t10\t12\t0\t0\n"
                + "chr1\t199\ts1\t0\t0\t30\t2\n"
                + "chr1\t299\ts1\t8\t0\t9\t0\n"
                + "chr1\t399\ts1\t0\t40\t0\t35\n"
                + "chr2\t99\ts1\t20\t0\t0\t21\n";
        final String sdB = "chr1\t99\ts2\t14\t9\t0\t0\n"
                + "chr1\t299\ts2\t11\t0\t12\t1\n"
                + "chr1\t399\ts2\t0\t20\t0\t60\n"
                + "chr2\t99\ts2\t30\t0\t0\t25\n";
        Files.writeString(dir.resolve("depth.sd.txt"), sdA, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("depth2.sd.txt"), sdB, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sd.list"),
                "/work/fixtures/depth.sd.txt\n/work/fixtures/depth2.sd.txt\n",
                StandardCharsets.UTF_8);
        final String[][] snps = {
                {"chr1", "100", "A", "C"}, {"chr1", "150", "AC", "A"}, {"chr1", "200", "G", "T"},
                {"chr1", "300", "A", "G"}, {"chr1", "400", "C", "T"}, {"chr2", "100", "A", "T"}};
        for (final boolean swapped : new boolean[] {false, true}) {
            final StringBuilder vcf = new StringBuilder("##fileformat=VCFv4.2\n"
                    + "##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
            for (final String[] snp : snps) {
                final boolean swap = swapped && snp[2].length() == 1;
                vcf.append(snp[0]).append('\t').append(snp[1]).append("\t.\t")
                        .append(swap ? snp[3] : snp[2]).append('\t')
                        .append(swap ? snp[2] : snp[3]).append("\t.\t.\t.\n");
            }
            Files.writeString(dir.resolve(swapped ? "baf_sites2.vcf" : "baf_sites.vcf"),
                    vcf.toString(), StandardCharsets.UTF_8);
        }
        Files.writeString(dir.resolve("sv.dict"),
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n",
                StandardCharsets.UTF_8);

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
        // The three Mutect gathers. `GatherPileupSummaries` gets the REFERENCE's own summaries of
        // one sample over three windows: two holding sites and one holding none, named by a
        // `.list` out of order, so the gather has to sort them by their first record and drop the
        // empty one. The population's sites sit below 5000, so the windows split there. A second
        // `.list` adds the first window again under another sample's name, which the gather
        // refuses: the corpus's two BAMs carry the same sample, so the name is rewritten here.
        for (final String[] window : new String[][] {
                {"chr1:2001-60000", "summaries_b.table"}, {"chr1:1-2000", "summaries_a.table"},
                {"chr1:90001-100000", "summaries_empty.table"}}) {
            new org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries()
                    .instanceMain(new String[] {
                            "--input", dir.resolve("reads.bam").toString(),
                            "--variant", population.toString(),
                            "--intervals", window[0],
                            "--output", dir.resolve(window[1]).toString(),
                    });
        }
        Files.writeString(dir.resolve("pileups.list"),
                "/work/fixtures/summaries_b.table\n/work/fixtures/summaries_empty.table\n"
                        + "/work/fixtures/summaries_a.table\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("summaries_other.table"),
                Files.readString(dir.resolve("summaries_a.table"))
                        .replace("SAMPLE=sample1", "SAMPLE=other"),
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("mixed_pileups.list"),
                "/work/fixtures/summaries_b.table\n/work/fixtures/summaries_other.table\n",
                StandardCharsets.UTF_8);
        // `GatherNormalArtifactData` gets the reference's own tables for each sample of
        // `tumor_normal.bam` taken as the normal, so the two shards hold different records.
        for (final String sample : new String[] {"normal", "tumor"}) {
            new org.broadinstitute.hellbender.tools.walkers.mutect.GetNormalArtifactData()
                    .instanceMain(new String[] {
                            "--input", dir.resolve("tumor_normal.bam").toString(),
                            "--reference", dir.resolve("reference.fasta").toString(),
                            "--normal-sample", sample,
                            "--output", dir.resolve("artifact_" + sample + ".table").toString(),
                    });
        }
        // `tumor` as the normal finds nothing, so its table is a header alone. A third table in the
        // writer's own layout (ints, a double as `Double.toString` writes it, the type's name)
        // gives the gather two shards with records, so their order shows in the output.
        Files.writeString(dir.resolve("artifact_extra.table"),
                "normal_alt\tnormal_dp\ttumor_alt\ttumor_dp\tdownsampling\ttype\n"
                        + "0\t10\t3\t12\t0.5\tSNV\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("artifacts.list"),
                "/work/fixtures/artifact_extra.table\n/work/fixtures/artifact_tumor.table\n"
                        + "/work/fixtures/artifact_normal.table\n",
                StandardCharsets.UTF_8);
        // `MergeMutectStats` reads the two-column table Mutect2 writes. The `.list` names one
        // shard twice, which the tool's `LinkedHashSet` reads once, and the second `.list` adds a
        // shard carrying a statistic the aggregation map does not hold.
        Files.writeString(dir.resolve("a.stats"), "statistic\tvalue\ncallable\t1000.0\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("b.stats"), "statistic\tvalue\ncallable\t2.5E7\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("odd.stats"),
                "statistic\tvalue\ncallable\t3.0\nrejected\t1.0\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("stats.list"),
                "/work/fixtures/a.stats\n/work/fixtures/b.stats\n/work/fixtures/a.stats\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("odd_stats.list"),
                "/work/fixtures/a.stats\n/work/fixtures/odd.stats\n", StandardCharsets.UTF_8);
        // Segment files for the copy-number utilities, which the codec reads only under a `.seg`,
        // `.maf` or `.maf.annotated` name: `regions.tsv` holds the same rows and is refused by
        // its name alone. `regions.seg` is unsorted, overlaps itself in a chain, abuts, and has a
        // chr2 row that a chr1-only dictionary sorts last; `regions2.seg` uses another spelling of
        // the locatable columns and carries comments. `tumour.seg` and `normal.seg` are called
        // segments: the normal amplification shares both breakpoints of one tumour segment within
        // a few bases and reciprocally overlaps another, and its deletion matches nothing.
        final String regions = "CONTIG\tSTART\tEND\tname\tvalue\tCALL\n"
                + "chr1\t300\t400\tc\t1.5\t-\n"
                + "chr1\t1\t100\ta\t0.5\t+\n"
                + "chr1\t50\t150\tb\t0.5\t+\n"
                + "chr1\t120\t200\tb\t0.7\t+\n"
                + "chr1\t401\t500\td\t1.5\t-\n"
                + "chr2\t1\t100\te\t2.0\t0\n";
        Files.writeString(dir.resolve("regions.seg"), regions, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("regions.tsv"), regions, StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("regions2.seg"),
                "#a note\nChromosome\tStart_Position\tEnd_Position\tname\tCALL\n"
                        + "chr1\t10\t60\tx\t+\n"
                        + "chr1\t61\t90\ty\t+\n"
                        + "chr1\t500\t900\tz\t-\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("annotations.list"), "CALL\nname\n", StandardCharsets.UTF_8);
        // A reference naming both contigs the segment files use, for the tools that require one:
        // under `reference.fasta` alone a file with a chr2 row is refused by the dictionary sort.
        try (final htsjdk.samtools.reference.FastaReferenceWriter both =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(dir.resolve("two_contigs.fasta"))
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            final StringBuilder bases = new StringBuilder();
            for (int i = 0; i < 30000; i++) {
                bases.append("ACGT".charAt(i % 4));
            }
            both.startSequence("chr1").appendBases(bases.toString());
            both.startSequence("chr2").appendBases(bases.toString());
        }
        Files.writeString(dir.resolve("tumour.seg"),
                "CONTIG\tSTART\tEND\tCALL\tMEAN\n"
                        + "chr1\t1\t1000\t+\t0.9\n"
                        + "chr1\t1001\t5000\t0\t0.0\n"
                        + "chr1\t5001\t8000\t+\t0.8\n"
                        + "chr1\t8001\t9000\t-\t-0.9\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("normal.seg"),
                "CONTIG\tSTART\tEND\tCALL\n"
                        + "chr1\t5\t995\t+\n"
                        + "chr1\t5200\t8100\t+\n"
                        + "chr1\t20000\t21000\t-\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("normal2.seg"),
                "CONTIG\tSTART\tEND\tCALL\n"
                        + "chr1\t8001\t9000\t-\n",
                StandardCharsets.UTF_8);
        // Tranche shards for `GatherTranches`, in the VQSLOD layout of version 6: two shards at
        // the same four levels, a third with a level the others lack, and one of version 5, which
        // the reader refuses by its version line. A `.list` names the first three, and another
        // names the requested sensitivities, out of order, since the tool sorts them in place.
        final String trancheHeader = "# Variant quality score tranches file\n# Version number 6\n"
                + "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,"
                + "accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";
        Files.writeString(dir.resolve("shard1.tranches"), trancheHeader
                + "4.0000,100,20,2.0000,1.5000,4.0000,VQSRTranche,SNP,1000,500,0.5000\n"
                + "2.0000,200,50,2.1000,1.6000,2.0000,VQSRTranche,SNP,1000,800,0.8000\n"
                + "0.0000,300,90,2.2000,1.7000,0.0000,VQSRTranche,SNP,1000,950,0.9500\n"
                + "-2.0000,400,150,2.3000,1.8000,-2.0000,VQSRTranche,SNP,1000,990,0.9900\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("shard2.tranches"), trancheHeader
                + "4.0000,60,10,1.8000,1.4000,4.0000,VQSRTranche,SNP,1000,450,0.4500\n"
                + "2.0000,130,30,1.9000,1.5000,2.0000,VQSRTranche,SNP,1000,780,0.7800\n"
                + "0.0000,220,70,2.0000,1.6000,0.0000,VQSRTranche,SNP,1000,940,0.9400\n"
                + "-2.0000,330,120,2.1000,1.7000,-2.0000,VQSRTranche,SNP,1000,985,0.9850\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("shard3.tranches"), trancheHeader
                + "4.0000,10,2,2.0000,1.5000,4.0000,VQSRTranche,SNP,1000,400,0.4000\n"
                + "1.0000,90,25,2.0000,1.5000,1.0000,VQSRTranche,SNP,1000,700,0.7000\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("old.tranches"),
                trancheHeader.replace("Version number 6", "Version number 5")
                        + "4.0000,10,2,2.0000,1.5000,4.0000,VQSRTranche,SNP,1000,400,0.4000\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("tranches.list"),
                "/work/fixtures/shard1.tranches\n/work/fixtures/shard2.tranches\n"
                        + "/work/fixtures/shard3.tranches\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("levels.list"), "90.0\n99.0\n95.0\n", StandardCharsets.UTF_8);
        // Mutect2-shaped calls for the mitochondrial filters: `AD` per allele, `AF` per alternate,
        // and the `AS_FilterStatus` a filtered Mutect2 VCF carries, one entry per alternate. Five
        // unfiltered sites sit below the low-heteroplasmy fraction of 0.1, which is past the
        // default allowance of three; one filtered low site does not count; the multi-allelic
        // record has one deep and one shallow alternate, so an allele filter splits it. The
        // second file is the same calls without `AS_FilterStatus`, which the filters' merge
        // refuses the moment anything is filtered.
        final String mtHeader = "##fileformat=VCFv4.2\n"
                + "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
                + "##FILTER=<ID=weak_evidence,Description=\"Mutation does not meet likelihood threshold\">\n"
                + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                + "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=\"Filter status for each allele\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n";
        final String[][] calls = {
                {"100", "A", "C", "PASS", "SITE", "0/1", "190,10", "0.05"},
                {"200", "C", "G", "PASS", "SITE", "0/1", "180,9", "0.048"},
                {"300", "G", "T", "PASS", "SITE", "0/1", "150,40", "0.21"},
                {"400", "T", "A", "PASS", "SITE", "0/1", "200,8", "0.038"},
                {"500", "A", "G", "weak_evidence", "weak_evidence", "0/1", "300,5", "0.016"},
                {"600", "C", "T", "PASS", "SITE", "0/1", "100,6", "0.057"},
                {"700", "G", "A", "PASS", "SITE", "0/1", "120,7", "0.055"},
                {"800", "T", "C,G", "PASS", "SITE|SITE", "0/1/2", "90,45,4", "0.33,0.03"}};
        final StringBuilder mt = new StringBuilder(mtHeader);
        final StringBuilder mtPlain = new StringBuilder(mtHeader);
        for (final String[] call : calls) {
            final String prefix = "chr1\t" + call[0] + "\t.\t" + call[1] + "\t" + call[2] + "\t.\t" + call[3] + "\t";
            final String sample = "\tGT:AD:AF\t" + call[5] + ":" + call[6] + ":" + call[7] + "\n";
            mt.append(prefix).append("AS_FilterStatus=").append(call[4]).append(sample);
            mtPlain.append(prefix).append(".").append(sample);
        }
        Files.writeString(dir.resolve("mito.vcf"), mt.toString(), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("mito_plain.vcf"), mtPlain.toString(), StandardCharsets.UTF_8);
        // For `CalculateAverageCombinedAnnotations`: sites as GenomicsDB leaves them, with
        // `RAW_GT_COUNT` (hom-ref, het, hom-var) and two summed annotations. One site's divisor is
        // zero, so it is written through untouched; one carries only one of the sums. The second
        // file drops `RAW_GT_COUNT` from its third site, which the tool refuses where it reaches it.
        final String sumsHeader = "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=MQ_SUM,Number=1,Type=Float,Description=\"Summed mapping quality\">\n"
                + "##INFO=<ID=QD_SUM,Number=1,Type=Float,Description=\"Summed quality by depth\">\n"
                + "##INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=\"Genotype counts\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        final String[] sumsRows = {
                "chr1\t100\t.\tA\tC\t50\tPASS\tMQ_SUM=180.0;QD_SUM=33.3;RAW_GT_COUNT=4,2,1\n",
                "chr1\t200\t.\tC\tG\t50\tPASS\tMQ_SUM=60.0;QD_SUM=10.0;RAW_GT_COUNT=7,0,0\n",
                "chr1\t300\t.\tG\tT\t50\tPASS\tQD_SUM=7.25;RAW_GT_COUNT=1,3,0\n",
                "chr1\t400\t.\tT\tA\t50\tPASS\tMQ_SUM=240.5;QD_SUM=41.0;RAW_GT_COUNT=0,1,3\n"};
        Files.writeString(dir.resolve("summed.vcf"), sumsHeader + String.join("", sumsRows),
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("summed_missing.vcf"), sumsHeader + sumsRows[0] + sumsRows[1]
                        + sumsRows[2].replace(";RAW_GT_COUNT=1,3,0", "") + sumsRows[3],
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sums.list"), "QD_SUM\nMQ_SUM\nABSENT_SUM\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sums_two.list"), "MQ_SUM\nQD_SUM\n", StandardCharsets.UTF_8);
        // For `FilterVariantTranches`: calls scored in `CNN_2D`, one unscored, one already filtered,
        // two indels, and resources that share their sites. The indel resource spells the deletion
        // at 600 with one more reference base, so the match goes through `isAlleleInList`'s
        // extension of the shorter reference; the second resource file carries the indels only.
        final String cnnHeader = "##fileformat=VCFv4.2\n"
                + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
                + "##INFO=<ID=CNN_2D,Number=1,Type=Float,Description=\"2D CNN score\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        Files.writeString(dir.resolve("cnn.vcf"), cnnHeader
                + "chr1\t100\t.\tA\tC\t50\t.\tCNN_2D=2.5\n"
                + "chr1\t200\t.\tC\tT\t50\t.\tCNN_2D=-1.0\n"
                + "chr1\t300\t.\tG\tA\t50\tPASS\tCNN_2D=0.5\n"
                + "chr1\t400\t.\tT\tG\t50\t.\tCNN_2D=1.8\n"
                + "chr1\t500\t.\tA\tAT\t50\t.\tCNN_2D=-0.5\n"
                + "chr1\t600\t.\tCTG\tC\t50\t.\tCNN_2D=1.2\n"
                + "chr1\t700\t.\tG\tC\t50\t.\t.\n"
                + "chr1\t800\t.\tT\tC\t50\tLowQual\tCNN_2D=-2.0\n",
                StandardCharsets.UTF_8);
        final String resourceHeader = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        Files.writeString(dir.resolve("snp_truth.vcf"), resourceHeader
                + "chr1\t100\t.\tA\tC\t.\t.\t.\n"
                + "chr1\t200\t.\tC\tT\t.\t.\t.\n"
                + "chr1\t300\t.\tG\tA,T\t.\t.\t.\n"
                + "chr1\t400\t.\tT\tG\t.\t.\t.\n"
                + "chr1\t500\t.\tA\tAT\t.\t.\t.\n"
                + "chr1\t800\t.\tT\tA\t.\t.\t.\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("indel_truth.vcf"), resourceHeader
                + "chr1\t500\t.\tA\tAT\t.\t.\t.\n"
                + "chr1\t600\t.\tCTGG\tCG\t.\t.\t.\n",
                StandardCharsets.UTF_8);
        // A resource is queried by interval, which needs an index.
        for (final String truth : new String[] {"snp_truth.vcf", "indel_truth.vcf"}) {
            new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                    .instanceMain(new String[] {"-I", dir.resolve(truth).toString()});
        }
        Files.writeString(dir.resolve("truth.list"),
                "/work/fixtures/snp_truth.vcf\n/work/fixtures/indel_truth.vcf\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("snp_tranches.list"), "99.0\n90.0\n99.0\n",
                StandardCharsets.UTF_8);
        vqsrFixtures(dir);
        svStratifyFixtures(dir);
        svConcordanceFixtures(dir);
        svClusterFixtures(dir);
        svAnnotateFixtures(dir);
        referenceBlockFixtures(dir);
        // The sites `ASEReadCounter` reads, with one sample's genotypes, indexed because the
        // walker queries them by locus; and the same file without its index, which it refuses.
        final String aseSites = "##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
                + "chr1\t1001\trs1\tA\tG\t.\t.\t.\tGT\t0/1\n"
                + "chr1\t1102\trs2\tC\tT\t.\t.\t.\tGT\t0/1\n"
                + "chr1\t1203\trs3\tG\tA,T\t.\t.\t.\tGT\t1/2\n"
                + "chr1\t1304\t.\tT\tC\t.\t.\t.\tGT\t1/1\n";
        Files.writeString(dir.resolve("ase_sites.vcf"), aseSites, StandardCharsets.UTF_8);
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", dir.resolve("ase_sites.vcf").toString()});
        Files.writeString(dir.resolve("ase_sites_unindexed.vcf"), aseSites, StandardCharsets.UTF_8);
        // A SAM TEXT donor for ConvertHeaderlessHadoopBamShardToBam, whose header is read by
        // SamReaderFactory rather than the engine. The version is 1.5 because that tool keeps it,
        // and the record after the header is where SAMTextReader stops reading it.
        Files.writeString(dir.resolve("donor.sam"),
                "@HD\tVN:1.5\tSO:queryname\n"
                + "@SQ\tSN:chr1\tLN:100000\n"
                + "@SQ\tSN:chr2\tLN:50000\n"
                + "@RG\tID:donor\tSM:donor\tPL:ILLUMINA\n"
                + "@CO\tdonated header\n"
                + "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t####\n"
                + "@CO\tafter the first record, so not header\n",
                StandardCharsets.UTF_8);
        pathSeqTaxonomy(dir);
        System.out.println("wrote " + dir);
    }

    /**
     * What `ReferenceBlockConcordance` compares: two single-sample GVCFs whose reference blocks cut
     * the same stretch differently and at different confidences, each with a variant site the
     * hom-ref filter drops, and a two-sample GVCF the length extraction refuses.
     */
    static void referenceBlockFixtures(final Path dir) throws Exception {
        final String header = "##fileformat=VCFv4.2\n"
                + "##ALT=<ID=NON_REF,Description=\"Any other allele\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                + "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum depth\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##contig=<ID=chr2,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        Files.writeString(dir.resolve("rbc_truth.g.vcf"), header + "\tS1\n"
                + "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:GQ:MIN_DP\t0/0:40:10\n"
                + "chr1\t101\t.\tA\tC,<NON_REF>\t50\t.\t.\tGT:GQ\t0/1:50\n"
                + "chr1\t102\t.\tC\t<NON_REF>\t.\t.\tEND=300\tGT:GQ:MIN_DP\t0/0:60:20\n"
                + "chr1\t301\t.\tG\t<NON_REF>\t.\t.\tEND=301\tGT:GQ:MIN_DP\t0/0:99:30\n"
                + "chr2\t10\t.\tT\t<NON_REF>\t.\t.\tEND=50\tGT:GQ:MIN_DP\t0/0:20:5\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("rbc_eval.g.vcf"), header + "\tS1\n"
                + "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=50\tGT:GQ:MIN_DP\t0/0:30:10\n"
                + "chr1\t51\t.\tT\t<NON_REF>\t.\t.\tEND=200\tGT:GQ:MIN_DP\t0/0:45:12\n"
                + "chr1\t201\t.\tA\tG,<NON_REF>\t60\t.\t.\tGT:GQ\t1/1:60\n"
                + "chr1\t202\t.\tC\t<NON_REF>\t.\t.\tEND=301\tGT:GQ:MIN_DP\t0/0:99:25\n"
                + "chr2\t20\t.\tT\t<NON_REF>\t.\t.\tEND=40\tGT:GQ:MIN_DP\t0/0:20:5\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("rbc_multi.g.vcf"), header + "\tS1\tS2\n"
                + "chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:GQ\t0/0:40\t0/0:30\n",
                StandardCharsets.UTF_8);
    }

    /**
     * What `SVAnnotate` reads: SVs of every type the rules tell apart, a protein-coding GTF and a
     * non-coding BED.
     *
     * `sv_genes.gtf` holds `alpha` on the plus strand, two exons with CDS, start and stop codons and
     * UTRs, and `beta` on the minus strand, one exon and a CDS, in the GENCODE form the codec
     * accepts. `sv_noncoding.bed` has an enhancer and a silencer. `sv_annotate.vcf` places a
     * deletion in a CDS, a duplication in a UTR, a deletion in the enhancer far from both genes, a
     * duplication spanning `beta`, an insertion in an intron, an inversion spanning `alpha`, a
     * breakend pair short enough for `--max-breakend-as-cnv-length`, a multiallelic CNV, a complex
     * `dupINV` and a translocation.
     */
    static void svAnnotateFixtures(final Path dir) throws Exception {
        Files.writeString(dir.resolve("sv_genes.gtf"), "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tgene_id \"GA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; level 2; havana_gene \"OTTHUMG1.1\";\n" + "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; tag \"basic\"; havana_gene \"OTTHUMG1.1\";\n" + "chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 1; exon_id \"TA.1.1\";\n" + "chr1\tHAVANA\texon\t2000\t3000\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 2; exon_id \"TA.1.2\";\n" + "chr1\tHAVANA\tCDS\t1100\t1200\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 1; exon_id \"TA.1.1\";\n" + "chr1\tHAVANA\tCDS\t2000\t2800\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 2; exon_id \"TA.1.2\";\n" + "chr1\tHAVANA\tstart_codon\t1100\t1102\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 1; exon_id \"TA.1.1\";\n" + "chr1\tHAVANA\tstop_codon\t2798\t2800\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 2; exon_id \"TA.1.2\";\n" + "chr1\tHAVANA\tUTR\t1000\t1099\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 1; exon_id \"TA.1.1\";\n" + "chr1\tHAVANA\tUTR\t2801\t3000\t.\t+\t.\tgene_id \"GA.1\"; transcript_id \"TA.1\"; gene_type \"protein_coding\"; gene_name \"alpha\"; transcript_type \"protein_coding\"; transcript_name \"TA.1\"; level 2; exon_number 2; exon_id \"TA.1.2\";\n" + "chr1\tHAVANA\tgene\t10000\t12000\t.\t-\t.\tgene_id \"GB.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; level 2; havana_gene \"OTTHUMG2.1\";\n" + "chr1\tHAVANA\ttranscript\t10000\t12000\t.\t-\t.\tgene_id \"GB.1\"; transcript_id \"TB.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TB.1\"; level 2; tag \"basic\"; havana_gene \"OTTHUMG2.1\";\n" + "chr1\tHAVANA\texon\t10000\t12000\t.\t-\t.\tgene_id \"GB.1\"; transcript_id \"TB.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TB.1\"; level 2; exon_number 1; exon_id \"TB.1.1\";\n" + "chr1\tHAVANA\tCDS\t10200\t11800\t.\t-\t.\tgene_id \"GB.1\"; transcript_id \"TB.1\"; gene_type \"protein_coding\"; gene_name \"beta\"; transcript_type \"protein_coding\"; transcript_name \"TB.1\"; level 2; exon_number 1; exon_id \"TB.1.1\";\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_noncoding.bed"),
                "chr1\t4000\t4500\tenhancer\t.\t+\nchr1\t20000\t20100\tsilencer\t.\t-\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_annotate.vcf"), "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">\n"
                + "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second end\">\n"
                + "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">\n"
                + "##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description=\"Complex type\">\n"
                + "##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description=\"Complex intervals\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##contig=<ID=chr2,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "chr1\t1150\tdel_cds\tN\t<DEL>\t.\tPASS\tEND=1180;SVTYPE=DEL;SVLEN=31\n"
                + "chr1\t1300\tins_intron\tN\t<INS:ME>\t.\tPASS\tEND=1301;SVTYPE=INS;SVLEN=300\n"
                + "chr1\t2850\tdup_utr\tN\t<DUP>\t.\tPASS\tEND=2900;SVTYPE=DUP;SVLEN=51\n"
                + "chr1\t4100\tdel_enh\tN\t<DEL>\t.\tPASS\tEND=4200;SVTYPE=DEL;SVLEN=101\n"
                + "chr1\t500\tinv_span\tN\t<INV>\t.\tPASS\tEND=3500;SVTYPE=INV;SVLEN=3001\n"
                + "chr1\t9000\tdup_span\tN\t<DUP>\t.\tPASS\tEND=13000;SVTYPE=DUP;SVLEN=4001\n"
                + "chr1\t10500\tcnv\tN\t<CNV>\t.\tPASS\tEND=10800;SVTYPE=CNV;SVLEN=301\n"
                + "chr1\t11500\tbnd_short\tN\tN[chr1:11900[\t.\tPASS\tSVTYPE=BND;CHR2=chr1;END2=11900;STRANDS=+-;SVLEN=400\n"
                + "chr1\t15000\tcpx\tN\t<CPX>\t.\tPASS\tEND=16000;SVTYPE=CPX;SVLEN=1000;CPX_TYPE=dupINV;CPX_INTERVALS=DUP_chr1:2500-2700,INV_chr1:2700-16000\n"
                + "chr1\t20050\tctx\tN\t<CTX>\t.\tPASS\tEND=20050;SVTYPE=CTX;CHR2=chr2;END2=500\n",
                StandardCharsets.UTF_8);
    }

    /**
     * What `SVCluster` reads: a reference the SV calls fit on, calls that group, and ploidies.
     *
     * `sv_reference.fasta` names chr1 and chr2 at 100 kb, the lengths the SV fixtures' headers
     * declare. `sv_cluster_input.vcf` holds three overlapping deletions (two with split-read or
     * paired-end evidence and one depth-only, which the breakpoints ignore), two duplications, two
     * insertions of different mobile-element subtypes, an inversion, two breakends to the same far
     * end, and a depth-only deletion and duplication that only `--enable-cnv` groups into a CNV.
     * Every genotype carries `ECN`, which carrier counting requires. `sv_ploidy.tsv` gives both
     * samples a ploidy of two; `sv_ploidy_haploid.tsv` gives `s2` one on chr1, which only a sample
     * `--fast-mode` dropped ever reads.
     */
    static void svClusterFixtures(final Path dir) throws Exception {
        try (final htsjdk.samtools.reference.FastaReferenceWriter writer =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(dir.resolve("sv_reference.fasta"))
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            final StringBuilder bases = new StringBuilder();
            for (int i = 0; i < 100000; i++) {
                bases.append("ACGT".charAt(i % 4));
            }
            writer.startSequence("chr1").appendBases(bases.toString());
            writer.startSequence("chr2").appendBases(bases.toString());
        }
        final String header = "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">\n"
                + "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence\">\n"
                + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">\n"
                + "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second end\">\n"
                + "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=ECN,Number=1,Type=Integer,Description=\"Expected copy number\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##contig=<ID=chr2,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";
        Files.writeString(dir.resolve("sv_cluster_input.vcf"), header
                + "chr1\t1000\td1\tN\t<DEL>\t.\tPASS\tEND=2000;SVTYPE=DEL;ALGORITHMS=manta;EVIDENCE=PE,SR\tGT:ECN\t0/1:2\t0/0:2\n"
                + "chr1\t1050\td2\tN\t<DEL>\t.\tPASS\tEND=2050;SVTYPE=DEL;ALGORITHMS=wham;EVIDENCE=PE\tGT:ECN\t0/0:2\t0/1:2\n"
                + "chr1\t1100\td3\tN\t<DEL>\t.\tLowQual\tEND=1990;SVTYPE=DEL;ALGORITHMS=depth;EVIDENCE=RD\tGT:ECN\t0/1:2\t0/1:2\n"
                + "chr1\t5000\tu1\tN\t<DUP>\t.\tPASS\tEND=8000;SVTYPE=DUP;ALGORITHMS=depth\tGT:ECN\t0/1:2\t0/0:2\n"
                + "chr1\t5100\tu2\tN\t<DUP>\t.\tPASS\tEND=8100;SVTYPE=DUP;ALGORITHMS=depth\tGT:ECN\t0/0:2\t0/1:2\n"
                + "chr1\t20000\ti1\tN\t<INS:ME:ALU>\t.\tPASS\tEND=20001;SVTYPE=INS;SVLEN=300;ALGORITHMS=manta;STRANDS=+-\tGT:ECN\t1/1:2\t0/0:2\n"
                + "chr1\t20020\ti2\tN\t<INS:ME:LINE1>\t.\tPASS\tEND=20021;SVTYPE=INS;SVLEN=310;ALGORITHMS=melt;STRANDS=+-\tGT:ECN\t0/0:2\t0/1:2\n"
                + "chr1\t40000\tv1\tN\t<INV>\t.\tPASS\tEND=42000;SVTYPE=INV;ALGORITHMS=manta;STRANDS=++\tGT:ECN\t0/1:2\t0/0:2\n"
                + "chr2\t3000\tb1\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=60000;STRANDS=+-;ALGORITHMS=manta\tGT:ECN\t0/1:2\t0/0:2\n"
                + "chr2\t3010\tb2\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=60010;STRANDS=+-;ALGORITHMS=wham\tGT:ECN\t0/0:2\t0/1:2\n"
                + "chr2\t9000\tc1\tN\t<DEL>\t.\tPASS\tEND=9500;SVTYPE=DEL;ALGORITHMS=depth\tGT:ECN\t0/1:2\t0/0:2\n"
                + "chr2\t9020\tc2\tN\t<DUP>\t.\tPASS\tEND=9510;SVTYPE=DUP;ALGORITHMS=depth\tGT:ECN\t0/0:2\t0/1:2\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_ploidy.tsv"),
                "SAMPLE\tchr1\tchr2\ns1\t2\t2\ns2\t2\t2\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_ploidy_haploid.tsv"),
                "SAMPLE\tchr1\tchr2\ns1\t2\t2\ns2\t1\t2\n", StandardCharsets.UTF_8);
        // What `GroupedSVCluster` adds to them: a track over the deletions and the insertions, the
        // strata to group by, and the thresholds per group, in full and with one group missing,
        // which `onTraversalStart` refuses after the writer exists.
        Files.writeString(dir.resolve("sv_grp_track.bed"),
                "chr1\t899\t2100\nchr1\t19989\t20030\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_grp_strata.tsv"),
                "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\n"
                + "DEL_small\tDEL\t50\t5000\tNA\n"
                + "DEL_large\tDEL\t5000\tNA\tNA\n"
                + "DUP_all\tDUP\tNA\tNA\tNA\n"
                + "INS_tr\tINS\tNA\tNA\tTR\n"
                + "BND_all\tBND\tNA\tNA\tNA\n", StandardCharsets.UTF_8);
        final String clustering = "NAME\tRECIPROCAL_OVERLAP\tSIZE_SIMILARITY\tBREAKEND_WINDOW\tSAMPLE_OVERLAP\n"
                + "DEL_small\t0.5\t0.0\t500\t0.0\n"
                + "DEL_large\t0.8\t0.0\t1000\t0.0\n"
                + "DUP_all\t0.8\t0.0\t1000\t0.0\n"
                + "INS_tr\t0.1\t0.0\t100\t0.0\n";
        Files.writeString(dir.resolve("sv_grp_clustering.tsv"),
                clustering + "BND_all\t0.0\t0.0\t50\t0.0\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_grp_clustering_short.tsv"), clustering,
                StandardCharsets.UTF_8);
    }

    /**
     * What `SVConcordance` compares: a truth and an evaluation SV callset sharing one sample.
     *
     * Each eval record sits near a truth record of its type, except one deletion with nothing near
     * it (a false positive) and a CNV, whose copy numbers are compared instead of its genotypes.
     * `s1` is in both files and `s2`, `s3` in one each, so only `s1` is annotated. The truth
     * insertion carries AC, AF and AN, which are copied rather than recounted. `sv_eval_sites.vcf`
     * is the eval file without samples, which gets no allele counts of its own.
     */
    static void svConcordanceFixtures(final Path dir) throws Exception {
        final String lines = "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">\n"
                + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">\n"
                + "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second end\">\n"
                + "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">\n"
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
                + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##contig=<ID=chr2,length=100000>\n";
        final String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        Files.writeString(dir.resolve("sv_truth.vcf"), lines + columns + "\ts1\ts2\n"
                + "chr1\t1000\tt_del\tN\t<DEL>\t.\tPASS\tEND=1300;SVTYPE=DEL;ALGORITHMS=manta\tGT\t0/1\t1/1\n"
                + "chr1\t5000\tt_dup\tN\t<DUP>\t.\tPASS\tEND=6000;SVTYPE=DUP;ALGORITHMS=depth\tGT\t0/1\t0/0\n"
                + "chr1\t30000\tt_ins\tN\t<INS>\t.\tPASS\tEND=30001;SVTYPE=INS;SVLEN=400;ALGORITHMS=manta;AC=3;AF=0.75;AN=4\tGT\t1/1\t0/1\n"
                + "chr1\t40000\tt_cnv\tN\t<DEL>,<DUP>\t.\tPASS\tEND=45000;SVTYPE=CNV;ALGORITHMS=depth\tGT:CN\t./.:1\t./.:3\n"
                + "chr2\t2000\tt_bnd\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=9000;STRANDS=+-;ALGORITHMS=manta\tGT\t0/1\t0/0\n",
                StandardCharsets.UTF_8);
        final String evalRecords =
                "chr1\t1010\te_del\tN\t<DEL>\t.\tPASS\tEND=1290;SVTYPE=DEL;SVLEN=-280;ALGORITHMS=manta"
                + "\tchr1\t5100\te_dup\tN\t<DUP>\t.\t.\tEND=6100;SVTYPE=DUP;ALGORITHMS=depth"
                + "\tchr1\t20000\te_fp\tN\t<DEL>\t.\tPASS\tEND=21000;SVTYPE=DEL;ALGORITHMS=manta"
                + "\tchr1\t30050\te_ins\tN\t<INS>\t.\tPASS\tEND=30051;SVTYPE=INS;SVLEN=380;ALGORITHMS=manta"
                + "\tchr1\t40100\te_cnv\tN\t<DEL>,<DUP>\t.\tPASS\tEND=45100;SVTYPE=CNV;ALGORITHMS=depth"
                + "\tchr2\t2010\te_bnd\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=9005;STRANDS=+-;ALGORITHMS=manta";
        final String[] records = evalRecords.split("\t(?=chr[12]\t)");
        final String[] genotypes = {"GT\t0/1\t0/0", "GT\t1/1\t./.", "GT\t0/1\t0/1", "GT\t1/1\t0/1",
                "GT:CN\t./.:1\t./.:2", "GT\t0/1\t0/0"};
        final StringBuilder eval = new StringBuilder(lines + columns + "\ts1\ts3\n");
        final StringBuilder sites = new StringBuilder(lines + columns.replace("\tFORMAT", "") + "\n");
        for (int i = 0; i < records.length; i++) {
            eval.append(records[i]).append('\t').append(genotypes[i]).append('\n');
            sites.append(records[i]).append('\n');
        }
        Files.writeString(dir.resolve("sv_eval.vcf"), eval.toString(), StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_eval_sites.vcf"), sites.toString(), StandardCharsets.UTF_8);
    }

    /**
     * What `SVStratify` reads: SV calls of every type the engine tells apart, two BED tracks, and
     * two stratification tables.
     *
     * `sv_calls.vcf` has a small and a large deletion, a duplication, an insertion with a length and
     * one without, an inversion and a breakend to chr2. `sv_bad.vcf` is the same header over a
     * deletion with no `ALGORITHMS`, which `SVCallRecordUtils.create` refuses inside `apply`.
     * `sv_rm.bed` covers the small deletion, the insertion and the breakend's far end; `sv_sd.bed`
     * part of the large deletion and the inversion. `sv_strata.tsv` names only the RM track and its
     * strata are mutually exclusive; `sv_overlap.tsv` names both tracks and has two deletion strata
     * one record can match at once, which only `--allow-multiple-matches` lets through. The two
     * `.list` files expand `--track-name` and `--track-intervals` to both tracks.
     */
    static void svStratifyFixtures(final Path dir) throws Exception {
        final String header = "##fileformat=VCFv4.2\n"
                + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                + "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Algorithms\">\n"
                + "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence\">\n"
                + "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">\n"
                + "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second end\">\n"
                + "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Strands\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "##contig=<ID=chr2,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n";
        Files.writeString(dir.resolve("sv_calls.vcf"), header
                + "chr1\t1000\tdel_small\tN\t<DEL>\t.\tPASS\tEND=1300;SVTYPE=DEL;SVLEN=-300;ALGORITHMS=manta;EVIDENCE=PE,SR\tGT\t0/1\n"
                + "chr1\t5000\tdel_large\tN\t<DEL>\t.\tPASS\tEND=15000;SVTYPE=DEL;ALGORITHMS=depth;EVIDENCE=RD\tGT\t0/1\n"
                + "chr1\t20000\tdup1\tN\t<DUP>\t.\tPASS\tEND=22000;SVTYPE=DUP;ALGORITHMS=depth\tGT\t0/1\n"
                + "chr1\t30000\tins1\tN\t<INS>\t.\tPASS\tEND=30001;SVTYPE=INS;SVLEN=400;ALGORITHMS=manta;STRANDS=+-\tGT\t1/1\n"
                + "chr1\t40000\tinv1\tN\t<INV>\t.\tPASS\tEND=41000;SVTYPE=INV;ALGORITHMS=manta;STRANDS=++\tGT\t0/1\n"
                + "chr1\t50000\tbnd1\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END2=7000;ALGORITHMS=manta;STRANDS=+-\tGT\t0/1\n"
                + "chr2\t8000\tins_undef\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=-1;ALGORITHMS=melt\tGT\t0/1\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_bad.vcf"), header
                + "chr1\t1000\tdel_small\tN\t<DEL>\t.\tPASS\tEND=1300;SVTYPE=DEL;ALGORITHMS=manta\tGT\t0/1\n"
                + "chr1\t5000\tno_algorithms\tN\t<DEL>\t.\tPASS\tEND=15000;SVTYPE=DEL\tGT\t0/1\n",
                StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_rm.bed"),
                "chr1\t899\t1400\nchr1\t29989\t30010\nchr2\t6989\t7010\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_sd.bed"),
                "chr1\t4999\t10000\nchr1\t39999\t41000\n", StandardCharsets.UTF_8);
        final String columns = "NAME\tSVTYPE\tMIN_SIZE\tMAX_SIZE\tTRACKS\n";
        Files.writeString(dir.resolve("sv_strata.tsv"), columns
                + "DEL_small\tDEL\t50\t5000\tNA\n"
                + "DEL_large\tDEL\t5000\tNA\tNA\n"
                + "DUP_any\tDUP\tNA\tNA\tNA\n"
                + "INS_rm\tINS\tNA\tNA\tRM\n"
                + "INV_any\tINV\t-1\t-1\t-1\n"
                + "BND_rm\tBND\tNA\tNA\tRM\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_overlap.tsv"), columns
                + "DEL_any\tDEL\tNA\tNA\tNA\n"
                + "DEL_rm\tDEL\tNA\tNA\tRM\n"
                + "INV_sd\tINV\tNA\tNA\tSD\n"
                + "INS_rm\tINS\tNA\tNA\tRM,SD\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_track_names.list"), "RM\nSD\n", StandardCharsets.UTF_8);
        Files.writeString(dir.resolve("sv_track_files.list"),
                "/work/fixtures/sv_rm.bed\n/work/fixtures/sv_sd.bed\n", StandardCharsets.UTF_8);
    }

    /**
     * What `ApplyVQSR` reads: variants, the recal file `VariantRecalibrator` would have written for
     * them, and a tranches file per mode.
     *
     * The variants cover every class the mode test and the allele-specific path tell apart: SNPs,
     * an insertion, a deletion, a MIXED site, a multiallelic SNP, a SNP beside a spanning deletion,
     * and a SNP already filtered `LowQual`. The SITE recal file carries one record per site, ending
     * where the variant ends; the ALLELE recal file one per alternate, which is what `-AS` matches
     * on, and a site-mode run over it takes the first record at the right end. The site file is
     * written a second time without its index, which the first query refuses.
     *
     * `vqsr_snp_applied.vcf` is the REFERENCE's own first pass, `-AS` in SNP mode at 99.0, so a row
     * reading it is a second run: its header carries SNP tranche filters, which
     * `checkForPreviousApplyRecalRun` finds, and its records carry the lists an INDEL run copies.
     */
    static void vqsrFixtures(final Path dir) throws Exception {
        final String contig = "##contig=<ID=chr1,length=100000>\n";
        final String info = "##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"Log odds\">\n"
                + "##INFO=<ID=culprit,Number=1,Type=String,Description=\"Culprit\">\n"
                + "##INFO=<ID=POSITIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"Positive\">\n"
                + "##INFO=<ID=NEGATIVE_TRAIN_SITE,Number=0,Type=Flag,Description=\"Negative\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n";
        final String input = "##fileformat=VCFv4.2\n"
                + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + contig
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
                + "chr1\t100\trs100\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n"
                + "chr1\t200\t.\tC\tT\t40\t.\tDP=8\tGT\t1/1\n"
                + "chr1\t300\t.\tG\tGA\t30\tPASS\tDP=5\tGT\t0/1\n"
                + "chr1\t400\t.\tT\tC,TA\t20\tPASS\tDP=7\tGT\t1/2\n"
                + "chr1\t500\t.\tA\tC\t10\tLowQual\tDP=3\tGT\t0/1\n"
                + "chr1\t600\t.\tCA\tC\t60\t.\tDP=12\tGT\t0/1\n"
                + "chr1\t700\t.\tG\tA,T\t45\tPASS\tDP=9\tGT\t1/2\n"
                + "chr1\t800\t.\tA\tC,*\t35\tPASS\tDP=6\tGT\t1/2\n";
        Files.writeString(dir.resolve("vqsr_input.vcf"), input, StandardCharsets.UTF_8);

        final String recalHeader = "##fileformat=VCFv4.2\n" + info + contig
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        final String sites = recalHeader
                + "chr1\t100\t.\tA\t<VQSR>\t.\t.\tEND=100;POSITIVE_TRAIN_SITE;VQSLOD=5.21;culprit=QD\n"
                + "chr1\t200\t.\tC\t<VQSR>\t.\t.\tEND=200;VQSLOD=-2.40;culprit=FS\n"
                + "chr1\t300\t.\tG\t<VQSR>\t.\t.\tEND=300;VQSLOD=1.50;culprit=MQ\n"
                + "chr1\t400\t.\tT\t<VQSR>\t.\t.\tEND=400;VQSLOD=0.80;culprit=QD\n"
                + "chr1\t500\t.\tA\t<VQSR>\t.\t.\tEND=500;VQSLOD=2.00;culprit=SOR\n"
                + "chr1\t600\t.\tCA\t<VQSR>\t.\t.\tEND=601;NEGATIVE_TRAIN_SITE;VQSLOD=-4.10;culprit=FS\n"
                + "chr1\t700\t.\tG\t<VQSR>\t.\t.\tEND=700;VQSLOD=3.10;culprit=MQRankSum\n"
                + "chr1\t800\t.\tA\t<VQSR>\t.\t.\tEND=800;VQSLOD=0.20;culprit=QD\n";
        final String alleles = recalHeader
                + "chr1\t100\t.\tA\tG\t.\t.\tEND=100;POSITIVE_TRAIN_SITE;VQSLOD=5.21;culprit=QD\n"
                + "chr1\t200\t.\tC\tT\t.\t.\tEND=200;VQSLOD=-2.40;culprit=FS\n"
                + "chr1\t300\t.\tG\tGA\t.\t.\tEND=300;VQSLOD=1.50;culprit=MQ\n"
                + "chr1\t400\t.\tT\tC\t.\t.\tEND=400;VQSLOD=0.80;culprit=QD\n"
                + "chr1\t400\t.\tT\tTA\t.\t.\tEND=400;VQSLOD=-0.50;culprit=FS\n"
                + "chr1\t500\t.\tA\tC\t.\t.\tEND=500;VQSLOD=2.00;culprit=SOR\n"
                + "chr1\t600\t.\tCA\tC\t.\t.\tEND=601;NEGATIVE_TRAIN_SITE;VQSLOD=-4.10;culprit=FS\n"
                + "chr1\t700\t.\tG\tA\t.\t.\tEND=700;VQSLOD=3.10;culprit=MQRankSum\n"
                + "chr1\t700\t.\tG\tT\t.\t.\tEND=700;VQSLOD=-1.50;culprit=QD\n"
                + "chr1\t800\t.\tA\tC\t.\t.\tEND=800;VQSLOD=0.20;culprit=QD\n";
        for (final String[] recal : new String[][] {
                {"vqsr_sites.recal.vcf", sites}, {"vqsr_alleles.recal.vcf", alleles}}) {
            Files.writeString(dir.resolve(recal[0]), recal[1], StandardCharsets.UTF_8);
            new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                    .instanceMain(new String[] {"-I", dir.resolve(recal[0]).toString()});
        }
        Files.writeString(dir.resolve("vqsr_unindexed.recal.vcf"), sites, StandardCharsets.UTF_8);

        final String columns = "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,"
                + "filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n";
        for (final String model : new String[] {"SNP", "INDEL"}) {
            final String name = "VQSRTranche" + model;
            Files.writeString(dir.resolve(model.toLowerCase() + ".tranches"),
                    "# Variant quality score tranches file\n# Version number 5\n" + columns
                    + "90.00,0,3,0.0000,2.0000,4.0000," + name + "0.00to90.00," + model + ",10,9,0.9000\n"
                    + "99.00,0,4,0.0000,2.0000,1.0000," + name + "90.00to99.00," + model + ",10,9,0.9900\n"
                    + "99.90,0,5,0.0000,2.0000,-1.0000," + name + "99.00to99.90," + model + ",10,9,0.9990\n"
                    + "100.00,0,6,0.0000,2.0000,-3.0000," + name + "99.90to100.00," + model + ",10,10,1.0000\n",
                    StandardCharsets.UTF_8);
        }

        new org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR().instanceMain(new String[] {
                "-V", dir.resolve("vqsr_input.vcf").toString(),
                "--recal-file", dir.resolve("vqsr_alleles.recal.vcf").toString(),
                "--tranches-file", dir.resolve("snp.tranches").toString(),
                "--truth-sensitivity-filter-level", "99.0",
                "--use-allele-specific-annotations", "--mode", "SNP",
                "--add-output-vcf-command-line", "false",
                "--create-output-variant-index", "false",
                "-O", dir.resolve("vqsr_snp_applied.vcf").toString()});
    }

    /**
     * What `PathSeqBuildReferenceTaxonomy` reads: a reference whose contig names carry every form
     * the tool parses, a taxonomy dump as a tar.gz, and the two catalogs.
     *
     * The RefSeq catalog is gzipped and the GenBank one is NOT, because the tool gunzips a catalog
     * by its NAME and not by its bytes: the pair covers both branches of `makeReaderMaybeGzipped`.
     * The second reference drops the contigs only a catalog can place, so a row holding it and no
     * GenBank catalog still has taxa, and the second dump renames a species and adds a genus, so
     * the two dumps write different trees over the same contigs.
     */
    static void pathSeqTaxonomy(final Path dir) throws Exception {
        final String[][] contigs = {
                {"ref|NC_VIRUS.1|", "300"}, {"ref|NC_BACT.1|", "1000"}, {"ref|NC_SHORT.1|", "100"},
                {"taxid|562|", "800"}, {"ACC_PLAIN.1", "900"}, {"gi|9|ref|NC_BOTH.1|taxid|11234|", "700"}};
        pathSeqReference(dir.resolve("pathseq.fasta"), contigs);
        pathSeqKmerReference(dir.resolve("pathseq_kmers.fasta"));
        pathSeqReference(dir.resolve("pathseq2.fasta"), new String[][] {
                contigs[3], contigs[5], {"taxid|9606|", "600"}});

        final String names = String.join("\n",
                "1\t|\troot\t|\t\t|\tscientific name\t|",
                "2\t|\tBacteria\t|\t\t|\tscientific name\t|",
                "10239\t|\tViruses\t|\t\t|\tscientific name\t|",
                "562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|",
                "11234\t|\tMeasles morbillivirus\t|\t\t|\tscientific name\t|",
                "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|",
                "40674\t|\tMammalia\t|\t\t|\tscientific name\t|") + "\n";
        final String nodes = String.join("\n",
                "1\t|\t1\t|\tno rank\t|",
                "2\t|\t1\t|\tsuperkingdom\t|",
                "10239\t|\t1\t|\tsuperkingdom\t|",
                "562\t|\t2\t|\tspecies\t|",
                "11234\t|\t10239\t|\tspecies\t|",
                "40674\t|\t1\t|\tclass\t|",
                "9606\t|\t40674\t|\tspecies\t|") + "\n";
        taxdump(dir.resolve("taxdump.tar.gz"), names, nodes);
        taxdump(dir.resolve("taxdump2.tar.gz"),
                names.replace("Escherichia coli", "Escherichia coli K-12")
                        + "561\t|\tEscherichia\t|\t\t|\tscientific name\t|\n",
                nodes.replace("562\t|\t2\t|", "562\t|\t561\t|")
                        + "561\t|\t2\t|\tgenus\t|\n");

        try (final OutputStream out = new java.util.zip.GZIPOutputStream(
                Files.newOutputStream(dir.resolve("refseq.catalog.gz")))) {
            out.write(String.join("\n",
                    "11234\tsomething\tNC_VIRUS.1\tmore",
                    "562\tsomething\tNC_BACT.1\tmore",
                    "562\tsomething\tNC_SHORT.1\tmore").concat("\n")
                    .getBytes(StandardCharsets.UTF_8));
        }
        Files.writeString(dir.resolve("genbank.catalog"),
                "a\tACC_PLAIN.1\tc\td\te\tf\t9606\th\n", StandardCharsets.UTF_8);
    }

    /**
     * A reference whose bases are not periodic, which `PathSeqBuildKmers` needs.
     *
     * The corpus's other references repeat `ACGT`, and a periodic sequence has four distinct
     * 31-mers whatever its length: every row of that tool's array wrote a set of two entries and
     * forty-two bytes, which compares a serializer rather than a k-mer set. These bases come from
     * a fixed linear congruential sequence, so they are varied and still the same on every run.
     * One contig also carries a run of `N` and a lower-case stretch: a bad base costs a whole
     * window, and a lower-case one is not upper-cased before it is read.
     */
    static void pathSeqKmerReference(final Path fasta) throws Exception {
        final StringBuilder first = new StringBuilder();
        long state = 12345L;
        for (int i = 0; i < 2000; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            first.append("ACGT".charAt((int) ((state >>> 33) & 3)));
        }
        final StringBuilder second = new StringBuilder();
        for (int i = 0; i < 400; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            second.append("ACGT".charAt((int) ((state >>> 33) & 3)));
        }
        second.append("NNNNNNNNNN");
        for (int i = 0; i < 200; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            second.append("acgt".charAt((int) ((state >>> 33) & 3)));
        }
        try (final htsjdk.samtools.reference.FastaReferenceWriter writer =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(fasta)
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            writer.startSequence("host_1").appendBases(first.toString());
            writer.startSequence("host_2").appendBases(second.toString());
        }
    }

    static void pathSeqReference(final Path fasta, final String[][] contigs) throws Exception {
        try (final htsjdk.samtools.reference.FastaReferenceWriter writer =
                     new htsjdk.samtools.reference.FastaReferenceWriterBuilder()
                             .setFastaFile(fasta)
                             .setMakeFaiOutput(true)
                             .setMakeDictOutput(true)
                             .build()) {
            for (final String[] contig : contigs) {
                final StringBuilder bases = new StringBuilder();
                for (int i = 0; i < Integer.parseInt(contig[1]); i++) {
                    bases.append("ACGT".charAt(i % 4));
                }
                writer.startSequence(contig[0]).appendBases(bases.toString());
            }
        }
    }

    static void taxdump(final Path path, final String names, final String nodes) throws Exception {
        try (final org.apache.commons.compress.archivers.tar.TarArchiveOutputStream tar =
                     new org.apache.commons.compress.archivers.tar.TarArchiveOutputStream(
                             new java.util.zip.GZIPOutputStream(Files.newOutputStream(path)))) {
            for (final String[] entry : new String[][] {{"names.dmp", names}, {"nodes.dmp", nodes}}) {
                final byte[] bytes = entry[1].getBytes(StandardCharsets.UTF_8);
                final org.apache.commons.compress.archivers.tar.TarArchiveEntry header =
                        new org.apache.commons.compress.archivers.tar.TarArchiveEntry(entry[0]);
                header.setSize(bytes.length);
                // A fixed time, so the tarball is the same bytes on every run.
                header.setModTime(0L);
                tar.putArchiveEntry(header);
                tar.write(bytes);
                tar.closeArchiveEntry();
            }
        }
    }
}
