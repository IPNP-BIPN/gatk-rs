/*
 * GetPileupSummaries, taken from the reference.
 *
 * The tool that produces what CalculateContamination consumes, so this closes that chain from both
 * ends: the table format and the model are already pinned, and this is where the table comes from.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE DEFAULT FILTER SET IS ELEVEN FILTERS, not the walker's usual two, and one of them is
 *     parameterised: MappingQualityReadFilter at 50, which is far above anything else in the port.
 *     A read at mapping quality 60 survives it and one at 30 does not;
 *   - A SITE IS SUMMARISED ONLY IF ITS FIRST VARIANT IS BIALLELIC AND A SNP, so a triallelic site
 *     and an indel are both skipped, and only the FIRST variant at a locus is looked at;
 *   - THE ALLELE FREQUENCY BOUNDS ARE STRICT AT BOTH ENDS: `min < af && af < max`, so a site whose
 *     AF is exactly the default 0.01 or exactly 0.2 is excluded;
 *   - A VARIANT WITH NO AF IS SKIPPED, and if NO variant had one the tool throws at the END of the
 *     traversal rather than at the start;
 *   - A HEADER WITHOUT AF IS REFUSED BEFORE ANY READ IS TOUCHED, which is a different exception
 *     from the one above and a different moment;
 *   - THE COUNTS COME FROM getBaseCounts, which counts A, C, G and T only: a deletion at the site
 *     is skipped and an N is not counted, so `otherAlts` is total minus ref minus alt over those
 *     four bases and not over the pileup's depth;
 *   - AND THE TABLE CARRIES THE SAMPLE AS METADATA, read from the read group of the FIRST sample
 *     the header names.
 *
 * Output:
 *
 *     table\t<label>\t<the whole table, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GetPileupSummariesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GetPileupSummariesDump {

    /**
     * The population VCF, against chr1 of ReadWalkerDump.FASTA. The fixture's reads are
     * `ACGTACGTAC` at 10, 65, 120, 140, 150, 160, 170 and 180 on chr1.
     *
     *  - 12 G>C AF=0.1, in range, and covered by the read at 10;
     *  - 13 T>A AF=0.01, exactly the lower bound, which is exclusive;
     *  - 14 A>C AF=0.2, exactly the upper bound, which is exclusive;
     *  - 15 C>G AF=0.005, below the range;
     *  - 16 G>T AF=0.5, above it;
     *  - 17 T>A,C AF=0.1, triallelic;
     *  - 18 AC>A AF=0.1, an indel;
     *  - 19 G>C with no AF at all;
     *  - 66 C>A AF=0.15, covered by a read in the soft-masked stretch.
     */
    static final String VARIANTS =
            "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "##contig=<ID=chr2,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t12\t.\tG\tC\t50\tPASS\tAF=0.1\n"
            + "chr1\t13\t.\tT\tA\t50\tPASS\tAF=0.01\n"
            + "chr1\t14\t.\tA\tC\t50\tPASS\tAF=0.2\n"
            + "chr1\t15\t.\tC\tG\t50\tPASS\tAF=0.005\n"
            + "chr1\t16\t.\tG\tT\t50\tPASS\tAF=0.5\n"
            + "chr1\t17\t.\tT\tA,C\t50\tPASS\tAF=0.1\n"
            + "chr1\t18\t.\tAC\tA\t50\tPASS\tAF=0.1\n"
            + "chr1\t19\t.\tG\tC\t50\tPASS\t.\n"
            + "chr1\t66\t.\tC\tA\t50\tPASS\tAF=0.15\n";

    /** The same VCF with no AF in the header at all. */
    static final String NO_AF_HEADER =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t12\t.\tG\tC\t50\tPASS\t.\n";

    /** A header that declares AF and records that never carry one. */
    static final String NO_AF_RECORDS =
            "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t12\t.\tG\tC\t50\tPASS\t.\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("getpileupsummaries-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        final Path variants = write(dir, "variants.vcf", VARIANTS);
        final Path noHeader = write(dir, "no-af-header.vcf", NO_AF_HEADER);
        final Path noRecords = write(dir, "no-af-records.vcf", NO_AF_RECORDS);

        System.out.println("# GetPileupSummariesDump: the table CalculateContamination reads");

        run("default", dir, fasta, bam, variants, "-L", "chr1");
        // The bounds opened up, so the two sites at the defaults' edges and the one above come in.
        run("wide-bounds", dir, fasta, bam, variants, "-L", "chr1",
                "--minimum-population-allele-frequency", "0.001",
                "--maximum-population-allele-frequency", "0.9");
        // A window holding only the site in the soft-masked stretch.
        run("one-site", dir, fasta, bam, variants, "-L", "chr1:66-66");
        // A window with reads but no variants at all, which writes a table of nothing.
        run("no-variants", dir, fasta, bam, variants, "-L", "chr1:150-160");
        // The two refusals, which happen at opposite ends of the run.
        run("header-without-af", dir, fasta, bam, noHeader, "-L", "chr1");
        run("records-without-af", dir, fasta, bam, noRecords, "-L", "chr1");

        // A second fixture stacking six reads over position 12, which is the only way the counting
        // rule is visible: the shared fixture covers every site once.
        final java.nio.file.Path stack = dir.resolve("stack.bam");
        buildStack(stack.toFile());
        run("stacked", dir, fasta, stack, variants, "-L", "chr1:12-12");
        // The same stack with the mapping quality filter's threshold in play: one read is at 30.
        run("stacked-window", dir, fasta, stack, variants, "-L", "chr1:12-16");
    }

    /**
     * Six reads over chr1:10-15, differing at position 12, which is `G` in the reference.
     *
     * Two carry `G`, one `C` (the alternate), one `A` (an other-alt), one `N` (which getBaseCounts
     * does not count at all) and one deletes the base (which it skips). A seventh sits at mapping
     * quality 30, below the filter's threshold of 50, and never reaches apply.
     */
    static void buildStack(final java.io.File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        dictionary.addSequence(new SAMSequenceRecord("chr2", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("stacked1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(stacked(header, "s001", "ACGTAC", "6M", 60));
        records.add(stacked(header, "s002", "ACGTAC", "6M", 60));
        records.add(stacked(header, "s003", "ACCTAC", "6M", 60));
        records.add(stacked(header, "s004", "ACATAC", "6M", 60));
        records.add(stacked(header, "s005", "ACNTAC", "6M", 60));
        // A deletion at 12: two matched bases, one deleted, then three more.
        records.add(stacked(header, "s006", "ACTAC", "2M1D3M", 60));
        // Below the mapping quality threshold, so it is filtered out entirely.
        records.add(stacked(header, "s007", "ACGTAC", "6M", 30));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord stacked(final SAMFileHeader header, final String name, final String bases,
                             final String cigar, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(10);
        record.setCigarString(cigar);
        record.setMappingQuality(mappingQuality);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, text.getBytes());
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", path.toString()});
        return path;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final Path variants, final String... extra) throws Exception {
        final Path out = dir.resolve("table-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-I", bam.toString(),
                "-V", variants.toString(),
                "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new GetPileupSummaries().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
