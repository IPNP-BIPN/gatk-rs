/*
 * MethylationTypeCaller's output, taken from the reference.
 *
 * A LocusWalker over a bisulfite-sequenced BAM: at every reference C it counts the reads that stayed
 * C (unconverted, methylated) against those that became T (converted), on the FORWARD strand only,
 * and at every reference G it counts G against A on the REVERSE strand only. Each covered site
 * becomes one VCF record.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE STRAND IS CHOSEN BY THE REFERENCE BASE, not by the reads: a C is counted from forward
 *     reads alone and a G from reverse reads alone, so a C covered only by reverse reads produces
 *     NO RECORD AT ALL while its coverage is still counted in DP;
 *   - A REFERENCE BASE THAT IS NEITHER C NOR G IS SKIPPED before anything is counted;
 *   - A SITE WITH NO METHYLATED COVERAGE IS SKIPPED, `unconverted + converted > 0`, so a C covered
 *     only by A and G reads writes nothing even though DP would have been positive;
 *   - DP IS THE WHOLE PILEUP, both strands and every base, while the two counts are one strand and
 *     two bases: the three numbers of one record do not add up, on purpose;
 *   - THE FORWARD CONTEXT IS getBases(0, 2), three bases starting at the site, AND THE REVERSE
 *     CONTEXT IS getBases(2, 0) REVERSE COMPLEMENTED, three bases ending at the site read backwards.
 *     A site within two bases of a contig edge gets a SHORTER context rather than an error;
 *   - THE ALT ALLELE IS THE CONVERTED BASE, T for a C and A for a G, whether or not any read shows
 *     it, so a fully unconverted site still writes an alt;
 *   - ITS FILTERS ARE THE LOCUS WALKER'S, which drop a read with no cigar or an unmapped one, and
 *     the traversal's own rules drop deletions and Ns from the pileup;
 *   - THE SAMPLES OF THE VCF HEADER ARE THE READ GROUPS' SAMPLES, sorted then deduplicated, and the
 *     records carry no genotypes at all, so a header with two samples has two empty columns;
 *   - AND THE DEFAULT HEADER CARRIES A ##GATKCommandLine LINE WITH THE RUN'S DATE, which is why the
 *     comparable run passes --add-output-vcf-command-line false. That line is dumped separately with
 *     its date masked, because its shape is measurable and its content is not.
 *
 * Output:
 *
 *     reference\t<the reference bases of chr1>
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MethylationTypeCallerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.walkers.MethylationTypeCaller;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MethylationTypeCallerDump {

    /**
     * Sixty bases with C and G in every context the tool distinguishes, including a C two bases from
     * the end so its context is truncated rather than padded.
     */
    static final String CHR1 = "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGCTCTAGCATCC";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("methylation-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MethylationTypeCallerDump: MethylationTypeCaller's output, from the reference");
        System.out.printf("reference\t%s%n", CHR1);

        final Path fasta = writeReference(dir);

        final Path bisulfite = dir.resolve("bisulfite.bam");
        buildBisulfite(bisulfite.toFile());
        emitFixture(dir, bisulfite, "bisulfite");

        // A second sample, so the header's sample list is two names rather than one.
        final Path samples = dir.resolve("samples.bam");
        buildTwoSamples(samples.toFile());
        emitFixture(dir, samples, "samples");

        run(dir, bisulfite, fasta, "plain",
                new String[] {"--add-output-vcf-command-line", "false"});
        run(dir, bisulfite, fasta, "interval",
                new String[] {"--add-output-vcf-command-line", "false", "-L", "chr1:18-24"});
        run(dir, samples, fasta, "two-samples",
                new String[] {"--add-output-vcf-command-line", "false"});
        // The default, whose header carries the date the run happened.
        run(dir, bisulfite, fasta, "with-command-line", new String[] {});
    }

    static void emitFixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    /**
     * Reads over the reference's C and G runs.
     *
     * The reference at 17..20 is GGGG and at 21..24 is CCCC, which is where the two strands are
     * exercised; the C at 59 is two bases from the end, where the context runs short.
     */
    static void buildBisulfite(final File file) {
        final SAMFileHeader header = header("s1", null);
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Forward reads over the CCCC run: one keeps its Cs, one has them converted to T.
            writer.addAlignment(read(header, "fwd-unconverted", 17, "10M",
                    CHR1.substring(16, 26), 0, "rg1"));
            writer.addAlignment(read(header, "fwd-converted", 17, "10M",
                    convert(CHR1.substring(16, 26), 'C', 'T'), 0, "rg1"));
            // A reverse read over the same run: it is counted at the Gs and not at the Cs.
            writer.addAlignment(read(header, "rev-unconverted", 17, "10M",
                    CHR1.substring(16, 26), 0x10, "rg1"));
            writer.addAlignment(read(header, "rev-converted", 17, "10M",
                    convert(CHR1.substring(16, 26), 'G', 'A'), 0x10, "rg1"));
            // A forward read whose bases at the Cs are neither C nor T, so those sites have DP but
            // no methylated coverage.
            writer.addAlignment(read(header, "fwd-other", 21, "4M",
                    "GGGG", 0, "rg1"));
            // An unmapped read, which the locus walker's filters drop; it still carries a position,
            // so it travels in coordinate order like any other record.
            writer.addAlignment(read(header, "unmapped", 21, "4M", "CCCC", 0x4, "rg1"));
            // A read with a deletion across a C, which the pileup does not count as a base.
            writer.addAlignment(read(header, "deletion", 29, "3M2D5M",
                    CHR1.substring(28, 31) + CHR1.substring(33, 38), 0, "rg1"));
            // A read carrying N over a C.
            writer.addAlignment(read(header, "with-n", 40, "5M",
                    "N" + CHR1.substring(40, 44), 0, "rg1"));
            // The C two bases from the end of the contig, where the context is truncated.
            writer.addAlignment(read(header, "contig-end", 55, "6M",
                    CHR1.substring(54, 60), 0, "rg1"));
        }
    }

    /** Two read groups with different samples, so the header lists two. */
    static void buildTwoSamples(final File file) {
        final SAMFileHeader header = header("s2", "s1");
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            writer.addAlignment(read(header, "one", 17, "10M", CHR1.substring(16, 26), 0, "rg1"));
            writer.addAlignment(read(header, "two", 17, "10M",
                    convert(CHR1.substring(16, 26), 'C', 'T'), 0, "rg2"));
        }
    }

    static SAMFileHeader header(final String first, final String second) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample(first);
        header.addReadGroup(group);
        if (second != null) {
            final SAMReadGroupRecord other = new SAMReadGroupRecord("rg2");
            other.setSample(second);
            header.addReadGroup(other);
        }
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases, final int flags,
                          final String readGroup) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, (byte) 35);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", readGroup);
        return record;
    }

    /** Every occurrence of one base replaced by another, which is what bisulfite conversion does. */
    static String convert(final String bases, final char from, final char to) {
        return bases.replace(from, to);
    }

    /** One run of the tool, with every line of the VCF it wrote. */
    static void run(final Path dir, final Path input, final Path fasta, final String label,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve("MethylationTypeCaller." + label + ".vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-R", fasta.toString(), "-O", output.toString(),
                "--use-jdk-inflater", "true", "--use-jdk-deflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new MethylationTypeCaller().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }

        for (final String line : Files.readAllLines(output)) {
            if (line.startsWith("##GATKCommandLine")) {
                // The date is the run's own, so only the shape can be compared.
                System.out.printf("commandline\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(line.replaceAll("Date=\"[^\"]*\"", "Date=\"MASKED\"")));
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + CHR1.length() + "\n",
                StandardCharsets.UTF_8);
        return fasta;
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
