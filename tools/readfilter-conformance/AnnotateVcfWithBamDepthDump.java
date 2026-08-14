/*
 * AnnotateVcfWithBamDepth's output VCF, taken from the reference.
 *
 * A VariantWalker that writes its input back out with one INFO field added, counting the reads of a
 * separate bam that cover each variant. The count is not a pileup: it is five conditions on the
 * read's own coordinates and flags.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A READ ONE BASE LONG IS NEVER COUNTED, because the test is `read.getEnd() > read.getStart()`
 *     rather than `>=`: a 1M read sitting exactly on the variant contributes nothing;
 *   - THE READ MUST CONTAIN THE VARIANT'S WHOLE SPAN, `new SimpleInterval(read).contains(vc)` being
 *     containment and not overlap, so a read covering the first base of a four-base deletion is not
 *     counted, and a record carrying END is asked for that whole block;
 *   - DUPLICATES, VENDOR-FAILED AND UNMAPPED READS ARE EXCLUDED BY THE TOOL, in its own condition
 *     rather than by a read filter, so the count differs from the traversal's read set;
 *   - THE ANNOTATION IS WRITTEN EVEN WHEN IT IS ZERO, and a record no read covers carries
 *     `BAM_DEPTH=0` rather than nothing;
 *   - AN EXISTING BAM_DEPTH IS OVERWRITTEN, `VariantContextBuilder.attribute` replacing the value
 *     the input carried;
 *   - THE INFO COLUMN IS SORTED BY THE WRITER, so BAM_DEPTH lands where its name sorts and not at
 *     the end of the line;
 *   - AND THE HEADER GOES THROUGH A HashSet BEFORE IT IS WRITTEN, `new HashSet<>(
 *     inputHeader.getMetaDataInSortedOrder())`, which the writer then sorts again: the golden is
 *     what says whether anything of the input's order survives.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one line of the output vcf, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AnnotateVcfWithBamDepthDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithBamDepth;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AnnotateVcfWithBamDepthDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=BAM_DEPTH,Number=1,Type=Integer,Description=\"was already here\">\n"
                    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the block\">\n"
                    + "##INFO=<ID=ZZ,Number=1,Type=Integer,Description=\"sorts after BAM_DEPTH\">\n"
                    + "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"sorts before BAM_DEPTH\">\n"
                    + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=200>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("annotatevcfwithbamdepth-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AnnotateVcfWithBamDepthDump: five conditions, and where the number lands");

        // Name, start, cigar, bases, flags, qualities.
        //  - at 20, three reads that contain the site, one duplicate and one vendor-failed;
        //  - at 40, a read one base long sitting exactly on the site;
        //  - at 60, a read that covers the first base of the deletion but not its last;
        //  - at 100, nothing at all.
        buildBam(dir, "reads", new String[][] {
            {"a1", "15", "20M", "AAAAAAAAAAAAAAAAAAAA", "0", "IIIIIIIIIIIIIIIIIIII"},
            {"a2", "15", "20M", "AAAAAAAAAAAAAAAAAAAA", "0", "IIIIIIIIIIIIIIIIIIII"},
            {"a3", "15", "20M", "AAAAAAAAAAAAAAAAAAAA", "1024", "IIIIIIIIIIIIIIIIIIII"},
            {"a4", "15", "20M", "AAAAAAAAAAAAAAAAAAAA", "512", "IIIIIIIIIIIIIIIIIIII"},
            {"b1", "40", "1M", "A", "0", "I"},
            {"c1", "58", "5M", "AAAAA", "0", "IIIII"},
            {"c2", "58", "20M", "AAAAAAAAAAAAAAAAAAAA", "0", "IIIIIIIIIIIIIIIIIIII"},
            {"d1", "78", "20M", "AAAAAAAAAAAAAAAAAAAA", "0", "IIIIIIIIIIIIIIIIIIII"},
        });

        final Path variants = writeVcf(dir, "variants",
                // Three reads contain it, of which one is a duplicate and one vendor-failed.
                "chr1\t20\t.\tA\tC\t50\tPASS\tAA=1;ZZ=9\tGT\t0/1",
                // Only a one-base read sits here, and the > test excludes it.
                "chr1\t40\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                // A four-base deletion: one read covers its first base, the other its whole span.
                "chr1\t60\t.\tACCC\tA\t50\tPASS\t.\tGT\t0/1",
                // A block carrying END, which is what containment is asked about.
                "chr1\t80\t.\tA\t<DEL>\t50\tPASS\tEND=110\tGT\t0/1",
                // Nothing covers it at all.
                "chr1\t150\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                // A record that already carries the annotation the tool writes.
                "chr1\t160\t.\tA\tC\t50\tPASS\tBAM_DEPTH=99\tGT\t0/1");

        final Path bam = dir.resolve("reads.bam");

        run(dir, "annotated", variants, bam, "annotated.vcf");
        // The same file with no reads at all, where every record is BAM_DEPTH=0.
        run(dir, "no-reads", variants, null, "no-reads.vcf");
        // An interval, which changes which records are written rather than what they carry.
        run(dir, "one-interval", variants, bam, "one-interval.vcf", "-L", "chr1:15-45");
        // The refusal.
        run(dir, "output-is-a-directory", variants, bam, ".");
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void buildBam(final Path dir, final String label, final String[][] reads)
            throws Exception {
        final Path bam = dir.resolve(label + ".bam");
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            for (final String[] spec : reads) {
                writer.addAlignment(read(header, spec));
            }
        }
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 200))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("pooled");
        header.addReadGroup(group);
        return header;
    }

    /** name, start, cigar, bases, flags, quality characters. */
    static SAMRecord read(final SAMFileHeader header, final String[] spec) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(spec[0]);
        record.setFlags(Integer.parseInt(spec[4]));
        record.setReferenceName("chr1");
        record.setAlignmentStart(Integer.parseInt(spec[1]));
        record.setCigarString(spec[2]);
        record.setReadBases(spec[3].getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[spec[5].length()];
        for (int i = 0; i < quals.length; i++) {
            quals[i] = (byte) (spec[5].charAt(i) - 33);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static void run(final Path dir, final String label, final Path input, final Path bam,
                    final String output, final String... arguments) throws Exception {
        final Path file = dir.resolve(output);
        final List<String> all = new ArrayList<>(List.of(
                "-V", input.toString(), "-O", file.toString()));
        if (bam != null) {
            all.addAll(List.of("-I", bam.toString()));
        }
        all.addAll(List.of(arguments));
        try {
            new AnnotateVcfWithBamDepth().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, file);
    }

    static void print(final String label, final Path output) {
        final List<String> lines;
        try {
            lines = Files.readAllLines(output, StandardCharsets.UTF_8);
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
            return;
        }
        for (final String line : lines) {
            if (line.startsWith("##GATKCommandLine")) {
                System.out.printf("commandline\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(line.replaceAll("Date=\"[^\"]*\"", "Date=\"MASKED\"")));
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
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
