/*
 * CalculateMixingFractions' table, taken from the reference.
 *
 * A VariantWalker that counts alt reads at singleton het SNPs, one bucket per sample, and divides
 * each sample's alt fraction by the sum of all of them.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE ROW ORDER IS A HASH BUCKET ORDER. `sampleCounts` is a `HashMap<String, ...>` and
 *     `onTraversalSuccess` walks its `entrySet()`, so the table's rows come out in the order
 *     `String.hashCode()` puts the sample names in a sixteen-bucket table, not in the header's
 *     order and not sorted;
 *   - A SAMPLE WITH NO COUNTED SITE POISONS EVERY ROW. Its alt fraction is 0/0, which is NaN, and
 *     the normalizer is the SUM of every sample's fraction, so one uncounted sample makes every
 *     mixing fraction NaN;
 *   - SINGLETON IS EITHER OF TWO TESTS, and the first one wins: `hasAttribute(AC) && AC[0] == 1`,
 *     or, only when that is false, exactly one het genotype. A record with AC=2 and one het is
 *     therefore still counted, and one with AC=1 and two hets is counted as well;
 *   - THE SAMPLE IS THE FIRST HET IN GENOTYPE ORDER, `findFirst()` over the genotypes, so a record
 *     with two hets is attributed entirely to whichever the header lists first;
 *   - AND A RECORD WITH NO HET AT ALL IS DROPPED after passing the singleton test, since
 *     `variantSample` is empty;
 *   - THE PILEUP IS BUILT BY WALKING EACH READ TO THE VARIANT'S START with an AlignmentStateMachine,
 *     so a read whose alignment passes over the position by a deletion contributes its deletion
 *     rather than nothing, and a read that never reaches the position contributes nothing at all;
 *   - VENDOR-FAILED READS ARE SKIPPED BY THE TOOL ITSELF, before the pileup rather than by a filter;
 *   - AND WITHOUT -I THERE IS NO PILEUP AT ALL, so every counted sample is 0/0 and every row is NaN.
 *
 * Output:
 *
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     input\t<label>\t<the whole input vcf, escaped>
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CalculateMixingFractionsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.CalculateMixingFractions;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CalculateMixingFractionsDump {

    /** Three samples whose names are not in hash order, which is the point of one of the rows. */
    static final String SAMPLES = "zebra\talpha\tmike";

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##contig=<ID=chr1,length=200>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + SAMPLES + "\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("calculatemixingfractions-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CalculateMixingFractionsDump: three buckets, and the order they print in");

        // The reads every run pileups over. Name, start, cigar, bases, flags, qualities.
        //  - at 20, four reads of which two carry the alt C;
        //  - at 40, three reads of which one is vendor-failed and carries the alt;
        //  - at 60, a read that spans the position with a deletion and one that ends before it;
        //  - at 80, one read carrying the alt.
        buildBam(dir, "reads", new String[][] {
            {"a1", "18", "10M", "AACAAAAAAA", "0", "IIIIIIIIII"},
            {"a2", "18", "10M", "AACAAAAAAA", "0", "IIIIIIIIII"},
            {"a3", "18", "10M", "AAAAAAAAAA", "0", "IIIIIIIIII"},
            {"a4", "18", "10M", "AAAAAAAAAA", "0", "IIIIIIIIII"},
            {"b1", "38", "10M", "AACAAAAAAA", "512", "IIIIIIIIII"},
            {"b2", "38", "10M", "AAAAAAAAAA", "0", "IIIIIIIIII"},
            {"b3", "38", "10M", "AACAAAAAAA", "0", "IIIIIIIIII"},
            {"c2", "50", "5M", "AAAAA", "0", "IIIII"},
            {"c1", "58", "2M1D7M", "AAAAAAAAA", "0", "IIIIIIIII"},
            {"d1", "78", "10M", "AACAAAAAAA", "0", "IIIIIIIIII"},
        });

        // Every shape the two tests can see.
        final Path everyShape = writeVcf(dir, "every-shape",
                // AC=1 and one het: counted, and attributed to zebra.
                "chr1\t20\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0",
                // No AC at all, exactly one het: counted through the second test, alpha's.
                "chr1\t40\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\t0/1\t0/0",
                // AC=2 with one het: the first test fails and the second passes, so it counts.
                "chr1\t60\t.\tA\tC\t50\tPASS\tAC=2\tGT\t0/0\t0/0\t0/1",
                // AC=1 with two hets: the first test passes, and the FIRST het takes it all.
                "chr1\t80\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t0/1\t0/0",
                // AC=1 with no het at all: singleton, but there is nobody to attribute it to.
                "chr1\t100\t.\tA\tC\t50\tPASS\tAC=1\tGT\t1/1\t0/0\t0/0",
                // Multi-allelic, so not biallelic.
                "chr1\t120\t.\tA\tC,G\t50\tPASS\tAC=1,0\tGT\t0/1\t0/0\t0/0",
                // An indel, so not a SNP.
                "chr1\t140\t.\tA\tACC\t50\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0");

        // One site only, so two of the three samples are never counted.
        final Path oneSite = writeVcf(dir, "one-site",
                "chr1\t20\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0");

        // Two sites that leave every sample with a count, which is the only way out of NaN.
        final Path everySampleCounted = writeVcf(dir, "every-sample-counted",
                "chr1\t20\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t0/0\t0/0",
                "chr1\t40\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/0\t0/1\t0/0",
                "chr1\t80\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/0\t0/0\t0/1");

        final Path bam = dir.resolve("reads.bam");

        run(dir, "every-shape", everyShape, bam, "every-shape.table");
        run(dir, "one-site", oneSite, bam, "one-site.table");
        run(dir, "every-sample-counted", everySampleCounted, bam, "every-sample-counted.table");
        // The same file with no reads behind it: every counted sample is 0/0.
        run(dir, "no-reads", everySampleCounted, null, "no-reads.table");
        // An interval that leaves only the third sample counted.
        run(dir, "one-interval", everySampleCounted, bam, "one-interval.table",
                "-L", "chr1:75-85");
        // The refusal.
        run(dir, "output-is-a-directory", everySampleCounted, bam, ".");
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
            new CalculateMixingFractions().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
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
