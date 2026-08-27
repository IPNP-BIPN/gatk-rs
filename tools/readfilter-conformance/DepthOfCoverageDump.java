/*
 * DepthOfCoverage's tables, taken from the reference.
 *
 * How deep the reads lie over a reference, counted per locus, per sample and per interval. The tool
 * writes a family of files rather than one, and which of them appear is decided by four `omit`
 * arguments rather than by the data.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE PER-LOCUS TABLE IS EVERY BASE OF EVERY INTERVAL, including those no read reaches, so its
 *     length is the intervals' and not the reads';
 *   - --min-base-quality FILTERS A BASE RATHER THAN A READ, so one low base lowers the depth at one
 *     locus and nowhere else;
 *   - --max-base-quality DOES THE SAME FROM ABOVE, which is a filter no other coverage tool has;
 *   - THE PARTITION IS THE SAMPLE, so two read groups of one sample are one column;
 *   - --print-base-counts ADDS A PER-BASE BREAKDOWN in a fixed `A: C: G: T: N:` order, each pair
 *     followed by a space, so every row ends with one;
 *   - THE CUMULATIVE TABLES COUNT LOCI AT OR ABOVE EACH DEPTH, so their first column is every
 *     locus and they fall monotonically;
 *   - EACH `omit` ARGUMENT REMOVES ITS OWN FILE and leaves the others, so the set of files written
 *     is a function of the arguments alone;
 *   - --omit-interval-statistics AND --calculate-coverage-over-genes ARE MUTUALLY EXCLUSIVE, and
 *     the refusal names both;
 *   - A QUALITY OUTSIDE 0..127 IS REFUSED BY THE ARGUMENT ITSELF, before a read is seen;
 *   - AND THE SUMMARY CARRIES A `Total` ROW whose quantiles are `N/A`, because a quantile over
 *     the partitions is not a quantile over anything the tool measured.
 *
 * Output:
 *
 *     files\t<label>=<the suffixes written, sorted, comma separated>
 *     out\t<label>.<suffix>=<that file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: DepthOfCoverageDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.coverage.DepthOfCoverage;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DepthOfCoverageDump {

    static final int CONTIG_LENGTH = 199980;

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // TWO read groups of ONE sample, and one of another: the partition is the SAMPLE.
        header.addReadGroup(group("rgA1", "sampleA"));
        header.addReadGroup(group("rgA2", "sampleA"));
        header.addReadGroup(group("rgB", "sampleB"));
        return header;
    }

    static SAMReadGroupRecord group(final String id, final String sample) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        group.setSample(sample);
        group.setPlatform("ILLUMINA");
        return group;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String readGroup,
                          final int start, final int length, final int baseQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(length + "M");
        record.setMappingQuality(60);
        final StringBuilder bases = new StringBuilder();
        while (bases.length() < length) {
            bases.append("ACGT");
        }
        record.setReadBases(bases.substring(0, length).getBytes(StandardCharsets.UTF_8));
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", readGroup);
        return record;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("depth-of-coverage-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# DepthOfCoverageDump: how deep the reads lie over a reference");

        final Path fasta = writeReference(dir);
        writeDictionary(dir);

        final SAMFileHeader header = header();
        final List<SAMRecord> records = new ArrayList<>();
        // sampleA, two read groups over the same ten bases: one column of depth two.
        records.add(read(header, "a1", "rgA1", 1000, 10, 30));
        records.add(read(header, "a2", "rgA2", 1000, 10, 30));
        // sampleB over the first five only, so the depth steps down inside the interval.
        records.add(read(header, "b1", "rgB", 1000, 5, 30));
        // A read whose bases are all at quality 5, below the floor the runs use.
        records.add(read(header, "lowq", "rgB", 1000, 10, 5));
        // A read whose bases are all at quality 60, above the ceiling one run sets.
        records.add(read(header, "highq", "rgB", 1000, 10, 60));
        records.sort((a, b) -> Integer.compare(a.getAlignmentStart(), b.getAlignmentStart()));

        final Path bam = dir.resolve("reads.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }

        // An interval that runs past the reads on both sides, so the table carries uncovered bases.
        final Path intervals = write(dir, "intervals.list", "chr1:995-1015\n");

        run(dir, "default", bam, fasta, intervals, List.of());
        run(dir, "base-counts", bam, fasta, intervals, List.of("--print-base-counts", "true"));
        run(dir, "min-baseq", bam, fasta, intervals, List.of("--min-base-quality", "10"));
        run(dir, "max-baseq", bam, fasta, intervals, List.of("--max-base-quality", "40"));
        run(dir, "omit-locus-table", bam, fasta, intervals, List.of("--omit-locus-table", "true"));
        run(dir, "omit-per-base", bam, fasta, intervals,
                List.of("--omit-depth-output-at-each-base", "true"));
        run(dir, "omit-per-sample", bam, fasta, intervals,
                List.of("--omit-per-sample-statistics", "true"));
        run(dir, "omit-intervals", bam, fasta, intervals,
                List.of("--omit-interval-statistics", "true"));

        // The two mutually exclusive arguments.
        run(dir, "both-interval-arguments", bam, fasta, intervals, List.of(
                "--omit-interval-statistics", "true",
                "--calculate-coverage-over-genes", intervals.toString()));
        // A quality the argument itself refuses.
        run(dir, "quality-too-high", bam, fasta, intervals,
                List.of("--min-base-quality", "200"));
        run(dir, "quality-negative", bam, fasta, intervals,
                List.of("--min-base-quality", "-1"));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final Path intervals, final List<String> extra) throws Exception {
        final Path base = dir.resolve("out-" + label);
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-R", fasta.toString(),
                "-L", intervals.toString(),
                "-O", base.toString()));
        argv.addAll(extra);
        try {
            new DepthOfCoverage().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        // Whatever the run wrote, by suffix.
        final List<String> suffixes = new ArrayList<>();
        try (final java.util.stream.Stream<Path> listing = Files.list(dir)) {
            for (final Path path : listing.sorted().toList()) {
                final String name = path.getFileName().toString();
                if (name.startsWith("out-" + label)) {
                    suffixes.add(name.substring(("out-" + label).length()));
                }
            }
        }
        java.util.Collections.sort(suffixes);
        System.out.printf("files\t%s=%s%n", label, String.join(",", suffixes));
        for (final String suffix : suffixes) {
            final Path path = dir.resolve("out-" + label + suffix);
            System.out.printf("out\t%s%s=%s%n", label, suffix.isEmpty() ? ".base" : suffix,
                    ReferenceQueryDump.escape(masked(Files.readString(path), dir)));
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        return fasta;
    }

    /** One contig, because the dictionary and the FASTA index must agree on their number. */
    static void writeDictionary(final Path dir) throws Exception {
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH)));
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
