/*
 * PileupSummary and its table, taken from the reference.
 *
 * The second shared piece of CalculateContamination, and the record every contamination tool
 * passes around: six columns, a sample carried as metadata, and a handful of quantities derived
 * from the counts.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE ALLELE FREQUENCY IS WRITTEN BY DataLine.set(int, double), WHOSE ROUNDING BRANCH IS DEAD.
 *     The method computes `Math.round(value)`, writes the short form when it matches, and then
 *     RETURNS `set(index, Double.toString(value))`, which overwrites what it just wrote. So 1.0 is
 *     written `1.0` and never `1`, and a port that implements the method as written would differ;
 *   - AND AN INTEGER COLUMN IS NOT THE SAME PATH, so position and the counts come out with no
 *     decimal point at all;
 *   - THE METADATA TAG IS `SAMPLE`, in upper case, and the whole comment line goes through the same
 *     csv writer as a row: a sample name holding a tab makes the writer QUOTE THE COMMENT LINE
 *     ITSELF, `"#<METADATA>SAMPLE=has\ta tab"`, which the reader then parses back intact;
 *   - THE TOTAL IS DERIVED, NOT STORED: `refCount + altCount + otherAltsCount` in this constructor,
 *     while the VariantContext one derives otherAlts from the total instead. The table has no
 *     total column, so a read record recomputes it;
 *   - THE ALT FRACTION GUARDS AGAINST ZERO and returns 0 for an empty pileup rather than NaN, so
 *     the minor allele fraction of an empty site is 0 and not NaN, while a frequency of 2 gives a
 *     reference frequency of -1 because nothing bounds it;
 *   - GATHERING TAKES THE SAMPLE FROM THE FIRST FILE, so a second sample is a UserException.BadInput
 *     naming both, A FILE WITH NO SAMPLE FIRST is an IllegalArgumentException from the writer
 *     ("Null object is not allowed here."), AND THE SAME FILE SECOND IS A NullPointerException,
 *     because the comparison calls equals on the value that is missing. Three refusals for one
 *     mistake, depending only on where the file sits in the list;
 *   - AND A READER TAKES THE MISSING SAMPLE IN ITS STRIDE, returning null beside the records;
 *   - THE COMPARATOR ORDERS BY THE DICTIONARY'S INDEX, which is -1 for a contig the dictionary
 *     does not have, so an unknown contig sorts BEFORE every known one;
 *   - AND getDouble's REFUSAL SAYS "expected int value". The message is the int one, copied twice
 *     into the double getter, so a malformed allele frequency is reported as a bad integer.
 *
 * Output:
 *
 *     written\t<label>\t<the whole file, escaped>
 *     derived\t<label>\t<row>\t<total>,<altFraction>,<minorAlleleFraction>,<refFrequency>
 *     read\t<label>\t<sample>\t<row>\t<contig>,<position>,<ref>,<alt>,<other>,<frequency>
 *     gathered\t<label>\t<the whole file, escaped>
 *     sorted\t<label>\t<contig:position, comma separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PileupSummaryDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class PileupSummaryDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("pileupsummary-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PileupSummaryDump: the pileup summary record and its table, from the reference");

        // Every kind of allele frequency the writer might spell differently, and one empty site.
        final List<PileupSummary> records = new ArrayList<>();
        records.add(new PileupSummary("chr1", 100, 10, 5, 0, 1.0));
        records.add(new PileupSummary("chr1", 200, 10, 5, 2, 0.5));
        records.add(new PileupSummary("chr1", 300, 0, 0, 0, 0.0));
        records.add(new PileupSummary("chr1", 400, 1, 3, 0, 1.0 / 3));
        records.add(new PileupSummary("chr1", 500, 7, 0, 0, 1e-7));
        records.add(new PileupSummary("chr1", 600, 2, 2, 0, 2.0));
        records.add(new PileupSummary("chr1", 700, 3, 1, 0, Double.NaN));
        records.add(new PileupSummary("chr1", 800, 3, 1, 0, Double.POSITIVE_INFINITY));
        write(dir, "frequencies", "sample-a", records);
        for (int i = 0; i < records.size(); i++) {
            final PileupSummary record = records.get(i);
            System.out.printf("derived\t%s\t%d\t%s,%s,%s,%s%n", "frequencies", i,
                    String.valueOf(record.getTotalCount()),
                    String.valueOf(record.getAltFraction()),
                    String.valueOf(record.getMinorAlleleFraction()),
                    String.valueOf(record.getRefFrequency()));
        }
        readFile(dir, "frequencies", dir.resolve("frequencies.table"));

        // A table with no records at all, whose header and metadata are still written.
        write(dir, "empty", "sample-a", List.of());
        readFile(dir, "empty", dir.resolve("empty.table"));

        // A sample name that would have to be quoted if it were a value rather than metadata.
        write(dir, "odd-sample", "has\ta tab", List.of(
                new PileupSummary("chr1", 1, 1, 1, 1, 0.25)));
        readFile(dir, "odd-sample", dir.resolve("odd-sample.table"));

        // Gathering: two files of the same sample, in the order they are given.
        write(dir, "part-one", "sample-a", List.of(
                new PileupSummary("chr1", 100, 10, 5, 0, 0.5),
                new PileupSummary("chr1", 200, 10, 5, 0, 0.5)));
        write(dir, "part-two", "sample-a", List.of(
                new PileupSummary("chr2", 100, 1, 1, 0, 0.25)));
        write(dir, "other-sample", "sample-b", List.of(
                new PileupSummary("chr1", 300, 1, 1, 0, 0.25)));
        gather(dir, "same-sample",
                List.of(dir.resolve("part-two.table"), dir.resolve("part-one.table")));
        gather(dir, "two-samples",
                List.of(dir.resolve("part-one.table"), dir.resolve("other-sample.table")));

        // A file with the columns and no sample metadata, first and then second in the list.
        final Path nameless = dir.resolve("nameless.table");
        Files.writeString(nameless,
                "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                        + "chr1\t100\t10\t5\t0\t0.5\n", StandardCharsets.UTF_8);
        System.out.printf("written\t%s\t%s%n", "nameless",
                ReferenceQueryDump.escape(Files.readString(nameless)));
        gather(dir, "nameless-first", List.of(nameless, dir.resolve("part-one.table")));
        gather(dir, "nameless-second", List.of(dir.resolve("part-one.table"), nameless));
        readFile(dir, "nameless", nameless);

        // The comparator, with a contig the dictionary does not have.
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 250_000_000),
                new SAMSequenceRecord("chr2", 240_000_000),
                new SAMSequenceRecord("chr10", 130_000_000)));
        final List<PileupSummary> unsorted = new ArrayList<>(List.of(
                new PileupSummary("chr10", 5, 1, 1, 0, 0.5),
                new PileupSummary("chr2", 50, 1, 1, 0, 0.5),
                new PileupSummary("chrUn", 1, 1, 1, 0, 0.5),
                new PileupSummary("chr1", 900, 1, 1, 0, 0.5),
                new PileupSummary("chr1", 100, 1, 1, 0, 0.5)));
        unsorted.sort(new PileupSummary.PileupSummaryComparator(dictionary));
        final List<String> places = new ArrayList<>();
        for (final PileupSummary record : unsorted) {
            places.add(record.getContig() + ":" + record.getStart());
        }
        System.out.printf("sorted\t%s\t%s%n", "dictionary-order", String.join(",", places));

        // A row whose position is not a number, which only the reader refuses.
        final Path broken = dir.resolve("broken.table");
        Files.writeString(broken,
                "#<METADATA>sample=sample-a\n"
                        + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                        + "chr1\tx\t10\t5\t0\t0.5\n", StandardCharsets.UTF_8);
        readFile(dir, "broken", broken);

        // And one whose allele frequency is not a number, which is the other getter.
        final Path brokenFrequency = dir.resolve("broken-frequency.table");
        Files.writeString(brokenFrequency,
                "#<METADATA>sample=sample-a\n"
                        + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
                        + "chr1\t100\t10\t5\t0\tnot-a-number\n", StandardCharsets.UTF_8);
        readFile(dir, "broken-frequency", brokenFrequency);
    }

    /** Write a table through PileupSummary.writeToFile, with its sample metadata. */
    static void write(final Path dir, final String label, final String sample,
                      final List<PileupSummary> records) throws Exception {
        final Path file = dir.resolve(label + ".table");
        PileupSummary.writeToFile(sample, records, file.toFile());
        System.out.printf("written\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
    }

    /** Read a table back, with the sample the metadata carried. */
    static void readFile(final Path dir, final String label, final Path file) {
        final ImmutablePair<String, List<PileupSummary>> pair;
        try {
            pair = PileupSummary.readFromFile(file.toFile());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        final List<PileupSummary> records = pair.getRight();
        if (records.isEmpty()) {
            System.out.printf("read\t%s\t%s\t\t%n", label,
                    ReferenceQueryDump.escape(String.valueOf(pair.getLeft())));
            return;
        }
        for (int i = 0; i < records.size(); i++) {
            final PileupSummary record = records.get(i);
            System.out.printf("read\t%s\t%s\t%d\t%s,%d,%d,%d,%d,%s%n", label,
                    ReferenceQueryDump.escape(String.valueOf(pair.getLeft())), i,
                    record.getContig(), record.getStart(), record.getRefCount(),
                    record.getAltCount(), record.getOtherAltCount(),
                    String.valueOf(record.getAlleleFrequency()));
        }
    }

    /** Gather several tables into one, which is what GatherPileupSummaries does. */
    static void gather(final Path dir, final String label, final List<Path> inputs) {
        final Path output = dir.resolve(label + "-gathered.table");
        final List<File> files = new ArrayList<>();
        for (final Path input : inputs) {
            files.add(input.toFile());
        }
        try {
            PileupSummary.writeToFile(files, output.toFile());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        try {
            System.out.printf("gathered\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(Files.readString(output)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label + "-read", e.getClass().getName(),
                    String.valueOf(e.getMessage()));
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
