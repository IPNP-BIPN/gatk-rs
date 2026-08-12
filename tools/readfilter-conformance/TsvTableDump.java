/*
 * The tsv TableWriter and TableReader, taken from the reference.
 *
 * GATK's own table format, under every `.table` file the contamination tools read and write. It is
 * not a plain TSV: it is opencsv configured with a tab separator, a quote character and an escape
 * character, plus a comment convention that carries metadata.
 *
 * Nine behaviours this is built to catch.
 *
 *   - A VALUE IS QUOTED ONLY WHEN IT HAS TO BE, and what forces a quote is measured rather than
 *     assumed: a tab, a quote, a newline, and whatever else opencsv decides;
 *   - THE ESCAPE CHARACTER IS THE BACKSLASH inside quotes, so a value containing a quote comes back
 *     out with the quote escaped and the whole field quoted;
 *   - METADATA IS A COMMENT LINE WITH A TAG, `#<METADATA>key=value`, written before the header. A
 *     comment without that tag is a plain comment and contributes nothing to the map, which is why
 *     a hand-written `#sample=s1` is NOT metadata however much it looks like it;
 *   - AND THE MAP IS FILLED BY processCommentLine ITSELF, so a subclass that overrides it without
 *     calling super LOSES EVERY PAIR. The same file read twice, once with an override and once
 *     without, gives an empty map and `sample=s1`;
 *   - THE HEADER IS WRITTEN LAZILY, on the first record or on an explicit call, so a table with no
 *     records at all still has its header if the writer was asked for it;
 *   - A ROW WITH THE WRONG NUMBER OF COLUMNS IS A UserException.BadInput naming the line and its
 *     number, and the message is the reader's rather than opencsv's;
 *   - A MISSING COLUMN IN THE HEADER IS A DIFFERENT REFUSAL from a bad row, raised before any row
 *     is read;
 *   - A DOUBLE IS WRITTEN BY String.valueOf, so an integral double keeps its `.0` and a very small
 *     one comes out in exponent form;
 *   - AND AN EMPTY FIELD READS BACK AS THE EMPTY STRING, which `getInt` then refuses with its own
 *     message naming the column.
 *
 * Output:
 *
 *     written\t<label>\t<the whole file, escaped>
 *     read\t<label>\t<row>\t<the values, comma separated>
 *     metadata\t<label>\t<key=value, comma separated>
 *     comment\t<label>\t<the comments the reader kept>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: TsvTableDump
 */

import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class TsvTableDump {

    static final TableColumnCollection COLUMNS =
            new TableColumnCollection("contig", "position", "value", "note");

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("tsvtable-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# TsvTableDump: GATK's own tsv table format, from the reference");

        // Plain values, then every character that might force a quote.
        write(dir, "plain", true, new String[][] {
                {"chr1", "100", "1.0", "ordinary"},
                {"chr2", "200", "0.5", "another"},
        });
        write(dir, "quoting", true, new String[][] {
                {"chr1", "1", "1.0", "has a space"},
                {"chr1", "2", "1.0", "has\ta tab"},
                {"chr1", "3", "1.0", "has a \"quote\""},
                {"chr1", "4", "1.0", "has a \\backslash"},
                {"chr1", "5", "1.0", "has a , comma"},
                {"chr1", "6", "1.0", ""},
                {"chr1", "7", "1.0", "#starts with a comment prefix"},
        });
        // Doubles, whose spelling is String.valueOf's and not a format's.
        write(dir, "doubles", true, new String[][] {
                {"chr1", "1", String.valueOf(1.0), "integral"},
                {"chr1", "2", String.valueOf(0.1), "tenth"},
                {"chr1", "3", String.valueOf(1.0 / 3), "third"},
                {"chr1", "4", String.valueOf(1e-7), "small"},
                {"chr1", "5", String.valueOf(1e21), "large"},
                {"chr1", "6", String.valueOf(Double.NaN), "nan"},
                {"chr1", "7", String.valueOf(Double.POSITIVE_INFINITY), "infinity"},
        });
        // A header and no records at all.
        write(dir, "empty", true, new String[][] {});
        // And one where the header was never asked for.
        write(dir, "no-header", false, new String[][] {});

        for (final String label : new String[] {"plain", "quoting", "doubles", "empty"}) {
            read(dir, label, dir.resolve(label + ".table"));
            // The same file through a reader that does NOT override processCommentLine, which is
            // what fills the metadata map: overriding it without calling super loses every pair.
            readWithDefaults(dir, label, dir.resolve(label + ".table"));
        }

        // A file whose row has one column too few, and one too many.
        handWritten(dir, "short-row",
                "#sample=s1\ncontig\tposition\tvalue\tnote\nchr1\t100\t1.0\n");
        handWritten(dir, "long-row",
                "#sample=s1\ncontig\tposition\tvalue\tnote\nchr1\t100\t1.0\tnote\textra\n");
        // One whose header is missing a column the reader asks for.
        handWritten(dir, "missing-column",
                "#sample=s1\ncontig\tposition\tvalue\nchr1\t100\t1.0\n");
        // One with two metadata lines of the same key, and a plain comment between them.
        handWritten(dir, "repeated-metadata",
                "#sample=first\n# just a comment\n#sample=second\ncontig\tposition\tvalue\tnote\nchr1\t1\t1.0\tx\n");
        // One whose numeric column is empty, which only fails when it is asked for as a number.
        handWritten(dir, "empty-number",
                "#sample=s1\ncontig\tposition\tvalue\tnote\nchr1\t\t1.0\tx\n");
    }

    /** Write a table through TableWriter, with its metadata line. */
    static void write(final Path dir, final String label, final boolean withHeader,
                      final String[][] rows) throws Exception {
        final Path file = dir.resolve(label + ".table");
        try (final TableWriter<String[]> writer = new TableWriter<String[]>(file, COLUMNS) {
            @Override
            protected void composeLine(final String[] record, final DataLine dataLine) {
                dataLine.set("contig", record[0])
                        .set("position", record[1])
                        .set("value", record[2])
                        .set("note", record[3]);
            }
        }) {
            writer.writeMetadata("sample", "s1");
            if (withHeader) {
                writer.writeHeaderIfApplies();
            }
            for (final String[] row : rows) {
                writer.writeRecord(row);
            }
        }
        System.out.printf("written\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
    }

    /** A file written by hand, so a malformed one can be measured. */
    static void handWritten(final Path dir, final String label, final String text) throws Exception {
        final Path file = dir.resolve(label + ".table");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("written\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
        read(dir, label, file);
    }

    /** Read a table back, with its metadata and its comments. */
    static void read(final Path dir, final String label, final Path file) {
        final List<String> comments = new ArrayList<>();
        final List<String> rows = new ArrayList<>();
        Map<String, String> metadata = Map.of();
        try (final TableReader<String[]> reader = new TableReader<String[]>(file) {
            @Override
            protected String[] createRecord(final DataLine dataLine) {
                return new String[] {
                        dataLine.get("contig"),
                        // Asked for as a number, which is where an empty field fails.
                        String.valueOf(dataLine.getInt("position")),
                        dataLine.get("value"),
                        dataLine.get("note"),
                };
            }

            @Override
            protected void processCommentLine(final String comment, final long lineNumber) {
                comments.add(comment);
            }
        }) {
            for (final String[] record : reader) {
                rows.add(String.join(",", record));
            }
            metadata = reader.getMetadata();
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }

        for (int i = 0; i < rows.size(); i++) {
            System.out.printf("read\t%s\t%d\t%s%n", label, i, ReferenceQueryDump.escape(rows.get(i)));
        }
        final List<String> pairs = new ArrayList<>();
        metadata.forEach((key, value) -> pairs.add(key + "=" + value));
        System.out.printf("metadata\t%s\t%s%n", label, String.join(",", pairs));
        System.out.printf("comment\t%s\t%s%n", label,
                ReferenceQueryDump.escape(String.join("|", comments)));
    }

    /** The same read, with the base class's own comment handling left in place. */
    static void readWithDefaults(final Path dir, final String label, final Path file) {
        try (final TableReader<String> reader = new TableReader<String>(file) {
            @Override
            protected String createRecord(final DataLine dataLine) {
                return dataLine.get("contig");
            }
        }) {
            for (final String ignored : reader) {
                // The records do not matter here; the metadata does.
            }
            final List<String> pairs = new ArrayList<>();
            reader.getMetadata().forEach((key, value) -> pairs.add(key + "=" + value));
            System.out.printf("metadatadefault\t%s\t%s%n", label, String.join(",", pairs));
        } catch (final Exception e) {
            System.out.printf("error\t%s-default\t%s:%s%n", label, e.getClass().getName(),
                    e.getMessage());
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
