/*
 * What a -L argument turns out to be, taken from the reference.
 *
 * IntervalWalkerDump measures what happens once the intervals exist. This measures the step before
 * it: IntervalUtils.parseIntervalArguments deciding whether one -L string is a Feature file, a GATK
 * interval file, a file it refuses, or a literal interval.
 *
 * The order of those tests is the whole behaviour, and it is observable:
 *
 *   - a Feature file is recognised first, by asking every registered codec whether it can decode
 *     the file, so a `.list` whose contents happen to be a Feature file goes down the Feature path
 *     and never reaches the interval-file reader;
 *   - the interval-file test is by *extension only*, lower-cased, and throws when the extension
 *     matches but the file is absent. It cannot test the contents, because a contig name may
 *     contain a period and the reference says so;
 *   - an existing file that is neither is an error naming both possibilities;
 *   - only a non-existent argument with neither extension is parsed as a literal interval, which
 *     is why a typo in a filename surfaces as "contig not in the dictionary".
 *
 * `.interval_list` and `.bed` are here because they are Feature files rather than interval files:
 * they never reach the interval-file branch at all, which is not what their names suggest.
 *
 * Output:
 *
 *     case\t<label>\t<ok|E:class>\t<n>\t<interval|...>
 *
 * Usage: IntervalFileDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class IntervalFileDump {

    static final int CONTIG_LENGTH = 200;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("intervalfile-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                new SAMSequenceRecord("chr2", CONTIG_LENGTH))));
        final GenomeLocParser parser = new GenomeLocParser(header.getSequenceDictionary());

        System.out.println("# IntervalFileDump: what one -L argument resolves to");

        // An ordinary GATK interval file, and the same contents under the other extension.
        final Path list = write(dir, "a.list", "chr1:1-10\nchr1:50-60\nchr2\n");
        final Path intervals = write(dir, "b.intervals", "chr1:1-10\nchr1:50-60\nchr2\n");
        probe(parser, "list", list.toString());
        probe(parser, "intervals", intervals.toString());

        // Whitespace and blank lines: each line is trimmed and empty ones are skipped.
        probe(parser, "whitespace",
                write(dir, "ws.list", "\n  chr1:1-10  \n\n\tchr2:5-6\n\n").toString());

        // A file whose only content is blank lines holds no intervals, which is an error rather
        // than an empty traversal.
        probe(parser, "blank-only", write(dir, "blank.list", "\n\n   \n").toString());

        // A completely empty file, which is a different path to the same claim.
        probe(parser, "empty", write(dir, "empty.list", "").toString());

        // The extension test is lower-cased, so an upper-case extension is still an interval file.
        probe(parser, "uppercase-extension",
                write(dir, "c.LIST", "chr1:1-5\n").toString());

        // A missing file with an interval-file extension throws, rather than being parsed as a
        // literal interval whose contig happens to end in `.list`.
        probe(parser, "missing-list", dir.resolve("absent.list").toString());

        // A file that exists with no recognised extension: neither Features nor intervals.
        probe(parser, "unknown-extension",
                write(dir, "d.txt", "chr1:1-10\n").toString());

        // A file that does not exist and has no recognised extension is a literal interval, so the
        // failure is about the contig rather than about the file.
        probe(parser, "missing-unknown-extension", dir.resolve("absent.txt").toString());

        // The removed multi-interval syntax.
        probe(parser, "semicolon", "chr1:1-10;chr2:1-10");

        // Feature files. These never reach the interval-file branch: `.interval_list` and `.bed`
        // are decoded by a codec, and their intervals come from the Features they contain.
        probe(parser, "picard-interval-list", writeIntervalList(dir, header).toString());
        probe(parser, "bed", write(dir, "f.bed", "chr1\t0\t10\nchr2\t4\t6\n").toString());

        // A `.list` file whose contents are a BED body. The extension says interval file, the
        // codec says Feature file, and the codec is asked first.
        probe(parser, "bed-contents-list-extension",
                write(dir, "g.list", "chr1\t0\t10\nchr2\t4\t6\n").toString());

        // Plain literal intervals, for the row that anchors the rest.
        probe(parser, "literal", "chr1:1-10");
        probe(parser, "literal-whole-contig", "chr2");
    }

    static void probe(final GenomeLocParser parser, final String label, final String argument) {
        String outcome = "ok";
        List<GenomeLoc> locs = List.of();
        try {
            locs = IntervalUtils.parseIntervalArguments(parser, argument);
        } catch (final Exception e) {
            // The class only: the message quotes the run's absolute paths.
            outcome = "E:" + e.getClass().getName();
        }
        final StringBuilder text = new StringBuilder();
        for (final GenomeLoc loc : locs) {
            if (text.length() > 0) {
                text.append('|');
            }
            text.append(loc.toString());
        }
        System.out.printf("case\t%s\t%s\t%d\t%s%n", label, outcome, locs.size(), text);
    }

    static Path write(final Path dir, final String name, final String contents) throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, contents.getBytes());
        return path;
    }

    /** A Picard-style `.interval_list`, which needs the SAM header its codec validates against. */
    static Path writeIntervalList(final Path dir, final SAMFileHeader header) throws Exception {
        final htsjdk.samtools.util.IntervalList list = new htsjdk.samtools.util.IntervalList(header);
        list.add(new htsjdk.samtools.util.Interval("chr1", 1, 10));
        list.add(new htsjdk.samtools.util.Interval("chr2", 5, 6));
        final Path path = dir.resolve("e.interval_list");
        list.write(path.toFile());
        return path;
    }
}
