/*
 * CountVariants' count, and where it does and does not go, taken from the reference.
 *
 * The tool's own body is three lines: a counter, an increment in `apply`, and `out.print(count)`.
 * Everything there is to be identical about is therefore the traversal that decides what `apply`
 * is called on, and the output collection that decides what reaches a file.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE COUNT REACHES NO STREAM WITHOUT -O. The class documentation says "The tool prints the
 *     count to standard output (and can optionally write it to a file)", and
 *     `OptionalTextOutputArgumentCollection.print` writes NOTHING AT ALL when `output` is null. The
 *     count is the traversal's return value and the log line's, and neither is stdout;
 *   - AND THE FILE HAS NO TRAILING NEWLINE, because `onTraversalSuccess` calls `print` rather than
 *     `println`: a count of 5 is one byte, not two;
 *   - THE FILE IS TRUNCATED RATHER THAN APPENDED, since `Files.write` with no options is CREATE,
 *     TRUNCATE_EXISTING and WRITE: a count of 5 over a ten-byte file leaves one byte;
 *   - EVERY ROW COUNTS, filtered or not: there is no variant filter on this walker, so a FILTER
 *     column, a symbolic allele and a duplicate record at the same position each count once;
 *   - A RECORD IS SELECTED BY ITS WHOLE SPAN AND NOT BY ITS POSITION, because the interval
 *     traversal is a Tribble query and the codec's stop comes from the END attribute or the length
 *     of REF: a record at chr1:100 with END=400 is counted by `-L chr1:300-310`, which its POS
 *     never reaches;
 *   - AND A RECORD SPANNING TWO INTERVALS IS COUNTED ONCE, because `FeatureIntervalIterator`
 *     drops a feature that overlaps the PREVIOUS interval;
 *   - -L WITHOUT AN INDEX IS A REFUSAL, thrown by `setIntervalsForTraversal` before any record is
 *     read, and -O onto a directory is a different one, thrown after the whole traversal has run.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     count\t<label>\t<the traversal's return value>\t<its class>
 *     file\t<label>\t<present|absent>\t<byte count>\t<content, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CountVariantsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.CountVariants;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CountVariantsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                    + "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n"
                    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End of the block\">\n"
                    + "##ALT=<ID=NON_REF,Description=\"Any other allele\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=LowQD,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "##contig=<ID=chr2,length=900>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("countvariants-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CountVariantsDump: what is counted, and where the count goes");

        // Five rows of five different shapes: a plain SNP, a filtered record, a multi-allelic, a
        // symbolic block with an END, and a second copy of the first record's position.
        final Path plain = writeVcf(dir, "plain",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=1;AN=2\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tLowQD\tAC=1;AN=2\tGT\t0/1",
                "chr1\t300\t.\tA\tC,G\t50\t.\tAC=1,1;AN=2\tGT\t1/2",
                "chr1\t400\t.\tA\t<NON_REF>\t50\t.\tEND=450;AN=2\tGT\t0/0",
                "chr1\t400\t.\tA\tG\t50\t.\tAC=1;AN=2\tGT\t0/1");

        // Nothing but records the FILTER column rejects.
        final Path filteredOnly = writeVcf(dir, "filtered-only",
                "chr1\t100\t.\tA\tC\t50\tLowQD\tAC=1;AN=2\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tLowQD\tAC=1;AN=2\tGT\t0/1");

        // A header and no records at all.
        final Path empty = writeVcf(dir, "empty");

        // One record whose POS is at 100 and whose END is at 400, and one whose span is the length
        // of its REF allele. Both are wider than the position an interval would match.
        final Path spanning = writeVcf(dir, "spanning",
                "chr1\t100\t.\tA\t<NON_REF>\t50\t.\tEND=400;AN=2\tGT\t0/0",
                "chr1\t600\t.\tAAAAAAAAAA\tA\t50\t.\tAC=1;AN=2\tGT\t0/1");

        // Two contigs, so an interval on one of them is a way of counting the other's rows out.
        final Path twoContigs = writeVcf(dir, "two-contigs",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=1;AN=2\tGT\t0/1",
                "chr2\t100\t.\tA\tC\t50\t.\tAC=1;AN=2\tGT\t0/1",
                "chr2\t200\t.\tA\tC\t50\t.\tAC=1;AN=2\tGT\t0/1");

        // The same records with no index beside them, which is what -L needs.
        final Path unindexed = writeUnindexed(dir, "unindexed",
                "chr1\t100\t.\tA\tC\t50\t.\tAC=1;AN=2\tGT\t0/1");

        // What is counted, and what the count is worth without -O.
        run(dir, "plain-no-output", plain, null);
        run(dir, "plain", plain, "plain.count");
        run(dir, "filtered-only", filteredOnly, "filtered-only.count");
        run(dir, "empty", empty, "empty.count");

        // The output file is overwritten rather than appended to: this one is ten bytes before the
        // run and one byte after it.
        final Path preexisting = dir.resolve("preexisting.count");
        Files.writeString(preexisting, "9999999999", StandardCharsets.UTF_8);
        System.out.printf("file\tbefore-overwrite\tpresent\t%d\t%s%n",
                Files.size(preexisting), ReferenceQueryDump.escape("9999999999"));
        run(dir, "overwrite", plain, "preexisting.count");

        // The span, not the position: chr1:300-310 is inside the END block that starts at 100, and
        // chr1:605-606 is inside the deletion whose REF is ten bases long.
        run(dir, "span-by-end", spanning, "span-by-end.count", "-L", "chr1:300-310");
        run(dir, "span-by-ref-length", spanning, "span-by-ref-length.count", "-L", "chr1:605-606");
        run(dir, "span-missed", spanning, "span-missed.count", "-L", "chr1:500-510");

        // One record over two intervals, counted once.
        run(dir, "two-intervals-one-record", spanning, "two-intervals.count",
                "-L", "chr1:150-160", "-L", "chr1:350-360");

        // An interval that matches nothing at all, and one that selects a contig.
        run(dir, "interval-matches-nothing", plain, "nothing.count", "-L", "chr1:900-950");
        run(dir, "interval-selects-contig", twoContigs, "contig.count", "-L", "chr2");

        // The two refusals.
        run(dir, "interval-without-index", unindexed, "unindexed.count", "-L", "chr1:100-200");
        run(dir, "output-is-a-directory", plain, ".");
        run(dir, "interval-off-the-dictionary", plain, "off.count", "-L", "chr3:1-10");
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        return write(dir, label, true, records);
    }

    /** The same file with no index beside it, which is what traversal by interval needs. */
    static Path writeUnindexed(final Path dir, final String label, final String... records)
            throws Exception {
        return write(dir, label, false, records);
    }

    static Path write(final Path dir, final String label, final boolean index,
                      final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        if (index) {
            new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        }
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    /**
     * One run, reported as its return value and as whatever the output file holds afterwards.
     *
     * `output` is resolved against the dump directory, and null means the run is made with no -O
     * at all, which is the case where the count reaches nothing that can be read back.
     */
    static void run(final Path dir, final String label, final Path input, final String output,
                    final String... arguments) throws Exception {
        final List<String> all = new ArrayList<>(List.of("-V", input.toString()));
        final Path file = output == null ? null : dir.resolve(output);
        if (file != null) {
            all.addAll(List.of("-O", file.toString()));
        }
        all.addAll(List.of(arguments));

        Object result;
        try {
            result = new CountVariants().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("count\t%s\t%s\t%s%n", label, result,
                result == null ? "null" : result.getClass().getName());

        if (file == null) {
            System.out.printf("file\t%s\tno-output-argument%n", label);
            return;
        }
        if (!Files.isRegularFile(file)) {
            System.out.printf("file\t%s\tabsent%n", label);
            return;
        }
        final byte[] bytes = Files.readAllBytes(file);
        System.out.printf("file\t%s\tpresent\t%d\t%s%n", label, bytes.length,
                ReferenceQueryDump.escape(new String(bytes, StandardCharsets.UTF_8)));
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
