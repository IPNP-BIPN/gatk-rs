/*
 * MergeAnnotatedRegions, taken from the reference.
 *
 * A segment file in, a segment file out, with touching or overlapping regions merged and their
 * annotations reconciled. The file format is the annotated-interval collection's, which four tools
 * share, so this run pins the reading and the writing as well as the merging.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE LOCATABLE COLUMNS ARE FOUND BY NAME FROM A FIXED LIST. `CONTIG`, `chrom`, `Chromosome`
 *     and nine others all name the contig; `Position` names BOTH the start and the end, being in
 *     both lists, so a file with a single `Position` column parses as one-base regions;
 *   - THE OUTPUT IS ALWAYS TAB SEPARATED AND ALWAYS CARRIES A SAM HEADER, whatever the input had:
 *     a file read with no `@` lines comes back with a bare `@HD`;
 *   - COMMENT LINES BECOME SAM HEADER COMMENTS, so a `#` preamble comes out as `@CO` lines;
 *   - THE ANNOTATION COLUMNS ARE SORTED ALPHABETICALLY, not kept in the input's order, and the
 *     three locatable columns are renamed to `CONTIG`, `START` and `END`;
 *   - MERGING IS ONE PASS WITH A PEEK OVER A SORTED LIST, and the merged region is what the next
 *     comparison uses, so a chain merges into one row however long it is. Sorting first is what
 *     makes that safe: a region written out of order between two that overlap it is pulled into
 *     the chain rather than left behind, which the `skipped-overlap` row shows;
 *   - ABUTTING REGIONS ARE NOT MERGED. `IntervalUtils.overlaps` is a real overlap, so `1-100` and
 *     `101-200` come out as two rows despite the tool's own summary saying "touching";
 *   - A CONFLICTING ANNOTATION IS SPLIT ON THE SEPARATOR, DEDUPLICATED, SORTED AND REJOINED, so
 *     merging three regions with values `b`, `a` and `b` gives `a__b` and not `b__a__b`;
 *   - AN ANNOTATION MISSING FROM ONE SIDE IS PASSED THROUGH UNCHANGED, not merged with an empty
 *     string;
 *   - THE ROWS ARE SORTED BY THE SEQUENCE DICTIONARY BEFORE ANY OF THIS, so an input in the wrong
 *     contig order is reordered rather than refused;
 *   - A CONTIG THE REFERENCE DOES NOT CARRY IS PASSED THROUGH, not refused: the sort tolerates it
 *     and the row is written back out as it came in;
 *   - AND A FILE OF NO ROWS THROWS. The collection reads its annotations off the first record, so
 *     a file with a column line and nothing under it dies with an `IndexOutOfBoundsException`
 *     rather than writing an empty result;
 *   - THE OUTPUT HEADER CARRIES THREE `@CO` LINES naming the locatable columns it renamed
 *     (`_ContigHeader`, `_StartHeader`, `_EndHeader`), after any comments the input had.
 *
 * Output:
 *
 *     merged\t<label>=<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MergeAnnotatedRegionsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.copynumber.utils.MergeAnnotatedRegions;

import java.nio.file.Files;
import java.nio.file.Path;

public class MergeAnnotatedRegionsDump {

    /** A plain segment file with three annotations, in the order a caller happened to write them. */
    static final String PLAIN =
            "CONTIG\tSTART\tEND\tname\tvalue\tcall\n"
            + "chr1\t1\t100\tone\t0.5\t+\n"
            + "chr1\t50\t150\ttwo\t0.5\t+\n"
            + "chr1\t300\t400\tthree\t1.5\t-\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("merge-annotated-regions-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, PreprocessIntervalsDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# MergeAnnotatedRegionsDump: the segments a merge leaves behind");

        run("plain", dir, fasta, PLAIN);
        // Abutting regions, which are NOT merged.
        run("abutting", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\ta\n"
                + "chr1\t101\t200\tb\n");
        // A chain of three overlaps, merged left to right, whose annotation values repeat.
        run("chain", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\tb\n"
                + "chr1\t50\t150\ta\n"
                + "chr1\t120\t200\tb\n");
        // A region written between two that overlap it, which the sort pulls into the chain.
        run("skipped-overlap", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\ta\n"
                + "chr1\t120\t130\tb\n"
                + "chr1\t90\t200\tc\n");
        // An annotation present in one row and absent from the other.
        run("missing-annotation", dir, fasta,
                "CONTIG\tSTART\tEND\tname\tvalue\n"
                + "chr1\t1\t100\ta\t\n"
                + "chr1\t50\t150\tb\t7\n");
        // Values that already carry the separator, which are split before they are rejoined.
        run("separator-in-value", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\tb__c\n"
                + "chr1\t50\t150\ta__b\n");
        // Rows out of dictionary order, which are sorted first.
        run("unsorted", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chr2\t1\t100\tz\n"
                + "chr1\t200\t240\ty\n"
                + "chr1\t1\t100\tx\n");
        // Other column names, including the `Position` that names both ends.
        run("other-column-names", dir, fasta,
                "Chromosome\tStart_Position\tEnd_Position\tname\n"
                + "chr1\t1\t100\ta\n"
                + "chr1\t50\t150\tb\n");
        run("position-column", dir, fasta,
                "chrom\tPosition\tname\n"
                + "chr1\t10\ta\n"
                + "chr1\t10\tb\n"
                + "chr1\t20\tc\n");
        // A comment preamble, which becomes SAM header comments.
        run("comments", dir, fasta,
                "#a note\n"
                + "#another\n"
                + "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\ta\n");
        // A SAM header on the input, whose sequence lines are replaced by the reference's.
        run("sam-header", dir, fasta,
                "@HD\tVN:1.6\n"
                + "@SQ\tSN:chr1\tLN:240\n"
                + "CONTIG\tSTART\tEND\tname\n"
                + "chr1\t1\t100\ta\n");
        // No rows at all, only a column line.
        run("no-rows", dir, fasta, "CONTIG\tSTART\tEND\tname\n");
        // A contig the reference does not carry, which is passed through rather than refused.
        run("unknown-contig", dir, fasta,
                "CONTIG\tSTART\tEND\tname\n"
                + "chrX\t1\t100\ta\n");
        // A file with no locatable columns at all.
        run("no-locatable-columns", dir, fasta, "name\tvalue\na\t1\n");
    }

    static void run(final String label, final Path dir, final Path fasta, final String input)
            throws Exception {
        final Path in = dir.resolve(label + ".seg");
        Files.write(in, input.getBytes());
        final Path out = dir.resolve("merged-" + label + ".seg");
        try {
            new MergeAnnotatedRegions().instanceMain(new String[] {
                    "-R", fasta.toString(), "--segments", in.toString(), "-O", out.toString()});
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("merged\t%s=%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
