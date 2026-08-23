/*
 * GtfToBed's output, taken from the reference.
 *
 * A Gencode GTF reduced to one row per gene, or one per transcript. The whole tool is a map keyed
 * by gene or transcript id, a comparator, and four columns.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE OUTPUT IS NOT ZERO-BASED: the start written is the GTF's own start, so what is called a
 *     BED file is one-based and half its coordinates are off by one against every other BED;
 *   - A GENE'S INTERVAL IS WIDENED BY ITS TRANSCRIPTS, the start taken down and the end taken up,
 *     so a gene row can be wider than the gene line the file carried;
 *   - --sort-by-transcript DOES NOT SORT, IT SELECTS: it decides which of the two kinds of row is
 *     written, and the rows are sorted the same way either way;
 *   - THE ORDER IS THE DICTIONARY'S CONTIG INDEX, THEN THE START, THEN THE KEY, so two features at
 *     one position are separated by their gene or transcript id as a STRING;
 *   - A GENE ROW IS FOUR COLUMNS AND A TRANSCRIPT ROW IS FOUR, the transcript's fourth being the
 *     gene name and the transcript id joined with a COMMA rather than a fifth column;
 *   - --use-basic-transcript KEEPS ONLY TRANSCRIPTS CARRYING A tag WHOSE VALUE IS basic, and it
 *     processes such a transcript ONCE PER MATCHING TAG;
 *   - AND A TRANSCRIPT IT DROPS NEVER WIDENS ITS GENE, so the gene rows change with the flag too;
 *   - THE DICTIONARY IS REQUIRED, and its absence is a UserException naming the argument;
 *   - AND A CONTIG THE DICTIONARY DOES NOT KNOW IS A refusal from the comparator rather than from
 *     the traversal.
 *
 * Output:
 *
 *     input\t<label>=<the whole GTF, escaped>
 *     bed\t<label>=<the whole output, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GtfToBedDump
 */

import org.broadinstitute.hellbender.tools.walkers.conversion.GtfToBed;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GtfToBedDump {

    /** One GTF line, whose attributes are the ones the codec needs and the tags asked for. */
    static String gene(final String contig, final int start, final int end, final String geneId,
                       final String geneName) {
        return contig + "\tHAVANA\tgene\t" + start + "\t" + end + "\t.\t+\t.\t"
                + "gene_id \"" + geneId + "\"; gene_type \"protein_coding\"; gene_name \"" + geneName
                + "\"; level 2; havana_gene \"OTTHUMG00000000001.1\";\n";
    }

    static String transcript(final String contig, final int start, final int end,
                             final String geneId, final String geneName, final String transcriptId,
                             final String... tags) {
        final StringBuilder line = new StringBuilder(contig + "\tHAVANA\ttranscript\t" + start + "\t"
                + end + "\t.\t+\t.\t"
                + "gene_id \"" + geneId + "\"; transcript_id \"" + transcriptId + "\"; "
                + "gene_type \"protein_coding\"; gene_name \"" + geneName + "\"; "
                + "transcript_type \"protein_coding\"; transcript_name \"" + transcriptId + "\"; "
                + "level 2;");
        for (final String tag : tags) {
            line.append(" tag \"").append(tag).append("\";");
        }
        return line.append(" havana_gene \"OTTHUMG00000000001.1\";\n").toString();
    }

    static String exon(final String contig, final int start, final int end, final String geneId,
                       final String geneName, final String transcriptId) {
        return contig + "\tHAVANA\texon\t" + start + "\t" + end + "\t.\t+\t.\t"
                + "gene_id \"" + geneId + "\"; transcript_id \"" + transcriptId + "\"; "
                + "gene_type \"protein_coding\"; gene_name \"" + geneName + "\"; "
                + "transcript_type \"protein_coding\"; transcript_name \"" + transcriptId + "\"; "
                + "exon_number 1; exon_id \"" + transcriptId + ".1\"; level 2;\n";
    }

    /**
     * Two contigs and four genes.
     *
     * The first gene's transcripts reach past it at both ends, so its row is wider than its line.
     * The second gene has one basic transcript and one that is not. The third and fourth start at
     * the same position, so their order is settled by their ids as strings.
     */
    static String gtf() {
        return gene("chr1", 100, 200, "GENE_B.1", "beta")
                + transcript("chr1", 50, 250, "GENE_B.1", "beta", "TX_B1.1", "basic")
                + exon("chr1", 50, 250, "GENE_B.1", "beta", "TX_B1.1")
                + transcript("chr1", 120, 180, "GENE_B.1", "beta", "TX_B2.1")
                + exon("chr1", 120, 180, "GENE_B.1", "beta", "TX_B2.1")
                + gene("chr1", 300, 400, "GENE_A.1", "alpha")
                + transcript("chr1", 300, 400, "GENE_A.1", "alpha", "TX_A1.1", "basic", "basic")
                + exon("chr1", 300, 400, "GENE_A.1", "alpha", "TX_A1.1")
                + gene("chr1", 300, 400, "GENE_C.1", "gamma")
                + transcript("chr1", 300, 500, "GENE_C.1", "gamma", "TX_C1.1")
                + exon("chr1", 300, 500, "GENE_C.1", "gamma", "TX_C1.1")
                + gene("chr2", 10, 20, "GENE_D.1", "delta")
                + transcript("chr2", 10, 20, "GENE_D.1", "delta", "TX_D1.1")
                + exon("chr2", 10, 20, "GENE_D.1", "delta", "TX_D1.1");
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("gtf-to-bed-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GtfToBedDump: a Gencode GTF reduced to one row per gene or transcript");

        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "gtf", List.of("chr1", "chr2"));
        final Path gtf = dir.resolve("annotation.gtf");
        Files.writeString(gtf, gtf(), StandardCharsets.UTF_8);
        System.out.printf("input\tannotation=%s%n", ReferenceQueryDump.escape(gtf()));

        run(dir, "genes", gtf, dict);
        run(dir, "transcripts", gtf, dict, "--sort-by-transcript", "true");
        run(dir, "genes-basic-only", gtf, dict, "--use-basic-transcript", "true");
        run(dir, "transcripts-basic-only", gtf, dict,
                "--sort-by-transcript", "true", "--use-basic-transcript", "true");
        // No dictionary at all.
        run(dir, "no-dictionary", gtf, null);
        // A contig the dictionary does not know.
        final Path narrow = MultiFeatureWalkerDump.writeDictionary(dir, "chr1only", List.of("chr1"));
        run(dir, "unknown-contig", gtf, narrow);
    }

    static void run(final Path dir, final String label, final Path gtf, final Path dictionary,
                    final String... extra) throws Exception {
        final Path out = dir.resolve(label + ".bed");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--gtf-path", gtf.toString(), "--output", out.toString()));
        if (dictionary != null) {
            argv.addAll(Arrays.asList("--sequence-dictionary", dictionary.toString()));
        }
        argv.addAll(Arrays.asList(extra));
        try {
            new GtfToBed().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("bed\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
