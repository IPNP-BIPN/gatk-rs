/*
 * ReferenceDataSource's answers, taken from the reference.
 *
 * A walker sees the reference through ReferenceDataSource.queryAndPrefetch(contig, start, stop),
 * which is htsjdk's indexed FASTA reader underneath. Everything that compares a read to the
 * reference (every annotation, every caller) reads its bases through this path, so a query that
 * returns different bytes, or that answers where the reference throws, is wrong everywhere at once.
 *
 * The fixture travels in the golden: the FASTA and its .fai are printed as escaped text, so the
 * port reads exactly the bytes the reference read rather than a reconstruction of them. The FASTA
 * is deliberately awkward: mixed case (soft-masked regions are lower case and stay lower case), a
 * line width that does not divide the sequence length, a contig whose length *is* a multiple of the
 * line width, an N run, and IUPAC codes.
 *
 * Output:
 *
 *     fasta\t<escaped FASTA text>
 *     fai\t<escaped .fai text>
 *     query\t<contig>\t<start>\t<stop>\t<bases or E>
 *
 * Usage: ReferenceQueryDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ReferenceQueryDump {

    /** Twelve bases a line, which divides neither contig evenly on purpose. */
    static final String FASTA =
            ">chr1 first contig\n"
            + "ACGTACGTACGT\n"
            + "acgtNNNNacgt\n"
            + "ACGTRYKMSWBD\n"
            + "HVNACGT\n"
            + ">chr2\n"
            + "TTTTGGGGCCCC\n"
            + "AAAATTTTGGGG\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("refquery");
        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        final Path fai = dir.resolve("ref.fasta.fai");
        // GATK refuses a reference with no sequence dictionary, so the harness makes one the same
        // way a user would. It carries no information the .fai does not, and does not travel.
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# ReferenceQueryDump: ReferenceDataSource's bases, from the reference");
        System.out.printf("fasta\t%s%n", escape(new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n", escape(new String(Files.readAllBytes(fai))));

        try (final ReferenceDataSource source = new ReferenceFileSource(fasta)) {
            final List<int[]> spans = new ArrayList<>();
            // Line boundaries, contig edges, single bases, inverted and out-of-range spans.
            for (final int[] span : new int[][] {
                    {1, 1}, {1, 12}, {1, 13}, {12, 13}, {13, 24}, {5, 20}, {24, 25},
                    {31, 43}, {43, 43}, {1, 43}, {40, 50}, {44, 44}, {0, 5}, {5, 4}, {-3, 2}}) {
                spans.add(span);
            }
            for (final String contig : new String[] {"chr1", "chr2", "chr3"}) {
                for (final int[] span : spans) {
                    final int start = span[0];
                    final int stop = span[1];
                    String bases;
                    try {
                        bases = new String(source
                                .queryAndPrefetch(new SimpleInterval(contig, start, stop))
                                .getBases());
                    } catch (final Exception | AssertionError e) {
                        bases = "E";
                    }
                    System.out.printf("query\t%s\t%d\t%d\t%s%n", contig, start, stop, bases);
                }
            }
        }
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}
