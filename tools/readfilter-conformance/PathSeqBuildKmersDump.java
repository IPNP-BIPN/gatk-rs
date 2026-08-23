/*
 * PathSeqBuildKmers's output, taken from the reference.
 *
 * A host reference turned into the k-mer set PathSeq subtracts against. The output is a Kryo
 * serialisation, so what is printed here is the set read back: every k-mer as the long it is stored
 * as, in order, with the bases that long spells.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - A K-MER IS TWO BITS PER BASE, A=0 C=1 G=2 T=3, the first base in the HIGH bits, and the long
 *     is the whole of it;
 *   - THE SET IS CANONICAL: a k-mer whose middle base is G or T is replaced by its reverse
 *     complement, which is why the size must be odd, and an even size is refused by canonical()
 *     itself rather than by the argument parser;
 *   - MASKING IS AN AND WITH A LONG BUILT FROM THE MASKED POSITIONS, and the positions are counted
 *     from the START of the k-mer, so masking position 0 clears the first base's bits;
 *   - MASKING HAPPENS AFTER CANONICALISATION, so which bases a mask clears depends on which strand
 *     won;
 *   - A MASK INDEX OUTSIDE THE K-MER IS REFUSED with the index it was given;
 *   - ANY BASE THAT IS NOT ACGT RESTARTS THE COUNT, the valid-base counter being set to -1 and
 *     then incremented, so a k-mer is emitted only after k good bases in a row and an N costs a
 *     whole window;
 *   - LOWER CASE COUNTS, a and A being the same base;
 *   - --kmer-spacing IS APPLIED BY REWINDING THE VALID-BASE COUNTER to kSize - spacing, so a
 *     spacing of one is every position and a spacing of k is non-overlapping windows;
 *   - THE SET HOLDS EACH DISTINCT LONG ONCE, however many places it came from;
 *   - A K LONGER THAN EVERY CONTIG PRODUCES NOTHING AND IS THEN REFUSED BY THE SET ITSELF, so no
 *     empty set is ever written;
 *   - AND A BLOOM FILTER ANSWERS THE SAME QUESTIONS with a theoretical false positive probability
 *     the tool computes and reports.
 *
 * Output:
 *
 *     fixture\t<name>=<the whole file, escaped>
 *     count\t<label>=<number of distinct k-mers>
 *     kmer\t<label>\t<long>\t<bases>
 *     bloom\t<label>\t<key>=<value>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PathSeqBuildKmersDump
 */

import org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerBloomFilter;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerCollection;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerSet;
import org.broadinstitute.hellbender.tools.spark.pathseq.PSKmerUtils;
import org.broadinstitute.hellbender.tools.spark.pathseq.PathSeqBuildKmers;
import org.broadinstitute.hellbender.tools.spark.utils.LongIterator;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVKmerShort;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PathSeqBuildKmersDump {

    /**
     * Two contigs: one plain, one carrying an N and a lower-case run, both short enough that every
     * k-mer can be printed.
     */
    static final String CONTIG_ONE = "ACGTTGCAAGGCTTACCATGG";
    static final String CONTIG_TWO = "ACGTNACGTacgtTTTTT";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("pathseq-kmers-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# PathSeqBuildKmersDump: a host reference turned into a k-mer set");

        final Path fasta = writeReference(dir);

        run(dir, "k5", fasta, 5, "", 1, 0);
        run(dir, "k5-mask-middle", fasta, 5, "2", 1, 0);
        run(dir, "k5-mask-first-and-last", fasta, 5, "0,4", 1, 0);
        run(dir, "k5-spacing-three", fasta, 5, "", 3, 0);
        run(dir, "k5-spacing-five", fasta, 5, "", 5, 0);
        run(dir, "k7", fasta, 7, "", 1, 0);
        run(dir, "k5-bloom", fasta, 5, "", 1, 0.01);
        // An even k, which canonical() refuses.
        run(dir, "even-k", fasta, 4, "", 1, 0);
        // A mask index the k-mer does not have.
        run(dir, "mask-out-of-range", fasta, 5, "9", 1, 0);
        // A k longer than either contig, which produces nothing at all.
        run(dir, "k-longer-than-reference", fasta, 31, "", 1, 0);
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("host.fasta");
        Files.writeString(fasta, ">one\n" + CONTIG_ONE + "\n>two\n" + CONTIG_TWO + "\n",
                StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("host.dict")});
        System.out.printf("fixture\tone=%s%n", CONTIG_ONE);
        System.out.printf("fixture\ttwo=%s%n", CONTIG_TWO);
        return fasta;
    }

    static void run(final Path dir, final String label, final Path fasta, final int kmerSize,
                    final String mask, final int spacing, final double bloomFpp) throws Exception {
        final Path out = dir.resolve(label + (bloomFpp > 0 ? ".bfi" : ".hss"));
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "--reference", fasta.toString(),
                "--output", out.toString(),
                "--kmer-size", Integer.toString(kmerSize),
                "--kmer-mask", mask,
                "--kmer-spacing", Integer.toString(spacing),
                "--bloom-false-positive-probability", Double.toString(bloomFpp)));
        try {
            new PathSeqBuildKmers().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        final PSKmerCollection collection = PSKmerUtils.readKmerFilter(out.toString());
        if (collection instanceof PSKmerBloomFilter) {
            final PSKmerBloomFilter filter = (PSKmerBloomFilter) collection;
            System.out.printf("bloom\t%s\tkmer-size=%d%n", label, filter.kmerSize());
            System.out.printf("bloom\t%s\tfalse-positive-probability=%s%n", label,
                    Double.toString(filter.getFalsePositiveProbability()));
            // Every k-mer the plain set holds must be in the filter, and a few that are not in the
            // reference at all are asked as well.
            for (final long value : sortedKmers(dir, label, fasta, filter.kmerSize())) {
                System.out.printf("bloom\t%s\tcontains-%d=%b%n", label, value,
                        filter.contains(new SVKmerShort(value)));
            }
            return;
        }
        final PSKmerSet set = (PSKmerSet) collection;
        final List<Long> values = new ArrayList<>();
        final LongIterator iterator = set.iterator();
        while (iterator.hasNext()) {
            values.add(iterator.next());
        }
        values.sort(Long::compareTo);
        System.out.printf("count\t%s=%d%n", label, values.size());
        for (final long value : values) {
            System.out.printf("kmer\t%s\t%d\t%s%n", label, value,
                    new SVKmerShort(value).toString(kmerSize));
        }
    }

    /** The k-mers of the plain run of the same size, which the Bloom filter is asked about. */
    static List<Long> sortedKmers(final Path dir, final String label, final Path fasta,
                                  final int kmerSize) throws Exception {
        final Path plain = dir.resolve("k" + kmerSize + ".hss");
        final List<Long> values = new ArrayList<>();
        if (!Files.exists(plain)) {
            return values;
        }
        final PSKmerSet set = (PSKmerSet) PSKmerUtils.readKmerFilter(plain.toString());
        final LongIterator iterator = set.iterator();
        while (iterator.hasNext()) {
            values.add(iterator.next());
        }
        values.sort(Long::compareTo);
        return values;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
