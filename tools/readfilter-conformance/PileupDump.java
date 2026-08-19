/*
 * Pileup, taken from the reference.
 *
 * A samtools-style pileup line per locus, and three optional columns that each have their own
 * shape. Small, and every part of it is a format rather than a computation.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE DEFAULT FILTERS ARE THE WALKER'S PLUS THREE. `LocusWalker` contributes wellformed and
 *     mapped; this tool adds NotDuplicate, PassesVendorQualityCheck and NotSecondaryAlignment, so a
 *     duplicate that reaches CountReads never reaches this one;
 *   - THE PILEUP IS FILTERED OF DELETIONS BEFORE ANYTHING IS PRINTED, so a deletion contributes to
 *     neither the base string nor the insert lengths -- and the VERBOSE column then reports a
 *     deletion count of zero for every locus, because it counts them in the already-filtered
 *     pileup;
 *   - WITHOUT A REFERENCE THE BASE IS 'N', which changes every base in the string that would have
 *     matched the reference: `.` and `,` become the base itself;
 *   - THE FEATURES COLUMN IS ALWAYS PRINTED, even when it is empty, so every line ends with a
 *     trailing space when there is no metadata;
 *   - THE INSERT LENGTHS ARE ONE PER READ, comma separated, in pileup order, and a read with no
 *     mate has a fragment length of zero;
 *   - AND THE VERBOSE COLUMN IS `name@offset@length@mappingQuality` per read, joined by commas,
 *     after the deletion count and a space.
 *
 * Output:
 *
 *     pileup\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PileupDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.qc.Pileup;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PileupDump {

    /** Metadata features: two records, one of which overlaps a locus the pileup prints. */
    static final String METADATA =
            "##fileformat=VCFv4.2\n"
            + "##contig=<ID=chr1,length=200>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t12\t.\tG\tC\t50\tPASS\t.\n"
            + "chr1\t13\t.\tT\tA\t50\tPASS\t.\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("pileup-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("reads.bam");
        ReadWalkerDump.buildFixture(bam.toFile());

        final Path metadata = dir.resolve("metadata.vcf");
        Files.write(metadata, METADATA.getBytes());
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", metadata.toString()});

        System.out.println("# PileupDump: a samtools-style pileup and its optional columns");

        // The read at chr1:10, which is ten bases over ten loci.
        run("plain", dir, fasta, bam, "-L", "chr1:10-19");
        run("verbose", dir, fasta, bam, "-L", "chr1:10-19", "--show-verbose", "true");
        run("insert-length", dir, fasta, bam, "-L", "chr1:10-19",
                "--output-insert-length", "true");
        run("both", dir, fasta, bam, "-L", "chr1:10-19", "--show-verbose", "true",
                "--output-insert-length", "true");
        // Metadata, whose features are printed in brackets at the loci they overlap.
        run("metadata", dir, fasta, bam, "-L", "chr1:10-19", "--metadata", metadata.toString());
        // No reference at all, so every reference base is N.
        run("no-reference", dir, null, bam, "-L", "chr1:10-19");
        // The deletion read, whose deleted bases are filtered out of the pileup entirely.
        run("deletion", dir, fasta, bam, "-L", "chr1:140-160", "--show-verbose", "true");
        // The duplicate at 170, which this tool's own filters remove.
        run("duplicate", dir, fasta, bam, "-L", "chr1:170-179");
        // The soft-masked stretch, where the reference base is upper-cased before it is compared.
        run("masked", dir, fasta, bam, "-L", "chr1:65-74");
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("pileup-" + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-O", out.toString()));
        if (fasta != null) {
            argv.add("-R");
            argv.add(fasta.toString());
        }
        argv.addAll(Arrays.asList(extra));
        try {
            new Pileup().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("pileup\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
