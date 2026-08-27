/*
 * GroundTruthReadsBuilder's CSV, taken from the reference.
 *
 * Every read scored against the haplotype its two ancestral references give it, written as one
 * CSV row. The flow-based scoring engine is not what this measures: what it measures is the
 * translation from the aligned contig to the two ancestral ones, which reads survive it, and what
 * the row carries.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE TOOL REFUSES ANY BUT A FLOW-BASED ENGINE outright, so every run has to name one;
 *   - IT ALSO NEEDS THE ALIGNED REFERENCE, not only the two ancestral ones: every read is scored
 *     against it too, and without one the haplotype is null and the engine refuses;
 *   - THE TRANSLATION IS A TABLE OF OFFSETS, one CSV per ancestor and contig, whose first line
 *     is IGNORED and whose rows are a position and the offset that applies from it on;
 *   - A POSITION BETWEEN TWO ROWS TAKES THE EARLIER ROW'S OFFSET, the search falling back on the
 *     insertion point less two;
 *   - THE TRANSLATED CONTIG IS THE READ'S OWN NAME WITH THE ANCESTOR APPENDED, so the two
 *     reference files have to be named for it;
 *   - A READ WHOSE TRANSLATED END IS NOT PAST ITS TRANSLATED START IS SKIPPED AND COUNTED, not
 *     refused: the exception is caught and the traversal carries on;
 *   - AN HAPLOTYPE IDENTICAL TO THE REFERENCE IS NOT RESCORED, taking the reference's own score,
 *     so a read away from either ancestral difference reports the two as equal;
 *   - A READ BELOW --min-mq IS DROPPED;
 *   - THE SOFT-CLIP FILTER LOOKS AT THE CLIP AT THE END OF THE READ rather than at either end,
 *     so a read clipped only at its front is kept whatever the argument says;
 *   - --min-haplotype-score-delta DROPS A READ WHOSE TWO HAPLOTYPES ARE TOO FAR APART, which is
 *     every read away from both ancestral differences;
 *   - --output-flow-length FIXES THE LENGTH OF THE TWO HAPLOTYPE KEYS IN BOTH DIRECTIONS, padding
 *     a short one and TRUNCATING a long one, and it does not touch the read's own key at all;
 *   - AND THE CSV IS QUOTED: the flow keys hold commas of their own, so a reader that splits on
 *     the comma alone reads the columns out of step.
 *
 * Output:
 *
 *     fasta\t<label>=<that reference's contig and length>
 *     csv\t<label>=<that translation table, escaped>
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     out\t<label>=<the whole output csv, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GroundTruthReadsBuilderDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthReadsBuilder;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GroundTruthReadsBuilderDump {

    static final int LENGTH = 2400;
    static final String FLOW_ORDER = "TGCA";

    static byte referenceBase(final int position) {
        return (byte) "TGCA".charAt((position - 1) % 4);
    }

    static String referenceBases(final int start, final int length) {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < length; i++) {
            bases.append((char) referenceBase(start + i));
        }
        return bases.toString();
    }

    /** One ancestral reference: the shared sequence with a handful of bases changed. */
    static String ancestralBases(final int[] substitutions) {
        final char[] bases = referenceBases(1, LENGTH).toCharArray();
        for (final int offset : substitutions) {
            bases[offset] = bases[offset] == 'T' ? 'A' : 'T';
        }
        return new String(bases);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("ground-truth-reads-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GroundTruthReadsBuilderDump: every read scored against the "
                + "haplotype its two ancestral references give it");

        // The two ancestors differ from one another at two positions, so the haplotypes differ.
        // The aligned reference too: the tool scores every read against it as well as against
        // the two ancestors, and without one the haplotype is null and the engine refuses.
        final Path reference = writeReference(dir, "reference", "chr1", referenceBases(1, LENGTH));
        final Path maternal = writeReference(dir, "maternal", "chr1_maternal",
                ancestralBases(new int[] {1049}));
        final Path paternal = writeReference(dir, "paternal", "chr1_paternal",
                ancestralBases(new int[] {1199}));

        // The translation tables: an identity for the maternal one and a shift for the paternal.
        final String maternalCsv = String.join("\n",
                "pos,offset", "1,0", "");
        final String paternalCsv = String.join("\n",
                "pos,offset", "1,0", "1500,10", "");
        write(dir, "maternal.chr1.csv", maternalCsv);
        write(dir, "paternal.chr1.csv", paternalCsv);
        System.out.printf("csv\tmaternal=%s%n", ReferenceQueryDump.escape(maternalCsv));
        System.out.printf("csv\tpaternal=%s%n", ReferenceQueryDump.escape(paternalCsv));

        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam);
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(bam)));

        run(dir, reference, "default", bam, maternal, paternal, List.of());
        // A mapping-quality floor, which drops the badly-mapped read.
        run(dir, reference, "min-mq-thirty", bam, maternal, paternal, List.of("--min-mq", "30"));
        // The soft-clipped reads kept rather than dropped.
        run(dir, reference, "keep-softclipped", bam, maternal, paternal,
                List.of("--discard-non-polyt-softclipped-reads", "false"));
        // A sequence added to the haplotype at each end.
        run(dir, reference, "prepend-append", bam, maternal, paternal,
                List.of("--prepend-sequence", "TTTT", "--append-sequence", "CCCC"));
        // A flow length the keys have to reach.
        run(dir, reference, "flow-length-ten", bam, maternal, paternal,
                List.of("--output-flow-length", "10"));
        run(dir, reference, "flow-length-large", bam, maternal, paternal,
                List.of("--output-flow-length", "400"));
        // The two score filters.
        run(dir, reference, "min-score", bam, maternal, paternal,
                List.of("--min-haplotype-score", "-1.0"));
        run(dir, reference, "min-score-delta", bam, maternal, paternal,
                List.of("--min-haplotype-score-delta", "1.0"));
        // A cap on how many reads reach the output.
        run(dir, reference, "max-two-reads", bam, maternal, paternal, List.of("--max-output-reads", "2"));
        // A read whose translated span collapses, which is a refusal rather than a filter.
        final String collapsing = String.join("\n", "pos,offset", "1,0", "1010,-200", "");
        write(dir, "paternal.chr1.csv", collapsing);
        System.out.printf("csv\tcollapsing=%s%n", ReferenceQueryDump.escape(collapsing));
        run(dir, reference, "collapsing-translation", bam, maternal, paternal, List.of());
    }

    /**
     * The reads. Six over the shared contig.
     *
     * One is badly mapped, one is soft-clipped with a poly-T clip, one is soft-clipped with
     * something else, and one sits past the paternal table's second row so its offset differs.
     */
    static void writeReads(final Path bam) {
        final SAMFileHeader header = readHeader();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            add(writer, header, "r-plain", 1000, "100M", 60, "");
            add(writer, header, "r-second", 1100, "100M", 60, "");
            add(writer, header, "r-low-mq", 1200, "100M", 10, "");
            add(writer, header, "r-polyt-clip", 1300, "4S96M", 60, "TTTT");
            add(writer, header, "r-other-clip", 1400, "4S96M", 60, "ACGT");
            // Past the paternal table's second row, so its paternal offset is ten.
            add(writer, header, "r-shifted", 1600, "100M", 60, "");
        }
    }

    static void add(final SAMFileWriter writer, final SAMFileHeader header, final String name,
                    final int start, final String cigar, final int mappingQuality,
                    final String clip) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        int aligned = 0;
        for (final htsjdk.samtools.CigarElement element : record.getCigar()) {
            if (element.getOperator() == htsjdk.samtools.CigarOperator.M) {
                aligned += element.getLength();
            }
        }
        final String bases = clip + referenceBases(start, aligned);
        record.setReadString(bases);
        final byte[] quality = new byte[bases.length()];
        Arrays.fill(quality, (byte) 40);
        record.setBaseQualities(quality);
        record.setMappingQuality(mappingQuality);
        record.setAttribute("RG", "rg1");
        record.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, new byte[bases.length()]);
        writer.addAlignment(record);
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ULTIMA");
        group.setFlowOrder(FLOW_ORDER);
        group.setSample("sample");
        header.addReadGroup(group);
        return header;
    }

    static Path writeReference(final Path dir, final String name, final String contig,
                               final String bases) throws Exception {
        final Path fasta = dir.resolve(name + ".fasta");
        final StringBuilder text = new StringBuilder(">" + contig + "\n");
        for (int i = 0; i < bases.length(); i += 60) {
            text.append(bases, i, Math.min(i + 60, bases.length())).append("\n");
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve(name + ".dict")});
        System.out.printf("fasta\t%s=%s%n", name, contig + ":" + bases.length());
        return fasta;
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static String asSam(final Path bam) throws Exception {
        final StringBuilder text = new StringBuilder();
        try (final htsjdk.samtools.SamReader reader =
                     htsjdk.samtools.SamReaderFactory.makeDefault().open(bam.toFile())) {
            for (final SAMRecord record : reader) {
                text.append(record.getSAMString());
            }
        }
        return text.toString();
    }

    static void run(final Path dir, final Path reference, final String label, final Path bam,
                    final Path maternal, final Path paternal, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".csv");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-R", reference.toString(),
                "--maternal-ref", maternal.toString(),
                "--paternal-ref", paternal.toString(),
                "--ancestral-translators-base-path", dir.toString() + "/",
                "--output-csv", out.toString(),
                // The tool refuses any other engine outright, so every run names this one.
                "--likelihood-calculation-engine", "FlowBased"));
        argv.addAll(extra);
        try {
            new GroundTruthReadsBuilder().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("none\t%s=no csv%n", label);
            return;
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
