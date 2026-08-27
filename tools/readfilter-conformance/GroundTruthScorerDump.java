/*
 * GroundTruthScorer's CSV and its report, taken from the reference.
 *
 * Every read scored against the reference it aligns to, written as one CSV row and summarised in
 * a GATK report of four tables. The flow-based scoring is not what this measures: what it
 * measures is which reads are scored at all, what the row carries, and how the report is cut.
 *
 * Ten behaviours this is built to catch.
 *
 *   - --omit-zeros-from-report IS NOT OPTIONAL IN PRACTICE: the four-level table is allocated as
 *     61 x 101 x 202 x 4 accumulators and writes every one of them out, so a report that keeps
 *     its zeros runs to five million rows. Every run below passes the option, and none of them
 *     measures what leaving it off would produce;
 *   - THE REPORT IS FOUR TABLES, named PhredBinAccumulator, qualReport, qual_hmerReport and
 *     qual_hmer_deviation_base_Report, in that order;
 *   - --quality-percentiles NAMES COLUMNS OF THE FIRST TABLE, whose width is seven plus however
 *     many percentiles were asked for: the default five make twelve columns and three make ten;
 *   - --exclude-zero-flows EMPTIES THE FLOWS WHOSE CALL IS ZERO rather than skipping their rows,
 *     so the row stays with a count of zero and every one of its statistics comes out `NaN`;
 *   - --add-mean-call ADDS TWO COLUMNS to the CSV, taking it from fifteen to seventeen, rather
 *     than changing the ones already there;
 *   - --normalized-score-threshold IS COMPARED AGAINST A NEGATIVE DEFAULT of -0.1, so a
 *     threshold of 0.1 is already above every score and empties the CSV;
 *   - AN EMPTY CSV IS STILL WRITTEN WITH ITS HEADER, and so is the report beside it;
 *   - --gt-no-output EMPTIES THE CSV THE SAME WAY while leaving the report untouched, so the two
 *     arguments are told apart by the report and not by the CSV;
 *   - --use-softclipped-bases CHANGES BOTH FILES, the clipped read's row and the histograms it
 *     feeds;
 *   - AND A READ THAT IS NOT FLOW-BASED IS REFUSED BY NAME AND POSITION, the tool having nothing
 *     to score it with.
 *
 * Output:
 *
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     out\t<label>=<the whole output csv, escaped>
 *     report\t<label>=<the whole report file, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GroundTruthScorerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthScorer;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GroundTruthScorerDump {

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

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("ground-truth-scorer-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GroundTruthScorerDump: every read scored against the reference it "
                + "aligns to, and the report that summarises them");

        final Path fasta = writeReference(dir, referenceBases(1, LENGTH));
        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam, false);
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(bam)));

        run(dir, "default", bam, fasta, List.of());
        // The percentile columns, named and defaulted.
        run(dir, "percentiles", bam, fasta,
                List.of("--quality-percentiles", "5,50,95"));
        run(dir, "exclude-zero-flows", bam, fasta,
                List.of("--exclude-zero-flows", "true"));
        // The two extra columns.
        run(dir, "mean-call", bam, fasta, List.of("--add-mean-call", "true"));
        // The score threshold, once mild and once above every score.
        run(dir, "threshold-mild", bam, fasta,
                List.of("--normalized-score-threshold", "0.1"));
        run(dir, "threshold-above-everything", bam, fasta,
                List.of("--normalized-score-threshold", "2.0"));
        // The soft-clipped bases scored rather than skipped.
        run(dir, "use-softclipped", bam, fasta, List.of("--use-softclipped-bases", "true"));
        // No CSV output at all, which still writes the report.
        run(dir, "no-output", bam, fasta, List.of("--gt-no-output", "true"));

        // A BAM whose reads carry no flow matrix at all.
        final Path plain = dir.resolve("plain.bam").toAbsolutePath();
        writeReads(plain, true);
        System.out.printf("sam\tplain=%s%n", ReferenceQueryDump.escape(asSam(plain)));
        run(dir, "not-flow-based", plain, fasta, List.of());
    }

    /** The reads: five over the reference, one of them soft-clipped. */
    static void writeReads(final Path bam, final boolean withoutFlowMatrix) {
        final SAMFileHeader header = readHeader(withoutFlowMatrix);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            add(writer, header, "r-plain", 1000, "100M", "", withoutFlowMatrix);
            add(writer, header, "r-second", 1100, "100M", "", withoutFlowMatrix);
            add(writer, header, "r-third", 1200, "100M", "", withoutFlowMatrix);
            add(writer, header, "r-clipped", 1300, "8S92M", "TTTTGGCC", withoutFlowMatrix);
            add(writer, header, "r-fifth", 1400, "100M", "", withoutFlowMatrix);
        }
    }

    static void add(final SAMFileWriter writer, final SAMFileHeader header, final String name,
                    final int start, final String cigar, final String clip,
                    final boolean withoutFlowMatrix) {
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
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        if (!withoutFlowMatrix) {
            record.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, new byte[bases.length()]);
        }
        writer.addAlignment(record);
    }

    static SAMFileHeader readHeader(final boolean withoutFlowOrder) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample");
        if (withoutFlowOrder) {
            // A read group with no flow order and no flow platform, which is what makes a read
            // not flow-based.
            group.setPlatform("ILLUMINA");
        } else {
            group.setPlatform("ULTIMA");
            group.setFlowOrder(FLOW_ORDER);
        }
        header.addReadGroup(group);
        return header;
    }

    static Path writeReference(final Path dir, final String bases) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder text = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            text.append(bases, i, Math.min(i + 60, bases.length())).append("\n");
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("reference.dict")});
        return fasta;
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

    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".csv");
        final Path report = dir.resolve("report-" + label + ".txt");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-R", fasta.toString(),
                "--output-csv", out.toString(),
                "--report-file", report.toString(),
                "--likelihood-calculation-engine", "FlowBased"));
        // EVERY run omits the zero rows. The four-level table is allocated as
        // 61 x 101 x 202 x 4 accumulators, so writing its zeros out means five million rows and
        // a report no golden can hold: the option is not a variation here, it is what makes the
        // tool usable at all.
        if (!extra.contains("--omit-zeros-from-report")) {
            argv.addAll(List.of("--omit-zeros-from-report", "true"));
        }
        argv.addAll(extra);
        try {
            new GroundTruthScorer().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        for (final String[] pair : new String[][] {{"out", out.toString()},
                {"report", report.toString()}}) {
            final Path file = Path.of(pair[1]);
            if (Files.exists(file)) {
                System.out.printf("%s\t%s=%s%n", pair[0], label,
                        ReferenceQueryDump.escape(masked(Files.readString(file), dir)));
            } else {
                System.out.printf("none\t%s=no %s%n", label, pair[0]);
            }
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
