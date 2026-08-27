/*
 * CalibrateDragstrModel's parameter table, taken from the reference.
 *
 * The DRAGstr model is a table of three parameters per period and repeat length, estimated from
 * the reads that pile up over a reference's short tandem repeats. The maximum-likelihood search
 * itself is not what this measures: what it measures is the table around it, which sites the tool
 * is willing to use, and what a period with no data gets instead.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE OUTPUT IS THREE BLOCKS, GOP, GCP and API, one row per period and one column per repeat
 *     length, under a header that repeats the whole command line;
 *   - THE TABLE'S SHAPE IS THE HYPER-PARAMETERS' AND NOT THE DATA'S: --max-period and
 *     --max-repeats decide how many rows and columns there are whatever the reference holds;
 *   - A PERIOD WITH NO DATA KEEPS THE DEFAULTS rather than being left out, so the table always
 *     has every cell filled and a reader cannot tell an estimate from a default by its shape;
 *   - GCP IS NEVER ESTIMATED: every row of it is ten over the period, repeated across the row;
 *   - GOP AND API ARE CONSTANT ACROSS A ROW too, the estimation being per period rather than per
 *     repeat length once the repeat lengths are grouped;
 *   - A SITE IS ONLY USED IF ITS DEPTH REACHES --minimum-depth, which is ten: a site with fewer
 *     reads is `skipped` rather than `used`, and the sites table is the only place that shows;
 *   - AND THE REST ARE `downsampled-out`, --down-sample-size being a cap per period and repeat
 *     length rather than a total;
 *   - AN HOMOPOLYMER WRITTEN AS A LONGER UNIT IS STILL PERIOD ONE, so a fixture that wants
 *     period two has to spell it with two different bases;
 *   - THE STR TABLE MUST HAVE BEEN COMPOSED FOR THE SAME REFERENCE, one composed for another
 *     being refused by its dictionary;
 *   - A --max-period BELOW THE TABLE'S OWN IS AN ArrayIndexOutOfBoundsException rather than a
 *     refusal, the estimator allocating from the argument and indexing from the file;
 *   - AND --down-sample-size AND --shard-size EACH HAVE A MINIMUM, refused by the parser.
 *
 * Output:
 *
 *     fixture\t<name>=<that number or sequence>
 *     out\t<label>=<the whole parameter table, escaped>
 *     rows\t<label>=<how many sites the debug table holds>
 *     status\t<label>\t<outcome>=<how many sites took it>
 *     shape\t<label>\t<period>,<repeats>,<outcome>=<how many>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CalibrateDragstrModelDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.dragstr.CalibrateDragstrModel;
import org.broadinstitute.hellbender.tools.dragstr.ComposeSTRTableFile;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CalibrateDragstrModelDump {

    /** How many repeat blocks the reference carries, and how far apart they sit. */
    static final int BLOCKS = 60;
    static final int SPACING = 300;
    static final int READ_LENGTH = 100;

    /**
     * The reference: a filler of `ACGT` with a repeat block every {@link #SPACING} bases.
     *
     * The blocks cycle through periods one to three and through repeat lengths, so the table has
     * data in more than one cell.
     */
    static String reference() {
        final StringBuilder bases = new StringBuilder();
        for (int block = 0; block < BLOCKS; block++) {
            while (bases.length() < block * SPACING) {
                bases.append("ACGT".charAt(bases.length() % 4));
            }
            final int period = (block % 3) + 1;
            final int repeats = 4 + (block % 7);
            // A period of two written `AA` is an homopolymer, which the scanner calls period one,
            // so each period gets a unit whose own bases differ.
            final String unit = period == 1 ? "A" : period == 2 ? "AC" : "AGT";
            for (int i = 0; i < repeats; i++) {
                bases.append(unit);
            }
        }
        while (bases.length() < BLOCKS * SPACING) {
            bases.append("ACGT".charAt(bases.length() % 4));
        }
        return bases.toString();
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("calibrate-dragstr-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CalibrateDragstrModelDump: the DRAGstr parameter table estimated "
                + "from the reads over a reference's repeats");

        final String bases = reference();
        System.out.printf("fixture\tlength=%d%n", bases.length());
        System.out.printf("fixture\tfirst-200=%s%n", bases.substring(0, 200));

        final Path fasta = writeReference(dir, "reference", bases);
        final Path table = composeTable(dir, "table", fasta, 3, 20);

        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam, bases);

        run(dir, "forced", table, bam, fasta, List.of("--force-estimation", "true"));
        // Without the force, the same data is too little.
        run(dir, "unforced", table, bam, fasta, List.of());
        // The table's shape is the hyper-parameters' and not the data's.
        run(dir, "period-two", table, bam, fasta,
                List.of("--force-estimation", "true", "--max-period", "2"));
        run(dir, "repeat-eight", table, bam, fasta,
                List.of("--force-estimation", "true", "--max-repeats", "8"));
        // The two hyper-parameters that decide how much data a repeat length needs.
        run(dir, "min-loci-one", table, bam, fasta,
                List.of("--force-estimation", "true", "--min-loci-count", "1"));
        run(dir, "min-depth-ten", table, bam, fasta,
                List.of("--force-estimation", "true", "--minimum-depth", "10"));
        // In parallel, which should change nothing about the table.
        run(dir, "parallel", table, bam, fasta,
                List.of("--force-estimation", "true", "--parallel", "true", "--threads", "2"));
        // The two arguments with a minimum.
        run(dir, "down-sample-too-small", table, bam, fasta,
                List.of("--force-estimation", "true", "--down-sample-size", "10"));
        run(dir, "shard-too-small", table, bam, fasta,
                List.of("--force-estimation", "true", "--shard-size", "10"));
        // A table composed for another reference entirely.
        final Path other = writeReference(dir, "other", bases.substring(0, 5000));
        final Path otherTable = composeTable(dir, "other-table", other, 3, 20);
        run(dir, "wrong-reference", otherTable, bam, fasta,
                List.of("--force-estimation", "true"));
    }

    /** The reads: one per repeat block, plus a second carrying an error in the repeat. */
    static void writeReads(final Path bam, final String bases) {
        final SAMFileHeader header = readHeader(bases.length());
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (int block = 0; block < BLOCKS; block++) {
                final int start = Math.max(1, block * SPACING - 20);
                if (start + READ_LENGTH > bases.length()) {
                    break;
                }
                // Twelve reads over each block: the estimator wants --minimum-depth of them,
                // which is ten by default, and a site with fewer is skipped rather than used.
                for (int copy = 0; copy < 12; copy++) {
                    final String read = bases.substring(start - 1, start - 1 + READ_LENGTH);
                    writer.addAlignment(record(header, "r" + block + "-" + copy, start, read,
                            copy % 4 == 3));
                }
            }
        }
    }

    static SAMRecord record(final SAMFileHeader header, final String name, final int start,
                            final String bases, final boolean reverse) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(bases.length() + "M");
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quality = new byte[bases.length()];
        Arrays.fill(quality, (byte) 35);
        record.setBaseQualities(quality);
        record.setMappingQuality(60);
        record.setReadNegativeStrandFlag(reverse);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static SAMFileHeader readHeader(final int length) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", length))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ILLUMINA");
        group.setSample("sample");
        header.addReadGroup(group);
        return header;
    }

    static Path writeReference(final Path dir, final String name, final String bases)
            throws Exception {
        final Path fasta = dir.resolve(name + ".fasta");
        final StringBuilder text = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            text.append(bases, i, Math.min(i + 60, bases.length())).append("\n");
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve(name + ".dict")});
        return fasta;
    }

    static Path composeTable(final Path dir, final String name, final Path fasta,
                             final int maxPeriod, final int maxRepeat) throws Exception {
        final Path out = dir.resolve(name + ".zip");
        new ComposeSTRTableFile().instanceMain(new String[] {
                "-R", fasta.toString(),
                "-O", out.toString(),
                "--decimation", "NONE",
                "--max-period", Integer.toString(maxPeriod),
                "--max-repeats", Integer.toString(maxRepeat)});
        return out;
    }

    static void run(final Path dir, final String label, final Path table, final Path bam,
                    final Path fasta, final List<String> extra) throws Exception {
        final Path out = dir.resolve("out-" + label + ".txt");
        final Path sites = dir.resolve("sites-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(List.of(
                "-R", fasta.toString(),
                "-I", bam.toString(),
                "--str-table-path", table.toString(),
                "-O", out.toString(),
                "--debug-sites-output", sites.toString()));
        argv.addAll(extra);
        try {
            new CalibrateDragstrModel().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("out\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        } else {
            System.out.printf("none\t%s=no parameter table%n", label);
        }
        if (Files.exists(sites)) {
            // The sites table has one row per site, tens of thousands of them, so what is
            // reported is a census: how many sites each outcome took, and how many of each
            // period and repeat length reached the estimator at all.
            final java.util.TreeMap<String, Integer> byStatus = new java.util.TreeMap<>();
            final java.util.TreeMap<String, Integer> byShape = new java.util.TreeMap<>();
            int rows = 0;
            for (final String line : Files.readString(sites).split("\n", -1)) {
                if (line.isEmpty()) {
                    continue;
                }
                rows++;
                final String[] columns = line.split("\t");
                final String status = columns[columns.length - 1];
                byStatus.merge(status, 1, Integer::sum);
                if (!status.equals("downsampled-out")) {
                    byShape.merge(columns[1] + "," + columns[2] + "," + status, 1, Integer::sum);
                }
            }
            System.out.printf("rows\t%s=%d%n", label, rows);
            for (final java.util.Map.Entry<String, Integer> entry : byStatus.entrySet()) {
                System.out.printf("status\t%s\t%s=%d%n", label, entry.getKey(),
                        entry.getValue());
            }
            for (final java.util.Map.Entry<String, Integer> entry : byShape.entrySet()) {
                System.out.printf("shape\t%s\t%s=%d%n", label, entry.getKey(), entry.getValue());
            }
        } else {
            System.out.printf("none\t%s=no sites table%n", label);
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
