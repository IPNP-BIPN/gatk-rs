/*
 * ApplyBQSR's output, taken from the reference.
 *
 * The eighth whole tool of the record-transform archetype, and the one every BQSR brick was built
 * for. Its body is four lines: make a BQSRReadTransformer from the recalibration file, write every
 * read the traversal hands it, close the writer. Everything interesting is in where that transformer
 * is hooked and in what the traversal hands over.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE TRANSFORMER IS A POST-READ-FILTER ONE. `makePostReadFilterTransformer` runs AFTER the
 *     read filters, so a read the filters drop is never recalibrated AND never written. The default
 *     filter is WellformedReadFilter, so a read with no read group or with mismatched qualities
 *     disappears from the output rather than passing through unrecalibrated;
 *   - THE RECALIBRATION FILE DECIDES THE READ GROUP KEYS, NOT THE BAM. The covariates come from the
 *     report's own RecalTable0, so a BAM whose read groups the report does not name is refused
 *     unless --allow-missing-read-group is set, and then those reads are quantized but not
 *     recalibrated;
 *   - THE WRITER IS CREATED PRESORTED. `createSAMWriter(output, true)` promises the reads are
 *     already in order, so the output keeps the traversal's order and an index is written beside it;
 *   - THE OUTPUT HEADER GAINS A @PG RECORD for this tool, and its command line is the one the
 *     argument parser reconstructed rather than the one given;
 *   - --emit-original-quals WRITES THE OQ TAG ONLY WHERE THERE IS NONE, so a read that already
 *     carries one keeps its own;
 *   - AND --quantize-quals AND --static-quantized-quals ARE SEPARATE MECHANISMS that both reach the
 *     output, the second applying on top of the first.
 *
 * Output, one row per (label, kind):
 *
 *     recal\t<label>\t<the recalibration table as text, escaped>
 *     header\t<label>\t<the output header, escaped>
 *     commandline\t<label>\t<the @PG command line>
 *     output\t<label>\t<the output BAM, base64>
 *     index\t<label>\t<the index, base64, or absent>
 *     error\t<label>\t<exception>\t<message>
 *
 * Usage: ApplyBqsrDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.tools.walkers.bqsr.ApplyBQSR;
import org.broadinstitute.hellbender.utils.recalibration.EventType;
import org.broadinstitute.hellbender.utils.recalibration.QuantizationInfo;
import org.broadinstitute.hellbender.utils.recalibration.RecalDatum;
import org.broadinstitute.hellbender.utils.recalibration.RecalUtils;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.utils.recalibration.covariates.ContextCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.CycleCovariate;
import org.broadinstitute.hellbender.utils.recalibration.covariates.StandardCovariateList;
import org.broadinstitute.hellbender.utils.report.GATKReport;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ApplyBqsrDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# ApplyBqsrDump: ApplyBQSR's output, from the reference");
        System.out.printf("deflater\t%s%n",
                htsjdk.samtools.util.BlockCompressedOutputStream.getDefaultDeflaterFactory()
                        .getClass().getName());

        // Before the fixture is written: the factory is static and first writer wins.
        htsjdk.samtools.util.BlockCompressedOutputStream.setDefaultDeflaterFactory(
                new htsjdk.samtools.util.zip.DeflaterFactory());

        // Relative on purpose, for the reason RecordTransformDump gives: the string handed to -I and
        // -O is the string recorded inside the output BAM's own @PG, so an absolute temporary path
        // would make every output byte unstable and canonicalisation cannot reach inside base64.
        final Path dir = Path.of("applybqsr-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);
        final Path bam = dir.resolve("input.bam");
        buildFixture(bam.toFile());

        // The recalibration table this tool is given, written by the reference from tables built
        // here. Its read groups are what the covariates are numbered from.
        final Path recal = dir.resolve("recal.table");
        writeRecalTable(recal, List.of("rg1"));
        System.out.printf("recal\tone-group\t%s%n",
                ReferenceQueryDump.escape(Files.readString(recal)));

        // A second table naming a read group the BAM does not use, so every read of the BAM is
        // outside it.
        final Path other = dir.resolve("other.table");
        writeRecalTable(other, List.of("elsewhere"));
        System.out.printf("recal\telsewhere\t%s%n",
                ReferenceQueryDump.escape(Files.readString(other)));

        run(dir, bam, recal, "plain", new String[] {});
        run(dir, bam, recal, "emit-original-quals",
                new String[] {"--emit-original-quals", "true"});
        run(dir, bam, recal, "quantize-4", new String[] {"--quantize-quals", "4"});
        run(dir, bam, recal, "no-quantization", new String[] {"--quantize-quals", "0"});
        run(dir, bam, recal, "static-quals",
                new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "30"});
        run(dir, bam, recal, "static-quals-round-down",
                new String[] {"--static-quantized-quals", "10", "--static-quantized-quals", "30",
                        "--round-down-quantized", "true"});
        run(dir, bam, recal, "preserve-25",
                new String[] {"--preserve-qscores-less-than", "25"});
        run(dir, bam, recal, "global-prior",
                new String[] {"--global-qscore-prior", "20"});

        // The read group the table does not name, refused and then allowed.
        run(dir, bam, other, "missing-read-group", new String[] {});
        run(dir, bam, other, "missing-read-group-allowed",
                new String[] {"--allow-missing-read-group", "true"});

        // A BAM with a read the default filters drop, to show the transformer never sees it.
        final Path filtered = dir.resolve("filtered.bam");
        buildFilteredFixture(filtered.toFile());
        // The inputs travel too, so the port opens the same bytes rather than a rebuilt lookalike.
        for (final String name : new String[] {"input", "filtered"}) {
            System.out.printf("fixture\t%s\t%s%n", name,
                    RecordTransformDump.base64(dir.resolve(name + ".bam")));
            System.out.printf("fixtureindex\t%s\t%s%n", name,
                    RecordTransformDump.base64(dir.resolve(name + ".bai")));
        }
        run(dir, filtered, recal, "filtered-read", new String[] {});
        run(dir, filtered, recal, "filters-disabled",
                new String[] {"--disable-read-filter", "WellformedReadFilter"});
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.delete(entry);
            }
        }
    }

    /** Six reads in one read group, at qualities the recalibration table has datums for. */
    static void buildFixture(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            for (int i = 0; i < 6; i++) {
                final SAMRecord record = read(header, "r" + i, 100 + i * 20, "rg1");
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) (i % 2 == 0 ? 20 : 30));
                record.setBaseQualities(quals);
                // One read already carries an OQ, so --emit-original-quals has one to leave alone.
                if (i == 2) {
                    record.setAttribute("OQ", "!!!!!!!!!!");
                }
                writer.addAlignment(record);
            }
        }
    }

    /** The same, plus one read WellformedReadFilter drops: its qualities are the wrong length. */
    static void buildFilteredFixture(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            final SAMRecord good = read(header, "good", 100, "rg1");
            final byte[] quals = new byte[10];
            Arrays.fill(quals, (byte) 30);
            good.setBaseQualities(quals);
            writer.addAlignment(good);

            // A cigar that does not match the read's length, which WellformedReadFilter's
            // READ_LENGTH_EQUALS_CIGAR_LENGTH drops. Qualities of the wrong length would do it too
            // and cannot be written: the BAM index builder refuses them before any filter runs.
            final SAMRecord malformed = read(header, "malformed", 200, "rg1");
            malformed.setCigarString("5M");
            final byte[] malformedQuals = new byte[10];
            Arrays.fill(malformedQuals, (byte) 30);
            malformed.setBaseQualities(malformedQuals);
            writer.addAlignment(malformed);
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String group) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes(StandardCharsets.UTF_8));
        record.setMappingQuality(60);
        record.setAttribute("RG", group);
        return record;
    }

    /** A recalibration table with a datum in every place the transformer looks. */
    static void writeRecalTable(final Path path, final List<String> readGroups) throws Exception {
        final RecalibrationArgumentCollection rac = new RecalibrationArgumentCollection();
        final StandardCovariateList covariates = new StandardCovariateList(rac, readGroups);
        final RecalibrationTables tables = new RecalibrationTables(covariates);

        for (int group = 0; group < readGroups.size(); group++) {
            for (final EventType event : EventType.values()) {
                final RecalDatum datum = new RecalDatum(100000L, 1000.0, (byte) 1);
                datum.setReportedQuality(30.0);
                tables.getReadGroupTable().put(datum, group, event.ordinal());
            }
            for (final int quality : new int[] {20, 30}) {
                for (final EventType event : EventType.values()) {
                    final RecalDatum datum = new RecalDatum(10000L, 50.0, (byte) 1);
                    datum.setReportedQuality(quality);
                    tables.getQualityScoreTable().put(datum, group, quality, event.ordinal());
                }
                // Every context this fixture's reads produce, and every cycle key of a ten-base
                // forward read.
                for (final String context : new String[] {"AC", "CG", "GT", "TA"}) {
                    final RecalDatum datum = new RecalDatum(1000L, 5.0, (byte) 1);
                    datum.setReportedQuality(quality);
                    tables.getTable(2).put(datum, group, quality,
                            ContextCovariate.keyFromContext(context), 0);
                }
                for (int cycle = 1; cycle <= 10; cycle++) {
                    final RecalDatum datum = new RecalDatum(1000L, 3.0, (byte) 1);
                    datum.setReportedQuality(quality);
                    tables.getTable(3).put(datum, group, quality,
                            CycleCovariate.keyFromCycle(cycle, 500), 0);
                }
            }
        }

        final QuantizationInfo quantization = new QuantizationInfo(tables, rac.QUANTIZING_LEVELS);
        final GATKReport report = RecalUtils.createRecalibrationGATKReport(
                rac.generateReportTable(covariates.covariateNames()), quantization, tables,
                covariates);
        try (final PrintStream out = new PrintStream(path.toFile(), StandardCharsets.UTF_8)) {
            report.print(out);
        }
    }

    /** One run of the tool, with its output BAM, its header and its index. */
    static void run(final Path dir, final Path input, final Path recal, final String label,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve("ApplyBQSR." + label + ".bam");
        // --use-jdk-deflater for the same reason every other record-transform dump names it: the
        // GKL deflater's bytes are not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "--bqsr-recal-file", recal.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new ApplyBQSR().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getSimpleName(),
                    e.getMessage());
            return;
        }

        String commandLine = "";
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(output.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMProgramRecord record : header.getProgramRecords()) {
                if (record.getCommandLine() != null) {
                    commandLine = record.getCommandLine();
                }
            }
            System.out.printf("header\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(header.getSAMString()));
        }
        System.out.printf("commandline\t%s\t%s%n", label, commandLine);
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }
}
