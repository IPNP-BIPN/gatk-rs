/*
 * CRAMIssue8768Detector's report, taken from the reference.
 *
 * A CRAM walked container by container, looking for the containers that GATK issue 8768 would have
 * corrupted. What the tool decides is a state machine over container metadata, and what it prints
 * is a report whose numbers do not all mean what their labels say.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE SUSPECT IS THE CONTAINER AFTER THE ONE STARTING AT POSITION 1, not that container
 *     itself. The test reads the PREVIOUS container's alignment start, so a contig whose first
 *     container starts at 1 gets exactly one bad container, the second one, whatever that one's
 *     own start is;
 *   - THE FIRST CONTAINER OF A REFERENCE CONTEXT IS NEVER BAD, by two separate branches: the first
 *     container of the file, and the first after any change of reference context;
 *   - THE PRINTED "Mismatch Rate/Count" COUNT IS THE TOTAL BASE COUNT, not the mismatch count.
 *     `analyzeContainerBaseMismatches` returns `Tuple<>(totalBases, misMatches/(double) totalBases)`
 *     and the caller stores `.a` into a field named misMatchCount, so the rate is a mismatch rate
 *     and the count beside it is the denominator;
 *   - A BAD CONTIG IS KEYED BY REFERENCE ID, NOT BY NAME: the map key is
 *     `ReferenceContext.toString()`, which reads `SINGLE_REFERENCE: 0`. Only the TSV ever resolves
 *     a contig name;
 *   - ONLY FOUR GOOD CONTAINERS PER CONTIG ARE REPORTED unless --verbose, and the counter is set to
 *     1 by the new-context branch rather than 0, so the fifth container of a contig is the first
 *     one dropped;
 *   - MULTI-REF AND UNMAPPED CONTAINERS ARE COUNTED BUT NEVER RECORDED: the ordinal advances for
 *     them and they become the previous context, but recordContainerStats returns without adding a
 *     row, so ordinals in the report skip;
 *   - AN EMPTY CRAM DIVIDES BY ZERO: with no good containers the average is 0.0/0, and the report
 *     says `NaN`;
 *   - THE THRESHOLD LINE IS PRINTED ONLY WHEN EXCEEDED, and it prints the same number twice in two
 *     formats, `%f` for the measured rate and `%1.2f` for the threshold;
 *   - AND A FOREIGN CRAM STOPS THE ANALYSIS DEAD: isForeignCRAM returns true from inside the
 *     container loop and doAnalysis RETURNS, so no report body, no averages, and a return code of
 *     0 as though nothing were wrong.
 *
 * Output:
 *
 *     fixture\t<label>\tcontainers=<count>
 *     container\t<label>\t<ordinal in file>\tcontext=<toString>\tstart=<n>\tspan=<n>\
 *         \tslices=<n>\treference-required=<bool>\tembedded=<id>\tbases=<n>\tmismatches=<n>
 *     report\t<label>=<the whole text report, escaped>
 *     tsv\t<label>=<the whole tsv, escaped>
 *     stdout\t<label>=<what the tool wrote to System.out, escaped>
 *     code\t<label>\t<the value doWork returned>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CRAMIssue8768DetectorDump
 */

import htsjdk.samtools.CRAMContainerStreamWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CRAMReferenceRegion;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.cram.structure.Slice;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.CRAMIssue8768Detector;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CRAMIssue8768DetectorDump {

    /**
     * Ten bases matching the reference only where a read starts on a base congruent to 1 mod 4.
     *
     * The reference is ACGT repeated, so a read of ACGTACGTAC placed at 1 matches and the same
     * read placed at 11 is out of phase. That is deliberate: containers hold three reads and only
     * one of the three is in phase, which gives the good containers a rate that is neither zero
     * nor the same as the bad ones.
     */
    static final String MATCHING = "ACGTACGTAC";
    /** Ten bases matching only where the reference holds a T, which is a quarter of the time. */
    static final String MISMATCHING = "TTTTTTTTTT";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("cram-issue-8768-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CRAMIssue8768DetectorDump: the containers issue 8768 would have corrupted");

        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "cram8768",
                List.of("chr1", "chr2", "chr3"));
        final Path fasta = dir.resolve("cram8768.fasta");

        // The measured input. chr1 and chr3 both start at position 1, so both are bad contigs;
        // chr2 starts at 5 and never is, and carries five containers so the fourth-good cap bites.
        final Path cram = dir.resolve("reads.cram");
        writeCram(cram, fasta, dict);
        describe("input", cram, fasta);

        run(dir, "default", cram, fasta, List.of());
        run(dir, "tsv", cram, fasta, List.of("--output-tsv", dir.resolve("report.tsv").toString()));
        run(dir, "verbose", cram, fasta, List.of("--verbose", "true"));
        // A threshold under the measured bad rate, which prints the extra line, and one over it,
        // which does not.
        run(dir, "threshold-low", cram, fasta, List.of("--mismatch-rate-threshold", "0.01"));
        run(dir, "threshold-high", cram, fasta, List.of("--mismatch-rate-threshold", "0.9"));
        run(dir, "echo", cram, fasta, List.of("--echo-to-stdout", "true"));

        // A CRAM no container of which starts at position 1, which is the other branch of the
        // summary and a return code of 0.
        final Path clean = dir.resolve("clean.cram");
        writeCleanCram(clean, fasta, dict);
        describe("clean", clean, fasta);
        run(dir, "clean", clean, fasta, List.of());

        // Two slices to a container, which the tool reads as a file it did not write.
        final Path foreign = dir.resolve("foreign.cram");
        writeMultiSliceCram(foreign, fasta, dict);
        describe("foreign", foreign, fasta);
        run(dir, "foreign", foreign, fasta, List.of());

        // Nothing but an EOF container, so the average of no good containers is 0.0/0.
        final Path empty = dir.resolve("empty.cram");
        writeEmptyCram(empty, fasta, dict);
        describe("empty", empty, fasta);
        run(dir, "empty", empty, fasta, List.of());
    }

    /**
     * Every container's metadata, and the two numbers the tool computes from its records.
     *
     * This walks the file the way the analyzer does, so what is reported here is what the analyzer
     * had in hand: the port reproduces the report from these rows rather than from CRAM bytes.
     */
    static void describe(final String label, final Path cram, final Path fasta) throws Exception {
        final ReferenceSource referenceSource = new ReferenceSource(fasta);
        final CompressorCache cache = new CompressorCache();
        int ordinal = 0;
        try (final SeekablePathStream stream = new SeekablePathStream(cram)) {
            final CramHeader cramHeader = CramIO.readCramHeader(stream);
            final SAMFileHeader samHeader = Container.readSAMFileHeaderContainer(
                    cramHeader.getCRAMVersion(), stream, cram.toString());
            final SAMSequenceDictionary dictionary = samHeader.getSequenceDictionary();
            for (boolean isEOF = false; !isEOF;) {
                final Container container = new Container(
                        cramHeader.getCRAMVersion(), stream, stream.position());
                ordinal++;
                // The EOF container carries no compression header at all, so everything read off
                // it here is guarded rather than assumed.
                final List<Slice> slices = container.isEOF() ? List.of() : container.getSlices();
                long bases = -1;
                long mismatches = -1;
                if (container.getAlignmentContext().getReferenceContext().isMappedSingleRef()) {
                    final long[] counted = count(container, referenceSource, dictionary, samHeader,
                            cache);
                    bases = counted[0];
                    mismatches = counted[1];
                }
                System.out.printf(
                        "container\t%s\t%d\tcontext=%s\tstart=%d\tspan=%d\tslices=%d"
                                + "\treference-required=%s\tembedded=%d\tbases=%d\tmismatches=%d%n",
                        label,
                        ordinal,
                        container.getAlignmentContext().getReferenceContext().toString(),
                        container.getAlignmentContext().getAlignmentStart(),
                        container.getAlignmentContext().getAlignmentSpan(),
                        slices.size(),
                        container.getCompressionHeader() == null
                                ? "none"
                                : Boolean.toString(
                                        container.getCompressionHeader().isReferenceRequired()),
                        slices.isEmpty()
                                ? Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID
                                : slices.get(0).getEmbeddedReferenceContentID(),
                        bases,
                        mismatches);
                isEOF = container.isEOF();
            }
        }
        System.out.printf("fixture\t%s\tcontainers=%d%n", label, ordinal);
    }

    /** `analyzeContainerBaseMismatches`, as totals rather than as the tuple it returns. */
    static long[] count(final Container container, final ReferenceSource referenceSource,
                        final SAMSequenceDictionary dictionary, final SAMFileHeader samHeader,
                        final CompressorCache cache) {
        final List<SAMRecord> records = container.getSAMRecords(
                ValidationStringency.LENIENT,
                new CRAMReferenceRegion(referenceSource, dictionary),
                cache,
                samHeader);
        final CRAMReferenceRegion region = new CRAMReferenceRegion(referenceSource, dictionary);
        region.fetchReferenceBases(
                container.getAlignmentContext().getReferenceContext().getReferenceContextID());
        final byte[] referenceBases = region.getCurrentReferenceBases();
        long bases = 0;
        long mismatches = 0;
        for (final SAMRecord record : records) {
            bases += record.getReadLength();
            mismatches += SequenceUtil.countMismatches(record, referenceBases);
        }
        return new long[] {bases, mismatches};
    }

    static void run(final Path dir, final String label, final Path cram, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path report = dir.resolve("report-" + label + ".txt");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", cram.toString(),
                "-O", report.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);

        // The tool writes its summary to System.out whether or not --echo-to-stdout was given, so
        // that stream is captured and reported rather than left to land in this dump's own output.
        final PrintStream previous = System.out;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        Object code = null;
        String failure = null;
        try {
            System.setOut(new PrintStream(captured, true, StandardCharsets.UTF_8));
            try {
                code = new CRAMIssue8768Detector().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception | AssertionError e) {
                failure = e.getClass().getName() + ":" + masked(String.valueOf(e.getMessage()), dir);
            }
        } finally {
            System.setOut(previous);
        }

        System.out.printf("stdout\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(captured.toString(StandardCharsets.UTF_8), dir)));
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(failure));
            return;
        }
        System.out.printf("code\t%s\t%s%n", label, code);
        if (Files.exists(report)) {
            System.out.printf("report\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(report), dir)));
        }
        final Path tsv = dir.resolve("report.tsv");
        if (extra.contains("--output-tsv") && Files.exists(tsv)) {
            System.out.printf("tsv\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(tsv), dir)));
        }
    }

    /**
     * Three reads to a container, and three contigs that differ in what the tool makes of them.
     *
     * chr1 and chr3 both open at position 1, so on both the SECOND container is the suspect. chr2
     * opens at 5 and has five containers, so it is the contig where the four-good-containers cap
     * decides what the report shows.
     */
    static void writeCram(final Path cram, final Path fasta, final Path dict) throws Exception {
        final SAMFileHeader header = readHeader(dict);
        try (final java.io.OutputStream out = Files.newOutputStream(cram)) {
            final CRAMContainerStreamWriter writer = writer(cram, fasta, header, out, 3, 1);
            writer.writeHeader();
            // chr1: three containers, the first one opening at position 1. The second holds reads
            // that do not match the reference, which is what a corrupt container looks like.
            write(writer, header, "chr1", 1, 3, MATCHING);
            write(writer, header, "chr1", 31, 3, MISMATCHING);
            write(writer, header, "chr1", 61, 3, MATCHING);
            // chr2: five containers, none opening at position 1, so none is ever suspected.
            write(writer, header, "chr2", 5, 15, MATCHING);
            // chr3: two containers, the first opening at position 1.
            write(writer, header, "chr3", 1, 3, MATCHING);
            write(writer, header, "chr3", 41, 3, MISMATCHING);
            // Unmapped reads, whose container advances the ordinal and records nothing.
            for (int index = 0; index < 3; index++) {
                writer.writeAlignment(read(header, "unmapped-" + index, null, 0, MATCHING));
            }
            writer.finish(true);
        }
    }

    /** The same shape, with nothing on the first base of any contig. */
    static void writeCleanCram(final Path cram, final Path fasta, final Path dict) throws Exception {
        final SAMFileHeader header = readHeader(dict);
        try (final java.io.OutputStream out = Files.newOutputStream(cram)) {
            final CRAMContainerStreamWriter writer = writer(cram, fasta, header, out, 3, 1);
            writer.writeHeader();
            write(writer, header, "chr1", 5, 6, MATCHING);
            write(writer, header, "chr2", 9, 6, MISMATCHING);
            writer.finish(true);
        }
    }

    /** Two slices to a container, which is the shape the tool calls foreign. */
    static void writeMultiSliceCram(final Path cram, final Path fasta, final Path dict)
            throws Exception {
        final SAMFileHeader header = readHeader(dict);
        try (final java.io.OutputStream out = Files.newOutputStream(cram)) {
            final CRAMContainerStreamWriter writer = writer(cram, fasta, header, out, 2, 2);
            writer.writeHeader();
            write(writer, header, "chr1", 1, 8, MATCHING);
            writer.finish(true);
        }
    }

    static void writeEmptyCram(final Path cram, final Path fasta, final Path dict) throws Exception {
        final SAMFileHeader header = readHeader(dict);
        try (final java.io.OutputStream out = Files.newOutputStream(cram)) {
            final CRAMContainerStreamWriter writer = new CRAMContainerStreamWriter(
                    new CRAMEncodingStrategy(), new ReferenceSource(fasta), header, out, null,
                    cram.getFileName().toString());
            writer.writeHeader();
            writer.finish(true);
        }
    }

    /**
     * A writer whose containers are small enough to be a fixture.
     *
     * setMinimumSingleReferenceSliceSize comes first because setReadsPerSlice validates against it,
     * and the default ten thousand reads to a slice would give one container for the whole file.
     */
    static CRAMContainerStreamWriter writer(final Path cram, final Path fasta,
                                            final SAMFileHeader header,
                                            final java.io.OutputStream out,
                                            final int readsPerSlice,
                                            final int slicesPerContainer) {
        final CRAMEncodingStrategy strategy = new CRAMEncodingStrategy()
                .setMinimumSingleReferenceSliceSize(1)
                .setReadsPerSlice(readsPerSlice)
                .setSlicesPerContainer(slicesPerContainer);
        return new CRAMContainerStreamWriter(strategy, new ReferenceSource(fasta), header, out,
                null, cram.getFileName().toString());
    }

    static void write(final CRAMContainerStreamWriter writer, final SAMFileHeader header,
                      final String contig, final int start, final int count, final String bases) {
        for (int index = 0; index < count; index++) {
            writer.writeAlignment(read(header, contig + "-" + (start + index * 10), contig,
                    start + index * 10, bases));
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        if (contig == null) {
            record.setReadUnmappedFlag(true);
        } else {
            record.setReferenceName(contig);
            record.setAlignmentStart(start);
            record.setCigarString("10M");
            record.setMappingQuality(60);
        }
        return record;
    }

    static SAMFileHeader readHeader(final Path dict) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(dict.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            return header;
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
