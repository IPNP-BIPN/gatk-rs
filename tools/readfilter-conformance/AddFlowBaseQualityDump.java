/*
 * AddFlowBaseQuality's output, taken from the reference.
 *
 * A flow-based read's indel error probabilities turned into a per-base mismatch quality. What the
 * tool computes is an enumeration over the ways the flow key could have been misread, and what it
 * writes back is a string whose middle characters are not computed at all.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE MIDDLE BASES OF AN HMER ARE NEVER COMPUTED: generateBaseErrorProbability writes the
 *     first base and the last base of each hmer and advances `base` past the rest, so a run of
 *     five bases keeps three error probabilities of exactly ZERO, which convertErrorProbToPhred
 *     turns into the MAXIMAL quality rather than the minimal one;
 *   - THE FIRST BASE OF THE READ IS OVERRIDDEN with the hmer's own key probability rather than
 *     with the computed one, and so is the LAST BASE OF THE READ, by a second override that fires
 *     whenever the write cursor has reached the end of the array;
 *   - AN HMER OF LENGTH ONE GETS ITS SECOND PROBABILITY LEFT AT ZERO and walks BOTH sides:
 *     generateHmerBaseErrorProbabilities skips the right-hand call, and
 *     generateSidedHmerBaseErrorProbability walks `sideIncr` and `-sideIncr` for that case alone;
 *   - THE SLICE WINDOW IS THE FLOW ORDER'S CYCLE LENGTH, computed by calcFlowOrderLength as the
 *     distance to the SECOND occurrence of the order's first base, so a flow order of TTGCA gives
 *     a cycle of ONE and a window of a single flow;
 *   - A SLICE IS INVALID WHEN IT HOLDS flowOrderLength - 1 CONSECUTIVE ZEROS, so with a cycle of
 *     one every zero invalidates, and the enumeration collapses to the key slice alone;
 *   - THE ERROR RATE FLOOR IS APPLIED TO EVERY BAND INCLUDING THE IMPOSSIBLE ONES: a key of 0 has
 *     no shorter neighbour and a key at maxHmer has no longer one, and both get the floor rather
 *     than zero;
 *   - THE PHRED CONVERSION TRUNCATES rather than rounds, `(int)(-10 * log10(p))`, and clamps to
 *     --maximal-quality-score, with a probability of exactly zero taking the clamp directly;
 *   - AND --replace-quality-mode MOVES THE OLD QUALITIES TO OQ and overwrites QUAL, where the
 *     default writes the new string to XQ and leaves QUAL alone.
 *
 * Output:
 *
 *     read\t<name>\tgroup=<rg>\tbases=<bases>\tqual=<phred string>\ttp=<comma separated>
 *     group\t<rg>\torder=<flow order>\tmaxClass=<n>
 *     flow\t<name>\tkey=<comma separated>\tmaxHmer=<n>
 *     prob\t<name>\t<flow>\tminus=<double>\tkey=<double>\tplus=<double>
 *     xq\t<label>\t<name>=<the XQ attribute>
 *     qual\t<label>\t<name>=<the quality string after the run>
 *     oq\t<label>\t<name>=<the OQ attribute>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AddFlowBaseQualityDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.FlowBasedArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.groundtruth.AddFlowBaseQuality;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AddFlowBaseQualityDump {

    /** The read group every read but one belongs to, whose flow order cycles every four flows. */
    static final String FOUR = "TGCA";
    /**
     * A flow order whose first base repeats immediately, so calcFlowOrderLength answers ONE and
     * both the slice window and the validity rule collapse.
     */
    static final String ONE = "TTGCA";

    record Read(String name, String group, String bases, String quals, byte[] tp) { }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("add-flow-base-quality-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AddFlowBaseQualityDump: flow indel probabilities as a base quality");

        final List<Read> reads = List.of(
                // Every hmer of length one, so every hmer walks both sides.
                new Read("singles", "rg4", "TGCATGCA", "IIIIIIII", tp(8)),
                // Hmers of length two, which is the shortest hmer with a computed last base and no
                // middle bases at all.
                new Read("doubles", "rg4", "TTGGCCAA", "IIIIIIII", tp(8)),
                // A five-base run, whose three middle bases are never computed and come out at the
                // maximal quality.
                new Read("long-hmer", "rg4", "TTTTTGCA", "IIIIIIII", tp(8)),
                // A read that opens on a zero flow, so the first key entry is 0 and its shorter
                // band is the floor rather than a probability.
                new Read("leading-zero", "rg4", "GCATGCA", "IIIIIII", tp(7)),
                // A read ending inside a long hmer, which is where the last-base override fires on
                // a base the hmer loop had already written.
                new Read("trailing-hmer", "rg4", "TGCAAAA", "IIIIIII", tp(7)),
                // Qualities that differ per base, so the bands are not all equal and the
                // enumeration has something to weigh.
                new Read("varied-quality", "rg4", "TTGCAA", "I5+I5+", tp(6)));

        final Path bam = dir.resolve("reads.bam");
        writeBam(bam, reads);
        describe(bam);

        run(dir, "default", bam, List.of());
        run(dir, "replace", bam, List.of("--replace-quality-mode", "true"));
        run(dir, "max-quality", bam, List.of("--maximal-quality-score", "10"));
        run(dir, "min-error-rate", bam, List.of("--minimal-error-rate", "0.1"));
        run(dir, "both", bam,
                List.of("--replace-quality-mode", "true", "--maximal-quality-score", "20"));

        // The same bases under a flow order whose first base repeats immediately. The cycle is
        // one, so the slice is a single flow and the enumeration indexes one past it. It is its
        // own input because the walker dies on the first read that reaches it.
        final Path cycle = dir.resolve("cycle-one.bam");
        writeBam(cycle, List.of(new Read("cycle-one", "rg1", "TGCATGCA", "IIIIIIII", tp(8))));
        describe(cycle);
        run(dir, "cycle-one", cycle, List.of());
    }

    /** A `tp` of all zeros: every base's quality is the probability of the hmer's own length. */
    static byte[] tp(final int length) {
        return new byte[length];
    }

    /**
     * The key and the three raw probabilities the tool bands, per read.
     *
     * These are read off a FlowBasedRead built exactly as the tool builds it, so the port can
     * compute the bands, the floor and the whole enumeration from them without a flow matrix of
     * its own.
     */
    static void describe(final Path bam) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(bam.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMReadGroupRecord group : header.getReadGroups()) {
                System.out.printf("group\t%s\torder=%s\tmaxClass=%s%n",
                        group.getId(),
                        group.getFlowOrder(),
                        new FlowBasedReadUtils.ReadGroupInfo(group).maxClass);
            }
            for (final SAMRecord record : reader) {
                final StringBuilder tpText = new StringBuilder();
                for (final byte value : record.getSignedByteArrayAttribute(
                        FlowBasedRead.FLOW_MATRIX_TAG_NAME)) {
                    tpText.append(tpText.length() == 0 ? "" : ",").append(value);
                }
                System.out.printf("read\t%s\tgroup=%s\tbases=%s\tqual=%s\ttp=%s%n",
                        record.getReadName(),
                        record.getReadGroup().getId(),
                        record.getReadString(),
                        record.getBaseQualityString(),
                        tpText);
                final FlowBasedReadUtils.ReadGroupInfo info = FlowBasedReadUtils.getReadGroupInfo(
                        header, new SAMRecordToGATKReadAdapter(record));
                // The tool sets keepBoundaryFlows before it walks, so the key measured here is
                // the key it computes rather than the default one.
                final FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();
                fbargs.keepBoundaryFlows = true;
                final FlowBasedRead flowRead = new FlowBasedRead(record, info.flowOrder,
                        info.maxClass, fbargs);
                final int[] key = flowRead.getKey();
                final StringBuilder keyText = new StringBuilder();
                for (final int hmer : key) {
                    keyText.append(keyText.length() == 0 ? "" : ",").append(hmer);
                }
                System.out.printf("flow\t%s\tkey=%s\tmaxHmer=%d%n",
                        record.getReadName(), keyText, flowRead.getMaxHmer());
                for (int flow = 0; flow < key.length; flow++) {
                    System.out.printf("prob\t%s\t%d\tminus=%s\tkey=%s\tplus=%s%n",
                            record.getReadName(),
                            flow,
                            key[flow] > 0
                                    ? Double.toString(flowRead.getProb(flow, key[flow] - 1))
                                    : "none",
                            Double.toString(flowRead.getProb(flow, key[flow])),
                            key[flow] < flowRead.getMaxHmer()
                                    ? Double.toString(flowRead.getProb(flow, key[flow] + 1))
                                    : "none");
                }
            }
        }
    }

    static void run(final Path dir, final String label, final Path bam, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".bam");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(), "-O", out.toString()));
        argv.addAll(extra);
        try {
            new AddFlowBaseQuality().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // The frame matters here: the index is computed from the flow order's cycle, so
            // naming where it lands is what makes the refusal a measurement rather than a crash.
            System.out.printf("error\t%s\t%s:%s\tat=%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)),
                    frame(e));
            return;
        }
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(out.toFile())) {
            for (final SAMRecord record : reader) {
                final String xq = record.getStringAttribute(
                        AddFlowBaseQuality.BASE_QUALITY_ATTRIBUTE_NAME);
                if (xq != null) {
                    System.out.printf("xq\t%s\t%s=%s%n", label, record.getReadName(),
                            ReferenceQueryDump.escape(xq));
                }
                final String oq = record.getStringAttribute(
                        AddFlowBaseQuality.OLD_QUALITY_ATTRIBUTE_NAME);
                if (oq != null) {
                    System.out.printf("oq\t%s\t%s=%s%n", label, record.getReadName(),
                            ReferenceQueryDump.escape(oq));
                }
                System.out.printf("qual\t%s\t%s=%s%n", label, record.getReadName(),
                        ReferenceQueryDump.escape(record.getBaseQualityString()));
            }
        }
    }

    static void writeBam(final Path file, final List<Read> reads) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addReadGroup(group("rg4", FOUR));
        header.addReadGroup(group("rg1", ONE));
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, file.toFile())) {
            int start = 100;
            for (final Read read : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(read.name());
                record.setReferenceName("chr1");
                record.setAlignmentStart(start);
                start += 100;
                record.setCigarString(read.bases().length() + "M");
                record.setReadBases(read.bases().getBytes(StandardCharsets.UTF_8));
                record.setBaseQualityString(read.quals());
                record.setMappingQuality(60);
                record.setAttribute("RG", read.group());
                record.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, read.tp());
                writer.addAlignment(record);
            }
        }
    }

    static SAMReadGroupRecord group(final String id, final String flowOrder) {
        final SAMReadGroupRecord group = new SAMReadGroupRecord(id);
        // ULTIMA is the platform that makes a read group a flow read group, and it is refused
        // without a flow order.
        group.setPlatform("ULTIMA");
        group.setFlowOrder(flowOrder);
        group.setSample("sample");
        return group;
    }

    /** The first frame of the throwable that belongs to the tool, method and line. */
    static String frame(final Throwable thrown) {
        Throwable cause = thrown;
        while (cause.getCause() != null) {
            cause = cause.getCause();
        }
        for (final StackTraceElement element : cause.getStackTrace()) {
            if (element.getClassName().startsWith("org.broadinstitute.hellbender")) {
                return element.getClassName() + "." + element.getMethodName()
                        + ":" + element.getLineNumber();
            }
        }
        return "none";
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
