/*
 * AddFlowSNVQuality's qualities, taken from the reference.
 *
 * The sibling of AddFlowBaseQuality: the same enumeration over the ways a flow key could have been
 * misread, but producing a probability for each of the three bases that were NOT called, and then
 * throwing away the base quality it just computed.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE COMPUTED BASE QUALITY IS DISCARDED. After the hmer loop fills baseProbs, the last loop
 *     overwrites EVERY entry with `1 - snvqProbs[calledIndex][ofs]`, which is itself
 *     `1 - sum(alt probabilities)`. The reference's own comment reads "at this point, bq becomes
 *     trivial (?)";
 *   - THE PHRED CONVERSION ROUNDS, `(byte)Math.round(-10 * log10(p))`, where the sibling tool
 *     TRUNCATES the same expression;
 *   - THE SIDE WALK IS BOUNDED BY THE SLICE, `if (sideFlow < minIndex || sideFlow > maxIndex)
 *     break`, which is exactly the bound AddFlowBaseQuality lacks, so a flow order whose cycle is
 *     ONE runs here and throws there;
 *   - THE SNVQ MODE IS A CHOICE OF FOUR FORMULAE over the same two probabilities: Legacy takes the
 *     slice probability, Optimistic their product, Pessimistic one minus the product of the
 *     complements, and Geometric the square root of the two, which is the default;
 *   - THE ALT PROBABILITIES ARE KEYED BY BASE, NOT BY FLOW, so two side flows carrying the same
 *     base overwrite each other in a LinkedHashMap and the LAST one wins;
 *   - EVERY MIDDLE BASE OF AN HMER TAKES THE MINIMUM RATE rather than a computed one, and the
 *     called base's own slot takes `max(0, 1 - altP)`;
 *   - --max-phred-score MOVES BOTH THE CLAMP AND THE FLOOR, since it sets maxQualityScore and
 *     minLikelihoodProbRate together;
 *   - THE OUTPUT ATTRIBUTES ARE NAMED FOR THE FLOW ORDER, lower-cased and prefixed with q, so a
 *     read group ordered TGCA writes qt, qg, qc and qa;
 *   - AND --output-quality-attribute LEAVES QUAL ALONE, writing the base quality into a tag
 *     instead, where the default overwrites QUAL.
 *
 * Output:
 *
 *     read\t<name>\tgroup=<rg>\tbases=<bases>\tqual=<phred string>\ttp=<comma separated>
 *     group\t<rg>\torder=<flow order>\tmaxClass=<n>
 *     flow\t<name>\tkey=<comma separated>\tmaxHmer=<n>
 *     prob\t<name>\t<flow>\tminus=<double>\tkey=<double>\tplus=<double>
 *     qual\t<label>\t<name>=<the quality string after the run>
 *     attr\t<label>\t<name>\t<tag>=<the attribute>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AddFlowSNVQualityDump
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
import org.broadinstitute.hellbender.tools.walkers.featuremapping.AddFlowSNVQuality;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class AddFlowSNVQualityDump {

    /** A four-flow order, whose four bases are the four output attributes. */
    static final String FOUR = "TGCA";
    /**
     * A flow order whose first base repeats immediately. calcFlowOrderLength answers ONE, and this
     * tool survives it where AddFlowBaseQuality throws, because its side walk stops at the slice.
     */
    static final String ONE = "TTGCA";

    record Read(String name, String group, String bases, String quals, byte[] tp) { }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("add-flow-snv-quality-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AddFlowSNVQualityDump: a probability for every base that was not called");

        final List<Read> reads = List.of(
                // Every hmer of length one, so every hmer walks both sides.
                new Read("singles", "rg4", "TGCATGCA", "IIIIIIII", tp(8)),
                // Hmers of length two, the shortest with a computed last base.
                new Read("doubles", "rg4", "TTGGCCAA", "IIIIIIII", tp(8)),
                // A five-base run, whose middle bases take the minimum rate.
                new Read("long-hmer", "rg4", "TTTTTGCA", "IIIIIIII", tp(8)),
                // A read opening on a zero flow.
                new Read("leading-zero", "rg4", "GCATGCA", "IIIIIII", tp(7)),
                // A read ending inside an hmer.
                new Read("trailing-hmer", "rg4", "TGCAAAA", "IIIIIII", tp(7)),
                // Qualities that differ per base, so the bands are not all equal.
                new Read("varied-quality", "rg4", "TTGCAA", "I5+I5+", tp(6)));

        final Path bam = dir.resolve("reads.bam");
        writeBam(bam, reads);
        describe(bam);

        run(dir, "default", bam, List.of());
        run(dir, "legacy", bam, List.of("--snvq-mode", "Legacy"));
        run(dir, "optimistic", bam, List.of("--snvq-mode", "Optimistic"));
        run(dir, "pessimistic", bam, List.of("--snvq-mode", "Pessimistic"));
        // The clamp and the floor move together.
        run(dir, "max-phred-20", bam, List.of("--max-phred-score", "20"));
        run(dir, "max-phred-10", bam, List.of("--max-phred-score", "10"));
        // The base quality written to a tag rather than over QUAL.
        run(dir, "attribute", bam, List.of("--output-quality-attribute", "BQ"));

        // A flow order whose cycle is ONE. The side walk survives it, unlike the sibling tool's,
        // but the normalisation does not: with a cycle of one only the first base of the order is
        // considered, so a read carrying any other base leaves calledIndex at -1 and the array is
        // indexed with it. The guard above it tests `calledBase < 0`, which an ASCII base never
        // is, so it never fires.
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
     * Read off a FlowBasedRead built as the tool builds it, so the port can compute the bands and
     * the whole enumeration from them without a flow matrix of its own.
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
                // This tool leaves the flow arguments at their defaults, unlike its sibling, so
                // the key measured here is the default one.
                final FlowBasedRead flowRead = new FlowBasedRead(record, info.flowOrder,
                        info.maxClass, new FlowBasedArgumentCollection());
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
            new AddFlowSNVQuality().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(out.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            for (final SAMRecord record : reader) {
                System.out.printf("qual\t%s\t%s=%s%n", label, record.getReadName(),
                        ReferenceQueryDump.escape(record.getBaseQualityString()));
                final String order = record.getReadGroup().getFlowOrder();
                final int cycle = FlowBasedReadUtils.calcFlowOrderLength(order);
                for (int i = 0; i < cycle; i++) {
                    final String tag = AddFlowSNVQuality.attrNameForNonCalledBase(order.charAt(i));
                    final String value = record.getStringAttribute(tag);
                    if (value != null) {
                        System.out.printf("attr\t%s\t%s\t%s=%s%n", label, record.getReadName(),
                                tag, ReferenceQueryDump.escape(value));
                    }
                }
                final String bq = record.getStringAttribute("BQ");
                if (bq != null) {
                    System.out.printf("attr\t%s\t%s\tBQ=%s%n", label, record.getReadName(),
                            ReferenceQueryDump.escape(bq));
                }
                // The flow order of the group, so the attribute names are derivable rather than
                // guessed.
                System.out.printf("order\t%s\t%s=%s%n", label, record.getReadName(), order);
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
        group.setPlatform("ULTIMA");
        group.setFlowOrder(flowOrder);
        group.setSample("sample");
        return group;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
