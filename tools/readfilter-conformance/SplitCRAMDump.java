/*
 * SplitCRAM's output, taken from the reference.
 *
 * A CRAM cut into shards at container boundaries, without decoding a single record. The tool reads
 * containers, counts the records each one declares, and starts a new output once the count reaches
 * the threshold, so what is measurable is where the cuts fall and what the outputs are called.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE CUT IS AT A CONTAINER BOUNDARY AND THE THRESHOLD IS A MINIMUM, NOT A MAXIMUM: the inner
 *     loop tests `records < shardRecords` BEFORE reading a container and adds its whole record
 *     count after, so a shard overshoots by up to one container. The test being strict, a threshold
 *     exactly the size of a container still gives one container per shard, and one above it gives
 *     two;
 *   - EVERY SHARD IS A WHOLE CRAM: the CRAM header, a SAM header container and an EOF container are
 *     written around each group, so a shard is readable on its own and every shard repeats the
 *     input's header;
 *   - THE OUTPUT NAME IS String.format ON A TEMPLATE, with a counter that starts at zero and is
 *     incremented whether or not the shard turns out to be written;
 *   - A TEMPLATE WITH NO %d IS REFUSED IN onStartup, by an IllegalArgumentException rather than a
 *     UserException, and the pattern it is checked against is `%[0-9]*d`, so `%04d` passes and a
 *     width with a flag such as `%-4d` does not;
 *   - --shard-max-output-count DOES NOT LIMIT ANYTHING ABOVE ONE: the counter it is compared with
 *     is declared inside the outer loop, so it is reset to zero for every shard and is always one
 *     when tested, and only the value 1 ever stops the run;
 *   - AND AN EMPTY CRAM PRODUCES NO SHARD AT ALL, the outer loop never running.
 *
 * Output:
 *
 *     fixture\t<label>\tcontainers=<records per container, comma separated>
 *     shard\t<label>\t<file name>\trecords=<n>\tnames=<read names, comma separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SplitCRAMDump
 */

import htsjdk.samtools.CRAMContainerStreamWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CramContainerIterator;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.SplitCRAM;

import java.io.BufferedInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class SplitCRAMDump {

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("split-cram-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SplitCRAMDump: a CRAM cut into shards at container boundaries");

        // Four contigs, so the CRAM writer closes a container at every reference change and the
        // input has containers small enough to shard.
        final Path dict = MultiFeatureWalkerDump.writeDictionary(dir, "split",
                List.of("chr1", "chr2", "chr3", "chr4"));
        final Path fasta = dir.resolve("split.fasta");
        final Path cram = dir.resolve("reads.cram");
        writeCram(cram, fasta, dict);
        System.out.printf("fixture\tinput\tcontainers=%s%n", containerRecords(cram));

        // One container per shard, which is what a threshold of one gives.
        run(dir, "one-per-shard", cram, fasta, "shard_%d.cram", "--shard-records", "1");
        // A threshold exactly the size of a container, which is still one container per shard:
        // the test is strict, and a container that reaches the threshold ends the shard.
        run(dir, "exact-threshold", cram, fasta, "exact_%d.cram", "--shard-records", "3");
        // A threshold one above it, which the first container does not reach, so the shard takes a
        // second and overshoots.
        run(dir, "overshoot", cram, fasta, "pair_%d.cram", "--shard-records", "4");
        // A threshold no shard reaches, which is one shard holding everything.
        run(dir, "one-shard", cram, fasta, "all_%d.cram", "--shard-records", "1000");
        // The maximum output count, at every value that could stop the run.
        run(dir, "max-one", cram, fasta, "max1_%d.cram",
                "--shard-records", "1", "--shard-max-output-count", "1");
        run(dir, "max-two", cram, fasta, "max2_%d.cram",
                "--shard-records", "1", "--shard-max-output-count", "2");
        run(dir, "max-three", cram, fasta, "max3_%d.cram",
                "--shard-records", "1", "--shard-max-output-count", "3");
        // A padded template, which the pattern accepts, and two it does not.
        run(dir, "padded-template", cram, fasta, "padded_%04d.cram", "--shard-records", "1000");
        run(dir, "flagged-template", cram, fasta, "flagged_%-4d.cram", "--shard-records", "1000");
        run(dir, "no-formatter", cram, fasta, "plain.cram", "--shard-records", "1000");

        // A CRAM with no records at all, whose outer loop never runs.
        final Path empty = dir.resolve("empty.cram");
        writeEmptyCram(empty, fasta, dict);
        System.out.printf("fixture\tempty\tcontainers=%s%n", containerRecords(empty));
        run(dir, "empty-input", empty, fasta, "empty_%d.cram", "--shard-records", "1");
    }

    /** The record count each container of a CRAM declares, which is what the tool counts with. */
    static String containerRecords(final Path cram) throws Exception {
        final List<String> counts = new ArrayList<>();
        try (final CramContainerIterator containers =
                     new CramContainerIterator(new BufferedInputStream(Files.newInputStream(cram)))) {
            while (containers.hasNext()) {
                final Container container = containers.next();
                counts.add(Integer.toString(container.getContainerHeader().getNumberOfRecords()));
            }
        }
        return String.join(",", counts);
    }

    /** Three reads on each of four contigs, plus three unmapped, coordinate sorted.
     *
     * Written through CRAMContainerStreamWriter with three reads to a slice and one slice to a
     * container, because the default is ten thousand and one container would not be a fixture. */
    static void writeCram(final Path cram, final Path fasta, final Path dict) throws Exception {
        final SAMFileHeader header = readHeader(dict);
        final CRAMEncodingStrategy strategy = new CRAMEncodingStrategy()
                .setMinimumSingleReferenceSliceSize(1)
                .setReadsPerSlice(3)
                .setSlicesPerContainer(1);
        try (final java.io.OutputStream out = Files.newOutputStream(cram)) {
            final CRAMContainerStreamWriter writer = new CRAMContainerStreamWriter(strategy,
                    new ReferenceSource(fasta), header, out, null, cram.getFileName().toString());
            writer.writeHeader();
            for (final String contig : List.of("chr1", "chr2", "chr3", "chr4")) {
                for (int index = 0; index < 3; index++) {
                    writer.writeAlignment(read(header, contig + "-" + index, contig, 1 + index * 10));
                }
            }
            for (int index = 0; index < 3; index++) {
                writer.writeAlignment(read(header, "unmapped-" + index, null, 0));
            }
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

    static SAMFileHeader readHeader(final Path dict) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(dict.toFile())) {
            final SAMFileHeader header = reader.getFileHeader();
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
            return header;
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadBases("ACGTACGTAC".getBytes());
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

    static void run(final Path dir, final String label, final Path cram, final Path fasta,
                    final String template, final String... extra) throws Exception {
        final Path work = dir.resolve(label);
        Files.createDirectories(work);
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", cram.toString(),
                "-O", work.resolve(template).toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new SplitCRAM().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
        }
        try (final Stream<Path> written = Files.list(work)) {
            final List<Path> shards = written
                    .sorted(Comparator.comparing(path -> path.getFileName().toString()))
                    .toList();
            for (final Path shard : shards) {
                System.out.printf("shard\t%s\t%s\trecords=%s\tnames=%s%n", label,
                        shard.getFileName(), containerRecords(shard), readNames(shard, fasta));
            }
        }
    }

    /** Every read name of a shard, in order, which says where the cut fell. */
    static String readNames(final Path shard, final Path fasta) throws Exception {
        final List<String> names = new ArrayList<>();
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(fasta)
                .validationStringency(ValidationStringency.SILENT)
                .open(shard.toFile())) {
            for (final SAMRecord record : reader) {
                names.add(record.getReadName());
            }
        }
        return String.join(",", names);
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
