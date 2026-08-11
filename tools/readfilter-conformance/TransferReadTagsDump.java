/*
 * TransferReadTags, taken from the reference.
 *
 * The tenth whole tool of the record-transform archetype, and the FIRST THAT IS NOT A WALKER. It
 * extends GATKTool and overrides traverse() itself, so there is no filter chain, no read
 * transformer and no interval bound: the reads it sees are the data source's, unfiltered. It is
 * also the first with TWO read inputs, walked in lockstep by query name.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE TRAVERSAL IS THE TOOL'S OWN, so none of the archetype's usual machinery runs. A port that
 *     reached for the read walker would apply WellformedReadFilter and drop reads this tool keeps.
 *     The fixture carries a read with no read group, which the wellformed filter rejects and this
 *     tool writes;
 *   - EVERY TAG IS TRANSFERRED AS A STRING, whatever it was. `getAttributeAsString` then
 *     `setAttribute(name, String)`, so an integer `XI:i:42` in the unmapped file arrives as
 *     `XI:Z:42` in the output and a float `XN:f:1.5` as `XN:Z:1.5`. The type is lost, and it is lost
 *     the same way every time;
 *   - THE ALIGNED FILE MUST BE A SUBSET OF THE UNMAPPED ONE, and the catch-up is one-directional.
 *     When the aligned read is ahead, the unmapped iterator is played forward; when it is behind,
 *     the tool throws. Both throw sites carry the same message and only one of them is reachable
 *     from the outer loop;
 *   - AN ALIGNED READ PAST THE END OF THE UNMAPPED FILE IS SILENTLY DROPPED. The catch-up loop is
 *     `while (unmappedSamIterator.hasNext())`, so when it runs out the loop simply ends, nothing is
 *     written, and the outer loop moves on. No exception, no warning, one fewer read in the output.
 *     This is the finding the dump exists for;
 *   - THE WRITER IS NOT TOLD THE READS ARE SORTED. `createSAMWriter(..., false)` on a queryname
 *     header means htsjdk sorts what it is handed with SAMRecordQueryNameComparator, which has six
 *     tie-breaks after the name. The fixture carries a name whose two records are in the file in
 *     the order that comparator would swap;
 *   - THE SORT ORDER IS CHECKED, AND ONLY ON THE ALIGNED FILE. A coordinate-sorted `-I` is refused;
 *     the unmapped file is not checked at all, because "the SortOrder field is often not
 *     populated", and the traversal is relied on to notice;
 *   - THREE MORE REFUSALS, EACH FROM A DIFFERENT LAYER: `--read-tags` omitted is refused by
 *     Barclay before the tool is built at all, and never reaches the `Utils.nonEmpty` in
 *     onTraversalStart, because a List argument with no `optional = true` is required; an unmapped
 *     file that is empty while the aligned one is not is a UserException; and a matched unmapped
 *     read that does not carry the tag asked for is an IllegalArgumentException out of
 *     `Utils.nonNull`, naming the unmapped read rather than the aligned one.
 *
 * Neither input is opened with an index: a queryname-sorted BAM cannot have one. NO OUTPUT CARRIES
 * ONE EITHER, whatever `--create-output-bam-index` says, for the same reason: measured absent on
 * all six runs that finish.
 *
 * Output:
 *
 *     deflater\t<class>
 *     fixture\t<label>\t<base64 bam>
 *     header\t<label>\t<escaped SAM header>
 *     commandline\t<label>\t<@PG command line>
 *     output\t<label>\t<base64 bam>
 *     index\t<label>\t<base64 bai or absent>
 *     reads\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<tags>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: TransferReadTagsDump
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
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import org.broadinstitute.hellbender.tools.walkers.qc.TransferReadTags;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class TransferReadTagsDump {

    public static void main(final String[] args) throws Exception {
        // The factory is static and the first writer wins. This dump calls no Picard entry point,
        // so nothing should replace it; the pin makes that a fact rather than a hope.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("transferreadtags-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# TransferReadTagsDump: TransferReadTags");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        // The unmapped side: six names, every tag present, three different tag types.
        final Path unmapped = dir.resolve("unmapped.bam");
        buildUnmapped(unmapped.toFile(), new String[] {"a1", "a2", "a3", "a4", "a5", "a6"}, null);
        fixture(unmapped, "unmapped");

        // The same, with a3 missing the RX tag the runs ask for.
        final Path unmappedMissing = dir.resolve("unmapped_missing.bam");
        buildUnmapped(unmappedMissing.toFile(),
                new String[] {"a1", "a2", "a3", "a4", "a5", "a6"}, "a3");
        fixture(unmappedMissing, "unmapped_missing");

        // The same, with a3 absent entirely, which is what makes the inner throw reachable.
        final Path unmappedGap = dir.resolve("unmapped_gap.bam");
        buildUnmapped(unmappedGap.toFile(), new String[] {"a1", "a2", "a4", "a5", "a6"}, null);
        fixture(unmappedGap, "unmapped_gap");

        final Path unmappedEmpty = dir.resolve("unmapped_empty.bam");
        buildUnmapped(unmappedEmpty.toFile(), new String[] {}, null);
        fixture(unmappedEmpty, "unmapped_empty");

        // The aligned side: a strict subset, so the catch-up loop runs on every read but the first.
        final Path aligned = dir.resolve("aligned.bam");
        buildAligned(aligned.toFile(), new String[] {"a1", "a3", "a5"},
                SAMFileHeader.SortOrder.queryname, false, true);
        fixture(aligned, "aligned");

        // A name past the end of the unmapped file: the catch-up loop runs out and the read is
        // silently dropped.
        final Path alignedTail = dir.resolve("aligned_tail.bam");
        buildAligned(alignedTail.toFile(), new String[] {"a1", "a9"},
                SAMFileHeader.SortOrder.queryname, false, true);
        fixture(alignedTail, "aligned_tail");

        // A name before the first unmapped one: the outer comparison is negative and the tool
        // throws before the catch-up loop is ever entered.
        final Path alignedBefore = dir.resolve("aligned_before.bam");
        buildAligned(alignedBefore.toFile(), new String[] {"a0", "a1"},
                SAMFileHeader.SortOrder.queryname, false, true);
        fixture(alignedBefore, "aligned_before");

        // Coordinate sorted, which is the one sort order that is checked.
        final Path alignedCoordinate = dir.resolve("aligned_coordinate.bam");
        buildAligned(alignedCoordinate.toFile(), new String[] {"a1", "a3", "a5"},
                SAMFileHeader.SortOrder.coordinate, false, true);
        fixture(alignedCoordinate, "aligned_coordinate");

        final Path alignedEmpty = dir.resolve("aligned_empty.bam");
        buildAligned(alignedEmpty.toFile(), new String[] {},
                SAMFileHeader.SortOrder.queryname, false, true);
        fixture(alignedEmpty, "aligned_empty");

        // One name whose two records are in the file in the order the queryname comparator would
        // swap, and one read with no read group at all.
        final Path alignedUnsorted = dir.resolve("aligned_unsorted.bam");
        buildAligned(alignedUnsorted.toFile(), new String[] {"a1", "a3", "a5"},
                SAMFileHeader.SortOrder.queryname, true, false);
        fixture(alignedUnsorted, "aligned_unsorted");

        transfer(dir, aligned, unmapped, "rx", new String[] {"--read-tags", "RX"});
        transfer(dir, aligned, unmapped, "alltypes",
                new String[] {"--read-tags", "RX", "--read-tags", "XI", "--read-tags", "XN"});
        transfer(dir, alignedTail, unmapped, "tail", new String[] {"--read-tags", "RX"});
        transfer(dir, alignedUnsorted, unmapped, "unsorted", new String[] {"--read-tags", "RX"});
        transfer(dir, alignedEmpty, unmapped, "emptyaligned", new String[] {"--read-tags", "RX"});
        transfer(dir, alignedEmpty, unmappedEmpty, "bothempty", new String[] {"--read-tags", "RX"});

        transfer(dir, aligned, unmappedGap, "gap", new String[] {"--read-tags", "RX"});
        transfer(dir, alignedBefore, unmapped, "before", new String[] {"--read-tags", "RX"});
        transfer(dir, aligned, unmappedMissing, "missingtag", new String[] {"--read-tags", "RX"});
        transfer(dir, alignedCoordinate, unmapped, "coordinate",
                new String[] {"--read-tags", "RX"});
        transfer(dir, aligned, unmappedEmpty, "emptyunmapped", new String[] {"--read-tags", "RX"});
        transfer(dir, aligned, unmapped, "notags", new String[] {});

        for (final String label : new String[] {
                "rx", "alltypes", "tail", "unsorted", "emptyaligned", "bothempty"}) {
            reads(dir, label);
        }
        // What the inputs carried, so the transfer can be read as a difference rather than a state.
        readsOf(unmapped, "in:unmapped");
        readsOf(aligned, "in:aligned");
        readsOf(alignedUnsorted, "in:aligned_unsorted");
    }

    /** Every read of an output, with the tags that are the whole point of the tool. */
    static void reads(final Path dir, final String label) throws Exception {
        readsOf(dir.resolve("TransferReadTags." + label + ".bam"), label);
    }

    static void readsOf(final Path bam, final String label) throws Exception {
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(bam.toFile())) {
            for (final SAMRecord record : reader) {
                final StringBuilder tags = new StringBuilder();
                for (final SAMRecord.SAMTagAndValue tag : record.getAttributes()) {
                    if (tags.length() > 0) {
                        tags.append(';');
                    }
                    // The class is printed with the value: an integer transferred as a string is
                    // the same characters and a different type, and only the type says so.
                    tags.append(tag.tag).append('=').append(tag.value)
                            .append(':').append(tag.value.getClass().getSimpleName());
                }
                System.out.printf("reads\t%s\t%s\t%d\t%s\t%d\t%s%n", label, record.getReadName(),
                        record.getFlags(), record.getReferenceName(), record.getAlignmentStart(),
                        tags.length() == 0 ? "-" : tags.toString());
            }
        }
    }

    /** A fixture. No index travels: a queryname-sorted BAM cannot have one. */
    static void fixture(final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, base64(bam));
    }

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 100))));
        header.setSortOrder(order);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);
        return header;
    }

    /**
     * The unmapped side: one record per name, three tag types on each.
     *
     * `withoutTag` names the one record that is missing `RX`, which is what the tool's
     * `Utils.nonNull` is measured on.
     */
    static void buildUnmapped(final File file, final String[] names, final String withoutTag) {
        final SAMFileHeader header = header(SAMFileHeader.SortOrder.queryname);
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, file)) {
            for (final String name : names) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(name);
                record.setReadUnmappedFlag(true);
                record.setReadBases("ACGTACGTAC".getBytes());
                final byte[] quals = new byte[10];
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                record.setAttribute("RG", "rg1");
                if (!name.equals(withoutTag)) {
                    record.setAttribute("RX", "AAA-CCC");
                }
                // An integer and a float, so the string conversion is visible on both.
                record.setAttribute("XI", 42);
                record.setAttribute("XN", 1.5f);
                writer.addAlignment(record);
            }
        }
    }

    /**
     * The aligned side.
     *
     * `swapPairOrder` writes the second-of-pair record of the first name before its first-of-pair
     * record, which is an order the queryname comparator swaps back; `readGroups` false drops the
     * `RG` tag, which the wellformed filter would reject and this tool does not run.
     */
    static void buildAligned(final File file, final String[] names,
                             final SAMFileHeader.SortOrder order, final boolean swapPairOrder,
                             final boolean readGroups) {
        final SAMFileHeader header = header(order);
        int start = 1;
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, file)) {
            for (final String name : names) {
                if (swapPairOrder && name.equals(names[0])) {
                    // Second of pair first, which is the order the comparator undoes.
                    writer.addAlignment(aligned(header, name, start, 0x1 | 0x80, readGroups));
                    writer.addAlignment(aligned(header, name, start, 0x1 | 0x40, readGroups));
                } else {
                    writer.addAlignment(aligned(header, name, start, 0, readGroups));
                }
                start += 10;
            }
        }
    }

    static SAMRecord aligned(final SAMFileHeader header, final String name, final int start,
                             final int flags, final boolean readGroup) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes());
        final byte[] quals = new byte[10];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setFlags(flags);
        if (readGroup) {
            record.setAttribute("RG", "rg1");
        }
        return record;
    }

    static void transfer(final Path dir, final Path aligned, final Path unmapped,
                         final String label, final String[] extra) throws Exception {
        final Path output = dir.resolve("TransferReadTags." + label + ".bam");
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", aligned.toString(), "--unmapped-sam", unmapped.toString(),
                "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new TransferReadTags().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            // A refusal is the observable behaviour, so it is dumped rather than swallowed.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()).replace('\n', ' '));
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
        System.out.printf("output\t%s\t%s%n", label, base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? base64(index) : "absent");
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }

    static String base64(final Path path) throws Exception {
        return Base64.getEncoder().encodeToString(Files.readAllBytes(path));
    }
}
