/*
 * SetNmMdAndUqTags's output, taken from the reference.
 *
 * NM, MD and UQ recalculated against the reference for every mapped record. `SetNmAndUqTags` is the
 * deprecated subclass of this and adds nothing, so measuring this one measures both.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE TAGS ARE REPLACED, NOT FILLED IN: a record arriving with a wrong NM, a wrong MD and a
 *     wrong UQ comes out with all three corrected, and one arriving with none gains all three;
 *   - AN UNMAPPED READ IS LEFT ALONE ENTIRELY, `fixRecord` returning before it touches anything, so
 *     a wrong tag on an unmapped read SURVIVES;
 *   - UQ IS ONLY SET WHEN THE READ HAS QUALITIES: `fixUq` checks `NULL_QUALS`, so a record whose
 *     qualities are `*` keeps whatever UQ it arrived with, while its NM and MD are still fixed;
 *   - SET_ONLY_UQ SKIPS NM AND MD, so a record with a wrong NM keeps it while its UQ is corrected;
 *   - UQ IS A SUM OF BASE QUALITIES AT MISMATCHES, so two mismatches of quality 30 make 60 and a
 *     read that matches perfectly gets UQ=0 rather than no tag at all;
 *   - THE MISMATCH TEST IS IUPAC-AWARE: an `N` in the reference against any read base is a
 *     mismatch, while an `N` in the READ against an `N` in the reference is a match;
 *   - AN INSERTION IS NOT COMPARED AND A DELETION IS, so a spliced or indel-carrying cigar changes
 *     NM by the indel lengths and MD by the deleted bases with their `^`;
 *   - IS_BISULFITE_SEQUENCE CHANGES BOTH TAGS: a C in the reference read as T is no longer a
 *     mismatch, so NM falls, and the NM is then recomputed by a DIFFERENT function than the one
 *     that wrote MD, which is why MD can disagree with NM on such a read;
 *   - AND A QUERYNAME SORTED INPUT IS REFUSED, the message naming the order it found.
 *
 * Output:
 *
 *     fasta\t<the reference, escaped>
 *     fixture\t<label>\t<the input BAM, base64>
 *     output\t<label>\t<the fixed BAM, base64>
 *     sam\t<label>=<the fixed BAM as text, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SetNmMdAndUqTagsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import picard.sam.SetNmMdAndUqTags;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SetNmMdAndUqTagsDump {

    /** Ten bases a line, with an N in the middle of chr1 and a C-rich run for the bisulfite run. */
    static final String FASTA =
            ">chr1\n"
            + "ACGTACGTAC\n"
            + "GTACGTNCGT\n"
            + "ACGTACGTAC\n"
            + ">chr2\n"
            + "CCCCCCCCCC\n"
            + "CCCCCCCCCC\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("set-nm-md-uq-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.writeString(fasta, FASTA, StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        System.out.println("# SetNmMdAndUqTagsDump: NM, MD and UQ recalculated against the reference");
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(FASTA));

        final Path sorted = dir.resolve("sorted.bam");
        buildBam(sorted, SAMFileHeader.SortOrder.coordinate);
        System.out.printf("fixture\tsorted\t%s%n", RecordTransformDump.base64(sorted));

        final Path queryname = dir.resolve("queryname.bam");
        buildBam(queryname, SAMFileHeader.SortOrder.queryname);
        System.out.printf("fixture\tqueryname\t%s%n", RecordTransformDump.base64(queryname));

        run(dir, "defaults", sorted, fasta);
        run(dir, "only-uq", sorted, fasta, "SET_ONLY_UQ=true");
        run(dir, "bisulfite", sorted, fasta, "IS_BISULFITE_SEQUENCE=true");
        run(dir, "queryname", queryname, fasta);
    }

    /**
     * One record per behaviour, in coordinate order.
     *
     * The reference is `ACGTACGTAC GTACGTNCGT ACGTACGTAC` on chr1 and twenty Cs on chr2.
     */
    static void buildBam(final Path file, final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 30),
                new SAMSequenceRecord("chr2", 20))));
        header.setSortOrder(order);

        // The queryname fixture holds the same records; the writer sorts them, since they are
        // written in coordinate order.
        final boolean presorted = order == SAMFileHeader.SortOrder.coordinate;
        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, presorted, file.toFile())) {
            // A perfect match arriving with no tags at all: gains NM=0, MD=8 and UQ=0.
            writer.addAlignment(read(header, "perfect", "chr1", 1, "8M", "ACGTACGT",
                    quals(8, 30), null, null, null));
            // Two mismatches, both quality 30, arriving with tags that are all wrong.
            writer.addAlignment(read(header, "wrong-tags", "chr1", 1, "8M", "AAGTACCT",
                    quals(8, 30), 99, "wrong", 12345));
            // A read over the reference's N, which is a mismatch whatever the read base is, and a
            // read whose own base is N there, which is a match.
            writer.addAlignment(read(header, "over-n", "chr1", 13, "6M", "ACGTAC",
                    quals(6, 20), null, null, null));
            writer.addAlignment(read(header, "n-in-read", "chr1", 13, "6M", "ACGTNC",
                    quals(6, 20), null, null, null));
            // An insertion and a deletion, which move NM by their lengths and MD by the deleted
            // bases.
            writer.addAlignment(read(header, "indels", "chr1", 21, "3M2I3M2D2M", "ACGTTACGTA",
                    quals(10, 25), null, null, null));
            // A read with no qualities at all, which keeps the UQ it arrived with.
            writer.addAlignment(read(header, "no-quals", "chr1", 21, "4M", "ACGA",
                    SAMRecord.NULL_QUALS, null, null, 777));
            // A C-rich read on chr2 where every C was read as a T, which bisulfite treatment
            // forgives and the default run does not.
            writer.addAlignment(read(header, "bisulfite", "chr2", 1, "8M", "TTTTTTTT",
                    quals(8, 40), null, null, null));
            // And an unmapped read carrying wrong tags, which nothing touches.
            final SAMRecord unmapped = new SAMRecord(header);
            unmapped.setReadName("unmapped");
            unmapped.setReadUnmappedFlag(true);
            unmapped.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            unmapped.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            unmapped.setReadBases("ACGTACGT".getBytes(StandardCharsets.UTF_8));
            unmapped.setBaseQualities(quals(8, 30));
            unmapped.setAttribute("NM", 42);
            unmapped.setAttribute("MD", "nonsense");
            unmapped.setAttribute("UQ", 4242);
            writer.addAlignment(unmapped);
        }
    }

    static byte[] quals(final int length, final int value) {
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) value);
        return qualities;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start, final String cigar, final String bases,
                          final byte[] qualities, final Integer nm, final String md,
                          final Integer uq) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        if (nm != null) {
            record.setAttribute("NM", nm);
        }
        if (md != null) {
            record.setAttribute("MD", md);
        }
        if (uq != null) {
            record.setAttribute("UQ", uq);
        }
        return record;
    }

    static void run(final Path dir, final String label, final Path input, final Path fasta,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("fixed-" + label + ".bam");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "I=" + input, "O=" + out, "R=" + fasta));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new SetNmMdAndUqTags().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s\t%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(out));
        System.out.printf("sam\t%s=%s%n", label, ReferenceQueryDump.escape(asText(out)));
    }

    /** The records as text, so a divergence reads as a tag rather than as a byte offset. */
    static String asText(final Path bam) {
        final StringBuilder text = new StringBuilder();
        try (SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(htsjdk.samtools.ValidationStringency.SILENT)
                .open(new File(bam.toString()))) {
            for (final SAMRecord record : reader) {
                text.append(PrintReadsDump.samLine(record));
            }
        } catch (final Exception e) {
            text.append("error: ").append(e);
        }
        return text.toString();
    }
}
