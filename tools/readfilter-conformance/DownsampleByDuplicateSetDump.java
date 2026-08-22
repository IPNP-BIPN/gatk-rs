/*
 * DownsampleByDuplicateSet's output, taken from the reference.
 *
 * Whole molecules dropped rather than reads, so that a mixture keeps its family-size distribution.
 * A hundred and five lines of tool over a hundred and fifty of `DuplicateSetWalker`, and the walker
 * is where the surprises are.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE LAST DUPLICATE SET OF THE FILE ESCAPES EVERY REJECTION RULE. `processLastReadSet` calls
 *     `apply` directly and never asks `rejectSet`, so a trailing set that is too small, or has an
 *     odd number of reads, or has too few reads on a strand, is offered to the tool anyway;
 *   - A SET WITH AN ODD NUMBER OF READS IS REJECTED AT THE DEFAULTS, `size() % 2 == 1` being one of
 *     the three rejection rules, so a three-read molecule is dropped whatever the minimums say;
 *   - A REJECTED SET DOES NOT CONSUME A RANDOM DRAW, the rejection happening before `apply`, so
 *     adding a rejectable molecule to the front of a file changes which of the LATER molecules
 *     survive;
 *   - THE SEED IS FIXED AT 142, so a run is reproducible but the outcome depends entirely on the
 *     ORDER the sets arrive in;
 *   - THE DRAW IS `rng.nextDouble() < fractionToKeep`, so a fraction of 1.0 keeps everything and a
 *     fraction of 0.0 keeps nothing, and both still draw;
 *   - THE SETS ARE CUT ON THE MOLECULE NUMBER OF THE MI TAG, the part before the slash, so
 *     `MI:Z:0/A` and `MI:Z:0/B` are ONE set and the strand suffix only feeds the per-strand
 *     minimum;
 *   - A MOLECULE NUMBER THAT GOES BACKWARDS IS A UserException naming the tag;
 *   - --min-reads AND --min-per-strand-reads EACH REJECT, and the per-strand one counts the two
 *     suffixes separately;
 *   - AND THE WALKER'S OWN READ FILTERS RUN FIRST, so an unmapped read never reaches a set at all
 *     and can change a set's parity.
 *
 * Output:
 *
 *     input\t<label>=<the whole input as text, escaped>
 *     sam\t<label>=<the whole output as text, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: DownsampleByDuplicateSetDump
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
import org.broadinstitute.hellbender.tools.walkers.consensus.DownsampleByDuplicateSet;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DownsampleByDuplicateSetDump {

    /** One read: the molecule id, the strand suffix, and whether it is mapped. */
    record Read(int molecule, String strand, boolean mapped) {}

    static Read read(final int molecule, final String strand) {
        return new Read(molecule, strand, true);
    }

    public static void main(final String[] args) throws Exception {
        // The deflater is pinned, and the class is printed, because every BAM below travels as
        // text and the pin is what makes that text re-derivable.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        System.out.printf("deflater\t%s%n",
                new DeflaterFactory().makeDeflater(5, true).getClass().getName());

        final Path dir = Path.of("downsample-by-duplicate-set-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# DownsampleByDuplicateSetDump: whole molecules dropped, not reads");

        // Ten molecules of two reads each, one read per strand, which is the shape the tool is
        // written for. Every fraction below sees the same order, so the draws line up.
        final List<Read> even = new ArrayList<>();
        for (int molecule = 0; molecule < 10; molecule++) {
            even.add(read(molecule, "A"));
            even.add(read(molecule, "B"));
        }
        run(dir, "keep-all", even, "1.0");
        run(dir, "keep-none", even, "0.0");
        run(dir, "keep-half", even, "0.5");
        run(dir, "keep-most", even, "0.95");

        // The same ten, with a three-read molecule in FRONT. It is rejected for being odd, which
        // costs no draw, so every later molecule keeps the decision it had before.
        final List<Read> odd = new ArrayList<>();
        odd.add(read(-1, "A"));
        odd.add(read(-1, "A"));
        odd.add(read(-1, "B"));
        odd.addAll(even);
        run(dir, "odd-set-in-front", odd, "0.5");

        // And a three-read molecule at the END, which no rejection rule ever sees.
        final List<Read> trailing = new ArrayList<>(even);
        trailing.add(read(10, "A"));
        trailing.add(read(10, "A"));
        trailing.add(read(10, "B"));
        run(dir, "odd-set-at-the-end", trailing, "1.0");

        // A single molecule of one read, which is both odd and below the minimum, and is the last
        // set of the file, so it is written anyway.
        run(dir, "one-read-file", List.of(read(0, "A")), "1.0");

        // The minimums, each rejecting a set the other would keep.
        run(dir, "min-reads-4", even, "1.0", "--min-reads", "4");
        run(dir, "min-per-strand-2", even, "1.0", "--min-per-strand-reads", "2");

        // A molecule number that goes backwards.
        final List<Read> backwards = List.of(
                read(5, "A"), read(5, "B"), read(1, "A"), read(1, "B"));
        run(dir, "unsorted-molecule-ids", backwards, "1.0");

        // An unmapped read inside a molecule, which the walker's filters drop, leaving the set odd.
        final List<Read> withUnmapped = List.of(
                read(0, "A"), read(0, "B"),
                new Read(1, "A", false), read(1, "A"), read(1, "B"),
                read(2, "A"), read(2, "B"));
        run(dir, "unmapped-read-inside", withUnmapped, "1.0");
    }

    static void run(final Path dir, final String label, final List<Read> reads,
                    final String fraction, final String... extra) throws Exception {
        final Path in = dir.resolve(label + ".bam");
        buildBam(in, reads);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(asText(in)));

        final Path out = dir.resolve(label + "-downsampled.bam");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", in.toString(), "-O", out.toString(), "--fraction-to-keep", fraction));
        argv.addAll(Arrays.asList(extra));
        try {
            new DownsampleByDuplicateSet().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("sam\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(asText(out), dir)));
        }
    }

    /** A coordinate-sorted BAM whose reads carry MI tags in the order given. */
    static void buildBam(final Path file, final List<Read> reads) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 10000))));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        group.setLibrary("lib1");
        group.setPlatform("illumina");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);

        try (SAMFileWriter writer =
                     new SAMFileWriterFactory().makeBAMWriter(header, true, file.toFile())) {
            int index = 0;
            for (final Read read : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("r" + index++);
                record.setReadBases("ACGTACGTAC".getBytes());
                record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
                record.setAttribute("RG", "rg1");
                record.setAttribute("MI", read.molecule() + "/" + read.strand());
                if (read.mapped()) {
                    record.setReferenceName("chr1");
                    record.setAlignmentStart(100 + 10 * (read.molecule() + 1));
                    record.setCigarString("10M");
                    record.setMappingQuality(60);
                } else {
                    record.setReadUnmappedFlag(true);
                    record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                    record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    record.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
                }
                writer.addAlignment(record);
            }
        }
    }

    static String asText(final Path bam) {
        final StringBuilder text = new StringBuilder();
        try (SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(new File(bam.toString()))) {
            text.append(reader.getFileHeader().getSAMString());
            for (final SAMRecord record : reader) {
                text.append(PrintReadsDump.samLine(record));
            }
        } catch (final Exception e) {
            text.append("error: ").append(e);
        }
        return text.toString();
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>")
                .replaceAll("@PG[^\n]*\n", "");
    }
}
