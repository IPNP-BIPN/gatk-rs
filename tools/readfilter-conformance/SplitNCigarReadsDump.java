/*
 * SplitNCigarReads' output, taken from the reference.
 *
 * The tool the two previous suites were built for: it splits a read at every N of its cigar, hands
 * the family to OverhangFixingManager, and writes what comes back. It is a MultiplePassReadWalker,
 * so the input is traversed TWICE: the first pass records which reads the clipper moved, the second
 * one writes them out with their mates' information repaired.
 *
 * Ten behaviours this is built to catch.
 *
 *   - A READ WITH k N ELEMENTS BECOMES k+1 READS, each keeping the bases of one section and soft
 *     clipping the rest, so a 3M2N5M read comes out as 3M5S and 3S5M. The bases and qualities are
 *     the whole read's in every piece: nothing is removed, only clipped;
 *   - A SECTION THAT ENDS OR STARTS WITH D IS TRIMMED FIRST, so the first section of 1M1D1N1M is
 *     1M and not 1M1D, and a read whose sections would be empty between two Ns is an
 *     IllegalArgumentException naming the cigar;
 *   - A CIGAR THAT ENDS IN N EMITS ONE PIECE AND LOSES THE N: `8M2N` comes out `8M`, not split and
 *     not left alone. A CIGAR THAT BEGINS WITH N IS PASSED THROUGH UNTOUCHED, `2N8M` and all,
 *     because the leading N produces no section and the tool returns the read it was given;
 *   - N-D-N IS REFUSED unless --refactor-cigar-string is given, which merges the three elements into
 *     one N of their total length BEFORE the read reaches the filters;
 *   - THE MAPPING QUALITY TRANSFORM IS ON BY DEFAULT: 255 becomes 60, and only 255. It runs AFTER
 *     the read filters, so a read the filters dropped never sees it;
 *   - ITS DEFAULT READ FILTER IS ALLOW_ALL_READS, so a malformed read is split and written like any
 *     other. That is the same pattern UnmarkDuplicates takes and not the engine's;
 *   - AN MC TAG IS REWRITTEN TO WHAT THE MATE'S CIGAR WILL BECOME, computed by running the split on
 *     an artificial read carrying that cigar;
 *   - A SECONDARY ALIGNMENT IS NOT SPLIT unless --process-secondary-alignments is given: it is
 *     passed to the manager whole and written as it arrived, with its mate information still
 *     repaired;
 *   - THE FAMILY IS MARKED SUPPLEMENTARY AND GIVEN SA TAGS on the way out, and NM, MD and NH are
 *     cleared from every piece;
 *   - AND ITS WRITER IS NOT PRESORTED, `createSAMWriter(OUTPUT, false)`, because the manager does
 *     not guarantee it hands reads over in coordinate order.
 *
 * Output, one row per (label, kind):
 *
 *     reference\t<the reference bases of chr1>
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     header\t<label>\t<the output header, escaped>
 *     commandline\t<label>\t<the @PG command line>
 *     output\t<label>\t<the output BAM, base64>
 *     index\t<label>\t<the index, base64, or absent>
 *     samline\t<label>\t<one output record, escaped>
 *     error\t<label>\t<exception class>:<message>
 *     ndn\t<label>\t<the cigar refactorNDNtoN returns>
 *
 * Usage: SplitNCigarReadsDump
 */

import htsjdk.samtools.Cigar;
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
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.hellbender.tools.walkers.rnaseq.SplitNCigarReads;
import org.broadinstitute.hellbender.transformers.NDNCigarReadTransformer;

import java.io.File;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SplitNCigarReadsDump {

    /** 120 bases, non-repeating, so a clip lands where the cigar says and not where a repeat does. */
    static final String CHR1 =
            "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGTGATTACAGGCTCTAGCATCGATCGATCGATTAGCTAGCTAGCTAACCGGTTACGTAGGCTTACCGGATCGATCGATCGAT";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("splitncigar-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SplitNCigarReadsDump: SplitNCigarReads' output, from the reference");
        System.out.printf("reference\t%s%n", CHR1);

        final Path fasta = writeReference(dir);

        // The cigar refactor on its own, which runs before the filters and can be measured alone.
        final Method refactor = NDNCigarReadTransformer.class.getDeclaredMethod(
                "refactorNDNtoN", Cigar.class);
        refactor.setAccessible(true);
        for (final String cigar : new String[] {
                "3M2N1D2N4M",   // the motif, merged into one N of five
                "3M2N4M",       // one N, left alone
                "3M2N1D2N1D2N2M", // two motifs in a row
                "3M2N1I2N4M",   // N-I-N, which is not the motif
                "2N3M",         // a leading N with nothing after it to merge
                "3M2N",         // a trailing N, which has fewer than two elements after it
                "10M"}) {
            System.out.printf("ndn\t%s\t%s%n", cigar,
                    refactor.invoke(null, TextCigarCodec.decode(cigar)));
        }

        // Every shape of split, in one file: no N, one N, two Ns, D beside N, and a leading N.
        final Path splits = dir.resolve("splits.bam");
        buildSplits(splits.toFile());
        emitFixture(dir, splits, "splits");

        // Mapping qualities, including the 255 the transform rewrites.
        final Path qualities = dir.resolve("qualities.bam");
        buildQualities(qualities.toFile());
        emitFixture(dir, qualities, "qualities");

        // A pair whose MC tags name each other's cigars, and a secondary alignment.
        final Path pairs = dir.resolve("pairs.bam");
        buildPairs(pairs.toFile());
        emitFixture(dir, pairs, "pairs");

        // One read whose N leaves a mismatching overhang across the splice of another.
        final Path overhangs = dir.resolve("overhangs.bam");
        buildOverhangs(overhangs.toFile());
        emitFixture(dir, overhangs, "overhangs");

        // A read whose cigar is N-D-N, which is refused without the refactor argument.
        final Path ndn = dir.resolve("ndn.bam");
        buildNdn(ndn.toFile());
        emitFixture(dir, ndn, "ndn");

        run(dir, splits, fasta, "splits", new String[] {});
        run(dir, qualities, fasta, "qualities", new String[] {});
        run(dir, qualities, fasta, "qualities-skip-mq",
                new String[] {"--skip-mapping-quality-transform", "true"});
        run(dir, pairs, fasta, "pairs", new String[] {});
        run(dir, pairs, fasta, "pairs-secondary",
                new String[] {"--process-secondary-alignments", "true"});
        run(dir, overhangs, fasta, "overhangs", new String[] {});
        run(dir, overhangs, fasta, "overhangs-not-fixed",
                new String[] {"--do-not-fix-overhangs", "true"});
        run(dir, overhangs, fasta, "overhangs-strict",
                new String[] {"--max-mismatches-in-overhang", "0"});
        run(dir, ndn, fasta, "ndn", new String[] {});
        run(dir, ndn, fasta, "ndn-refactored",
                new String[] {"--refactor-cigar-string", "true"});
    }

    static void emitFixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    /** Every shape of split the tool has a branch for. */
    static void buildSplits(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // No N at all: the read goes through untouched.
            writer.addAlignment(read(header, "no-n", 10, "10M", 0, 60));
            // One N: two pieces, 3M5S and 3S5M.
            writer.addAlignment(read(header, "one-n", 20, "3M2N5M", 0, 60));
            // Two Ns: three pieces.
            writer.addAlignment(read(header, "two-n", 30, "2M2N3M2N3M", 0, 60));
            // A deletion beside an N, whose section is trimmed back before the clip.
            writer.addAlignment(read(header, "d-before-n", 40, "3M1D2N5M", 0, 60));
            writer.addAlignment(read(header, "d-after-n", 50, "3M2N1D5M", 0, 60));
            // A cigar that ends in N: the last section is not emitted.
            writer.addAlignment(read(header, "trailing-n", 60, "8M2N", 0, 60));
            // A cigar that begins with N: nothing before it to emit.
            writer.addAlignment(read(header, "leading-n", 70, "2N8M", 0, 60));
            // Soft clips around an N, which travel with their piece.
            writer.addAlignment(read(header, "soft-clipped", 80, "2S3M2N3M2S", 0, 60));
            // An insertion inside a section, which consumes read bases and no reference.
            writer.addAlignment(read(header, "insertion", 90, "3M2I3M2N2M", 0, 60));
        }
    }

    /** Mapping qualities around the 255 the transform rewrites. */
    static void buildQualities(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            writer.addAlignment(read(header, "mq-255", 10, "3M2N5M", 0, 255));
            writer.addAlignment(read(header, "mq-254", 20, "3M2N5M", 0, 254));
            writer.addAlignment(read(header, "mq-60", 30, "3M2N5M", 0, 60));
            writer.addAlignment(read(header, "mq-0", 40, "3M2N5M", 0, 0));
            // A read the engine's filters would drop, which ALLOW_ALL_READS keeps: its cigar is
            // longer than its bases.
            writer.addAlignment(read(header, "malformed", 50, "20M", 0, 60));
        }
    }

    /** A pair carrying MC tags, and a secondary alignment. */
    static void buildPairs(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            final SAMRecord first = read(header, "pair", 10, "3M2N5M", 0x1 | 0x40, 60);
            first.setMateReferenceIndex(0);
            first.setMateAlignmentStart(40);
            first.setAttribute("MC", "4M2N4M");
            writer.addAlignment(first);

            final SAMRecord second = read(header, "pair", 40, "4M2N4M", 0x1 | 0x80, 60);
            second.setMateReferenceIndex(0);
            second.setMateAlignmentStart(10);
            second.setAttribute("MC", "3M2N5M");
            writer.addAlignment(second);

            // A secondary alignment with an N, which is only split when asked for.
            writer.addAlignment(read(header, "secondary", 70, "3M2N5M", 0x100, 60));
        }
    }

    /** Two reads sharing a splice, one of which leaves a mismatching overhang across it. */
    static void buildOverhangs(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // The splice this read declares runs from 48 to 57.
            writer.addAlignment(bases(header, "splicer", 45, "3M10N7M",
                    CHR1.substring(44, 47) + CHR1.substring(57, 64), 0, 60));
            // A read whose bases across that splice disagree with the reference: the overhang goes.
            writer.addAlignment(bases(header, "overhanging", 50,
                    "20M", mutate(CHR1.substring(49, 69), 0, 1, 2), 0, 60));
            // And one whose bases agree, which is left alone.
            writer.addAlignment(bases(header, "matching", 50, "20M", CHR1.substring(49, 69), 0, 60));
        }
    }

    /** The N-D-N motif, which is refused unless the refactor argument is given. */
    static void buildNdn(final File file) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            writer.addAlignment(read(header, "ndn", 10, "3M2N1D2N5M", 0, 60));
        }
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    /** A ten-base read, whose bases are the reference's at its start. */
    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final int flags, final int mapq) {
        return bases(header, name, start, cigar, CHR1.substring(start - 1, start + 9), flags, mapq);
    }

    static SAMRecord bases(final SAMFileHeader header, final String name, final int start,
                           final String cigar, final String bases, final int flags, final int mapq) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[bases.length()];
        for (int i = 0; i < quals.length; i++) {
            // A gradient, so a clip is visible in the qualities as well as in the cigar.
            quals[i] = (byte) (20 + i);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(mapq);
        record.setAttribute("RG", "rg1");
        // The three tags the family repair clears.
        record.setAttribute("NM", 1);
        record.setAttribute("MD", "10");
        record.setAttribute("NH", 2);
        return record;
    }

    static String mutate(final String bases, final int... offsets) {
        final char[] chars = bases.toCharArray();
        for (final int offset : offsets) {
            chars[offset] = chars[offset] == 'A' ? 'C' : 'A';
        }
        return new String(chars);
    }

    /** One run of the tool, with its header, its bytes and every record it wrote. */
    static void run(final Path dir, final Path input, final Path fasta, final String label,
                    final String[] extra) throws Exception {
        final Path output = dir.resolve("SplitNCigarReads." + label + ".bam");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-R", fasta.toString(), "-O", output.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new SplitNCigarReads().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            return;
        }

        String commandLine = "";
        final List<String> lines = new ArrayList<>();
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
            for (final SAMRecord record : reader) {
                lines.add(record.getSAMString().trim());
            }
        }
        System.out.printf("commandline\t%s\t%s%n", label, commandLine);
        System.out.printf("output\t%s\t%s%n", label, RecordTransformDump.base64(output));

        final Path index = dir.resolve(output.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("index\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");

        for (final String line : lines) {
            System.out.printf("samline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + CHR1.length() + "\n",
                StandardCharsets.UTF_8);
        return fasta;
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
}
