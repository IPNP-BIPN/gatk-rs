/*
 * ClipReads, taken from the reference.
 *
 * The ninth whole tool of the record-transform archetype, and the FIRST THAT WRITES A SECOND
 * OUTPUT THAT IS NOT A BAM. Everything the tool does to a read goes through ReadClipper, which is
 * already ported and already measured; what is left for this dump is the three clippers that build
 * the ops, the representation that decides how the ops are applied, and the statistics file, which
 * is Java text formatting rather than htsjdk bytes.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE QUALITY CLIPPER READS THE READ IN MACHINE-CYCLE ORDER, NOT IN REFERENCE ORDER. It walks
 *     i from the last base to the first but indexes with `isReverseStrand() ? readLen - i - 1 : i`,
 *     so a reverse-strand read is scanned from its front, and the clip it emits is `0..clipPoint`
 *     rather than `clipPoint..readLen-1`. Both strands are in the fixture with the qualities
 *     mirrored, so a port that dropped the strand test would clip the wrong end of one of them;
 *   - THE REPRESENTATION DECIDES WHETHER THE WRITER SORTS. `presorted` is true only for WRITE_NS,
 *     WRITE_NS_Q0S and WRITE_Q0S; the three clipping representations that can move a read get a
 *     sorting writer. The fixture carries a soft-clipped read whose start moves backwards under
 *     REVERT_SOFTCLIPPED_BASES and HARDCLIP_BASES, so the output order is not the input order;
 *   - HARDCLIP_BASES AND REVERT_SOFTCLIPPED_BASES REVERT BEFORE THEY CLIP, in apply, before the
 *     ReadClipper is even constructed. The read every clipper then sees has its soft clip turned
 *     back into a match, so its start has already moved and its cigar has already changed by the
 *     time an op is built. Measured on the fixture's 3S7M read, which arrives at position 6 and
 *     leaves at position 3;
 *   - THE SEQUENCE CLIPPER IS A java.util.regex MATCH, case-insensitive, against the read's bases
 *     as a String, with the pattern reverse-complemented for a reverse-strand read. It loops until
 *     find() fails, so one sequence can clip a read more than once. The counts are keyed by the
 *     sequence string in a TreeMap, so they print in ASCII order of the argument as typed, which
 *     puts an upper-case argument before a lower-case one;
 *   - THE CYCLE CLIPPER IS 1-BASED AND INCLUSIVE, tolerates a stop past the end of the read by
 *     clamping it, drops a range whose start is past the end, and goes through the same strand-
 *     aware flip as the quality clipper;
 *   - THE ADAPTER CLIPPER DOES NOT FLIP FOR STRAND. XF and XT are read as they are, XT is the
 *     first base to clip and XF is the first base NOT to clip, both 1-based, and both zero means
 *     the whole read. It also writes the tf and tm tags onto the read, which is the one place the
 *     tool adds a tag rather than removing bases;
 *   - --read DROPS EVERY READ IT DOES NOT NAME, rather than passing it through, because the whole
 *     of apply including the write is inside the name test. A --read that names nothing therefore
 *     produces an empty BAM and a statistics file whose percentages are `NaN`, which is what
 *     `String.format("%.2f", 0.0 / 0)` prints in Java.
 *
 * The output BAM travels in the golden in full, base64, index included, and the statistics file
 * travels as escaped text. The deflater is pinned and recorded for the reason #160 gives, even
 * though this dump makes no Picard call: the pin costs one line and the audit is what makes the
 * absence of a call checkable rather than assumed.
 *
 * Output:
 *
 *     deflater\t<class>
 *     clipfasta\t<escaped FASTA text>
 *     fixture\t<label>\t<base64 bam>
 *     fixtureindex\t<label>\t<base64 bai>
 *     header\t<label>\t<escaped SAM header>
 *     commandline\t<label>\t<@PG command line>
 *     output\t<label>\t<base64 bam>
 *     index\t<label>\t<base64 bai or absent>
 *     stats\t<label>\t<escaped statistics text>
 *     reads\t<label>\t<name>\t<flags>\t<contig>\t<start>\t<cigar>\t<bases>\t<quals>\t<tags>
 *     error\t<label>\t<class>:<message>
 *
 * Usage: ClipReadsDump
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

import org.broadinstitute.hellbender.tools.ClipReads;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class ClipReadsDump {

    /**
     * The sequences handed to -XF, which reach the tool through
     * ReferenceSequenceFileFactory rather than through the reference machinery.
     *
     * Two of them, so the statistics file has two rows to order, and the second is a
     * reverse-complement palindrome so it matches on either strand.
     */
    static final String CLIP_FASTA =
            ">clipA\n"
            + "GGGGG\n"
            + ">clipB\n"
            + "GGTACC\n";

    public static void main(final String[] args) throws Exception {
        // The factory is static and the first writer wins. This dump calls no Picard entry point,
        // so nothing should replace it; the pin makes that a fact rather than a hope.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Relative on purpose: the string handed to -I and -O is the string recorded inside the
        // output BAM's own @PG, so an absolute temporary path would make every output byte
        // unstable and canonicalization cannot reach inside base64.
        final Path dir = Path.of("clipreads-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ClipReadsDump: ClipReads");
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        final Path clipFasta = dir.resolve("clips.fasta");
        Files.write(clipFasta, CLIP_FASTA.getBytes());
        System.out.printf("clipfasta\t%s%n",
                ReferenceQueryDump.escape(new String(Files.readAllBytes(clipFasta))));

        final Path plain = dir.resolve("plain.bam");
        buildFixture(plain.toFile());
        fixture(dir, plain, "plain");

        // The three clippers on their own, at the default representation.
        clip(dir, plain, "qt", new String[] {"-QT", "10"});
        clip(dir, plain, "cycles", new String[] {"-CT", "1-3,8-12"});
        // Upper case before lower case in the statistics, and the lower-case argument proves the
        // pattern is compiled CASE_INSENSITIVE.
        clip(dir, plain, "seq", new String[] {"-X", "GGGGG", "-X", "acgt"});
        clip(dir, plain, "seqfile", new String[] {"-XF", clipFasta.toString()});
        // All three at once, which is the combination the tool's own documentation shows.
        clip(dir, plain, "combo",
                new String[] {"-QT", "10", "-CT", "1-2", "-X", "GGGGG"});

        // The same clip under each representation. The last three get a sorting writer.
        clip(dir, plain, "q0s", new String[] {"-QT", "10", "-CR", "WRITE_Q0S"});
        clip(dir, plain, "nsq0s", new String[] {"-QT", "10", "-CR", "WRITE_NS_Q0S"});
        clip(dir, plain, "soft", new String[] {"-QT", "10", "-CR", "SOFTCLIP_BASES"});
        clip(dir, plain, "hard", new String[] {"-QT", "10", "-CR", "HARDCLIP_BASES"});
        clip(dir, plain, "revert",
                new String[] {"-QT", "10", "-CR", "REVERT_SOFTCLIPPED_BASES"});

        // The tags, which are the one thing the tool adds rather than removes.
        clip(dir, plain, "adapter", new String[] {"-CA", "true"});
        // Everything except the one named read is dropped, not passed through.
        clip(dir, plain, "onlyread", new String[] {"-QT", "10", "--read", "r0"});
        // Nothing is named, so nothing is examined, so every percentage is NaN.
        clip(dir, plain, "noread", new String[] {"-QT", "10", "--read", "nosuchread"});
        // A length filter that only the short read fails.
        clip(dir, plain, "minlen",
                new String[] {"-QT", "10", "--min-read-length-to-output", "8"});
        // No clipping argument at all: the tool runs, the quality clipper still runs with its
        // default threshold of -1, and nothing is clipped.
        clip(dir, plain, "noclip", new String[] {});

        for (final String label : new String[] {
                "qt", "cycles", "seq", "seqfile", "combo", "q0s", "nsq0s", "soft", "hard",
                "revert", "adapter", "onlyread", "noread", "minlen", "noclip"}) {
            reads(dir, label);
        }
    }

    /** Every read of an output, with everything the three clippers can have changed. */
    static void reads(final Path dir, final String label) throws Exception {
        final Path output = dir.resolve("ClipReads." + label + ".bam");
        try (final SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(output.toFile())) {
            for (final SAMRecord record : reader) {
                final StringBuilder tags = new StringBuilder();
                for (final SAMRecord.SAMTagAndValue tag : record.getAttributes()) {
                    if (tags.length() > 0) {
                        tags.append(';');
                    }
                    tags.append(tag.tag).append('=').append(tag.value);
                }
                System.out.printf("reads\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s%n", label,
                        record.getReadName(), record.getFlags(), record.getReferenceName(),
                        record.getAlignmentStart(), record.getCigarString(),
                        record.getReadString(), record.getBaseQualityString(),
                        tags.length() == 0 ? "-" : tags.toString());
            }
        }
    }

    /**
     * A fixture and the index written beside it.
     *
     * The index travels too because the port's reader needs one to open the file at all, and a
     * test that built its own would be inventing part of the input rather than reading it.
     */
    static void fixture(final Path dir, final Path bam, final String label) throws Exception {
        System.out.printf("fixture\t%s\t%s%n", label, base64(bam));
        final Path index = dir.resolve(bam.getFileName().toString().replace(".bam", ".bai"));
        System.out.printf("fixtureindex\t%s\t%s%n", label, base64(index));
    }

    /**
     * Seven reads, one per thing a clipper can do.
     *
     * The two quality reads carry mirrored qualities on opposite strands, so the machine-cycle
     * scan has to reach opposite ends of them, and a port that ignored the strand would clip
     * nothing at all on the reverse one. The soft-clipped read is what makes the sorting writer
     * observable: reverting its soft clip moves its start from 6 back to 3, in front of a read
     * that was ahead of it in the input.
     */
    static void buildFixture(final File file) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", 60))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        final SAMProgramRecord existing = new SAMProgramRecord("upstream");
        existing.setProgramVersion("1.0");
        header.addProgramRecord(existing);

        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, file)) {
            // Forward, good qualities then bad: the scan clips the tail.
            writer.addAlignment(read(header, "r0", 1, "10M", "ACGTAGGTAC",
                    "IIIII#####", 0));
            // Reverse, with the bad qualities first in ARRAY order, which is last in cycle
            // order: the scan reaches them at the start of its walk and clips the front. A port
            // that walked the array backwards instead would find no clip point at all.
            writer.addAlignment(read(header, "r1", 5, "10M", "GGGGGACGTA",
                    "#####IIIII", 16));
            // Soft-clipped, so REVERT and HARDCLIP move it before anything else happens. Also
            // carries the GGTACC palindrome for -XF.
            writer.addAlignment(read(header, "r2", 6, "3S7M", "TTTGGTACCA",
                    "IIIIIIIIII", 0));
            // Reverse strand with adapter tags, which are NOT flipped for strand.
            writer.addAlignment(read(header, "r3", 15, "10M", "ACGTACGTAC",
                    "IIIIIIIIII", 16, 3, 8));
            // Both adapter tags zero: the whole read is clipped.
            writer.addAlignment(read(header, "r4", 21, "10M", "TTTTTTTTTT",
                    "IIIIIIIIII", 0, 0, 0));
            // Two copies of GGGGG, so the sequence clipper's find() loop has to go round twice.
            writer.addAlignment(read(header, "r5", 25, "10M", "GGGGGGGGGG",
                    "IIIIIIIIII", 0));
            // Five bases, which --min-read-length-to-output 8 rejects and every other run keeps.
            writer.addAlignment(read(header, "r6", 35, "5M", "ACGTA", "IIIII", 0));
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases, final String quals,
                          final int flags) {
        return read(header, name, start, cigar, bases, quals, flags, -1, -1);
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final String bases, final String quals,
                          final int flags, final int xf, final int xt) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[quals.length()];
        for (int i = 0; i < quals.length(); i++) {
            qualities[i] = (byte) (quals.charAt(i) - 33);
        }
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        record.setFlags(flags);
        record.setAttribute("RG", "rg1");
        if (xf >= 0) {
            record.setAttribute(ClipReads.FIVE_PRIME_ADAPTER_LOCATION_TAG, xf);
        }
        if (xt >= 0) {
            record.setAttribute(ClipReads.THREE_PRIME_ADAPTER_LOCATION_TAG, xt);
        }
        return record;
    }

    static void clip(final Path dir, final Path input, final String label, final String[] extra)
            throws Exception {
        final Path output = dir.resolve("ClipReads." + label + ".bam");
        final Path stats = dir.resolve("ClipReads." + label + ".stats");
        // --use-jdk-deflater is the knob that decides which bytes come out, for the same reason
        // PrintReadsDump names it: the GKL deflater's output is not yet reproduced.
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", input.toString(), "-O", output.toString(),
                "-os", stats.toString(),
                "--use-jdk-deflater", "true", "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        try {
            new ClipReads().instanceMain(argv.toArray(new String[0]));
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
        System.out.printf("stats\t%s\t%s%n", label,
                Files.exists(stats)
                        ? ReferenceQueryDump.escape(new String(Files.readAllBytes(stats)))
                        : "absent");
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
