/*
 * FindBreakpointEvidenceSpark and StructuralVariationDiscoveryPipelineSpark through Main, on one
 * partition, with every debug file they can write.
 *
 * The fixture is a reference and the reads a donor genome would give it: the donor lacks a stretch
 * of chr1 and has a stretch of chr2 inserted backwards further along, so the evidence is real, the
 * reads overlap, and the assembler has something to assemble. Everything is drawn from seeded
 * generators, and nothing printed holds a time, a host address or a temporary path:
 *
 *   - an image holds BWA's pointers, so its row has them zeroed and the runs use the file as made;
 *   - `produced N contigs in S secs.` holds a duration, so its number is masked;
 *   - a directory's listing skips Hadoop's hidden checksums.
 *
 * Output:
 *
 *     input\t<name>\t<base64 of a fixture file>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     stdout\t<label>\t<what the tool printed, escaped>
 *     error\t<label>\t<exception class>:<root cause's message, escaped, cwd masked>
 *     dir\t<label>\t<a directory under the run's output>
 *     file\t<label>\t<path under the run's output>\t<base64>
 *
 * Usage: SvEvidenceDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SvEvidenceDump {

    static final Path DIR = Path.of("sv-evidence-dump");
    static final int READ_LENGTH = 150;
    static final String[] CONTIGS = {"chr1", "chr2"};
    static final int[] LENGTHS = {16000, 9000};

    /** A piece of the donor: where it is on the reference, and which way round. */
    record Piece(int contig, int start, int length, boolean reverse) {}

    /** Where a stretch of a read lands. */
    record Hit(Piece piece, int readStart, int refStart, int length) {}

    static String[] reference;

    public static void main(final String[] args) throws Exception {
        PrintReadsDump.emptyDirectory(DIR);
        Files.createDirectories(DIR);
        System.out.println("# SvEvidenceDump: the last two Spark tools through Main, on one partition");
        fixtures();
        for (final String name : new String[] {"sv.fasta", "sv.fasta.fai", "sv.dict", "sv.fasta.img", "reads.bam",
                "reads.bai", "shallow.bam", "shallow.bai", "other.bam", "other.bai", "not-sorted.bam",
                "kill.txt", "kill-wrong-length.txt", "kill-not-bases.txt", "external.bed", "external-unsorted.bed",
                "external-off.bed", "gaps.txt", "gaps-columns.txt", "gaps-contig.txt", "cross.txt", "cross-off.txt"}) {
            final byte[] bytes = Files.readAllBytes(DIR.resolve(name));
            System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(
                    name.endsWith(".img") ? withoutAddresses(bytes) : bytes));
        }
        findBreakpointRows();
        pipelineRows();
        System.out.flush();
        // Spark's threads would keep the container up.
        System.exit(0);
    }

    static String at(final String name) {
        return DIR.resolve(name).toString();
    }

    static String out(final String label, final String name) throws Exception {
        final Path dir = DIR.resolve("out").resolve(label);
        Files.createDirectories(dir);
        return dir.resolve(name).toString();
    }

    /** The bases of a stretch of the donor, on the reference's forward strand or against it. */
    static byte[] bases(final Piece piece, final int from, final int length) {
        final byte[] copied = new byte[length];
        for (int at = 0; at < length; at++) {
            final int offset = piece.reverse ? piece.start + piece.length - 1 - from - at : piece.start + from + at;
            final char base = reference[piece.contig].charAt(offset);
            copied[at] = (byte) (piece.reverse ? "TGCA".charAt("ACGT".indexOf(base)) : base);
        }
        return copied;
    }

    static byte[] revcomp(final byte[] bases) {
        final byte[] flipped = new byte[bases.length];
        for (int at = 0; at < bases.length; at++) {
            flipped[bases.length - 1 - at] = (byte) "TGCA".charAt("ACGT".indexOf(bases[at]));
        }
        return flipped;
    }

    /** A stretch of the donor, cut at the joins of its pieces. */
    static List<Hit> hits(final List<Piece> donor, final int start, final int length) {
        final List<Hit> found = new ArrayList<>();
        int pieceStart = 0;
        for (final Piece piece : donor) {
            final int from = Math.max(start, pieceStart);
            final int to = Math.min(start + length, pieceStart + piece.length);
            if (from < to) {
                final int inside = from - pieceStart;
                final int refStart = piece.reverse ? piece.start + piece.length - inside - (to - from) : piece.start + inside;
                found.add(new Hit(piece, from - start, refStart, to - from));
            }
            pieceStart += piece.length;
        }
        return found;
    }

    static List<Piece> donor() {
        return List.of(new Piece(0, 0, 6000, false), new Piece(0, 6800, 4200, false),
                new Piece(1, 3000, 1500, true), new Piece(0, 11000, 5000, false),
                new Piece(1, 0, 3000, false), new Piece(1, 4500, 4500, false));
    }

    /** One end of a template, as records of its own reading direction. */
    static List<SAMRecord> placed(final SAMFileHeader header, final String name, final int ordinalFlag,
                                  final byte[] read, final byte[] quals, final List<Hit> pieces,
                                  final boolean readAgainstDonor, final String group, final Random random) {
        // A read kept in the order it was sequenced, whose `Hit`s run along the donor.
        final List<SAMRecord> records = new ArrayList<>();
        Hit primary = pieces.get(0);
        for (final Hit hit : pieces) {
            if (hit.length > primary.length) primary = hit;
        }
        final List<Hit> order = new ArrayList<>();
        order.add(primary);
        for (final Hit hit : pieces) {
            if (hit != primary) order.add(hit);
        }
        final List<String> tags = new ArrayList<>();
        for (final Hit hit : order) {
            final boolean reverse = hit.piece.reverse != readAgainstDonor;
            final int before = hit.readStart;
            final int after = READ_LENGTH - hit.readStart - hit.length;
            // Along the reference, a piece that runs backwards turns the clips round.
            final int first = hit.piece.reverse ? after : before;
            final int last = hit.piece.reverse ? before : after;
            final char clip = hit == primary ? 'S' : 'H';
            tags.add(CONTIGS[hit.piece.contig] + "," + (hit.refStart + 1) + "," + (reverse ? '-' : '+') + ","
                    + cigar(first, hit.length, last, 'S') + ",60,0;");
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(name);
            int flags = 0x1 | ordinalFlag | (reverse ? 0x10 : 0);
            if (hit != primary) flags |= 0x800;
            record.setFlags(flags);
            record.setReferenceName(CONTIGS[hit.piece.contig]);
            record.setAlignmentStart(hit.refStart + 1);
            record.setMappingQuality(60);
            record.setCigarString(cigar(first, hit.length, last, clip));
            final byte[] stored = reverse ? revcomp(read) : read;
            final byte[] storedQuals = quals.clone();
            if (reverse) {
                for (int at = 0; at < quals.length; at++) storedQuals[at] = quals[quals.length - 1 - at];
            }
            if (hit == primary) {
                record.setReadBases(stored);
                record.setBaseQualities(storedQuals);
            } else {
                record.setReadBases(Arrays.copyOfRange(stored, first, first + hit.length));
                record.setBaseQualities(Arrays.copyOfRange(storedQuals, first, first + hit.length));
            }
            record.setAttribute("RG", group);
            record.setAttribute("NM", 0);
            record.setAttribute("AS", hit.length);
            records.add(record);
        }
        if (records.size() > 1) {
            for (int at = 0; at < records.size(); at++) {
                final StringBuilder others = new StringBuilder();
                for (int other = 0; other < tags.size(); other++) {
                    if (other != at) others.append(tags.get(other));
                }
                records.get(at).setAttribute("SA", others.toString());
            }
        }
        return records;
    }

    static String cigar(final int first, final int length, final int last, final char clip) {
        return (first > 0 ? first + "" + clip : "") + length + "M" + (last > 0 ? last + "" + clip : "");
    }

    /** Both ends' mate fields, from what each ended up as. */
    static void mates(final List<SAMRecord> ends, final SAMRecord otherPrimary, final int size) {
        for (final SAMRecord record : ends) {
            record.setFlags(record.getFlags() | (otherPrimary.getReadNegativeStrandFlag() ? 0x20 : 0));
            record.setMateReferenceName(otherPrimary.getReferenceName());
            record.setMateAlignmentStart(otherPrimary.getAlignmentStart());
            record.setAttribute("MC", otherPrimary.getCigarString());
            if (record.getReferenceName().equals(otherPrimary.getReferenceName())) {
                final int span = Math.abs(record.getAlignmentStart() - otherPrimary.getAlignmentStart());
                record.setInferredInsertSize(record.getAlignmentStart() <= otherPrimary.getAlignmentStart() ? span + READ_LENGTH : -span - READ_LENGTH);
            }
            if (size > 0 && Math.abs(record.getInferredInsertSize()) < 1000
                    && record.getReferenceName().equals(otherPrimary.getReferenceName())
                    && record.getReadNegativeStrandFlag() != otherPrimary.getReadNegativeStrandFlag()) {
                record.setFlags(record.getFlags() | 0x2);
            }
        }
    }

    static List<SAMRecord> simulate(final SAMFileHeader header, final List<Piece> donor, final double coverage,
                                    final long seed) {
        final int donorLength = donor.stream().mapToInt(Piece::length).sum();
        final Random random = new Random(seed);
        final int pairs = (int) (coverage * donorLength / (2 * READ_LENGTH));
        final List<SAMRecord> reads = new ArrayList<>();
        for (int pair = 0; pair < pairs; pair++) {
            final int size = 380 + random.nextInt(41);
            final int start = random.nextInt(donorLength - size);
            final String name = "t" + pair;
            final String group = "rg" + (1 + random.nextInt(2));
            final byte[] forward = new byte[READ_LENGTH];
            final byte[] backward = new byte[READ_LENGTH];
            int filled = 0;
            for (final Hit hit : hits(donor, start, READ_LENGTH)) {
                System.arraycopy(bases(hit.piece, hit.piece.reverse
                        ? hit.piece.length - (hit.refStart - hit.piece.start) - hit.length : hit.refStart - hit.piece.start,
                        hit.length), 0, forward, filled, hit.length);
                filled += hit.length;
            }
            filled = 0;
            for (final Hit hit : hits(donor, start + size - READ_LENGTH, READ_LENGTH)) {
                System.arraycopy(bases(hit.piece, hit.piece.reverse
                        ? hit.piece.length - (hit.refStart - hit.piece.start) - hit.length : hit.refStart - hit.piece.start,
                        hit.length), 0, backward, filled, hit.length);
                filled += hit.length;
            }
            final byte[] quals1 = qualities(random);
            final byte[] quals2 = qualities(random);
            noise(forward, random);
            noise(backward, random);
            // The second end is read against the donor: its `Hit`s are numbered from the far end.
            final List<Hit> far = new ArrayList<>();
            for (final Hit hit : hits(donor, start + size - READ_LENGTH, READ_LENGTH)) {
                far.add(new Hit(hit.piece, READ_LENGTH - hit.readStart - hit.length, hit.refStart, hit.length));
            }
            final boolean firstIsFirst = random.nextBoolean();
            final List<SAMRecord> one = placed(header, name, firstIsFirst ? 0x40 : 0x80, forward, quals1,
                    hits(donor, start, READ_LENGTH), false, group, random);
            final List<SAMRecord> two = placed(header, name, firstIsFirst ? 0x80 : 0x40, revcomp(backward), quals2,
                    far, true, group, random);
            mates(one, two.get(0), size);
            mates(two, one.get(0), size);
            reads.addAll(one);
            reads.addAll(two);
        }
        // Some of what a run has to skip or carry: a duplicate, a failed read, unmapped mates.
        final List<SAMRecord> extra = new ArrayList<>();
        for (int at = 0; at < reads.size(); at += 97) {
            final SAMRecord record = reads.get(at);
            if ((record.getFlags() & 0x800) != 0) continue;
            record.setFlags(record.getFlags() | (at % 2 == 0 ? 0x400 : 0x200));
        }
        reads.addAll(extra);
        return reads;
    }

    static byte[] qualities(final Random random) {
        final byte[] quals = new byte[READ_LENGTH];
        for (int at = 0; at < READ_LENGTH; at++) {
            quals[at] = (byte) (25 + random.nextInt(11));
        }
        return quals;
    }

    static void noise(final byte[] bases, final Random random) {
        for (int at = 0; at < bases.length; at++) {
            if (random.nextInt(500) == 0) {
                bases[at] = (byte) "ACGT".charAt(random.nextInt(4));
            }
        }
    }

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        for (int contig = 0; contig < CONTIGS.length; contig++) {
            header.addSequence(new SAMSequenceRecord(CONTIGS[contig], LENGTHS[contig]));
        }
        for (final String[] group : new String[][] {{"rg1", "libA"}, {"rg2", "libA"}}) {
            final SAMReadGroupRecord record = new SAMReadGroupRecord(group[0]);
            record.setSample("sampleS");
            record.setLibrary(group[1]);
            record.setPlatform("ILLUMINA");
            header.addReadGroup(record);
        }
        header.setSortOrder(order);
        return header;
    }

    static void write(final Path path, final SAMFileHeader header, final List<SAMRecord> records) {
        final boolean sorted = header.getSortOrder() == SAMFileHeader.SortOrder.coordinate;
        try (SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(sorted)
                .makeBAMWriter(header, false, path)) {
            records.forEach(writer::addAlignment);
        }
    }

    static void fixtures() throws Exception {
        final Random random = new Random(20261102L);
        reference = new String[CONTIGS.length];
        final StringBuilder fasta = new StringBuilder();
        final StringBuilder dict = new StringBuilder("@HD\tVN:1.6\n");
        for (int contig = 0; contig < CONTIGS.length; contig++) {
            final StringBuilder sequence = new StringBuilder(LENGTHS[contig]);
            for (int at = 0; at < LENGTHS[contig]; at++) {
                sequence.append("ACGT".charAt(random.nextInt(4)));
            }
            reference[contig] = sequence.toString();
            fasta.append('>').append(CONTIGS[contig]).append('\n');
            for (int at = 0; at < LENGTHS[contig]; at += 60) {
                fasta.append(reference[contig], at, Math.min(at + 60, LENGTHS[contig])).append('\n');
            }
            dict.append("@SQ\tSN:").append(CONTIGS[contig]).append("\tLN:").append(LENGTHS[contig]).append('\n');
        }
        final Path fastaPath = DIR.resolve("sv.fasta");
        Files.writeString(fastaPath, fasta.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fastaPath, true);
        Files.writeString(DIR.resolve("sv.dict"), dict.toString(), StandardCharsets.UTF_8);
        new Main().instanceMain(new String[] {"BwaMemIndexImageCreator",
                "-I", fastaPath.toString(), "-O", at("sv.fasta.img")});

        final SAMFileHeader sorted = header(SAMFileHeader.SortOrder.coordinate);
        write(DIR.resolve("reads.bam"), sorted, simulate(sorted, donor(), 40, 7L));
        write(DIR.resolve("shallow.bam"), sorted, simulate(sorted, donor(), 8, 8L));
        // The same genome with no deletion on it, whose evidence is only what chance gives.
        write(DIR.resolve("other.bam"), sorted, simulate(sorted, List.of(new Piece(0, 0, 16000, false),
                new Piece(1, 0, 9000, false)), 30, 9L));
        final SAMFileHeader unsorted = header(SAMFileHeader.SortOrder.unsorted);
        write(DIR.resolve("not-sorted.bam"), unsorted, simulate(unsorted, donor(), 2, 10L));

        final String kmer = reference[0].substring(1000, 1051);
        Files.writeString(DIR.resolve("kill.txt"), kmer + "\n" + reference[1].substring(7000, 7051) + "\n");
        Files.writeString(DIR.resolve("kill-wrong-length.txt"), kmer + "\nACGT\n");
        Files.writeString(DIR.resolve("kill-not-bases.txt"), kmer + "\n" + "N".repeat(51) + "\n");
        Files.writeString(DIR.resolve("external.bed"), "chr1\t2000\t2100\nchr1\t3000\t4000\nchr2\t1000\t1100\n");
        Files.writeString(DIR.resolve("external-unsorted.bed"), "chr1\t3000\t3100\nchr1\t2000\t2100\n");
        Files.writeString(DIR.resolve("external-off.bed"), "chrOther\t3000\t3100\n");
        Files.writeString(DIR.resolve("gaps.txt"), "# a comment\nchr1\t5900\t5950\n");
        Files.writeString(DIR.resolve("gaps-columns.txt"), "chr1\t5900\n");
        Files.writeString(DIR.resolve("gaps-contig.txt"), "chrOther\t5900\t5950\n");
        Files.writeString(DIR.resolve("cross.txt"), "chr2\n");
        Files.writeString(DIR.resolve("cross-off.txt"), "chrOther\n");
    }

    /** BWA's pointers, which nothing reads back; see SparkWalkersDump. */
    static byte[] withoutAddresses(final byte[] image) {
        final byte[] cleaned = image.clone();
        final java.nio.ByteBuffer read = java.nio.ByteBuffer.wrap(image).order(java.nio.ByteOrder.LITTLE_ENDIAN);
        final int bwt = Math.toIntExact(1120 + read.getLong(56) * 4 + read.getLong(1104) * 8);
        final List<Integer> addresses = new ArrayList<>(List.of(1112, bwt + 16, bwt + 32, bwt + 40));
        final int annotations = bwt + 48 + 16 * read.getInt(bwt + 24);
        for (int sequence = 0; sequence < read.getInt(bwt + 8); sequence++) {
            addresses.add(annotations + 40 * sequence + 24);
            addresses.add(annotations + 40 * sequence + 32);
        }
        addresses.forEach(address -> Arrays.fill(cleaned, address, address + 8, (byte) 0));
        return cleaned;
    }

    static void run(final String label, final String tool, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>();
        argv.add(tool);
        argv.addAll(Arrays.asList(rest));
        argv.addAll(Arrays.asList("--spark-master", "local[1]"));
        System.out.printf("argv\t%s\t%s%n", label, ReferenceQueryDump.escape(String.join("\t", argv)));
        final Path outDir = DIR.resolve("out").resolve(label);
        Files.createDirectories(outDir);
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        String failure = null;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            try {
                new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception e) {
                Throwable cause = e;
                while (cause.getCause() != null && cause.getCause() != cause) {
                    cause = cause.getCause();
                }
                failure = e.getClass().getName() + ":" + String.valueOf(cause.getMessage());
            }
        } finally {
            System.out.flush();
            System.setOut(realOut);
        }
        System.out.printf("stdout\t%s\t%s%n", label,
                ReferenceQueryDump.escape(masked(outBytes.toString(StandardCharsets.UTF_8))));
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(masked(failure)));
        }
        final List<Path> produced;
        try (Stream<Path> walk = Files.walk(outDir)) {
            produced = walk.filter(path -> !path.equals(outDir)).sorted().collect(Collectors.toList());
        }
        for (final Path path : produced) {
            final Path relative = outDir.relativize(path);
            final String name = path.getFileName().toString();
            if (Files.isDirectory(path)) {
                System.out.printf("dir\t%s\t%s%n", label, relative);
            } else if (!name.endsWith(".crc")) {
                byte[] bytes = Files.readAllBytes(path);
                if (name.equals("intervals.txt")) {
                    bytes = new String(bytes, StandardCharsets.ISO_8859_1)
                            .replaceAll("contigs in \\d+ secs\\.", "contigs in <S> secs.")
                            .getBytes(StandardCharsets.ISO_8859_1);
                }
                System.out.printf("file\t%s\t%s\t%s%n", label, relative, Base64.getEncoder().encodeToString(bytes));
            }
        }
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath().toString(), "<cwd>");
    }

    /** The tool with every debug file it has, under a directory of the run's. */
    static String[] debug(final String label, final String input, final String... rest) throws Exception {
        final String fastq = out(label, "fastq");
        Files.createDirectories(Path.of(fastq));
        final List<String> argv = new ArrayList<>(List.of(
                "-I", at(input),
                "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"),
                "--read-metadata", out(label, "metadata.txt"),
                "--high-coverage-intervals", out(label, "high-coverage.bed"),
                "--breakpoint-evidence-dir", out(label, "evidence"),
                "--unfiltered-breakpoint-evidence-dir", out(label, "unfiltered"),
                "--breakpoint-intervals", out(label, "intervals.txt"),
                "--qname-intervals-mapped", out(label, "qnames-mapped.txt"),
                "--kmer-intervals", out(label, "kmers.txt"),
                "--qname-intervals-for-assembly", out(label, "qnames-assembly.txt"),
                "--target-link-file", out(label, "links.bedpe"),
                "--fastq-dir", fastq));
        argv.addAll(Arrays.asList(rest));
        return argv.toArray(new String[0]);
    }

    static void find(final String label, final String input, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList(debug(label, input, rest)));
        argv.addAll(List.of("-O", out(label, "contigs.sam")));
        run(label, "FindBreakpointEvidenceSpark", argv.toArray(new String[0]));
    }

    static void findBreakpointRows() throws Exception {
        find("find-defaults", "reads.bam");
    }

    static void pipelineRows() throws Exception {
    }
}
