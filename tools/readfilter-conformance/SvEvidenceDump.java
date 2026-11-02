/*
 * FindBreakpointEvidenceSpark and StructuralVariationDiscoveryPipelineSpark through Main, with
 * every debug file they can write: on one partition, and on as many as a small
 * `--bam-partition-size` cuts the same fixtures into.
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
 *     crc\t<label>\t<path under the run's output>\t<length>\t<CRC-32>, for a file over 16 KiB
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
    /** Above this a file's row is its length and CRC: the FASTQs and k-mers of 94 runs would be 100 MB. */
    static final int WHOLE_FILE_LIMIT = 16384;
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
        System.out.println("# SvEvidenceDump: the last two Spark tools through Main, on one partition and on several");
        fixtures();
        for (final String name : new String[] {"sv.fasta", "sv.fasta.fai", "sv.dict", "sv.fasta.img", "reads.bam",
                "reads.bai", "shallow.bam", "shallow.bai", "other.bam", "other.bai", "not-sorted.bam",
                "kill.txt", "kill31.txt", "kill32.txt", "kill-wrong-length.txt", "kill-not-bases.txt", "cnv.vcf", "external.bed", "external-unsorted.bed",
                "external-off.bed", "gaps.txt", "gaps-columns.txt", "gaps-contig.txt", "cross.txt", "cross-off.txt", "external-long.bed",
                "sbi.bam", "sbi.bam.sbi", "wide.bam", "wide.bai", "wrapped.bam", "wrapped.bai"}) {
            final byte[] bytes = Files.readAllBytes(DIR.resolve(name));
            System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(
                    name.endsWith(".img") ? withoutAddresses(bytes) : bytes));
        }
        jdkDeflaterRows();
        findBreakpointRows();
        pipelineRows();
        partitionRows();
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
        new Main().instanceMain(new String[] {"BwaMemIndexImageCreator", "--use-jdk-deflater",
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
        Files.writeString(DIR.resolve("kill31.txt"), reference[0].substring(1000, 1031) + "\n");
        Files.writeString(DIR.resolve("kill32.txt"), reference[0].substring(1000, 1032) + "\n");
        Files.writeString(DIR.resolve("cnv.vcf"), "##fileformat=VCFv4.2\n"
                + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
                + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">\n"
                + "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copies\">\n"
                + "##contig=<ID=chr1,length=16000>\n##contig=<ID=chr2,length=9000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleS\n"
                + "chr1\t6001\tcnvDel\t" + reference[0].charAt(6000) + "\t<DEL>\t60\t.\tEND=6800;SVTYPE=DEL;SVLEN=-799\tGT:CN\t0/1:1\n");
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

        // Features across the whole of the contigs, long ones reaching over partitions.
        Files.writeString(DIR.resolve("external-long.bed"), "chr1\t500\t600\nchr1\t1000\t12000\n"
                + "chr1\t13000\t13100\nchr2\t500\t8000\nchr2\t8500\t8600\n");
        Files.copy(DIR.resolve("reads.bam"), DIR.resolve("sbi.bam"));
        htsjdk.samtools.BAMSBIIndexer.createIndex(DIR.resolve("sbi.bam"), 4096L);
        // The reads under a dictionary that makes KmerCleaner's table over 2^27 entries.
        reheadered("wide.bam", 70_000_000);
        // And one whose table size narrows to a negative int.
        reheadered("wrapped.bam", 1_200_000_000);
    }

    static void reheadered(final String name, final int firstContigLength) throws Exception {
        final List<SAMRecord> records = new ArrayList<>();
        final SAMFileHeader header;
        try (htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault()
                .open(DIR.resolve("reads.bam"))) {
            header = reader.getFileHeader().clone();
            header.getSequence(0).setSequenceLength(firstContigLength);
            for (final SAMRecord record : reader) {
                record.setHeader(header);
                records.add(record);
            }
        }
        write(DIR.resolve(name), header, records);
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
                if (name.endsWith(".vcf")) {
                    bytes = new String(bytes, StandardCharsets.ISO_8859_1)
                            .replaceAll("Date=\"[^\"]*\"", "Date=\"<D>\"")
                            .getBytes(StandardCharsets.ISO_8859_1);
                }
                if (bytes.length > WHOLE_FILE_LIMIT) {
                    final java.util.zip.CRC32 crc = new java.util.zip.CRC32();
                    crc.update(bytes);
                    System.out.printf("crc\t%s\t%s\t%d\t%d%n", label, relative, bytes.length, crc.getValue());
                } else {
                    System.out.printf("file\t%s\t%s\t%s%n", label, relative, Base64.getEncoder().encodeToString(bytes));
                }
            }
        }
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath().toString(), "<cwd>");
    }

    /** The tool with every debug file it has, under a directory of the run's. `rest` may name other k-mers or another image. */
    static String[] debug(final String label, final String input, final String... rest) throws Exception {
        final String fastq = out(label, "fastq");
        Files.createDirectories(Path.of(fastq));
        final List<String> argv = new ArrayList<>(List.of(
                "-I", at(input),
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
        final List<String> named = Arrays.asList(rest);
        if (!named.contains("--kmers-to-ignore")) {
            argv.addAll(List.of("--kmers-to-ignore", at("kill.txt")));
        }
        if (!named.contains("--aligner-index-image")) {
            argv.addAll(List.of("--aligner-index-image", at("sv.fasta.img")));
        }
        argv.addAll(named);
        return argv.toArray(new String[0]);
    }

    static void find(final String label, final String input, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList(debug(label, input, rest)));
        argv.addAll(List.of("-O", out(label, "contigs.sam")));
        run(label, "FindBreakpointEvidenceSpark", argv.toArray(new String[0]));
    }

    /** One k-mer-less refusal: the parser and `validate` come before anything is read. */
    static void refusal(final String label, final String tool, final String... rest) throws Exception {
        run(label, tool, rest);
    }

    /**
     * BAMs through the JDK's deflater, for a host where GKL's level 2 is out of reach. They run
     * first, and the image is made with the flag, because `CommandLineProgram` only ever sets the
     * Intel factory and the JDK's is what a JVM starts with: after any run without the flag, the
     * flag would change nothing here.
     */
    static void jdkDeflaterRows() throws Exception {
        run("find-bam-jdk", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "--use-jdk-deflater", "-O", out("find-bam-jdk", "contigs.bam"));
        pipeline("pipeline-jdk-deflater", "reads.bam", "--use-jdk-deflater");
    }

    static void findBreakpointRows() throws Exception {
        find("find-defaults", "reads.bam");
        find("find-shallow", "shallow.bam");
        find("find-other", "other.bam");
        find("find-queryname", "reads.bam", "--assembled-contigs-output-order", "queryname");
        find("find-gfa", "reads.bam", "--write-gfas");
        find("find-no-mapping", "reads.bam", "--include-mapping-location", "false");
        find("find-pop", "reads.bam", "--pop-variant-bubbles");
        find("find-keep-shadowed", "reads.bam", "--remove-shadowed-contigs", "false");
        find("find-no-expand", "reads.bam", "--expand-assembly-graph", "false");
        find("find-neither", "reads.bam", "--remove-shadowed-contigs", "false", "--expand-assembly-graph", "false");
        find("find-interval-only", "reads.bam", "--interval-only-assembly");
        find("find-adapter", "reads.bam", "--adapter-sequence",
                "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG");
        find("find-k31", "reads.bam", "--k-size", "31", "--kmers-to-ignore", at("kill31.txt"));
        find("find-k-even", "reads.bam", "--k-size", "32", "--kmers-to-ignore", at("kill32.txt"));
        find("find-kill-wrong-length", "reads.bam", "--kmers-to-ignore", at("kill-wrong-length.txt"));
        find("find-kill-not-bases", "reads.bam", "--kmers-to-ignore", at("kill-not-bases.txt"));
        find("find-kill-missing", "reads.bam", "--kmers-to-ignore", at("nowhere.txt"));
        find("find-min-kmers", "reads.bam", "--min-kmers-per-interval", "1000000");
        find("find-cleaner", "reads.bam", "--cleaner-max-intervals", "1", "--cleaner-min-kmer-count", "6",
                "--cleaner-max-copy-number", "2");
        find("find-dust", "reads.bam", "--kmer-max-dust-score", "20");
        find("find-small-fastq", "reads.bam", "--max-fastq-size", "150000");
        find("find-external", "reads.bam", "--external-evidence", at("external.bed"));
        find("find-external-weight", "reads.bam", "--external-evidence", at("external.bed"),
                "--external-evidence-weight", "3", "--external-evidence-uncertainty", "400");
        find("find-external-unsorted", "reads.bam", "--external-evidence", at("external-unsorted.bed"));
        find("find-external-off", "reads.bam", "--external-evidence", at("external-off.bed"));
        find("find-external-missing", "reads.bam", "--external-evidence", at("nowhere.bed"));
        find("find-gaps", "reads.bam", "--exclusion-intervals", at("gaps.txt"), "--exclusion-interval-padding", "200");
        find("find-gaps-columns", "reads.bam", "--exclusion-intervals", at("gaps-columns.txt"));
        find("find-gaps-contig", "reads.bam", "--exclusion-intervals", at("gaps-contig.txt"));
        find("find-gaps-missing", "reads.bam", "--exclusion-intervals", at("nowhere.txt"));
        find("find-cross", "reads.bam", "--cross-contigs-to-ignore", at("cross.txt"));
        find("find-cross-off", "reads.bam", "--cross-contigs-to-ignore", at("cross-off.txt"));
        find("find-evidence-weight", "reads.bam", "--min-evidence-coverage-ratio", "5.0",
                "--min-coherent-evidence-coverage-ratio", "5.0");
        find("find-evidence-thin", "reads.bam", "--min-evidence-coverage-ratio", "0.05",
                "--min-coherent-evidence-coverage-ratio", "0.02");
        find("find-mapq", "reads.bam", "--min-evidence-mapq", "61");
        find("find-match-length", "reads.bam", "--min-evidence-match-length", "140");
        find("find-overhang", "reads.bam", "--allowed-short-fragment-overhang", "0");
        find("find-fragment", "reads.bam", "--max-tracked-fragment-length", "300");
        find("find-depth", "reads.bam", "--high-depth-coverage-peak-factor", "1", "--high-depth-coverage-factor", "1");
        find("find-ratio", "reads.bam", "--assembly-to-mapped-size-ratio-guess", "1");
        find("find-z-drop", "reads.bam", "--z-dropoff", "1");
        find("find-not-sorted", "not-sorted.bam");
        // Without FASTQs: local[1] launches the next task before the failure of the one before has
        // reached the scheduler, so how many more it writes is a race.
        run("find-image-missing", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("nowhere.img"), "-O", out("find-image-missing", "contigs.sam"));
        find("find-input-missing", "nowhere.bam");

        final String[][] invalid = {
                {"k-size", "0"}, {"kmer-max-dust-score", "0"}, {"min-evidence-mapq", "-1"},
                {"min-evidence-match-length", "0"}, {"allowed-short-fragment-overhang", "-1"},
                {"max-tracked-fragment-length", "0"}, {"high-depth-coverage-peak-factor", "0"},
                {"min-evidence-coverage-ratio", "0"}, {"min-coherent-evidence-coverage-ratio", "-1.5"},
                {"min-kmers-per-interval", "0"}, {"cleaner-max-intervals", "0"}, {"cleaner-min-kmer-count", "0"},
                {"cleaner-max-copy-number", "0"}, {"assembly-to-mapped-size-ratio-guess", "0"},
                {"max-fastq-size", "0"}, {"exclusion-interval-padding", "-1"}, {"external-evidence-weight", "0"},
                {"external-evidence-uncertainty", "0"}, {"assembled-contigs-output-order", "unsorted"},
                {"assembled-contigs-output-order", "bogus"}};
        for (final String[] one : invalid) {
            find("find-invalid-" + one[0] + "-" + one[1], "reads.bam", "--" + one[0], one[1]);
        }
        // The one run whose FASTQ directory nothing made first.
        run("find-fastq-no-dir", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "--fastq-dir", out("find-fastq-no-dir", "fastq"),
                "-O", out("find-fastq-no-dir", "contigs.sam"));
        final String label = "find-bad-output";
        run(label, "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "-O", out(label, "contigs.txt"));
        run("find-no-debug-sam", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "-O", out("find-no-debug-sam", "contigs.sam"));
        run("find-no-debug-bam", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "-O", out("find-no-debug-bam", "contigs.bam"));
        run("find-no-debug-bam-queryname", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore",
                at("kill.txt"), "--aligner-index-image", at("sv.fasta.img"),
                "--assembled-contigs-output-order", "queryname", "-O", out("find-no-debug-bam-queryname", "contigs.bam"));
        run("find-no-kmers", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"),
                "--aligner-index-image", at("sv.fasta.img"), "-O", out("find-no-kmers", "contigs.sam"));
        run("find-no-image", "FindBreakpointEvidenceSpark", "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"),
                "-O", out("find-no-image", "contigs.sam"));
    }

    /**
     * Forced partitions: a local file is cut at this size or at 32 MiB, whichever is less, the last
     * split taking a tenth over. reads.bam is some 300 KiB, so these make a handful to a hundred.
     */
    static void partitionRows() throws Exception {
        final String pipelineTool = "StructuralVariationDiscoveryPipelineSpark";
        for (final String size : new String[] {"131072", "65536", "20000", "7000", "1500"}) {
            find("split-" + size, "reads.bam", "--bam-partition-size", size);
        }
        // The slop: a tenth over the size is one split, and more than that is two.
        final long bytes = Files.size(DIR.resolve("reads.bam"));
        find("split-slop-under", "reads.bam", "--bam-partition-size", String.valueOf(bytes * 100 / 105));
        find("split-slop-over", "reads.bam", "--bam-partition-size", String.valueOf(bytes * 100 / 115));
        // Every 64 KiB block its own split or none, and splits with no block starting in them.
        find("split-300", "reads.bam", "--bam-partition-size", "300");
        // `(int) splitSize`: negative and zero leave the default, and the rest is what is left.
        find("split-2g", "reads.bam", "--bam-partition-size", "2147483648");
        find("split-4g", "reads.bam", "--bam-partition-size", "4294967296");
        find("split-4g-and-64k", "reads.bam", "--bam-partition-size", "4295032832");
        find("split-over-block", "reads.bam", "--bam-partition-size", "1073741824");
        find("split-nio", "reads.bam", "--bam-partition-size", "65536", "--use-nio");
        find("split-nio-small", "reads.bam", "--bam-partition-size", "7000", "--use-nio");
        find("split-nio-default", "reads.bam", "--use-nio");
        find("split-sbi", "sbi.bam", "--bam-partition-size", "65536");
        find("split-sbi-small", "sbi.bam", "--bam-partition-size", "7000");
        find("split-intervals", "reads.bam", "--bam-partition-size", "65536", "-L", "chr1:2000-12000");
        find("split-intervals-two", "reads.bam", "--bam-partition-size", "20000", "-L", "chr1:1-9000", "-L", "chr2");
        find("split-shallow", "shallow.bam", "--bam-partition-size", "20000");
        find("split-other", "other.bam", "--bam-partition-size", "65536");
        find("split-k31", "reads.bam", "--bam-partition-size", "65536", "--k-size", "31",
                "--kmers-to-ignore", at("kill31.txt"));
        find("split-cleaner", "reads.bam", "--bam-partition-size", "65536", "--cleaner-max-intervals", "1",
                "--cleaner-min-kmer-count", "6", "--cleaner-max-copy-number", "2");
        find("split-min-kmers", "reads.bam", "--bam-partition-size", "65536", "--min-kmers-per-interval", "1000000");
        find("split-thin", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "0.05",
                "--min-coherent-evidence-coverage-ratio", "0.02");
        find("split-coherent", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "5.0",
                "--min-coherent-evidence-coverage-ratio", "0.1");
        find("split-weight", "reads.bam", "--bam-partition-size", "20000", "--min-evidence-coverage-ratio", "5.0",
                "--min-coherent-evidence-coverage-ratio", "5.0");
        find("split-external", "reads.bam", "--bam-partition-size", "65536", "--external-evidence", at("external.bed"));
        find("split-external-long", "reads.bam", "--bam-partition-size", "20000",
                "--external-evidence", at("external-long.bed"));
        find("split-external-empty-partition", "reads.bam", "--bam-partition-size", "300",
                "--external-evidence", at("external.bed"));
        find("split-gaps", "reads.bam", "--bam-partition-size", "65536", "--exclusion-intervals", at("gaps.txt"),
                "--exclusion-interval-padding", "200");
        find("split-cross", "reads.bam", "--bam-partition-size", "65536", "--cross-contigs-to-ignore", at("cross.txt"));
        find("split-depth", "reads.bam", "--bam-partition-size", "65536", "--high-depth-coverage-peak-factor", "1",
                "--high-depth-coverage-factor", "1");
        find("split-interval-only", "reads.bam", "--bam-partition-size", "65536", "--interval-only-assembly");
        find("split-queryname", "reads.bam", "--bam-partition-size", "65536",
                "--assembled-contigs-output-order", "queryname");
        find("split-not-sorted", "not-sorted.bam", "--bam-partition-size", "7000");
        pipeline("split-pipeline", "reads.bam", "--bam-partition-size", "65536");
        pipeline("split-pipeline-small", "reads.bam", "--bam-partition-size", "7000", "--exp-interpret");
        pipeline("split-pipeline-threshold", "reads.bam", "--bam-partition-size", "20000",
                "--imprecise-variant-evidence-threshold", "1");
        pipeline("split-pipeline-other", "other.bam", "--bam-partition-size", "65536", "--exp-interpret");
        pipeline("split-pipeline-external", "reads.bam", "--bam-partition-size", "65536",
                "--external-evidence", at("external-long.bed"));
        pipeline("split-pipeline-cnv", "reads.bam", "--bam-partition-size", "65536", "--cnv-calls", at("cnv.vcf"));

        // KmerCleaner's table, from twice the dictionary over the partitions, narrowed to an int.
        find("table-wide", "wide.bam");
        find("table-wide-split", "wide.bam", "--bam-partition-size", "65536");
        find("table-wide-k31", "wide.bam", "--k-size", "31", "--kmers-to-ignore", at("kill31.txt"));
        find("table-wrapped", "wrapped.bam");
        find("table-wrapped-split", "wrapped.bam", "--bam-partition-size", "65536");
    }

    static void pipeline(final String label, final String input, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList(debug(label, input, rest)));
        argv.addAll(List.of("-R", at("sv.fasta"), "--contig-sam-file", out(label, "contigs.sam"),
                "-O", out(label, "variants")));
        run(label, "StructuralVariationDiscoveryPipelineSpark", argv.toArray(new String[0]));
    }

    static void pipelineRows() throws Exception {
        final String tool = "StructuralVariationDiscoveryPipelineSpark";
        pipeline("pipeline-defaults", "reads.bam");
        pipeline("pipeline-exp", "reads.bam", "--exp-interpret");
        pipeline("pipeline-other", "other.bam");
        pipeline("pipeline-shallow", "shallow.bam", "--exp-interpret");
        pipeline("pipeline-threshold", "reads.bam", "--imprecise-variant-evidence-threshold", "1");
        pipeline("pipeline-threshold-size", "reads.bam", "--imprecise-variant-evidence-threshold", "1",
                "--max-callable-imprecise-deletion-size", "100");
        pipeline("pipeline-external", "reads.bam", "--external-evidence", at("external.bed"));
        pipeline("pipeline-cross", "reads.bam", "--cross-contigs-to-ignore", at("cross.txt"), "--exp-interpret");
        pipeline("pipeline-both-external", "reads.bam", "--external-evidence", at("external.bed"),
                "--cnv-calls", at("cnv.vcf"));
        pipeline("pipeline-cnv", "reads.bam", "--cnv-calls", at("cnv.vcf"));
        pipeline("pipeline-truth", "reads.bam", "--truth-vcf", at("cnv.vcf"));
        pipeline("pipeline-min-mq", "reads.bam", "--exp-interpret", "--min-mq", "61");
        pipeline("pipeline-min-align", "reads.bam", "--exp-interpret", "--min-align-length", "500");
        pipeline("pipeline-too-few-kmers", "reads.bam", "--min-kmers-per-interval", "1000000");
        pipeline("pipeline-invalid-evidence", "reads.bam", "--k-size", "0");
        pipeline("pipeline-invalid-discover", "reads.bam", "--min-align-length", "0");
        final String label = "pipeline-no-reference";
        run(label, tool, "-I", at("reads.bam"), "--kmers-to-ignore", at("kill.txt"), "--aligner-index-image",
                at("sv.fasta.img"), "--contig-sam-file", out(label, "contigs.sam"), "-O", out(label, "variants"));
        run("pipeline-no-contig-file", tool, "-I", at("reads.bam"), "-R", at("sv.fasta"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "-O", out("pipeline-no-contig-file", "variants"));
        run("pipeline-contigs-bam", tool, "-I", at("reads.bam"), "-R", at("sv.fasta"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "--contig-sam-file", out("pipeline-contigs-bam", "contigs.bam"),
                "-O", out("pipeline-contigs-bam", "variants"));
        run("pipeline-contigs-other-name", tool, "-I", at("reads.bam"), "-R", at("sv.fasta"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "--contig-sam-file", out("pipeline-contigs-other-name", "contigs.txt"),
                "-O", out("pipeline-contigs-other-name", "variants"));
        run("pipeline-nested-out", tool, "-I", at("reads.bam"), "-R", at("sv.fasta"), "--kmers-to-ignore", at("kill.txt"),
                "--aligner-index-image", at("sv.fasta.img"), "--contig-sam-file", out("pipeline-nested-out", "contigs.sam"),
                "-O", out("pipeline-nested-out", "a/b/variants"));
    }
}
