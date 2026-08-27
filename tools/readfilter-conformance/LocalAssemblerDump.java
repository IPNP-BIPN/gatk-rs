/*
 * LocalAssembler's GFA and its scaffolds, taken from the reference.
 *
 * A de Bruijn assembly of the reads over one interval, written as a graph and as a FASTA of the
 * paths through it. What is measured is the graph the reads produce and what each argument
 * changes about it.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE GRAPH IS GFA 2.0, not 1.0: its header reads `VN:Z:2.0`, its segments are `S` lines
 *     carrying a length and three offsets, and its edges are `E` lines and not `L` ones;
 *   - THE KMER SIZE IS 31 AND FIXED, so a read of twenty bases contributes nothing and assembles
 *     into an empty graph;
 *   - A CONTIG MUST BE SEEN --min-thin-observations TIMES, which is four, so a fixture of one
 *     read per sequence assembles into nothing at all: every case here writes five copies;
 *   - RAISING THAT FLOOR TO TEN CHANGES NOTHING when the reads already meet it, so the argument
 *     is a floor on OBSERVATIONS and not on reads;
 *   - EVERY BASE UNDER --q-min BREAKS THE KMER RUN, so a floor above the fixture's own quality
 *     empties the graph outright;
 *   - AN `N` BREAKS THE RUN TOO, BUT THE GAP-FILLING PUTS IT BACK: a read with one `N` in the
 *     middle still assembles into a single contig of its full length, the two halves being
 *     rejoined by --min-gapfill-count observations of the missing kmers;
 *   - TWO READS THAT OVERLAP BY MORE THAN A KMER ASSEMBLE INTO ONE CONTIG of their union, and
 *     two that do not stay as two contigs with no edge between them;
 *   - A SUBSTITUTION IN THE MIDDLE OF AN OTHERWISE SHARED SEQUENCE MAKES A BUBBLE: four segments
 *     and four edges, the trunk either side and one branch per allele;
 *   - --no-scaffolding LEAVES THE GRAPH ALONE AND CHANGES THE FASTA, so the two runs agree on
 *     the GFA byte for byte and differ in the paths;
 *   - THE FASTA NAMES EACH PATH `<assembly>_t<n>` AND LISTS ITS SEGMENTS, with `RC` on a segment
 *     traversed backwards;
 *   - AND AN INTERVAL WITH NO READ IN IT WRITES BOTH FILES ALL THE SAME, the GFA holding its
 *     header alone and the FASTA nothing.
 *
 * The reference is drawn from a linear congruential generator whose state carries forward. Two
 * earlier fixtures were periodic and the assembly collapsed into a single forty-base contig
 * looping on itself: one repeated `ATGCAGCATG`, and one hashed the POSITION, which still left
 * runs like `ATTACC` every six bases.
 *
 * Output:
 *
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     gfa\t<label>=<the whole gfa, escaped>
 *     fasta\t<label>=<the whole fasta, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: LocalAssemblerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.LocalAssembler;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LocalAssemblerDump {

    static final int CONTIG_LENGTH = 6000;

    /**
     * The whole reference, drawn once from a linear congruential generator.
     *
     * The sequence has to be genuinely aperiodic or every 31-mer repeats and the assembly
     * collapses. Two earlier fixtures failed that way: one used `ATGCAGCATG` repeated, and one
     * used a hash OF THE POSITION, which still produced runs like `ATTACC` every six bases. Only
     * a generator whose state carries forward gives a sequence the assembler treats as unique.
     */
    static final String REFERENCE = drawReference();

    static String drawReference() {
        final StringBuilder bases = new StringBuilder(CONTIG_LENGTH);
        long state = 20260828L;
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            state = state * 6364136223846793005L + 1442695040888963407L;
            bases.append("ACGT".charAt((int) ((state >>> 33) & 3L)));
        }
        return bases.toString();
    }

    /** The reference from `start`, one-based, for `length` bases. */
    static String referenceBases(final int start, final int length) {
        return REFERENCE.substring(start - 1, start - 1 + length);
    }

    /** name, start, sequence, qualities. */
    record Read(String name, int start, String bases, byte quality) { }

    /**
     * The same read `COPIES` times over.
     *
     * A contig has to be seen --min-thin-observations times to survive, which is four by default,
     * so a fixture of one read per sequence assembles into nothing at all.
     */
    static final int COPIES = 5;

    static List<Read> copies(final String name, final int start, final String bases,
                             final byte quality) {
        final List<Read> reads = new ArrayList<>();
        for (int i = 0; i < COPIES; i++) {
            reads.add(new Read(name + "-" + i, start, bases, quality));
        }
        return reads;
    }

    static List<Read> concat(final List<Read>... groups) {
        final List<Read> reads = new ArrayList<>();
        for (final List<Read> group : groups) {
            reads.addAll(group);
        }
        return reads;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("local-assembler-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# LocalAssemblerDump: the graph the reads over one interval produce");

        final Path fasta = writeReference(dir);

        // One read, which assembles into one contig.
        run(dir, "one-read", fasta, copies("r1", 1000, referenceBases(1000, 100), (byte) 40));

        // Two reads overlapping by more than a kmer, which join.
        run(dir, "overlapping-reads", fasta, concat(
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40),
                copies("r2", 1050, referenceBases(1050, 100), (byte) 40)));

        // Two reads that do not overlap at all, which stay apart.
        run(dir, "disjoint-reads", fasta, concat(
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40),
                copies("r2", 1500, referenceBases(1500, 100), (byte) 40)));

        // A read too short to kmerize at all.
        run(dir, "read-shorter-than-k", fasta,
                copies("r1", 1000, referenceBases(1000, 20), (byte) 40));

        // A substitution in the middle of an otherwise shared sequence: a bubble.
        final String shared = referenceBases(1000, 120);
        final char[] mutated = shared.toCharArray();
        mutated[60] = mutated[60] == 'A' ? 'T' : 'A';
        run(dir, "bubble", fasta, concat(
                copies("r1", 1000, shared, (byte) 40),
                copies("r2", 1000, new String(mutated), (byte) 40)));

        // A low-quality base in the middle, which breaks the kmer run.
        run(dir, "low-quality-base", fasta,
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40), "--q-min", "45");

        // An `N` in the middle, which breaks it too.
        final char[] withN = referenceBases(1000, 100).toCharArray();
        withN[50] = 'N';
        run(dir, "n-in-the-middle", fasta,
                copies("r1", 1000, new String(withN), (byte) 40));

        // The thin-observation floor raised past what the reads support.
        run(dir, "thin-observations-ten", fasta, concat(
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40),
                copies("r2", 1050, referenceBases(1050, 100), (byte) 40)),
                "--min-thin-observations", "10");

        // The traversals instead of the scaffolds.
        run(dir, "no-scaffolding", fasta, concat(
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40),
                copies("r2", 1050, referenceBases(1050, 100), (byte) 40)),
                "--no-scaffolding", "true");

        // An interval that holds no read at all.
        run(dir, "empty-interval", fasta,
                copies("r1", 1000, referenceBases(1000, 100), (byte) 40),
                "-L", "chr1:4000-4100");
    }

    static void run(final Path dir, final Path fasta, final String label, final List<Read> reads,
                    final String... extra) throws Exception {
        run(dir, label, fasta, reads, extra);
    }

    static void run(final Path dir, final String label, final Path fasta, final List<Read> reads,
                    final String... extra) throws Exception {
        final Path bam = dir.resolve(label + ".bam");
        writeReads(bam, reads);
        System.out.printf("sam\t%s=%s%n", label, ReferenceQueryDump.escape(asSam(bam)));
        final Path gfa = dir.resolve("out-" + label + ".gfa");
        final Path out = dir.resolve("out-" + label + ".fa");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-R", fasta.toString(),
                "--assembly-name", label,
                "--gfa-file", gfa.toString(),
                "--fasta-file", out.toString()));
        final List<String> extras = Arrays.asList(extra);
        if (!extras.contains("-L")) {
            argv.addAll(List.of("-L", "chr1:900-1700"));
        }
        argv.addAll(extras);
        try {
            new LocalAssembler().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        for (final Path file : List.of(gfa, out)) {
            final String kind = file.toString().endsWith(".gfa") ? "gfa" : "fasta";
            if (Files.exists(file)) {
                System.out.printf("%s\t%s=%s%n", kind, label,
                        ReferenceQueryDump.escape(masked(Files.readString(file), dir)));
            } else {
                System.out.printf("none\t%s=no %s%n", label, kind);
            }
        }
    }

    static void writeReads(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = readHeader();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read read : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(read.name());
                record.setReferenceName("chr1");
                record.setAlignmentStart(read.start());
                record.setCigarString(read.bases().length() + "M");
                record.setReadString(read.bases());
                final byte[] quality = new byte[read.bases().length()];
                Arrays.fill(quality, read.quality());
                record.setBaseQualities(quality);
                record.setMappingQuality(60);
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        return header;
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final String bases = REFERENCE;
        final StringBuilder text = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            text.append(bases, i, Math.min(i + 60, bases.length())).append("\n");
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static String asSam(final Path bam) throws Exception {
        final StringBuilder text = new StringBuilder();
        try (final htsjdk.samtools.SamReader reader =
                     htsjdk.samtools.SamReaderFactory.makeDefault().open(bam.toFile())) {
            for (final SAMRecord record : reader) {
                text.append(record.getSAMString());
            }
        }
        return text.toString();
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
