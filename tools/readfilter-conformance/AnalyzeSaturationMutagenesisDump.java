/*
 * AnalyzeSaturationMutagenesis' three reports, taken from the reference.
 *
 * Reads over a saturation-mutagenesis ORF, counted by what each one turned out to be and by the
 * variants it carried. One prefix produces three files and none is a subset of the others: a
 * census of the reads, one of the codons, and one of the variants.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - PAIRED MODE REFUSES A COORDINATE-SORTED BAM OUTRIGHT, the mates having to be adjacent,
 *     which is what queryname order gives;
 *   - THE TRIMMING READS THE FRAGMENT LENGTH, so a properly-paired read whose TLEN is zero is
 *     trimmed away to nothing and counted LOW_QUALITY however good its bases are;
 *   - EVERY READ LANDS IN EXACTLY ONE REPORT TYPE, and the census is a tree of percentages whose
 *     denominators change from level to level: the top three are over all reads, the three
 *     categories under them over the evaluable ones, and each category's own rows over itself;
 *   - THE OVERLAPPING LINE COUNTS READS AND NOT PAIRS, being twice the pair count, while the
 *     rows beneath it count pairs;
 *   - A READ BELOW --min-mapq IS COUNTED UNMAPPED rather than badly mapped;
 *   - A READ WHOSE HIGH-QUALITY RUN IS SHORTER THAN --min-length IS LOW_QUALITY;
 *   - A VARIANT NEEDS --min-flanking-length WILD-TYPE CALLS ON EACH SIDE, and a read whose
 *     variant sits at its own edge is counted `Insufficient flank` rather than wild type;
 *   - A VARIANT IS ONLY REPORTED ONCE IT HAS --min-variant-obs OBSERVATIONS, so the same reads
 *     under a lower threshold report more rows;
 *   - A VARIANT ROW NAMES ITS BASE CHANGE, ITS CODON CHANGE AND ITS AMINO-ACID CHANGE, the last
 *     translated through --codon-translation;
 *   - THE ORF IS ONE-BASED AND INCLUSIVE and may be several intervals, spliced before
 *     translation; its total length must divide by three and the translation must be exactly
 *     sixty-four characters, both refused by name;
 *   - A REFUSAL REACHES THE CONSOLE AS `A USER ERROR has occurred: <message>`, the exception's
 *     own class never appearing;
 *   - AND THE COUNTERS ARE STATIC FIELDS THAT NOTHING RESETS, so two invocations in one JVM add
 *     up: the same BAM run twice reports thirty-four reads rather than seventeen.
 *
 * Every configuration below runs in a JVM of its own for that last reason.
 *
 * Output:
 *
 *     fixture\t<name>=<that sequence>
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     reads\t<label>=<the whole readCounts file, escaped>
 *     variants\t<label>=<the whole variantCounts file, escaped>
 *     codons\t<label>=<the codonCounts file, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<the tool's own message>
 *
 * Usage: AnalyzeSaturationMutagenesisDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.AnalyzeSaturationMutagenesis;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AnalyzeSaturationMutagenesisDump {

    /** The amplicon: 300 bases, of which 1 to 300 are the ORF, so a hundred codons. */
    static final int LENGTH = 300;
    static final int READ_LENGTH = 120;

    static byte referenceBase(final int position) {
        // A repeating unit of six, so every codon of the ORF differs from its neighbours.
        return (byte) "ACGTCA".charAt((position - 1) % 6);
    }

    static String referenceBases(final int start, final int length) {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < length; i++) {
            bases.append((char) referenceBase(start + i));
        }
        return bases.toString();
    }

    /** name, start, cigar, the offsets to substitute, the offsets to lower in quality, mapq. */
    record Read(String name, int start, String cigar, int[] substitutions, int[] lowQuality,
                int mappingQuality, boolean paired, int mateStart) { }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("analyze-saturation-mutagenesis-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# AnalyzeSaturationMutagenesisDump: reads over an ORF, counted by "
                + "what each one turned out to be");

        final String bases = referenceBases(1, LENGTH);
        System.out.printf("fixture\treference=%s%n", bases);

        final Path fasta = writeReference(dir, bases);
        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam);
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(bam)));

        run(dir, "default", bam, fasta, List.of("--orf", "1-300"));
        // One observation is enough, which is what puts the rarer variants in the report.
        run(dir, "min-obs-one", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1"));
        // A wider flank, which drops the variants near the edges of their reads.
        run(dir, "flank-ten", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1",
                        "--min-flanking-length", "10"));
        // A higher mapping-quality floor, which turns a mapped read into an unmapped one.
        run(dir, "mapq-thirty", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1", "--min-mapq", "30"));
        // A longer minimum length, which turns evaluable reads into low-quality ones.
        run(dir, "min-length-hundred", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1", "--min-length", "100"));
        // Unpaired mode, which evaluates the two mates of every pair on their own.
        run(dir, "unpaired-mode", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1", "--paired-mode", "false"));
        // The disjoint pairs combined rather than evaluated apart.
        run(dir, "combine-disjoint", bam, fasta,
                List.of("--orf", "1-300", "--min-variant-obs", "1",
                        "--dont-ignore-disjoint-pairs", "true"));
        // An ORF of two intervals, spliced before translation.
        run(dir, "two-intervals", bam, fasta,
                List.of("--orf", "1-147,151-300", "--min-variant-obs", "1"));
        // An ORF whose length is not a multiple of three.
        run(dir, "orf-not-a-codon", bam, fasta, List.of("--orf", "1-299"));
        // A translation table of the wrong length.
        run(dir, "short-translation", bam, fasta,
                List.of("--orf", "1-300", "--codon-translation", "KNKN"));

        // The counters are static fields that nothing resets, so running the tool twice in ONE
        // JVM adds the second run's reads to the first's.
        final Path twice = dir.resolve("out-twice-in-one-jvm");
        final String[] argv = {
                "-I", bam.toString(), "-R", fasta.toString(), "-O", twice.toString(),
                "--orf", "1-300", "--min-variant-obs", "1"};
        new AnalyzeSaturationMutagenesis().instanceMain(argv);
        System.out.printf("reads\tonce-in-one-jvm=%s%n", ReferenceQueryDump.escape(
                masked(Files.readString(Path.of(twice + ".readCounts")), dir)));
        new AnalyzeSaturationMutagenesis().instanceMain(argv);
        System.out.printf("reads\ttwice-in-one-jvm=%s%n", ReferenceQueryDump.escape(
                masked(Files.readString(Path.of(twice + ".readCounts")), dir)));
    }

    /**
     * The reads. Nine pairs and three singletons over the amplicon.
     *
     * Three pairs carry the same substitution so it clears the default observation threshold; one
     * carries a rarer one; one has a variant flush against its own edge; one is soft-clipped; one
     * has a run of low-quality bases; one is mapped below the default floor; and two are a
     * disjoint pair whose mates do not overlap.
     */
    static void writeReads(final Path bam) {
        final SAMFileHeader header = readHeader();
        final List<Read> reads = new ArrayList<>();
        // Three pairs with the same substitution at reference position 61, which is codon 21.
        for (int i = 0; i < 3; i++) {
            reads.add(new Read("common-" + i, 1, "120M", new int[] {60}, new int[0], 60, true,
                    100));
            reads.add(new Read("common-" + i, 100, "120M", new int[0], new int[0], 60, true, 1));
        }
        // A rarer substitution, seen once, at position 121.
        reads.add(new Read("rare", 1, "120M", new int[0], new int[0], 60, true, 100));
        reads.add(new Read("rare", 100, "120M", new int[] {21}, new int[0], 60, true, 1));
        // A substitution flush against the read's own start, which the flank test refuses.
        reads.add(new Read("edge", 40, "120M", new int[] {0}, new int[0], 60, true, 100));
        reads.add(new Read("edge", 100, "120M", new int[0], new int[0], 60, true, 40));
        // A run of low-quality bases in the middle, which the trimming cuts around.
        reads.add(new Read("lowq", 1, "120M", new int[0], range(40, 60), 60, true, 100));
        reads.add(new Read("lowq", 100, "120M", new int[0], new int[0], 60, true, 1));
        // Mapped below the default floor of four.
        reads.add(new Read("badmapq", 1, "120M", new int[0], new int[0], 2, false, 0));
        // A disjoint pair: the two mates do not overlap at all.
        reads.add(new Read("disjoint", 1, "60M", new int[0], new int[0], 60, true, 200));
        reads.add(new Read("disjoint", 200, "60M", new int[] {10}, new int[0], 60, true, 1));
        // Two singletons, one soft-clipped.
        reads.add(new Read("single", 90, "120M", new int[0], new int[0], 60, false, 0));
        reads.add(new Read("clipped", 90, "10S110M", new int[] {50}, new int[0], 60, false, 0));

        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            boolean first = true;
            String previous = "";
            for (final Read read : reads) {
                final boolean isFirst = !read.name().equals(previous) || first;
                writer.addAlignment(record(header, read, isFirst));
                previous = read.name();
                first = false;
            }
        }
    }

    static int[] range(final int from, final int to) {
        final int[] offsets = new int[to - from];
        for (int i = 0; i < offsets.length; i++) {
            offsets[i] = from + i;
        }
        return offsets;
    }

    static SAMRecord record(final SAMFileHeader header, final Read read, final boolean isFirst) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(read.name());
        record.setReferenceName("amplicon");
        record.setAlignmentStart(read.start());
        record.setCigarString(read.cigar());
        int length = 0;
        for (final htsjdk.samtools.CigarElement element : record.getCigar()) {
            if (element.getOperator().consumesReadBases()) {
                length += element.getLength();
            }
        }
        final char[] bases = referenceBases(read.start(), length).toCharArray();
        for (final int offset : read.substitutions()) {
            bases[offset] = bases[offset] == 'A' ? 'T' : 'A';
        }
        record.setReadString(new String(bases));
        final byte[] quality = new byte[length];
        Arrays.fill(quality, (byte) 40);
        for (final int offset : read.lowQuality()) {
            quality[offset] = 2;
        }
        record.setBaseQualities(quality);
        record.setMappingQuality(read.mappingQuality());
        record.setAttribute("RG", "rg1");
        if (read.paired()) {
            record.setReadPairedFlag(true);
            record.setProperPairFlag(true);
            record.setFirstOfPairFlag(isFirst);
            record.setSecondOfPairFlag(!isFirst);
            record.setMateReferenceName("amplicon");
            record.setMateAlignmentStart(read.mateStart());
            record.setMateNegativeStrandFlag(!isFirst);
            record.setReadNegativeStrandFlag(!isFirst);
            // The fragment length has to be set: the trimming reads it to keep from running past
            // the end of a fragment shorter than the read, and a TLEN of zero trims everything
            // away, which makes every properly-paired read LOW_QUALITY.
            final int span = Math.abs(read.mateStart() - read.start()) + length;
            record.setInferredInsertSize(isFirst ? span : -span);
        }
        return record;
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("amplicon", LENGTH))));
        // Paired mode refuses a coordinate-sorted BAM outright: the mates have to be adjacent,
        // which is what queryname order gives.
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ILLUMINA");
        group.setSample("sample");
        header.addReadGroup(group);
        return header;
    }

    static Path writeReference(final Path dir, final String bases) throws Exception {
        final Path fasta = dir.resolve("amplicon.fasta");
        final StringBuilder text = new StringBuilder(">amplicon\n");
        for (int i = 0; i < bases.length(); i += 60) {
            text.append(bases, i, Math.min(i + 60, bases.length())).append("\n");
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("amplicon.dict")});
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

    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final List<String> extra) throws Exception {
        final Path prefix = dir.resolve("out-" + label);
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-R", fasta.toString(),
                "-O", prefix.toString()));
        argv.addAll(extra);
        // Each configuration runs in a JVM OF ITS OWN. The tool's counters are static fields
        // that nothing resets, so two invocations in one JVM add up: the run below named
        // `twice-in-one-jvm` is what shows that, and every other run is kept clear of it.
        final List<String> command = new ArrayList<>(List.of(
                "java", "-cp", System.getenv("ORACLE_CP"),
                "org.broadinstitute.hellbender.Main", "AnalyzeSaturationMutagenesis"));
        command.addAll(argv);
        final Process process = new ProcessBuilder(command)
                .redirectErrorStream(true)
                .start();
        final String log = new String(process.getInputStream().readAllBytes(),
                StandardCharsets.UTF_8);
        final int status = process.waitFor();
        if (status != 0) {
            // The tool's own message is the line naming the exception, which is what the golden
            // keeps: the stack trace below it is not stable enough to compare.
            // A refusal reaches the console as `A USER ERROR has occurred: <message>`, the
            // exception's own class never appearing: the wording is all the caller gets.
            String message = "exit " + status;
            for (final String line : log.split("\n")) {
                final String trimmed = line.trim();
                if (trimmed.startsWith("A USER ERROR has occurred:")
                        || trimmed.startsWith("Exception in thread")) {
                    message = trimmed;
                    break;
                }
            }
            System.out.printf("error\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(masked(message, dir)));
            return;
        }
        for (final String[] pair : new String[][] {
                {"reads", ".readCounts"}, {"variants", ".variantCounts"},
                {"codons", ".codonCounts"}}) {
            final Path file = Path.of(prefix + pair[1]);
            if (Files.exists(file)) {
                System.out.printf("%s\t%s=%s%n", pair[0], label,
                        ReferenceQueryDump.escape(masked(Files.readString(file), dir)));
            } else {
                System.out.printf("none\t%s=no %s%n", label, pair[1]);
            }
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
