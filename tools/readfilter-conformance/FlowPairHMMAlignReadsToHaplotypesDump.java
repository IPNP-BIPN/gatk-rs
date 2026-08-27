/*
 * FlowPairHMMAlignReadsToHaplotypes' read-by-haplotype matrix, taken from the reference.
 *
 * Every read scored against every haplotype of a FASTA, written either as the whole matrix or as
 * one line per read naming its best haplotype. The alignment engine's arithmetic is not what this
 * measures: what it measures is the two output formats, and above all what the concise one makes
 * of the reference haplotype.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE EXPANDED FORMAT IS ONE COLUMN PER HAPLOTYPE, headed by the FASTA's own names in the
 *     FASTA's own order, and one row per read;
 *   - THE CONCISE FORMAT IS FIVE COLUMNS: the read, its best haplotype, that score, and the two
 *     differences;
 *   - THE REFERENCE SCORE IS ONLY RECORDED WHILE THE REFERENCE IS THE BEST HAPLOTYPE SO FAR,
 *     the assignment sitting inside the branch that raises the best score, so a reference that
 *     comes after a better haplotype is never recorded at all;
 *   - WHICH MAKES THE COLUMN DEPEND ON THE FASTA'S ORDER: the same haplotypes with the reference
 *     first report a real difference and with the reference last report `Infinity`;
 *   - AND WITH NO --ref-haplotype, OR ONE THE FASTA DOES NOT NAME, every row reports `Infinity`;
 *   - A READ THAT MATCHES NO HAPLOTYPE HAS NO BEST ONE: its name column is EMPTY, its score is
 *     `-Infinity`, and both differences are `NaN`, being an infinity less itself;
 *   - THE TWO ENGINES DISAGREE ABOUT THAT READ, FlowBased leaving it unmatched where
 *     FlowBasedHMM gives it the reference;
 *   - THE BEST HAPLOTYPE IS THE FIRST OF AN EXACT TIE, the comparison being strict, though two
 *     scores that PRINT the same need not be tied: a `Diff_from_second` of 0.000 is a rounded
 *     difference and not an equality;
 *   - EVERY SCORE IS PRINTED WITH THREE DECIMALS, and an infinity prints as `Infinity` rather
 *     than as a number;
 *   - AND AN ENGINE THE TOOL DOES NOT HAVE IS A BARE RuntimeException naming the two it does.
 *
 * Output:
 *
 *     fasta\t<label>=<that fasta, escaped>
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     out\t<label>=<the whole output file, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FlowPairHMMAlignReadsToHaplotypesDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowPairHMMAlignReadsToHaplotypes;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FlowPairHMMAlignReadsToHaplotypesDump {

    static final int CONTIG_LENGTH = 9960;
    static final String FLOW_ORDER = "TGCA";
    static final String UNIT = "TGCA";
    /** Where every haplotype and every read starts. */
    static final int START = 1000;
    static final int LENGTH = 40;

    static byte referenceBase(final int position) {
        return (byte) UNIT.charAt((position - 1) % UNIT.length());
    }

    static byte[] referenceBases(final int start, final int length) {
        final byte[] bases = new byte[length];
        for (int i = 0; i < length; i++) {
            bases[i] = referenceBase(start + i);
        }
        return bases;
    }

    /** The reference bases with the offsets given substituted. */
    static String withChanges(final int... offsets) {
        final byte[] bases = referenceBases(START, LENGTH);
        for (final int offset : offsets) {
            bases[offset] = bases[offset] == 'T' ? (byte) 'A' : (byte) 'T';
        }
        return new String(bases, StandardCharsets.UTF_8);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("flow-pairhmm-align-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FlowPairHMMAlignReadsToHaplotypesDump: every read scored against "
                + "every haplotype, and what the concise format makes of it");

        // The haplotypes: the reference itself and three that differ from it in one, two and
        // three places.
        final String fasta = String.join("\n",
                ">hap_ref", withChanges(),
                ">hap_one", withChanges(10),
                ">hap_two", withChanges(10, 20),
                ">hap_three", withChanges(10, 20, 30),
                "");
        final Path fastaPath = write(dir, "haplotypes.fasta", fasta);
        System.out.printf("fasta\thaplotypes=%s%n", ReferenceQueryDump.escape(fasta));

        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam);
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(bam)));

        // The expanded format, with and without a named reference haplotype: the reference is
        // named in the FASTA either way, so only the concise format can tell the difference.
        run(dir, "expanded", bam, fastaPath, List.of());
        run(dir, "expanded-ref", bam, fastaPath, List.of("--ref-haplotype", "hap_ref"));
        // The concise format, which is where the reference score is read.
        run(dir, "concise-ref", bam, fastaPath,
                List.of("--concise-output-format", "true", "--ref-haplotype", "hap_ref"));
        // The same with no reference haplotype at all.
        run(dir, "concise-no-ref", bam, fastaPath,
                List.of("--concise-output-format", "true"));
        // A reference haplotype that the FASTA does not name, which is the same as none.
        run(dir, "concise-unknown-ref", bam, fastaPath,
                List.of("--concise-output-format", "true", "--ref-haplotype", "nothere"));
        // The other engine, on both formats.
        run(dir, "hmm-expanded", bam, fastaPath, List.of("-E", "FlowBasedHMM"));
        run(dir, "hmm-concise", bam, fastaPath,
                List.of("-E", "FlowBasedHMM", "--concise-output-format", "true",
                        "--ref-haplotype", "hap_ref"));
        // An engine the tool does not have.
        run(dir, "unknown-engine", bam, fastaPath, List.of("-E", "PairHMM"));

        // The SAME haplotypes with the reference LAST, which is what the reference score's
        // recording condition turns on: it is only ever set while the reference is the best
        // haplotype so far, so a reference that comes after a better one is never recorded.
        final String reordered = String.join("\n",
                ">hap_one", withChanges(10),
                ">hap_two", withChanges(10, 20),
                ">hap_three", withChanges(10, 20, 30),
                ">hap_ref", withChanges(),
                "");
        final Path reorderedPath = write(dir, "reordered.fasta", reordered);
        System.out.printf("fasta\treordered=%s%n", ReferenceQueryDump.escape(reordered));
        run(dir, "concise-ref-last", bam, reorderedPath,
                List.of("--concise-output-format", "true", "--ref-haplotype", "hap_ref"));
        run(dir, "expanded-ref-last", bam, reorderedPath, List.of());

        // Two haplotypes equally far from the reference read, so its two best scores tie.
        final String tied = String.join("\n",
                ">hap_ref", withChanges(),
                ">hap_left", withChanges(10),
                ">hap_right", withChanges(20),
                "");
        final Path tiedPath = write(dir, "tied.fasta", tied);
        System.out.printf("fasta\ttied=%s%n", ReferenceQueryDump.escape(tied));
        run(dir, "concise-tied", bam, tiedPath,
                List.of("--concise-output-format", "true", "--ref-haplotype", "hap_ref"));
        run(dir, "expanded-tied", bam, tiedPath, List.of());
    }

    /**
     * The reads. One matching each haplotype exactly, and one matching none.
     *
     * The read that matches the reference exactly is the one whose concise line names the
     * reference as best; the others are the ones whose reference difference is an infinity.
     */
    static void writeReads(final Path bam) {
        final SAMFileHeader header = readHeader();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            writer.addAlignment(record(header, "r-like-ref", withChanges()));
            writer.addAlignment(record(header, "r-like-one", withChanges(10)));
            writer.addAlignment(record(header, "r-like-two", withChanges(10, 20)));
            writer.addAlignment(record(header, "r-like-three", withChanges(10, 20, 30)));
            writer.addAlignment(record(header, "r-like-none", withChanges(5, 15, 25, 35)));
        }
    }

    static SAMRecord record(final SAMFileHeader header, final String name, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(START);
        record.setCigarString(bases.length() + "M");
        record.setReadBases(bases.getBytes(StandardCharsets.UTF_8));
        final byte[] quality = new byte[bases.length()];
        Arrays.fill(quality, (byte) 40);
        record.setBaseQualities(quality);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        // A `tp` of all zeros: every base's quality is the probability of its own hmer length.
        record.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, new byte[bases.length()]);
        return record;
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ULTIMA");
        group.setFlowOrder(FLOW_ORDER);
        group.setSample("sm1");
        header.addReadGroup(group);
        return header;
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
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
        final Path out = dir.resolve("out-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-H", fasta.toString(),
                "-O", out.toString()));
        argv.addAll(extra);
        try {
            new FlowPairHMMAlignReadsToHaplotypes().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("none\t%s=no matrix file%n", label);
            return;
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
