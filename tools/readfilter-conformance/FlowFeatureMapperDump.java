/*
 * FlowFeatureMapper's mapped features, taken from the reference.
 *
 * Every base of a read that differs from the reference, written out as a VCF record of its own.
 * The score is the flow matrix's, which is not what this measures: what it measures is which
 * bases become records at all, what each record carries, and what takes a record away again.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - A FEATURE IS A SINGLE BASE THAT DIFFERS FROM THE REFERENCE AND IS SURROUNDED BY BASES THAT
 *     DO NOT, so a mismatch on the first base of its cigar element has nothing before it and is
 *     never written;
 *   - HOW MANY BASES MUST MATCH IS AN ARGUMENT: two mismatches two apart both survive a surround
 *     of one, because the base between them matches, and both go when the surround is two;
 *   - --snv-identical-bases-after TAKES THE SAME VALUE UNLESS GIVEN ONE OF ITS OWN, so asking
 *     for three before and one after keeps the FIRST of a pair two apart and drops the second;
 *   - THE SURROUND MUST FIT INSIDE THE CIGAR ELEMENT: an element shorter than the surround plus
 *     one is skipped whole, so a mismatch inside a 3M run is seen at a surround of one and not
 *     at a surround of two;
 *   - AN `N` IN THE REFERENCE IS NOT A MISMATCH, though it does count towards the edit distance;
 *   - X_FC1 IS THE READ'S MISMATCH COUNT AND X_FC2 ITS FEATURE COUNT, so the two differ whenever
 *     a mismatch failed the surround test, and X_FC2 moves with the surround arguments;
 *   - X_EDIST IS THE LEVENSHTEIN DISTANCE between the read's aligned bases and the reference the
 *     walker handed it, which is not its mismatch count;
 *   - X_INDEX IS THE OFFSET IN THE WHOLE READ, soft clip included, and X_LENGTH THE UNCLIPPED
 *     LENGTH;
 *   - AN INTERVAL SELECTS READS AND NOT FEATURES, so a read that starts inside it contributes
 *     every feature it carries, including those past the interval's end;
 *   - --min-score AND --max-score BOUND THE SCORE and drop the records outside;
 *   - A DUPLICATE READ IS DROPPED unless --include-dup-reads asks for it, and its records then
 *     carry the flag in X_FLAGS;
 *   - AND --copy-attr COPIES A BAM TAG ONTO EVERY RECORD, under a prefix, with a type and a
 *     description read out of the argument itself.
 *
 * Output:
 *
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     out\t<label>=<the whole output vcf without its `##` lines, escaped>
 *     header\t<label>\t<one `##INFO` or `##FILTER` line>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FlowFeatureMapperDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class FlowFeatureMapperDump {

    static final int CONTIG_LENGTH = 9960;
    static final String FLOW_ORDER = "TGCA";
    /** The reference is `TGCA` repeated, which is also the flow order. */
    static final String UNIT = "TGCA";

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

    /** name, start, cigar, the offsets in the read to change, duplicate, soft-clip. */
    record Read(String name, int start, String cigar, int[] mismatches, boolean duplicate) { }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("flow-feature-mapper-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FlowFeatureMapperDump: which bases of a read become features, and "
                + "what each feature's record carries");

        final Path fasta = writeReference(dir);
        final Path bam = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(bam);
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(bam)));

        run(dir, "default", bam, fasta, List.of());
        // Two identical bases either side rather than one, which drops the features whose
        // neighbours are themselves mismatched.
        run(dir, "identical-two", bam, fasta, List.of("--snv-identical-bases", "2"));
        // A different count after, which is what the second argument is for.
        run(dir, "identical-three-after-one", bam, fasta,
                List.of("--snv-identical-bases", "3", "--snv-identical-bases-after", "1"));
        // The duplicate read, in and out.
        run(dir, "include-duplicates", bam, fasta, List.of("--include-dup-reads", "true"));
        // The score bounds.
        // The scores this fixture produces are 4.95249, 5.98951 and 11.97903, so each bound cuts
        // the set somewhere rather than emptying it.
        run(dir, "max-score-6", bam, fasta, List.of("--max-score", "6.0"));
        run(dir, "min-score-5-5", bam, fasta, List.of("--min-score", "5.5"));
        run(dir, "max-score-0", bam, fasta, List.of("--max-score", "0.0"));
        // A copied attribute, under a prefix and with a type of its own.
        run(dir, "copy-attr", bam, fasta,
                List.of("--copy-attr", "za,Integer,a number", "--copy-attr-prefix", "P_"));
        // An interval, which limits the walk rather than the reads.
        run(dir, "one-interval", bam, fasta, List.of("-L", "chr1:1000-1100"));
        // Every base reported rather than only the mismatched ones, which also turns the
        // surround off.
        run(dir, "report-all-alts", bam, fasta, List.of("--report-all-alts", "true"));
    }

    /**
     * The reads. Each is the reference over its span with a handful of bases changed.
     *
     * One has two mismatches one apart, one has a mismatch two from the end of its match run, one
     * has an `N` in the reference beneath it, one is a duplicate, and one is soft-clipped.
     */
    static void writeReads(final Path bam) {
        final SAMFileHeader header = readHeader();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            final List<Read> reads = List.of(
                    // Three well-separated mismatches.
                    new Read("r-three", 1000, "40M", new int[] {10, 20, 30}, false),
                    // Two mismatches one apart: neither is surrounded, so neither is written.
                    new Read("r-adjacent", 1100, "40M", new int[] {10, 12}, false),
                    // A mismatch on the very first base of the element, which has no base before.
                    new Read("r-edge", 1200, "40M", new int[] {0, 20}, false),
                    // A read whose match runs are short: 3M1I3M1D30M, so the two short elements
                    // are skipped whole.
                    new Read("r-short-elements", 1300, "3M1I3M1D33M", new int[] {1, 20}, false),
                    // Over the reference's `N` run at 1500.
                    new Read("r-over-n", 1480, "40M", new int[] {10, 25}, false),
                    new Read("r-duplicate", 1600, "40M", new int[] {10, 20, 30}, true),
                    new Read("r-clipped", 1700, "5S35M", new int[] {10, 20}, false));
            for (final Read read : reads) {
                writer.addAlignment(record(header, read));
            }
        }
    }

    static SAMRecord record(final SAMFileHeader header, final Read read) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(read.name());
        record.setReferenceName("chr1");
        record.setAlignmentStart(read.start());
        record.setCigarString(read.cigar());
        final byte[] bases = walk(record, read.start());
        for (final int offset : read.mismatches()) {
            bases[offset] = bases[offset] == 'T' ? (byte) 'A' : (byte) 'T';
        }
        record.setReadBases(bases);
        final byte[] quality = new byte[bases.length];
        Arrays.fill(quality, (byte) 40);
        record.setBaseQualities(quality);
        record.setMappingQuality(60);
        record.setDuplicateReadFlag(read.duplicate());
        record.setAttribute("RG", "rg1");
        record.setAttribute("za", read.start());
        // A `tp` of all zeros: every base's quality is the probability of its own hmer length.
        final byte[] tp = new byte[bases.length];
        record.setAttribute(FlowBasedRead.FLOW_MATRIX_TAG_NAME, tp);
        return record;
    }

    /** The reference bases a cigar reads, with `A` for every inserted base. */
    static byte[] walk(final SAMRecord record, final int start) {
        final java.io.ByteArrayOutputStream bases = new java.io.ByteArrayOutputStream();
        int position = start;
        for (final htsjdk.samtools.CigarElement element : record.getCigar()) {
            final htsjdk.samtools.CigarOperator operator = element.getOperator();
            if (operator.consumesReadBases() && operator.consumesReferenceBases()) {
                bases.write(referenceBases(position, element.getLength()), 0, element.getLength());
                position += element.getLength();
            } else if (operator.consumesReadBases()) {
                // A soft clip reads the bases BEFORE the alignment start, so it is written from
                // the reference there rather than as a run of one base.
                for (int i = 0; i < element.getLength(); i++) {
                    bases.write('A');
                }
            } else if (operator.consumesReferenceBases()) {
                position += element.getLength();
            }
        }
        return bases.toByteArray();
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ULTIMA");
        group.setFlowOrder(FLOW_ORDER);
        group.setSample("sample");
        header.addReadGroup(group);
        return header;
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
        final Path out = dir.resolve("out-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "-O", out.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new FlowFeatureMapper().instanceMain(argv.toArray(new String[0]));
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
            System.out.printf("none\t%s=no vcf%n", label);
            return;
        }
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(out).split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            if (line.startsWith("##INFO") || line.startsWith("##FILTER")) {
                // The header lines are their own rows, so a new INFO field is visible without
                // reading the whole header back.
                System.out.printf("header\t%s\t%s%n", label, line);
            } else if (!line.startsWith("##")) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(body.toString(), dir)));
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            final byte[] line = referenceBases(i * 60 + 1, 60);
            // A run of `N` at 1500, which is not a mismatch on either side of the comparison.
            for (int j = 0; j < line.length; j++) {
                final int position = i * 60 + 1 + j;
                if (position >= 1500 && position < 1510) {
                    line[j] = 'N';
                }
            }
            bases.append(new String(line, StandardCharsets.UTF_8)).append("\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
