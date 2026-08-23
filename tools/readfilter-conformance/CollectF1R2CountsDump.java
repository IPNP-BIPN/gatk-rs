/*
 * CollectF1R2Counts's output, taken from the reference.
 *
 * The counts Mutect2's orientation-bias filter is trained on: for every three-base reference
 * context, how deep the reference sites were, and what the alt sites looked like. The traversal is
 * a LocusWalker; what is measured here is F1R2CountsCollector, which decides which pileups count
 * at all and which of three places a site is written to.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - A SITE GOES TO ONE OF THREE PLACES: a reference site increments a depth histogram, an alt
 *     site with exactly one alt read increments a depth-one histogram keyed by context, alt base
 *     and orientation, and any other alt site becomes a row of the alt table;
 *   - THE ALT BASE IS maxElementIndex OVER THE BASE COUNTS WITH THE REFERENCE SET TO -1, so a tie
 *     between two alt bases goes to the lower base index and a site with no alt read at all reads
 *     as a reference site;
 *   - AND THE ALT TABLE'S DEPTH COLUMN IS refCount + altCount, NOT THE PILEUP'S DEPTH, so that
 *     tied site reports four over a pileup of six, the second alt base counting nowhere;
 *   - ORIENTATION IS isReverseStrand() != isFirstOfPair(), which for an UNPAIRED read makes a
 *     forward read F2R1 and a reverse read F1R2;
 *   - THE BASE QUALITY TEST IS STRICT: a base at exactly --f1r2-min-bq is excluded from the
 *     pileup, so a lone alt base at that quality leaves a reference site behind, and one below it
 *     turns the same site into a depth-one alt;
 *   - THE SITE IS SKIPPED WHOLE when its median mapping quality is below --f1r2-median-mq, when
 *     more than a hundredth of the pileup is indel, or when the three-base context runs off the
 *     reference or contains an N;
 *   - AND WHEN IT IS SKIPPED FOR ONE SAMPLE IT IS SKIPPED FOR EVERY SAMPLE AFTER IT: the loop over
 *     the split pileup RETURNS rather than continuing, and so do the reference-site and
 *     one-alt-read branches, so at any site only the samples up to the first one that returns are
 *     ever counted, in the order the split map hands them over. That order is a HashMap order over
 *     the sample names and not the header's: over alpha and bravo it is bravo first, and every one
 *     of alpha's alt sites is lost because bravo is reference at each of them;
 *   - DEPTH IS CAPPED AT --f1r2-max-depth before it is counted, and the histograms are prefilled
 *     with every bin from one to that depth;
 *   - EVERY CONTEXT IS PRESENT WHETHER OR NOT IT WAS SEEN, all 64 of them, so the shape of the
 *     output does not depend on the data;
 *   - THE OUTPUT IS A TAR.GZ of one alt table, one reference histogram and one alt histogram per
 *     sample, named by the URL-encoded sample name;
 *   - AND THE ORDER OF THE REFERENCE HISTOGRAMS IS A java.util.HashMap ORDER OVER THE 64 CONTEXT
 *     STRINGS, which is reproducible, while the alt histograms are keyed by a pair holding an ENUM,
 *     whose hashCode is an identity hash, so their order is not. The alt histogram file is
 *     therefore reported here with its columns sorted, and the reference histogram as it is.
 *
 * Output:
 *
 *     tar\t<label>=<entry names, comma separated>
 *     file\t<label>\t<name>=<the whole file, escaped>
 *     sorted\t<label>\t<name>=<the whole file with its histogram columns sorted, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CollectF1R2CountsDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.broadinstitute.hellbender.tools.walkers.readorientation.CollectF1R2Counts;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class CollectF1R2CountsDump {

    static final int LENGTH = 200;
    static final int READ_LENGTH = 30;

    /**
     * The reference: a forty-base motif repeated, which gives many different three-base contexts
     * where a four-base cycle would give four, with one N to make a context that is skipped.
     */
    static final String MOTIF = "ACGTTGCAAGGCTTACCATGGACTTCAGATCCGTAACGGT";

    static String reference() {
        final StringBuilder bases = new StringBuilder();
        for (int index = 0; index < LENGTH; index++) {
            bases.append(index == 99 ? 'N' : MOTIF.charAt(index % MOTIF.length()));
        }
        return bases.toString();
    }

    static char refBase(final int oneBased) {
        return reference().charAt(oneBased - 1);
    }

    /** The rank-th base of ACGT that is not the reference base, counted from zero. */
    static char altBase(final int oneBased, final int rank) {
        int seen = 0;
        for (final char base : new char[] {'A', 'C', 'G', 'T'}) {
            if (base != refBase(oneBased) && seen++ == rank) {
                return base;
            }
        }
        throw new IllegalStateException("unreachable");
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("collect-f1r2-counts-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CollectF1R2CountsDump: the F1R2 counts the orientation-bias filter is trained on");

        final Path fasta = writeReference(dir);

        // One sample, where every branch of the collector is reachable.
        final Path single = dir.resolve("single.bam");
        writeBam(single, List.of("alpha"), true);
        run(dir, "single-sample", single, fasta, 4);

        // The same reads under two samples, where the returns stop the loop.
        final Path pair = dir.resolve("pair.bam");
        writeBam(pair, List.of("alpha", "bravo"), true);
        run(dir, "two-samples", pair, fasta, 4);

        // A depth cap low enough to be reached, and one high enough not to be.
        run(dir, "max-depth-two", single, fasta, 2);
        // A median mapping quality low enough to let the low-quality block through.
        run(dir, "low-median-mq", single, fasta, 4, "--f1r2-median-mq", "20");
        // A base quality one below the default, which lets the excluded base back in.
        run(dir, "min-bq-nineteen", single, fasta, 4, "--f1r2-min-bq", "19");
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("f1r2.fasta");
        final StringBuilder text = new StringBuilder(">chr1\n");
        final String bases = reference();
        for (int index = 0; index < bases.length(); index += 50) {
            text.append(bases, index, Math.min(index + 50, bases.length())).append('\n');
        }
        Files.writeString(fasta, text.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("f1r2.dict")});
        System.out.printf("reference\tchr1=%s%n", bases);
        return fasta;
    }

    /**
     * Four blocks of reads, each thirty bases long and every one unpaired.
     *
     * The first block is where the counting happens: three forward and three reverse reads per
     * sample, with alt bases planted at five offsets. The second spans the N. The third carries a
     * mapping quality the default median test refuses. The fourth carries a deletion, and its
     * cigar consumes exactly the thirty bases the read holds: a cigar that does not is dropped
     * whole by WellformedReadFilter's read-length rule long before the collector sees it, and the
     * block would then measure nothing at all.
     */
    static void writeBam(final Path bam, final List<String> samples, final boolean index)
            throws Exception {
        final SAMFileHeader header = new SAMFileHeader(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String sample : samples) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg-" + sample);
            group.setSample(sample);
            header.addReadGroup(group);
        }
        final List<SAMRecord> records = new ArrayList<>();
        for (final String sample : samples) {
            records.addAll(block(header, sample, 41, 60, false));
            records.addAll(block(header, sample, 91, 60, false));
            records.addAll(block(header, sample, 131, 30, false));
            records.addAll(block(header, sample, 161, 60, true));
        }
        records.sort((left, right) -> Integer.compare(left.getAlignmentStart(), right.getAlignmentStart()));
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(index)
                .makeBAMWriter(header, true, bam)) {
            records.forEach(writer::addAlignment);
        }
    }

    /** Six reads over one window, three forward and three reverse. */
    static List<SAMRecord> block(final SAMFileHeader header, final String sample, final int start,
                                 final int mappingQuality, final boolean deletion) {
        final List<SAMRecord> records = new ArrayList<>();
        for (int index = 0; index < 6; index++) {
            final boolean reverse = index >= 3;
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(sample + "-" + start + "-" + index);
            record.setReferenceName("chr1");
            record.setAlignmentStart(start);
            record.setMappingQuality(mappingQuality);
            record.setReadNegativeStrandFlag(reverse);
            record.setAttribute("RG", "rg-" + sample);
            final StringBuilder bases = new StringBuilder();
            final byte[] quals = new byte[READ_LENGTH];
            Arrays.fill(quals, (byte) 40);
            for (int offset = 0; offset < READ_LENGTH; offset++) {
                final int locus = start + offset + (deletion && offset >= 12 ? 2 : 0);
                char base = refBase(locus);
                // The first block carries the plants; every other block is pure reference.
                if (start == 41) {
                    // Two alt reads of the first sample, which is an alt table row.
                    if (offset == 3 && sample.equals("alpha") && index < 2) {
                        base = altBase(locus, 0);
                    }
                    // One alt read, on the reverse strand, which is a depth-one histogram entry.
                    if (offset == 6 && sample.equals("alpha") && index == 3) {
                        base = altBase(locus, 0);
                    }
                    // Two alt reads of the SECOND sample, at a locus where the first sample is
                    // reference and returns before this one is reached.
                    if (offset == 9 && sample.equals("bravo") && index < 2) {
                        base = altBase(locus, 0);
                    }
                    // One alt read whose base quality is exactly the default threshold.
                    if (offset == 14 && sample.equals("alpha") && index == 0) {
                        base = altBase(locus, 0);
                        quals[offset] = 20;
                    }
                    // Two alt reads of one base and two of another, which is a tie that
                    // maxElementIndex settles on the lower base index.
                    if (offset == 19 && sample.equals("alpha")) {
                        if (index < 2) {
                            base = altBase(locus, 0);
                        } else if (index < 4) {
                            base = altBase(locus, 1);
                        }
                    }
                }
                bases.append(base);
            }
            record.setReadString(bases.toString());
            record.setBaseQualities(quals);
            record.setCigarString(deletion ? "12M2D18M" : READ_LENGTH + "M");
            records.add(record);
        }
        return records;
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta,
                    final int maxDepth, final String... extra) throws Exception {
        final Path work = dir.resolve(label);
        Files.createDirectories(work);
        final Path out = work.resolve("counts.tar.gz");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", bam.toString(), "-R", fasta.toString(), "-O", out.toString(),
                "--f1r2-max-depth", Integer.toString(maxDepth)));
        argv.addAll(Arrays.asList(extra));
        try {
            new CollectF1R2Counts().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        final Map<String, String> entries = untar(out);
        System.out.printf("tar\t%s=%s%n", label, String.join(",", entries.keySet()));
        for (final Map.Entry<String, String> entry : entries.entrySet()) {
            final String name = entry.getKey();
            if (name.endsWith(".alt_histogram")) {
                // The column order of this one is an identity hash order, so it is not stable
                // from run to run and only the sorted form can be a golden.
                System.out.printf("sorted\t%s\t%s=%s%n", label, name,
                        ReferenceQueryDump.escape(sortColumns(entry.getValue())));
            } else {
                System.out.printf("file\t%s\t%s=%s%n", label, name,
                        ReferenceQueryDump.escape(masked(entry.getValue(), dir)));
            }
        }
    }

    /**
     * Every entry of the tar.gz, keyed by name and reported in name order.
     *
     * The order inside the archive is File.listFiles over a temporary directory, which is the
     * filesystem's and not the tool's, so it is not what this compares.
     */
    static Map<String, String> untar(final Path archive) throws Exception {
        final Map<String, String> entries = new TreeMap<>();
        try (final InputStream raw = Files.newInputStream(archive);
             final TarArchiveInputStream tar =
                     new TarArchiveInputStream(new GzipCompressorInputStream(raw))) {
            TarArchiveEntry entry;
            while ((entry = tar.getNextEntry()) != null) {
                if (entry.isDirectory()) {
                    continue;
                }
                final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
                tar.transferTo(bytes);
                entries.put(entry.getName(), bytes.toString(StandardCharsets.UTF_8));
            }
        }
        return entries;
    }

    /**
     * A metrics file whose histogram table has its columns put in order.
     *
     * The bin column stays first; the rest are sorted by name, and every row is reordered with the
     * header. Everything outside the table is left alone.
     */
    static String sortColumns(final String metrics) {
        final List<String> lines = new ArrayList<>(Arrays.asList(metrics.split("\n", -1)));
        int header = -1;
        for (int index = 0; index < lines.size(); index++) {
            if (lines.get(index).startsWith("## HISTOGRAM")) {
                header = index + 1;
                break;
            }
        }
        if (header < 0 || header >= lines.size()) {
            return metrics;
        }
        final String[] columns = lines.get(header).split("\t", -1);
        final Integer[] order = new Integer[columns.length - 1];
        for (int index = 0; index < order.length; index++) {
            order[index] = index + 1;
        }
        Arrays.sort(order, (left, right) -> columns[left].compareTo(columns[right]));
        for (int index = header; index < lines.size(); index++) {
            final String line = lines.get(index);
            if (line.isEmpty()) {
                break;
            }
            final String[] fields = line.split("\t", -1);
            final StringBuilder rebuilt = new StringBuilder(fields[0]);
            for (final int column : order) {
                rebuilt.append('\t').append(column < fields.length ? fields[column] : "");
            }
            lines.set(index, rebuilt.toString());
        }
        return String.join("\n", lines);
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
