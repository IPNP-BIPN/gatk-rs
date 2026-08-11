/*
 * What a ReadWalker hands to apply(), taken from the reference.
 *
 * Every read-based tool in GATK is this traversal: the reads data source is iterated, each read
 * goes through the pre-filter transformer, the tool's read filters, then the post-filter
 * transformer, and what survives reaches apply() together with a ReferenceContext and a
 * FeatureContext built from *that read's* interval. The order of those three steps, and which
 * reads never arrive, is the tool's output.
 *
 * Two of its decisions are not visible from any tool's own output:
 *
 *   - getReadInterval returns **null** both for an unmapped read and for a mapped read whose
 *     coordinates do not form a valid SimpleInterval. A mapped read with an empty cigar has an
 *     alignment end one before its start, so it reaches apply() with an *empty* ReferenceContext:
 *     the walker is handed a read and no reference under it;
 *   - the filter runs between the two transformers, so the filter judges the pre-transformed read
 *     while the walker receives the post-transformed one.
 *
 * The probe walker below is a real ReadWalker run through the real command line, so the traversal
 * measured is the one a tool gets rather than a reconstruction of it. The fixture BAM, its index
 * and the FASTA all travel in the golden.
 *
 * Output:
 *
 *     bam\t<base64>            fai\t<escaped>
 *     bai\t<base64>            fasta\t<escaped>
 *     apply\t<label>\t<index>\t<name>|<start>|<cigar>|<flags>|<window>|<bases>
 *     summary\t<label>\t<CountingReadFilter summary line>
 *     count\t<label>\t<number of apply calls>
 *
 * Usage: ReadWalkerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class ReadWalkerDump {

    static final int CONTIG_LENGTH = 200;

    /**
     * A 200-base contig per line-wrapped chunk of 60, with a soft-masked stretch so the
     * upper-casing every ReferenceContext applies stays visible in the window.
     */
    static final String FASTA = ">chr1\n"
            + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            + "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt\n"
            + "ACGTNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            + "ACGTACGTACGTACGTACGT\n"
            + ">chr2\n"
            + "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC\n"
            + "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC\n"
            + "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC\n"
            + "TTTTGGGGCCCCAAAATTTT\n";

    /** Every apply() call of the current traversal, filled by the probe walker. */
    static final List<String> APPLIED = new ArrayList<>();

    @CommandLineProgramProperties(
            summary = "Records what a ReadWalker hands to apply()",
            oneLineSummary = "ReadWalker traversal probe",
            programGroup = ReadDataManipulationProgramGroup.class)
    public static final class ProbeWalker extends ReadWalker {
        @Override
        public void apply(final GATKRead read, final ReferenceContext reference,
                          final FeatureContext features) {
            // The window is null when the read has no valid interval, and getBases() is then
            // empty rather than an error: the walker gets a read with no reference under it.
            APPLIED.add(String.format("%s|%d|%s|%d|%s|%s",
                    read.getName(),
                    read.getStart(),
                    read.getCigar().toString(),
                    read.getFlags(),
                    reference.getWindow() == null ? "null" : reference.getWindow().toString(),
                    new String(reference.getBases())));
        }
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("readwalker");
        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        // The Picard call above replaced the deflater factory. It is a static on
        // BlockCompressedOutputStream, and picard's CommandLineProgram installs
        // com.intel.gkl.compression.IntelDeflaterFactory into it, so every BAM written after that
        // call carries GKL bytes rather than the JDK deflater's. The fixture below is only ever
        // read, never reproduced, so no claim in this golden was wrong; what was wrong is that the
        // golden did not say which deflater wrote its input. It says so now, and the factory is
        // pinned before the fixture rather than after it.
        System.out.printf("deflaterafterdict\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        System.out.printf("deflater\t%s%n",
                BlockCompressedOutputStream.getDefaultDeflaterFactory().getClass().getName());

        final Path bam = dir.resolve("reads.bam");
        buildFixture(bam.toFile());
        final Path bai = dir.resolve("reads.bai");

        System.out.println("# ReadWalkerDump: what a ReadWalker hands to apply()");
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(fasta))));
        System.out.printf("fai\t%s%n", ReferenceQueryDump.escape(
                new String(Files.readAllBytes(dir.resolve("ref.fasta.fai")))));
        System.out.printf("bam\t%s%n", Base64.getEncoder().encodeToString(
                Files.readAllBytes(bam)));
        System.out.printf("bai\t%s%n", Base64.getEncoder().encodeToString(
                Files.readAllBytes(bai)));

        // No intervals: the whole file, unplaced reads included.
        traverse("all", bam, fasta, new String[] {});
        traverse("chr1", bam, fasta, new String[] {"-L", "chr1"});
        traverse("chr1:1-60", bam, fasta, new String[] {"-L", "chr1:1-60"});
        traverse("chr1:100-160", bam, fasta, new String[] {"-L", "chr1:100-160"});
        // Two intervals that abut, which the interval parser merges before the query runs.
        traverse("chr1:1-100+101-200", bam, fasta,
                new String[] {"-L", "chr1:1-100", "-L", "chr1:101-200"});
        traverse("chr2", bam, fasta, new String[] {"-L", "chr2"});
        // Every default filter disabled: the reads the WellformedReadFilter was dropping arrive.
        traverse("all-nofilter", bam, fasta,
                new String[] {"--disable-tool-default-read-filters", "true"});
        // One filter added on top of the default.
        traverse("all-nodup", bam, fasta,
                new String[] {"--read-filter", "NotDuplicateReadFilter"});
        // No reference at all: the walker is handed empty ReferenceContexts throughout.
        traverse("all-noref", bam, null, new String[] {});
        // -L unmapped. setTraversalBounds makes a traversal bounded when it has intervals *or*
        // when unmapped reads were asked for, and loadNextIterator runs the interval query first
        // and the unmapped query second, so the unplaced reads are a tail rather than an
        // interleaving. On its own it is a bounded traversal of nothing but that tail, which is a
        // different answer from an unbounded traversal that happens to include the same reads.
        traverse("unmapped", bam, fasta, new String[] {"-L", "unmapped"});
        traverse("unmapped-and-chr1", bam, fasta,
                new String[] {"-L", "unmapped", "-L", "chr1"});
        traverse("unmapped-and-chr2", bam, fasta,
                new String[] {"-L", "unmapped", "-L", "chr2"});
        // The order the arguments are given in cannot matter: the unmapped request is separated
        // out of the interval list before the list is sorted.
        traverse("chr1-and-unmapped", bam, fasta,
                new String[] {"-L", "chr1", "-L", "unmapped"});
        // An unmapped read carrying its mate's position is *not* in the unmapped tail: it is
        // returned by an interval query overlapping that position, which is why this differs from
        // the unmapped-only run by more than the placed reads.
        traverse("unmapped-narrow", bam, fasta,
                new String[] {"-L", "unmapped", "-L", "chr1:1-1"});
    }

    static void traverse(final String label, final Path bam, final Path fasta,
                         final String[] extra) {
        APPLIED.clear();
        final List<String> argv = new ArrayList<>(Arrays.asList("-I", bam.toString()));
        if (fasta != null) {
            argv.add("-R");
            argv.add(fasta.toString());
        }
        argv.addAll(Arrays.asList(extra));

        String summary;
        try {
            final ProbeWalker walker = new ProbeWalker();
            walker.instanceMain(argv.toArray(new String[0]));
            summary = "ok";
        } catch (final Exception | AssertionError e) {
            summary = "E";
        }
        for (int i = 0; i < APPLIED.size(); i++) {
            System.out.printf("apply\t%s\t%d\t%s%n", label, i, APPLIED.get(i));
        }
        System.out.printf("summary\t%s\t%s%n", label, summary);
        System.out.printf("count\t%s\t%d%n", label, APPLIED.size());
    }

    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        for (final String contig : new String[] {"chr1", "chr2"}) {
            dictionary.addSequence(new SAMSequenceRecord(contig, CONTIG_LENGTH));
        }
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);

        final List<SAMRecord> records = new ArrayList<>();
        records.add(mapped(header, "r001", "chr1", 10, "10M", true));
        // Inside the soft-masked stretch: the window comes back upper-cased.
        records.add(mapped(header, "r002", "chr1", 65, "10M", true));
        // Spans the N run at 125-128.
        records.add(mapped(header, "r003", "chr1", 120, "10M", true));
        // A deletion, so the read's span is longer than its bases.
        records.add(mapped(header, "r004", "chr1", 140, "5M10D5M", true));
        // No read group: WellformedReadFilter drops it, and only that filter does.
        records.add(mapped(header, "r005", "chr1", 150, "10M", false));
        // An N operator in the cigar, which WellformedReadFilter also drops.
        records.add(mapped(header, "r006", "chr1", 160, "4M2N4M", true));
        // A duplicate, which the default filters keep and NotDuplicateReadFilter drops.
        final SAMRecord duplicate = mapped(header, "r007", "chr1", 170, "10M", true);
        duplicate.setDuplicateReadFlag(true);
        records.add(duplicate);
        // Mapped with an empty cigar: its alignment end is its start minus one, so the interval
        // is invalid and the ReferenceContext arrives empty.
        records.add(mapped(header, "m001", "chr1", 180, "*", true));
        records.add(mapped(header, "r101", "chr2", 10, "10M", true));
        // Unmapped but parked at its mate's coordinate.
        records.add(unmappedAtMate(header, "u001", "chr2", 20));
        records.add(unplaced(header, "x001"));

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord mapped(final SAMFileHeader header, final String name, final String contig,
                            final int start, final String cigar, final boolean readGroup) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        if (readGroup) {
            record.setAttribute("RG", "rg1");
        }
        return record;
    }

    static SAMRecord unmappedAtMate(final SAMFileHeader header, final String name,
                                    final String contig, final int start) {
        final SAMRecord record = mapped(header, name, contig, start, "*", true);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReadUnmappedFlag(true);
        record.setMateReferenceName(contig);
        record.setMateAlignmentStart(start);
        record.setMappingQuality(0);
        return record;
    }

    static SAMRecord unplaced(final SAMFileHeader header, final String name) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReadUnmappedFlag(true);
        record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setMappingQuality(0);
        record.setReadBases("ACGTACGTAC".getBytes());
        record.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        record.setAttribute("RG", "rg1");
        return record;
    }
}
