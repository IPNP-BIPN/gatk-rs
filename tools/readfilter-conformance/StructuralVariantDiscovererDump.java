/*
 * StructuralVariantDiscoverer's calls, taken from the reference.
 *
 * How the split alignments of an assembled contig become a structural variant. A contig that aligns
 * in two pieces with a gap between them is a deletion; the same two pieces in the other order or on
 * the other strand are something else; and which it is comes from the SIGNATURE of the pair rather
 * than from any argument.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE READS MUST BE QUERYNAME-SORTED, refused outright otherwise, because the tool gathers a
 *     contig's alignments by walking consecutive records of the same name;
 *   - A CONTIG ALIGNING IN TWO PIECES WITH A REFERENCE GAP IS A DELETION, whose length is the gap;
 *   - THE SAME TWO PIECES WITH THE SECOND ON THE OTHER STRAND PRODUCE NOTHING: a strand flip
 *     alone is not an inversion signature, and the contig is dropped in silence;
 *   - A CONTIG WITH ONE ALIGNMENT PRODUCES NOTHING, whatever else is in the file;
 *   - A SECONDARY ALIGNMENT IS FILTERED and an UNMAPPED one too, before any signature is read;
 *   - THE SAMPLE ID COMES FROM THE READ GROUP, so the output's one sample is the header's;
 *   - THE OUTPUT IS A VCF WHOSE ALTERNATES ARE SYMBOLIC, one record per adjacency rather than one
 *     per contig;
 *   - AND A CONTIG WHOSE PIECES OVERLAP ON THE REFERENCE IS A TANDEM DUPLICATION rather than
 *     nothing: the overlap becomes the repeat unit, and the call is <DUP> with a duplication
 *     number of 1,2.
 *
 * Output:
 *
 *     bam\treads=<one line per record: name flags contig start mapq cigar>
 *     out\t<label>=<the whole output vcf without its header, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: StructuralVariantDiscovererDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.StructuralVariantDiscoverer;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class StructuralVariantDiscovererDump {

    static final int CONTIG_LENGTH = 199980;

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(order);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sampleA");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        return header;
    }

    /** One alignment of an assembled contig. */
    static SAMRecord piece(final SAMFileHeader header, final String name, final int start,
                           final String cigar, final boolean reverse, final boolean supplementary,
                           final int length) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        final StringBuilder bases = new StringBuilder();
        while (bases.length() < length) {
            bases.append("ACGT");
        }
        record.setReadBases(bases.substring(0, length).getBytes(StandardCharsets.UTF_8));
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setReadNegativeStrandFlag(reverse);
        record.setSupplementaryAlignmentFlag(supplementary);
        record.setAttribute("RG", "rg1");
        return record;
    }

    static String describe(final List<SAMRecord> records) {
        final StringBuilder text = new StringBuilder();
        for (final SAMRecord record : records) {
            text.append(record.getReadName()).append('\t')
                    .append(record.getFlags()).append('\t')
                    .append(record.getReferenceName()).append('\t')
                    .append(record.getAlignmentStart()).append('\t')
                    .append(record.getMappingQuality()).append('\t')
                    .append(record.getCigarString()).append('\n');
        }
        return text.toString();
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sv-discoverer-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# StructuralVariantDiscovererDump: how split alignments become a "
                + "structural variant");

        final Path fasta = writeReference(dir);
        final SAMFileHeader queryname = header(SAMFileHeader.SortOrder.queryname);

        final List<SAMRecord> records = new ArrayList<>();
        // A deletion: two pieces of one contig with a 500-base gap on the reference.
        records.add(piece(queryname, "ctg-del", 10000, "100M100S", false, false, 200));
        records.add(piece(queryname, "ctg-del", 10600, "100S100M", false, true, 200));
        // The second piece on the other strand, which turns out to produce no call at all.
        records.add(piece(queryname, "ctg-inv", 20000, "100M100S", false, false, 200));
        records.add(piece(queryname, "ctg-inv", 20600, "100S100M", true, true, 200));
        // One alignment only, which is no adjacency at all.
        records.add(piece(queryname, "ctg-single", 30000, "200M", false, false, 200));
        // Two pieces that OVERLAP on the reference rather than leaving a gap, which is read as a
        // tandem duplication.
        records.add(piece(queryname, "ctg-overlap", 40000, "100M100S", false, false, 200));
        records.add(piece(queryname, "ctg-overlap", 40050, "100S100M", false, true, 200));
        // A secondary alignment and an unmapped record, both filtered before any signature.
        final SAMRecord secondary = piece(queryname, "ctg-secondary", 50000, "200M", false, false,
                200);
        secondary.setSecondaryAlignment(true);
        records.add(secondary);
        final SAMRecord unmapped = piece(queryname, "ctg-unmapped", 60000, "200M", false, false,
                200);
        unmapped.setReadUnmappedFlag(true);
        records.add(unmapped);

        final Path bam = dir.resolve("contigs.bam");
        // presorted=false, so htsjdk puts the records in queryname order itself: the fixture is
        // written in the order the behaviours read best, which is not alphabetical.
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(queryname, false, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("bam\treads=%s%n", ReferenceQueryDump.escape(describe(records)));

        run(dir, "default", bam, fasta, List.of());

        // The same records in coordinate order, which the tool refuses.
        final SAMFileHeader coordinate = header(SAMFileHeader.SortOrder.coordinate);
        final Path sorted = dir.resolve("coordinate.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(coordinate, false, sorted.toFile())) {
            for (final SAMRecord record : records) {
                final SAMRecord copy = record.deepCopy();
                copy.setHeader(coordinate);
                writer.addAlignment(copy);
            }
        }
        run(dir, "coordinate-sorted", sorted, fasta, List.of());
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
            new StructuralVariantDiscoverer().instanceMain(argv.toArray(new String[0]));
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
            System.out.printf("none\t%s=no output file%n", label);
            return;
        }
        final StringBuilder body = new StringBuilder();
        for (final String line : Files.readString(out).split("\n", -1)) {
            if (!line.startsWith("##") && !line.isEmpty()) {
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
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
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
