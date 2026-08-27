/*
 * HaplotypeBasedVariantRecaller's likelihood matrix, taken from the reference.
 *
 * Every allele a haplotype carries, scored against every read that haplotype spans. The tool's
 * own arithmetic is the PairHMM's, which is not what this measures: what it measures is
 * everything around it, which is a haplotype group chosen by a fitness score, a matrix line built
 * per read, a sort, and the ways a line is dropped or comes out wrong.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - A HAPLOTYPE IS A READ WHOSE NAME STARTS `HC_`, and any other record in the haplotype BAM
 *     is passed over however well it fits;
 *   - HAPLOTYPES ARE GROUPED BY IDENTICAL SPAN, a record whose span differs from the group's
 *     first closing the group and opening a new one;
 *   - THE GROUP CHOSEN IS THE ONE THAT CENTRES THE VARIANT BEST, scored as one less twice the
 *     distance from the halfway point, which is why the variant at 1050 is scored against the
 *     group spanning 1000 to 1099 and not against the wider one that also contains it;
 *   - THE ALLELES WRITTEN ARE THE HAPLOTYPES' AND NOT THE VCF'S: the VCF decides only WHICH
 *     positions are scored, so a site written `T -> TCG` comes out `G* GTT`;
 *   - EACH LINE NAMES ITS READ, ITS KEY LENGTH, ITS DUPLICATE AND REVERSE FLAGS AND ITS MAPPING
 *     QUALITY before the likelihoods, and the key length is 0 for a read that is not flow-based;
 *   - THE LINES ARE SORTED BY THE LAST ALLELE'S LIKELIHOOD, descending, AND NOT BY THE BEST OF
 *     THEM: the sort key is overwritten by each allele in turn, so a read with the best first
 *     column can still sort last;
 *   - A READ THAT DOES NOT SPAN THE WHOLE VARIANT IS DROPPED, silently, the bases column coming
 *     out empty and an empty bases column ending the line before it is added;
 *   - A VARIANT INSIDE A READ'S DELETION IS NOT DROPPED, THOUGH, AND COMES OUT WRONG: the cigar
 *     walk subtracts the deletion from the offset without ever refusing, so the offset goes
 *     NEGATIVE inside the element and the base reported is the one that many positions early;
 *   - THE READS ARE TRIMMED TO THE HAPLOTYPE BEFORE ANY OF THIS, so a soft-clipped read's
 *     unclipped offset is its plain offset: the clip is gone by the time the line is built;
 *   - THE HEADER LINE OMITS THE END POSITION for a variant of one base and for a MIXED one,
 *     whatever its length;
 *   - AND AN INTERVAL WITH NO ALLELE IN IT WRITES AN EMPTY FILE rather than none.
 *
 * Output:
 *
 *     vcf\t<label>=<that vcf, escaped>
 *     sam\t<label>=<that bam as sam, without its header, escaped>
 *     out\t<label>=<the whole matrix csv, escaped>
 *     none\t<label>=<what was not written>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HaplotypeBasedVariantRecallerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.variantrecalling.HaplotypeBasedVariantRecaller;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HaplotypeBasedVariantRecallerDump {

    static final int CONTIG_LENGTH = 19980;
    /** The reference is `ACGT` repeated, so a base at any position is known from the position. */
    static final String UNIT = "ACGT";

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

    static List<String> header() {
        return new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsm1"));
    }

    static String site(final int position, final String reference, final String alternate) {
        return "chr1\t" + position + "\t.\t" + reference + "\t" + alternate
                + "\t100.00\tPASS\t.\tGT\t0/1";
    }

    static String vcf(final List<String> sites) {
        final List<String> lines = header();
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("haplotype-based-variant-recaller-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HaplotypeBasedVariantRecallerDump: every allele of a VCF scored "
                + "against every read a haplotype spans");

        final Path fasta = writeReference(dir);

        // The alleles. A SNP, an insertion, a deletion whose END is past its START, and a mixed
        // site, which is the one whose header line omits its end.
        // The reference at these four positions is C, T, T and T, the contig being `ACGT`
        // repeated. The alleles written into the matrix come from the HAPLOTYPES rather than
        // from here, so what this file decides is only WHICH positions are scored.
        final String alleles = vcf(List.of(
                site(1050, "C", "A"),
                // An insertion and a deletion are anchored on the base BEFORE them, so the
                // positions here are one less than the bases the haplotypes change.
                site(1199, "T", "TCG"),
                site(1299, "GTAC", "G"),
                // A mixed site, whose header line omits its end however long it is.
                site(1400, "TA", "GC,TACG")));
        final Path vcfPath = index(write(dir, "alleles.vcf", alleles));
        System.out.printf("vcf\talleles=%s%n", ReferenceQueryDump.escape(alleles));

        final Path reads = dir.resolve("reads.bam").toAbsolutePath();
        writeReads(reads.toFile());
        System.out.printf("sam\treads=%s%n", ReferenceQueryDump.escape(asSam(reads)));

        final Path haplotypes = dir.resolve("haplotypes.bam").toAbsolutePath();
        writeHaplotypes(haplotypes.toFile());
        System.out.printf("sam\thaplotypes=%s%n",
                ReferenceQueryDump.escape(asSam(haplotypes)));

        run(dir, "whole-reference", vcfPath, reads, haplotypes, fasta, List.of());
        // One interval, which is what limits the walk rather than the VCF.
        run(dir, "one-interval", vcfPath, reads, haplotypes, fasta,
                List.of("-L", "chr1:1000-1250"));
        // An interval no allele falls in, which writes an empty file rather than none.
        run(dir, "empty-interval", vcfPath, reads, haplotypes, fasta,
                List.of("-L", "chr1:5000-5100"));
        // The mapping-quality filter, which is what takes reads out of the matrix.
        run(dir, "mapping-quality-30", vcfPath, reads, haplotypes, fasta,
                List.of("--minimum-mapping-quality", "30"));
    }

    /**
     * The reads. Six over the first haplotype's span and three over the second's.
     *
     * One does not span its variant, one has the variant inside a deletion, one is on the reverse
     * strand, one is a duplicate, and one is soft-clipped so its unclipped offset differs from its
     * offset.
     */
    static void writeReads(final File file) {
        final SAMFileHeader header = readHeader();
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(
                             header, true, file)) {
            add(writer, header, "r-spans-all", 1000, "120M", false, false, 60);
            // Ends before the insertion at 1200, so its bases column comes out empty there.
            add(writer, header, "r-short", 1000, "60M", false, false, 60);
            // A deletion over the SNP at 1050, which is the other way a line loses its bases.
            add(writer, header, "r-deletion", 1000, "45M10D65M", false, false, 60);
            add(writer, header, "r-reverse", 1000, "120M", true, false, 60);
            add(writer, header, "r-duplicate", 1000, "120M", false, true, 60);
            // Soft-clipped at the front, so its unclipped start is ten before its start.
            add(writer, header, "r-clipped", 1000, "10S110M", false, false, 20);
            // The only read over the insertion the far group carries.
            add(writer, header, "r-far-region", 1150, "120M", false, false, 60);
            add(writer, header, "r-second-region", 1250, "220M", false, false, 60);
            add(writer, header, "r-second-reverse", 1250, "220M", true, false, 60);
            add(writer, header, "r-second-low-mq", 1250, "220M", false, false, 25);
        }
    }

    /**
     * The haplotypes. Two groups over the first region, one centred on the SNP and one not, and
     * one group over the second.
     *
     * A tenth record does NOT start with `HC_`, so it is passed over however well it fits.
     */
    static void writeHaplotypes(final File file) {
        final SAMFileHeader header = readHeader();
        // Written UNSORTED and sorted on the way out: the two groups are added in fitness order
        // rather than in coordinate order, and a presorted writer refuses that.
        try (final SAMFileWriter writer =
                     new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(
                             header, false, file)) {
            // The off-centre group: the variant at 1050 sits near its start, and one of its
            // haplotypes carries the insertion at 1200.
            haplotype(writer, header, "HC_far_ref", 1040, "200M", -1);
            haplotype(writer, header, "HC_far_alt", 1040, "200M", 20);
            haplotype(writer, header, "HC_far_ins", 1040, "160M2I40M", -1);
            // The centred group, which the fitness score prefers for the variant at 1050.
            haplotype(writer, header, "HC_near_ref", 1000, "100M", -1);
            haplotype(writer, header, "HC_near_alt", 1000, "100M", 50);
            // A record in the haplotype BAM that is not a haplotype, however well it fits.
            haplotype(writer, header, "not_a_haplotype", 1000, "100M", 10);
            haplotype(writer, header, "HC_second_ref", 1250, "250M", -1);
            haplotype(writer, header, "HC_second_alt", 1250, "250M", 150);
            haplotype(writer, header, "HC_second_del", 1250, "50M3D197M", -1);
        }
    }

    static SAMFileHeader readHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sm1");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        return header;
    }

    static void add(final SAMFileWriter writer, final SAMFileHeader header, final String name,
                    final int start, final String cigar, final boolean reverse,
                    final boolean duplicate, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        // The bases are the reference's, so every likelihood is a reference match except where a
        // haplotype carries a substitution.
        int length = 0;
        for (final htsjdk.samtools.CigarElement element : record.getCigar()) {
            if (element.getOperator().consumesReadBases()) {
                length += element.getLength();
            }
        }
        record.setReadBases(referenceBases(start, length));
        final byte[] quality = new byte[length];
        Arrays.fill(quality, (byte) 35);
        record.setBaseQualities(quality);
        record.setReadNegativeStrandFlag(reverse);
        record.setDuplicateReadFlag(duplicate);
        record.setMappingQuality(mappingQuality);
        record.setAttribute("RG", "rg1");
        writer.addAlignment(record);
    }

    /**
     * One haplotype: the reference over its span, with one base substituted when asked for.
     */
    static void haplotype(final SAMFileWriter writer, final SAMFileHeader header,
                          final String name, final int start, final String cigar,
                          final int substitutionOffset) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        // The haplotype's bases are the reference's, walked through its own cigar, so an
        // insertion is two bases the reference does not have and a deletion is bases left out.
        final byte[] bases = walk(record, start);
        final int length = bases.length;
        if (substitutionOffset >= 0) {
            bases[substitutionOffset] = bases[substitutionOffset] == 'A' ? (byte) 'G' : (byte) 'A';
        }
        record.setReadBases(bases);
        final byte[] quality = new byte[length];
        Arrays.fill(quality, (byte) 40);
        record.setBaseQualities(quality);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        writer.addAlignment(record);
    }

    /** The reference bases a cigar reads, with `T` for every inserted base. */
    static byte[] walk(final SAMRecord record, final int start) {
        final java.io.ByteArrayOutputStream bases = new java.io.ByteArrayOutputStream();
        int position = start;
        for (final htsjdk.samtools.CigarElement element : record.getCigar()) {
            final htsjdk.samtools.CigarOperator operator = element.getOperator();
            if (operator.consumesReadBases() && operator.consumesReferenceBases()) {
                bases.write(referenceBases(position, element.getLength()), 0, element.getLength());
                position += element.getLength();
            } else if (operator.consumesReadBases()) {
                for (int i = 0; i < element.getLength(); i++) {
                    bases.write('T');
                }
            } else if (operator.consumesReferenceBases()) {
                position += element.getLength();
            }
        }
        return bases.toByteArray();
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static Path index(final Path path) throws Exception {
        htsjdk.tribble.index.IndexFactory.createLinearIndex(path.toFile(),
                new htsjdk.variant.vcf.VCFCodec()).writeBasedOnFeatureFile(path.toFile());
        return path;
    }

    /** A BAM read back as SAM text, without its header, which is what the golden carries. */
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

    static void run(final Path dir, final String label, final Path alleles, final Path reads,
                    final Path haplotypes, final Path fasta, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".csv");
        final List<String> argv = new ArrayList<>(List.of(
                "--alleles-file-vcf", alleles.toString(),
                "--haplotypes-file-bam", haplotypes.toString(),
                "--matrix-file-csv", out.toString(),
                "-I", reads.toString(),
                "-R", fasta.toString()));
        argv.addAll(extra);
        try {
            new HaplotypeBasedVariantRecaller().instanceMain(argv.toArray(new String[0]));
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

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            bases.append(new String(referenceBases(i * 60 + 1, 60), StandardCharsets.UTF_8))
                    .append("\n");
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
