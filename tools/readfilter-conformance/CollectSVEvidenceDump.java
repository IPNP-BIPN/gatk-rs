/*
 * CollectSVEvidence's four evidence files, taken from the reference.
 *
 * What one BAM contributes to structural-variant calling: discordant pairs, split-read positions,
 * per-site allele depths and per-interval read counts. One traversal writes all four, and each has
 * its own rule for which reads it will look at.
 *
 * Ten behaviours this is built to catch.
 *
 *   - A DISCORDANT PAIR IS WRITTEN ONCE, BY ONE OF ITS TWO READS: the one whose contig index is
 *     smaller, or whose start is smaller, and at an EQUAL start the FIRST of the two seen, tracked
 *     by name within the locus;
 *   - THE STRAND FIELD HOLDS `!isReverseStrand` AND THE CODEC PRINTS `+` FOR TRUE, so the two
 *     negations cancel and the output shows each read's own strand: the record written for the
 *     same-start pair prints `+` for its forward read and `-` for its reverse mate;
 *   - A READ IS SPLIT ONLY IF EXACTLY ONE END IS SOFT-CLIPPED: clipped at both ends it is not
 *     counted at all;
 *   - WHICH END DECIDES THE DIRECTION AND THE POSITION: a leading match gives a RIGHT position at
 *     start plus every reference-consuming length, and a leading clip gives a LEFT position at the
 *     start itself;
 *   - THE MATCH LENGTH COUNTS DELETIONS, because it sums every operator that consumes reference;
 *   - SPLIT POSITIONS ARE COUNTED, so two reads clipped at the same place print a count of two,
 *     and the same position in the two directions stays two records;
 *   - THE SITE DEPTH READS ONLY BIALLELIC SNPS at new loci, so an indel, a triallelic site and a
 *     repeated position are all skipped;
 *   - --site-depth-min-mapq AND --site-depth-min-baseq FILTER DIFFERENT THINGS, one the whole read
 *     and one a single base;
 *   - AN INTERVAL WITH NO READS IS STILL WRITTEN, with a count of zero;
 *   - AND EVERY OUTPUT REFUSES A FILE NAME IT COULD NOT READ BACK, each with its own message
 *     naming its own three extensions.
 *
 * Output:
 *
 *     bam\treads=<one line per record: name flags contig start mapq cigar mate-contig mate-start>
 *     vcf\tsites=<the site-depth vcf, escaped>
 *     bed\tintervals=<the depth-evidence intervals, escaped>
 *     out\t<label>.<kind>=<that evidence file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CollectSVEvidenceDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.tools.walkers.sv.CollectSVEvidence;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectSVEvidenceDump {

    static final int CONTIG_LENGTH = 199980;
    static final String SAMPLE = "NA1";

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                new SAMSequenceRecord("chr2", CONTIG_LENGTH))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample(SAMPLE);
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        return header;
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String contig,
                          final int start, final String cigar, final int mappingQuality,
                          final int baseQuality, final int length) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setMappingQuality(mappingQuality);
        final StringBuilder bases = new StringBuilder();
        while (bases.length() < length) {
            bases.append("ACGT");
        }
        record.setReadBases(bases.substring(0, length).getBytes());
        final byte[] qualities = new byte[length];
        Arrays.fill(qualities, (byte) baseQuality);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", "rg1");
        return record;
    }

    /** A read of one mate of a pair. `proper` decides whether the PE writer will look at it. */
    static SAMRecord paired(final SAMFileHeader header, final String name, final String contig,
                            final int start, final String cigar, final boolean first,
                            final String mateContig, final int mateStart, final boolean proper,
                            final boolean reverse, final boolean mateReverse) {
        final SAMRecord record = read(header, name, contig, start, cigar, 60, 30, 20);
        record.setReadPairedFlag(true);
        record.setProperPairFlag(proper);
        record.setFirstOfPairFlag(first);
        record.setSecondOfPairFlag(!first);
        record.setMateReferenceName(mateContig);
        record.setMateAlignmentStart(mateStart);
        record.setReadNegativeStrandFlag(reverse);
        record.setMateNegativeStrandFlag(mateReverse);
        return record;
    }

    /**
     * A read that is PROPERLY PAIRED and so never reaches the discordant writer.
     *
     * Every read in the fixture has to be one, because an UNPAIRED read reports
     * `isProperlyPaired() == false` and the discordant writer then asks it for its mate: see the
     * `unpaired` run.
     */
    static SAMRecord solo(final SAMFileHeader header, final String name, final int start,
                          final String cigar, final int mappingQuality, final int baseQuality,
                          final int length) {
        final SAMRecord record = read(header, name, "chr1", start, cigar, mappingQuality,
                baseQuality, length);
        record.setReadPairedFlag(true);
        record.setProperPairFlag(true);
        record.setFirstOfPairFlag(true);
        record.setMateReferenceName("chr1");
        record.setMateAlignmentStart(start + 200);
        return record;
    }

    static List<SAMRecord> buildReads(final SAMFileHeader header) {
        final List<SAMRecord> records = new ArrayList<>();

        // A properly paired pair, which the discordant writer never looks at.
        records.add(paired(header, "ok", "chr1", 1000, "20M", true, "chr1", 1200, true,
                false, true));
        records.add(paired(header, "ok", "chr1", 1200, "20M", false, "chr1", 1000, true,
                true, false));

        // A discordant pair on one contig: only the read with the smaller start writes it, and the
        // strands are written inverted.
        records.add(paired(header, "disc", "chr1", 2000, "20M", true, "chr1", 2500, false,
                false, false));
        records.add(paired(header, "disc", "chr1", 2500, "20M", false, "chr1", 2000, false,
                false, false));

        // A discordant pair at the SAME start, where the first one seen writes it and the second
        // is recognised by name and dropped.
        records.add(paired(header, "same", "chr1", 3000, "20M", true, "chr1", 3000, false,
                false, true));
        records.add(paired(header, "same", "chr1", 3000, "20M", false, "chr1", 3000, false,
                true, false));

        // Two DIFFERENT pairs at the same start, which are two names and so two records.
        records.add(paired(header, "sameA", "chr1", 4000, "20M", true, "chr1", 4000, false,
                false, false));
        records.add(paired(header, "sameB", "chr1", 4000, "20M", true, "chr1", 4000, false,
                false, false));

        // A pair across contigs, written by the read on the SMALLER contig index only.
        records.add(paired(header, "cross", "chr1", 5000, "20M", true, "chr2", 5000, false,
                false, false));

        // Split reads. A leading match gives a RIGHT position past the alignment; a leading clip
        // gives a LEFT position at the start.
        records.add(solo(header, "clipRight", 10000, "20M5S", 60, 30, 25));
        records.add(solo(header, "clipLeft", 11000, "5S20M", 60, 30, 25));
        // Two reads clipped at the same place, which is a count of two.
        records.add(solo(header, "twoA", 12000, "20M5S", 60, 30, 25));
        records.add(solo(header, "twoB", 12000, "20M5S", 60, 30, 25));
        // The same position reached from the two directions, which stays two records.
        records.add(solo(header, "bothA", 13000, "20M5S", 60, 30, 25));
        records.add(solo(header, "bothB", 13020, "5S20M", 60, 30, 25));
        // Clipped at BOTH ends, which is not a split read at all.
        records.add(solo(header, "clipBoth", 14000, "5S20M5S", 60, 30, 30));
        // A deletion inside the alignment, which the match length counts.
        records.add(solo(header, "withDeletion", 15000, "10M5D10M5S", 60, 30, 25));
        // A supplementary alignment, which neither writer looks at.
        final SAMRecord supplementary = solo(header, "supp", 16000, "20M5S", 60, 30, 25);
        supplementary.setSupplementaryAlignmentFlag(true);
        records.add(supplementary);
        // A duplicate, which the default read filters remove.
        final SAMRecord duplicate = solo(header, "dup", 17000, "20M5S", 60, 30, 25);
        duplicate.setDuplicateReadFlag(true);
        records.add(duplicate);

        // Reads over the site-depth loci and the depth intervals.
        records.add(solo(header, "depth1", 20000, "20M", 60, 30, 20));
        records.add(solo(header, "depth2", 20000, "20M", 60, 30, 20));
        // One below the mapping-quality floor the runs use, and one whose bases are all below the
        // base-quality floor.
        records.add(solo(header, "lowMapq", 20000, "20M", 5, 30, 20));
        records.add(solo(header, "lowBaseq", 20000, "20M", 60, 5, 20));
        // A read over the second depth interval and none over the third.
        records.add(solo(header, "depth3", 30000, "20M", 60, 30, 20));

        records.sort((a, b) -> {
            final int contig = Integer.compare(a.getReferenceIndex(), b.getReferenceIndex());
            return contig != 0 ? contig : Integer.compare(a.getAlignmentStart(),
                    b.getAlignmentStart());
        });
        return records;
    }

    /**
     * The reads as text, so the golden says what was fed in.
     *
     * The BASES and the BASE QUALITIES are reported too: without them the site-depth counts cannot
     * be recomputed from the golden, only from the source of this file.
     */
    static String describe(final List<SAMRecord> records) {
        final StringBuilder text = new StringBuilder();
        for (final SAMRecord record : records) {
            final StringBuilder qualities = new StringBuilder();
            for (final byte quality : record.getBaseQualities()) {
                if (qualities.length() > 0) {
                    qualities.append(',');
                }
                qualities.append(quality);
            }
            text.append(record.getReadName()).append('\t')
                    .append(record.getFlags()).append('\t')
                    .append(record.getReferenceName()).append('\t')
                    .append(record.getAlignmentStart()).append('\t')
                    .append(record.getMappingQuality()).append('\t')
                    .append(record.getCigarString()).append('\t')
                    .append(record.getReadPairedFlag() ? record.getMateReferenceName() : ".")
                    .append('\t')
                    .append(record.getReadPairedFlag() ? record.getMateAlignmentStart() : 0)
                    .append('\t')
                    .append(record.getReadString()).append('\t')
                    .append(qualities)
                    .append('\n');
        }
        return text.toString();
    }

    /** Biallelic SNPs, plus the three kinds of site the iterator walks past. */
    static String buildSites() {
        return String.join("\n",
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##contig=<ID=chr2,length=" + CONTIG_LENGTH + ">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                // Taken: a biallelic SNP.
                "chr1\t20005\t.\tA\tC\t.\t.\t.",
                // Skipped: the same position again.
                "chr1\t20005\t.\tA\tG\t.\t.\t.",
                // Skipped: triallelic.
                "chr1\t20008\t.\tA\tC,G\t.\t.\t.",
                // Skipped: an indel.
                "chr1\t20010\t.\tAC\tA\t.\t.\t.",
                // Taken.
                "chr1\t20012\t.\tA\tT\t.\t.\t.",
                // Taken, and no read covers it.
                "chr1\t50000\t.\tA\tT\t.\t.\t.",
                "");
    }

    /** Three intervals, the last of which no read reaches. */
    static String buildIntervals() {
        return String.join("\n",
                "chr1\t19999\t20100",
                "chr1\t29999\t30100",
                "chr1\t99999\t100100",
                "");
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("collect-sv-evidence-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CollectSVEvidenceDump: what one BAM contributes to structural-variant"
                + " calling");

        final SAMFileHeader header = header();
        final List<SAMRecord> records = buildReads(header);
        final Path bam = dir.resolve("reads.bam");
        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("bam\treads=%s%n", ReferenceQueryDump.escape(describe(records)));

        final String sites = buildSites();
        final String intervals = buildIntervals();
        final Path sitesPath = write(dir, "sites.vcf", sites);
        final Path intervalsPath = write(dir, "intervals.bed", intervals);
        System.out.printf("vcf\tsites=%s%n", ReferenceQueryDump.escape(sites));
        System.out.printf("bed\tintervals=%s%n", ReferenceQueryDump.escape(intervals));

        run(dir, "default", bam, List.of(
                "-PE", dir.resolve("out-default.pe.txt").toString(),
                "-SR", dir.resolve("out-default.sr.txt").toString(),
                "-SD", dir.resolve("out-default.sd.txt").toString(),
                "-F", sitesPath.toString(),
                "-RD", dir.resolve("out-default.rd.txt").toString(),
                "-DI", intervalsPath.toString()));
        // A mapping-quality floor above the low read's, for both counters.
        run(dir, "high-mapq", bam, List.of(
                "-SD", dir.resolve("out-high-mapq.sd.txt").toString(),
                "-F", sitesPath.toString(),
                "-RD", dir.resolve("out-high-mapq.rd.txt").toString(),
                "-DI", intervalsPath.toString(),
                "--site-depth-min-mapq", "30",
                "--depth-evidence-min-mapq", "30"));
        // A base-quality floor of zero. The default is 20 and the low read's bases are 5, so it is
        // ALREADY excluded by default: lowering the floor is what makes its base appear, and it
        // filters a base rather than a read, which is what separates it from the mapping-quality
        // floor.
        run(dir, "low-baseq", bam, List.of(
                "-SD", dir.resolve("out-low-baseq.sd.txt").toString(),
                "-F", sitesPath.toString(),
                "--site-depth-min-baseq", "0"));
        // Each output on its own, which is how the four are shown to be independent.
        run(dir, "pe-only", bam, List.of(
                "-PE", dir.resolve("out-pe-only.pe.txt").toString()));
        run(dir, "sr-only", bam, List.of(
                "-SR", dir.resolve("out-sr-only.sr.txt").toString()));

        // A BAM holding one UNPAIRED read, with the discordant writer asked for. An unpaired read
        // reports isProperlyPaired() == false, so the writer asks it for a mate it does not have.
        final SAMFileHeader loneHeader = header();
        final SAMRecord lone = read(loneHeader, "lone", "chr1", 1000, "20M", 60, 30, 20);
        final Path loneBam = dir.resolve("unpaired.bam");
        try (final SAMFileWriter writer = factory.makeBAMWriter(loneHeader, true,
                loneBam.toFile())) {
            writer.addAlignment(lone);
        }
        run(dir, "unpaired", loneBam, List.of(
                "-PE", dir.resolve("out-unpaired.pe.txt").toString()));

        // No output at all.
        run(dir, "no-output", bam, List.of());
        // A file name each writer could not read back.
        run(dir, "bad-pe-name", bam, List.of("-PE", dir.resolve("wrong.txt").toString()));
        run(dir, "bad-sr-name", bam, List.of("-SR", dir.resolve("wrong.txt").toString()));
        run(dir, "bad-sd-name", bam, List.of(
                "-SD", dir.resolve("wrong.txt").toString(), "-F", sitesPath.toString()));
        run(dir, "bad-rd-name", bam, List.of(
                "-RD", dir.resolve("wrong.txt").toString(), "-DI", intervalsPath.toString()));
        // An intervals file with nothing in it.
        final Path empty = write(dir, "empty.bed", "");
        run(dir, "empty-intervals", bam, List.of(
                "-RD", dir.resolve("out-empty.rd.txt").toString(), "-DI", empty.toString()));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path bam, final List<String> extra)
            throws Exception {
        final List<String> argv = new ArrayList<>(List.of(
                "-I", bam.toString(),
                "--sample-name", SAMPLE));
        argv.addAll(extra);
        try {
            new CollectSVEvidence().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        for (final String kind : new String[] {"pe", "sr", "sd", "rd"}) {
            final Path out = dir.resolve("out-" + label + "." + kind + ".txt");
            if (!Files.exists(out)) {
                continue;
            }
            System.out.printf("out\t%s.%s=%s%n", label, kind,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
