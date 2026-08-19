/*
 * GetNormalArtifactData, taken from the reference.
 *
 * The training data for Mutect2's normal-artifact filter, which this port already has. A locus
 * walker with Mutect2's own filters, two rejection rules, a binomial p-value, and a RANDOM
 * DOWNSAMPLE -- seeded, and therefore reproducible.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE ALLELE IS CHOSEN IN THE NORMAL AND COUNTED IN THE TUMOUR. The best non-reference allele
 *     is the normal's most common one, and the tumour count is that allele's, whatever the tumour's
 *     own most common is;
 *   - THE COUNTS HAVE SIX SLOTS, NOT FOUR. Indices 4 and 5 are "before an insertion" and "before a
 *     deletion start", which is how an indel becomes an allele, and index >= 4 is what makes the
 *     record's type INDEL rather than SNV;
 *   - A SITE WITH NO NORMAL ALTERNATE IS SKIPPED, and so is one whose normal alternate is more than
 *     a fifth of the normal pileup: the tool is looking for artefacts, not variants;
 *   - THE DOWNSAMPLE IS RANDOM BUT SEEDED. `Utils.getRandomGenerator()` is `new Random(47382911)`,
 *     shared across the whole run, so the sites that survive depend on how many `nextDouble` calls
 *     came before -- a port that drew at a different moment keeps a different set of sites;
 *   - THE KEEP PROBABILITY HAS A FLOOR of 0.05 and comes from the tumour's binomial p-value, so a
 *     site the tumour supports strongly is always kept and one it does not is usually dropped;
 *   - A TUMOUR ALTERNATE ABOVE HALF THE TUMOUR PILEUP IS SKIPPED, and that test happens AFTER the
 *     random draw, so it consumes a number from the generator either way;
 *   - AND THE TABLE IS SIX COLUMNS with the downsample probability written as a double.
 *
 * Output:
 *
 *     table\t<label>\t<the whole output file, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GetNormalArtifactDataDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.walkers.mutect.GetNormalArtifactData;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GetNormalArtifactDataDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("normalartifact-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        final Path fasta = dir.resolve("ref.fasta");
        Files.write(fasta, ReadWalkerDump.FASTA.getBytes());
        FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("ref.dict")});

        final Path bam = dir.resolve("artifact.bam");
        buildFixture(bam.toFile());

        System.out.println("# GetNormalArtifactDataDump: the normal-artifact training data");

        run("default", dir, fasta, bam, "-normal", "NORMAL", "-L", "chr1:10-19");
        // A higher error probability, which raises every p-value and therefore every keep
        // probability, so more sites survive the draw.
        run("high-error", dir, fasta, bam, "-normal", "NORMAL", "-L", "chr1:10-19",
                "--error-prob", "0.1");
        // The whole forty bases, so the generator is drawn on many more times.
        run("wide", dir, fasta, bam, "-normal", "NORMAL", "-L", "chr1:10-49");
        // A single locus, so the generator is at its first draw when it matters.
        run("one-locus", dir, fasta, bam, "-normal", "NORMAL", "-L", "chr1:12-12");
        // A window whose normal has no alternate at all: every site is skipped.
        // Position 13's reference is `A`, which is what the alternate read carries, so the normal
        // has no alternate there and the site is skipped.
        run("no-alternate", dir, fasta, bam, "-normal", "NORMAL", "-L", "chr1:13-13");
        // The normal named as a sample that is not in the BAM.
        run("unknown-normal", dir, fasta, bam, "-normal", "NOBODY", "-L", "chr1:10-19");

        // A second fixture whose TUMOUR carries no alternate at all. The p-value is then 1, the
        // keep probability falls to its floor of 0.05, and which sites survive is decided entirely
        // by the seeded generator -- which is the only way that draw is observable.
        final Path floorBam = dir.resolve("floor.bam");
        buildFloorFixture(floorBam.toFile());
        run("floor", dir, fasta, floorBam, "-normal", "NORMAL", "-L", "chr1:10-49");
    }

    /** The same normal, and a tumour that carries the reference everywhere. */
    static void buildFloorFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String sample : new String[] {"NORMAL", "TUMOR"}) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg" + sample);
            group.setSample(sample);
            group.setPlatform("ILLUMINA");
            header.addReadGroup(group);
        }

        final StringBuilder matching = new StringBuilder();
        for (int i = 0; i < 10; i++) {
            matching.append("CGTA");
        }
        final StringBuilder alternate = new StringBuilder();
        for (int i = 0; i < 40; i++) {
            alternate.append('A');
        }

        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < 9; i++) {
            records.add(read(header, "n" + i, "rgNORMAL", 10, matching.toString()));
        }
        records.add(read(header, "n9", "rgNORMAL", 10, alternate.toString()));
        for (int i = 0; i < 10; i++) {
            records.add(read(header, "t" + i, "rgTUMOR", 10, matching.toString()));
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /**
     * Two samples over chr1:10-49, whose reference bases repeat `ACGT`.
     *
     * The reads are FORTY bases long, because Mutect2's standard filters include
     * `ReadLengthReadFilter(30, MAX_VALUE)` and a ten-base read never reaches apply at all.
     *
     * Nine normal reads and seven tumour reads carry the reference bases; one normal and three
     * tumour reads carry `A` everywhere. So the normal has one alternate in ten -- under the
     * one-fifth rule -- except at the positions where the reference is itself `A`, which have none
     * and are skipped.
     */
    static void buildFixture(final File bam) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 200));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String sample : new String[] {"NORMAL", "TUMOR"}) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg" + sample);
            group.setSample(sample);
            group.setPlatform("ILLUMINA");
            header.addReadGroup(group);
        }

        // The reference at position 10 is `C`, so a matching read starts there.
        final StringBuilder matching = new StringBuilder();
        for (int i = 0; i < 10; i++) {
            matching.append("CGTA");
        }
        final StringBuilder alternate = new StringBuilder();
        for (int i = 0; i < 40; i++) {
            alternate.append('A');
        }

        final List<SAMRecord> records = new ArrayList<>();
        for (int i = 0; i < 9; i++) {
            records.add(read(header, "n" + i, "rgNORMAL", 10, matching.toString()));
        }
        records.add(read(header, "n9", "rgNORMAL", 10, alternate.toString()));
        for (int i = 0; i < 7; i++) {
            records.add(read(header, "t" + i, "rgTUMOR", 10, matching.toString()));
        }
        for (int i = 7; i < 10; i++) {
            records.add(read(header, "t" + i, "rgTUMOR", 10, alternate.toString()));
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);
        try (final SAMFileWriter writer = factory.makeBAMWriter(header, true, bam)) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    static SAMRecord read(final SAMFileHeader header, final String name, final String group,
                          final int start, final String bases) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString(bases.length() + "M");
        record.setMappingQuality(60);
        record.setReadBases(bases.getBytes());
        final byte[] qualities = new byte[bases.length()];
        java.util.Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setAttribute("RG", group);
        return record;
    }

    static void run(final String label, final Path dir, final Path fasta, final Path bam,
                    final String... extra) throws Exception {
        final Path out = dir.resolve("data-" + label + ".tsv");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(), "-I", bam.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new GetNormalArtifactData().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("table\t%s\t%s%n", label,
                ReferenceQueryDump.escape(new String(Files.readAllBytes(out))));
    }
}
