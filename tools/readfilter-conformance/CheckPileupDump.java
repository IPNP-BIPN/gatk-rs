/*
 * CheckPileup's output, taken from the reference.
 *
 * A LocusWalker that compares GATK's own pileup against a samtools mpileup file, locus by locus, and
 * prints every disagreement. The truth file goes through SAMPileupCodec, which the previous suite
 * measured.
 *
 * Nine behaviours this is built to catch.
 *
 *   - ITS FILTERS ARE THE LOCUS WALKER'S PLUS THREE OF SAMTOOLS', NOT_DUPLICATE,
 *     PASSES_VENDOR_QUALITY_CHECK and NOT_SECONDARY_ALIGNMENT, so a duplicate is absent from GATK's
 *     pileup and its absence is a size disagreement rather than a base one;
 *   - IT FIXES OVERLAPPING PAIRS BY DEFAULT, the way samtools does: where two mates cover one locus
 *     with the same base one quality is raised and the other set to ZERO, so the pileup GATK reports
 *     is not the qualities the reads carry. --ignore-overlaps turns that off;
 *   - A LOCUS THE TRUTH FILE DOES NOT COVER IS A REFUSAL, and the line it prints first is GATK's own
 *     pileup, so the file is written before the exception is thrown;
 *   - A DISAGREEMENT PRINTS `<gatk> vs. <truth>` and then refuses, with a message naming WHICH of
 *     the four comparisons failed: size, location, bases or quals, in that order;
 *   - THE FOUR COMPARISONS STOP AT THE FIRST FAILURE, so a locus that differs in size never reports
 *     its bases;
 *   - THE BASES ARE COMPARED CASE-INSENSITIVELY AND THE QUALITIES ARE NOT, because the bases go
 *     through toUpperCase and the qualities are compared as the characters they print as;
 *   - --continue-after-error PRINTS THE SAME LINES AND CARRIES ON, so the file of a failing run and
 *     the file of a continuing one differ by everything after the first failure;
 *   - THE SUMMARY IS RETURNED, NOT PRINTED TO THE FILE: "Validated %d sites covered by %d bases",
 *     where the bases are counted AFTER the filters and the overlap fixing;
 *   - AND THE COUNTERS ADVANCE EVEN FOR A LOCUS THAT DISAGREED, because `nLoci++` sits after the
 *     comparison rather than inside its success branch.
 *
 * Output:
 *
 *     reference\t<the reference bases of chr1>
 *     fixture\t<label>\t<the input BAM, base64>
 *     fixtureindex\t<label>\t<the index, base64>
 *     truth\t<label>\t<the pileup file, escaped>
 *     report\t<label>\t<the whole output file, escaped>
 *     summary\t<label>\t<what onTraversalSuccess returned>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CheckPileupDump
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
import org.broadinstitute.hellbender.tools.walkers.qc.CheckPileup;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CheckPileupDump {

    /** Forty bases, with a run of As where the reads sit. */
    static final String CHR1 = "ACGTACGTACGTTTTTGGGGCCCCAAAAACGTACGTACGT";

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final Path dir = Path.of("checkpileup-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CheckPileupDump: CheckPileup's output, from the reference");
        System.out.printf("reference\t%s%n", CHR1);

        final Path fasta = writeReference(dir);

        // Two reads over 21..24, whose bases are the reference's.
        buildBam(dir, "plain", new String[][] {
                {"r1", "21", "4M", "CCCC", "0", "IIII"},
                {"r2", "21", "4M", "CCCC", "0", "JJJJ"},
        });
        // The same, with the second read a duplicate: the filters drop it.
        buildBam(dir, "duplicate", new String[][] {
                {"r1", "21", "4M", "CCCC", "0", "IIII"},
                {"r2", "21", "4M", "CCCC", "1024", "JJJJ"},
        });
        // An overlapping pair: two mates of one fragment over the same four bases.
        buildBam(dir, "overlapping", new String[][] {
                {"pair", "21", "4M", "CCCC", "99", "IIII"},
                {"pair", "21", "4M", "CCCC", "147", "JJJJ"},
        });

        // A truth file that agrees with the plain fixture.
        truth(dir, "agrees", new String[] {
                "chr1\t21\tC\t2\t..\tIJ",
                "chr1\t22\tC\t2\t..\tIJ",
                "chr1\t23\tC\t2\t..\tIJ",
                "chr1\t24\tC\t2\t..\tIJ",
        });
        // One whose first locus has one read too few.
        truth(dir, "wrong-size", new String[] {
                "chr1\t21\tC\t1\t.\tI",
                "chr1\t22\tC\t2\t..\tIJ",
                "chr1\t23\tC\t2\t..\tIJ",
                "chr1\t24\tC\t2\t..\tIJ",
        });
        // One whose bases differ at the first locus, in lower case at the second.
        truth(dir, "wrong-bases", new String[] {
                "chr1\t21\tC\t2\tAA\tIJ",
                "chr1\t22\tC\t2\tcc\tIJ",
                "chr1\t23\tC\t2\t..\tIJ",
                "chr1\t24\tC\t2\t..\tIJ",
        });
        // One whose qualities differ at the first locus.
        truth(dir, "wrong-quals", new String[] {
                "chr1\t21\tC\t2\t..\tII",
                "chr1\t22\tC\t2\t..\tIJ",
                "chr1\t23\tC\t2\t..\tIJ",
                "chr1\t24\tC\t2\t..\tIJ",
        });
        // One that stops before the last locus.
        truth(dir, "incomplete", new String[] {
                "chr1\t21\tC\t2\t..\tIJ",
                "chr1\t22\tC\t2\t..\tIJ",
                "chr1\t23\tC\t2\t..\tIJ",
        });
        // What the overlapping pair looks like once samtools has tweaked it.
        truth(dir, "overlapping-fixed", new String[] {
                "chr1\t21\tC\t2\t..\t^!",
                "chr1\t22\tC\t2\t..\t^!",
                "chr1\t23\tC\t2\t..\t^!",
                "chr1\t24\tC\t2\t..\t^!",
        });

        run(dir, fasta, "agrees", "plain", "agrees", new String[] {});
        run(dir, fasta, "wrong-size", "plain", "wrong-size", new String[] {});
        run(dir, fasta, "wrong-bases", "plain", "wrong-bases", new String[] {});
        run(dir, fasta, "wrong-quals", "plain", "wrong-quals", new String[] {});
        run(dir, fasta, "incomplete", "plain", "incomplete", new String[] {});
        // The same failures with the tool told to carry on, where the file holds every one of them.
        run(dir, fasta, "wrong-bases-continue", "plain", "wrong-bases",
                new String[] {"--continue-after-error", "true"});
        run(dir, fasta, "incomplete-continue", "plain", "incomplete",
                new String[] {"--continue-after-error", "true"});
        // A duplicate, which the filters drop and the truth file therefore over-counts.
        run(dir, fasta, "duplicate", "duplicate", "agrees",
                new String[] {"--continue-after-error", "true"});
        // The overlapping pair, with and without the fixing.
        run(dir, fasta, "overlapping", "overlapping", "agrees",
                new String[] {"--continue-after-error", "true"});
        run(dir, fasta, "overlapping-ignored", "overlapping", "agrees",
                new String[] {"--ignore-overlaps", "true", "--continue-after-error", "true"});
    }

    static void buildBam(final Path dir, final String label, final String[][] reads)
            throws Exception {
        final Path bam = dir.resolve(label + ".bam");
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            for (final String[] spec : reads) {
                writer.addAlignment(read(header, spec));
            }
        }
        System.out.printf("fixture\t%s\t%s%n", label, RecordTransformDump.base64(bam));
        final Path index = dir.resolve(label + ".bai");
        System.out.printf("fixtureindex\t%s\t%s%n", label,
                Files.exists(index) ? RecordTransformDump.base64(index) : "absent");
    }

    static void truth(final Path dir, final String label, final String[] lines) throws Exception {
        final Path file = dir.resolve(label + ".pileup");
        Files.writeString(file, String.join("\n", lines) + "\n", StandardCharsets.UTF_8);
        // A feature file must support random access, so it is indexed the way the tool's own
        // message tells the user to. Without this every run dies before it reads a locus.
        new org.broadinstitute.hellbender.tools.IndexFeatureFile()
                .instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("truth\t%s\t%s%n", label,
                ReferenceQueryDump.escape(Files.readString(file)));
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1.length()))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("s1");
        header.addReadGroup(group);
        return header;
    }

    /** name, start, cigar, bases, flags, quality characters. */
    static SAMRecord read(final SAMFileHeader header, final String[] spec) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(spec[0]);
        record.setFlags(Integer.parseInt(spec[4]));
        record.setReferenceName("chr1");
        record.setAlignmentStart(Integer.parseInt(spec[1]));
        record.setCigarString(spec[2]);
        record.setReadBases(spec[3].getBytes(StandardCharsets.UTF_8));
        final byte[] quals = new byte[spec[5].length()];
        for (int i = 0; i < quals.length; i++) {
            quals[i] = (byte) (spec[5].charAt(i) - 33);
        }
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        record.setAttribute("RG", "rg1");
        if ((Integer.parseInt(spec[4]) & 0x1) != 0) {
            record.setMateReferenceName("chr1");
            record.setMateAlignmentStart(Integer.parseInt(spec[1]));
            record.setInferredInsertSize(4);
        }
        return record;
    }

    /** One run of the tool, with its file and the summary it returned. */
    static void run(final Path dir, final Path fasta, final String label, final String bam,
                    final String pileup, final String[] extra) throws Exception {
        final Path output = dir.resolve("CheckPileup." + label + ".txt");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-I", dir.resolve(bam + ".bam").toString(),
                "-R", fasta.toString(),
                "--pileup", dir.resolve(pileup + ".pileup").toString(),
                "-O", output.toString(),
                "--use-jdk-inflater", "true"));
        argv.addAll(Arrays.asList(extra));

        final Object summary;
        try {
            summary = new CheckPileup().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
            emitReport(dir, label, output);
            return;
        }
        System.out.printf("summary\t%s\t%s%n", label,
                ReferenceQueryDump.escape(String.valueOf(summary)));
        emitReport(dir, label, output);
    }

    /** The file the run left behind, which a failing run writes before it throws. */
    static void emitReport(final Path dir, final String label, final Path output) throws Exception {
        if (Files.exists(output)) {
            System.out.printf("report\t%s\t%s%n", label,
                    ReferenceQueryDump.escape(Files.readString(output)));
        } else {
            System.out.printf("report\t%s\tabsent%n", label);
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + CHR1.length() + "\n",
                StandardCharsets.UTF_8);
        return fasta;
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
