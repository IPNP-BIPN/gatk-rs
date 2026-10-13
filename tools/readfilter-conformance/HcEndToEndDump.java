/*
 * HaplotypeCaller end to end, taken from the reference.
 *
 * Two samples' reads over two random contigs, written as an indexed BAM, called by the tool
 * itself the way a user runs it. Everything the earlier HaplotypeCaller suites measured piece by
 * piece is chained here: the read filters, the activity profile and its soft-clip spreading, the
 * region boundaries, the assembly, the likelihoods, the genotyping, the default annotations, the
 * physical phasing and the VCF writer.
 *
 * Each sample's reads come from two haplotypes laid over the reference: SNPs het and hom, a
 * two-base deletion, a three-base insertion, two SNPs thirty bases apart in cis for one sample
 * and in trans for the other, and a SNP on the second contig. One read starts every three bases
 * per sample, alternating haplotype and strand; every tenth is mapped at quality 10 and every
 * fifteenth is a duplicate; reads starting over chr1:1290-1330 carry twenty high-quality
 * soft-clipped bases; one read carries an N and one an IUPAC R; and sixty reads start at
 * chr1:250, more than the positional downsampler keeps, so it draws from the generator.
 *
 * Runs, in this order and in ONE JVM: the random generator is never reset between them, so a QD
 * above 35 draws from the same stream the port's process-wide generator replays.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output VCF from its #CHROM line on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcEndToEndDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcEndToEndDump {

    static final int READ_LENGTH = 100;
    static final String[] CONTIGS = {"chr1", "chr2"};
    static final int[] LENGTHS = {1500, 800};

    /** A variant: a SNP (len 0), a deletion of len bases after pos, or an insertion after pos. */
    record Variant(int contig, int pos, int deleted, String inserted) {}

    static final Variant[] VARIANTS = {
            new Variant(0, 300, 0, null),     // 0: SNP
            new Variant(0, 520, 0, null),     // 1: SNP
            new Variant(0, 760, 2, null),     // 2: deletion of 761-762
            new Variant(0, 1000, 0, "TTG"),   // 3: insertion after 1000
            new Variant(0, 1200, 0, null),    // 4: SNP
            new Variant(0, 1230, 0, null),    // 5: SNP
            new Variant(1, 400, 0, null),     // 6: SNP on chr2
    };

    /** Each sample's two haplotypes, as the variants they carry. */
    static final int[][][] HAPLOTYPES = {
            // s1: 300 het, 520 hom, the indels het, 1200 and 1230 het in cis, chr2:400 het.
            {{1}, {0, 1, 2, 3, 4, 5, 6}},
            // s2: 300 hom, 520 het, the deletion het, 1200 and 1230 het in trans.
            {{0, 1, 4}, {0, 2, 5}},
    };

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-end-to-end-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcEndToEndDump: HaplotypeCaller end to end");

        final Random random = new Random(20261010L);
        final String[] references = new String[CONTIGS.length];
        for (int c = 0; c < CONTIGS.length; c++) {
            final StringBuilder ref = new StringBuilder();
            for (int i = 0; i < LENGTHS[c]; i++) {
                ref.append("ACGT".charAt(random.nextInt(4)));
            }
            references[c] = ref.toString();
            System.out.printf("fasta\t%s=%s%n", CONTIGS[c], references[c]);
        }
        final Path fasta = writeReference(dir, references);

        final String sam = sam(references, random);
        final Path samPath = dir.resolve("input.sam");
        Files.writeString(samPath, sam, StandardCharsets.UTF_8);
        final Path bam = dir.resolve("input.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(samPath);
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam));

        run(dir, "default", bam, fasta, List.of());
        run(dir, "interval", bam, fasta, List.of("-L", "chr1:200-600", "-L", "chr2"));
        run(dir, "no-phasing", bam, fasta, List.of("--do-not-run-physical-phasing", "true"));
        run(dir, "call-threshold-500", bam, fasta,
                List.of("--standard-min-confidence-threshold-for-calling", "500"));
        run(dir, "no-downsampling", bam, fasta, List.of("--max-reads-per-alignment-start", "0"));
        run(dir, "low-mapq-kept", bam, fasta,
                List.of("--disable-read-filter", "MappingQualityReadFilter"));
        run(dir, "min-base-quality-30", bam, fasta, List.of("--min-base-quality-score", "30"));
    }

    /** One haplotype as (reference position, base) pairs; position 0 marks an inserted base. */
    static List<int[]> haplotype(final String reference, final int contig, final int[] carried) {
        final List<int[]> elements = new ArrayList<>();
        int skip = 0;
        for (int position = 1; position <= reference.length(); position++) {
            if (skip > 0) {
                skip--;
                continue;
            }
            char base = reference.charAt(position - 1);
            String inserted = null;
            for (final int v : carried) {
                final Variant variant = VARIANTS[v];
                if (variant.contig() != contig || variant.pos() != position) {
                    continue;
                }
                if (variant.inserted() != null) {
                    inserted = variant.inserted();
                } else if (variant.deleted() > 0) {
                    skip = variant.deleted();
                } else {
                    base = other(base);
                }
            }
            elements.add(new int[]{position, base});
            if (inserted != null) {
                for (final char b : inserted.toCharArray()) {
                    elements.add(new int[]{0, b});
                }
            }
        }
        return elements;
    }

    static char other(final char base) {
        return "CGTA".charAt("ACGT".indexOf(base));
    }

    /** A read off a haplotype from its first aligned base at or after start. */
    static String[] read(final List<int[]> hap, final int start, final int count) {
        int first = 0;
        while (hap.get(first)[0] == 0 || hap.get(first)[0] < start) {
            first++;
        }
        int last = Math.min(first + READ_LENGTH, hap.size());
        while (hap.get(last - 1)[0] == 0) {
            last--;
        }
        final StringBuilder bases = new StringBuilder();
        final StringBuilder quals = new StringBuilder();
        final List<String> ops = new ArrayList<>();
        int previous = -1;
        for (int i = first; i < last; i++) {
            final int[] element = hap.get(i);
            bases.append((char) element[1]);
            quals.append((char) (33 + 25 + (count * 7 + i) % 15));
            if (element[0] == 0) {
                ops.add("I");
            } else {
                if (previous != -1 && element[0] > previous + 1) {
                    for (int d = previous + 1; d < element[0]; d++) {
                        ops.add("D");
                    }
                }
                ops.add("M");
                previous = element[0];
            }
        }
        return new String[]{String.valueOf(hap.get(first)[0]), bases.toString(), quals.toString(),
                String.join("", ops)};
    }

    static String sam(final String[] references, final Random random) {
        final List<String[]> rows = new ArrayList<>();
        int count = 0;
        for (int c = 0; c < CONTIGS.length; c++) {
            for (int s = 0; s < HAPLOTYPES.length; s++) {
                final List<List<int[]>> haps = List.of(haplotype(references[c], c, HAPLOTYPES[s][0]),
                        haplotype(references[c], c, HAPLOTYPES[s][1]));
                for (int start = 1 + s; start + READ_LENGTH + 10 <= LENGTHS[c]; start += 3, count++) {
                    final String[] r = read(haps.get(count % 2), start, count);
                    int pos = Integer.parseInt(r[0]);
                    String bases = r[1];
                    String ops = r[3];
                    int flag = (count / 2) % 2 == 0 ? 0 : 16;
                    int mapq = 60;
                    // Every tenth read is mapped at quality 10, every fifteenth a duplicate.
                    if (count % 10 == 7) {
                        mapq = 10;
                    }
                    if (count % 15 == 4) {
                        flag |= 1024;
                    }
                    // High-quality soft clips over chr1:1300-1330: the first twenty bases clipped
                    // and replaced by random ones, so the pileup sees them as clipped evidence.
                    if (c == 0 && s == 0 && pos >= 1290 && pos <= 1330 && ops.startsWith("MMMMMMMMMMMMMMMMMMMM")) {
                        final StringBuilder clipped = new StringBuilder();
                        for (int i = 0; i < 20; i++) {
                            clipped.append("ACGT".charAt(random.nextInt(4)));
                        }
                        bases = clipped + bases.substring(20);
                        ops = "SSSSSSSSSSSSSSSSSSSS" + ops.substring(20);
                        pos += 20;
                    }
                    // One read with an N and one with an IUPAC code.
                    if (count == 40) {
                        bases = bases.substring(0, 50) + "N" + bases.substring(51);
                    }
                    if (count == 41) {
                        bases = bases.substring(0, 50) + "R" + bases.substring(51);
                    }
                    rows.add(new String[]{String.valueOf(c), String.valueOf(pos),
                            "r" + count + "\t" + flag + "\t" + CONTIGS[c] + "\t" + pos + "\t" + mapq
                                    + "\t" + cigar(ops) + "\t*\t0\t0\t" + bases + "\t" + r[2]
                                    + "\tRG:Z:rg" + (s + 1)});
                }
            }
        }
        // Sixty reads of s2 starting at chr1:250, alternately with and without the SNP at 300:
        // more than the positional downsampler keeps, so which ten it drops moves s2's depths.
        final List<List<int[]>> stack = List.of(haplotype(references[0], 0, new int[]{}),
                haplotype(references[0], 0, new int[]{0}));
        for (int i = 0; i < 60; i++, count++) {
            final String[] r = read(stack.get(i % 2), 250, count);
            rows.add(new String[]{"0", r[0], "s" + i + "\t" + (i % 2 == 0 ? 0 : 16) + "\tchr1\t" + r[0]
                    + "\t60\t" + cigar(r[3]) + "\t*\t0\t0\t" + r[1] + "\t" + r[2] + "\tRG:Z:rg2"});
        }
        rows.sort((a, b) -> a[0].equals(b[0]) ? Integer.compare(Integer.parseInt(a[1]), Integer.parseInt(b[1]))
                : a[0].compareTo(b[0]));
        final StringBuilder out = new StringBuilder();
        out.append("@HD\tVN:1.6\tSO:coordinate\n");
        for (int c = 0; c < CONTIGS.length; c++) {
            out.append("@SQ\tSN:").append(CONTIGS[c]).append("\tLN:").append(LENGTHS[c]).append("\n");
        }
        out.append("@RG\tID:rg1\tSM:s1\n");
        out.append("@RG\tID:rg2\tSM:s2\n");
        for (final String[] row : rows) {
            out.append(row[2]).append("\n");
        }
        return out.toString();
    }

    static String cigar(final String ops) {
        final StringBuilder out = new StringBuilder();
        int i = 0;
        while (i < ops.length()) {
            int j = i;
            while (j < ops.length() && ops.charAt(j) == ops.charAt(i)) {
                j++;
            }
            out.append(j - i).append(ops.charAt(i));
            i = j;
        }
        return out.toString();
    }

    static Path writeReference(final Path dir, final String[] references) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder();
        final List<htsjdk.samtools.SAMSequenceRecord> records = new ArrayList<>();
        for (int c = 0; c < CONTIGS.length; c++) {
            bases.append(">").append(CONTIGS[c]).append("\n");
            for (int i = 0; i < LENGTHS[c]; i += 60) {
                bases.append(references[c], i, Math.min(i + 60, LENGTHS[c])).append("\n");
            }
            records.add(new htsjdk.samtools.SAMSequenceRecord(CONTIGS[c], LENGTHS[c]));
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(records));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
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
            new HaplotypeCaller().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        final StringBuilder body = new StringBuilder();
        boolean inBody = false;
        for (final String line : Files.readAllLines(out)) {
            inBody |= line.startsWith("#CHROM");
            if (inBody && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
