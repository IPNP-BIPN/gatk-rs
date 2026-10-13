/*
 * HaplotypeCaller in DRAGEN mode end to end, taken from the reference (G4.1, #1330).
 *
 * Two samples' reads over one random contig with three short tandem repeats spliced in, written as
 * an indexed BAM, called by the tool itself. The data are built for what DRAGEN mode changes:
 *
 *  - a one-base deletion inside a ten-base A homopolymer and a CA insertion inside a seven-unit
 *    CA repeat, where DRAGstr's priors and gap penalties differ from the flat ones;
 *  - a site where s2 is homozygous reference but six reads mapped at quality 6 to 14 carry an
 *    alternate base, which FRD explains as foreign reads;
 *  - a site where s1 is homozygous reference but its reads carry an alternate base only in their
 *    last twenty-five bases, at base quality 22, which BQD explains as a dropout;
 *  - an ordinary het SNP and a hom SNP for contrast;
 *  - an XQ tag on every read (two and a half times its mapping quality, at most 250), which
 *    --transform-dragen-mapping-quality swaps in as the mapping quality.
 *
 * Runs, in this order and in ONE JVM, so the generator carries over as the earlier suite's does:
 * the default, --dragen-mode, --dragen-mode with a DRAGstr table, --apply-bqd, --apply-frd,
 * --apply-frd with a depth cap, --transform-dragen-mapping-quality, the DRAGstr table alone, the
 * table without its priors, the table without its PairHMM scores, the table with a het:hom ratio
 * of 4, --dragen-mode in GVCF mode on s1, --dragen-mode with an explicit calling threshold (the mode
 * does not override what the command line set), and --dragen-mode beside --flow-mode (refused).
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     dragstr\tparams=<the DRAGstr table file, escaped>
 *     out\t<label>=<the output VCF's ##INFO lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcDragenEndToEndDump
 */

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;

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

public class HcDragenEndToEndDump {

    static final int READ_LENGTH = 100;
    static final String CONTIG = "chr1";
    static final int LENGTH = 1400;

    /** A variant: a SNP (deleted 0, inserted null), a deletion after pos, or an insertion after pos. */
    record Variant(int pos, int deleted, String inserted) {}

    // The repeats are spliced in at 401-410 (A x10), 701-714 (CA x7) and 1001-1008 (T x8).
    static final Variant[] VARIANTS = {
            new Variant(300, 0, null),      // 0: het SNP in s1, hom in s2
            new Variant(405, 1, null),      // 1: one A deleted inside the homopolymer
            new Variant(712, 0, "CA"),      // 2: one CA unit inserted inside the repeat
            new Variant(900, 0, null),      // 3: het SNP in s2
    };

    /** Each sample's two haplotypes, as the variants they carry. */
    static final int[][][] HAPLOTYPES = {
            {{0, 1}, {1, 2}},
            {{0, 2}, {0, 3}},
    };

    /** FRD's site: s2's reads at low mapping quality carry an alternate here. */
    static final int FOREIGN_SITE = 560;
    /** BQD's site: s1's reads carry an alternate here only in their last five bases. */
    static final int DROPOUT_SITE = 1180;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-dragen-end-to-end-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcDragenEndToEndDump: HaplotypeCaller in DRAGEN mode end to end");

        final Random random = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        splice(ref, 401, "AAAAAAAAAA");
        splice(ref, 701, "CACACACACACACA");
        splice(ref, 1001, "TTTTTTTT");
        // Keep each repeat bounded by a different base, so its length is what was spliced in.
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String reference = ref.toString();
        System.out.printf("fasta\t%s=%s%n", CONTIG, reference);
        final Path fasta = writeReference(dir, reference);

        final String sam = sam(reference);
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

        final Path params = dir.resolve("dragstr.txt");
        DragstrParamUtils.print(DragstrParams.DEFAULT, new GATKPath(params.toString()), "origin", "dump");
        System.out.printf("dragstr\tparams=%s%n", ReferenceQueryDump.escape(Files.readString(params)));
        final String table = params.toString();

        run(dir, "default", bam, fasta, List.of());
        run(dir, "dragen-mode", bam, fasta, List.of("--dragen-mode"));
        run(dir, "dragen-mode-dragstr", bam, fasta, List.of("--dragen-mode", "--dragstr-params-path", table));
        run(dir, "apply-bqd", bam, fasta, List.of("--apply-bqd"));
        run(dir, "apply-frd", bam, fasta, List.of("--apply-frd"));
        run(dir, "apply-frd-depth-cap", bam, fasta, List.of("--apply-frd", "--max-effective-depth-adjustment-for-frd", "4"));
        run(dir, "transform-mq", bam, fasta, List.of("--transform-dragen-mapping-quality"));
        run(dir, "dragstr", bam, fasta, List.of("--dragstr-params-path", table));
        run(dir, "dragstr-no-priors", bam, fasta, List.of("--dragstr-params-path", table, "--dont-use-dragstr-priors"));
        run(dir, "dragstr-no-hmm-scores", bam, fasta, List.of("--dragstr-params-path", table, "--dont-use-dragstr-pair-hmm-scores"));
        run(dir, "dragstr-het-hom-ratio", bam, fasta, List.of("--dragstr-params-path", table, "--dragstr-het-hom-ratio", "4"));
        run(dir, "dragen-mode-gvcf", bam, fasta, List.of("--dragen-mode", "-ERC", "GVCF", "--sample-name", "s1"));
        run(dir, "dragen-mode-explicit-threshold", bam, fasta, List.of("--dragen-mode", "--standard-min-confidence-threshold-for-calling", "60"));
        run(dir, "dragen-mode-and-flow-mode", bam, fasta, List.of("--dragen-mode", "--flow-mode", "STANDARD"));
    }

    static void splice(final StringBuilder ref, final int position, final String repeat) {
        ref.replace(position - 1, position - 1 + repeat.length(), repeat);
    }

    /** One haplotype as (reference position, base) pairs; position 0 marks an inserted base. */
    static List<int[]> haplotype(final String reference, final int[] carried) {
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
                if (variant.pos() != position) {
                    continue;
                }
                if (variant.inserted() != null) {
                    inserted = variant.inserted();
                } else if (variant.deleted() > 0) {
                    skip = variant.deleted();
                } else {
                    base = HcEndToEndDump.other(base);
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

    static String sam(final String reference) {
        final List<String[]> rows = new ArrayList<>();
        int count = 0;
        for (int s = 0; s < HAPLOTYPES.length; s++) {
            final List<List<int[]>> haps = List.of(haplotype(reference, HAPLOTYPES[s][0]),
                    haplotype(reference, HAPLOTYPES[s][1]));
            for (int start = 1 + s; start + READ_LENGTH + 10 <= LENGTH; start += 4, count++) {
                final String[] r = HcEndToEndDump.read(haps.get(count % 2), start, count);
                final int pos = Integer.parseInt(r[0]);
                final int end = pos + READ_LENGTH - 1;
                final char[] bases = r[1].toCharArray();
                final char[] quals = r[2].toCharArray();
                final boolean plain = !r[3].contains("I") && !r[3].contains("D");
                int mapq = count % 11 == 5 ? 20 : 60;
                // FRD: s2's reads over the foreign site at low mapping quality carry an alternate.
                if (s == 1 && plain && pos <= FOREIGN_SITE && FOREIGN_SITE <= end && count % 3 == 0) {
                    mapq = 6 + (count % 9);
                    bases[FOREIGN_SITE - pos] = HcEndToEndDump.other(bases[FOREIGN_SITE - pos]);
                }
                // BQD: s1's reads ending within twenty-five bases past the dropout site carry an
                // alternate there at base quality 22; those on the reverse strand begin there instead.
                final int flag = (count / 2) % 2 == 0 ? 0 : 16;
                if (s == 0 && plain && pos <= DROPOUT_SITE && DROPOUT_SITE <= end) {
                    final int fromEnd = flag == 0 ? end - DROPOUT_SITE : DROPOUT_SITE - pos;
                    if (fromEnd < 25) {
                        bases[DROPOUT_SITE - pos] = HcEndToEndDump.other(bases[DROPOUT_SITE - pos]);
                        quals[DROPOUT_SITE - pos] = (char) (33 + 22);
                    }
                }
                final int xq = Math.min(250, mapq * 5 / 2);
                rows.add(new String[]{String.valueOf(pos),
                        "r" + count + "\t" + flag + "\t" + CONTIG + "\t" + pos + "\t" + mapq
                                + "\t" + HcEndToEndDump.cigar(r[3]) + "\t*\t0\t0\t" + new String(bases) + "\t" + new String(quals)
                                + "\tRG:Z:rg" + (s + 1) + "\tXQ:i:" + xq});
            }
        }
        rows.sort((a, b) -> Integer.compare(Integer.parseInt(a[0]), Integer.parseInt(b[0])));
        final StringBuilder out = new StringBuilder();
        out.append("@HD\tVN:1.6\tSO:coordinate\n");
        out.append("@SQ\tSN:").append(CONTIG).append("\tLN:").append(LENGTH).append("\n");
        out.append("@RG\tID:rg1\tSM:s1\n");
        out.append("@RG\tID:rg2\tSM:s2\n");
        for (final String[] row : rows) {
            out.append(row[1]).append("\n");
        }
        return out.toString();
    }

    static Path writeReference(final Path dir, final String reference) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder();
        bases.append(">").append(CONTIG).append("\n");
        for (int i = 0; i < LENGTH; i += 60) {
            bases.append(reference, i, Math.min(i + 60, LENGTH)).append("\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(
                List.of(new htsjdk.samtools.SAMSequenceRecord(CONTIG, LENGTH))));
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
                    ReferenceQueryDump.escape(HcEndToEndDump.masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        final StringBuilder body = new StringBuilder();
        boolean inBody = false;
        for (final String line : Files.readAllLines(out)) {
            inBody |= line.startsWith("#CHROM");
            if ((inBody || line.startsWith("##INFO")) && !line.isEmpty()) {
                body.append(line).append("\n");
            }
        }
        System.out.printf("out\t%s=%s%n", label, ReferenceQueryDump.escape(body.toString()));
    }
}
