/*
 * The callers with the flow-based likelihood engine, end to end (G3.2f, #1335).
 *
 * --likelihood-calculation-engine FlowBased: FlowBasedAlignmentLikelihoodEngine over Ultima reads.
 * Each sample's reads become flow reads under the flow order and maximal hmer of its first read's
 * read group, each haplotype a key in that order; a read is trimmed in bases to the first
 * haplotype's span (invalid below ten bases), the haplotype clipped to the flows the read covers
 * on it widened by four, and the read's matrix read along it from the first flow reading its first
 * base, stepping four flows, the best sum kept. The matrix is normalized and the poorly modelled
 * reads dropped against the engine's own threshold: three errors at the expected rate per base and
 * two at the filling value. Under a flow engine the reference model skips low-quality deletions
 * and gives a deletion the quality of the base after it when that base is the reference's.
 *
 * One contig of 1400 random bases with a het SNP, a het one-base deletion in an A homopolymer with
 * a het SNP twenty bases on (so a collapsed haplotype carries a substitution), and a het one-base
 * insertion in a C homopolymer; one sample's reads, 100 bases every three, from
 * the two haplotypes alternately, each with a tp array (mostly zero, some +-1) and qualities that
 * follow it, a t0 string on every other read, under flow order TGCA.
 *
 * Output:
 *
 *     fasta\tchr1=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcFlowDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

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

public class HcFlowDump {

    static final int LENGTH = 1400;
    static final int READ_LENGTH = 100;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-flow-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcFlowDump: the callers with the flow-based likelihood engine");

        final Random random = new Random(20261014L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        // An A homopolymer at 601-606 and a C homopolymer at 901-905, fenced by other bases.
        HcDragenEndToEndDump.splice(ref, 601, "AAAAAA");
        ref.setCharAt(599, 'G');
        ref.setCharAt(606, 'T');
        HcDragenEndToEndDump.splice(ref, 901, "CCCCC");
        ref.setCharAt(899, 'A');
        ref.setCharAt(905, 'G');
        if (ref.charAt(299) == 'G') {
            ref.setCharAt(299, 'C');
        }
        final String reference = ref.toString();
        System.out.printf("fasta\tchr1=%s%n", reference);
        final Path fasta = HcDragenEndToEndDump.writeReference(dir, reference);

        // The alternate haplotype as (reference position or 0 for an inserted base, base).
        final List<int[]> alt = new ArrayList<>();
        for (int pos = 1; pos <= LENGTH; pos++) {
            if (pos == 606) {
                continue; // the last A of the homopolymer deleted
            }
            final char base = pos == 300 ? 'G' : pos == 625 ? HcEndToEndDump.other(reference.charAt(pos - 1)) : reference.charAt(pos - 1);
            alt.add(new int[]{pos, base});
            if (pos == 905) {
                alt.add(new int[]{0, 'C'}); // one more C
            }
        }
        final List<int[]> ref0 = new ArrayList<>();
        for (int pos = 1; pos <= LENGTH; pos++) {
            ref0.add(new int[]{pos, reference.charAt(pos - 1)});
        }

        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:").append(LENGTH).append("\n");
        sam.append("@RG\tID:rgF\tSM:f1\tPL:ULTIMA\tFO:TGCATGCATGCA\n");
        final List<String[]> rows = new ArrayList<>();
        int count = 0;
        for (int start = 1; start + READ_LENGTH + 2 <= LENGTH; start += 3, count++) {
            final List<int[]> hap = count % 2 == 0 ? ref0 : alt;
            int index = 0;
            while (hap.get(index)[0] < start) {
                index++;
            }
            if (hap.get(index)[0] == 0) {
                index++;
            }
            final StringBuilder bases = new StringBuilder();
            final StringBuilder cigar = new StringBuilder();
            char op = 0;
            int run = 0;
            int lastRef = hap.get(index)[0] - 1;
            for (int i = index; i < index + READ_LENGTH; i++) {
                final int[] entry = hap.get(i);
                bases.append((char) entry[1]);
                if (entry[0] != 0 && entry[0] > lastRef + 1) {
                    cigar.append(run).append(op);
                    cigar.append(entry[0] - lastRef - 1).append('D');
                    run = 0;
                    op = 'M';
                }
                final char next = entry[0] == 0 ? 'I' : 'M';
                if (next != op && run > 0) {
                    cigar.append(run).append(op);
                    run = 0;
                }
                op = next;
                run++;
                if (entry[0] != 0) {
                    lastRef = entry[0];
                }
            }
            cigar.append(run).append(op);
            final StringBuilder quals = new StringBuilder();
            final StringBuilder tp = new StringBuilder("tp:B:c");
            final StringBuilder t0 = new StringBuilder();
            for (int i = 0; i < READ_LENGTH; i++) {
                final int roll = random.nextInt(10);
                final int t = roll < 7 ? 0 : roll < 9 ? 1 : -1;
                tp.append(',').append(t);
                quals.append((char) (33 + (t == 0 ? 30 + random.nextInt(11) : 8 + random.nextInt(15))));
                t0.append((char) ('!' + 10 + random.nextInt(31)));
            }
            final int flag = (count / 2) % 2 == 0 ? 0 : 16;
            String row = "r" + count + "\t" + flag + "\tchr1\t" + hap.get(index)[0] + "\t60\t" + cigar
                    + "\t*\t0\t0\t" + bases + "\t" + quals + "\tRG:Z:rgF\t" + tp;
            if (count % 2 == 0) {
                row += "\tt0:Z:" + t0;
            }
            rows.add(new String[]{String.valueOf(hap.get(index)[0]), row});
        }
        rows.sort((a, b) -> Integer.compare(Integer.parseInt(a[0]), Integer.parseInt(b[0])));
        for (final String[] row : rows) {
            sam.append(row[1]).append("\n");
        }
        final Path samPath = dir.resolve("input.sam");
        Files.writeString(samPath, sam.toString(), StandardCharsets.UTF_8);
        final Path bam = dir.resolve("input.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(samPath);
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam.toString()));

        final List<String> flow = List.of("--likelihood-calculation-engine", "FlowBased");
        HcAllAnnotationsDump.run(dir, "hc-flow", HaplotypeCaller::new, bam, fasta, flow);
        HcAllAnnotationsDump.run(dir, "hc-flow-optimized", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("--flow-likelihood-optimized-comp")));
        HcAllAnnotationsDump.run(dir, "hc-flow-exact", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("--exact-matching", "--trim-to-haplotype", "false")));
        HcAllAnnotationsDump.run(dir, "hc-flow-untrimmed", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("--trim-to-haplotype", "false")));
        HcAllAnnotationsDump.run(dir, "hc-flow-matrix", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("--flow-use-t0-tag", "--flow-fill-empty-bins-value", "0.0001",
                        "--keep-boundary-flows", "--flow-symmetric-indel-probs", "--flow-lump-probs")));
        HcAllAnnotationsDump.run(dir, "hc-flow-dynamic", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("--enable-dynamic-read-disqualification-for-genotyping",
                        "--dynamic-read-disqualification-threshold", "10")));
        HcAllAnnotationsDump.run(dir, "hc-flow-gvcf", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of("-ERC", "GVCF", "-L", "chr1:550-700")));
        HcAllAnnotationsDump.run(dir, "mutect2-flow", Mutect2::new, bam, fasta, flow);
        // --flow-assembly-collapse-hmer-size: the A homopolymer of six and the C one of five are
        // longer than these thresholds, and the haplotypes' deletions beside them are put back.
        final String collapse = "--flow-assembly-collapse-hmer-size";
        HcAllAnnotationsDump.run(dir, "hc-collapse-4", HaplotypeCaller::new, bam, fasta, List.of(collapse, "4"));
        HcAllAnnotationsDump.run(dir, "hc-collapse-3-partial", HaplotypeCaller::new, bam, fasta,
                List.of(collapse, "3", "--flow-assembly-collapse-partial-mode"));
        HcAllAnnotationsDump.run(dir, "hc-collapse-determined", HaplotypeCaller::new, bam, fasta, List.of(collapse, "-1"));
        HcAllAnnotationsDump.run(dir, "hc-flow-collapse-4", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(flow, List.of(collapse, "4")));
        HcAllAnnotationsDump.run(dir, "ramped-collapse-4", org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller::new,
                bam, fasta, List.of(collapse, "4"));
        HcAllAnnotationsDump.run(dir, "mutect2-collapse-4", Mutect2::new, bam, fasta, List.of(collapse, "4"));
        // --flow-mode: its presets, set where the command line did not set them.
        HcAllAnnotationsDump.run(dir, "hc-flow-mode-standard", HaplotypeCaller::new, bam, fasta, List.of("--flow-mode", "STANDARD"));
        HcAllAnnotationsDump.run(dir, "hc-flow-mode-advanced", HaplotypeCaller::new, bam, fasta, List.of("--flow-mode", "ADVANCED"));
        HcAllAnnotationsDump.run(dir, "hc-flow-mode-standard-sor-1", HaplotypeCaller::new, bam, fasta,
                List.of("--flow-mode", "STANDARD", "--flow-filter-alleles-sor-threshold", "1"));
        HcAllAnnotationsDump.run(dir, "ramped-flow-mode-advanced", org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller::new,
                bam, fasta, List.of("--flow-mode", "ADVANCED"));
        HcAllAnnotationsDump.run(dir, "mutect2-flow-mode-standard", Mutect2::new, bam, fasta, List.of("--flow-mode", "STANDARD"));
        HcAllAnnotationsDump.run(dir, "mutect2-flow-mode-advanced", Mutect2::new, bam, fasta, List.of("--flow-mode", "ADVANCED"));
        // --likelihood-calculation-engine FlowBasedHMM: the PairHMM over flows.
        final List<String> hmm = List.of("--likelihood-calculation-engine", "FlowBasedHMM");
        HcAllAnnotationsDump.run(dir, "hc-flow-hmm", HaplotypeCaller::new, bam, fasta, hmm);
        HcAllAnnotationsDump.run(dir, "hc-flow-hmm-penalties", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(hmm, List.of("--flow-hmm-engine-min-indel-adjust", "20",
                        "--flow-hmm-engine-flat-insertion-penatly", "30", "--flow-hmm-engine-flat-deletion-penatly", "35",
                        "--pcr-indel-model", "AGGRESSIVE")));
        HcAllAnnotationsDump.run(dir, "hc-flow-hmm-no-pcr", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(hmm, List.of("--pcr-indel-model", "NONE", "--enable-dynamic-read-disqualification-for-genotyping")));
        HcAllAnnotationsDump.run(dir, "hc-flow-hmm-cycle-8", HaplotypeCaller::new, bam, fasta,
                HcPedigreeDump.concat(hmm, List.of("--flow-order-cycle-length", "8")));
        HcAllAnnotationsDump.run(dir, "mutect2-flow-hmm", Mutect2::new, bam, fasta, hmm);
    }
}
