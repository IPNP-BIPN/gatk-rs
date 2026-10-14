/*
 * Mutect2's Permutect datasets, end to end (G3.2t, #1354).
 *
 * --permutect-training-dataset and --permutect-test-dataset, each with its truth VCF, and
 * --permutect-ref-downsample, --permutect-alt-downsample, --permutect-non-artifact-ratio and
 * --permutect-dataset-mode: every annotated call adds one datum per alternate its engine labels,
 * and closing the engine writes contigs.table and read-groups.table beside the dataset. A training
 * dataset makes every candidate locus active.
 *
 * A read set above its downsampling count is cut by Collections.shuffle without a generator, so the
 * reference's output is then not reproducible; the runs here downsample to counts no read set
 * reaches. ILLUMINA mode reads each read's mate, which an unpaired read refuses.
 *
 * Read set: HcDragenEndToEndDump's, as Mutect2MitochondriaDump builds it, unpaired as it is, paired
 * (consecutive reads as mates), and with Ultima's si and rq tags, paired and not.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\t<name>=<a SAM text, escaped>
 *     vcf\t<name>=<a truth VCF, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     file\t<label>\t<name>=<a dataset or table, escaped>
 *     error\t<label>\t<exception class>:<message, the dump's directory as <dir>>
 *
 * Usage: Mutect2PermutectDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Mutect2PermutectDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-permutect-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2PermutectDump: Mutect2's Permutect datasets");

        // HcDragenEndToEndDump's contig, built as Mutect2MitochondriaDump builds it.
        final Random random = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        HcDragenEndToEndDump.splice(ref, 401, "AAAAAAAAAA");
        HcDragenEndToEndDump.splice(ref, 701, "CACACACACACACA");
        HcDragenEndToEndDump.splice(ref, 1001, "TTTTTTTT");
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String reference = ref.toString();
        final String contig = HcDragenEndToEndDump.CONTIG;
        System.out.printf("fasta\t%s=%s%n", contig, reference);
        final Path fasta = HcDragenEndToEndDump.writeReference(dir, reference);

        final String both = HcDragenEndToEndDump.sam(reference);
        final String one = java.util.Arrays.stream(both.split("\n"))
                .filter(line -> !line.contains("rg2"))
                .collect(java.util.stream.Collectors.joining("\n")) + "\n";
        final String pairedOne = paired(one, false);
        final String pairedBoth = paired(both, false);
        final String ultima = paired(one, true);
        final String ultimaUnpaired = tagged(one);
        final Path unpairedBam = bam(dir, "unpaired", one);
        final Path pairedBam = bam(dir, "paired", pairedOne);
        final Path pairedBothBam = bam(dir, "paired-both", pairedBoth);
        final Path ultimaBam = bam(dir, "ultima", ultima);
        final Path ultimaUnpairedBam = bam(dir, "ultima-unpaired", ultimaUnpaired);

        // The calls of a first plain run over the paired reads, for the truth VCFs.
        final Path plain = dir.resolve("plain.vcf");
        new Mutect2().instanceMain(new String[] {"-I", pairedBam.toString(), "-O", plain.toString(), "-R", fasta.toString()});
        final List<String[]> calls = new ArrayList<>();
        for (final String line : Files.readAllLines(plain)) {
            if (!line.startsWith("#")) {
                calls.add(line.split("\t"));
            }
        }
        final StringBuilder truth = new StringBuilder();
        for (int i = 0; i < calls.size(); i++) {
            final String[] call = calls.get(i);
            if (i % 2 == 0) {
                truth.append(contig).append('\t').append(call[1]).append("\t.\t").append(call[3]).append('\t')
                        .append(call[4].split(",")[0]).append(i == 2 ? "\t.\tLowQual\t.\n" : "\t.\tPASS\t.\n");
            }
        }
        final String truthVcf = Mutect2ResourcesDump.vcf(dir, "truth", truth.toString()).toString();

        final List<String> big = List.of("--permutect-ref-downsample", "1000", "--permutect-alt-downsample", "1000");
        run(dir, "test", pairedBam, fasta, List.of("test"), big, List.of());
        run(dir, "test-truth", pairedBam, fasta, List.of("test"), big, List.of("--permutect-test-truth", truthVcf));
        run(dir, "training", pairedBam, fasta, List.of("training"), big, List.of());
        run(dir, "training-truth", pairedBam, fasta, List.of("training"), big,
                List.of("--permutect-training-truth", truthVcf));
        run(dir, "training-ratio", pairedBam, fasta, List.of("training"), big,
                List.of("--permutect-training-truth", truthVcf, "--permutect-non-artifact-ratio", "3"));
        run(dir, "both-engines", pairedBam, fasta, List.of("training", "test"), big,
                List.of("--permutect-test-truth", truthVcf));
        run(dir, "tumor-normal", pairedBothBam, fasta, List.of("training", "test"), big,
                List.of("--normal-sample", "s2", "--permutect-training-truth", truthVcf));
        run(dir, "two-tumors", pairedBothBam, fasta, List.of("test"), big, List.of());
        run(dir, "ultima", ultimaBam, fasta, List.of("test"), big, List.of("--permutect-dataset-mode", "ULTIMA"));
        run(dir, "ultima-unpaired", ultimaUnpairedBam, fasta, List.of("test"), big,
                List.of("--permutect-dataset-mode", "ULTIMA"));
        run(dir, "unpaired", unpairedBam, fasta, List.of("test"), big, List.of());
    }

    /**
     * The SAM with consecutive reads as mates, both mapped, the first of the two as first of pair;
     * with Ultima, each read gains an si and an rq tag.
     */
    static String paired(final String sam, final boolean ultima) {
        final StringBuilder out = new StringBuilder();
        final List<String[]> reads = new ArrayList<>();
        for (final String line : sam.split("\n")) {
            if (line.startsWith("@")) {
                out.append(line).append('\n');
            } else if (!line.isEmpty()) {
                reads.add(line.split("\t"));
            }
        }
        for (int i = 0; i + 1 < reads.size(); i += 2) {
            final String[] a = reads.get(i);
            final String[] b = reads.get(i + 1);
            final int aStart = Integer.parseInt(a[3]);
            final int bStart = Integer.parseInt(b[3]);
            final int bEnd = bStart + b[9].length() - 1;
            final int aEnd = aStart + a[9].length() - 1;
            final int insert = Math.max(aEnd, bEnd) - Math.min(aStart, bStart) + 1;
            final boolean aReverse = (Integer.parseInt(a[1]) & 16) != 0;
            final boolean bReverse = (Integer.parseInt(b[1]) & 16) != 0;
            a[1] = String.valueOf(1 | 2 | 64 | (aReverse ? 16 : 0) | (bReverse ? 32 : 0));
            b[1] = String.valueOf(1 | 2 | 128 | (bReverse ? 16 : 0) | (aReverse ? 32 : 0));
            b[0] = a[0];
            a[6] = "=";
            b[6] = "=";
            a[7] = String.valueOf(bStart);
            b[7] = String.valueOf(aStart);
            a[8] = String.valueOf(insert);
            b[8] = String.valueOf(-insert);
            for (final String[] read : new String[][] {a, b}) {
                out.append(String.join("\t", read));
                if (ultima) {
                    out.append("\tsi:i:").append(100 + i).append("\trq:f:0.").append(9 - (i % 7));
                }
                out.append('\n');
            }
        }
        return out.toString();
    }

    /** The SAM as it is, each read with Ultima's si and rq tags. */
    static String tagged(final String sam) {
        final StringBuilder out = new StringBuilder();
        int i = 0;
        for (final String line : sam.split("\n")) {
            if (line.isEmpty()) {
                continue;
            }
            out.append(line);
            if (!line.startsWith("@")) {
                out.append("\tsi:i:").append(100 + i).append("\trq:f:0.").append(9 - (i % 7));
                i++;
            }
            out.append('\n');
        }
        return out.toString();
    }

    static Path bam(final Path dir, final String name, final String sam) throws Exception {
        System.out.printf("sam\t%s=%s%n", name, ReferenceQueryDump.escape(sam));
        return Mutect2ResourcesDump.bam(dir, name, sam);
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta, final List<String> engines,
                    final List<String> downsampling, final List<String> extra) throws Exception {
        final Path caseDir = dir.resolve(label);
        Files.createDirectories(caseDir);
        final List<String> argv = new ArrayList<>();
        for (final String engine : engines) {
            argv.add("--permutect-" + engine + "-dataset");
            argv.add(caseDir.resolve(engine + ".dataset").toString());
        }
        argv.addAll(downsampling);
        argv.addAll(extra);
        Mutect2MitochondriaDump.run(dir, label, Mutect2::new, bam, fasta, argv);
        for (final String name : List.of("training.dataset", "test.dataset", "contigs.table", "read-groups.table")) {
            final Path file = caseDir.resolve(name);
            if (Files.exists(file)) {
                System.out.printf("file\t%s\t%s=%s%n", label, name,
                        ReferenceQueryDump.escape(Files.readString(file, StandardCharsets.UTF_8)));
            }
        }
    }
}
