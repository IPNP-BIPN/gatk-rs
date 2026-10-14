/*
 * --flow-filter-alleles-debug-graphs in the callers, end to end (G3.2u, #1355).
 *
 * With --flow-filter-alleles, AlleleFiltering writes, for every linked group of two or more alleles
 * no haplotype carries together, the group's interaction graph to allele.interaction.<contig>.
 * <first start>-<last start>.dot in the working directory: a vertex per allele labelled with its
 * quality against its inverse, and an edge from each allele to every other whose quality drops when
 * the first one's haplotypes are disabled, weighted by the drop. jgrapht's DOTExporter numbers the
 * vertices in the graph's order.
 *
 * Read set: one contig, reads tiled every three bases over four haplotypes in turn: SNVs at 300 and
 * 450, SNVs at 302 and 452, a SNV at 454, and the reference, so 300 and 302 form a group of two and
 * 450, 452 and 454 a group of three. Runs in ONE JVM, in this order; the graphs each run writes are
 * printed and removed before the next.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     graph\t<label>\t<file name>=<the DOT file, escaped>, files sorted by name
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcAlleleGraphsDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.RampedHaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;
import java.util.function.Supplier;

public class HcAlleleGraphsDump {

    static final String CONTIG = "chr1";
    static final int LENGTH = 600;
    static final int READ_LENGTH = 100;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-allele-graphs-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcAlleleGraphsDump: --flow-filter-alleles-debug-graphs in the callers");

        final Random random = new Random(20261014L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        final String reference = ref.toString();
        System.out.printf("fasta\t%s=%s%n", CONTIG, reference);
        final Path fasta = writeReference(dir, reference);

        final int[][] snvs = {{300, 450}, {302, 452}, {454}, {}};
        final String[] haplotypes = new String[snvs.length];
        for (int h = 0; h < snvs.length; h++) {
            final StringBuilder bases = new StringBuilder(reference);
            for (final int pos : snvs[h]) {
                final char base = bases.charAt(pos - 1);
                bases.setCharAt(pos - 1, base == 'A' ? 'C' : base == 'C' ? 'G' : base == 'G' ? 'T' : 'A');
            }
            haplotypes[h] = bases.toString();
        }
        final StringBuilder sam = new StringBuilder("@HD\tVN:1.6\tSO:coordinate\n")
                .append("@SQ\tSN:").append(CONTIG).append("\tLN:").append(LENGTH).append('\n')
                .append("@RG\tID:rg1\tSM:s1\n");
        final String quals = "?".repeat(READ_LENGTH);
        int count = 0;
        for (int start = 1; start + READ_LENGTH - 1 <= LENGTH; start += 3, count++) {
            final String bases = haplotypes[count % haplotypes.length].substring(start - 1, start - 1 + READ_LENGTH);
            sam.append('r').append(count).append('\t').append((count / 4) % 2 == 0 ? 0 : 16).append('\t')
                    .append(CONTIG).append('\t').append(start).append("\t60\t").append(READ_LENGTH).append("M\t*\t0\t0\t")
                    .append(bases).append('\t').append(quals).append("\tRG:Z:rg1\n");
        }
        final Path bam = HcLinkedGraphDump.bam(dir, "input", sam.toString());

        final String filter = "--flow-filter-alleles";
        final String graphs = "--flow-filter-alleles-debug-graphs";
        run(dir, "hc-graphs", HaplotypeCaller::new, bam, fasta, List.of(filter, graphs));
        run(dir, "hc-graphs-qual-300", HaplotypeCaller::new, bam, fasta,
                List.of(filter, graphs, "--flow-filter-alleles-qual-threshold", "300"));
        run(dir, "hc-graphs-unfiltered", HaplotypeCaller::new, bam, fasta, List.of(graphs));
        run(dir, "ramped-graphs", RampedHaplotypeCaller::new, bam, fasta, List.of(filter, graphs));
        run(dir, "mutect2-graphs", Mutect2::new, bam, fasta, List.of(filter, graphs));
    }

    /** The contig in lines of sixty, with its index and dictionary. */
    static Path writeReference(final Path dir, final String reference) throws Exception {
        final StringBuilder fasta = new StringBuilder(">").append(CONTIG).append('\n');
        final int offset = fasta.length();
        for (int i = 0; i < reference.length(); i += 60) {
            fasta.append(reference, i, Math.min(reference.length(), i + 60)).append('\n');
        }
        final Path path = dir.resolve("reference.fasta");
        Files.writeString(path, fasta.toString(), StandardCharsets.US_ASCII);
        Files.writeString(dir.resolve("reference.fasta.fai"),
                CONTIG + "\t" + reference.length() + "\t" + offset + "\t60\t61\n", StandardCharsets.US_ASCII);
        Files.writeString(dir.resolve("reference.dict"),
                "@HD\tVN:1.6\n@SQ\tSN:" + CONTIG + "\tLN:" + reference.length() + "\n", StandardCharsets.US_ASCII);
        return path;
    }

    static void run(final Path dir, final String label, final Supplier<org.broadinstitute.hellbender.cmdline.CommandLineProgram> tool,
                    final Path bam, final Path fasta, final List<String> extra) throws Exception {
        HcAllAnnotationsDump.run(dir, label, tool, bam, fasta, extra);
        final TreeMap<String, String> files = new TreeMap<>();
        try (var listing = Files.list(Path.of("."))) {
            for (final Path path : (Iterable<Path>) listing::iterator) {
                final String name = path.getFileName().toString();
                if (name.startsWith("allele.interaction.") && name.endsWith(".dot")) {
                    files.put(name, Files.readString(path, StandardCharsets.UTF_8));
                    Files.delete(path);
                }
            }
        }
        for (final var file : files.entrySet()) {
            System.out.printf("graph\t%s\t%s=%s%n", label, file.getKey(), ReferenceQueryDump.escape(file.getValue()));
        }
    }
}
