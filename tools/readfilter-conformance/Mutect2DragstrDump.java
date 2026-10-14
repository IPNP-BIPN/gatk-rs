/*
 * Mutect2 with DRAGstr's PairHMM gap penalties, end to end (G3.2u, #1355).
 *
 * --dragstr-params-path reaches Mutect2 through createLikelihoodCalculationEngine alone: the
 * PairHMM's gap open and continuation penalties at each read's short tandem repeats, unless
 * --dont-use-dragstr-pair-hmm-scores, in which case the table is not even read.
 * --dragstr-het-hom-ratio feeds HaplotypeCaller's genotype priors and changes nothing here.
 *
 * Read set: HcDragenEndToEndDump's (an A homopolymer at 401 and a CA repeat at 701), as
 * Mutect2MitochondriaDump builds it, and DragstrParams.DEFAULT printed to a file.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text of the one-sample BAM, escaped>
 *     sam\tboth=<the SAM text of the two-sample BAM, escaped>
 *     dragstr\tparams=<the DRAGstr table file, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     error\t<label>\t<exception class>:<message, the dump's directory as <dir>>
 *
 * Usage: Mutect2DragstrDump
 */

import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Random;

public class Mutect2DragstrDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-dragstr-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2DragstrDump: Mutect2 with DRAGstr's PairHMM gap penalties");

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
        System.out.printf("fasta\t%s=%s%n", HcDragenEndToEndDump.CONTIG, reference);
        final Path fasta = HcDragenEndToEndDump.writeReference(dir, reference);

        final String both = HcDragenEndToEndDump.sam(reference);
        final String one = java.util.Arrays.stream(both.split("\n"))
                .filter(line -> !line.contains("rg2"))
                .collect(java.util.stream.Collectors.joining("\n")) + "\n";
        final Path twoSamples = Mutect2ResourcesDump.bam(dir, "both", both);
        final Path bam = Mutect2ResourcesDump.bam(dir, "input", one);
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(one));
        System.out.printf("sam\tboth=%s%n", ReferenceQueryDump.escape(both));

        final Path params = dir.resolve("dragstr.txt");
        DragstrParamUtils.print(DragstrParams.DEFAULT, new GATKPath(params.toString()), "origin", "dump");
        System.out.printf("dragstr\tparams=%s%n", ReferenceQueryDump.escape(Files.readString(params)));
        final String path = params.toString();
        final String missing = dir.resolve("absent.txt").toString();

        run(dir, "plain", bam, fasta, List.of());
        run(dir, "dragstr", bam, fasta, List.of("--dragstr-params-path", path));
        run(dir, "dragstr-ratio", bam, fasta, List.of("--dragstr-params-path", path, "--dragstr-het-hom-ratio", "3"));
        run(dir, "dragstr-off", bam, fasta, List.of("--dragstr-params-path", path, "--dont-use-dragstr-pair-hmm-scores"));
        run(dir, "dragstr-off-missing", bam, fasta, List.of("--dragstr-params-path", missing, "--dont-use-dragstr-pair-hmm-scores"));
        run(dir, "dragstr-tumor-normal", twoSamples, fasta,
                List.of("--dragstr-params-path", path, "--normal-sample", "s2"));
        run(dir, "dragstr-missing", bam, fasta, List.of("--dragstr-params-path", missing));
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta, final List<String> extra)
            throws Exception {
        Mutect2MitochondriaDump.run(dir, label, Mutect2::new, bam, fasta, extra);
    }
}
