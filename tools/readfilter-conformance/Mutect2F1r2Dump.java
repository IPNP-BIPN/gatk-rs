/*
 * Mutect2's F1R2 counts, end to end (G3.2s, #1353).
 *
 * --f1r2-tar-gz with --f1r2-median-mq, --f1r2-min-bq and --f1r2-max-depth: an F1R2CountsCollector
 * over the tumour samples counts the tumour pileup of every locus isActive examines (a forced
 * locus and an empty pileup are not), and onTraversalSuccess writes each sample's alt table, its
 * reference histograms and its depth-one alt histograms into the archive. The name must end in
 * .tar.gz, checked in onTraversalStart.
 *
 * Read set: HcDragenEndToEndDump's, as Mutect2MitochondriaDump builds it.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text of the one-sample BAM, escaped>
 *     sam\tboth=<the SAM text of the two-sample BAM, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     member\t<label>\t<member name>=<the member's text, escaped>, members sorted by name and an
 *         alt histogram's columns sorted
 *     error\t<label>\t<exception class>:<message, the dump's directory as <dir>>
 *
 * Usage: Mutect2F1r2Dump
 */

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

public class Mutect2F1r2Dump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-f1r2-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2F1r2Dump: Mutect2's F1R2 counts");

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
        final String sam = java.util.Arrays.stream(both.split("\n"))
                .filter(line -> !line.contains("rg2"))
                .collect(java.util.stream.Collectors.joining("\n")) + "\n";
        final Path twoSamples = Mutect2ResourcesDump.bam(dir, "both", both);
        final Path bam = Mutect2ResourcesDump.bam(dir, "input", sam);
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam));
        System.out.printf("sam\tboth=%s%n", ReferenceQueryDump.escape(both));

        run(dir, "f1r2", bam, fasta, "f1r2.tar.gz", List.of());
        run(dir, "f1r2-thresholds", bam, fasta, "f1r2-thresholds.tar.gz",
                List.of("--f1r2-median-mq", "30", "--f1r2-min-bq", "10", "--f1r2-max-depth", "20"));
        run(dir, "f1r2-strict", bam, fasta, "f1r2-strict.tar.gz", List.of("--f1r2-median-mq", "61"));
        run(dir, "f1r2-tumor-normal", twoSamples, fasta, "f1r2-tumor-normal.tar.gz", List.of("--normal-sample", "s2"));
        run(dir, "f1r2-both-tumors", twoSamples, fasta, "f1r2-both-tumors.tar.gz", List.of());
        run(dir, "f1r2-gvcf", bam, fasta, "f1r2-gvcf.tar.gz", List.of("-ERC", "GVCF", "-L", HcDragenEndToEndDump.CONTIG + ":380-460"));
        run(dir, "f1r2-interval", bam, fasta, "f1r2-interval.tar.gz", List.of("-L", HcDragenEndToEndDump.CONTIG + ":1-5"));
        run(dir, "f1r2-bad-name", bam, fasta, "f1r2.tgz", List.of());
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta, final String archive,
                    final List<String> extra) throws Exception {
        final Path tarGz = dir.resolve(archive);
        final List<String> argv = new ArrayList<>(List.of("--f1r2-tar-gz", tarGz.toString()));
        argv.addAll(extra);
        Mutect2MitochondriaDump.run(dir, label, Mutect2::new, bam, fasta, argv);
        if (!Files.exists(tarGz)) {
            return;
        }
        final TreeMap<String, String> members = new TreeMap<>();
        try (InputStream in = Files.newInputStream(tarGz);
             TarArchiveInputStream tar = new TarArchiveInputStream(new GZIPInputStream(in))) {
            TarArchiveEntry entry;
            while ((entry = tar.getNextTarEntry()) != null) {
                if (!entry.isDirectory()) {
                    members.put(entry.getName(), new String(tar.readAllBytes(), StandardCharsets.UTF_8));
                }
            }
        }
        for (final var member : members.entrySet()) {
            // An alt histogram's columns come in an identity hash order, not stable from run to
            // run, so only the sorted form can be a golden.
            final String content = member.getKey().endsWith(".alt_histogram")
                    ? CollectF1R2CountsDump.sortColumns(member.getValue()) : member.getValue();
            System.out.printf("member\t%s\t%s=%s%n", label, member.getKey(), ReferenceQueryDump.escape(content));
        }
    }
}
