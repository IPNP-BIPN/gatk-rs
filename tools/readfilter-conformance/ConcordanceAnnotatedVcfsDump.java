/*
 * Concordance's three annotated VCFs, taken from the reference.
 *
 * The tool's optional record outputs: `-tpfn`, `-tpfp` and `-ftnfn`, each record written out with a
 * STATUS attribute naming the state it was in. Five behaviours this is built to catch.
 *
 *   - THREE OF THE FIVE STATES WRITE TO TWO FILES, AND NOT THE SAME RECORD. A true positive puts
 *     the TRUTH record in -tpfn and the EVAL record in -tpfp, and a filtered false negative puts the
 *     TRUTH record in -tpfn and the EVAL record in -ftnfn. The two files therefore disagree on
 *     everything the two records disagree on while agreeing on the STATUS;
 *   - -tpfn IS WRITTEN AGAINST THE TRUTH HEADER and the other two against the eval header, so the
 *     sample column and the INFO declarations of -tpfn are the truth file's;
 *   - A FILTERED FALSE NEGATIVE IS LABELLED FFN IN BOTH FILES IT REACHES, not FN in one of them,
 *     even though -tpfn is documented as true positives and false negatives;
 *   - THE STATUS LINE AND THE DEFAULT TOOL LINES ARE ADDED IN OPPOSITE ORDERS for -tpfp and -ftnfn,
 *     to the SAME VCFHeader object:
 *
 *         defaultToolHeaderLines.forEach(evalHeader::addMetaDataLine);
 *         evalHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);      // -tpfp
 *         ...
 *         evalHeader.addMetaDataLine(TRUTH_STATUS_HEADER_LINE);
 *         defaultToolHeaderLines.forEach(evalHeader::addMetaDataLine);   // -ftnfn
 *
 *     and the measurement is that IT MAKES NO DIFFERENCE: the writer emits
 *     `getMetaDataInSortedOrder`, so the two insertion orders produce the same lines and the two
 *     runs asking for one output each differ only in their command line. The three runs here are
 *     what says so, rather than a reading of the code;
 *   - AND THE DEFAULT TOOL LINES DO REACH THE FILE, unlike the sibling tool
 *     AnnotateVcfWithExpectedAlleleFraction, which builds its header before adding them: `##source=`
 *     and `##GATKCommandLine` are in all three outputs;
 *   - AND THE ANNOTATION IS A BUILDER COPY, `new VariantContextBuilder(vc).attribute(...)`, so
 *     everything else the record carried survives into the output, filters included: the -ftnfn file
 *     holds records whose FILTER column is not PASS.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     vcfline\t<run>-<which>\t<one line of the output vcf, escaped>
 *     commandline\t<run>-<which>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ConcordanceAnnotatedVcfsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.Concordance;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ConcordanceAnnotatedVcfsDump {

    /** The truth header: its own sample, its own INFO key, and no FILTER line it does not need. */
    static final String TRUTH_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=TRUTHONLY,Number=1,Type=Integer,Description=\"declared by the truth file alone\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruthsample\n";

    /** The eval header: a different sample and a different INFO key. */
    static final String EVAL_HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##INFO=<ID=EVALONLY,Number=1,Type=Integer,Description=\"declared by the eval file alone\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tevalsample\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("concordance-annotated-vcfs-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConcordanceAnnotatedVcfsDump: five states, three files, two records each");

        // One of each state that writes anything, and every pair of records deliberately different
        // in its ID, its QUAL and its INFO so that which side was written is visible.
        final Path truth = writeVcf(dir, "truth", TRUTH_HEADER,
                "chr1\t100\ttruth_tp\tA\tC\t11\tPASS\tTRUTHONLY=1\tGT\t0/1",
                "chr1\t200\ttruth_ffn\tA\tC\t22\tPASS\tTRUTHONLY=2\tGT\t0/1",
                "chr1\t300\ttruth_fn\tA\tC\t33\tPASS\tTRUTHONLY=3\tGT\t0/1");
        final Path eval = writeVcf(dir, "eval", EVAL_HEADER,
                "chr1\t100\teval_tp\tA\tC\t44\tPASS\tEVALONLY=11\tGT\t0/1",
                "chr1\t200\teval_ffn\tA\tC\t55\tweak\tEVALONLY=22\tGT\t0/1",
                "chr1\t400\teval_fp\tA\tC\t66\tPASS\tEVALONLY=44\tGT\t0/1",
                "chr1\t500\teval_ftn\tA\tC\t77\tweak\tEVALONLY=55\tGT\t0/1");

        run(dir, "all-three", truth, eval, true, true, true);
        run(dir, "tpfp-alone", truth, eval, false, true, false);
        run(dir, "ftnfn-alone", truth, eval, false, false, true);
    }

    static Path writeVcf(final Path dir, final String label, final String header,
                         final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path truth, final Path eval,
                    final boolean tpfn, final boolean tpfp, final boolean ftnfn) {
        final Path summary = dir.resolve(label + ".summary");
        final Path tpfnFile = dir.resolve(label + "-tpfn.vcf");
        final Path tpfpFile = dir.resolve(label + "-tpfp.vcf");
        final Path ftnfnFile = dir.resolve(label + "-ftnfn.vcf");
        final List<String> all = new ArrayList<>(List.of(
                "--truth", truth.toString(),
                "--evaluation", eval.toString(),
                "--summary", summary.toString()));
        if (tpfn) {
            all.add("-tpfn");
            all.add(tpfnFile.toString());
        }
        if (tpfp) {
            all.add("-tpfp");
            all.add(tpfpFile.toString());
        }
        if (ftnfn) {
            all.add("-ftnfn");
            all.add(ftnfnFile.toString());
        }
        try {
            new Concordance().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        if (tpfn) {
            print(label + "-tpfn", tpfnFile);
        }
        if (tpfp) {
            print(label + "-tpfp", tpfpFile);
        }
        if (ftnfn) {
            print(label + "-ftnfn", ftnfnFile);
        }
    }

    static void print(final String label, final Path output) {
        final List<String> lines;
        try {
            lines = Files.readAllLines(output, StandardCharsets.UTF_8);
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
            return;
        }
        for (final String line : lines) {
            if (line.startsWith("##GATKCommandLine")) {
                System.out.printf("commandline\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(line.replaceAll("Date=\"[^\"]*\"", "Date=\"MASKED\"")));
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
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
