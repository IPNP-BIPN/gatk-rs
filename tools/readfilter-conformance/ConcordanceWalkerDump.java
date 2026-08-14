/*
 * AbstractConcordanceWalker's ConcordanceIterator, taken from the reference.
 *
 * The engine slice under the two remaining validation tools: it walks a truth VCF and an eval VCF in
 * lockstep and labels each step with one of five states. The labels are not the interesting part;
 * which iterator each step consumes is.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A SAME-LOCUS DISAGREEMENT ADVANCES TRUTH ALONE. `return TruthVersusEval.falseNegative(
 *     truthIterator.next())` consumes truth and leaves the eval record in place, so that record is
 *     compared against the NEXT truth record: one disagreement moves the whole eval side out of
 *     step with the truth side rather than pairing locally;
 *   - AND THE COMMENT SAYS SO: "we could equally well advance eval". The choice decides which
 *     records the tool ever sees;
 *   - A FILTERED EVAL RECORD AT A TRUTH LOCUS CONSUMES BOTH, as a filtered false negative, without
 *     ever asking whether the two agree;
 *   - AN EVAL-ONLY FILTERED RECORD IS A FILTERED TRUE NEGATIVE, so a filtered record is labelled
 *     differently depending on whether truth has anything at that position;
 *   - THE DEFAULT TRUTH FILTER DROPS FILTERED TRUTH RECORDS, so a truth record with a FILTER column
 *     never reaches the iterator, and the eval record at that locus becomes a FALSE POSITIVE;
 *   - THE DEFAULT EVAL FILTER KEEPS EVERYTHING, which is what makes the two filtered states
 *     reachable at all: a walker filtering filtered records on both sides, as
 *     EvaluateInfoFieldConcordance does, can never produce them;
 *   - AND THE ORDER IS THE DICTIONARY'S, VariantContextComparator comparing contig index and then
 *     start, so nothing looks at the end of a record.
 *
 * Output:
 *
 *     input\t<label>\t<the whole vcf, escaped>
 *     state\t<run>\t<index>\t<STATE>\t<truth or ->\t<eval or ->
 *     error\t<run>\t<exception class>:<message>
 *
 * Usage: ConcordanceWalkerDump
 */

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.AbstractConcordanceWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

public class ConcordanceWalkerDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FILTER=<ID=weak,Description=\"Was already there\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "##contig=<ID=chr2,length=900>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    /** The run being printed, since the walkers are constructed by the command line. */
    static String run = "";

    /** The step counter, reset per run. */
    static int step = 0;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("concordancewalker-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ConcordanceWalkerDump: five states, and which iterator each one consumes");

        // Truth and eval built so that every state is reachable and the cascade is visible.
        //
        //   100  both agree                      -> true positive
        //   110  truth only                      -> false negative
        //   120  eval only                       -> false positive
        //   130  eval only, filtered             -> filtered true negative
        //   140  both, eval filtered             -> filtered false negative
        //   200  both, different reference base  -> disagreement: truth advances alone
        //   210  truth only, and it inherits the eval record left behind at 200
        //   300  truth filtered, eval present    -> false positive, the truth record never arriving
        final Path truth = writeVcf(dir, "truth",
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t110\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t140\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tAT\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t210\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tweak\t.\tGT\t0/1");
        final Path eval = writeVcf(dir, "eval",
                "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t120\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t130\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t140\t.\tA\tC\t50\tweak\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");

        // A second contig, to show the order is the dictionary's rather than the position's.
        final Path acrossContigs = writeVcf(dir, "across-contigs",
                "chr1\t500\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path secondContigOnly = writeVcf(dir, "second-contig-only",
                "chr2\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");

        run("default-filters", new DefaultFilterProbe(), truth, eval);
        run("filtered-dropped-both-sides", new DroppingProbe(), truth, eval);
        run("across-contigs", new DefaultFilterProbe(), acrossContigs, secondContigOnly);
    }

    static void run(final String label, final AbstractConcordanceWalker walker,
                    final Path truth, final Path eval) {
        run = label;
        step = 0;
        try {
            walker.instanceMain(new String[] {
                "--truth", truth.toString(), "--evaluation", eval.toString()});
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** `<contig>:<start> <ref>/<alts> <filters>`, or `-` when the state carries no such record. */
    static String render(final VariantContext variant) {
        if (variant == null) {
            return "-";
        }
        final StringBuilder alternates = new StringBuilder();
        for (int i = 0; i < variant.getAlternateAlleles().size(); i++) {
            if (i > 0) {
                alternates.append(',');
            }
            alternates.append(variant.getAlternateAllele(i).getDisplayString());
        }
        return String.format("%s:%d %s/%s %s", variant.getContig(), variant.getStart(),
                variant.getReference().getDisplayString(), alternates,
                variant.isFiltered() ? String.join(";", variant.getFilters()) : "PASS");
    }

    /** The base class's own filters: filtered truth records dropped, every eval record kept. */
    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class DefaultFilterProbe extends AbstractConcordanceWalker {
        @Override
        protected void apply(final TruthVersusEval truthVersusEval, final ReadsContext readsContext,
                             final ReferenceContext refContext) {
            System.out.printf("state\t%s\t%d\t%s\t%s\t%s%n", run, step++,
                    truthVersusEval.getConcordance().name(),
                    render(truthVersusEval.hasTruth() ? truthVersusEval.getTruth() : null),
                    render(truthVersusEval.hasEval() ? truthVersusEval.getEval() : null));
        }

        @Override
        protected boolean areVariantsAtSameLocusConcordant(final VariantContext truth,
                                                           final VariantContext eval) {
            return truth.getReference().equals(eval.getReference())
                    && eval.getAlternateAlleles().contains(truth.getAlternateAllele(0));
        }
    }

    /** What EvaluateInfoFieldConcordance does: filtered records dropped on both sides. */
    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class DroppingProbe extends DefaultFilterProbe {
        @Override
        protected org.apache.commons.collections4.Predicate<VariantContext> makeTruthVariantFilter() {
            return vc -> !vc.isFiltered() && !vc.isSymbolicOrSV();
        }

        @Override
        protected org.apache.commons.collections4.Predicate<VariantContext> makeEvalVariantFilter() {
            return vc -> !vc.isFiltered() && !vc.isSymbolicOrSV();
        }
    }

    static Path writeVcf(final Path dir, final String label, final String... records)
            throws Exception {
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
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
