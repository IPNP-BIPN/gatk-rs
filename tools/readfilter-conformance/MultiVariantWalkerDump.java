/*
 * MultiVariantWalker, MultiVariantWalkerGroupedOnStart and MultiVariantWalkerGroupedByOverlap,
 * taken from the reference through three probe walkers written for this measurement.
 *
 * The three classes are the engine under CombineGVCFs, ReblockGVCF, VariantEval, ApplyVQSR,
 * VariantRecalibrator, FilterAlignmentArtifacts, VCFComparator, SVStratify, SVCluster,
 * GroupedSVCluster and JointGermlineCNVSegmentation. Everything those tools share is what happens
 * BEFORE their apply sees a record: which inputs are accepted, whose dictionary wins, what the
 * merged header and sample list are, in which order records from several files arrive, and how
 * the grouped walkers cut that stream into groups. So the probes print exactly that.
 *
 * Output:
 *
 *     fasta\t<the reference FASTA, escaped>
 *     input\t<label>\t<indexed|unindexed>\t<the whole vcf, escaped>
 *     dict\t<run>\t<name:length,...>             the driving variants' merged dictionary
 *     samples\t<run>\t<s1,s2,...>                getSamplesForVariants()
 *     header\t<run>\t<the merged header's lines in sorted order, escaped>
 *     record\t<run>\t<n>\t<contig>:<start>-<end>\t<ID>
 *     group\t<run>\t<n>\t<ID,ID,...>\t<window contig:start-end>\t<window bases>
 *     error\t<run>\t<exception class>:<message, the work directory masked as <dir>>
 *
 * Usage: MultiVariantWalkerDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.MultiVariantWalker;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedByOverlap;
import org.broadinstitute.hellbender.engine.MultiVariantWalkerGroupedOnStart;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class MultiVariantWalkerDump {

    static final int CHR1 = 1000;
    static final int CHR2 = 900;

    static final String FORMAT_LINES =
            "##fileformat=VCFv4.2\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n";

    static final String CONTIGS =
            "##contig=<ID=chr1,length=1000>\n"
                    + "##contig=<ID=chr2,length=900>\n";

    static final String COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n";

    /** The run being printed, since the walkers are constructed by the command line. */
    static String run = "";

    /** The record or group counter, reset per run. */
    static int step = 0;

    static Path dir;

    public static void main(final String[] args) throws Exception {
        dir = Path.of("multivariantwalker-dump").toAbsolutePath();
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MultiVariantWalkerDump: the data source, the merge and the two groupings");
        final Path fasta = writeReference();

        // ---- the merge order -------------------------------------------------------------------
        final Path a = vcf("a", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\ta100\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\ta300\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t50\ta50\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path b = vcf("b", true, FORMAT_LINES + CONTIGS,
                "chr1\t200\tb200\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t10\tb10\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path t1 = vcf("t1", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\tt1-100\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t101\tt1-101\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path t2 = vcf("t2", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\tt2-100\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t100\tt2-100b\tA\tG\t50\tPASS\t.\tGT\t0/1");
        final Path t3 = vcf("t3", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\tt3-100\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t101\tt3-101\tA\tC\t50\tPASS\t.\tGT\t0/1");
        probe("interleaved", new Probe(), List.of("-V", a, "-V", b));
        probe("one-input", new Probe(), List.of("-V", b));
        probe("ties-two", new Probe(), List.of("-V", t1, "-V", t2));
        probe("ties-three", new Probe(), List.of("-V", t1, "-V", t2, "-V", t3));
        probe("ties-three-reversed", new Probe(), List.of("-V", t3, "-V", t2, "-V", t1));
        probe("intervals", new Probe(), List.of("-V", a, "-V", b, "-L", "chr1:150-350"));

        // ---- the dictionaries -------------------------------------------------------------------
        final Path d1 = vcf("d1", true, FORMAT_LINES
                        + "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=900>\n",
                "chr1\t10\td1-10\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t10\td1-chr2\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path d2 = vcf("d2", true, FORMAT_LINES
                        + "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr3,length=800>\n",
                "chr1\t20\td2-20\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr3\t5\td2-chr3\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path reordered = vcf("reordered", true, FORMAT_LINES
                        + "##contig=<ID=chr2,length=900>\n##contig=<ID=chr1,length=1000>\n",
                "chr2\t5\tr-chr2\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t5\tr-chr1\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path longer = vcf("longer", true, FORMAT_LINES + "##contig=<ID=chr1,length=2000>\n",
                "chr1\t30\tlonger-30\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path md5a = vcf("md5a", true, FORMAT_LINES + "##contig=<ID=chr1,length=1000,md5=aaaa>\n",
                "chr1\t40\tmd5a-40\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path md5b = vcf("md5b", true, FORMAT_LINES + "##contig=<ID=chr1,length=1000,md5=bbbb>\n",
                "chr1\t50\tmd5b-50\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path noContigsIndexed = vcf("no-contigs-indexed", true, FORMAT_LINES,
                "chr1\t60\tnci-60\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t60\tnci-chr2\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path noContigs = vcf("no-contigs", false, FORMAT_LINES,
                "chr1\t70\tnc-70\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path noContigsToo = vcf("no-contigs-too", false, FORMAT_LINES,
                "chr1\t80\tnct-80\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path unknownContig = vcf("unknown-contig", false, FORMAT_LINES + CONTIGS,
                "chr9\t5\tu-chr9\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path unsorted = vcf("unsorted", false, FORMAT_LINES + CONTIGS,
                "chr1\t500\tun-500\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t150\tun-150\tA\tC\t50\tPASS\t.\tGT\t0/1");
        probe("dictionaries-union", new Probe(), List.of("-V", d1, "-V", d2));
        probe("dictionaries-union-reversed", new Probe(), List.of("-V", d2, "-V", d1));
        probe("dictionaries-reordered", new Probe(), List.of("-V", d1, "-V", reordered));
        probe("length-mismatch", new Probe(), List.of("-V", a, "-V", longer));
        probe("length-mismatch-not-validated", new Probe(),
                List.of("-V", a, "-V", longer, "--disable-sequence-dictionary-validation"));
        probe("md5-mismatch", new Probe(), List.of("-V", md5a, "-V", md5b));
        probe("md5-one-side", new Probe(), List.of("-V", md5a, "-V", a));
        probe("no-contigs-indexed", new Probe(), List.of("-V", noContigsIndexed, "-V", a));
        probe("no-contigs-indexed-alone", new Probe(), List.of("-V", noContigsIndexed));
        probe("no-contigs-beside-one", new Probe(), List.of("-V", noContigs, "-V", a));
        probe("no-dictionary-anywhere", new Probe(), List.of("-V", noContigs, "-V", noContigsToo));
        probe("unknown-contig", new Probe(), List.of("-V", a, "-V", unknownContig));
        probe("unsorted", new Probe(), List.of("-V", unsorted, "-V", b));

        // ---- the inputs and the samples -----------------------------------------------------------
        final Path samplesA = vcf("samples-a", true, FORMAT_LINES + CONTIGS,
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts2\ts1",
                "chr1\t10\tsa-10\tA\tC\t50\tPASS\t.\tGT\t0/1\t0/0");
        final Path samplesB = vcf("samples-b", true, FORMAT_LINES + CONTIGS
                        + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts3",
                "chr1\t20\tsb-20\tA\tC\t50\tPASS\tDP=3\tGT\t0/1\t0/0");
        final Path samplesC = vcf("samples-c", true, FORMAT_LINES + CONTIGS
                        + "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth as a float\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts4",
                "chr1\t30\tsc-30\tA\tC\t50\tPASS\tDP=3.5\tGT\t0/1");
        probe("samples", new Probe(), List.of("-V", samplesA, "-V", samplesB));
        probe("header-type-conflict", new Probe(), List.of("-V", samplesB, "-V", samplesC));
        probe("duplicate-input", new Probe(), List.of("-V", a, "-V", a));
        probe("same-file-two-tags", new Probe(), List.of("-V:one", a, "-V:two", a));
        probe("throws-in-apply", new Throwing(), List.of("-V", a, "-V", b));

        // ---- grouped on start ------------------------------------------------------------------
        final Path g1 = vcf("g1", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\tg1-100\tACGTA\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t105\tg1-105\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t130\tg1-130\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t5\tg1-chr2\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path g2 = vcf("g2", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\tg2-100\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t101\tg2-101\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t998\tg2-998\tACG\tA\t50\tPASS\t.\tGT\t0/1");
        final List<Object> grouped = List.of("-V", g1, "-V", g2, "-R", fasta);
        onStart("on-start-default", grouped, List.of());
        onStart("on-start-distance-5", grouped, List.of("--combine-variants-distance", "5"));
        onStart("on-start-distance-40", grouped, List.of("--combine-variants-distance", "40"));
        onStart("on-start-max-distance-3", grouped,
                List.of("--combine-variants-distance", "40", "--max-distance", "3"));
        onStart("on-start-ref-padding-3", grouped, List.of("--ref-padding", "3"));
        onStart("on-start-intervals", grouped, List.of("-L", "chr1:101-200"));
        onStart("on-start-ignore-outside", grouped,
                List.of("-L", "chr1:101-200", "--ignore-variants-starting-outside-interval"));
        onStart("on-start-empty", grouped, List.of("-L", "chr1:400-500"));

        // ---- grouped by overlap ----------------------------------------------------------------
        final Path o1 = vcf("o1", true, FORMAT_LINES + CONTIGS,
                "chr1\t100\to1-100\tACGTACGTAC\tA\t50\tPASS\t.\tGT\t0/1",
                "chr1\t120\to1-120\tA\t<NON_REF>\t.\tPASS\tEND=140\tGT\t0/0",
                "chr1\t200\to1-200\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t1\to1-chr2\tA\tC\t50\tPASS\t.\tGT\t0/1");
        final Path o2 = vcf("o2", true, FORMAT_LINES + CONTIGS,
                "chr1\t105\to2-105\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t109\to2-109\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t110\to2-110\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr1\t140\to2-140\tA\tC\t50\tPASS\t.\tGT\t0/1",
                "chr2\t1\to2-chr2\tAC\tA\t50\tPASS\t.\tGT\t0/1");
        final List<Object> overlapping = List.of("-V", o1, "-V", o2, "-R", fasta);
        byOverlap("overlap-default", overlapping, List.of());
        byOverlap("overlap-ignore-reference-blocks", overlapping, List.of("--ignore-reference-blocks"));
        byOverlap("overlap-ignore-outside", overlapping,
                List.of("-L", "chr1:105-300", "--ignore-variants-starting-outside-interval"));
    }

    static void probe(final String label, final MultiVariantWalker walker, final List<Object> argv) {
        invoke(label, walker, argv);
    }

    static void onStart(final String label, final List<Object> inputs, final List<Object> extra) {
        final List<Object> argv = new ArrayList<>(inputs);
        argv.addAll(extra);
        invoke(label, new OnStartProbe(), argv);
    }

    static void byOverlap(final String label, final List<Object> inputs, final List<Object> extra) {
        final List<Object> argv = new ArrayList<>(inputs);
        argv.addAll(extra);
        invoke(label, new OverlapProbe(), argv);
    }

    static void invoke(final String label, final MultiVariantWalker walker, final List<Object> argv) {
        run = label;
        step = 0;
        final String[] strings = argv.stream().map(Object::toString).toArray(String[]::new);
        try {
            walker.instanceMain(strings);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()))));
        }
    }

    /** What the data source reports before the first record. */
    static void describe(final VCFHeader header, final SAMSequenceDictionary dictionary,
                         final java.util.SortedSet<String> samples) {
        System.out.printf("dict\t%s\t%s%n", run, dictionary == null ? "null"
                : dictionary.getSequences().stream()
                    .map(r -> r.getSequenceName() + ":" + r.getSequenceLength())
                    .collect(Collectors.joining(",")));
        System.out.printf("samples\t%s\t%s%n", run, String.join(",", samples));
        final StringBuilder lines = new StringBuilder();
        for (final VCFHeaderLine line : header.getMetaDataInSortedOrder()) {
            lines.append("##").append(line).append("\n");
        }
        System.out.printf("header\t%s\t%s%n", run, ReferenceQueryDump.escape(masked(lines.toString())));
    }

    static String render(final VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd();
    }

    /** Every record in the order the merge hands it over. */
    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class Probe extends MultiVariantWalker {
        @Override
        public void onTraversalStart() {
            describe(getHeaderForVariants(), getSequenceDictionaryForDrivingVariants(),
                    getSamplesForVariants());
        }

        @Override
        public void apply(final VariantContext variant, final ReadsContext readsContext,
                          final ReferenceContext referenceContext, final FeatureContext featureContext) {
            System.out.printf("record\t%s\t%d\t%s\t%s%n", run, step++, render(variant), variant.getID());
        }
    }

    /** A record whose ID is b200 throws inside apply, which the traversal wraps. */
    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class Throwing extends Probe {
        @Override
        public void apply(final VariantContext variant, final ReadsContext readsContext,
                          final ReferenceContext referenceContext, final FeatureContext featureContext) {
            if (variant.getID().equals("b200")) {
                throw new IllegalStateException("probe refused " + variant.getID());
            }
            super.apply(variant, readsContext, referenceContext, featureContext);
        }
    }

    static void group(final List<VariantContext> variants, final ReferenceContext reference) {
        System.out.printf("group\t%s\t%d\t%s\t%s\t%s%n", run, step++,
                variants.stream().map(VariantContext::getID).collect(Collectors.joining(",")),
                reference.getWindow().getContig() + ":" + reference.getWindow().getStart()
                        + "-" + reference.getWindow().getEnd(),
                new String(reference.getBases(), StandardCharsets.US_ASCII));
    }

    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class OnStartProbe extends MultiVariantWalkerGroupedOnStart {
        @Override
        public void onTraversalStart() {
            describe(getHeaderForVariants(), getSequenceDictionaryForDrivingVariants(),
                    getSamplesForVariants());
        }

        @Override
        public void apply(final List<VariantContext> variantContexts,
                          final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
            group(variantContexts, referenceContext);
        }
    }

    @CommandLineProgramProperties(summary = "probe", oneLineSummary = "probe",
            programGroup = VariantEvaluationProgramGroup.class)
    public static class OverlapProbe extends MultiVariantWalkerGroupedByOverlap {
        @Override
        public void onTraversalStart() {
            describe(getHeaderForVariants(), getSequenceDictionaryForDrivingVariants(),
                    getSamplesForVariants());
        }

        @Override
        public void apply(final List<VariantContext> variantContexts,
                          final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
            group(variantContexts, referenceContext);
        }
    }

    /** A VCF: the header lines, then either a #CHROM line among the records or the default one. */
    static Path vcf(final String label, final boolean indexed, final String headerLines,
                    final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(headerLines);
        if (records.length == 0 || !records[0].startsWith("#CHROM")) {
            text.append(COLUMNS);
        }
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        if (indexed) {
            new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        }
        System.out.printf("input\t%s\t%s\t%s%n", label, indexed ? "indexed" : "unindexed",
                ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    /** chr1 of 1000 bases and chr2 of 900, each a repeat of a period that is not a divisor of 60. */
    static Path writeReference() throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder();
        final String[] periods = {"ACGTTGCAAC", "GATTACA"};
        final int[] lengths = {CHR1, CHR2};
        final String[] names = {"chr1", "chr2"};
        for (int c = 0; c < 2; c++) {
            bases.append('>').append(names[c]).append('\n');
            for (int i = 0; i < lengths[c]; i++) {
                bases.append(periods[c].charAt(i % periods[c].length()));
                if (i % 60 == 59 || i == lengths[c] - 1) {
                    bases.append('\n');
                }
            }
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("chr1", CHR1), new SAMSequenceRecord("chr2", CHR2))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        System.out.printf("fasta\t%s%n", ReferenceQueryDump.escape(bases.toString()));
        return fasta;
    }

    static String masked(final String text) {
        return text.replace(dir.toString(), "<dir>");
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
