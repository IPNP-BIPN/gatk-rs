/*
 * FilterMutectCalls run end to end, taken from the reference.
 *
 * The tool the whole Mutect filtering stack exists for: a VCF of unfiltered calls and its `.stats`
 * table in, a filtered VCF and a filtering-stats file out. Six behaviours this is built to catch,
 * none of which the engine-level suites could reach.
 *
 *   - SIX FILTERS ARE BUILT ONLY WHEN THE RUN IS NOT MITOCHONDRIAL. ClusteredEventsFilter,
 *     MultiallelicFilter, FragmentLengthFilter, PolymeraseSlippageFilter, FilteredHaplotypeFilter
 *     and GermlineFilter sit inside an `if (!MTFAC.mitochondria)`, so the same input filtered twice
 *     in the two modes gives two different FILTER columns AND two different stats files. This is the
 *     only way the guard is visible: `buildFiltersList` is private and the stats file names only the
 *     filters that FIRED;
 *   - THE TOOL MAKES FOUR PASSES OVER THE INPUT, two of them learning, one for the threshold and one
 *     for calling, so the filters applied to the FIRST record come from a model that has already
 *     seen the LAST. A one-record input and a many-record input filter the same record differently;
 *   - THE THRESHOLD IS LEARNED IN A PASS OF ITS OWN, after the parameters have stopped moving, so
 *     the threshold in the stats file corresponds exactly to the filters that were applied;
 *   - THE HEADER IS REWRITTEN: Mutect2's `filtering_status` line is dropped and replaced under the
 *     same key, every ##FILTER line is added whether or not its filter runs, and AS_FilterStatus and
 *     STRQ arrive as ##INFO;
 *   - THE STATS FILE IS WRITTEN AFTER THE CALLING PASS and its rows are the filters that fired, with
 *     the false-positive and false-negative counts a real run produces rather than the hand-built
 *     ones the `filter-stats` suite pins;
 *   - AND A MISSING STATS TABLE IS A UserException NAMING THE FILE, not a silent default.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     stats\t<label>\t<the whole mutect stats table, escaped>
 *     vcfline\t<run>\t<one record line of the output vcf, escaped>
 *     header\t<run>\t<one ##FILTER, ##INFO or ##filtering_status line of the output vcf, escaped>
 *     filtering\t<run>\t<one line of the filtering-stats file, escaped>
 *     error\t<run>\t<exception class>:<message>
 *
 * Usage: FilterMutectCallsDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FilterMutectCallsDump {

    /** One kilobase of chr1, which is all the reference any of these records needs. */
    static final String CHR1 = "AC".repeat(500);

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
                    + "##filtering_status=Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove false positives.\n"
                    + "##INFO=<ID=CONTQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to contamination\">\n"
                    + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n"
                    + "##INFO=<ID=GERMQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not germline variants\">\n"
                    + "##INFO=<ID=ROQ,Number=1,Type=Float,Description=\"Phred-scaled qualities that alt allele are not due to read orientation artifact\">\n"
                    + "##INFO=<ID=SEQQ,Number=1,Type=Integer,Description=\"Phred-scaled quality that alt alleles are not sequencing errors\">\n"
                    + "##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description=\"Phred-scaled quality of strand bias artifact\">\n"
                    + "##INFO=<ID=ECNT,Number=1,Type=Integer,Description=\"Number of events in this haplotype\">\n"
                    + "##INFO=<ID=MBQ,Number=R,Type=Integer,Description=\"median base quality\">\n"
                    + "##INFO=<ID=MFRL,Number=R,Type=Integer,Description=\"median fragment length\">\n"
                    + "##INFO=<ID=MMQ,Number=R,Type=Integer,Description=\"median mapping quality\">\n"
                    + "##INFO=<ID=MPOS,Number=A,Type=Integer,Description=\"median distance from end of read\">\n"
                    + "##INFO=<ID=NALOD,Number=A,Type=Float,Description=\"Negative log 10 odds of artifact in normal\">\n"
                    + "##INFO=<ID=NLOD,Number=A,Type=Float,Description=\"Normal log 10 likelihood ratio of diploid het or hom alt genotypes\">\n"
                    + "##INFO=<ID=POPAF,Number=A,Type=Float,Description=\"negative log 10 population allele frequencies of alt alleles\">\n"
                    + "##INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Log 10 likelihood ratio score of variant existing versus not existing\">\n"
                    + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                    + "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fractions\">\n"
                    + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n"
                    + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    + "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"strand bias counts\">\n"
                    + "##contig=<ID=chr1,length=1000>\n"
                    + "##normal_sample=N1\n"
                    + "##tumor_sample=T1\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT1\tN1\n";

    /** The FORMAT of every record here, and one tumour and one normal genotype. */
    static String genotypes(final int tumorRef, final int tumorAlt, final int normalRef,
                            final int normalAlt, final String strandBias) {
        final double fraction = (double) tumorAlt / (tumorRef + tumorAlt);
        return "GT:AD:AF:DP:SB\t0/1:" + tumorRef + "," + tumorAlt + ":"
                + String.format("%.3f", fraction) + ":" + (tumorRef + tumorAlt) + ":" + strandBias
                + "\t0/0:" + normalRef + "," + normalAlt + ":0.010:" + (normalRef + normalAlt)
                + ":" + strandBias;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("filter-mutect-calls-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FilterMutectCallsDump: the tool end to end, in two modes");

        final Path fasta = writeReference(dir);

        // Eight records: three clean and clonal, one weak, one with a poor base quality, one with a
        // poor read position, one clustered, and one germline-looking.
        final Path calls = writeVcf(dir, "calls",
                "chr1\t100\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=30.0\t" + genotypes(80, 20, 99, 1, "20,20,10,10"),
                "chr1\t200\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=40.0\t" + genotypes(78, 22, 99, 1, "20,20,10,10"),
                "chr1\t300\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=35.0\t" + genotypes(79, 21, 99, 1, "20,20,10,10"),
                // Weak evidence: a TLOD barely above nothing.
                "chr1\t400\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=3.0\t" + genotypes(97, 3, 99, 1, "40,40,1,2"),
                // A poor median base quality on the alternate.
                "chr1\t500\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,5;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=30.0\t" + genotypes(80, 20, 99, 1, "20,20,10,10"),
                // A median read position at the very end of the reads.
                "chr1\t600\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=0;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=30.0\t" + genotypes(80, 20, 99, 1, "20,20,10,10"),
                // Nine events in one haplotype, which only a non-mitochondrial run can filter.
                "chr1\t700\t.\tA\tC\t.\t.\tDP=100;ECNT=9;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=30.0\t" + genotypes(80, 20, 99, 1, "20,20,10,10"),
                // Germline-looking: a common allele and a normal that carries it, which again only a
                // non-mitochondrial run has a filter for.
                "chr1\t800\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=-3.0;NLOD=-5.0;POPAF=0.1;TLOD=30.0\t" + genotypes(50, 50, 60, 40, "20,20,25,25"));

        // The same input cut to one record, which is what shows the learning passes mattering.
        final Path single = writeVcf(dir, "single",
                "chr1\t100\t.\tA\tC\t.\t.\tDP=100;ECNT=1;MBQ=30,30;MFRL=300,300;MMQ=60,60;MPOS=25;"
                        + "NALOD=2.0;NLOD=5.0;POPAF=6.0;TLOD=30.0\t" + genotypes(80, 20, 99, 1, "20,20,10,10"));

        final Path stats = writeStats(dir, "calls", 1000000.0);
        writeStats(dir, "single", 1000000.0);
        // A stats table saying there were no callable sites at all, which switches the empirical
        // priors off and warns rather than refusing.
        final Path noCallable = writeStats(dir, "no-callable", 0.0);

        run(dir, "default", calls, fasta, stats, List.of());
        run(dir, "mitochondria", calls, fasta, stats, List.of("--mitochondria-mode", "true"));
        run(dir, "single-record", single, fasta, dir.resolve("single.stats"), List.of());
        run(dir, "no-callable-sites", calls, fasta, noCallable, List.of());
        // A threshold strategy other than the default, which changes the threshold and therefore
        // every FILTER column.
        run(dir, "constant-threshold", calls, fasta, stats,
                List.of("--threshold-strategy", "CONSTANT", "--initial-threshold", "0.5"));
        // A stats table that is not there.
        run(dir, "missing-stats", calls, fasta, dir.resolve("no-such.stats"), List.of());
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

    /** The `.stats` table Mutect2 writes beside its VCF: two columns and one row. */
    static Path writeStats(final Path dir, final String label, final double callableSites)
            throws Exception {
        final String text = "statistic\tvalue\ncallable\t" + Double.toString(callableSites) + "\n";
        final Path file = dir.resolve(label + ".stats");
        Files.writeString(file, text, StandardCharsets.UTF_8);
        System.out.printf("stats\t%s\t%s%n", label, ReferenceQueryDump.escape(text));
        return file;
    }

    static void run(final Path dir, final String label, final Path variants, final Path fasta,
                    final Path stats, final List<String> extra) {
        final Path output = dir.resolve(label + ".out.vcf");
        final Path filtering = dir.resolve(label + ".filtering.tsv");
        final List<String> all = new ArrayList<>(List.of(
                "-V", variants.toString(),
                "-R", fasta.toString(),
                "--stats", stats.toString(),
                "--filtering-stats", filtering.toString(),
                "-O", output.toString()));
        all.addAll(extra);
        try {
            new FilterMutectCalls().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            for (Throwable cause = e.getCause(); cause != null; cause = cause.getCause()) {
                System.out.printf("cause\t%s\t%s:%s%n", label, cause.getClass().getName(),
                        ReferenceQueryDump.escape(String.valueOf(cause.getMessage())));
            }
            return;
        }
        try {
            for (final String line : Files.readAllLines(output, StandardCharsets.UTF_8)) {
                if (line.startsWith("##FILTER=") || line.startsWith("##INFO=<ID=AS_FilterStatus")
                        || line.startsWith("##INFO=<ID=STRQ") || line.startsWith("##filtering_status")) {
                    System.out.printf("header\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                } else if (!line.startsWith("#")) {
                    System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
                }
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s-read-vcf\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
        }
        try {
            for (final String line : Files.readAllLines(filtering, StandardCharsets.UTF_8)) {
                System.out.printf("filtering\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s-read-stats\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        Files.writeString(dir.resolve("reference.dict"),
                "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + CHR1.length() + "\n",
                StandardCharsets.UTF_8);
        return fasta;
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
