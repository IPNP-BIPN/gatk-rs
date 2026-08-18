/*
 * `FilteredHaplotypeFilter`, taken from the reference.
 *
 * The ninth of the ten filters the `filter-mutect-calls` golden needs, and the first with state
 * that outlives a record: it accumulates one artifact probability per phased haplotype across the
 * callset and answers from what the PREVIOUS pass learned. Seven behaviours this is built to catch.
 *
 *   - THE FIRST PASS ANSWERS ZERO TO EVERYTHING. `phasedProbabilities` is read and
 *     `accumulatingPhasedProbabilities` is written, and `learnParameters` moves one to the other
 *     between passes, so this filter's output depends on how many passes have run;
 *   - TWO OF THE THREE NAME COMPARISONS ARE `==` ON `String`. Germline and normal-artifact are
 *     matched by reference, the filter's own name by `.equals`. It works only because both sides
 *     are compile-time constants and therefore interned; a filter returning an EQUAL but
 *     non-interned name is silently not matched, which this dump measures directly;
 *   - `.get()` ON AN Optional WITH NO GUARD: a record with no tumour sample is a
 *     `NoSuchElementException` out of a filter;
 *   - THE RECORD FILTERS ITSELF. The distance test includes zero, so a site some other filter
 *     called an artifact in the first pass filters itself in the second. The filter excludes its
 *     own FILTER from the accumulated probability, as its comment says, but not its own LOCUS;
 *   - ONE ENTRY PER TUMOUR GENOTYPE, NOT PER RECORD: two tumour samples sharing a phasing string
 *     accumulate the same locus twice;
 *   - THE READING GENOTYPE IS THE ONE WITH THE GREATEST `AF`, AND A TIE KEEPS THE FIRST,
 *     `Stream.max` reducing with `BinaryOperator.maxBy`;
 *   - AND `learnParameters` MOVES THE REFERENCE rather than copying it.
 *
 * Output:
 *
 *     default\tmaxDistanceToFilteredCallOnSameHaplotype\t<value>
 *     name\thaplotype\t<filterName>,<errorType>,<annotation>,<required annotations>
 *     accumulated\t<label>\t<phasing string>=<[(locus, probability), ...]>
 *     learned\t<label>\t<phasing string>=<[(locus, probability), ...]>
 *     prob\t<label>\t<one error probability per alternate allele>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FilteredHaplotypeFilterDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorProbabilities;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorType;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilteredHaplotypeFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.io.File;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

public class FilteredHaplotypeFilterDump {

    /** A filter whose answer, name and type the dump chooses outright. */
    static class Fixed extends Mutect2Filter {
        private final String name;
        private final ErrorType type;
        private final double probability;

        Fixed(final String name, final ErrorType type, final double probability) {
            this.name = name;
            this.type = type;
            this.probability = probability;
        }

        @Override
        public ErrorType errorType() {
            return type;
        }

        @Override
        public String filterName() {
            return name;
        }

        @Override
        public Optional<String> phredScaledPosteriorAnnotationName() {
            return Optional.empty();
        }

        @Override
        protected List<String> requiredInfoAnnotations() {
            return List.of();
        }

        @Override
        public List<Double> errorProbabilities(final VariantContext vc,
                                               final Mutect2FilteringEngine engine,
                                               final ReferenceContext reference) {
            return List.of(probability);
        }
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# FilteredHaplotypeFilterDump: the two-pass haplotype filter");

        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        System.out.printf("default\tmaxDistanceToFilteredCallOnSameHaplotype\t%d%n",
                arguments.maxDistanceToFilteredCallOnSameHaplotype);

        final FilteredHaplotypeFilter identity =
                new FilteredHaplotypeFilter(arguments.maxDistanceToFilteredCallOnSameHaplotype);
        System.out.printf("name\thaplotype\t%s,%s,%s,%s%n", identity.filterName(), identity.errorType(),
                identity.phredScaledPosteriorAnnotationName().orElse("none"), "none");

        // Before anything is learned, every record is zero.
        final FilteredHaplotypeFilter fresh = new FilteredHaplotypeFilter(100);
        prob("first-pass", fresh, record(100, "0|1", "100_A_C", 0.3));

        // Accumulate three loci on one haplotype, then learn.
        final FilteredHaplotypeFilter filter = new FilteredHaplotypeFilter(100);
        accumulate(filter, record(100, "0|1", "100_A_C", 0.3), artifact(0.8));
        accumulate(filter, record(150, "0|1", "100_A_C", 0.3), artifact(0.4));
        accumulate(filter, record(500, "0|1", "100_A_C", 0.3), artifact(0.9));
        // A different haplotype at the same locus.
        accumulate(filter, record(120, "1|0", "100_A_C", 0.3), artifact(0.7));
        show("accumulated", "three-loci", filter, "accumulatingPhasedProbabilities");
        learn(filter);
        show("learned", "three-loci", filter, "phasedProbabilities");
        show("accumulated", "after-learning", filter, "accumulatingPhasedProbabilities");

        // The distance test, which includes zero: the record filters itself.
        prob("at-its-own-locus", filter, record(100, "0|1", "100_A_C", 0.3));
        // 150 is within 100 of both 100 and 150, so the maximum of 0.8 and 0.4 wins.
        prob("within-of-two", filter, record(150, "0|1", "100_A_C", 0.3));
        // 300 is within 100 of nothing on this haplotype.
        prob("out-of-range", filter, record(300, "0|1", "100_A_C", 0.3));
        // Exactly at the maximum distance, which `<=` accepts.
        prob("exactly-at-the-distance", filter, record(200, "0|1", "100_A_C", 0.3));
        prob("one-past-the-distance", filter, record(201, "0|1", "100_A_C", 0.3));
        // The other haplotype's accumulation is not visible from this one.
        prob("other-haplotype", filter, record(120, "1|0", "100_A_C", 0.3));
        // A phasing string nothing accumulated.
        prob("unknown-haplotype", filter, record(100, "1|1", "999_A_C", 0.3));
        // No PGT, and no PID: the phasing string is absent and the answer is zero.
        prob("no-pgt", filter, unphased(100, null, "100_A_C"));
        prob("no-pid", filter, unphased(100, "0|1", null));
        // No tumour sample at all.
        prob("no-tumour-sample", filter, normalOnly(100));

        // Which tumour genotype is read: the one with the greatest AF, and a tie keeps the first.
        // The two haplotypes carry DIFFERENT accumulated probabilities, so which genotype is read
        // is visible in the answer.
        final FilteredHaplotypeFilter twoSamples = new FilteredHaplotypeFilter(100);
        accumulate(twoSamples, record(100, "0|1", "100_A_C", 0.1), artifact(0.2));
        accumulate(twoSamples, record(100, "1|0", "100_A_C", 0.9), artifact(0.8));
        show("accumulated", "two-tumours", twoSamples, "accumulatingPhasedProbabilities");
        learn(twoSamples);
        prob("greatest-af-wins", twoSamples, twoTumours(100, 0.1, "0|1", 0.9, "1|0"));
        prob("tie-keeps-the-first", twoSamples, twoTumours(100, 0.5, "0|1", 0.5, "1|0"));

        // Two tumour samples on ONE haplotype accumulate the same locus twice.
        final FilteredHaplotypeFilter shared = new FilteredHaplotypeFilter(100);
        accumulate(shared, twoTumours(100, 0.1, "0|1", 0.9, "0|1"), artifact(0.6));
        show("accumulated", "shared-haplotype", shared, "accumulatingPhasedProbabilities");

        // What accumulateDataForLearning takes the maximum over.
        final VariantContext one = record(100, "0|1", "100_A_C", 0.3);
        accumulated("artifact-only", one, List.of(artifactFilter("base_qual", 0.3),
                artifactFilter("map_qual", 0.7)));
        // A NON_SOMATIC filter is excluded whatever its probability.
        accumulated("non-somatic-excluded", one, List.of(artifactFilter("base_qual", 0.3),
                new Fixed(GATKVCFConstants.CONTAMINATION_FILTER_NAME, ErrorType.NON_SOMATIC, 0.99)));
        // The haplotype filter's own name is excluded, by `.equals`.
        accumulated("self-excluded", one, List.of(artifactFilter("base_qual", 0.3),
                artifactFilter(GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME, 0.99)));
        // A germline probability above 0.25 drops the normal-artifact filter.
        accumulated("germline-drops-normal-artifact", one, List.of(
                new Fixed(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, ErrorType.NON_SOMATIC, 0.5),
                artifactFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
                artifactFilter("base_qual", 0.1)));
        // Below the threshold it keeps it.
        accumulated("germline-below-the-threshold", one, List.of(
                new Fixed(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, ErrorType.NON_SOMATIC, 0.1),
                artifactFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
                artifactFilter("base_qual", 0.1)));
        // Exactly at the threshold, which the strict `>` keeps.
        accumulated("germline-at-the-threshold", one, List.of(
                new Fixed(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, ErrorType.NON_SOMATIC, 0.25),
                artifactFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
                artifactFilter("base_qual", 0.1)));
        // AND THE NAME COMPARISON IS `==`: an EQUAL but non-interned germline name is not matched,
        // so the normal-artifact filter is kept and the answer changes.
        accumulated("germline-name-not-interned", one, List.of(
                new Fixed(new String(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME),
                        ErrorType.NON_SOMATIC, 0.5),
                artifactFilter(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME, 0.9),
                artifactFilter("base_qual", 0.1)));
        // The same for the normal-artifact name itself.
        accumulated("normal-artifact-name-not-interned", one, List.of(
                new Fixed(GATKVCFConstants.GERMLINE_RISK_FILTER_NAME, ErrorType.NON_SOMATIC, 0.5),
                artifactFilter(new String(GATKVCFConstants.ARTIFACT_IN_NORMAL_FILTER_NAME), 0.9),
                artifactFilter("base_qual", 0.1)));
        // And the filter's own name, which is compared with `.equals` and IS matched.
        accumulated("self-name-not-interned", one, List.of(
                artifactFilter(new String(GATKVCFConstants.BAD_HAPLOTYPE_FILTER_NAME), 0.99),
                artifactFilter("base_qual", 0.1)));
    }

    static Mutect2Filter artifactFilter(final String name, final double probability) {
        return new Fixed(name, ErrorType.ARTIFACT, probability);
    }

    /** The one filter every accumulation case shares, so its probability is the only variable. */
    static List<Mutect2Filter> artifact(final double probability) {
        return List.of(artifactFilter("base_qual", probability));
    }

    static Mutect2FilteringEngine engine() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "T2", "N1"));
        return new Mutect2FilteringEngine(new M2FiltersArgumentCollection(), header,
                new File("no-such-stats-file.tsv"));
    }

    static Genotype tumour(final String sample, final String pgt, final String pid, final double af) {
        final GenotypeBuilder builder = new GenotypeBuilder(sample, List.of(Allele.REF_A, Allele.ALT_C))
                .AD(new int[] {80, 20}).attribute("AF", List.of(af));
        if (pgt != null) {
            builder.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY, pgt);
        }
        if (pid != null) {
            builder.attribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY, pid);
        }
        return builder.make();
    }

    static VariantContext record(final int start, final String pgt, final String pid, final double af) {
        return new VariantContextBuilder("dump", "chr1", start, start,
                List.of(Allele.REF_A, Allele.ALT_C))
                .genotypes(List.of(tumour("T1", pgt, pid, af),
                        new GenotypeBuilder("N1", List.of(Allele.REF_A, Allele.REF_A))
                                .AD(new int[] {90, 1}).make()))
                .make();
    }

    static VariantContext unphased(final int start, final String pgt, final String pid) {
        return record(start, pgt, pid, 0.3);
    }

    static VariantContext normalOnly(final int start) {
        return new VariantContextBuilder("dump", "chr1", start, start,
                List.of(Allele.REF_A, Allele.ALT_C))
                .genotypes(List.of(new GenotypeBuilder("N1", List.of(Allele.REF_A, Allele.REF_A))
                        .AD(new int[] {90, 1}).make()))
                .make();
    }

    static VariantContext twoTumours(final int start, final double firstAf, final String firstPgt,
                                     final double secondAf, final String secondPgt) {
        return new VariantContextBuilder("dump", "chr1", start, start,
                List.of(Allele.REF_A, Allele.ALT_C))
                .genotypes(List.of(tumour("T1", firstPgt, "100_A_C", firstAf),
                        tumour("T2", secondPgt, "100_A_C", secondAf)))
                .make();
    }

    static void prob(final String label, final FilteredHaplotypeFilter filter, final VariantContext vc) {
        try {
            System.out.printf("prob\t%s\t%s%n", label, filter.errorProbabilities(vc, engine(), null));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    /** `accumulateDataForLearning`, which is protected. */
    static void accumulate(final FilteredHaplotypeFilter filter, final VariantContext vc,
                           final List<Mutect2Filter> filters) throws Exception {
        final Mutect2FilteringEngine engine = engine();
        final ErrorProbabilities probabilities = new ErrorProbabilities(filters, vc, engine, null);
        final Method method = Mutect2Filter.class.getDeclaredMethod("accumulateDataForLearning",
                VariantContext.class, ErrorProbabilities.class, Mutect2FilteringEngine.class);
        method.setAccessible(true);
        method.invoke(filter, vc, probabilities, engine);
    }

    /** One accumulation, printed: what the maximum over the other filters came to. */
    static void accumulated(final String label, final VariantContext vc,
                            final List<Mutect2Filter> filters) throws Exception {
        final FilteredHaplotypeFilter filter = new FilteredHaplotypeFilter(100);
        accumulate(filter, vc, filters);
        show("accumulated", label, filter, "accumulatingPhasedProbabilities");
    }

    /** `learnParametersAndClearAccumulatedData`, which is protected. */
    static void learn(final FilteredHaplotypeFilter filter) throws Exception {
        final Method method = Mutect2Filter.class.getDeclaredMethod("learnParametersAndClearAccumulatedData");
        method.setAccessible(true);
        method.invoke(filter);
    }

    @SuppressWarnings("unchecked")
    static void show(final String kind, final String label, final FilteredHaplotypeFilter filter,
                     final String fieldName) throws Exception {
        final Field field = FilteredHaplotypeFilter.class.getDeclaredField(fieldName);
        field.setAccessible(true);
        final Map<String, List<?>> map = (Map<String, List<?>>) field.get(filter);
        if (map.isEmpty()) {
            System.out.printf("%s\t%s\t(empty)%n", kind, label);
            return;
        }
        final List<String> keys = new ArrayList<>(map.keySet());
        keys.sort(String::compareTo);
        for (final String key : keys) {
            System.out.printf("%s\t%s\t%s=%s%n", kind, label, key, map.get(key));
        }
    }
}
