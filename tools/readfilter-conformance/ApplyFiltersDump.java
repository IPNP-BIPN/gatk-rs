/*
 * `Mutect2FilteringEngine.applyFiltersAndAccumulateOutputStats`, taken from the reference.
 *
 * The step that turns a record's per-filter probabilities into the FILTER column and the
 * `AS_FilterStatus` annotation. Seven behaviours this is built to catch.
 *
 *   - A FILTER CAN FIRE WITHOUT BEING NAMED IN THE FILTER COLUMN. Only the entries whose probability
 *     reaches `min(maxErrorProb, MIN_REPORTABLE_ERROR_PROBABILITY = 0.1)` are written out;
 *   - `SITE` IS A PLACEHOLDER, NOT A FILTER. Every allele that passes a filter is recorded as
 *     `SITE`, and `getDistinctFiltersForAllele` removes it when the allele has any real filter and
 *     adds it back when it has none;
 *   - A SYMBOLIC ALTERNATE DOES NOT CONSUME FROM THE ITERATOR, so where the symbolic alleles are
 *     decides which filter string each real allele gets;
 *   - A SITE-LEVEL FILTER IS DERIVED FROM THE ALLELE-LEVEL ONES ONLY WHEN EVERY ALLELE AGREES;
 *   - `FAIL` IS THE SITE'S ANSWER WHEN EVERY ALLELE IS FILTERED FOR A DIFFERENT REASON, and only
 *     when no site filter has already been recorded;
 *   - THE PHRED-SCALED ANNOTATION IS WRITTEN ONLY WHEN THE FILTER'S REQUIRED ANNOTATIONS ARE
 *     PRESENT, a second and independent check of what `ErrorProbabilities` already applied;
 *   - AND THE THRESHOLD IS CLAMPED INTO [1e-10, 1 - 1e-10] before use.
 *
 * The engine's filter list is replaced by reflection with filters whose name, type, arity and
 * probabilities this dump chooses outright, so the input to the step is exact.
 *
 * Output:
 *
 *     applied\t<label>\t<FILTER column>|<AS_FilterStatus>
 *     attribute\t<label>\t<key>=<value>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ApplyFiltersDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ErrorType;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilteringOutputStats;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2AlleleFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2VariantFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;

public class ApplyFiltersDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT_C = Allele.create("C", false);
    static final Allele ALT_G = Allele.create("G", false);
    static final Allele SYMBOLIC = Allele.create("<NON_REF>", false);

    /**
     * A per-allele filter whose answers this dump chooses. It must extend `Mutect2AlleleFilter`
     * rather than `Mutect2Filter`, because `ErrorProbabilities` partitions on
     * `Mutect2VariantFilter.class.isAssignableFrom(...)` and would otherwise call it a site filter.
     */
    static class FixedAllele extends Mutect2AlleleFilter {
        private final String name;
        private final List<Double> probabilities;

        FixedAllele(final String name, final List<Double> probabilities) {
            this.name = name;
            this.probabilities = probabilities;
        }

        @Override
        public ErrorType errorType() {
            return ErrorType.ARTIFACT;
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
        public List<String> requiredInfoAnnotations() {
            return List.of();
        }

        @Override
        protected List<Double> calculateErrorProbabilityForAlleles(final VariantContext vc,
                                                                   final Mutect2FilteringEngine engine,
                                                                   final ReferenceContext reference) {
            return probabilities;
        }
    }

    /** A per-site filter, whose one probability the base class copies to every alternate allele. */
    static class FixedSite extends Mutect2VariantFilter {
        private final String name;
        private final double probability;
        private final String annotation;
        private final List<String> required;

        FixedSite(final String name, final double probability, final String annotation,
                  final List<String> required) {
            this.name = name;
            this.probability = probability;
            this.annotation = annotation;
            this.required = required;
        }

        @Override
        public ErrorType errorType() {
            return ErrorType.NON_SOMATIC;
        }

        @Override
        public String filterName() {
            return name;
        }

        @Override
        public Optional<String> phredScaledPosteriorAnnotationName() {
            return Optional.ofNullable(annotation);
        }

        @Override
        public List<String> requiredInfoAnnotations() {
            return required;
        }

        @Override
        protected double calculateErrorProbability(final VariantContext vc,
                                                   final Mutect2FilteringEngine engine,
                                                   final ReferenceContext reference) {
            return probability;
        }
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# ApplyFiltersDump: turning probabilities into the FILTER column");

        final VariantContext biallelic = record(List.of(REF, ALT_C));
        final VariantContext triallelic = record(List.of(REF, ALT_C, ALT_G));
        final VariantContext withSymbolic = record(List.of(REF, ALT_C, SYMBOLIC));
        final VariantContext symbolicFirst = record(List.of(REF, SYMBOLIC, ALT_C));

        // Nothing above the threshold: the record passes and every allele is SITE.
        apply("everything-passes", biallelic, 0.5,
                List.of(allele("base_qual", List.of(0.1, 0.1))));

        // One allele filtered on a triallelic record: the other stays SITE.
        apply("one-allele-filtered", triallelic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1))));

        // Every allele filtered by the same filter, which becomes a site-level filter.
        apply("every-allele-same-filter", triallelic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.9))));

        // Every allele filtered, by DIFFERENT filters: no single site filter, so FAIL.
        apply("every-allele-different-filters", triallelic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1)),
                        allele("map_qual", List.of(0.1, 0.9))));

        // One allele carries two filters at once.
        apply("one-allele-two-filters", triallelic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1)),
                        allele("map_qual", List.of(0.9, 0.1))));

        // A site filter alongside allele filters, which suppresses the FAIL check.
        apply("site-and-allele-filters", triallelic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1)),
                        allele("map_qual", List.of(0.1, 0.9)),
                        site("germline", 0.95, null)));

        // A RECORD WITH NO PER-ALLELE FILTER AT ALL CRASHES: `orderedASFilterStrings` walks an empty
        // iterator. Every per-allele filter answers an empty list when its annotations are missing,
        // so this is a record's worth of missing annotations away.
        apply("only-site-filters", biallelic, 0.5, List.of(site("germline", 0.9, null)));

        // THE REPORTING FLOOR: a filter above the threshold but below 0.1 is applied and NOT named
        // when a stronger one is present. A passing allele filter is present so the walk succeeds.
        final Mutect2Filter passing = allele("base_qual", List.of(0.0));
        apply("below-the-reporting-floor", biallelic, 0.01,
                List.of(passing, site("germline", 0.99, null), site("contamination", 0.05, null)));
        // And alone, where the floor is the maximum itself.
        apply("below-the-floor-alone", biallelic, 0.01,
                List.of(passing, site("contamination", 0.05, null)));

        // The phred-scaled annotation, written only when the required annotations are there.
        apply("annotation-written", biallelic, 0.5,
                List.of(passing, site("germline", 0.9, "GERMQ", List.of())));
        apply("annotation-required-annotation-missing", biallelic, 0.5,
                List.of(passing, site("germline", 0.9, "GERMQ", List.of("NOT_THERE"))));
        // The annotation is written even when the filter did not fire.
        apply("annotation-without-the-filter", biallelic, 0.5,
                List.of(passing, site("germline", 0.1, "GERMQ", List.of())));

        // A symbolic alternate takes SITE without consuming a filter string.
        apply("symbolic-last", withSymbolic, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1))));
        apply("symbolic-first", symbolicFirst, 0.5,
                List.of(allele("base_qual", List.of(0.9, 0.1))));

        // An empty list from a filter is dropped rather than counted.
        apply("empty-list", triallelic, 0.5,
                List.of(allele("base_qual", List.of()), allele("map_qual", List.of(0.9, 0.9))));

        // No filters at all.
        apply("no-filters", triallelic, 0.5, List.of());

        // The threshold's clamp, at each end.
        apply("threshold-zero", biallelic, 0.0, List.of(allele("base_qual", List.of(0.0))));
        apply("threshold-one", biallelic, 1.0, List.of(allele("base_qual", List.of(1.0))));
    }

    static VariantContext record(final List<Allele> alleles) {
        return new VariantContextBuilder("dump", "chr1", 100, 100, alleles).make();
    }

    static Mutect2Filter allele(final String name, final List<Double> probabilities) {
        return new FixedAllele(name, probabilities);
    }

    static Mutect2Filter site(final String name, final double probability, final String annotation) {
        return site(name, probability, annotation, List.of());
    }

    static Mutect2Filter site(final String name, final double probability, final String annotation,
                              final List<String> required) {
        return new FixedSite(name, probability, annotation, required);
    }

    @SuppressWarnings("unchecked")
    static void apply(final String label, final VariantContext vc, final double threshold,
                      final List<Mutect2Filter> filters) throws Exception {
        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        arguments.initialPosteriorThreshold = threshold;
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine =
                new Mutect2FilteringEngine(arguments, header, new File("no-such-stats-file.tsv"));
        final Field field = Mutect2FilteringEngine.class.getDeclaredField("filters");
        field.setAccessible(true);
        final List<Mutect2Filter> engineFilters = (List<Mutect2Filter>) field.get(engine);
        engineFilters.clear();
        engineFilters.addAll(filters);
        // `FilteringOutputStats` counts by filter OBJECT, and its map is built at construction, so
        // it has to be rebuilt over the list that was just replaced.
        final Field stats = Mutect2FilteringEngine.class.getDeclaredField("filteringOutputStats");
        stats.setAccessible(true);
        stats.set(engine, new FilteringOutputStats(engineFilters));

        try {
            final VariantContext filtered = engine.applyFiltersAndAccumulateOutputStats(vc, null);
            final String column = filtered.isNotFiltered() ? "PASS"
                    : String.join(";", new TreeSet<>(filtered.getFilters()));
            System.out.printf("applied\t%s\t%s|%s%n", label, column,
                    filtered.getAttributeAsString("AS_FilterStatus", "(none)"));
            final List<String> keys = new ArrayList<>(filtered.getAttributes().keySet());
            keys.sort(String::compareTo);
            for (final String key : keys) {
                if (key.equals("AS_FilterStatus")) {
                    continue;
                }
                System.out.printf("attribute\t%s\t%s=%s%n", label, key,
                        filtered.getAttributeAsString(key, ""));
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }
}
