/*
 * `Mutect2FilteringEngine.buildFiltersList` and the mode-dependent getters on
 * `M2FiltersArgumentCollection`, taken from the reference.
 *
 * Which filters a run builds is decided once, at construction, from two boolean flags, and neither
 * the engine nor its stats file will say which. Five behaviours this is built to catch.
 *
 *   - THE LIST IS BUILT IN ONE ORDER AND THAT ORDER IS THE ENGINE'S. `ErrorProbabilities` collects
 *     it into a `LinkedHashMap`, so the construction order is the iteration order downstream;
 *   - MITOCHONDRIAL MODE DROPS SIX FILTERS AND ADDS NONE; MICROBIAL MODE DROPS THE SAME SIX AND ADDS
 *     ONE BACK, so `PolymeraseSlippageFilter` sits in a DIFFERENT POSITION in microbial mode than in
 *     the default one;
 *   - THE COMMENT CLAIMS A CONDITION THE CODE DOES NOT HAVE. "Normal Artifact Filter doesn't apply
 *     to mitochondria because we are not comparing tumor and normal" sits directly above an
 *     unguarded `filters.add(new NormalArtifactFilter(...))`;
 *   - `getMinMedianMappingQuality()` WRITES TO THE FIELD IT READS. It memoises on first call, so a
 *     collection asked once before `microbial` is set keeps the non-microbial default for ever;
 *   - AND `getLogSnvPrior()` COMPARES AGAINST THE DEFAULT BY VALUE, so passing the default
 *     explicitly is indistinguishable from not passing it and silently yields the mitochondrial
 *     prior instead.
 *
 * `ReadOrientationFilter` needs a tar.gz of artifact priors, which is a fixture rather than an
 * argument; the branch is recorded here and not exercised.
 *
 * Output:
 *
 *     list\t<label>\t<count>=<class,class,...>
 *     argument\t<label>\t<name>=<value>
 *
 * Usage: FilterListByModeDump
 */

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;

import java.io.File;
import java.lang.reflect.Field;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class FilterListByModeDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# FilterListByModeDump: which filters a mode builds");

        list("default", false, false);
        list("mitochondria", true, false);
        list("microbial", false, true);
        list("both-modes", true, true);

        // The mode-dependent getters, each on its own fresh collection.
        argument("mapping-quality-default", fresh(false, false), "minMedianMappingQuality");
        argument("mapping-quality-mitochondria", fresh(true, false), "minMedianMappingQuality");
        argument("mapping-quality-microbial", fresh(false, true), "minMedianMappingQuality");
        // Set explicitly: the memoising getter leaves it alone.
        final M2FiltersArgumentCollection explicit = fresh(false, true);
        explicit.minMedianMappingQuality = 42;
        argument("mapping-quality-explicit", explicit, "minMedianMappingQuality");

        // THE GETTER REMEMBERS: asked before the flag is set, it keeps the first answer.
        final M2FiltersArgumentCollection remembered = fresh(false, false);
        System.out.printf("argument\tmapping-quality-asked-first\tminMedianMappingQuality=%d%n",
                remembered.getMinMedianMappingQuality());
        remembered.microbial = true;
        System.out.printf("argument\tmapping-quality-asked-again\tminMedianMappingQuality=%d%n",
                remembered.getMinMedianMappingQuality());
        // And the other way round, which does take the microbial default.
        final M2FiltersArgumentCollection setFirst = fresh(false, false);
        setFirst.microbial = true;
        System.out.printf("argument\tmapping-quality-flag-set-first\tminMedianMappingQuality=%d%n",
                setFirst.getMinMedianMappingQuality());

        // The priors, whose mitochondrial values are chosen by comparing against the default.
        priors("priors-default", fresh(false, false));
        priors("priors-mitochondria", fresh(true, false));
        priors("priors-microbial", fresh(false, true));
        // Explicitly the default value, under mitochondrial mode: indistinguishable from unset.
        final M2FiltersArgumentCollection explicitPriors = fresh(true, false);
        explicitPriors.logSNVPrior = fresh(false, false).logSNVPrior;
        explicitPriors.logIndelPrior = fresh(false, false).logIndelPrior;
        priors("priors-explicitly-the-default", explicitPriors);
        // And explicitly something else, which is kept.
        final M2FiltersArgumentCollection otherPriors = fresh(true, false);
        otherPriors.logSNVPrior = -12.0;
        otherPriors.logIndelPrior = -13.0;
        priors("priors-explicitly-other", otherPriors);
    }

    static M2FiltersArgumentCollection fresh(final boolean mitochondria, final boolean microbial) {
        final M2FiltersArgumentCollection arguments = new M2FiltersArgumentCollection();
        arguments.mitochondria = mitochondria;
        arguments.microbial = microbial;
        return arguments;
    }

    @SuppressWarnings("unchecked")
    static void list(final String label, final boolean mitochondria, final boolean microbial)
            throws Exception {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFHeaderLine("normal_sample", "N1"));
        final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
        final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(
                fresh(mitochondria, microbial), header, new File("no-such-stats-file.tsv"));
        final Field field = Mutect2FilteringEngine.class.getDeclaredField("filters");
        field.setAccessible(true);
        final List<Mutect2Filter> filters = (List<Mutect2Filter>) field.get(engine);
        final String names = filters.stream().map(f -> f.getClass().getSimpleName())
                .collect(Collectors.joining(","));
        System.out.printf("list\t%s\t%d=%s%n", label, filters.size(), names);
    }

    static void argument(final String label, final M2FiltersArgumentCollection arguments,
                         final String name) {
        System.out.printf("argument\t%s\t%s=%d%n", label, name,
                arguments.getMinMedianMappingQuality());
    }

    static void priors(final String label, final M2FiltersArgumentCollection arguments) {
        System.out.printf("argument\t%s\tlogSnvPrior=%s%n", label,
                Double.toString(arguments.getLogSnvPrior()));
        System.out.printf("argument\t%s\tlogIndelPrior=%s%n", label,
                Double.toString(arguments.getLogIndelPrior()));
    }
}
