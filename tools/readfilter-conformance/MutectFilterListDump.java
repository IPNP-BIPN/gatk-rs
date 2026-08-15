/*
 * The engine's filter list and the tool's header lines, taken from the reference.
 *
 * `Mutect2FilteringEngine.buildFiltersList` decides which filters exist for a run, and
 * `FilterMutectCalls.onTraversalStart` decides which header lines describe them. The two lists are
 * not the same list. Six behaviours this is built to catch.
 *
 *   - SIX FILTERS ARE BUILT ONLY WHEN THE RUN IS NOT MITOCHONDRIAL. ClusteredEventsFilter,
 *     MultiallelicFilter, FragmentLengthFilter, PolymeraseSlippageFilter, FilteredHaplotypeFilter
 *     and GermlineFilter sit inside an `if (!MTFAC.mitochondria)`, so a mitochondrial run has a
 *     strictly shorter list;
 *   - `ReadOrientationFilter` IS BUILT ONLY WHEN A PRIOR TAR.GZ WAS GIVEN, so its absence is an
 *     argument's absence rather than a default;
 *   - THE STATS FILE NAMES ONLY THE FILTERS THAT FIRED. An engine that has filtered nothing writes
 *     its metadata and its header row and NO filter rows at all, in every mode, so a missing row is
 *     a filter that found nothing rather than a filter that does not exist;
 *   - `MUTECT_FILTER_NAMES` IS A DIFFERENT LIST FROM THE FILTERS and it includes PASS. Every entry
 *     becomes a ##FILTER header line whether or not the filter runs;
 *   - `MUTECT_AS_FILTER_NAMES` IS ONE ENTRY, AS_FilterStatus, and it becomes an ##INFO line rather
 *     than a ##FILTER one;
 *   - AND THE FILTERS DISAGREE ON BOTH AXES THE ENGINE SORTS THEM BY: nine answer per allele and
 *     nine per site, and three error types are represented -- SEQUENCING for the tumour-evidence
 *     filter alone, NON_SOMATIC for contamination and germline, ARTIFACT for the rest -- which is
 *     what decides which probabilities are combined by a maximum and which by independence.
 *
 * The engine's list is private and its stats file names only the filters that FIRED, so what is
 * measured here is every filter's identity, constructed directly, beside the stats file of an engine
 * that has filtered nothing. WHICH filters a mode builds needs the tool run end to end, which is its
 * own slice.
 *
 * Output:
 *
 *     filters\t<label>\t<one filter name per row of the stats file, comma separated>
 *     count\t<label>\t<how many rows the stats file has>
 *     stats\t<label>\t<a whole line of the stats file>
 *     filterline\t<name>\t<the ##FILTER header line it produces>
 *     infoline\t<name>\t<the ##INFO header line it produces>
 *     names\t<label>\t<the constant's entries, comma separated>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MutectFilterListDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.BaseQualityFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ClusteredEventsFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ContaminationFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.DuplicatedAltReadFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilteredHaplotypeFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FragmentLengthFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.GermlineFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MappingQualityFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MinAlleleFractionFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MultiallelicFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2AlleleFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2Filter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NRatioFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.NormalArtifactFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.PanelOfNormalsFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.PolymeraseSlippageFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.ReadPositionFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrandArtifactFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.StrictStrandBiasFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.TumorEvidenceFilter;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.M2FiltersArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.Mutect2FilteringEngine;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public class MutectFilterListDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# MutectFilterListDump: which filters a run builds, and what describes them");

        // The default list.
        engine("default", new M2FiltersArgumentCollection());

        // The mitochondrial list, which is six filters shorter.
        final M2FiltersArgumentCollection mito = new M2FiltersArgumentCollection();
        mito.mitochondria = true;
        engine("mitochondria", mito);

        // A filter whose threshold is set does not change the list, only what it does.
        final M2FiltersArgumentCollection tuned = new M2FiltersArgumentCollection();
        tuned.minMedianBaseQuality = 30;
        tuned.minReadsOnEachStrand = 2;
        tuned.minAf = 0.05;
        engine("tuned", tuned);

        // Every filter the engine can build, constructed directly: its name, its error type and
        // whether it answers per allele. TWO FILTERS SHARE ONE NAME, and the error type is what
        // decides which probabilities are combined with which.
        filter("TumorEvidenceFilter", new TumorEvidenceFilter());
        filter("BaseQualityFilter", new BaseQualityFilter(20));
        filter("MappingQualityFilter", new MappingQualityFilter(30, 5));
        filter("DuplicatedAltReadFilter", new DuplicatedAltReadFilter(0));
        filter("StrandArtifactFilter", new StrandArtifactFilter());
        filter("ContaminationFilter", new ContaminationFilter(List.of(), 0.0));
        filter("StrictStrandBiasFilter", new StrictStrandBiasFilter(0));
        filter("ReadPositionFilter", new ReadPositionFilter(5));
        filter("MinAlleleFractionFilter", new MinAlleleFractionFilter(0.0));
        filter("NormalArtifactFilter", new NormalArtifactFilter(0.001));
        filter("NRatioFilter", new NRatioFilter(Double.POSITIVE_INFINITY));
        filter("PanelOfNormalsFilter", new PanelOfNormalsFilter());
        filter("ClusteredEventsFilter", new ClusteredEventsFilter(2, 2));
        filter("MultiallelicFilter", new MultiallelicFilter(1));
        filter("FragmentLengthFilter", new FragmentLengthFilter(10000));
        filter("PolymeraseSlippageFilter", new PolymeraseSlippageFilter(8, 0.1));
        filter("FilteredHaplotypeFilter", new FilteredHaplotypeFilter(100));
        filter("GermlineFilter", new GermlineFilter(List.of()));

        // Every filter name the header declares, whether or not its filter runs.
        System.out.printf("names\tMUTECT_FILTER_NAMES\t%s%n",
                String.join(",", GATKVCFConstants.MUTECT_FILTER_NAMES));
        System.out.printf("names\tMUTECT_AS_FILTER_NAMES\t%s%n",
                String.join(",", GATKVCFConstants.MUTECT_AS_FILTER_NAMES));

        for (final String name : GATKVCFConstants.MUTECT_FILTER_NAMES) {
            line("filterline", name, GATKVCFHeaderLines.getFilterLine(name));
        }
        for (final String name : GATKVCFConstants.MUTECT_AS_FILTER_NAMES) {
            line("infoline", name, GATKVCFHeaderLines.getInfoLine(name));
        }
        // The INFO line the tool adds beside them, which is not in either list.
        line("infoline", GATKVCFConstants.POLYMERASE_SLIPPAGE_QUAL_KEY,
                GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.POLYMERASE_SLIPPAGE_QUAL_KEY));

        // Mutect2's own line, which the tool strips before adding its replacement under the same key.
        System.out.printf("stats\tfiltering-status-mutect2\t%s%n", ReferenceQueryDump.escape(
                new VCFHeaderLine(FilterMutectCalls.FILTERING_STATUS_VCF_KEY,
                        "Warning: unfiltered Mutect 2 calls.  Please run FilterMutectCalls to remove"
                                + " false positives.").toString()));

        // The filtering-status line the tool writes over Mutect2's own.
        System.out.printf("stats\tfiltering-status-key\t%s%n",
                ReferenceQueryDump.escape(FilterMutectCalls.FILTERING_STATUS_VCF_KEY));
        System.out.printf("stats\tfiltering-status-line\t%s%n", ReferenceQueryDump.escape(
                new VCFHeaderLine(FilterMutectCalls.FILTERING_STATUS_VCF_KEY,
                        "These calls have been filtered by FilterMutectCalls to label false positives"
                                + " with a list of failed filters and true positives with PASS.").toString()));
    }

    /** One engine, written out through the stats file that names its filters. */
    static void engine(final String label, final M2FiltersArgumentCollection arguments) {
        try {
            final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
            lines.add(new VCFHeaderLine("normal_sample", "N1"));
            final VCFHeader header = new VCFHeader(lines, List.of("T1", "N1"));
            final Mutect2FilteringEngine engine = new Mutect2FilteringEngine(arguments, header,
                    new File("no-such-stats-file.tsv"));
            // Nothing has been filtered yet: the stats file has metadata and a header row and no
            // filter rows at all, which is what shows that a row is a firing rather than a filter.
            final Path empty = Files.createTempFile("filtering-stats-empty", ".tsv");
            engine.writeFilteringStats(empty);
            for (final String row : Files.readAllLines(empty)) {
                if (!row.startsWith("#")) {
                    System.out.printf("stats\t%s-empty\t%s%n", label, ReferenceQueryDump.escape(row));
                }
            }
            Files.deleteIfExists(empty);

            final Path stats = Files.createTempFile("filtering-stats", ".tsv");
            engine.writeFilteringStats(stats);
            final List<String> rows = Files.readAllLines(stats);
            Files.deleteIfExists(stats);

            final StringBuilder names = new StringBuilder();
            int count = 0;
            for (final String row : rows) {
                // The metadata rows start with a '#'; the header row names the columns.
                if (row.startsWith("#") || row.startsWith("filter\t")) {
                    continue;
                }
                if (names.length() > 0) {
                    names.append(',');
                }
                names.append(row.split("\t")[0]);
                count++;
            }
            System.out.printf("filters\t%s\t%s%n", label, names);
            System.out.printf("count\t%s\t%d%n", label, count);
            for (final String row : rows) {
                System.out.printf("stats\t%s\t%s%n", label, ReferenceQueryDump.escape(row));
            }
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void filter(final String label, final Mutect2Filter filter) {
        try {
            System.out.printf("filter\t%s\t%s,%s,%s%n", label,
                    ReferenceQueryDump.escape(filter.filterName()),
                    filter.errorType(),
                    filter instanceof Mutect2AlleleFilter ? "per-allele" : "per-site");
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
        }
    }

    static void line(final String kind, final String name, final Object headerLine) {
        System.out.printf("%s\t%s\t%s%n", kind, name,
                ReferenceQueryDump.escape(String.valueOf(headerLine)));
    }
}
