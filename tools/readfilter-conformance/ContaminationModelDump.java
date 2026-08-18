/*
 * ContaminationModel and the tool around it, taken from the reference.
 *
 * The last unblocked member of the reporting-walker archetype. The decomposition under it is
 * already ported and pinned by `kernel-segmentation`; this measures everything above it, from the
 * coverage filter the tool applies to the contamination estimate and its standard error.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE COVERAGE FILTER USES A MEDIAN AND A MEAN, not one or the other: the low threshold is a
 *     ratio of the MEDIAN and the high threshold a ratio of the MEAN, and both are computed AFTER
 *     dropping sites at or below MIN_COVERAGE, so an uncovered site that would have moved the
 *     median is already gone;
 *   - THE MODEL ITERATES THREE TIMES between minor allele fractions and contamination, each round
 *     feeding the previous contamination back in, so a port that solved once lands somewhere else;
 *   - THE LOH SEARCH WALKS THE THRESHOLD DOWN in steps of 0.04 and stops at the FIRST threshold
 *     keeping more than a quarter of the sites, so the answer depends on the order of the walk and
 *     not only on the sites;
 *   - THE STRATEGY CASCADE IS THREE STRATEGIES OVER ONE LOOP: hom alt above 0.25, hom ref down to
 *     0.20, then the unscrupulous hom ref, and it returns the first estimate whose standard error
 *     is small enough relative to the estimate itself;
 *   - THE STANDARD ERROR IS A BINARY SEARCH, not the closed formula: the formula is only the
 *     fallback for when the search fails to bracket a zero;
 *   - THE UNSCRUPULOUS STRATEGY GOES THROUGH A PERCENTILE, so its threshold is commons-math's
 *     interpolation of the 90th percentile and not a sorted index;
 *   - AND THE SEGMENTATION IS PER CONTIG, so a contig with too few het sites for the window is one
 *     segment rather than none.
 *
 * The sites are built from a formula rather than read from a file, so a port can rebuild exactly
 * the same input: no rounding, no RNG, integer counts throughout. The `sites` rows carry them so a
 * port that built a different corpus fails on the corpus rather than on the model.
 *
 * Output:
 *
 *     sites\t<label>\t<index>=<contig>,<start>,<ref>,<alt>,<other>,<allele frequency bits>
 *     coverage\t<label>\t<the starts that survive, comma separated>
 *     segments\t<label>\t<index>=<contig>,<start>,<end>,<site count>
 *     maf\t<label>\t<index>=<contig>,<start>,<end>,<minor allele fraction bits>
 *     contamination\t<label>\t<estimate bits>,<standard error bits>
 *
 * Usage: ContaminationModelDump
 */

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationModel;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationSegmenter;
import org.broadinstitute.hellbender.tools.walkers.contamination.MinorAlleleFractionRecord;
import org.broadinstitute.hellbender.tools.walkers.contamination.PileupSummary;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

public class ContaminationModelDump {

    /** `CalculateContamination.MIN_COVERAGE`. */
    static final int MIN_COVERAGE = 10;
    /** `CalculateContamination.DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD`. */
    static final double LOW_COVERAGE_RATIO_THRESHOLD = 1.0 / 2;
    /** `CalculateContamination.DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD`. */
    static final double HIGH_COVERAGE_RATIO_THRESHOLD = 3.0;

    public static void main(final String[] args) throws Exception {
        System.out.println("# ContaminationModelDump: the contamination model and the tool around it");

        final List<PileupSummary> tumor = new ArrayList<>();
        tumor.addAll(sites("chr1", 200, true));
        tumor.addAll(sites("chr2", 150, false));
        final List<PileupSummary> normal = new ArrayList<>();
        normal.addAll(sites("chr1", 200, false));

        print("tumor", tumor);
        print("normal", normal);

        // The coverage filter, which is the first thing the tool does to either table.
        final List<PileupSummary> filteredTumor = filterSitesByCoverage(tumor);
        final List<PileupSummary> filteredNormal = filterSitesByCoverage(normal);
        coverage("tumor", filteredTumor);
        coverage("normal", filteredNormal);

        // The segmentation, which is where the decomposition is reached.
        segments("tumor", ContaminationSegmenter.findSegments(filteredTumor));
        segments("normal", ContaminationSegmenter.findSegments(filteredNormal));

        // Tumour only: one model, used both to genotype and to segment.
        final ContaminationModel tumorModel = new ContaminationModel(filteredTumor);
        maf("tumor-only", tumorModel.segmentationRecords());
        contamination("tumor-only", tumorModel.calculateContaminationFromHoms(filteredTumor));

        // Matched normal: the normal genotypes, the tumour segments.
        final ContaminationModel normalModel = new ContaminationModel(filteredNormal);
        maf("matched-normal", normalModel.segmentationRecords());
        contamination("matched-normal", normalModel.calculateContaminationFromHoms(filteredTumor));

        // A short table, below the window the segmenter needs, so the whole contig is one segment.
        final List<PileupSummary> short_ = sites("chr3", 40, false);
        print("short", short_);
        final List<PileupSummary> filteredShort = filterSitesByCoverage(short_);
        coverage("short", filteredShort);
        segments("short", ContaminationSegmenter.findSegments(filteredShort));
        final ContaminationModel shortModel = new ContaminationModel(filteredShort);
        maf("short", shortModel.segmentationRecords());
        contamination("short", shortModel.calculateContaminationFromHoms(filteredShort));
    }

    /**
     * The corpus, from a formula. Every count is an integer computed by integer arithmetic, so
     * there is no rounding rule for a port to disagree with.
     *
     * The genotype cycles hom ref, het, het, hom alt, het. When `loh` is set the second half of the
     * contig has its hets at a minor allele fraction of a quarter rather than a half, which is the
     * loss of heterozygosity the segmenter is supposed to find.
     */
    static List<PileupSummary> sites(final String contig, final int count, final boolean loh) {
        final List<PileupSummary> sites = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            final int position = 1000 + 100 * i;
            final double alleleFrequency = 0.05 + (i % 19) * 0.05;
            final int depth = 50 + (i % 7) * 5;
            final int other = i % 4 == 0 ? 1 : 0;
            final int alt;
            switch (i % 5) {
                case 0:
                    alt = 1 + i % 2;                       // hom ref, with a little error
                    break;
                case 3:
                    alt = depth - other - (i % 2);         // hom alt
                    break;
                default:
                    alt = (loh && i >= count / 2) ? (depth - other) / 4 : (depth - other) / 2;
                    break;
            }
            final int ref = depth - alt - other;
            sites.add(new PileupSummary(contig, position, ref, alt, other, alleleFrequency));
        }
        return sites;
    }

    /** `CalculateContamination.filterSitesByCoverage`, reached by reflection so it is the reference's. */
    @SuppressWarnings("unchecked")
    static List<PileupSummary> filterSitesByCoverage(final List<PileupSummary> sites) throws Exception {
        final Class<?> tool = Class.forName(
                "org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination");
        final Object instance = tool.getDeclaredConstructor().newInstance();
        final Method method = tool.getDeclaredMethod("filterSitesByCoverage", List.class);
        method.setAccessible(true);
        return (List<PileupSummary>) method.invoke(instance, sites);
    }

    static void print(final String label, final List<PileupSummary> sites) {
        for (int i = 0; i < sites.size(); i++) {
            final PileupSummary site = sites.get(i);
            System.out.printf("sites\t%s\t%d=%s,%d,%d,%d,%d,%016x%n", label, i, site.getContig(),
                    site.getStart(), site.getRefCount(), site.getAltCount(), site.getOtherAltCount(),
                    Double.doubleToRawLongBits(site.getAlleleFrequency()));
        }
    }

    static void coverage(final String label, final List<PileupSummary> sites) {
        final StringBuilder text = new StringBuilder();
        for (final PileupSummary site : sites) {
            if (text.length() > 0) {
                text.append(',');
            }
            text.append(site.getContig()).append(':').append(site.getStart());
        }
        System.out.printf("coverage\t%s\t%s%n", label, text.length() == 0 ? "(none)" : text.toString());
    }

    static void segments(final String label, final List<List<PileupSummary>> segments) {
        for (int i = 0; i < segments.size(); i++) {
            final List<PileupSummary> segment = segments.get(i);
            System.out.printf("segments\t%s\t%d=%s,%d,%d,%d%n", label, i,
                    segment.get(0).getContig(), segment.get(0).getStart(),
                    segment.get(segment.size() - 1).getEnd(), segment.size());
        }
    }

    static void maf(final String label, final List<MinorAlleleFractionRecord> records) {
        for (int i = 0; i < records.size(); i++) {
            final MinorAlleleFractionRecord record = records.get(i);
            System.out.printf("maf\t%s\t%d=%s,%d,%d,%016x%n", label, i,
                    record.getContig(), record.getStart(), record.getEnd(),
                    Double.doubleToRawLongBits(record.getMinorAlleleFraction()));
        }
    }

    static void contamination(final String label, final Pair<Double, Double> answer) {
        System.out.printf("contamination\t%s\t%016x,%016x%n", label,
                Double.doubleToRawLongBits(answer.getLeft()),
                Double.doubleToRawLongBits(answer.getRight()));
    }
}
