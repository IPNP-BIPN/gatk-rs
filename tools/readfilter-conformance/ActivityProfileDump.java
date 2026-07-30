/*
 * Where an assembly region starts and stops, taken from the reference.
 *
 * This is what HaplotypeCaller and Mutect2 assemble over: a boundary that moves by one base changes
 * the haplotypes and therefore the calls, so it is upstream of every variant those tools emit.
 *
 * Four behaviours decide the boundaries and none follows from "cut where the activity changes".
 *
 *   - a probability is spread over a Gaussian and the spread is ADDED, not assigned.
 *     BandPassActivityProfile.processState turns one state into 2 * filterSize + 1 states carrying
 *     prob * kernel[i], and incorporateSingleState adds each into whatever is already there. So the
 *     probability at a site is the sum of its neighbours' tails, and a site never reported active
 *     can end up above the threshold. A probability of exactly 0.0 skips the filter entirely;
 *   - the filter size is chosen from the kernel's VALUES, not from sigma: determineFilterSize walks
 *     in from the edge while the values are at least 1e-5. The width then feeds back into
 *     getMaxProbPropagationDistance, which decides when a region may be popped at all;
 *   - the cut site is a local minimum searched backwards, and isMinimum is asymmetric:
 *     p[i] <= p[i+1] but p[i] < p[i-1]. On a plateau that picks the left-hand end;
 *   - nothing can be popped until the profile holds maxRegionSize + maxProbPropagationDistance
 *     states, so regions appear in bursts. forceConversion bypasses that and first trims every
 *     state the filter wrote past the last position added.
 *
 * The kernels are dumped as raw bits, because the whole point is the last bits: Math.exp is the
 * platform's and a port's exp may differ in the final ulp, which would move a cut site.
 *
 * Output:
 *
 *     kernel\t<maxFilterSize>\t<sigma>\t<filterSize>\t<bandSize>\t<raw bits, comma-separated>
 *     probs\t<label>\t<raw bits of every state probability, comma-separated>
 *     region\t<label>\t<n>\t<contig>:<start>-<end>\t<isActive>
 *     summary\t<label>\t<regions>\t<states left>\t<maxProbPropagationDistance>
 *     error\t<label>\t<class>\t<message>
 *
 * Usage: ActivityProfileDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.activityprofile.BandPassActivityProfile;

import java.lang.reflect.Method;
import java.util.List;
import java.util.StringJoiner;

public class ActivityProfileDump {

    static final int CONTIG_LENGTH = 10000;

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        return header;
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# ActivityProfileDump: where an assembly region starts and stops");

        // The kernels, at the sizes and sigmas the callers use and either side of them.
        for (final int maxFilterSize : new int[] {0, 1, 5, 25, 50}) {
            for (final double sigma : new double[] {0.5, 1.0, 2.0, 17.0}) {
                kernel(maxFilterSize, sigma);
            }
        }

        // A plain profile: no filter, so the probabilities are what was added.
        profile("plain-active-then-inactive", null, 0.002, 50,
                probs(20, 0.9), probs(20, 0.0));
        // The same input through the band pass filter, where the spread crosses the threshold on
        // sites that were added as zero.
        profile("bandpass-active-then-inactive", 17.0, 0.002, 50,
                probs(20, 0.9), probs(20, 0.0));
        // An active stretch longer than maxRegionSize, which forces a cut at a local minimum.
        profile("long-active", null, 0.002, 10,
                probs(5, 0.9), probs(1, 0.1), probs(5, 0.9), probs(1, 0.05), probs(30, 0.9));
        // A plateau at the minimum, where the asymmetric comparison decides which end is cut.
        profile("plateau", null, 0.002, 10,
                probs(4, 0.9), probs(4, 0.1), probs(30, 0.9));
        // Everything inactive, so the first region is an inactive one.
        profile("all-inactive", null, 0.002, 20, probs(60, 0.0));
        // A profile too short to pop anything without forcing.
        profile("too-short", null, 0.002, 50, probs(5, 0.9));
        // The zero-probability case, which skips the filter.
        profile("bandpass-zeros", 17.0, 0.002, 20, probs(40, 0.0));
        // A single state, forced.
        profile("single-state", null, 0.002, 20, probs(1, 0.9));
    }

    static double[] probs(final int count, final double value) {
        final double[] values = new double[count];
        java.util.Arrays.fill(values, value);
        return values;
    }

    static void kernel(final int maxFilterSize, final double sigma) throws Exception {
        try {
            final BandPassActivityProfile profile = new BandPassActivityProfile(
                    50, 0.002, maxFilterSize, sigma, header());
            final Method getKernel = BandPassActivityProfile.class.getDeclaredMethod("getKernel");
            getKernel.setAccessible(true);
            final double[] kernel = (double[]) getKernel.invoke(profile);
            final StringJoiner bits = new StringJoiner(",");
            for (final double value : kernel) {
                bits.add(Long.toString(Double.doubleToRawLongBits(value)));
            }
            System.out.printf("kernel\t%d\t%s\t%d\t%d\t%s%n", maxFilterSize, sigma,
                    profile.getFilteredSize(), profile.getBandSize(), bits);
        } catch (final Throwable t) {
            System.out.printf("error\tkernel-%d-%s\t%s\t%s%n", maxFilterSize, sigma,
                    t.getClass().getName(), oneLine(t.getMessage()));
        }
    }

    /**
     * Feed the probabilities in one at a time from chr1:100, popping after every add and once more
     * with forceConversion at the end. Popping as we go is what the walker does, and it matters:
     * the same probabilities popped only at the end give different regions.
     */
    static void profile(final String label, final Double sigma, final double threshold,
                        final int maxRegionSize, final double[]... blocks) {
        try {
            final ActivityProfile profile = sigma == null
                    ? new ActivityProfile(50, threshold, header())
                    : new BandPassActivityProfile(50, threshold, 50, sigma, header());

            int position = 100;
            int index = 0;
            for (final double[] block : blocks) {
                for (final double value : block) {
                    profile.add(new ActivityProfileState(
                            new SimpleInterval("chr1", position, position), value));
                    position++;
                    for (final AssemblyRegion region
                            : profile.popReadyAssemblyRegions(0, 1, maxRegionSize, false)) {
                        System.out.printf("region\t%s\t%d\t%s\t%b%n", label, index++,
                                region.getSpan(), region.isActive());
                    }
                }
            }

            // The probabilities left in the profile, before the flush: this is where a divergent
            // kernel shows up as a number rather than as a moved boundary.
            final StringJoiner bits = new StringJoiner(",");
            for (final double value : probabilities(profile)) {
                bits.add(Long.toString(Double.doubleToRawLongBits(value)));
            }
            System.out.printf("probs\t%s\t%s%n", label, bits);

            for (final AssemblyRegion region
                    : profile.popReadyAssemblyRegions(0, 1, maxRegionSize, true)) {
                System.out.printf("region\t%s\t%d\t%s\t%b%n", label, index++,
                        region.getSpan(), region.isActive());
            }

            System.out.printf("summary\t%s\t%d\t%d\t%d%n", label, index, profile.size(),
                    profile.getMaxProbPropagationDistance());
        } catch (final Throwable t) {
            System.out.printf("error\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static double[] probabilities(final ActivityProfile profile) throws Exception {
        final Method method = ActivityProfile.class.getDeclaredMethod("getProbabilitiesAsArray");
        method.setAccessible(true);
        return (double[]) method.invoke(profile);
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
