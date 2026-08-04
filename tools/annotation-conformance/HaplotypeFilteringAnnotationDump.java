/*
 * HaplotypeFilteringAnnotation, taken from the reference.
 *
 * Two counts, ASSEMBLED_HAPS and FILTERED_HAPS, read off a likelihood matrix whose allele axis is
 * a Haplotype rather than an Allele. The arithmetic is nothing; what the cases below pin is what
 * the two numbers actually count:
 *
 *   - ASSEMBLED_HAPS is alleles().size(), which is the size of an IndexedSet under
 *     Haplotype.equals. That equality is the uniqueness value, the reference flag and the bases,
 *     while Haplotype.hashCode is the bases alone. So two haplotypes with identical bases are one
 *     entry or two depending on a field that does not appear in the output, and the cases
 *     duplicate-bases-same-uniqueness and duplicate-bases-different-uniqueness are the pair that
 *     shows it;
 *   - FILTERED_HAPS is a getter over an int field that only AlleleFiltering ever writes, so a
 *     matrix nobody filtered reports 0 rather than no key at all. unfiltered is that case, and
 *     negative-filtered-count shows the getter does not clamp;
 *   - the engine hands this annotation either an AlleleLikelihoods<Fragment, Haplotype> or an
 *     AlleleLikelihoods<GATKRead, Haplotype>, chosen by a ternary in
 *     VariantAnnotatorEngine.addInfoAnnotations. by-fragment and by-read are the same matrix
 *     through both branches.
 *
 * Output:
 *
 *     anno\t<label>\t<key>=<value>[<class>];...
 *     keys\t<key>,<key>
 *
 * Usage: HaplotypeFilteringAnnotationDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import org.broadinstitute.hellbender.tools.walkers.annotator.HaplotypeFilteringAnnotation;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedAlleleList;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.Fragment;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class HaplotypeFilteringAnnotationDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT = Allele.create("C", false);
    static final SAMFileHeader HEADER = makeHeader();
    static final int VARIANT_START = 105;

    public static void main(final String[] args) {
        System.out.println("# HaplotypeFilteringAnnotationDump: ASSEMBLED_HAPS and FILTERED_HAPS");

        // The keys, in getKeyNames order, which is the declaration order and not the map's.
        final StringJoiner keys = new StringJoiner(",");
        new HaplotypeFilteringAnnotation().getKeyNames().forEach(keys::add);
        System.out.printf("keys\t%s%n", keys);

        // The count of haplotypes, over lists an IndexedSet treats differently.
        emit("no-haplotypes", byRead(List.of(), 0));
        emit("one-haplotype", byRead(List.of(hap("ACGT", false, 0)), 0));
        emit("three-haplotypes", byRead(List.of(
                hap("ACGT", true, 0), hap("ACGA", false, 0), hap("ACGC", false, 0)), 0));
        emit("duplicate-bases-same-uniqueness", byRead(List.of(
                hap("ACGT", false, 0), hap("ACGT", false, 0)), 0));
        emit("duplicate-bases-different-uniqueness", byRead(List.of(
                hap("ACGT", false, 0), hap("ACGT", false, 1)), 0));
        emit("duplicate-bases-different-ref-flag", byRead(List.of(
                hap("ACGT", true, 0), hap("ACGT", false, 0)), 0));
        emit("different-lengths", byRead(List.of(
                hap("ACGT", false, 0), hap("ACGTACGT", false, 0)), 0));

        // The filtered count, which is carried rather than computed.
        emit("unfiltered", byRead(List.of(
                hap("ACGT", true, 0), hap("ACGA", false, 0), hap("ACGC", false, 0)), 0));
        emit("two-filtered", byRead(List.of(
                hap("ACGT", true, 0), hap("ACGA", false, 0), hap("ACGC", false, 0)), 2));
        emit("filtered-exceeds-remaining", byRead(List.of(hap("ACGT", true, 0)), 5));
        emit("negative-filtered-count", byRead(List.of(
                hap("ACGT", true, 0), hap("ACGA", false, 0)), -1));

        // The two instantiations the engine's ternary chooses between, over the same haplotypes.
        final List<Haplotype> haplotypes = List.of(
                hap("ACGT", true, 0), hap("ACGA", false, 0), hap("ACGC", false, 0));
        emit("by-read", byRead(haplotypes, 1));
        emit("by-fragment", byFragment(haplotypes, 1));

        // No evidence at all, with haplotypes present: the counts do not depend on the evidence.
        emit("no-evidence", emptyEvidence(haplotypes, 1));
    }

    static SAMFileHeader makeHeader() {
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1000)));
        final SAMFileHeader header = new SAMFileHeader(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    static Haplotype hap(final String bases, final boolean isRef, final int uniqueness) {
        final Haplotype haplotype = new Haplotype(bases.getBytes(), isRef);
        haplotype.setUniquenessValue(uniqueness);
        return haplotype;
    }

    static VariantContext site() {
        return new VariantContextBuilder().chr("chr1").start(VARIANT_START).stop(VARIANT_START)
                .alleles(List.of(REF, ALT)).make();
    }

    static GATKRead read(final String name, final int start) {
        final SAMRecord record = new SAMRecord(HEADER);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(start);
        record.setCigarString("20M");
        final byte[] bases = new byte[20];
        Arrays.fill(bases, (byte) 'A');
        record.setReadBases(bases);
        final byte[] qualities = new byte[20];
        Arrays.fill(qualities, (byte) 30);
        record.setBaseQualities(qualities);
        record.setMappingQuality(60);
        record.setInferredInsertSize(300);
        return new SAMRecordToGATKReadAdapter(record);
    }

    /** AlleleLikelihoods<GATKRead, Haplotype>, the branch taken when the fragment one is absent. */
    static AlleleLikelihoods<GATKRead, Haplotype> byRead(final List<Haplotype> haplotypes,
                                                        final int filtered) {
        final List<GATKRead> reads = List.of(read("r0", 100), read("r1", 101));
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(reads));
        final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(haplotypes), bySample);
        fill(likelihoods, reads.size());
        likelihoods.setFilteredHaplotypeCount(filtered);
        return likelihoods;
    }

    /** AlleleLikelihoods<Fragment, Haplotype>, the branch preferred when it is present. */
    static AlleleLikelihoods<Fragment, Haplotype> byFragment(final List<Haplotype> haplotypes,
                                                             final int filtered) {
        final List<Fragment> fragments = List.of(
                Fragment.createAndAvoidFailure(List.of(read("r0", 100))),
                Fragment.createAndAvoidFailure(List.of(read("r1", 101))));
        final Map<String, List<Fragment>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>(fragments));
        final AlleleLikelihoods<Fragment, Haplotype> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(haplotypes), bySample);
        fill(likelihoods, fragments.size());
        likelihoods.setFilteredHaplotypeCount(filtered);
        return likelihoods;
    }

    static AlleleLikelihoods<GATKRead, Haplotype> emptyEvidence(final List<Haplotype> haplotypes,
                                                                final int filtered) {
        final Map<String, List<GATKRead>> bySample = new LinkedHashMap<>();
        bySample.put("s1", new ArrayList<>());
        final AlleleLikelihoods<GATKRead, Haplotype> likelihoods = new AlleleLikelihoods<>(
                new IndexedSampleList(List.of("s1")), new IndexedAlleleList<>(haplotypes), bySample);
        likelihoods.setFilteredHaplotypeCount(filtered);
        return likelihoods;
    }

    static <E extends htsjdk.samtools.util.Locatable> void fill(
            final AlleleLikelihoods<E, Haplotype> likelihoods, final int evidenceCount) {
        if (likelihoods.numberOfAlleles() == 0) {
            return;
        }
        final LikelihoodMatrix<E, Haplotype> m = likelihoods.sampleMatrix(0);
        for (int a = 0; a < likelihoods.numberOfAlleles(); a++) {
            for (int e = 0; e < evidenceCount; e++) {
                m.set(a, e, a == 0 ? -1 : -10);
            }
        }
    }

    static void emit(final String label,
                     final AlleleLikelihoods<? extends htsjdk.samtools.util.Locatable, Haplotype> haplotypeLikelihoods) {
        try {
            // The first two arguments are the reference and feature contexts, and the read and
            // fragment matrices; this annotation reads none of them, and the engine is entitled to
            // pass null for the fragment one, so nulls are what the dump passes.
            final Map<String, Object> result = new HaplotypeFilteringAnnotation()
                    .annotate(null, null, site(), null, null, haplotypeLikelihoods);
            final StringJoiner joiner = new StringJoiner(";");
            // Sorted, because the reference builds a HashMap and its iteration order over these two
            // keys is String.hashCode's. The engine copies the result into a LinkedHashMap the
            // encoder then sorts, so nothing observable depends on that order and a dump that
            // reported it would be pinning the JDK's hash rather than this annotation.
            result.entrySet().stream()
                    .sorted(Map.Entry.comparingByKey())
                    .forEach(entry -> {
                        final Object value = entry.getValue();
                        joiner.add(String.format("%s=%s[%s]", entry.getKey(), value,
                                value == null ? "null" : value.getClass().getName()));
                    });
            System.out.printf("anno\t%s\t%s%n", label, joiner);
        } catch (final Exception | AssertionError e) {
            System.out.printf("anno\t%s\tE:%s%n", label, e.getClass().getName());
        }
    }
}
