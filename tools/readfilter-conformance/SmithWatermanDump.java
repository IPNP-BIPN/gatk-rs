/*
 * The Smith-Waterman aligners over the same reference and alternate pairs, taken from the
 * reference.
 *
 * `HaplotypeCaller` and `Mutect2` realign every haplotype to the reference, and reads to
 * haplotypes, with Smith-Waterman, and `--smith-waterman FASTEST_AVAILABLE` resolves to
 * `SmithWatermanIntelAligner` (Intel GKL, through JNI) whenever the native library loads. Nobody
 * had checked whether it returns the same offset and CIGAR as the pure-Java
 * `SmithWatermanJavaAligner` a port would target, and a golden taken without checking would pin
 * whichever one happened to run: the trap the PairHMM measurement found for the likelihood kernel.
 *
 * What this decides:
 *
 *   - IF THE IMPLEMENTATIONS AGREE on every pair, parameter set and overhang strategy, the choice
 *     is free and the port targets the readable one;
 *   - IF THEY DISAGREE, the oracle contract has to name one, and every dump that reaches the
 *     aligner has to force it with `--smith-waterman JAVA`.
 *
 * WHICH IMPLEMENTATIONS RAN IS PART OF THE ANSWER. The native library loads on some hosts and may
 * not load on others (the pinned image carries no `libgomp.so.1`, which stopped the vectorised
 * PairHMM), so the dump prints what it could build and what it could not, and a `loaded` row of
 * `no` is a result rather than a failure.
 *
 * The aligner is run over every (pair, SWParameters, SWOverhangStrategy) triple, and
 * `CigarUtils.calculateCigar` (pad with Ns, align, trim the padding, left-align) over every
 * pair with the two strategies GATK passes it, `SOFTCLIP` and `INDEL`.
 *
 * Output:
 *
 *     loaded\t<implementation>\t<yes|no:reason>
 *     fastest\t<simple class name FASTEST_AVAILABLE resolves to>
 *     align\t<pair>\t<params>\t<strategy>\t<implementation>\t<offset>\t<cigar|exception:name>
 *     alignagree\t<pair>\t<params>\t<strategy>\t<yes|no|single>
 *     cigar\t<pair>\t<params>\t<strategy>\t<implementation>\t<cigar|null|exception:name>
 *     cigaragree\t<pair>\t<params>\t<strategy>\t<yes|no|single>
 *
 * Usage: SmithWatermanDump
 */

import htsjdk.samtools.Cigar;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignment;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanIntelAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanJavaAligner;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class SmithWatermanDump {

    record Pair(String label, String ref, String alt) {
    }

    // A 40-base unique stretch, so that an edit in the middle has no ambiguous placement.
    static final String U = "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT";

    public static void main(final String[] args) throws Exception {
        System.out.println("# SmithWatermanDump: SmithWatermanJavaAligner and SmithWatermanIntelAligner");

        final List<Pair> pairs = List.of(
                new Pair("identical", U, U),
                new Pair("sub-middle", U, "GATTACAGCCTAGGCTTAAAGTCCAGTTGACCATGCAAGT"),
                new Pair("sub-first-base", U, "CATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT"),
                new Pair("sub-last-base", U, "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGA"),
                new Pair("two-adjacent-subs", U, "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGGTAGT"),
                new Pair("two-distant-subs", U, "GATTACTGCCTAGGCTTAACGTCCAGTTGACCATGCAAGC"),
                new Pair("insertion-middle", U, "GATTACAGCCTAGGCTTAACGTCCGGGAGTTGACCATGCAAGT"),
                new Pair("deletion-middle", U, "GATTACAGCCTAGGCTTAACGCAGTTGACCATGCAAGT"),
                new Pair("insertion-near-start", U, "GATTTTTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT"),
                new Pair("insertion-near-end", U, "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGTTTTGT"),
                new Pair("deletion-near-start", U, "GATAGCCTAGGCTTAACGTCCAGTTGACCATGCAAGT"),
                new Pair("deletion-near-end", U, "GATTACAGCCTAGGCTTAACGTCCAGTTGACCATGCA"),
                new Pair("insertion-and-deletion", U, "GATTACAGCCTAGGGGCTTAACGTCCAGTGACCATGCAAGT"),
                new Pair("sub-then-deletion", U, "GATTACAGCCTAGGCTTAACGTCGAGTGACCATGCAAGT"),
                // A gap in a homopolymer or a tandem repeat can sit at several positions: the
                // placement is the tie-break.
                new Pair("homopolymer-insertion", "CGTCAAAAAAAGTCGATCGTAGCTA", "CGTCAAAAAAAAGTCGATCGTAGCTA"),
                new Pair("homopolymer-deletion", "CGTCAAAAAAAGTCGATCGTAGCTA", "CGTCAAAAAAGTCGATCGTAGCTA"),
                new Pair("homopolymer-deletion-2", "CGTCAAAAAAAGTCGATCGTAGCTA", "CGTCAAAAAGTCGATCGTAGCTA"),
                new Pair("tandem-insertion", "TTGCACACACACACGGATCCAGTT", "TTGCACACACACACACGGATCCAGTT"),
                new Pair("tandem-deletion", "TTGCACACACACACGGATCCAGTT", "TTGCACACACGGATCCAGTT"),
                new Pair("tandem-trinucleotide-deletion", "ATGCAGCAGCAGCAGCAGTTGCA", "ATGCAGCAGCAGTTGCA"),
                new Pair("tandem-trinucleotide-insertion", "ATGCAGCAGCAGTTGCA", "ATGCAGCAGCAGCAGCAGTTGCA"),
                new Pair("homopolymer-at-start-deletion", "AAAAAGTCCGATTGCA", "AAAAGTCCGATTGCA"),
                new Pair("homopolymer-at-end-insertion", "GTCCGATTGCAAAAA", "GTCCGATTGCAAAAAA"),
                // The alternate is longer than the reference, or only partly inside it.
                new Pair("alt-longer-both-ends", "ACGTCCAGTTGACCAT", "TTTTTACGTCCAGTTGACCATGGGGG"),
                new Pair("alt-overhangs-left", "ACGTCCAGTTGACCAT", "GGGGGACGTCCAGTTGACCAT"),
                new Pair("alt-overhangs-right", "ACGTCCAGTTGACCAT", "ACGTCCAGTTGACCATGGGGG"),
                new Pair("alt-overhangs-left-mismatch", "ACGTCCAGTTGACCAT", "GGGGGACGTCCAGATGACCAT"),
                new Pair("alt-overhangs-right-mismatch", "ACGTCCAGTTGACCAT", "ACGTCCAGTTGACGATGGGGG"),
                // The alternate is a part of the reference: prefix, suffix, interior, and an
                // interior stretch that occurs twice (the exact-match shortcut takes the last).
                new Pair("alt-is-prefix", U, U.substring(0, 22)),
                new Pair("alt-is-suffix", U, U.substring(15)),
                new Pair("alt-is-interior", U, U.substring(8, 30)),
                new Pair("alt-occurs-twice", "TTGACCAGTCATTGACCAGTCATTTGAC", "TTGACCAGTCA"),
                new Pair("alt-occurs-twice-with-mismatch", "TTGACCAGTCATTGACCAGTCATTTGAC", "TTGACCAGACA"),
                // Short and degenerate shapes the Java accepts.
                new Pair("one-base-same", "A", "A"),
                new Pair("one-base-different", "A", "C"),
                new Pair("one-base-ref-long-alt", "A", "ACGTACGT"),
                new Pair("long-ref-one-base-alt", U, "T"),
                new Pair("two-bases-swapped", "ACGT", "AGCT"),
                new Pair("nothing-in-common", "AAAAAAAAAAAA", "CCCCCCCCCCCC"),
                new Pair("poly-n", "NNNNNNNNNN", "NNNNNNNN"),
                // The Java refuses these, and the refusal is part of the contract.
                new Pair("empty-reference", "", "ACGT"),
                new Pair("empty-alternate", "ACGT", ""),
                // Haplotype-shaped: a 60-base reference with a 5-base tandem deletion and an SNV
                // elsewhere, the shape `calculateCigar` is written for.
                new Pair("haplotype-deletion-and-snv",
                        "CTGAACGTTAGCCATGCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGATCGTAC",
                        "CTGAACGTTAGCCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGTTCGTAC"),
                new Pair("haplotype-insertion-and-snv",
                        "CTGAACGTTAGCCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGATCGTAC",
                        "CTGAACGTTAGCCATGCATGCATGCATGCAGTCAGGATCCAATTGGCTAGCTTAGGTTCGTAC"),
                new Pair("haplotype-leading-deletion", "AAAAAAGTCCGATTGCAGGCTAAC", "AAAAGTCCGATTGCAGGCTAAC"),
                new Pair("haplotype-trailing-deletion", "GTCCGATTGCAGGCTAACTTTTTT", "GTCCGATTGCAGGCTAACTTTT"));

        final Map<String, SWParameters> parameterSets = new LinkedHashMap<>();
        parameterSets.put("ORIGINAL_DEFAULT", SmithWatermanAlignmentConstants.ORIGINAL_DEFAULT);
        parameterSets.put("STANDARD_NGS", SmithWatermanAlignmentConstants.STANDARD_NGS);
        parameterSets.put("NEW_SW_PARAMETERS", SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS);
        parameterSets.put("ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS",
                SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS);

        // What could be built.
        final Map<String, SmithWatermanAligner> aligners = new LinkedHashMap<>();
        aligners.put("SmithWatermanJavaAligner", SmithWatermanJavaAligner.getInstance());
        System.out.println("loaded\tSmithWatermanJavaAligner\tyes");
        try {
            aligners.put("SmithWatermanIntelAligner", new SmithWatermanIntelAligner());
            System.out.println("loaded\tSmithWatermanIntelAligner\tyes");
        } catch (final Throwable failure) {
            System.out.printf("loaded\tSmithWatermanIntelAligner\tno:%s%n",
                    ReferenceQueryDump.escape(failure.getClass().getSimpleName()));
        }
        try {
            final SmithWatermanAligner fastest = SmithWatermanAligner.getAligner(
                    SmithWatermanAligner.Implementation.FASTEST_AVAILABLE);
            System.out.printf("fastest\t%s%n", fastest.getClass().getSimpleName());
        } catch (final Throwable failure) {
            System.out.printf("fastest\tno:%s%n",
                    ReferenceQueryDump.escape(failure.getClass().getSimpleName()));
        }

        for (final Pair pair : pairs) {
            final byte[] ref = pair.ref().getBytes(StandardCharsets.US_ASCII);
            final byte[] alt = pair.alt().getBytes(StandardCharsets.US_ASCII);
            for (final Map.Entry<String, SWParameters> params : parameterSets.entrySet()) {
                for (final SWOverhangStrategy strategy : SWOverhangStrategy.values()) {
                    final List<String> renderings = new ArrayList<>();
                    for (final Map.Entry<String, SmithWatermanAligner> aligner : aligners.entrySet()) {
                        String offset;
                        String cigar;
                        if (aligner.getKey().equals("SmithWatermanIntelAligner")
                                && (ref.length == 0 || alt.length == 0)) {
                            // An empty array goes straight to native code, where a refusal could be
                            // a crash of the whole JVM, which no handler here can survive and which
                            // would take every other row with it.
                            offset = "-";
                            cigar = "not-run:empty-input-to-native";
                            renderings.add(offset + "," + cigar);
                            System.out.printf("align\t%s\t%s\t%s\t%s\t%s\t%s%n", pair.label(),
                                    params.getKey(), strategy, aligner.getKey(), offset, cigar);
                            continue;
                        }
                        try {
                            final SmithWatermanAlignment alignment =
                                    aligner.getValue().align(ref, alt, params.getValue(), strategy);
                            offset = Integer.toString(alignment.getAlignmentOffset());
                            cigar = alignment.getCigar().toString();
                        } catch (final Throwable failure) {
                            offset = "-";
                            cigar = "exception:" + failure.getClass().getSimpleName();
                        }
                        renderings.add(offset + "," + cigar);
                        System.out.printf("align\t%s\t%s\t%s\t%s\t%s\t%s%n", pair.label(),
                                params.getKey(), strategy, aligner.getKey(), offset,
                                ReferenceQueryDump.escape(cigar));
                    }
                    System.out.printf("alignagree\t%s\t%s\t%s\t%s%n", pair.label(), params.getKey(),
                            strategy, agreement(renderings));
                }
                for (final SWOverhangStrategy strategy : new SWOverhangStrategy[] {
                        SWOverhangStrategy.SOFTCLIP, SWOverhangStrategy.INDEL}) {
                    final List<String> renderings = new ArrayList<>();
                    for (final Map.Entry<String, SmithWatermanAligner> aligner : aligners.entrySet()) {
                        String rendering;
                        try {
                            final Cigar cigar = CigarUtils.calculateCigar(ref, alt,
                                    aligner.getValue(), params.getValue(), strategy);
                            rendering = cigar == null ? "null" : cigar.toString();
                        } catch (final Throwable failure) {
                            rendering = "exception:" + failure.getClass().getSimpleName();
                        }
                        renderings.add(rendering);
                        System.out.printf("cigar\t%s\t%s\t%s\t%s\t%s%n", pair.label(),
                                params.getKey(), strategy, aligner.getKey(),
                                ReferenceQueryDump.escape(rendering));
                    }
                    System.out.printf("cigaragree\t%s\t%s\t%s\t%s%n", pair.label(), params.getKey(),
                            strategy, agreement(renderings));
                }
            }
        }

        for (final SmithWatermanAligner aligner : aligners.values()) {
            aligner.close();
        }
    }

    static String agreement(final List<String> renderings) {
        if (renderings.size() == 1) {
            return "single";
        }
        return renderings.stream().distinct().count() == 1 ? "yes" : "no";
    }
}
