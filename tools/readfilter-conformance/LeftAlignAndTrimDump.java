/*
 * GATKVariantContextUtils.leftAlignAndTrim, taken from the reference.
 *
 * The workhorse under LeftAlignAndTrimVariants, and the third brick of the variant-transform
 * archetype. It slides an indel as far left as the reference lets it, widening its own window
 * until the shift stops hitting the edge.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE WINDOW DOUBLES, `for (leadingBases = min(maxLeadingBases, 10); ...; leadingBases =
 *     min(2 * leadingBases, maxLeadingBases))`, so an indel that can move more than 10 bases costs
 *     several passes over widening slices of reference;
 *   - AND THE LOOP ONLY CONTINUES WHEN THE SHIFT REACHED THE EDGE, `shifts.getLeft() ==
 *     variantOffsetInRef`, so an alignment that stops one base short of the window is FINAL even
 *     though a wider window was allowed;
 *   - A SHIFT OF ZERO RETURNS THE SAME OBJECT, before any allele is rebuilt;
 *   - maxLeadingBases <= 0 RETURNS AT ONCE, which is what the caller's
 *     `min(maxLeadingBases, distanceToLastVariant - 1)` produces for two adjacent variants;
 *   - AND A NON-INDEL RETURNS AT ONCE TOO, `isIndel()` being the type test, so a MIXED site is
 *     never aligned;
 *   - THE NEW ALLELES COME OUT OF A HashMap, `alleleMap.values()`, so their order is the map's and
 *     not the record's. On the three-allele case measured here the two happen to agree, which is
 *     what the golden pins: a port that rebuilds the list in record order matches this file, and
 *     nothing here says it would match every input;
 *   - A GENOTYPE ALLELE THAT IS NOT IN THE MAP WOULD BECOME NO_CALL, and no valid record can carry
 *     one: htsjdk refuses to build a VariantContext whose genotype calls an allele it does not
 *     list, "Allele in genotype * not in the variant context", so that branch is unreachable from
 *     anything this tool can be given;
 *   - THE START AND THE END MOVE BY DIFFERENT AMOUNTS, `start - shifts.getLeft()` and
 *     `stop - shifts.getRight()`, and normalizeAlleles can shift RIGHT when trimming is on;
 *   - AND WITH trim OFF THE ALLELES KEEP THEIR LENGTH, so the record moves without shrinking.
 *
 * Output:
 *
 *     reference\t<the reference bases>
 *     in\t<label>\t<contig>:<start>-<end>\t<alleles>\t<genotypes>
 *     out\t<label>\t<contig>:<start>-<end>\t<alleles>\t<genotypes>
 *     same\t<label>\t<true when the very same object came back>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: LeftAlignAndTrimDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.ReferenceFileSource;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LeftAlignAndTrimDump {

    /**
     * A reference with two homopolymers of different lengths and one dinucleotide repeat, so an
     * indel can walk a little, a lot, or not at all.
     *
     * Positions are one-based: the A run is 11..18, the T run 31..60, and the CA repeat 71..90.
     */
    static final String REFERENCE =
            "GGGGGGGGGG"                      // 1..10
            + "AAAAAAAA"                      // 11..18
            + "GGGGGGGGGGGG"                  // 19..30
            + "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" // 31..60
            + "GGGGGGGGGG"                    // 61..70
            + "CACACACACACACACACACA"          // 71..90
            + "GGGGGGGGGG";                   // 91..100

    static ReferenceDataSource reference;

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("leftalignandtrim-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# LeftAlignAndTrimDump: leftAlignAndTrim, from the reference");
        System.out.printf("reference\t%s%n", REFERENCE);

        final Path fasta = writeReference(dir);
        try (final ReferenceDataSource source = new ReferenceFileSource(fasta)) {
            reference = source;

            // A deletion of one A at the right end of the A run, which can walk to the run's start.
            run("deletion-in-short-run", deletion(18, 1), 1000, true);
            // The same deletion with only two bases of window, which stops at the edge and the loop
            // gives up when the window is already the maximum.
            run("deletion-narrow-window", deletion(18, 1), 2, true);
            // And with a window of exactly the distance it wants to move.
            run("deletion-exact-window", deletion(18, 1), 7, true);

            // A deletion in the long T run, which needs the doubling loop: 10, then 20, then 30.
            run("deletion-in-long-run", deletion(60, 1), 1000, true);
            run("deletion-in-long-run-narrow", deletion(60, 1), 10, true);
            run("deletion-in-long-run-twenty", deletion(60, 1), 20, true);

            // An insertion of one T at the end of the same run.
            run("insertion-in-long-run", insertion(60, "T"), 1000, true);
            // An insertion of the repeat unit in the CA repeat, which walks by twos.
            run("insertion-in-repeat", insertion(89, "AC"), 1000, true);

            // Already left aligned: the FIRST base of the A run, anchored on the G before it.
            run("already-aligned", deletion(11, 1), 1000, true);
            // A deleted base that differs from the one before it, so nothing can move.
            run("no-repeat", deletion(71, 1), 1000, true);

            // The two early exits.
            run("zero-window", deletion(18, 1), 0, true);
            run("negative-window", deletion(18, 1), -5, true);
            run("snp", snp(18), 1000, true);
            run("mixed", mixed(18), 1000, true);

            // Trimming off, so the alleles keep their length while the record moves.
            run("no-trim", deletion(18, 1), 1000, false);
            run("no-trim-untrimmed-alleles", untrimmed(18), 1000, false);
            run("trim-untrimmed-alleles", untrimmed(18), 1000, true);

            // Genotypes, which are remapped through the same map the alleles came out of.
            run("with-genotypes", withGenotypes(deletion(18, 1)), 1000, true);

            // Three alleles, which the doc says cannot happen and the code does not check: the
            // order the alleles come back in is the hash order of the map they went through.
            run("multiallelic", multiallelic(18), 1000, true);
        }
    }

    /**
     * A deletion of the `length` bases starting AT `position`, anchored on the base before it, so
     * the record starts one base earlier than what it removes.
     */
    static VariantContext deletion(final int position, final int length) {
        final int anchor = position - 1;
        final String reference = REFERENCE.substring(anchor - 1, position - 1 + length);
        final String alternate = REFERENCE.substring(anchor - 1, anchor);
        return new VariantContextBuilder("x", "chr1", anchor, anchor + length,
                List.of(Allele.create(reference, true), Allele.create(alternate, false))).make();
    }

    /** An insertion of `bases` after `position`. */
    static VariantContext insertion(final int position, final String bases) {
        final String reference = REFERENCE.substring(position - 1, position);
        return new VariantContextBuilder("x", "chr1", position, position,
                List.of(Allele.create(reference, true), Allele.create(reference + bases, false)))
                .make();
    }

    static VariantContext snp(final int position) {
        final String reference = REFERENCE.substring(position - 1, position);
        final String alternate = reference.equals("A") ? "C" : "A";
        return new VariantContextBuilder("x", "chr1", position, position,
                List.of(Allele.create(reference, true), Allele.create(alternate, false))).make();
    }

    /** A snp and a deletion at one site, which is MIXED and therefore not an indel. */
    static VariantContext mixed(final int position) {
        final String reference = REFERENCE.substring(position - 1, position + 1);
        return new VariantContextBuilder("x", "chr1", position, position + 1,
                List.of(Allele.create(reference, true),
                        Allele.create(REFERENCE.substring(position - 1, position), false),
                        Allele.create("CC", false))).make();
    }

    /** A deletion whose alleles share a trailing base as well as the leading one. */
    static VariantContext untrimmed(final int position) {
        final String reference = REFERENCE.substring(position - 1, position + 2);
        final String alternate = REFERENCE.substring(position - 1, position)
                + REFERENCE.substring(position + 1, position + 2);
        return new VariantContextBuilder("x", "chr1", position, position + 2,
                List.of(Allele.create(reference, true), Allele.create(alternate, false))).make();
    }

    /**
     * A deletion and an insertion of the same base at one site, both of which can walk left, so the
     * order the alleles come back in is the order the map hands them over.
     */
    static VariantContext multiallelic(final int position) {
        final int anchor = position - 1;
        final String reference = REFERENCE.substring(anchor - 1, position);
        final String deleted = REFERENCE.substring(anchor - 1, anchor);
        return new VariantContextBuilder("x", "chr1", anchor, position,
                List.of(Allele.create(reference, true),
                        Allele.create(deleted, false),
                        Allele.create(reference + REFERENCE.substring(anchor - 1, anchor), false)))
                .make();
    }

    static VariantContext withGenotypes(final VariantContext variant) {
        final List<Genotype> genotypes = new ArrayList<>();
        genotypes.add(new GenotypeBuilder("s1",
                List.of(variant.getReference(), variant.getAlternateAllele(0))).make());
        genotypes.add(new GenotypeBuilder("s2",
                List.of(variant.getAlternateAllele(0), variant.getAlternateAllele(0))).make());
        return new VariantContextBuilder(variant).genotypes(genotypes).make();
    }

    static void run(final String label, final VariantContext variant, final int maxLeadingBases,
                    final boolean trim) {
        System.out.printf("in\t%s\t%s\t%s\t%s%n", label, place(variant), alleles(variant),
                genotypes(variant));
        final VariantContext result;
        try {
            final ReferenceContext context = new ReferenceContext(reference,
                    new SimpleInterval(variant.getContig(), variant.getStart(), variant.getEnd()));
            result = GATKVariantContextUtils.leftAlignAndTrim(variant, context, maxLeadingBases, trim);
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("out\t%s\t%s\t%s\t%s%n", label, place(result), alleles(result),
                genotypes(result));
        System.out.printf("same\t%s\t%s%n", label, result == variant);
    }

    static String place(final VariantContext variant) {
        return variant.getContig() + ":" + variant.getStart() + "-" + variant.getEnd();
    }

    /** The alleles in the order the record carries them, the reference marked with a `*`. */
    static String alleles(final VariantContext variant) {
        final List<String> out = new ArrayList<>();
        for (final Allele allele : variant.getAlleles()) {
            out.add(allele.getBaseString() + (allele.isReference() ? "*" : ""));
        }
        return String.join(",", out);
    }

    static String genotypes(final VariantContext variant) {
        final List<String> out = new ArrayList<>();
        for (final Genotype genotype : variant.getGenotypes()) {
            final List<String> called = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                called.add(allele.getDisplayString());
            }
            out.add(genotype.getSampleName() + "=" + String.join("/", called));
        }
        return String.join(";", out);
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        Files.writeString(fasta, ">chr1\n" + REFERENCE + "\n", StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(fasta, true);
        final Path dict = dir.resolve("reference.dict");
        Files.writeString(dict, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:" + REFERENCE.length()
                + "\tM5:0\tUR:file:" + fasta + "\n", StandardCharsets.UTF_8);
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

    static {
        // Keeps the unused import warning away from Arrays, which the allele helpers use.
        Arrays.hashCode(new int[0]);
    }
}
