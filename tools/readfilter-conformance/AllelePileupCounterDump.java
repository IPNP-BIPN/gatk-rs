/*
 * What a validation pileup says about an allele, taken from the reference.
 *
 * `ValidateBasicSomaticShortMutations` reads a discovery call, counts the validation pileups under
 * it and decides whether the call is out of the noise floor. `SomaticValidationPowerDump` measured
 * the arithmetic of that decision; this measures everything that feeds it, which is a chain of
 * allele choices over pileup elements.
 *
 * Nine behaviours this is built to catch.
 *
 *   - `typeOfVariant` CALLS EQUAL-LENGTH ALLELES DIFFERING IN ONE BASE A SNP, not an MNP, however
 *     long they are, and calls a spanning-deletion alternate NO_VARIATION;
 *   - `isComplexIndel` IS A PREFIX TEST, so a deletion whose remaining base disagrees with the
 *     reference's first is complex and a genotype carrying it cannot be validated at all;
 *   - `chooseAlleleForRead` DOES NOT STOP AT THE FIRST MATCH. The loop assigns without breaking,
 *     so the LAST matching alternate wins;
 *   - THE REFERENCE TEST INCLUDES THE INDEL LOOKAHEAD. A read whose bases match the reference is
 *     still not reference if it is before an insertion or a deletion start, which is what lets an
 *     indel-carrying read be counted as its alternate;
 *   - A READ THAT ENDS INSIDE THE ALLELE ANSWERS UNKNOWN rather than false, and `basesMatch`
 *     against the short array is false, so such a read is neither reference nor alternate;
 *   - THE QUALITY CUTOFF IS APPLIED LAST, to the chosen allele's own bases, and the minimum over
 *     an empty array is -1, so a read that ends exactly at the allele is refused by quality;
 *   - THE DELETION MATCH COMPARES LENGTHS ONLY: any alternate whose deletion length matches is the
 *     pileup's allele regardless of its bases;
 *   - `AllelePileupCounter` DROPS MAPPING QUALITY 0 AND 255, and drops nothing else: a base under
 *     the cutoff is chosen as no allele rather than filtered out of the pileup;
 *   - A HAPLOID GENOTYPE THROWS rather than being refused, because `isAbleToValidateGenotype`
 *     types the variant from the second allele before the ploidy it just computed is consulted;
 *   - AND `calculateMaxAltRatio` COUNTS DELETIONS OUT BUT INDEL LOOKAHEAD IN. Its two filters are
 *     not complements: an element that is UNKNOWN is in neither the numerator nor the denominator,
 *     so a pileup of nothing but short reads is a ratio of zero rather than a NaN.
 *
 * Output:
 *
 *     type\t<ref>,<alt>=<type>,<complex>
 *     choose\t<case>\t<read>=<allele or none>
 *     count\t<case>\t<allele>=<count>...
 *     ratio\t<case>,<minqual>=<bits>
 *     support\t<case>,<alt>,<minqual>=<count>
 *     validatable\t<case>=<true|false>
 *     result\t<case>=<none> or the sixteen fields
 *     table\t<line>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: AllelePileupCounterDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import org.apache.commons.lang3.mutable.MutableInt;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.AllelePileupCounter;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BasicSomaticShortMutationValidator;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BasicValidationResult;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.PowerCalculationUtils;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class AllelePileupCounterDump {

    static final int CONTIG_LENGTH = 200;
    /** The locus every pileup is taken at, five bases into each read. */
    static final int LOCUS = 105;
    static final int START = 100;

    /** A read: name, cigar, bases, base qualities, mapping quality. */
    static String[] read(final String name, final String cigar, final String bases,
                         final String quals, final String mapq) {
        return new String[] {name, cigar, bases, quals, mapq};
    }

    /** Every pileup this dump measures, by case name. */
    static final Map<String, String[][]> CASES = Map.ofEntries(
            // A plain SNP: two reference reads, two alternate, one alternate under any cutoff.
            Map.entry("snp", new String[][] {
                read("ref1", "10M", "AAAAAAAAAA", "30", "60"),
                read("ref2", "10M", "AAAAAAAAAA", "30", "60"),
                read("alt1", "10M", "AAAAACAAAA", "30", "60"),
                read("alt2", "10M", "AAAAACAAAA", "30", "60"),
                read("altq", "10M", "AAAAACAAAA", "5", "60"),
            }),
            // Mapping quality 0 and the unavailable 255, which the counter drops and the ratio
            // does not.
            Map.entry("mapq", new String[][] {
                read("ref1", "10M", "AAAAAAAAAA", "30", "60"),
                read("zero", "10M", "AAAAACAAAA", "30", "0"),
                read("none", "10M", "AAAAACAAAA", "30", "255"),
            }),
            // An insertion of three bases immediately after the locus, and a read that carries a
            // deletion in the same place instead.
            Map.entry("insertion", new String[][] {
                read("ref1", "10M", "AAAAAAAAAA", "30", "60"),
                read("ins1", "6M3I4M", "AAAAAATTTAAAA", "30", "60"),
                read("ins2", "6M3I4M", "AAAAAATTTAAAA", "30", "60"),
                read("del1", "6M2D4M", "AAAAAAAAAA", "30", "60"),
            }),
            // A deletion of two bases starting after the locus, plus the deleted position itself
            // seen from a read that spans it.
            Map.entry("deletion", new String[][] {
                read("ref1", "10M", "AAAAAAAAAA", "30", "60"),
                read("del1", "6M2D4M", "AAAAAAAAAA", "30", "60"),
                read("del2", "6M2D4M", "AAAAAAAAAA", "30", "60"),
                read("del3", "6M5D4M", "AAAAAAAAAA", "30", "60"),
            }),
            // Reads that end inside a two-base allele, which is the UNKNOWN answer.
            Map.entry("short", new String[][] {
                read("stop", "6M", "AAAAAA", "30", "60"),
                read("stop2", "6M", "AAAAAC", "30", "60"),
            }),
            // Nothing at the locus at all.
            Map.entry("empty", new String[][] {}));

    public static void main(final String[] args) throws Exception {
        System.out.println("# AllelePileupCounterDump: what a validation pileup says about an allele");
        final SAMFileHeader header = header();

        final Allele refA = Allele.create("A", true);
        final Allele refAA = Allele.create("AA", true);
        final Allele refAAA = Allele.create("AAA", true);
        final Allele altC = Allele.create("C", false);
        final Allele altG = Allele.create("G", false);
        final Allele altGT = Allele.create("GT", false);
        final Allele insertion = Allele.create("ATTT", false);
        final Allele deletion = Allele.create("A", false);
        final Allele other = Allele.create("T", false);

        // The type of every pairing worth a row, including the two the reference refuses.
        final String[][] pairs = {
                {"A", "C"}, {"A", "A"}, {"AC", "AG"}, {"AC", "GT"}, {"ACGT", "ACGA"},
                {"ACGT", "AGGA"}, {"A", "ATTT"}, {"ATTT", "A"}, {"AAA", "A"}, {"AAA", "AA"},
                {"AAA", "TA"}, {"AAAA", "TTT"}, {"A", "*"}, {"A", "<NON_REF>"}, {"AC", "A"},
        };
        for (final String[] pair : pairs) {
            final Allele reference = Allele.create(pair[0], true);
            final Allele alternate = Allele.create(pair[1], false);
            System.out.printf("type\t%s,%s=%s,%s%n", pair[0], pair[1],
                    GATKVariantContextUtils.typeOfVariant(reference, alternate),
                    GATKVariantContextUtils.isComplexIndel(reference, alternate));
        }

        // The chosen allele, read by read, over the alternates each case is about.
        chooseCase(header, "snp", refA, List.of(altC), 0);
        chooseCase(header, "snp-q20", refA, List.of(altC), 20);
        chooseCase(header, "snp-two", refA, List.of(altC, altG), 0);
        // Two alternates that both match, to show the last one winning.
        chooseCase(header, "snp-both", refA, List.of(altC, altC), 0);
        chooseCase(header, "mapq", refA, List.of(altC), 0);
        chooseCase(header, "insertion", refA, List.of(insertion), 0);
        // A deletion alternate whose bases are nothing like the pileup's, matched on length alone.
        chooseCase(header, "deletion", refAAA, List.of(deletion), 0);
        chooseCase(header, "deletion-other", refAAA, List.of(other), 0);
        chooseCase(header, "short", refAA, List.of(altGT), 0);
        chooseCase(header, "short-q1", refAA, List.of(altGT), 1);

        // The counter's map, which is the choice tallied and the unusable reads dropped.
        countCase(header, "snp", refA, List.of(altC), 0);
        countCase(header, "snp-q20", refA, List.of(altC), 20);
        countCase(header, "mapq", refA, List.of(altC), 0);
        countCase(header, "insertion", refA, List.of(insertion), 0);
        countCase(header, "deletion", refAAA, List.of(deletion), 0);
        countCase(header, "short", refAA, List.of(altGT), 0);
        countCase(header, "empty", refA, List.of(altC), 0);
        // A pileup the counter is never given, which is the null branch.
        System.out.printf("count\tnull\t%s%n", counted(new AllelePileupCounter(refA, List.of(altC), 0)));

        for (final String name : new String[] {"snp", "mapq", "insertion", "deletion", "short", "empty"}) {
            for (final int minQuality : new int[] {0, 20}) {
                System.out.printf("ratio\t%s,%d=%016x%n", name, minQuality,
                        Double.doubleToRawLongBits(PowerCalculationUtils.calculateMaxAltRatio(
                                pileup(header, name), name.equals("short") ? refAA : refA, minQuality)));
            }
        }

        support(header, "snp", refA, altC, 0);
        support(header, "snp", refA, altC, 20);
        support(header, "snp", refA, altG, 0);
        support(header, "mapq", refA, altC, 0);
        support(header, "insertion", refA, insertion, 0);
        support(header, "deletion", refAAA, deletion, 0);
        support(header, "short", refAA, altGT, 0);
        support(header, "empty", refA, altC, 0);

        // The genotypes, from the plainly validatable to each way of being refused.
        final List<BasicValidationResult> results = new ArrayList<>();
        validate(header, "plain", genotype(refA, altC, new int[] {40, 10}, null), refA, "snp", 8, 30, 0, results);
        validate(header, "no-alt-reads", genotype(refA, altC, new int[] {40, 10}, null), refA, "snp", 0, 30, 0, results);
        validate(header, "one-alt-read", genotype(refA, altC, new int[] {40, 10}, null), refA, "snp", 1, 30, 0, results);
        validate(header, "filtered", genotype(refA, altC, new int[] {40, 10}, "weak_evidence"), refA, "snp", 8, 30, 0, results);
        validate(header, "zero-discovery", genotype(refA, altC, new int[] {0, 0}, null), refA, "snp", 8, 30, 0, results);
        validate(header, "no-ad", genotype(refA, altC, null, null), refA, "snp", 8, 30, 0, results);
        validate(header, "alt-first", genotype(altC, refA, new int[] {40, 10}, null), refA, "snp", 8, 30, 0, results);
        validate(header, "complex", genotype(refAAA, Allele.create("TA", false), new int[] {40, 10}, null), refAAA, "snp", 8, 30, 0, results);
        validate(header, "insertion", genotype(refA, insertion, new int[] {40, 10}, null), refA, "insertion", 8, 30, 0, results);
        validate(header, "deletion", genotype(refAAA, deletion, new int[] {40, 10}, null), refAAA, "deletion", 8, 30, 0, results);
        validate(header, "empty-normal", genotype(refA, altC, new int[] {40, 10}, null), refA, "empty", 8, 30, 0, results);
        validate(header, "min-quality", genotype(refA, altC, new int[] {40, 10}, null), refA, "snp", 8, 30, 20, results);

        // The table those results are written as, which is the tool's own output format.
        final File file = File.createTempFile("basic-validation", ".tsv");
        file.deleteOnExit();
        BasicValidationResult.write(results, new GATKPath(file.getAbsolutePath()));
        for (final String line : Files.readAllLines(file.toPath())) {
            System.out.printf("table\t%s%n", line);
        }

        // A haploid genotype does not come back refused: `isAbleToValidateGenotype` reads the
        // second allele to type the variant before the ploidy it just computed is consulted, so the
        // list access throws first.
        error("haploid-genotype", () -> BasicSomaticShortMutationValidator
                .isAbleToValidateGenotype(haploid(refA), refA));
        error("symbolic-reference", () -> new AllelePileupCounter(
                Allele.create("<NON_REF>", true), List.of(altC), 0));
        error("non-reference-reference", () -> new AllelePileupCounter(altC, List.of(altG), 0));
        error("reference-alternate", () -> new AllelePileupCounter(refA, List.of(Allele.create("C", true)), 0));
        error("negative-quality", () -> new AllelePileupCounter(refA, List.of(altC), -1));
        error("negative-quality-ratio", () -> PowerCalculationUtils.calculateMaxAltRatio(
                pileup(header(), "snp"), refA, -1));
        error("negative-quality-support", () -> PowerCalculationUtils.calculateNumReadsSupportingAllele(
                pileup(header(), "snp"), refA, altC, -1));
        error("symbolic-type", () -> GATKVariantContextUtils.typeOfVariant(
                Allele.create("<NON_REF>", true), altC));
        error("span-del-reference", () -> GATKVariantContextUtils.typeOfVariant(Allele.SPAN_DEL, altC));
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(
                List.of(new SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample");
        group.setPlatform("ILLUMINA");
        header.addReadGroup(group);
        return header;
    }

    static ReadPileup pileup(final SAMFileHeader header, final String name) {
        final List<GATKRead> reads = new ArrayList<>();
        for (final String[] spec : CASES.get(name)) {
            reads.add(makeRead(header, spec[0], spec[1], spec[2], Byte.parseByte(spec[3]),
                    Integer.parseInt(spec[4])));
        }
        return new ReadPileup(new SimpleInterval("chr1", LOCUS, LOCUS), reads);
    }

    /** The case name without the suffix that says which alleles it was asked about. */
    static String pileupName(final String label) {
        final int dash = label.indexOf('-');
        return dash < 0 ? label : label.substring(0, dash);
    }

    static void chooseCase(final SAMFileHeader header, final String label, final Allele reference,
                           final List<Allele> alternates, final int minQuality) {
        for (final PileupElement element : pileup(header, pileupName(label))) {
            final Allele chosen = GATKVariantContextUtils.chooseAlleleForRead(
                    element, reference, alternates, minQuality);
            System.out.printf("choose\t%s\t%s=%s%n", label, element.getRead().getName(),
                    chosen == null ? "none" : chosen.getBaseString());
        }
    }

    static void countCase(final SAMFileHeader header, final String label, final Allele reference,
                          final List<Allele> alternates, final int minQuality) {
        final AllelePileupCounter counter = new AllelePileupCounter(reference, alternates,
                minQuality, pileup(header, pileupName(label)));
        System.out.printf("count\t%s\t%s%n", label, counted(counter));
    }

    /** The count map, ordered by allele so the row does not depend on the hash order. */
    static String counted(final AllelePileupCounter counter) {
        final Map<Allele, MutableInt> counts = counter.getCountMap();
        final List<Allele> alleles = new ArrayList<>(counts.keySet());
        alleles.sort((left, right) -> left.getBaseString().compareTo(right.getBaseString()));
        final StringBuilder text = new StringBuilder();
        for (final Allele allele : alleles) {
            if (text.length() > 0) {
                text.append(',');
            }
            text.append(allele.getBaseString()).append('=').append(counts.get(allele).intValue());
        }
        return text.toString();
    }

    static void support(final SAMFileHeader header, final String name, final Allele reference,
                        final Allele alternate, final int minQuality) {
        System.out.printf("support\t%s,%s,%d=%d%n", name, alternate.getBaseString(), minQuality,
                PowerCalculationUtils.calculateNumReadsSupportingAllele(
                        pileup(header, name), reference, alternate, minQuality));
    }

    static Genotype genotype(final Allele first, final Allele second, final int[] depths,
                             final String filters) {
        final GenotypeBuilder builder = new GenotypeBuilder("sample", Arrays.asList(first, second));
        if (depths != null) {
            builder.AD(depths);
        }
        if (filters != null) {
            builder.filters(filters);
        }
        return builder.make();
    }

    static Genotype haploid(final Allele allele) {
        return new GenotypeBuilder("sample", List.of(allele)).AD(new int[] {40, 10}).make();
    }

    static void validate(final SAMFileHeader header, final String label, final Genotype genotype,
                         final Allele reference, final String pileupCase, final int altCount,
                         final int totalCount, final int minQuality,
                         final List<BasicValidationResult> results) {
        System.out.printf("validatable\t%s=%s%n", label,
                BasicSomaticShortMutationValidator.isAbleToValidateGenotype(genotype, reference));
        final BasicValidationResult result = BasicSomaticShortMutationValidator
                .calculateBasicValidationResult(genotype, reference, pileup(header, pileupCase),
                        altCount, totalCount, minQuality, new SimpleInterval("chr1", LOCUS, LOCUS),
                        "PASS");
        if (result == null) {
            System.out.printf("result\t%s=none%n", label);
            return;
        }
        results.add(result);
        System.out.printf("result\t%s=%s,%d,%d,%s,%s,%d,%b,%b,%016x,%d,%d,%d,%d,%s,%d%n", label,
                result.getContig(), result.getStart(), result.getEnd(),
                result.getReference().getBaseString(), result.getAlternate().getBaseString(),
                result.getMinValidationReadCount(), result.isEnoughValidationReads(),
                result.isOutOfNoiseFloor(), Double.doubleToRawLongBits(result.getPower()),
                result.getValidationAltCount(), result.getValidationRefCount(),
                result.getDiscoveryAltCount(), result.getDiscoveryRefCount(),
                result.getFilters(), result.getNumAltSupportingReadsInNormal());
    }

    static GATKRead makeRead(final SAMFileHeader header, final String name, final String cigar,
                             final String bases, final byte quality, final int mappingQuality) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceName("chr1");
        record.setAlignmentStart(START);
        record.setCigar(TextCigarCodec.decode(cigar));
        record.setReadBases(bases.getBytes());
        final byte[] quals = new byte[bases.length()];
        Arrays.fill(quals, quality);
        record.setBaseQualities(quals);
        record.setMappingQuality(mappingQuality);
        record.setAttribute("RG", "rg1");
        return new SAMRecordToGATKReadAdapter(record);
    }

    interface Body {
        Object run();
    }

    static void error(final String label, final Body body) {
        try {
            System.out.printf("unexpected\t%s\t%s%n", label, body.run());
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }
}
