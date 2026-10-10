/*
 * AssemblyRegionTrimmer: the part of an assembled region that is genotyped, and its flanks.
 *
 * The default trim spans the events overlapping the active span and pads each one, a SNP by 20,
 * an indel by 75, and an indel in a tandem repeat by 75 plus the longest run of the repeat either
 * allele carries in the reference after it (`TandemRepeat.getNumTandemRepeatUnits` on a reference
 * context over the padded span, as every caller builds it). The padded span is cut to the region's
 * padded span. The legacy trim pads the variant span once and caps it at the active span plus 25.
 * The result then trims the region (`getVariantRegion`) and gives its two non-variant flanks,
 * padded by `assemblyRegionPadding`.
 *
 * Output:
 *
 *     contig\t<bases>                                  (chr1, once)
 *     read\t<name>\t<start>\t<cigar>                    (the reads every region is given, once)
 *     case\t<label>\t<active start>\t<active end>\t<padding>\t<legacy>\t<snp>\t<indel>\t<str>\t<max extension>\t<region padding>
 *     event\t<label>\t<start>\t<ref>\t<alt>
 *     result\t<label>\t<variant span or null>\t<padded span or null>
 *     region\t<label>\t<which>\t<span>\t<padded span>\t<read names,...>   (variant, left, right; or "empty")
 *     error\t<label>\t<which>\t<exception class>: <message>
 *
 * Usage: AssemblyRegionTrimmerDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.ReferenceMemorySource;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.Random;
import java.util.TreeSet;
import java.util.function.Supplier;
import java.util.stream.Collectors;

public class AssemblyRegionTrimmerDump {

    static String contig;
    static SAMFileHeader header;
    static SAMSequenceDictionary dict;
    static final List<GATKRead> READS = new ArrayList<>();

    public static void main(final String[] args) {
        System.out.println("# AssemblyRegionTrimmerDump: AssemblyRegionTrimmer.trim");
        final Random random = new Random(20261011L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 1000; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
            // A dinucleotide repeat at 401-430 and a trinucleotide one at 601-636.
            if (i == 399) {
                sb.append("AC".repeat(15));
                i += 30;
            }
            if (i == 599) {
                sb.append("GTT".repeat(12));
                i += 36;
            }
        }
        contig = sb.substring(0, 1000);
        System.out.println("contig\t" + contig);
        dict = new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", contig.length())));
        header = new SAMFileHeader();
        header.setSequenceDictionary(dict);
        for (int start = 1; start <= 950; start += 37) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "r" + start, 0, start, 50);
            READS.add(read);
            System.out.println("read\tr" + start + "\t" + start + "\t50M");
        }

        final int[] hc = {20, 75, 75, 25, 100};
        run("no-events", 300, 500, 100, false, hc);
        run("one-snp", 300, 500, 100, false, hc, ev(350, 1, "T"));
        run("snp-at-active-start", 300, 500, 100, false, hc, ev(300, 1, "T"));
        run("snp-in-padding-only", 300, 500, 100, false, hc, ev(250, 1, "T"));
        run("snps-apart", 300, 500, 100, false, hc, ev(310, 1, "T"), ev(480, 1, "G"));
        run("deletion", 300, 500, 100, false, hc, ev(350, 4, null));
        run("insertion", 300, 500, 100, false, hc, ev(350, 1, "+GGA"));
        run("str-insertion", 300, 500, 100, false, hc, ev(400, 1, "+AC"));
        run("str-deletion", 300, 500, 100, false, hc, ev(400, 5, null));
        run("str-tri-deletion", 500, 700, 100, false, hc, ev(600, 4, null));
        run("str-tri-insertion", 500, 700, 100, false, hc, ev(600, 1, "+GTTGTT"));
        run("not-a-repeat-insertion", 300, 500, 100, false, hc, ev(400, 1, "+GG"));
        run("mnp", 300, 500, 100, false, hc, ev(350, 3, "x"));
        run("snp-and-str", 300, 500, 100, false, hc, ev(320, 1, "T"), ev(400, 1, "+AC"));
        run("event-spanning-active-end", 300, 500, 100, false, hc, ev(495, 10, null));
        run("near-contig-start", 10, 60, 20, false, hc, ev(15, 1, "T"), ev(30, 3, null));
        run("near-contig-end", 950, 990, 20, false, hc, ev(985, 1, "T"));
        run("small-padding", 300, 500, 10, false, hc, ev(350, 4, null));
        run("zero-snp-padding", 300, 500, 100, false, new int[]{0, 75, 75, 25, 100}, ev(350, 1, "T"));
        run("region-padding-0", 300, 500, 100, false, new int[]{20, 75, 75, 25, 0}, ev(350, 1, "T"));
        run("legacy-no-events", 300, 500, 100, true, hc);
        run("legacy-snp", 300, 500, 100, true, hc, ev(350, 1, "T"));
        run("legacy-indel", 300, 500, 100, true, hc, ev(350, 4, null));
        run("legacy-snp-at-edge", 300, 500, 100, true, hc, ev(498, 1, "T"));
        run("legacy-outside-only", 300, 500, 100, true, hc, ev(250, 1, "T"));
        run("legacy-snps-apart", 300, 500, 100, true, hc, ev(301, 1, "T"), ev(499, 1, "G"));
        run("legacy-str", 300, 500, 100, true, hc, ev(400, 1, "+AC"));
        run("legacy-near-contig-start", 10, 60, 20, true, hc, ev(12, 1, "T"));
        run("negative-snp-padding", 300, 500, 100, false, new int[]{-1, 75, 75, 25, 100}, ev(350, 1, "T"));
    }

    /**
     * An event at `start` whose reference allele is `refLength` contig bases: `null` alt is the
     * deletion of all but the first, `+bases` an insertion after the first, `x` each base shifted
     * round ACGT, anything else that literal alt.
     */
    static Event ev(final int start, final int refLength, final String alt) {
        final String ref = contig.substring(start - 1, start - 1 + refLength);
        final String a;
        if (alt == null) {
            a = ref.substring(0, 1);
        } else if (alt.startsWith("+")) {
            a = ref + alt.substring(1);
        } else if (alt.equals("x")) {
            final StringBuilder b = new StringBuilder();
            for (final char c : ref.toCharArray()) {
                b.append("ACGT".charAt(("ACGT".indexOf(c) + 1) % 4));
            }
            a = b.toString();
        } else {
            a = alt;
        }
        return new Event("chr1", start, Allele.create(ref, true), Allele.create(a, false));
    }

    static String render(final SimpleInterval i) {
        return i == null ? "null" : i.getContig() + ":" + i.getStart() + "-" + i.getEnd();
    }

    static void region(final String label, final String which, final Supplier<Optional<AssemblyRegion>> f) {
        try {
            final Optional<AssemblyRegion> r = f.get();
            if (r.isEmpty()) {
                System.out.println("region\t" + label + "\t" + which + "\tempty");
            } else {
                System.out.println("region\t" + label + "\t" + which + "\t" + render(r.get().getSpan()) + "\t" + render(r.get().getPaddedSpan())
                        + "\t" + r.get().getReads().stream().map(GATKRead::getName).collect(Collectors.joining(",")));
            }
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + which + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }

    static void run(final String label, final int activeStart, final int activeEnd, final int padding, final boolean legacy,
                    final int[] p, final Event... events) {
        System.out.println("case\t" + label + "\t" + activeStart + "\t" + activeEnd + "\t" + padding + "\t" + legacy
                + "\t" + p[0] + "\t" + p[1] + "\t" + p[2] + "\t" + p[3] + "\t" + p[4]);
        for (final Event e : events) {
            System.out.println("event\t" + label + "\t" + e.getStart() + "\t" + e.refAllele().getBaseString() + "\t" + e.altAllele().getBaseString());
        }
        final AssemblyRegion region = new AssemblyRegion(new SimpleInterval("chr1", activeStart, activeEnd), padding, header);
        for (final GATKRead read : READS) {
            if (region.getPaddedSpan().overlaps(read)) {
                region.add(read);
            }
        }
        final AssemblyRegionArgumentCollection a = new AssemblyRegionArgumentCollection();
        a.snpPaddingForGenotyping = p[0];
        a.indelPaddingForGenotyping = p[1];
        a.strPaddingForGenotyping = p[2];
        a.maxExtensionIntoRegionPadding = p[3];
        a.assemblyRegionPadding = p[4];
        a.enableLegacyAssemblyRegionTrimming = legacy;
        final AssemblyRegionTrimmer.Result result;
        try {
            final AssemblyRegionTrimmer trimmer = new AssemblyRegionTrimmer(a, dict);
            final ReferenceContext ref = new ReferenceContext(new ReferenceMemorySource(
                    new ReferenceBases(contig.getBytes(StandardCharsets.US_ASCII), new SimpleInterval("chr1", 1, contig.length())), dict),
                    region.getPaddedSpan());
            final TreeSet<Event> sorted = new TreeSet<>(org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet.HAPLOTYPE_EVENT_COMPARATOR);
            sorted.addAll(List.of(events));
            result = trimmer.trim(region, sorted, ref);
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\ttrim\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
            return;
        }
        SimpleInterval variant = null;
        SimpleInterval padded = null;
        try {
            final java.lang.reflect.Field vf = AssemblyRegionTrimmer.Result.class.getDeclaredField("variantSpan");
            final java.lang.reflect.Field pf = AssemblyRegionTrimmer.Result.class.getDeclaredField("paddedSpan");
            vf.setAccessible(true);
            pf.setAccessible(true);
            variant = (SimpleInterval) vf.get(result);
            padded = (SimpleInterval) pf.get(result);
        } catch (final ReflectiveOperationException e) {
            throw new RuntimeException(e);
        }
        System.out.println("result\t" + label + "\t" + render(variant) + "\t" + render(padded));
        region(label, "variant", () -> Optional.of(result.getVariantRegion()));
        region(label, "left", result::nonVariantLeftFlankRegion);
        region(label, "right", result::nonVariantRightFlankRegion);
    }
}
