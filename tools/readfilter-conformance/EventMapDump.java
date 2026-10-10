/*
 * EventMap: the events a haplotype's CIGAR makes against the reference.
 *
 * `EventMap.fromHaplotype` walks the CIGAR from `alignmentStartHapwrtRef`: an insertion becomes an
 * event padded with the base before it, unless it is the first or last element or starts the
 * contig; a deletion likewise unless it starts the contig; mismatches in an M, = or X element are
 * grouped into MNPs while the gap to the next is at most `maxMnpDistance`; any base that is not
 * A, C, G or T (either case) drops the event; S is skipped and N, H and P throw. Events sharing a
 * start are compounded (a SNP with an insertion or a deletion, an insertion with a deletion, and
 * two SNPs refused), and every `Event` is cut to its minimal representation by dropping the
 * trailing bases its alleles share.
 *
 * Output:
 *
 *     case\t<label>\t<ref>\t<ref start>\t<haplotype>\t<cigar>\t<alignment start>\t<max mnp distance>
 *     event\t<label>\t<index>\t<start>\t<end>\t<ref bases>\t<alt bases>\t<snp>\t<indel>\t<insertion>\t<deletion>\t<mnp>
 *     overlap\t<label>\t<locus>\t<index of each overlapping event, in order,...>
 *     error\t<label>\t<exception class>: <message>
 *     starts\t<label>\t<positions,...>          (getEventStartPositions over every case built)
 *     direct\t<label>\t<start>\t<end>\t<ref bases>\t<alt bases>   (new Event, or makeCompoundEvents)
 *
 * Usage: EventMapDump
 */

import htsjdk.samtools.TextCigarCodec;
import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.haplotype.Event;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class EventMapDump {

    static final String REF = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGT";
    static final int REF_START = 1000;
    static final List<Haplotype> BUILT = new ArrayList<>();

    public static void main(final String[] args) {
        System.out.println("# EventMapDump: EventMap.fromHaplotype");
        // Substitutions and their grouping into MNPs.
        run("ref", REF, REF.length() + "M", 0, 0);
        run("snp", mutate(REF, 10), REF.length() + "M", 0, 0);
        run("snp-first-base", mutate(REF, 0), REF.length() + "M", 0, 0);
        run("adjacent-snps-mnp0", mutate(mutate(REF, 10), 11), REF.length() + "M", 0, 0);
        run("adjacent-snps-mnp1", mutate(mutate(REF, 10), 11), REF.length() + "M", 0, 1);
        run("gap1-mnp0", mutate(mutate(REF, 10), 12), REF.length() + "M", 0, 0);
        run("gap1-mnp1", mutate(mutate(REF, 10), 12), REF.length() + "M", 0, 1);
        run("gap1-mnp2", mutate(mutate(REF, 10), 12), REF.length() + "M", 0, 2);
        run("chain-mnp1", mutate(mutate(mutate(mutate(mutate(mutate(REF, 10), 11), 12), 14), 15), 17), REF.length() + "M", 0, 1);
        run("eq-and-x", mutate(REF, 10), "10=1X" + (REF.length() - 11) + "=", 0, 0);
        run("lowercase", REF.substring(0, 10) + "t" + REF.substring(11), REF.length() + "M", 0, 0);
        run("star-in-haplotype", REF.substring(0, 10) + "*" + REF.substring(11), REF.length() + "M", 0, 0);
        run("star-in-mnp", REF.substring(0, 10) + "*C" + REF.substring(12), REF.length() + "M", 0, 1);
        run("n-in-haplotype", REF.substring(0, 10) + "N" + REF.substring(11), REF.length() + "M", 0, 0);

        // Indels.
        final String ins = REF.substring(0, 20) + "GG" + REF.substring(20);
        run("insertion", ins, "20M2I" + (REF.length() - 20) + "M", 0, 0);
        run("insertion-first-element", "GG" + REF, "2I" + REF.length() + "M", 0, 0);
        run("insertion-last-element", REF + "GG", REF.length() + "M2I", 0, 0);
        run("insertion-at-contig-start", "TTGG" + REF, "2S2I" + REF.length() + "M", 0, 0);
        run("insertion-with-n", REF.substring(0, 20) + "GN" + REF.substring(20), "20M2I" + (REF.length() - 20) + "M", 0, 0);
        final String del = REF.substring(0, 20) + REF.substring(23);
        run("deletion", del, "20M3D" + (REF.length() - 23) + "M", 0, 0);
        run("deletion-at-contig-start", REF.substring(3), "3D" + (REF.length() - 3) + "M", 0, 0);
        run("leading-deletion-in-window", REF.substring(8), "3D" + (REF.length() - 8) + "M", 5, 0);

        // Events sharing a start.
        final String snpIns = mutate(REF, 19).substring(0, 20) + "GG" + REF.substring(20);
        run("snp-plus-insertion", snpIns, "20M2I" + (REF.length() - 20) + "M", 0, 0);
        final String snpDel = mutate(REF, 19).substring(0, 20) + REF.substring(23);
        run("snp-plus-deletion", snpDel, "20M3D" + (REF.length() - 23) + "M", 0, 0);
        final String insDel = REF.substring(0, 20) + "TA" + REF.substring(22);
        run("insertion-plus-deletion", insDel, "20M2I2D" + (REF.length() - 22) + "M", 0, 0);
        // The deletion's last base is the insertion's: the compound is trimmed.
        final String trimmed = REF.substring(0, 20) + "C" + REF.charAt(21) + REF.substring(22);
        run("insertion-plus-deletion-trimmed", trimmed, "20M2I2D" + (REF.length() - 22) + "M", 0, 0);
        final String three = mutate(REF, 19).substring(0, 20) + "TA" + REF.substring(22);
        run("snp-insertion-deletion", three, "20M2I2D" + (REF.length() - 22) + "M", 0, 0);
        run("deletion-then-insertion", insDel, "20M2D2I" + (REF.length() - 22) + "M", 0, 0);

        // Offsets into the reference, clips, and refusals.
        run("window-offset", REF.substring(5, 40), "35M", 5, 0);
        run("window-offset-snp", mutate(REF, 12).substring(5, 40), "35M", 5, 0);
        run("negative-start", REF, REF.length() + "M", -1, 0);
        run("soft-clip", "TTTT" + mutate(REF, 10).substring(0, 30), "4S30M", 0, 0);
        run("hard-clip", mutate(REF, 10), "2H" + REF.length() + "M", 0, 0);
        run("negative-mnp", mutate(REF, 10), REF.length() + "M", 0, -1);

        // Events built directly: the minimal representation, and compounds no canonical CIGAR makes.
        direct("trim-two", "ACGT", "AGT");
        direct("trim-to-mnp", "ACGTT", "AGCTT");
        direct("no-trim-length-one", "A", "AGT");
        direct("no-trim-last-base-differs", "ACG", "AT");
        direct("identical", "ACG", "ACG");
        direct("lowercase-shared", "ACgt", "AGgt");
        compound("snp+snp", event(1019, "G", "T"), event(1019, "G", "C"));
        compound("snp+insertion", event(1019, "G", "T"), event(1019, "G", "GCC"));
        compound("insertion+snp", event(1019, "G", "GCC"), event(1019, "G", "T"));
        compound("snp+deletion", event(1019, "G", "T"), event(1019, "GAT", "G"));
        compound("deletion+insertion", event(1019, "GAT", "G"), event(1019, "G", "GCT"));
        compound("insertion+deletion", event(1019, "G", "GCC"), event(1019, "GAT", "G"));
        compound("insertion+insertion", event(1019, "G", "GCC"), event(1019, "G", "GA"));
        compound("different-starts", event(1019, "G", "T"), event(1020, "A", "AC"));
        try {
            new Event("chr1", 1019, Allele.create("G", false), Allele.create("T", false));
        } catch (final Exception e) {
            System.out.println("error\tnon-reference-ref\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }

        // getEventStartPositions over every haplotype built.
        final java.util.TreeSet<Integer> starts = EventMap.getEventStartPositions(BUILT);
        System.out.println("starts\tall\t" + starts.stream().map(String::valueOf).collect(Collectors.joining(",")));
    }

    static Event event(final int start, final String ref, final String alt) {
        return new Event("chr1", start, Allele.create(ref, true), Allele.create(alt, false));
    }

    static void print(final String label, final Event e) {
        System.out.println("direct\t" + label + "\t" + e.getStart() + "\t" + e.getEnd() + "\t"
                + e.refAllele().getBaseString() + "\t" + e.altAllele().getBaseString());
    }

    static void direct(final String label, final String ref, final String alt) {
        try {
            print(label, event(1019, ref, alt));
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }

    static void compound(final String label, final Event e1, final Event e2) {
        try {
            final java.lang.reflect.Method m = EventMap.class.getDeclaredMethod("makeCompoundEvents", Event.class, Event.class);
            m.setAccessible(true);
            print(label, (Event) m.invoke(null, e1, e2));
        } catch (final java.lang.reflect.InvocationTargetException e) {
            System.out.println("error\t" + label + "\t" + e.getCause().getClass().getSimpleName() + ": " + e.getCause().getMessage());
        } catch (final ReflectiveOperationException e) {
            throw new RuntimeException(e);
        }
    }

    static String mutate(final String s, final int at) {
        final char c = s.charAt(at);
        final char alt = c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A';
        return s.substring(0, at) + alt + s.substring(at + 1);
    }

    static void run(final String label, final String bases, final String cigar, final int alignmentStart, final int maxMnp) {
        System.out.println("case\t" + label + "\t" + REF + "\t" + REF_START + "\t" + bases + "\t" + cigar + "\t" + alignmentStart + "\t" + maxMnp);
        try {
            final Haplotype h = new Haplotype(bases.getBytes(StandardCharsets.US_ASCII), false);
            h.setCigar(TextCigarCodec.decode(cigar));
            h.setAlignmentStartHapwrtRef(alignmentStart);
            final SimpleInterval refLoc = new SimpleInterval("chr1", REF_START, REF_START + REF.length() - 1);
            final EventMap map = EventMap.fromHaplotype(h, REF.getBytes(StandardCharsets.US_ASCII), refLoc, maxMnp);
            h.setEventMap(map);
            BUILT.add(h);
            final List<Event> events = new ArrayList<>(map.getEvents());
            for (int i = 0; i < events.size(); i++) {
                final Event e = events.get(i);
                System.out.println("event\t" + label + "\t" + i + "\t" + e.getStart() + "\t" + e.getEnd()
                        + "\t" + e.refAllele().getBaseString() + "\t" + e.altAllele().getBaseString()
                        + "\t" + e.isSNP() + "\t" + e.isIndel() + "\t" + e.isSimpleInsertion() + "\t" + e.isSimpleDeletion() + "\t" + e.isMNP());
            }
            for (int locus = REF_START - 1; locus <= REF_START + REF.length(); locus++) {
                final List<Event> overlapping = map.getOverlappingEvents(locus);
                if (!overlapping.isEmpty()) {
                    System.out.println("overlap\t" + label + "\t" + locus + "\t"
                            + overlapping.stream().map(e -> String.valueOf(events.indexOf(e))).collect(Collectors.joining(",")));
                }
            }
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }
}
