/*
 * ReadThreadingAssembler.runLocalAssembly over fixed regions, in the default sequence-graph mode.
 *
 * Every graph stage the assembler runs is measured on its own by an earlier dump; this one measures
 * the orchestration: the k-mer sizes tried (sorted, then up to six retries ten larger when none
 * gives a graph), the refusals of createGraph (a reference shorter than k, non-unique reference
 * k-mers, cycles after pruning, a low-complexity graph), the order of the stages, the k best
 * haplotypes of every graph, their CIGARs against the reference haplotype
 * (`CigarUtils.calculateCigar` with NEW_SW_PARAMETERS, SOFTCLIP then INDEL), the haplotypes
 * dropped (an N, or under 30 reference bases), and the order of the result set, a LinkedHashSet in
 * which the reference haplotype is removed as an alternative and re-added at the end.
 *
 * Each case lays its inputs out first, so that the Rust test rebuilds the same region from the
 * golden: the contig, the active span and padding, the reference padding, the assembler's
 * settings, then every read with its sample, start, cigar, bases and qualities (Phred + 33).
 * Reads are soft-clip-hard-clipped by the assembler itself.
 *
 * Output:
 *
 *     case\t<label>\t<contig bases>\t<active start>\t<active end>\t<padding>\t<reference padding>
 *     settings\t<label>\t<key>=<value>\t...
 *     read\t<label>\t<name>\t<sample>\t<start>\t<cigar>\t<bases>\t<qualities>
 *     result\t<label>\t<index>\t<kmer or ->\t<status>\t<vertices>\t<edges>
 *     haplotype\t<label>\t<index>\t<ref>\t<bases>\t<cigar>\t<alignment start>\t<score bits>,<score>\t<kmer>\t<location>
 *     set\t<label>\t<count>\t<variation present>\t<reference index>
 *     error\t<label>\t<exception class>: <message>
 *
 * `result` rows are the package-private `assemble` on the same reads (what each k-mer size gave);
 * `haplotype` and `set` rows are `runLocalAssembly`'s result set.
 *
 * Usage: ReadThreadingAssemblerDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResult;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class ReadThreadingAssemblerDump {

    static final String CONTIG = "chr1";

    /** One read: its sample, start, cigar, bases and qualities. */
    record Read(String sample, int start, String cigar, String bases, byte[] quals) {
    }

    /** The assembler's settings, HaplotypeCaller's defaults unless a case changes them. */
    static final class Settings {
        List<Integer> kmers = List.of(10, 25);
        boolean dontIncrease = false;
        boolean allowNonUnique = false;
        int pruningSamples = 1;
        int pruneFactor = 2;
        boolean adaptive = false;
        double initialErrorRate = 0.001;
        double pruningLogOdds = MathUtils.log10ToLog(1.0);
        double seedingLogOdds = MathUtils.log10ToLog(4.0);
        int maxUnprunedVariants = 100;
        int maxPaths = 128;
        boolean legacyCycles = false;
        int minMatching = -1;
        boolean recoverDangling = true;
        boolean recoverAll = false;
        int minDanglingLength = 4;
        int minBaseQuality = 10;

        Settings copy() {
            final Settings s = new Settings();
            s.kmers = kmers;
            s.dontIncrease = dontIncrease;
            s.allowNonUnique = allowNonUnique;
            s.pruningSamples = pruningSamples;
            s.pruneFactor = pruneFactor;
            s.adaptive = adaptive;
            s.initialErrorRate = initialErrorRate;
            s.pruningLogOdds = pruningLogOdds;
            s.seedingLogOdds = seedingLogOdds;
            s.maxUnprunedVariants = maxUnprunedVariants;
            s.maxPaths = maxPaths;
            s.legacyCycles = legacyCycles;
            s.minMatching = minMatching;
            s.recoverDangling = recoverDangling;
            s.recoverAll = recoverAll;
            s.minDanglingLength = minDanglingLength;
            s.minBaseQuality = minBaseQuality;
            return s;
        }

        String render() {
            return "kmers=" + kmers.stream().map(String::valueOf).collect(Collectors.joining(","))
                    + "\tdontIncrease=" + dontIncrease
                    + "\tallowNonUnique=" + allowNonUnique
                    + "\tpruningSamples=" + pruningSamples
                    + "\tpruneFactor=" + pruneFactor
                    + "\tadaptive=" + adaptive
                    + "\tinitialErrorRate=" + Long.toHexString(Double.doubleToRawLongBits(initialErrorRate))
                    + "\tpruningLogOdds=" + Long.toHexString(Double.doubleToRawLongBits(pruningLogOdds))
                    + "\tseedingLogOdds=" + Long.toHexString(Double.doubleToRawLongBits(seedingLogOdds))
                    + "\tmaxUnprunedVariants=" + maxUnprunedVariants
                    + "\tmaxPaths=" + maxPaths
                    + "\tlegacyCycles=" + legacyCycles
                    + "\tminMatching=" + minMatching
                    + "\trecoverDangling=" + recoverDangling
                    + "\trecoverAll=" + recoverAll
                    + "\tminDanglingLength=" + minDanglingLength
                    + "\tminBaseQuality=" + minBaseQuality;
        }
    }

    static String contig;

    public static void main(final String[] args) throws Exception {
        System.out.println("# ReadThreadingAssemblerDump: ReadThreadingAssembler.runLocalAssembly");
        final Random random = new Random(20261010L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 400; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        contig = sb.toString();
        final Settings hc = new Settings();

        // The active span 150-250 padded by 50: the reference haplotype is 101-300 (1-based), and
        // the reference with 5 more bases each side starts at 96.
        run("ref-only", contig, 150, 250, 50, 5, hc, List.of());
        run("ref-reads", contig, 150, 250, 50, 5, hc, tile(contig, "s1", 101, 300, 4, 30));
        run("snp", contig, 150, 250, 50, 5, hc, het(contig, snp(contig, 200, 'x')));
        run("snp-minor", contig, 150, 250, 50, 5, hc,
                concat(tile(contig, "s1", 101, 300, 3, 30), tile(snp(contig, 200, 'x'), "s1", 101, 300, 11, 30)));
        run("insertion", contig, 150, 250, 50, 5, hc, het(contig, insert(contig, 190, "GATTACA")));
        run("deletion", contig, 150, 250, 50, 5, hc, het(contig, delete(contig, 210, 6)));
        run("two-snps", contig, 150, 250, 50, 5, hc, het(contig, snp(snp(contig, 170, 'x'), 230, 'x')));
        run("two-alts", contig, 150, 250, 50, 5, hc,
                concat(tile(snp(contig, 180, 'x'), "s1", 101, 300, 4, 30), tile(snp(contig, 220, 'x'), "s1", 101, 300, 4, 30)));
        run("mnp", contig, 150, 250, 50, 5, hc, het(contig, snp(snp(snp(contig, 200, 'x'), 201, 'x'), 202, 'x')));
        run("snp-near-edge", contig, 150, 250, 50, 5, hc, het(contig, snp(contig, 106, 'x')));
        run("errors-pruned", contig, 150, 250, 50, 5, hc, withErrors(het(contig, snp(contig, 200, 'x'))));

        Settings s = hc.copy();
        s.pruneFactor = 0;
        run("errors-unpruned", contig, 150, 250, 50, 5, s, withErrors(het(contig, snp(contig, 200, 'x'))));

        s = hc.copy();
        s.adaptive = true;
        s.pruneFactor = 0;
        run("adaptive", contig, 150, 250, 50, 5, s, withErrors(het(contig, snp(contig, 200, 'x'))));

        s = hc.copy();
        s.maxPaths = 2;
        run("max-paths-2", contig, 150, 250, 50, 5, s,
                concat(tile(snp(contig, 160, 'x'), "s1", 101, 300, 4, 30),
                        tile(snp(contig, 200, 'x'), "s1", 101, 300, 4, 30),
                        tile(snp(contig, 240, 'x'), "s1", 101, 300, 4, 30)));

        s = hc.copy();
        s.kmers = List.of(25, 10, 15);
        run("kmers-unsorted", contig, 150, 250, 50, 5, s, het(contig, snp(contig, 200, 'x')));

        s = hc.copy();
        s.kmers = List.of(25);
        run("kmer-25", contig, 150, 250, 50, 5, s, het(contig, insert(contig, 190, "GATTACA")));

        // A 30-base segment copied in place: the reference's k-mers of 10 and 25 repeat, 35 is the
        // first retry and does not.
        final String repeated = contig.substring(0, 180) + contig.substring(150, 180) + contig.substring(180);
        run("repeat-retry", repeated, 150, 250, 50, 5, hc, het(repeated, snp(repeated, 230, 'x')));
        s = hc.copy();
        s.dontIncrease = true;
        run("repeat-dont-increase", repeated, 150, 250, 50, 5, s, het(repeated, snp(repeated, 230, 'x')));
        s = hc.copy();
        s.allowNonUnique = true;
        run("repeat-allow-non-unique", repeated, 150, 250, 50, 5, s, het(repeated, snp(repeated, 230, 'x')));

        // A short tandem repeat whose expansion the reads carry: candidate cycles.
        final String str = contig.substring(0, 190) + "ACACACACACACACACACAC" + contig.substring(190);
        run("str-expansion", str, 150, 250, 50, 5, hc, het(str, insert(str, 195, "ACACACAC")));
        s = hc.copy();
        s.allowNonUnique = true;
        run("str-allow-non-unique", str, 150, 250, 50, 5, s, het(str, insert(str, 195, "ACACACAC")));
        s = hc.copy();
        s.legacyCycles = true;
        run("str-legacy-cycles", str, 150, 250, 50, 5, s, het(str, insert(str, 195, "ACACACAC")));

        // A homopolymer-heavy region: a low-complexity graph at 10.
        final String poly = contig.substring(0, 140) + "AC".repeat(50) + contig.substring(140);
        s = hc.copy();
        s.allowNonUnique = true;
        run("low-complexity", poly, 150, 250, 50, 5, s, het(poly, snp(poly, 260, 'x')));
        s.dontIncrease = true;
        run("low-complexity-allowed", poly, 150, 250, 50, 5, s, het(poly, snp(poly, 260, 'x')));

        // A region too short for any haplotype to keep 30 reference bases.
        run("short-region", contig, 200, 210, 5, 2, hc, het(contig, snp(contig, 205, 'x')));
        s = hc.copy();
        s.kmers = List.of(25);
        run("kmer-longer-than-reference", contig, 200, 210, 5, 2, s, het(contig, snp(contig, 205, 'x')));

        // Soft clips the assembler hard-clips away, low qualities and Ns it skips.
        run("soft-clips", contig, 150, 250, 50, 5, hc, softClipped(het(contig, snp(contig, 200, 'x'))));
        run("low-quality", contig, 150, 250, 50, 5, hc, lowQuality(het(contig, snp(contig, 200, 'x'))));

        // Dangling ends: reads that stop at a mismatch near their end.
        run("dangling", contig, 150, 250, 50, 5, hc, concat(tile(contig, "s1", 101, 300, 4, 30), dangling(contig, 190)));
        s = hc.copy();
        s.recoverDangling = false;
        run("dangling-off", contig, 150, 250, 50, 5, s, concat(tile(contig, "s1", 101, 300, 4, 30), dangling(contig, 190)));
        s = hc.copy();
        s.recoverAll = true;
        s.minDanglingLength = 1;
        run("dangling-all", contig, 150, 250, 50, 5, s, concat(tile(contig, "s1", 101, 300, 4, 30), dangling(contig, 190)));
        s = hc.copy();
        s.minMatching = 2;
        run("dangling-min-matching", contig, 150, 250, 50, 5, s, concat(tile(contig, "s1", 101, 300, 4, 30), dangling(contig, 190)));

        run("dangling-tail", contig, 150, 250, 50, 5, hc, concat(tile(contig, "s1", 101, 300, 4, 30), danglingTail(contig, 190, 44)));
        s = hc.copy();
        s.recoverDangling = false;
        run("dangling-tail-off", contig, 150, 250, 50, 5, s, concat(tile(contig, "s1", 101, 300, 4, 30), danglingTail(contig, 190, 44)));
        run("dangling-head", contig, 150, 250, 50, 5, hc, concat(tile(contig, "s1", 101, 300, 4, 30), danglingHead(contig, 210, 6)));
        s = hc.copy();
        s.minMatching = 0;
        run("dangling-head-min-matching-0", contig, 150, 250, 50, 5, s, concat(tile(contig, "s1", 101, 300, 4, 30), danglingHead(contig, 210, 6)));

        // Two samples: one carrier, then both.
        s = hc.copy();
        s.pruningSamples = 2;
        run("two-samples-one-carrier", contig, 150, 250, 50, 5, s,
                concat(tile(contig, "s1", 101, 300, 4, 30), tile(contig, "s2", 102, 300, 4, 30),
                        tile(snp(contig, 200, 'x'), "s1", 103, 300, 4, 30)));
        run("two-samples-both-carry", contig, 150, 250, 50, 5, s,
                concat(tile(contig, "s1", 101, 300, 4, 30), tile(contig, "s2", 102, 300, 4, 30),
                        tile(snp(contig, 200, 'x'), "s1", 103, 300, 8, 30), tile(snp(contig, 200, 'x'), "s2", 105, 300, 8, 30)));

        // A 30-base segment from later in the reference copied earlier in the reads: a cycle at 10
        // and at 25, none at 35.
        final String moved = contig.substring(0, 170) + contig.substring(220, 250) + contig.substring(170);
        run("cycle", contig, 150, 250, 50, 5, hc, het(contig, moved, 101, 330));
        s = hc.copy();
        s.dontIncrease = true;
        run("cycle-dont-increase", contig, 150, 250, 50, 5, s, het(contig, moved, 101, 330));
        s = hc.copy();
        s.legacyCycles = true;
        run("cycle-legacy", contig, 150, 250, 50, 5, s, het(contig, moved, 101, 330));

        // Near the contig start: the reference padding is cut at base 1.
        run("contig-start", contig, 20, 80, 15, 10, hc, het(contig, snp(contig, 50, 'x'), 1, 120));
    }

    // ---- read construction ----

    /** `base` replaced at 1-based `pos` by the next base round ACGT (`x`). */
    static String snp(final String seq, final int pos, final char x) {
        final char old = seq.charAt(pos - 1);
        final char alt = "ACGT".charAt(("ACGT".indexOf(old) + 1) % 4);
        return seq.substring(0, pos - 1) + alt + seq.substring(pos);
    }

    static String insert(final String seq, final int afterPos, final String bases) {
        return seq.substring(0, afterPos) + bases + seq.substring(afterPos);
    }

    static String delete(final String seq, final int fromPos, final int length) {
        return seq.substring(0, fromPos - 1) + seq.substring(fromPos - 1 + length);
    }

    /**
     * Reads of length 50 every `step` bases over [from, to] of `hap`, each labelled with the start
     * it would have on the reference (the haplotype coordinate, which is close enough for a region
     * that only checks overlap and order).
     */
    static List<Read> tile(final String hap, final String sample, final int from, final int to, final int step, final int qual) {
        final List<Read> reads = new ArrayList<>();
        for (int start = from; start + 50 - 1 <= Math.min(to, hap.length()); start += step) {
            final String bases = hap.substring(start - 1, start - 1 + 50);
            final byte[] q = new byte[50];
            java.util.Arrays.fill(q, (byte) qual);
            reads.add(new Read(sample, start, "50M", bases, q));
        }
        return reads;
    }

    static List<Read> het(final String ref, final String alt) {
        return het(ref, alt, 101, 300);
    }

    static List<Read> het(final String ref, final String alt, final int from, final int to) {
        return concat(tile(ref, "s1", from, to, 4, 30), tile(alt, "s1", from + 2, to, 4, 30));
    }

    @SafeVarargs
    static List<Read> concat(final List<Read>... lists) {
        final List<Read> all = new ArrayList<>();
        for (final List<Read> l : lists) {
            all.addAll(l);
        }
        return all;
    }

    /** Every seventh read gets one base changed: errors seen once. */
    static List<Read> withErrors(final List<Read> reads) {
        final List<Read> out = new ArrayList<>();
        for (int i = 0; i < reads.size(); i++) {
            final Read r = reads.get(i);
            if (i % 7 == 3) {
                final int at = 10 + (i % 30);
                out.add(new Read(r.sample(), r.start(), r.cigar(), snp(r.bases(), at + 1, 'x'), r.quals()));
            } else {
                out.add(r);
            }
        }
        return out;
    }

    /** Every third read's first five bases become soft-clipped junk. */
    static List<Read> softClipped(final List<Read> reads) {
        final List<Read> out = new ArrayList<>();
        for (int i = 0; i < reads.size(); i++) {
            final Read r = reads.get(i);
            if (i % 3 == 0) {
                out.add(new Read(r.sample(), r.start() + 5, "5S45M", "TTTTT" + r.bases().substring(5), r.quals()));
            } else if (i % 3 == 1) {
                out.add(new Read(r.sample(), r.start(), "46M4S", r.bases().substring(0, 46) + "GGGG", r.quals()));
            } else {
                out.add(r);
            }
        }
        return out;
    }

    /** Every other read gets a low-quality stretch and an N. */
    static List<Read> lowQuality(final List<Read> reads) {
        final List<Read> out = new ArrayList<>();
        for (int i = 0; i < reads.size(); i++) {
            final Read r = reads.get(i);
            if (i % 2 == 0) {
                final byte[] q = r.quals().clone();
                for (int j = 20; j < 24; j++) {
                    q[j] = 5;
                }
                final String b = r.bases().substring(0, 35) + 'N' + r.bases().substring(36);
                out.add(new Read(r.sample(), r.start(), r.cigar(), b, q));
            } else {
                out.add(r);
            }
        }
        return out;
    }

    /** Reads ending at `pos` whose last three bases are changed, five copies. */
    static List<Read> dangling(final String ref, final int pos) {
        final List<Read> out = new ArrayList<>();
        for (int copy = 0; copy < 5; copy++) {
            final int start = pos - 49;
            String b = ref.substring(start - 1, pos);
            b = snp(snp(b, 48, 'x'), 50, 'x');
            final byte[] q = new byte[50];
            java.util.Arrays.fill(q, (byte) 30);
            out.add(new Read("s1", start, "50M", b, q));
        }
        return out;
    }

    /** Five reads ending at `pos` with one base changed at `at` (1-based in the read). */
    static List<Read> danglingTail(final String ref, final int pos, final int at) {
        final List<Read> out = new ArrayList<>();
        for (int copy = 0; copy < 5; copy++) {
            final int start = pos - 49;
            final String b = snp(ref.substring(start - 1, pos), at, 'x');
            final byte[] q = new byte[50];
            java.util.Arrays.fill(q, (byte) 30);
            out.add(new Read("s1", start, "50M", b, q));
        }
        return out;
    }

    /** Five reads starting at `pos` with one base changed at `at` (1-based in the read). */
    static List<Read> danglingHead(final String ref, final int pos, final int at) {
        final List<Read> out = new ArrayList<>();
        for (int copy = 0; copy < 5; copy++) {
            final String b = snp(ref.substring(pos - 1, pos - 1 + 50), at, 'x');
            final byte[] q = new byte[50];
            java.util.Arrays.fill(q, (byte) 30);
            out.add(new Read("s1", pos, "50M", b, q));
        }
        return out;
    }

    // ---- the run ----

    static void run(final String label, final String contigBases, final int activeStart, final int activeEnd,
                    final int padding, final int refPadding, final Settings s, final List<Read> unsorted) throws Exception {
        System.out.println("case\t" + label + "\t" + contigBases + "\t" + activeStart + "\t" + activeEnd + "\t" + padding + "\t" + refPadding);
        System.out.println("settings\t" + label + "\t" + s.render());
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord(CONTIG, contigBases.length()))));
        for (final String sample : List.of("s1", "s2")) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg-" + sample);
            group.setSample(sample);
            header.addReadGroup(group);
        }
        final List<Read> reads = new ArrayList<>(unsorted);
        reads.sort(Comparator.comparingInt(Read::start));
        final AssemblyRegion region = new AssemblyRegion(new SimpleInterval(CONTIG, activeStart, activeEnd), true, padding, header);
        int n = 0;
        for (final Read r : reads) {
            final String name = "r" + n++;
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, name, 0, r.start(),
                    r.bases().getBytes(StandardCharsets.US_ASCII), r.quals(), r.cigar());
            read.setReadGroup("rg-" + r.sample());
            if (!region.getPaddedSpan().overlaps(read)) {
                continue;
            }
            final StringBuilder q = new StringBuilder();
            for (final byte b : r.quals()) {
                q.append((char) (b + 33));
            }
            System.out.println("read\t" + label + "\t" + name + "\t" + r.sample() + "\t" + r.start() + "\t" + r.cigar() + "\t" + r.bases() + "\t" + q);
            region.add(read);
        }
        try {
            final SimpleInterval padded = region.getPaddedSpan();
            final SimpleInterval refLoc = new SimpleInterval(CONTIG, Math.max(padded.getStart() - refPadding, 1),
                    Math.min(padded.getEnd() + refPadding, contigBases.length()));
            final byte[] fullRef = contigBases.substring(refLoc.getStart() - 1, refLoc.getEnd()).getBytes(StandardCharsets.US_ASCII);
            final byte[] refBases = contigBases.substring(padded.getStart() - 1, padded.getEnd()).getBytes(StandardCharsets.US_ASCII);
            final ReadThreadingAssembler assembler = assembler(s);
            final SmithWatermanAligner aligner = SmithWatermanAligner.getAligner(SmithWatermanAligner.Implementation.JAVA);
            final SWParameters dangling = SmithWatermanAlignmentConstants.STANDARD_NGS;
            final SWParameters toRef = SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS;

            // What each k-mer size gives, through the package-private `assemble`.
            final List<GATKRead> clipped = region.getReads().stream().map(ReadClipper::hardClipSoftClippedBases).collect(Collectors.toList());
            final Method assemble = ReadThreadingAssembler.class.getDeclaredMethod("assemble", List.class, Haplotype.class,
                    SAMFileHeader.class, SmithWatermanAligner.class, SWParameters.class);
            assemble.setAccessible(true);
            try {
                @SuppressWarnings("unchecked")
                final List<AssemblyResult> results = (List<AssemblyResult>) assemble.invoke(assembler, clipped,
                        referenceHaplotype(region, refBases, refLoc), header, aligner, dangling);
                int i = 0;
                for (final AssemblyResult r : results) {
                    final boolean failed = r.getStatus() == AssemblyResult.Status.FAILED;
                    System.out.println("result\t" + label + "\t" + i++ + "\t" + (failed ? "-" : String.valueOf(r.getKmerSize()))
                            + "\t" + r.getStatus()
                            + "\t" + (r.getSeqGraph() == null ? "-" : String.valueOf(r.getSeqGraph().vertexSet().size()))
                            + "\t" + (r.getSeqGraph() == null ? "-" : String.valueOf(r.getSeqGraph().edgeSet().size())));
                }
            } catch (final java.lang.reflect.InvocationTargetException e) {
                System.out.println("error\t" + label + "\tassemble\t" + e.getCause().getClass().getSimpleName() + ": " + e.getCause().getMessage());
            }

            final AssemblyResultSet set = assembler.runLocalAssembly(region, referenceHaplotype(region, refBases, refLoc), fullRef, refLoc,
                    null, header, aligner, null, dangling, toRef);
            final List<Haplotype> haplotypes = set.getHaplotypeList();
            int refIndex = -1;
            for (int i = 0; i < haplotypes.size(); i++) {
                final Haplotype h = haplotypes.get(i);
                if (h == set.getReferenceHaplotype()) {
                    refIndex = i;
                }
                System.out.println("haplotype\t" + label + "\t" + i + "\t" + h.isReference() + "\t" + h.getBaseString()
                        + "\t" + h.getCigar() + "\t" + h.getAlignmentStartHapwrtRef()
                        + "\t" + Long.toHexString(Double.doubleToRawLongBits(h.getScore())) + "," + h.getScore()
                        + "\t" + h.getKmerSize() + "\t" + h.getGenomeLocation());
            }
            System.out.println("set\t" + label + "\t" + haplotypes.size() + "\t" + set.isVariationPresent() + "\t" + refIndex);
        } catch (final Exception e) {
            System.out.println("error\t" + label + "\trun\t" + e.getClass().getSimpleName() + ": " + e.getMessage());
        }
    }

    /** `ReferenceConfidenceModel.createReferenceHaplotype`. */
    static Haplotype referenceHaplotype(final AssemblyRegion region, final byte[] refBases, final SimpleInterval refLoc) {
        final Haplotype h = new Haplotype(refBases, true);
        h.setGenomeLocation(region.getPaddedSpan());
        h.setAlignmentStartHapwrtRef(region.getPaddedSpan().getStart() - refLoc.getStart());
        final Cigar c = TextCigarCodec.decode(refBases.length + "M");
        h.setCigar(c);
        return h;
    }

    /** `HaplotypeCallerReadThreadingAssemblerArgumentCollection.makeReadThreadingAssembler`. */
    static ReadThreadingAssembler assembler(final Settings s) {
        final ReadThreadingAssembler a = new ReadThreadingAssembler(s.maxPaths, s.kmers, s.dontIncrease, s.allowNonUnique,
                s.pruningSamples, s.pruneFactor, s.adaptive, s.initialErrorRate, s.pruningLogOdds, s.seedingLogOdds,
                s.maxUnprunedVariants, false, s.legacyCycles, s.minMatching);
        a.setRecoverDanglingBranches(s.recoverDangling);
        a.setRecoverAllDanglingBranches(s.recoverAll);
        a.setMinDanglingBranchLength(s.minDanglingLength);
        a.setMinBaseQualityToUseInAssembly((byte) s.minBaseQuality);
        return a;
    }
}
