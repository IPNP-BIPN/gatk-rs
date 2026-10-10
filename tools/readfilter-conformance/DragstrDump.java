/*
 * DRAGstr in the read likelihoods: the read STR analyzer, the parameter tables, and the two
 * PairHMM input score imputators the default one does not cover.
 *
 * `DragstrReadSTRAnalyzer` gives, at each read position, the number of repeats of every period up
 * to the maximum and the period with the most (the smallest on a tie). `DragstrParams` looks GOP,
 * GCP and API up by period and repeat count, clamping both to the table; `DragstrParamUtils.parse`
 * reads the table file `CalibrateDragstrModel` writes. `DragstrPairHMMInputScoreImputator` gives
 * each base but the last the GOP (capped at 40) and GCP of its most repeated period, rounded, and
 * the last 45 and 10; `NonSymmetricalPairHMMInputScoreImputator` gives flat insertion, deletion
 * and continuation penalties. The last section runs `PairHMMLikelihoodCalculationEngine` with the
 * default DRAGstr parameters.
 *
 * Output:
 *
 *     str\t<maxPeriod>\t<bases>\t<pos>\t<most period>\t<most repeats>\t<repeats for period 1..max>
 *     param\t<table label>\t<gop|gcp|api>\t<period>\t<repeats>\t<bits>
 *     paramtext\t<label>\t<file text, lines joined by |>
 *     parsed\t<label>\t<maxPeriod>\t<maxRepeats>          (then param rows under that label)
 *     impute\t<table label>\t<bases>\t<gop,...>\t<gcp,...>
 *     nonsym\t<gcp>\t<ins>\t<del>\t<quals length>\t<read length>\t<del,...>\t<ins,...>\t<gcp,...>
 *     haplotype\t<index>\t<bases>
 *     read\t<name>\t<bases>\t<quals phred+33>\t<mapq>
 *     lk\t<allele>\t<read>\t<bits>
 *     error\t<label>\t<exception class>: <message>
 *
 * Usage: DragstrDump
 */

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.DragstrPairHMMInputScoreImputator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.NonSymmetricalPairHMMInputScoreImputator;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.PairHMMLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParamUtils;
import org.broadinstitute.hellbender.utils.dragstr.DragstrParams;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.pairhmm.DragstrReadSTRAnalyzer;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMMInputScoreImputation;
import org.broadinstitute.hellbender.utils.read.ArtificialReadUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class DragstrDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# DragstrDump: DragstrReadSTRAnalyzer, DragstrParams, the DRAGstr and non-symmetrical imputators");
        final Random random = new Random(20261013L);
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 60; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        final List<String> sequences = List.of("A", "AC", "AAAA", "ACACACAC", "AACAACAAC", "ACGTACGTACGTT",
                "TTTTTACACACGGGGCAGCAGCAGT", "ACGACGACGAC", "G".repeat(27), "ATATATATATCGCGCGCGAAAAAAAAAAAAAAAAAAAAAAAA",
                "AAGAAGAAGAAGTTTTGTTTTGTTTTG", sb.toString());
        final Method of = DragstrReadSTRAnalyzer.class.getDeclaredMethod("of", byte[].class, int.class);
        of.setAccessible(true);
        for (final int maxPeriod : new int[]{1, 2, 4, 8}) {
            for (final String s : sequences) {
                final DragstrReadSTRAnalyzer a = (DragstrReadSTRAnalyzer) of.invoke(null, s.getBytes(StandardCharsets.US_ASCII), maxPeriod);
                for (int pos = 0; pos < s.length(); pos++) {
                    final StringBuilder r = new StringBuilder();
                    for (int p = 1; p <= maxPeriod; p++) {
                        r.append(p == 1 ? "" : ",").append(a.numberOfRepeats(pos, p));
                    }
                    System.out.println("str\t" + maxPeriod + "\t" + s + "\t" + pos + "\t" + a.mostRepeatedPeriod(pos) + "\t"
                            + a.numberOfMostRepeats(pos) + "\t" + r);
                }
            }
        }

        lookups("default", DragstrParams.DEFAULT);

        final Path dir = Files.createTempDirectory("dragstr");
        final Path printed = dir.resolve("default.txt");
        DragstrParamUtils.print(DragstrParams.DEFAULT, new GATKPath(printed.toString()), "origin", "dump");
        final String defaultText = Files.readString(printed);
        final Map<String, String> texts = new LinkedHashMap<>();
        texts.put("printed-default", defaultText);
        texts.put("small", "# a comment\n    1      2      3\nGOP:\n10.00  20.00  30.00\n11.00  21.00  31.50\nGCP:\n5.00  5.00  5.00\n2.50  2.50  2.50\nAPI:\n1.00 2.00 3.00\n4.00 5.00 6.00\n");
        texts.put("extra-table", "1 2\nGOP:\n10 20\nXYZ:\n1 1\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("bad-header", "1 3\nGOP:\n10 20\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("missing-api", "1 2\nGOP:\n10 20\nGCP:\n5 5\n");
        texts.put("too-few-columns", "1 2\nGOP:\n10\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("too-many-columns", "1 2\nGOP:\n10 20 30\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("negative", "1 2\nGOP:\n10 -20\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("not-a-number", "1 2\nGOP:\n10 x\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("rows-differ", "1 2\nGOP:\n10 20\n11 21\nGCP:\n5 5\nAPI:\n1 2\n");
        texts.put("only-comments", "# nothing\n# here\n");
        texts.put("header-only", "1 2\n");
        texts.put("leading-spaces", "  1  2\nGOP:\n  10  20\nGCP:\n  5  5\nAPI:\n  1  2\n");
        final Map<String, DragstrParams> parsed = new LinkedHashMap<>();
        for (final Map.Entry<String, String> e : texts.entrySet()) {
            System.out.println("paramtext\t" + e.getKey() + "\t" + e.getValue().replace("\n", "|"));
            final Path f = dir.resolve(e.getKey() + ".txt");
            Files.writeString(f, e.getValue());
            try {
                final DragstrParams p = DragstrParamUtils.parse(new GATKPath(f.toString()));
                System.out.println("parsed\t" + e.getKey() + "\t" + p.maximumPeriod() + "\t" + p.maximumRepeats());
                lookups(e.getKey(), p);
                parsed.put(e.getKey(), p);
            } catch (final Exception ex) {
                // The temporary directory is not reproducible; the file name is.
                System.out.println("error\t" + e.getKey() + "\t" + ex.getClass().getSimpleName() + ": "
                        + String.valueOf(ex.getMessage()).replace(dir.toString() + "/", ""));
            }
        }

        final List<String> reads = List.of("ACGTTGCA", "AAAAAAAAAAAAAAAAAAAAAAAAACG", "CACACACACACACAGT", "T",
                "GTTGTTGTTGTTGTTGTTGTTA", sb.toString());
        for (final String label : List.of("default", "small")) {
            final DragstrParams p = label.equals("default") ? DragstrParams.DEFAULT : parsed.get(label);
            for (final String bases : reads) {
                final GATKRead read = ArtificialReadUtils.createArtificialRead(bases.getBytes(StandardCharsets.US_ASCII),
                        new byte[bases.length()], bases.length() + "M");
                final PairHMMInputScoreImputation imp = DragstrPairHMMInputScoreImputator.of(p).impute(read);
                System.out.println("impute\t" + label + "\t" + bases + "\t" + show(imp.insOpenPenalties()) + "\t" + show(imp.gapContinuationPenalties()));
            }
        }
        for (final String bases : List.of("ACGT", "A")) {
            final GATKRead read = ArtificialReadUtils.createArtificialRead(bases.getBytes(StandardCharsets.US_ASCII),
                    new byte[bases.length()], bases.length() + "M");
            final PairHMMInputScoreImputation imp = NonSymmetricalPairHMMInputScoreImputator.newInstance((byte) 10, (byte) 30, (byte) 35).impute(read);
            System.out.println("nonsym\t10\t30\t35\t" + read.getBaseQualityCount() + "\t" + read.getLength() + "\t"
                    + show(imp.delOpenPenalties()) + "\t" + show(imp.insOpenPenalties()) + "\t" + show(imp.gapContinuationPenalties()));
        }

        // The likelihood engine with the default DRAGstr parameters.
        final String ref = "ACGTACGTTTTTTTTTTGCACACACACAGGCTAGCTAGGATC";
        final List<String> haps = List.of(ref, ref.substring(0, 10) + "TT" + ref.substring(10), ref.substring(0, 20) + ref.substring(24),
                ref.substring(0, 30) + "A" + ref.substring(31));
        final List<Haplotype> haplotypes = new ArrayList<>();
        for (int i = 0; i < haps.size(); i++) {
            System.out.println("haplotype\t" + i + "\t" + haps.get(i));
            haplotypes.add(new Haplotype(haps.get(i).getBytes(StandardCharsets.US_ASCII), i == 0));
        }
        final SAMFileHeader header = ArtificialReadUtils.createArtificialSamHeader();
        final List<GATKRead> engineReads = new ArrayList<>();
        final String[] readBases = {ref.substring(2, 32), haps.get(1).substring(4, 34), haps.get(2).substring(0, 30),
                haps.get(3).substring(10, 40), ref.substring(5, 25)};
        for (int i = 0; i < readBases.length; i++) {
            final byte[] q = new byte[readBases[i].length()];
            for (int j = 0; j < q.length; j++) {
                q[j] = (byte) (20 + (i * 7 + j * 3) % 20);
            }
            final GATKRead read = ArtificialReadUtils.createArtificialRead(header, "e" + i, 0, 1,
                    readBases[i].getBytes(StandardCharsets.US_ASCII), q, readBases[i].length() + "M");
            read.setMappingQuality(60);
            engineReads.add(read);
            final StringBuilder qs = new StringBuilder();
            for (final byte b : q) {
                qs.append((char) (b + 33));
            }
            System.out.println("read\te" + i + "\t" + readBases[i] + "\t" + qs + "\t60");
        }
        final PairHMMLikelihoodCalculationEngine engine = new PairHMMLikelihoodCalculationEngine((byte) 10, DragstrParams.DEFAULT,
                null, PairHMM.Implementation.LOGLESS_CACHING, -4.5, PairHMMLikelihoodCalculationEngine.PCRErrorModel.NONE);
        final Map<String, List<GATKRead>> perSample = new LinkedHashMap<>();
        perSample.put("s1", engineReads);
        final AlleleLikelihoods<GATKRead, Haplotype> result = engine.computeReadLikelihoods(haplotypes, header,
                new IndexedSampleList(List.of("s1")), perSample, true);
        for (int a = 0; a < result.numberOfAlleles(); a++) {
            for (int r = 0; r < result.sampleEvidenceCount(0); r++) {
                System.out.println("lk\t" + a + "\t" + result.sampleEvidence(0).get(r).getName() + "\t"
                        + Long.toHexString(Double.doubleToRawLongBits(result.sampleMatrix(0).get(a, r))));
            }
        }
    }

    static String show(final byte[] b) {
        final StringBuilder s = new StringBuilder();
        for (int i = 0; i < b.length; i++) {
            s.append(i == 0 ? "" : ",").append(b[i]);
        }
        return s.toString();
    }

    static void lookups(final String label, final DragstrParams p) {
        for (int period = 1; period <= p.maximumPeriod() + 1; period++) {
            for (int repeats = 1; repeats <= p.maximumRepeats() + 1; repeats++) {
                System.out.println("param\t" + label + "\tgop\t" + period + "\t" + repeats + "\t" + Long.toHexString(Double.doubleToRawLongBits(p.gop(period, repeats))));
                System.out.println("param\t" + label + "\tgcp\t" + period + "\t" + repeats + "\t" + Long.toHexString(Double.doubleToRawLongBits(p.gcp(period, repeats))));
                System.out.println("param\t" + label + "\tapi\t" + period + "\t" + repeats + "\t" + Long.toHexString(Double.doubleToRawLongBits(p.api(period, repeats))));
            }
        }
    }
}
