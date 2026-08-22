/*
 * MTLowHeteroplasmyFilterTool's output, taken from the reference.
 *
 * Low-heteroplasmy mitochondrial calls filtered, but only once there are too many of them. It is a
 * TwoPassVariantWalker, and the two passes do entirely different things.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE FILTER IS ALL OR NOTHING ACROSS THE WHOLE FILE. The first pass counts the UNFILTERED low
 *     heteroplasmy sites; the second filters every low allele it can find, but ONLY if that count
 *     exceeded --max-allowed-low-hets. So a file with three such sites comes out untouched and a
 *     file with four comes out with all four filtered, and one record's fate is decided by the
 *     others;
 *   - THE COUNT IS OF SITES, NOT OF ALLELES, so a multiallelic record with two low alternates
 *     counts once;
 *   - A SITE THAT IS ALREADY FILTERED DOES NOT COUNT, `variant.isNotFiltered()` being the guard, so
 *     adding a filter upstream can keep the whole file from being filtered here;
 *   - AND `PASS` IS NOT FILTERED EITHER, since htsjdk's `isFiltered` asks whether the filter SET is
 *     non-empty and `PASS` leaves it empty, so a record marked PASS counts exactly like a record
 *     marked `.`;
 *   - THE THRESHOLD IS STRICT, `x < lowHetThreshold`, so an allele fraction of exactly 0.1 is not
 *     low at the default;
 *   - AND --low-het-threshold DOES NOTHING. The field is declared `private final double
 *     lowHetThreshold = 0.1`, which is a compile-time constant, so javac replaces every read of it
 *     with the literal and whatever Barclay writes into the field is never looked at again. The
 *     control is `maxAllowedLowHets`, declared without `final` in the same class, which works:
 *     a run at --low-het-threshold 0.6 over fractions of 0.05 and 0.5 filters only the 0.05 ones,
 *     which is what the default threshold does;
 *   - THE FIRST PASS READS EVERY GENOTYPE WITHOUT A PRECONDITION while the second reads only those
 *     that HAVE an allele fraction, so a genotype with no AF at all is asked for one anyway;
 *   - `AF=.` IS NOT A MISSING VALUE INSIDE THE ARRAY, IT IS AN ABSENT ATTRIBUTE, so it takes the
 *     `() -> null` default and the first pass throws a NullPointerException on it, exactly as it
 *     does for a genotype with no AF field at all. The Double.MAX_VALUE substitution is for a
 *     missing entry WITHIN a multi-valued array, `AF=0.5,.`, which is measured separately and is
 *     never below any threshold;
 *   - THE SECOND PASS TAKES THE MAXIMUM ACROSS SAMPLES per alternate allele;
 *   - AND WRITING THE AS_FilterStatus ATTRIBUTE THROWS when the record has none to merge into,
 *     exactly as NuMTFilterTool does.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     filtered\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MTLowHeteroplasmyDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.MTLowHeteroplasmyFilterTool;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MTLowHeteroplasmyDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##FILTER=<ID=weak_evidence,Description=\"Weak evidence\">\n"
            + "##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=\"Filter status for each allele\">\n"
            + "##contig=<ID=chrM,length=16569>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n";

    /** One biallelic record: position, FILTER column, and the two samples' allele fractions. */
    static String record(final int position, final String filter, final String first,
                         final String second) {
        return "chrM\t" + position + "\t.\tA\tC\t50\t" + filter
                + "\tAS_FilterStatus=SITE\tGT:AF\t0/1:" + first + "\t0/1:" + second + "\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mt-low-heteroplasmy-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MTLowHeteroplasmyDump: low-heteroplasmy calls filtered, but only in bulk");

        // Exactly three unfiltered low sites, which is the default allowance and not one more.
        final String three =
                record(100, ".", "0.05", "0.05")
                + record(200, ".", "0.05", "0.05")
                + record(300, ".", "0.05", "0.05")
                + record(400, ".", "0.5", "0.5");
        // A fourth, which trips the whole file.
        final String four = three + record(500, ".", "0.05", "0.05");
        run(dir, "three-low-hets", HEADER + three);
        run(dir, "four-low-hets", HEADER + four);
        // The same four, one of them already filtered, so only three count.
        run(dir, "one-already-filtered", HEADER + three + record(500, "weak_evidence", "0.05", "0.05"));
        // And the same four with PASS rather than a dot, which counts exactly the same.
        run(dir, "pass-not-dot", HEADER
                + record(100, "PASS", "0.05", "0.05") + record(200, "PASS", "0.05", "0.05")
                + record(300, "PASS", "0.05", "0.05") + record(400, "PASS", "0.5", "0.5")
                + record(500, "PASS", "0.05", "0.05"));
        // An allowance of zero, where a single low site is enough.
        run(dir, "allow-none", HEADER + record(100, ".", "0.05", "0.05") + record(200, ".", "0.5", "0.5"),
                "--max-allowed-low-hets", "0");
        // A threshold that moves which sites are low, over the same file.
        run(dir, "threshold-0.6", HEADER + three, "--low-het-threshold", "0.6",
                "--max-allowed-low-hets", "0");
        // Exactly the threshold, which is not below it.
        run(dir, "exactly-the-threshold", HEADER
                + record(100, ".", "0.1", "0.1") + record(200, ".", "0.1", "0.1")
                + record(300, ".", "0.1", "0.1") + record(400, ".", "0.1", "0.1"));
        // An AF that is a dot, which is an ABSENT attribute rather than a missing value.
        run(dir, "af-is-a-dot", HEADER
                + record(100, ".", ".", ".") + record(200, ".", ".", ".")
                + record(300, ".", ".", ".") + record(400, ".", ".", "."),
                "--max-allowed-low-hets", "0");
        // A missing value INSIDE a multi-valued array, which is the Double.MAX_VALUE substitution.
        run(dir, "af-entry-is-a-dot", HEADER
                + "chrM\t100\t.\tA\tC,G\t50\t.\tAS_FilterStatus=SITE|SITE\tGT:AF\t0/1:0.05,.\t0/1:0.05,.\n",
                "--max-allowed-low-hets", "0");
        // A multiallelic record, whose two alternates are one low and one not.
        run(dir, "multiallelic", HEADER
                + "chrM\t100\t.\tA\tC,G\t50\t.\tAS_FilterStatus=SITE|SITE\tGT:AF\t0/1:0.5,0.05\t0/1:0.5,0.05\n",
                "--max-allowed-low-hets", "0");
        // A genotype with no allele fraction at all, which the first pass asks for anyway.
        run(dir, "genotype-without-af", HEADER
                + "chrM\t100\t.\tA\tC\t50\t.\tAS_FilterStatus=SITE\tGT\t0/1\t0/1\n");
        // And a record with nothing to merge the attribute into.
        run(dir, "no-as-filter-status", HEADER
                + "chrM\t100\t.\tA\tC\t50\t.\t.\tGT:AF\t0/1:0.05\t0/1:0.05\n",
                "--max-allowed-low-hets", "0");
    }

    static void run(final Path dir, final String label, final String vcf, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, vcf, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(vcf));

        final Path out = dir.resolve(label + "-filtered.vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-V", in.toString(), "-O", out.toString()));
        argv.addAll(Arrays.asList(extra));
        try {
            new MTLowHeteroplasmyFilterTool().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        if (Files.exists(out)) {
            System.out.printf("filtered\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>")
                .replaceAll("##GATKCommandLine=<[^\n]*>\n", "");
    }
}
