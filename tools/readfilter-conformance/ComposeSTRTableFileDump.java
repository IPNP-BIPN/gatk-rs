/*
 * ComposeSTRTableFile's output, taken from the reference.
 *
 * A reference scanned for short tandem repeats, every site reported with the period and the number
 * of repeats that fit it best. The output is a zip; what is printed here is its sites table and its
 * summary.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - EVERY POSITION IS TRIED AS THE START OF A REPEAT UNIT, and when one is found THE SCAN JUMPS
 *     TO ITS END, so no position is looked at twice;
 *   - BUT THE SITES STILL OVERLAP, because the search reaches BACKWARDS from the position it
 *     started at: a homopolymer ending at 9 and a dinucleotide repeat beginning at 9 are both
 *     reported, sharing that base;
 *   - THE BEST PERIOD IS THE ONE WITH THE MOST REPEATS, ties going to the SHORTER period, and the
 *     repeat count is an INTEGER DIVISION of the span by the period, so trailing bases that do not
 *     complete a unit are inside the interval and count for nothing;
 *   - PERIOD ONE IS SPECIAL-CASED and is the only period tried when the base at the position is
 *     not ACGT, in which case NOTHING is emitted at all, which is why an N leaves a gap;
 *   - A PERIOD IS ABANDONED AS SOON AS ITS UNIT CONTAINS A NON-STANDARD BASE;
 *   - THE MASK IS A PER-CONTIG COUNTER PER PERIOD AND REPEAT, initialised to the CONTIG'S INDEX
 *     rather than to zero, so the first site of the second contig carries mask 1;
 *   - DECIMATION IS A BIT TEST between that mask and the table's entry for the period and repeat,
 *     so it keeps one site in every 2^n rather than a fraction, and under the default table the
 *     second contig's homopolymer is the one that goes;
 *   - --max-repeat CAPS THE INDEX THE MASK COUNTER IS KEPT UNDER, not the repeat reported, so
 *     lowering it makes distinct sites share a counter and CHANGES THE MASKS, and therefore which
 *     sites decimation removes;
 *   - --max-period CAPS THE SEARCH, so a region whose true period is longer breaks into whatever
 *     shorter period won, one site per base when that is period one with a single repeat;
 *   - -L RESTRICTS WHERE THE SCAN STARTS BUT NOT HOW FAR A SITE REACHES, so a site can begin before
 *     the interval and end after it;
 *   - AND THE SUMMARY COUNTS BOTH WHAT WAS EMITTED AND WHAT WAS DECIMATED, per period and repeat.
 *
 * Output:
 *
 *     fixture\tchr1=<the sequence>
 *     sites\t<label>=<the whole sites table, escaped>
 *     summary\t<label>=<the whole summary, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: ComposeSTRTableFileDump
 */

import org.broadinstitute.hellbender.tools.dragstr.ComposeSTRTableFile;

import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeMap;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class ComposeSTRTableFileDump {

    /**
     * One contig carrying a homopolymer, a dinucleotide repeat, a trinucleotide repeat, an N, a
     * second homopolymer, a four-base repeat and a short tail.
     */
    static final String CHR1 = "AAAAAAAA"
            + "ACACACACACACACAC"
            + "GTCGTCGTCGTC"
            + "N"
            + "TTTTTTTT"
            + "ACGTACGTACGTACG"
            + "CCCCC";

    /** A second contig, so the mask's per-contig initial value is visible. */
    static final String CHR2 = "GGGGGGGGTATATATATA";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("compose-str-table-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# ComposeSTRTableFileDump: a reference scanned for short tandem repeats");
        System.out.printf("fixture\tchr1=%s%n", CHR1);
        System.out.printf("fixture\tchr2=%s%n", CHR2);

        final Path fasta = writeReference(dir);

        run(dir, "default", fasta, "DEFAULT", 8, 20);
        run(dir, "no-decimation", fasta, "NONE", 8, 20);
        run(dir, "max-period-two", fasta, "NONE", 2, 20);
        run(dir, "max-repeat-three", fasta, "NONE", 8, 3);
        run(dir, "interval", fasta, "NONE", 8, 20, "-L", "chr1:20-40");
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("str.fasta");
        Files.writeString(fasta, ">chr1\n" + CHR1 + "\n>chr2\n" + CHR2 + "\n",
                StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[] {
                "R=" + fasta, "O=" + dir.resolve("str.dict")});
        return fasta;
    }

    static void run(final Path dir, final String label, final Path fasta, final String decimation,
                    final int maxPeriod, final int maxRepeat, final String... extra)
            throws Exception {
        final Path out = dir.resolve(label + ".zip");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "-R", fasta.toString(),
                "-O", out.toString(),
                "--decimation", decimation,
                "--max-period", Integer.toString(maxPeriod),
                "--max-repeats", Integer.toString(maxRepeat),
                "--generate-sites-text-output", "true"));
        argv.addAll(Arrays.asList(extra));
        try {
            new ComposeSTRTableFile().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
            return;
        }
        final Map<String, String> entries = unzip(out);
        for (final String name : List.of("sites.txt", "summary.txt")) {
            final String content = entries.get(name);
            if (content != null) {
                System.out.printf("%s\t%s=%s%n", name.replace(".txt", ""), label,
                        ReferenceQueryDump.escape(masked(content, dir)));
            }
        }
    }

    /** Every entry of the zip, keyed by name. */
    static Map<String, String> unzip(final Path archive) throws Exception {
        final Map<String, String> entries = new TreeMap<>();
        try (final InputStream raw = Files.newInputStream(archive);
             final ZipInputStream zip = new ZipInputStream(raw)) {
            ZipEntry entry;
            while ((entry = zip.getNextEntry()) != null) {
                if (entry.isDirectory()) {
                    continue;
                }
                final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
                zip.transferTo(bytes);
                entries.put(entry.getName(), bytes.toString(StandardCharsets.UTF_8));
            }
        }
        return entries;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
