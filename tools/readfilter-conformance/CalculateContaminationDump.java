/*
 * CalculateContamination's two files, taken from the reference.
 *
 * The model underneath is measured already (`contamination-model`); what this dump is for is the
 * tool around it: which sites survive the coverage filter, which table the genotyping is done
 * from, and what reaches the two output files.
 *
 * Nine behaviours this is built to catch.
 *
 *   - --low-coverage-ratio-threshold AND --high-coverage-ratio-threshold ARE ACCEPTED AND IGNORED.
 *     Barclay does set the fields, and a probe reading them back after the parse shows the value
 *     that was asked for; the filter never sees it, because both fields are `private final double`
 *     initialised from a constant expression, which makes them constant variables under JLS 4.12.4
 *     and their one read a compile-time constant. A low ratio of ten, which would drop every site
 *     in the table, changes nothing;
 *   - THE COVERAGE FILTER STILL CUTS, at a ceiling of three times the MEAN: six homozygous-
 *     alternate sites at depth four hundred are dropped, and the same six at depth sixty are kept
 *     and move the answer, which is how the golden shows the cut happened rather than assuming it;
 *   - A SITE AT OR UNDER `MIN_COVERAGE` IS DROPPED BEFORE EITHER STATISTIC IS TAKEN, so it moves
 *     neither the median nor the mean;
 *   - WITHOUT A MATCHED NORMAL THE TUMOUR IS GENOTYPED FROM ITSELF;
 *   - WITH ONE THE NORMAL SAYS WHICH SITES ARE HOMOZYGOUS AND THE TUMOUR'S OWN COUNTS ARE READ
 *     THERE: a normal whose homozygous-alternate sites sit at other positions answers 0.915 where
 *     the tumour alone answers 0.044;
 *   - A NORMAL WITH NO HOMOZYGOUS-ALTERNATE SITE AT ALL FALLS THROUGH to the hom-ref strategies,
 *     which read the tumour and not the model, so the answer returns to the tumour-only one;
 *   - --tumor-segmentation IS OPTIONAL and changes neither number in the contamination table;
 *   - THE SEGMENTATION IS THE TUMOUR'S OWN even when a matched normal did the genotyping, which is
 *     the one place the two models are built from different site lists;
 *   - AND THE SAMPLE NAME IS THE TUMOUR TABLE'S METADATA LINE, copied into both outputs, while the
 *     tool returns the string `SUCCESS`.
 *
 * Output:
 *
 *     table\t<label>=<that input file, escaped>
 *     contamination\t<label>\t<the contamination table, escaped>
 *     segments\t<label>\t<the segmentation table, escaped, or `absent`>
 *     returned\t<label>\t<what doWork returned>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: CalculateContaminationDump
 */

import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CalculateContaminationDump {

    /** One row of a pileup summary table. */
    record Site(int position, int ref, int alt, int other, double frequency) {}

    static String table(final String sample, final List<Site> sites) {
        final StringBuilder text = new StringBuilder("#<METADATA>SAMPLE=" + sample + "\n"
                + "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n");
        for (final Site site : sites) {
            text.append("chr1\t").append(site.position()).append('\t').append(site.ref())
                    .append('\t').append(site.alt()).append('\t').append(site.other())
                    .append('\t').append(String.format("%.6f", site.frequency())).append('\n');
        }
        return text.toString();
    }

    /**
     * A stretch of sites in the shape the model can fit.
     *
     * Every fifth site is homozygous reference with a little error, every fifth homozygous
     * alternate, and the rest heterozygous, which is what `ContaminationModel` needs to fit an
     * allele fraction before it can read a contamination off the homozygous-alternate ones. The
     * shape is `ContaminationModelDump`'s own, so the two suites are looking at the same fixture
     * from either side of the tool.
     */
    static List<Site> sites(final int count, final int depthOffset, final boolean loh,
                            final int shift) {
        final List<Site> sites = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            final int position = 1000 + 100 * i;
            final double frequency = 0.05 + (i % 19) * 0.05;
            final int depth = 50 + depthOffset + (i % 7) * 5;
            final int other = i % 4 == 0 ? 1 : 0;
            final int alt;
            switch ((i + shift) % 5) {
                case 0:
                    alt = 1 + i % 2;
                    break;
                case 3:
                    alt = depth - other - (i % 2);
                    break;
                default:
                    alt = (loh && i >= count / 2) ? (depth - other) / 4 : (depth - other) / 2;
                    break;
            }
            sites.add(new Site(position, depth - alt - other, alt, other, frequency));
        }
        return sites;
    }

    /** A stretch of sites at ONE depth, which is what moves the coverage statistics. */
    static List<Site> flat(final int from, final int count, final int depth) {
        final List<Site> sites = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            sites.add(new Site(from + i * 100, depth / 2, depth / 2, 0, 0.5));
        }
        return sites;
    }

    /**
     * A stretch of homozygous-alternate sites carrying a fifth of their reads as reference.
     *
     * A hom-alt site is what the contamination is read off, so these WOULD move the answer if the
     * coverage filter kept them: putting them at a depth the high threshold cuts is how the golden
     * shows the cut happened at all.
     */
    static List<Site> deepHomAlt(final int from, final int count, final int depth) {
        final List<Site> sites = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            sites.add(new Site(from + i * 100, depth / 5, depth - depth / 5, 0, 0.5));
        }
        return sites;
    }

    static List<Site> concat(final List<Site>... parts) {
        final List<Site> all = new ArrayList<>();
        for (final List<Site> part : parts) {
            all.addAll(part);
        }
        return all;
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("calculate-contamination-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# CalculateContaminationDump: the coverage filter, the genotyping "
                + "table, and the two files");

        // A hundred sites in the shape the model fits: hom ref, hom alt and het in turn.
        final List<Site> tumour = sites(100, 0, false, 0);
        final Path tumourPath = write(dir, "tumour.table", table("tumour", tumour), "tumour");

        // The same sites with a handful under the minimum coverage, which the filter drops before
        // it computes either statistic.
        final List<Site> withUncovered = concat(tumour, flat(90000, 5, 8));
        final Path uncoveredPath = write(dir, "uncovered.table",
                table("tumour", withUncovered), "with-uncovered-sites");

        // The same sites with a tail of very deep ones, which moves the mean and so the top.
        final List<Site> withDeep = concat(tumour, deepHomAlt(95000, 6, 400));
        final Path deepPath = write(dir, "deep.table", table("tumour", withDeep),
                "with-a-deep-tail");

        // The same six sites at a depth the filter keeps, which is what says the answer moves when
        // they are not cut.
        final List<Site> withShallow = concat(tumour, deepHomAlt(95000, 6, 60));
        final Path shallowPath = write(dir, "shallow.table", table("tumour", withShallow),
                "with-a-shallow-tail");

        // A matched normal, deeper than the tumour and with a loss of heterozygosity in its second
        // half, so the two models see different site lists.
        // Its homozygous-alternate sites sit at DIFFERENT positions, which is what makes the
        // genotyping source visible: the contamination is still read off the tumour's counts, but
        // at the positions the normal called homozygous.
        final List<Site> normal = sites(100, 10, true, 2);
        final Path normalPath = write(dir, "normal.table", table("normal", normal), "normal");

        // A normal in which nothing is homozygous alternate, so the genotyping model has nothing
        // to read a contamination off.
        final List<Site> homRef = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            homRef.add(new Site(1000 + 100 * i, 59, 1, 0, 0.05 + (i % 19) * 0.05));
        }
        final Path homRefPath = write(dir, "homref.table", table("normal", homRef),
                "normal-all-hom-ref");

        run(dir, "plain", tumourPath, null, false, List.of());
        run(dir, "with-segmentation", tumourPath, null, true, List.of());
        run(dir, "uncovered-sites", uncoveredPath, null, false, List.of());
        run(dir, "a-deep-tail", deepPath, null, false, List.of());
        run(dir, "a-shallow-tail", shallowPath, null, false, List.of());
        // The two ratios, including a low one of nought.
        run(dir, "low-ratio-zero", deepPath, null, false,
                List.of("--low-coverage-ratio-threshold", "0.0"));
        run(dir, "high-ratio-one", tumourPath, null, false,
                List.of("--high-coverage-ratio-threshold", "1.0"));
        run(dir, "low-ratio-one", tumourPath, null, false,
                List.of("--low-coverage-ratio-threshold", "1.0"));
        run(dir, "low-ratio-ten", tumourPath, null, false,
                List.of("--low-coverage-ratio-threshold", "10.0"));
        // The matched normal, with and without the segmentation the tumour's own model writes.
        run(dir, "matched-normal", tumourPath, normalPath, false, List.of());
        run(dir, "matched-normal-homref", tumourPath, homRefPath, false, List.of());
        run(dir, "matched-normal-segmented", tumourPath, normalPath, true, List.of());
    }

    static Path write(final Path dir, final String name, final String text, final String label)
            throws Exception {
        System.out.printf("table\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path input, final Path matched,
                    final boolean segmentation, final List<String> extra) throws Exception {
        final Path out = dir.resolve("contamination-" + label + ".table");
        final Path segments = dir.resolve("segments-" + label + ".table");
        final List<String> argv = new ArrayList<>(List.of("-I", input.toString(),
                "-O", out.toString()));
        if (matched != null) {
            argv.add("--matched-normal");
            argv.add(matched.toString());
        }
        if (segmentation) {
            argv.add("--tumor-segmentation");
            argv.add(segments.toString());
        }
        argv.addAll(extra);
        final Object returned;
        try {
            returned = new CalculateContamination().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(cause.getMessage())
                            .replace(dir.toString(), "<dir>")));
            return;
        }
        System.out.printf("returned\t%s\t%s%n", label, returned);
        System.out.printf("contamination\t%s\t%s%n", label,
                ReferenceQueryDump.escape(withoutComments(Files.readString(out,
                        StandardCharsets.UTF_8))));
        System.out.printf("segments\t%s\t%s%n", label, Files.exists(segments)
                ? ReferenceQueryDump.escape(withoutComments(Files.readString(segments,
                        StandardCharsets.UTF_8)))
                : "absent");
    }

    /** The comment lines carry the command line and the run's own clock. */
    static String withoutComments(final String text) {
        final List<String> kept = new ArrayList<>();
        for (final String line : text.split("\n", -1)) {
            if (!line.startsWith("#") || line.startsWith("#<METADATA>")) {
                if (!line.isEmpty()) {
                    kept.add(line);
                }
            }
        }
        return String.join("\n", kept);
    }
}
