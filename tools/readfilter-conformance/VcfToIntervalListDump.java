/*
 * VcfToIntervalList, taken from the reference.
 *
 * A VCF read as a stream of intervals and written back out as a Picard interval list, merged on the
 * way through. The merging is `IntervalList.IntervalMergerIterator`, which no other ported tool
 * reaches, and the naming is the part that carries the surprises.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE MERGER IS A STREAM AND NEVER SORTS. It compares each interval with the one before it, so
 *     an input whose records are out of order comes out out of order, and two overlapping records
 *     separated by a third do not merge;
 *   - ABUTTING INTERVALS MERGE, `combineAbuttingIntervals` being true, so 50-50 and 51-51 are one
 *     interval while 60-60 and 62-62 are two;
 *   - AN UNNAMED RECORD IS `interval-<n>` AND THE COUNTER ONLY COUNTS UNNAMED ONES, so it is not
 *     the record's position in the file;
 *   - AND THE COUNTER RUNS AFTER THE FILTERING, the stream filtering before it maps, so dropping a
 *     filtered record renumbers every unnamed interval after it;
 *   - THE MERGED NAME IS A LinkedHashSet JOINED WITH A PIPE, so a name shared by two merged records
 *     appears once, in the order first met;
 *   - `INCLUDE_FILTERED` CHANGES WHICH INTERVALS MERGE, not merely which are present: a filtered
 *     record between two others is the bridge that joins them;
 *   - A FILTER OF `PASS` IS NOT A FILTER, so a record carrying it survives the default run;
 *   - THE END IS `getAttributeAsInt(END, vc.getEnd())`, so a symbolic ALT with END=200 is a 101
 *     base interval and a deletion's interval is as long as its reference allele;
 *   - EVERY INTERVAL IS BUILT WITH `negative=false`, so the strand column is always `+`;
 *   - AND THE HEADER IS `new SAMFileHeader(vcfHeader.getSequenceDictionary())`, WHICH IS NULL WHEN
 *     THE FILE DECLARES NO CONTIGS: the codec then fails on the null dictionary rather than the
 *     tool refusing anything, and the output file is left behind.
 *
 * VARIANT_ID_METHOD is a STATIC field, so a run that sets it leaves it set for every later run in
 * the same JVM. The USE_FIRST run is therefore last here, and nothing may be added after it.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     list\t<label>=<the whole output interval list, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: VcfToIntervalListDump
 */

import picard.vcf.VcfToIntervalList;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class VcfToIntervalListDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
            + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
            + "##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
            + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
            + "##contig=<ID=chr1,length=1000>\n"
            + "##contig=<ID=chr2,length=1000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    /** One record per behaviour, in position order. */
    static final String RECORDS =
            // A named record on its own.
            "chr1\t10\trs1\tA\tC\t.\t.\t.\n"
            // An unnamed one, which is the first interval-<n>.
            + "chr1\t20\t.\tA\tC\t.\t.\t.\n"
            // A deletion whose interval is as long as its reference allele, and a named record
            // inside it: they overlap and merge, and the merged name has both parts.
            + "chr1\t30\t.\tACGT\tA\t.\t.\t.\n"
            + "chr1\t33\trs2\tA\tC\t.\t.\t.\n"
            // Abutting, which merges.
            + "chr1\t50\trs3\tA\tC\t.\t.\t.\n"
            + "chr1\t51\trs4\tA\tC\t.\t.\t.\n"
            // One base apart, which does not.
            + "chr1\t60\trs5\tA\tC\t.\t.\t.\n"
            + "chr1\t62\trs6\tA\tC\t.\t.\t.\n"
            // An UNNAMED filtered record, and an unnamed one after it: the counter runs after
            // the filtering, so including it renumbers every unnamed interval that follows.
            + "chr1\t70\t.\tA\tC\t.\tq10\t.\n"
            + "chr1\t71\t.\tA\tC\t.\t.\t.\n"
            // A filtered record BETWEEN two others, which bridges them when it is included.
            + "chr1\t80\trs7\tA\tC\t.\t.\t.\n"
            + "chr1\t81\tbridge\tA\tC\t.\tq10\t.\n"
            + "chr1\t82\trs8\tA\tC\t.\t.\t.\n"
            // A record whose FILTER is PASS, which is not filtered.
            + "chr1\t90\tpassing\tA\tC\t.\tPASS\t.\n"
            // A symbolic ALT whose END attribute is the interval's end.
            + "chr1\t100\tsymbolic\tA\t<DEL>\t.\t.\tEND=200\n"
            // Two records sharing one name, merged into one interval with one name.
            + "chr1\t300\tdup\tA\tC\t.\t.\t.\n"
            + "chr1\t301\tdup\tA\tC\t.\t.\t.\n"
            // An unnamed record merged with a named one, and a record whose ID holds two names.
            + "chr1\t400\t.\tA\tC\t.\t.\t.\n"
            + "chr1\t401\trs9;rs10\tA\tC\t.\t.\t.\n"
            // The contig changes, which never merges however close the positions are.
            + "chr2\t10\trs11\tA\tC\t.\t.\t.\n";

    /** The same records with three of them out of order, which the merger does not sort. */
    static final String UNSORTED =
            "chr1\t500\tlate\tA\tC\t.\t.\t.\n"
            + "chr1\t480\tearly\tA\tC\t.\t.\t.\n"
            + "chr1\t501\tafter-late\tA\tC\t.\t.\t.\n"
            + "chr1\t600\tlast\tA\tC\t.\t.\t.\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("vcf-to-interval-list-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# VcfToIntervalListDump: a VCF turned into a Picard interval list");

        final String input = HEADER + RECORDS;
        run("defaults", dir, input);
        // The two filtered records kept, which changes both the numbering and the merging.
        run("include-filtered", dir, input, "INCLUDE_FILTERED=true");
        // Records out of order, which come out out of order.
        run("unsorted", dir, HEADER + UNSORTED);
        // A file with no records, which is a header and nothing else.
        run("no-records", dir, HEADER);
        // A header with no contig lines, whose dictionary is null.
        run("no-contigs", dir,
                HEADER.replace("##contig=<ID=chr1,length=1000>\n", "")
                        .replace("##contig=<ID=chr2,length=1000>\n", "")
                + RECORDS);
        // VARIANT_ID_METHOD is static and this run leaves it set: it must stay last.
        run("use-first", dir, input, "VARIANT_ID_METHOD=USE_FIRST");
    }

    static void run(final String label, final Path dir, final String input, final String... extra)
            throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        final Path out = dir.resolve(label + ".interval_list");
        final List<String> argv = new ArrayList<>(Arrays.asList("I=" + in, "O=" + out));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new VcfToIntervalList().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("list\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
