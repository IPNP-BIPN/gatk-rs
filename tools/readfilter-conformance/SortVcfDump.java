/*
 * SortVcf, taken from the reference.
 *
 * One or more VCFs read into a sorting collection and written back out in dictionary order, under a
 * header merged from all of them.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE ORDER IS THE DICTIONARY'S, NOT THE ALPHABET'S. Records on `chr10` sort after `chr2` when
 *     the header declares them in that order, and a file whose contigs are declared in an unusual
 *     order sorts into that order rather than into a natural one;
 *   - THE SORT IS BY CONTIG AND POSITION ONLY, so two records at the same locus keep the order they
 *     were read in, which for several inputs is the order the inputs were given;
 *   - THE FIRST INPUT'S DICTIONARY BECOMES THE RUN'S, and a second input whose dictionary differs
 *     is a refusal rather than a merge;
 *   - THE SAMPLE CHECK IS ON THE SORTED NAMES, `getSampleNamesInOrder()`, so two files listing the
 *     same samples in different column orders are NOT refused: the output's columns are the sorted
 *     names and each file's genotypes are placed by name;
 *   - THE HEADER IS THE SMART MERGE of every input's, so a line only one file carried is in the
 *     output and two INFO lines with the same ID and different descriptions are a conflict;
 *   - A FILE WITH NO DICTIONARY AT ALL IS A REFUSAL unless SEQUENCE_DICTIONARY is given;
 *   - A RECORD ON A CONTIG THE DICTIONARY DOES NOT DECLARE stops the run with a NULL POINTER
 *     rather than a message: the comparator looks the contig up in a map and unboxes what it finds;
 *   - THE OUTPUT'S SAMPLE COLUMNS ARE THE SORTED NAMES, whatever order any input had;
 *   - AND THE RECORDS THEMSELVES GO BACK OUT THROUGH THE WRITER, so their QUAL and INFO spellings
 *     are the writer's rather than the inputs'.
 *
 * Output:
 *
 *     input\t<label>/<n>=<the whole input vcf, escaped>
 *     sorted\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: SortVcfDump
 */

import picard.vcf.SortVcf;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SortVcfDump {

    /** A header whose contigs are declared chr10 before chr2, which no alphabet would produce. */
    static String header(final String contigs, final String extraInfo, final String samples) {
        return "##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + extraInfo
                + contigs
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples + "\n";
    }

    static final String CONTIGS =
            "##contig=<ID=chr10,length=240>\n##contig=<ID=chr2,length=240>\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("sort-vcf-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# SortVcfDump: a VCF sorted into its dictionary's order");

        // Records given in the wrong order, on contigs whose declared order is not alphabetical.
        final String unsorted = header(CONTIGS, "", "NA1")
                + "chr2\t30\t.\tA\tC\t50.00\tPASS\tAC=1\tGT\t0/1\n"
                + "chr10\t20\t.\tA\tG\t.\tPASS\tAC=1\tGT\t0/1\n"
                + "chr2\t10\t.\tA\tT\t10.5\tPASS\tAC=1\tGT\t1/1\n"
                + "chr10\t5\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n";
        run("one-file", dir, new String[] {unsorted});

        // Two files whose records interleave, and two at the same locus.
        final String first = header(CONTIGS, "", "NA1")
                + "chr2\t10\tfirst\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"
                + "chr10\t20\tfirst\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n";
        final String second = header(CONTIGS, "", "NA1")
                + "chr2\t10\tsecond\tA\tG\t50\tPASS\tAC=1\tGT\t0/1\n"
                + "chr2\t20\tsecond\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n";
        run("two-files", dir, new String[] {first, second});
        // The same two the other way round, which is what says whether the tie is stable.
        run("two-files-reversed", dir, new String[] {second, first});

        // A line only the second file declares, which the merge keeps.
        final String withExtra = header(CONTIGS,
                "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">\n", "NA1")
                + "chr2\t15\t.\tA\tC\t50\tPASS\tDB\tGT\t0/1\n";
        run("merged-header", dir, new String[] {first, withExtra});

        // The same ID declared differently in the two files.
        final String conflicting = header(CONTIGS, "", "NA1")
                .replace("Description=\"Allele count\"", "Description=\"Something else\"")
                + "chr2\t15\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n";
        run("conflicting-header", dir, new String[] {first, conflicting});

        // Two samples, whose column order the output keeps.
        run("two-samples", dir, new String[] {
                header(CONTIGS, "", "zeta\talpha")
                        + "chr2\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t1/1\n"});
        // The same samples listed the other way round in the second file.
        run("sample-order-differs", dir, new String[] {
                header(CONTIGS, "", "zeta\talpha")
                        + "chr2\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t1/1\n",
                header(CONTIGS, "", "alpha\tzeta")
                        + "chr2\t20\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t1/1\n"});

        // A second file declaring a different dictionary.
        run("different-dictionaries", dir, new String[] {first,
                header("##contig=<ID=chr2,length=999>\n", "", "NA1")
                        + "chr2\t15\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"});

        // No dictionary at all.
        run("no-dictionary", dir, new String[] {
                header("", "", "NA1") + "chr2\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"});

        // A record on a contig the dictionary does not declare.
        run("undeclared-contig", dir, new String[] {header(CONTIGS, "", "NA1")
                + "chr2\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"
                + "chrX\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"});

        // A file with no records at all.
        run("no-records", dir, new String[] {header(CONTIGS, "", "NA1")});
    }

    static void run(final String label, final Path dir, final String[] inputs) throws Exception {
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < inputs.length; i++) {
            final Path in = dir.resolve(label + "-" + i + ".vcf");
            Files.writeString(in, inputs[i], StandardCharsets.UTF_8);
            System.out.printf("input\t%s/%d=%s%n", label, i,
                    ReferenceQueryDump.escape(inputs[i]));
            argv.add("I=" + in);
        }
        final Path out = dir.resolve("sorted-" + label + ".vcf");
        argv.addAll(Arrays.asList("O=" + out, "CREATE_INDEX=false"));
        try {
            final Object code = new SortVcf().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            // The dictionary refusal names the input's path, which is where the run happened.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())
                            .replaceAll("/work/sort-vcf-dump/[^ ]+", "<masked>")));
            return;
        }
        System.out.printf("sorted\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
