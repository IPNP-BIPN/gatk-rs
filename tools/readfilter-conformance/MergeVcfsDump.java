/*
 * MergeVcfs, taken from the reference.
 *
 * Several already-sorted VCFs merged into one by a merging iterator, under a header smart-merged
 * from all of them. The sibling of SortVcf, and the differences between the two are the point.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE MERGE IS A HEAP OVER SORTED INPUTS, not a sort of everything: a file whose own records
 *     are out of order is a REFUSAL from the merging iterator, naming the comparator class and no
 *     file at all;
 *   - A TIE BETWEEN TWO INPUTS IS NOT DECIDED BY THE ORDER THEY WERE GIVEN. Both orders of the
 *     same pair write the same two records in the same order, which is the opposite of SortVcf,
 *     where the input order decides;
 *   - THE CONTIG CHECK IS `isCompatible`, WHICH IS ABOUT INDICES: a second file must declare each
 *     shared contig at the SAME index, so both a subset that shifts an index and a reordering are
 *     refused, and the refusal names the file;
 *   - THE SAMPLE CHECK IS ON THE SORTED NAMES, so two files whose sample columns are ordered
 *     differently are accepted and the output's columns are the sorted names;
 *   - EVERY COMMENT IS WRITTEN AS A LINE KEYED `MergeVcfs.comment`, and the smart merge keys
 *     unstructured lines by their key, so TWO COMMENTS COLLAPSE INTO ONE: asking for two notes
 *     writes the first and silently drops the second;
 *   - THE HEADER IS A LinkedHashSet OF HEADERS, so two identical input headers collapse to one
 *     before the merge ever sees them;
 *   - A FILE WITH NO CONTIG LINES IS A REFUSAL unless a dictionary is supplied;
 *   - THE INDEX CHECK RUNS AFTER EVERY FILE IS READ, so its refusal names no file at all;
 *   - AND THE RECORDS GO BACK OUT THROUGH THE WRITER, so their QUAL spellings are the writer's.
 *
 * Output:
 *
 *     input\t<label>/<n>=<the whole input vcf, escaped>
 *     merged\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: MergeVcfsDump
 */

import picard.vcf.MergeVcfs;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MergeVcfsDump {

    static String header(final String contigs, final String samples) {
        return "##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + contigs
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples + "\n";
    }

    /** chr10 before chr2, so the order is visibly the dictionary's. */
    static final String CONTIGS =
            "##contig=<ID=chr10,length=240>\n##contig=<ID=chr2,length=240>\n";

    static String record(final String contig, final int position, final String id,
                         final String qual) {
        return contig + "\t" + position + "\t" + id + "\tA\tC\t" + qual + "\tPASS\tAC=1\tGT\t0/1\n";
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("merge-vcfs-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# MergeVcfsDump: several sorted VCFs merged into one");

        final String first = header(CONTIGS, "NA1")
                + record("chr10", 5, "first-a", "50.00")
                + record("chr2", 10, "first-b", "10.5")
                + record("chr2", 30, "first-c", "50");
        final String second = header(CONTIGS, "NA1")
                + record("chr10", 20, "second-a", "50")
                + record("chr2", 10, "second-b", "50")
                + record("chr2", 20, "second-c", "50");

        run("two-files", dir, new String[] {first, second});
        // The same two the other way round, which is what says whether the tie at chr2:10 is
        // decided by the input order. It is not.
        run("reversed", dir, new String[] {second, first});
        // One file alone.
        run("one-file", dir, new String[] {first});
        // Three files, the third holding one record between the others.
        run("three-files", dir, new String[] {first, second,
                header(CONTIGS, "NA1") + record("chr2", 15, "third-a", "50")});
        // Two identical files, which the LinkedHashSet of headers collapses.
        run("identical-files", dir, new String[] {first, first});
        // A file whose own records are out of order, which the merging iterator refuses.
        run("unsorted-input", dir, new String[] {
                header(CONTIGS, "NA1")
                        + record("chr2", 30, "out-a", "50")
                        + record("chr2", 10, "out-b", "50"),
                second});
        // A second file declaring only one of the first's contigs, which shifts that contig's
        // index and is therefore not compatible.
        run("subset-contigs", dir, new String[] {first,
                header("##contig=<ID=chr2,length=240>\n", "NA1") + record("chr2", 25, "sub", "50")});
        // A second file declaring the same contigs in the other order, which it does not.
        run("reordered-contigs", dir, new String[] {first,
                header("##contig=<ID=chr2,length=240>\n##contig=<ID=chr10,length=240>\n", "NA1")
                        + record("chr2", 25, "reordered", "50")});
        // Two samples listed in different orders.
        run("sample-order-differs", dir, new String[] {
                header(CONTIGS, "zeta\talpha")
                        + "chr2\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t1/1\n",
                header(CONTIGS, "alpha\tzeta")
                        + "chr2\t20\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\t1/1\n"});
        // A file with different samples, which is a refusal.
        run("different-samples", dir, new String[] {first,
                header(CONTIGS, "NA2") + record("chr2", 25, "other", "50")});
        // Two comments, of which only the first survives the merge.
        run("comments", dir, new String[] {first, second}, "CO=first note", "CO=second note");
        // No contig lines at all.
        run("no-contigs", dir, new String[] {
                header("", "NA1") + record("chr2", 10, "nodict", "50")});
        // A file with no records, merged with one that has them.
        run("empty-file", dir, new String[] {header(CONTIGS, "NA1"), second});
    }

    static void run(final String label, final Path dir, final String[] inputs,
                    final String... extra) throws Exception {
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < inputs.length; i++) {
            final Path in = dir.resolve(label + "-" + i + ".vcf");
            Files.writeString(in, inputs[i], StandardCharsets.UTF_8);
            System.out.printf("input\t%s/%d=%s%n", label, i, ReferenceQueryDump.escape(inputs[i]));
            argv.add("I=" + in);
        }
        final Path out = dir.resolve("merged-" + label + ".vcf");
        argv.addAll(Arrays.asList("O=" + out, "CREATE_INDEX=false"));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code = new MergeVcfs().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            // Two refusals name the input's path, which is where the run happened.
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())
                            .replaceAll("/work/merge-vcfs-dump/[^ ]+", "<masked>")));
            return;
        }
        System.out.printf("merged\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
