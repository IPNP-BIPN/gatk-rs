/*
 * GatherVcfsCloud's output and its refusals, taken from the reference.
 *
 * Concatenating vcfs that are already in order, which sounds like `cat` and is not: the tool has two
 * gathering modes with different code paths, and it validates before it writes.
 *
 * Nine behaviours this is built to catch, three of which are worse than they look.
 *
 *   - AUTOMATIC PICKS BLOCK ONLY WHEN EVERY INPUT AND THE OUTPUT ARE BLOCK COMPRESSED, so the same
 *     files gathered to a `.vcf` and to a `.vcf.gz` take two different code paths, and asking for
 *     BLOCK on a plain vcf is a BadInput rather than a fallback;
 *   - THERE ARE TWO ORDER CHECKS, NOT ONE, and they throw different classes. The validation before
 *     writing compares the FIRST record of each file, and refuses with an
 *     IllegalArgumentException; the conventional writer then compares the LAST record written
 *     against the next file's first, and refuses with an IllegalStateException naming both
 *     positions. Files that overlap in the middle pass the first and are caught by the second;
 *   - --disable-contig-ordering-check DOES NOT DISABLE THE CHECK, IT WEAKENS IT, and the result is
 *     an INVALID FILE: it compares positions only within a contig, so gathering a chr2 shard before
 *     a chr1 shard is accepted and the output holds chr2:100 followed by chr1:100. The tool writes
 *     a vcf whose records are not in dictionary order;
 *   - --ignore-safety-checks SKIPS THE SAMPLE-LIST COMPARISON AND WRITES THE RESULT ANYWAY, so a
 *     record belonging to sample s1 is written under a header that declares only s0. The genotype
 *     is not dropped, it is RELABELLED;
 *   - A MISSING SEQUENCE DICTIONARY IS REFUSED BY THE INDEXER, not by the validation: the message
 *     is "In order to index the resulting VCF, the input VCFs must contain ##contig lines", a plain
 *     UserException, and it arrives after the validation has already passed;
 *   - DIFFERING SAMPLE LISTS ARE AN IllegalArgumentException whose message lists the difference in
 *     BOTH directions, each as a sorted set, and names the file by URI;
 *   - THE ORDER COMPARISON IS THE DICTIONARY'S, through a VariantContextComparator, so contig order
 *     is the header's and not alphabetical;
 *   - GATHERING ONE FILE IS STILL A GATHER, with the same header rewriting as for many;
 *   - AND THE HEADER OF THE OUTPUT IS THE FIRST FILE'S, whichever mode ran.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: GatherVcfsDump
 */

import org.broadinstitute.hellbender.tools.GatherVcfsCloud;
import org.broadinstitute.hellbender.tools.IndexFeatureFile;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class GatherVcfsDump {

    static String header(final String samples, final boolean withDictionary) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        if (withDictionary) {
            text.append("##contig=<ID=chr1,length=100000>\n");
            text.append("##contig=<ID=chr2,length=90000>\n");
        }
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t").append(samples).append("\n");
        return text.toString();
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("gathervcfs-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# GatherVcfsDump: gathering vcfs, from the reference");

        // Three shards in order, one sample.
        final Path first = write(dir, "first", header("s0", true),
                "chr1\t100\t.\tA\tC\t50\t.\t.\tGT\t0/1",
                "chr1\t200\t.\tA\tG\t50\t.\t.\tGT\t0/1");
        final Path second = write(dir, "second", header("s0", true),
                "chr1\t300\t.\tA\tT\t50\t.\t.\tGT\t1/1");
        final Path third = write(dir, "third", header("s0", true),
                "chr2\t100\t.\tA\tC\t50\t.\t.\tGT\t0/1");

        // The same shards out of order, a second sample list, one without a dictionary, and one
        // whose records overlap the previous file in the middle but not at its first record.
        final Path otherSample = write(dir, "other-sample", header("s1", true),
                "chr1\t400\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        final Path noDictionary = write(dir, "no-dictionary", header("s0", false),
                "chr1\t100\t.\tA\tC\t50\t.\t.\tGT\t0/1");
        // Sorted within itself, but starting inside the first file's range: the check compares
        // first records only, so this is accepted while the files genuinely overlap.
        final Path overlapping = write(dir, "overlapping", header("s0", true),
                "chr1\t150\t.\tA\tG\t50\t.\t.\tGT\t0/1",
                "chr1\t250\t.\tA\tC\t50\t.\t.\tGT\t0/1");

        // The gathers that work, to a plain vcf and to a bgzipped one.
        run(dir, "three-shards", "out.vcf", List.of(first, second, third));
        run(dir, "three-shards-gz", "out.vcf.gz", List.of(first, second, third));
        run(dir, "one-shard", "one.vcf", List.of(first));
        // The same three with the mode forced either way.
        run(dir, "conventional", "conventional.vcf", List.of(first, second, third),
                "--gather-type", "CONVENTIONAL");
        run(dir, "block-refused", "block.vcf", List.of(first, second, third),
                "--gather-type", "BLOCK");

        // The refusals.
        run(dir, "out-of-order", "order.vcf", List.of(second, first));
        run(dir, "out-of-order-check-disabled", "order2.vcf", List.of(second, first),
                "--disable-contig-ordering-check", "true");
        run(dir, "contig-out-of-order", "contig.vcf", List.of(third, first));
        run(dir, "contig-out-of-order-check-disabled", "contig2.vcf", List.of(third, first),
                "--disable-contig-ordering-check", "true");
        run(dir, "different-samples", "samples.vcf", List.of(first, otherSample));
        run(dir, "different-samples-ignored", "samples2.vcf", List.of(first, otherSample),
                "--ignore-safety-checks", "true");
        run(dir, "no-dictionary", "nodict.vcf", List.of(noDictionary, second));
        // Overlapping records whose FIRST records are still in order.
        run(dir, "overlapping-records", "overlap.vcf", List.of(first, overlapping));
    }

    static Path write(final Path dir, final String label, final String header,
                      final String... records) throws Exception {
        final StringBuilder text = new StringBuilder(header);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final String outputName,
                    final List<Path> inputs, final String... arguments) {
        final Path output = dir.resolve(outputName);
        final List<String> all = new ArrayList<>();
        for (final Path input : inputs) {
            all.add("-I");
            all.add(input.toString());
        }
        all.add("-O");
        all.add(output.toString());
        all.addAll(List.of(arguments));
        try {
            new GatherVcfsCloud().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, output);
    }

    static void print(final String label, final Path output) {
        if (output.toString().endsWith(".gz")) {
            // A bgzipped output is bytes rather than lines; its size and its first block say which
            // path was taken without reproducing a compressor here.
            try {
                System.out.printf("gzsize\t%s\t%d%n", label, Files.size(output));
            } catch (final Exception e) {
                System.out.printf("error\t%s-size\t%s:%s%n", label, e.getClass().getName(),
                        String.valueOf(e.getMessage()));
            }
            return;
        }
        final List<String> lines;
        try {
            lines = Files.readAllLines(output, StandardCharsets.UTF_8);
        } catch (final Exception e) {
            System.out.printf("error\t%s-read\t%s:%s%n", label, e.getClass().getName(),
                    String.valueOf(e.getMessage()));
            return;
        }
        for (final String line : lines) {
            if (line.startsWith("##GATKCommandLine")) {
                System.out.printf("commandline\t%s\t%s%n", label,
                        ReferenceQueryDump.escape(line.replaceAll("Date=\"[^\"]*\"", "Date=\"MASKED\"")));
                continue;
            }
            System.out.printf("vcfline\t%s\t%s%n", label, ReferenceQueryDump.escape(line));
        }
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
