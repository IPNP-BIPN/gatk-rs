/*
 * RemoveNearbyIndels' output, taken from the reference.
 *
 * The first member of the variant-transform archetype: a VariantWalker that buffers one indel at a
 * time and drops any PAIR of indels closer than a given spacing, keeping whatever non-indels sat
 * between them.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE DISTANCE IS `right.getStart() - left.getEnd()`, END to START, so a deletion's length
 *     counts against it and two indels exactly `minIndelSpacing` apart are KEPT while one base
 *     closer is dropped;
 *   - AND lastIndel IS UPDATED EVEN WHEN THE INDEL WAS JUST THROWN AWAY, `lastIndel = vc.isIndel()
 *     ? vc : lastIndel`, so a third indel near the second is measured against an indel that never
 *     reached the output. A chain of three close indels loses all three;
 *   - THE FIRST BRANCH ONLY EVER FIRES ONCE, because lastIndel never returns to null, so the whole
 *     file after the first indel goes through `nearby` or the flush;
 *   - THROWING OUT A PAIR EMITS THE NON-INDELS BETWEEN THEM, which is what the doc means by
 *     "regardless of any intervening non-indel variants";
 *   - `nearby` REQUIRES THE SAME CONTIG, so a new contig flushes whatever was held;
 *   - THE LAST INDEL OF A FILE IS KEPT BY A REFERENCE COMPARISON, `vc == lastIndel` in
 *     emitRemaining, and not by an equality: without it the buffered indel would be measured
 *     against ITSELF, `start - end < spacing` being true for any single indel, and dropped;
 *   - A SPACING OF ZERO KEEPS EVERY INDEL THAT DOES NOT OVERLAP the one before it, since the test
 *     is strict;
 *   - THE HEADER IS REBUILT FROM THE INPUT'S SORTED METADATA plus the tool's own lines, so a
 *     ##GATKCommandLine line with the run's DATE is added, which is why it is masked here;
 *   - AND `isIndel` IS `getType() == INDEL`, WHICH A MIXED SITE IS NOT. A record whose alternates
 *     are all shorter than the reference is an INDEL and buffers like one; a record with a SNP and
 *     an insertion side by side is MIXED, so it is not buffered and cannot pair with anything, and
 *     a record with two SNPs is a plain SNP.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: RemoveNearbyIndelsDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.validation.RemoveNearbyIndels;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public class RemoveNearbyIndelsDump {

    static final String HEADER =
            "##fileformat=VCFv4.2\n"
                    + "##contig=<ID=chr1,length=250000000>\n"
                    + "##contig=<ID=chr2,length=240000000>\n"
                    + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("removenearbyindels-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# RemoveNearbyIndelsDump: RemoveNearbyIndels' output, from the reference");

        // Two indels ten apart, with a snp between them, and one isolated indel far away.
        final Path spaced = write(dir, "spaced",
                snp("chr1", 100),
                deletion("chr1", 200, 3),
                snp("chr1", 205),
                insertion("chr1", 210),
                snp("chr1", 300),
                insertion("chr1", 1000));

        // Three indels in a row, each close to the one before: the third is measured against an
        // indel that was already thrown away.
        final Path chain = write(dir, "chain",
                insertion("chr1", 100),
                insertion("chr1", 105),
                insertion("chr1", 110),
                snp("chr1", 500));

        // A file ending on an indel, which only emitRemaining can keep.
        final Path trailing = write(dir, "trailing",
                snp("chr1", 100),
                insertion("chr1", 500));

        // Two indels whose distance is measured from the END of the first: the deletion is long
        // enough that a following indel is nearby though its start is far.
        final Path longDeletion = write(dir, "long-deletion",
                deletion("chr1", 100, 40),
                insertion("chr1", 150));

        // Indels on either side of a contig boundary, which `nearby` refuses to pair.
        final Path contigs = write(dir, "contigs",
                insertion("chr1", 1000),
                insertion("chr2", 1),
                snp("chr2", 5));

        // Two multi-allelic sites. The first has two alleles that are BOTH shorter than the
        // reference, so its type is INDEL and it buffers like one. The second has a snp and an
        // insertion, so its type is MIXED and `isIndel` is FALSE: it is not buffered at all.
        final Path mixed = write(dir, "mixed",
                "chr1\t100\t.\tACG\tA,T\t.\t.\t.",
                snp("chr1", 104),
                insertion("chr1", 108));
        final Path multiAllelic = write(dir, "multi-allelic",
                "chr1\t100\t.\tA\tC,AGG\t.\t.\t.",
                insertion("chr1", 104),
                "chr1\t500\t.\tA\tC,G\t.\t.\t.",
                insertion("chr1", 504));

        // Two indels at exactly the spacing under test, and two one base closer.
        final Path boundary = write(dir, "boundary",
                deletion("chr1", 100, 2),
                insertion("chr1", 106),
                deletion("chr1", 200, 2),
                insertion("chr1", 205));

        for (final Path input : List.of(spaced, chain, trailing, longDeletion, contigs, mixed,
                multiAllelic)) {
            run(dir, input, 20);
        }
        // The boundary file at the two spacings that straddle it, and at zero.
        run(dir, boundary, 5);
        run(dir, boundary, 4);
        run(dir, boundary, 0);
        // And the spaced file with a spacing wide enough to swallow everything.
        run(dir, spaced, 1000);
    }

    static String snp(final String contig, final int position) {
        return contig + "\t" + position + "\t.\tA\tC\t.\t.\t.";
    }

    /** A deletion whose reference allele is `length` bases long, so its end is start + length - 1. */
    static String deletion(final String contig, final int position, final int length) {
        final StringBuilder reference = new StringBuilder("A");
        for (int i = 1; i < length; i++) {
            reference.append("C");
        }
        return contig + "\t" + position + "\t.\t" + reference + "\tA\t.\t.\t.";
    }

    static String insertion(final String contig, final int position) {
        return contig + "\t" + position + "\t.\tA\tAGG\t.\t.\t.";
    }

    /** A vcf written by hand and indexed, which is what a VariantWalker is given. */
    static Path write(final Path dir, final String label, final String... records) throws Exception {
        final Path file = dir.resolve(label + ".vcf");
        final StringBuilder text = new StringBuilder(HEADER);
        for (final String record : records) {
            text.append(record).append("\n");
        }
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label,
                ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final Path input, final int spacing) {
        final String name = input.getFileName().toString().replace(".vcf", "");
        final String label = name + "-at-" + spacing;
        final Path output = dir.resolve(label + ".vcf");
        try {
            new RemoveNearbyIndels().instanceMain(new String[] {
                    "-V", input.toString(),
                    "-O", output.toString(),
                    "--min-indel-spacing", String.valueOf(spacing),
            });
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        print(label, output);
    }

    /** Every line of the output, with the one line that carries the run's date masked. */
    static void print(final String label, final Path output) {
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
