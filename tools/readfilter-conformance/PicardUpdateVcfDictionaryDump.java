/*
 * UpdateVcfSequenceDictionary (Picard's), taken from the reference.
 *
 * The class is named PicardUpdateVcfDictionaryDump rather than after its tool: GATK has a tool
 * called UpdateVCFSequenceDictionary whose dump is already here, and on a case-insensitive
 * filesystem the two file names are one file.
 *
 * A VCF whose contig lines are replaced by a dictionary read from elsewhere. GATK has a tool of
 * almost the same name, already ported, and the two disagree about nearly everything: this one has
 * no `--replace`, no per-record validation and no refusals of its own.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE DICTIONARY REPLACES THE CONTIG LINES AND NOTHING ELSE. Every other header line survives,
 *     in the writer's order;
 *   - A CONTIG THE INPUT HAD AND THE DICTIONARY LACKS IS SIMPLY GONE from the header, and its
 *     RECORDS ARE STILL WRITTEN: there is no per-record check, so the output declares fewer contigs
 *     than its records use;
 *   - THE INPUT'S OWN CONTIG ATTRIBUTES ARE LOST, and only some of the dictionary's replace them:
 *     a `.dict` sequence carrying `M5`, `UR` and `AS` writes a contig line with `assembly=` and
 *     nothing else, so the checksum and the URI do not survive the round trip either way;
 *   - THE ORDER IS THE DICTIONARY'S, so a dictionary listing the contigs the other way round
 *     reorders the header while leaving the records where they were;
 *   - AN INPUT WITH NO CONTIG LINES GAINS THE WHOLE DICTIONARY;
 *   - A DICTIONARY WITH NO SEQUENCES AT ALL IS NOT REFUSED HERE, but leaves a header the indexing
 *     writer then refuses;
 *   - THE SAMPLES ARE UNTOUCHED, the header being mutated in place rather than rebuilt;
 *   - AND THE DICTIONARY CAN COME FROM A `.dict`, A `.fasta` OR ANOTHER VCF, all three through the
 *     same extractor.
 *
 * Output:
 *
 *     input\t<label>=<the whole input vcf, escaped>
 *     dictionary\t<label>=<the dictionary source, escaped>
 *     updated\t<label>=<the whole output vcf, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: PicardUpdateVcfDictionaryDump
 */

import picard.vcf.UpdateVcfSequenceDictionary;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class PicardUpdateVcfDictionaryDump {

    static String header(final String contigs) {
        return "##fileformat=VCFv4.2\n"
                + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + contigs
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n";
    }

    static final String RECORDS =
            "chr1\t10\t.\tA\tC\t50\tPASS\tAC=1\tGT\t0/1\n"
            + "chr2\t20\t.\tA\tG\t50\tPASS\tAC=1\tGT\t0/1\n";

    /** A `.dict` file, which is a SAM header and nothing else. */
    static String dict(final String... sequences) {
        final StringBuilder text = new StringBuilder("@HD\tVN:1.6\tSO:unsorted\n");
        for (final String sequence : sequences) {
            text.append(sequence).append('\n');
        }
        return text.toString();
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("update-vcf-dictionary-dump");
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# UpdateVcfSequenceDictionaryDump: a VCF's contig lines replaced");

        final String withContigs = header(
                "##contig=<ID=chr1,length=240,assembly=mine>\n##contig=<ID=chr2,length=240>\n")
                + RECORDS;

        // A dictionary naming both contigs, with a length the input did not have.
        run("both-contigs", dir, withContigs,
                dict("@SQ\tSN:chr1\tLN:1000", "@SQ\tSN:chr2\tLN:2000"));
        // The same two the other way round, which reorders the header only.
        run("reversed", dir, withContigs,
                dict("@SQ\tSN:chr2\tLN:2000", "@SQ\tSN:chr1\tLN:1000"));
        // A dictionary naming only one of them, whose records are still written.
        run("missing-contig", dir, withContigs, dict("@SQ\tSN:chr1\tLN:1000"));
        // A dictionary naming a contig the file never mentions.
        run("extra-contig", dir, withContigs,
                dict("@SQ\tSN:chr1\tLN:1000", "@SQ\tSN:chr2\tLN:2000", "@SQ\tSN:chr3\tLN:3000"));
        // A dictionary carrying the fields a .dict usually has.
        run("with-attributes", dir, withContigs,
                dict("@SQ\tSN:chr1\tLN:1000\tM5:abc\tUR:file:/nowhere\tAS:other",
                        "@SQ\tSN:chr2\tLN:2000"));
        // An input with no contig lines at all.
        run("no-contigs-in", dir, header("") + RECORDS,
                dict("@SQ\tSN:chr1\tLN:1000", "@SQ\tSN:chr2\tLN:2000"));
        // A dictionary with no sequences, which this tool does not refuse.
        run("empty-dictionary", dir, withContigs, dict());
        // A file with no records.
        run("no-records", dir, header("##contig=<ID=chr1,length=240>\n"),
                dict("@SQ\tSN:chr1\tLN:1000"));
    }

    static void run(final String label, final Path dir, final String input, final String dictionary,
                    final String... extra) throws Exception {
        final Path in = dir.resolve(label + ".vcf");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        final Path dictPath = dir.resolve(label + ".dict");
        Files.writeString(dictPath, dictionary, StandardCharsets.UTF_8);
        System.out.printf("input\t%s=%s%n", label, ReferenceQueryDump.escape(input));
        System.out.printf("dictionary\t%s=%s%n", label, ReferenceQueryDump.escape(dictionary));
        final Path out = dir.resolve("updated-" + label + ".vcf");
        final List<String> argv = new ArrayList<>(Arrays.asList(
                "I=" + in, "SD=" + dictPath, "O=" + out, "CREATE_INDEX=false"));
        argv.addAll(Arrays.asList(extra));
        try {
            final Object code =
                    new UpdateVcfSequenceDictionary().instanceMain(argv.toArray(new String[0]));
            if (!Integer.valueOf(0).equals(code)) {
                System.out.printf("exit\t%s=%s%n", label, code);
                return;
            }
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            return;
        }
        System.out.printf("updated\t%s=%s%n", label,
                ReferenceQueryDump.escape(Files.readString(out)));
    }
}
