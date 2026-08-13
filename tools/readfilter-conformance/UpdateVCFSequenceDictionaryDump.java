/*
 * UpdateVCFSequenceDictionary's output, taken from the reference.
 *
 * The second member of the variant-transform archetype, and the calibration point: the same walker
 * shape as RemoveNearbyIndels, doing something else entirely. It replaces the contig lines of a vcf
 * with a dictionary taken from somewhere else, and refuses in five different ways.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE HEADER IS REBUILT FROM `getMetaDataInInputOrder()`, where RemoveNearbyIndels uses
 *     `getMetaDataInSortedOrder()`, AND THE DIFFERENCE NEVER REACHES THE FILE: the writer emits the
 *     lines in its own order, so an input carrying INFO before ALT comes out ALT before INFO either
 *     way. Two tools of one archetype differ on that line and produce the same order;
 *   - AND setSequenceDictionary REPLACES THE CONTIG LINES rather than adding to them, so a contig
 *     the input had and the dictionary does not is gone from the header while its records stay;
 *   - WITHOUT --replace, AN INPUT THAT ALREADY HAS A DICTIONARY IS REFUSED, and the message names
 *     the FEATURE INPUT's name, which is the argument's value and not the file name;
 *   - THE CHECK IS ON THE HEADER READ FROM THE FILE, not on the engine's best dictionary, which
 *     the comment says is because the engine might dig one up from an index;
 *   - A VARIANT ON A CONTIG THE DICTIONARY DOES NOT HAVE IS REFUSED, and the message quotes the
 *     variant's ID, which is "." when the record has none;
 *   - A VARIANT THAT ENDS PAST THE SEQUENCE LENGTH IS REFUSED TOO, and END IS `vc.getEnd()`, which
 *     the INFO field END overrides: a one-base record with END=250000001 is refused though its
 *     reference allele is one base long;
 *   - THE VALIDATION IS PER RECORD AND THE OUTPUT IS ALREADY OPEN, so the refusal comes after a
 *     header and any earlier records have been written;
 *   - AN EMPTY DICTIONARY SOURCE IS A BadArgumentValue naming the source;
 *   - A DICTIONARY WITH A SEQUENCE OF UNKNOWN LENGTH IS REFUSED by a different exception, and
 *     UNKNOWN_SEQUENCE_LENGTH IS 0, so `LN:0` is what triggers it;
 *   - AND GIVING BOTH --source-dictionary AND --sequence-dictionary IS A PLAIN CommandLineException,
 *     while giving NEITHER is a MissingArgument naming only the first of the two.
 *
 * Output:
 *
 *     input\t<label>\t<the whole input vcf, escaped>
 *     dict\t<label>\t<the whole dictionary file, escaped>
 *     vcfline\t<label>\t<one line of the output VCF, escaped>
 *     commandline\t<label>\t<the ##GATKCommandLine line with its date masked>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: UpdateVCFSequenceDictionaryDump
 */

import org.broadinstitute.hellbender.tools.IndexFeatureFile;
import org.broadinstitute.hellbender.tools.walkers.variantutils.UpdateVCFSequenceDictionary;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class UpdateVCFSequenceDictionaryDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("updatevcfdictionary-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# UpdateVCFSequenceDictionaryDump: the dictionary replacement, from the reference");

        // A vcf with no contig lines at all, and one that has them, with the same records.
        final Path bare = writeVcf(dir, "bare", false,
                "chr1\t100\t.\tA\tC\t.\t.\t.",
                "chr2\t200\trs1\tA\tC\t.\t.\t.");
        final Path withContigs = writeVcf(dir, "with-contigs", true,
                "chr1\t100\t.\tA\tC\t.\t.\t.",
                "chr2\t200\trs1\tA\tC\t.\t.\t.");
        // One whose second record is on a contig no dictionary here has.
        final Path unknownContig = writeVcf(dir, "unknown-contig", false,
                "chr1\t100\t.\tA\tC\t.\t.\t.",
                "chrUn\t1\tbad\tA\tC\t.\t.\t.");
        // One whose record runs past the end of its contig, by its reference allele.
        final Path pastEnd = writeVcf(dir, "past-end", false,
                "chr2\t239999999\tlong\tACCCC\tA\t.\t.\t.");
        // And one whose END attribute is what runs past it, though the record is one base.
        final Path endAttribute = writeVcf(dir, "end-attribute", false,
                "chr1\t100\tsymbolic\tA\t<DEL>\t.\t.\tEND=250000001");

        // The dictionaries: one good, one empty, one whose second sequence has no length.
        final Path dictionary = writeDictionary(dir, "dictionary",
                "@SQ\tSN:chr1\tLN:250000000",
                "@SQ\tSN:chr2\tLN:240000000");
        final Path empty = writeDictionary(dir, "empty-dictionary");
        final Path noLength = writeDictionary(dir, "no-length",
                "@SQ\tSN:chr1\tLN:250000000",
                "@SQ\tSN:chr2\tLN:0");
        // A dictionary holding only one of the two contigs the input uses.
        final Path partial = writeDictionary(dir, "partial-dictionary",
                "@SQ\tSN:chr1\tLN:250000000");

        run(dir, "bare", bare, "--source-dictionary", dictionary.toString());
        run(dir, "with-contigs-refused", withContigs, "--source-dictionary", dictionary.toString());
        run(dir, "with-contigs-replaced", withContigs,
                "--source-dictionary", dictionary.toString(), "--replace", "true");
        // A dictionary missing a contig the input's HEADER has, which is dropped from the output
        // while its records are refused when they are reached.
        run(dir, "partial", bare, "--source-dictionary", partial.toString());
        run(dir, "unknown-contig", unknownContig, "--source-dictionary", dictionary.toString());
        run(dir, "past-end", pastEnd, "--source-dictionary", dictionary.toString());
        run(dir, "end-attribute", endAttribute, "--source-dictionary", dictionary.toString());
        run(dir, "empty-dictionary", bare, "--source-dictionary", empty.toString());
        run(dir, "no-length", bare, "--source-dictionary", noLength.toString());
        run(dir, "both-dictionaries", bare,
                "--source-dictionary", dictionary.toString(),
                "--sequence-dictionary", dictionary.toString());
        run(dir, "no-dictionary", bare);
    }

    /** A vcf written by hand and indexed, with or without contig lines of its own. */
    static Path writeVcf(final Path dir, final String label, final boolean withContigs,
                         final String... records) throws Exception {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n");
        text.append("##ALT=<ID=DEL,Description=\"Deletion\">\n");
        if (withContigs) {
            text.append("##contig=<ID=chr1,length=250000000>\n");
            text.append("##contig=<ID=chr2,length=240000000>\n");
        }
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final String record : records) {
            text.append(record).append("\n");
        }
        final Path file = dir.resolve(label + ".vcf");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        new IndexFeatureFile().instanceMain(new String[] {"-I", file.toString()});
        System.out.printf("input\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    /** A .dict file, which is a SAM header and nothing else. */
    static Path writeDictionary(final Path dir, final String label, final String... sequences)
            throws Exception {
        final StringBuilder text = new StringBuilder("@HD\tVN:1.6\n");
        for (final String sequence : sequences) {
            text.append(sequence).append("\n");
        }
        final Path file = dir.resolve(label + ".dict");
        Files.writeString(file, text.toString(), StandardCharsets.UTF_8);
        System.out.printf("dict\t%s\t%s%n", label, ReferenceQueryDump.escape(text.toString()));
        return file;
    }

    static void run(final Path dir, final String label, final Path input, final String... arguments) {
        final Path output = dir.resolve(label + "-out.vcf");
        final List<String> all = new ArrayList<>(List.of("-V", input.toString(), "-O", output.toString()));
        all.addAll(List.of(arguments));
        try {
            new UpdateVCFSequenceDictionary().instanceMain(all.toArray(new String[0]));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(String.valueOf(e.getMessage())));
            // The output may exist anyway: the writer is opened before the records are checked.
            print(label + "-partial", output);
            return;
        }
        print(label, output);
    }

    /** Every line of the output, with the one line that carries the run's date masked. */
    static void print(final String label, final Path output) {
        if (!Files.exists(output)) {
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
