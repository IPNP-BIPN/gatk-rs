/*
 * FuncotateSegments' annotated segments, taken from the reference.
 *
 * How a copy-number segment file is annotated from a folder of data sources. The folder's SHAPE is
 * as much of the tool as the annotation is: a manifest whose version must fall in a range, a
 * directory per source, a directory per reference version inside it, and a properties file that
 * says how to read the data.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE DATA SOURCE FOLDER IS THREE LEVELS DEEP: root, source name, reference version;
 *   - THE MANIFEST'S VERSION MUST FALL IN A RANGE, and a folder outside it is REFUSED, not
 *     skipped: a typo in the version stops the run;
 *   - A MISSING MANIFEST IS NOT REFUSED, though: the version check simply never happens and the
 *     sources load, so the file that guards the range is the one file that is optional;
 *   - A REFERENCE VERSION WITH NO DIRECTORY UNDER ANY SOURCE IS ALSO REFUSED, with a different
 *     message, and the same folder answers normally for the version it does have;
 *   - THE CONFIG FILE'S REQUIRED KEYS DEPEND ON ITS `type`, and a locatable XSV without its
 *     `end_column` is refused by name;
 *   - A GENCODE SOURCE IS REQUIRED, whatever else the folder holds: the renderers look for
 *     `Gencode_<version>_genes`, and the prefix is the source's own `name` key, capital G
 *     included, so a source named `gencode` in lower case annotates nothing findable;
 *   - A SEGMENT UNDER 150 BASES IS REFUSED rather than written through unannotated, and the
 *     complaint is about the variant context and not about the segment file;
 *   - THE SEG OUTPUT'S COLUMNS ARE A FIXED SET plus one group per source, so the locatable XSV
 *     source beside GENCODE contributes nothing visible to it;
 *   - A SEGMENT OVERLAPPING NO GENE IS STILL WRITTEN, with its Gencode columns empty, while the
 *     absent columns of a segment file are `__UNKNOWN__` and not empty;
 *   - AND EVERY RUN WRITES A SECOND FILE, `<output>.gene_list.txt`, one row per gene and exon,
 *     which holds ONLY the genes a segment covers: the segments covering nothing are absent from
 *     it, and a gene gets one row for the whole-segment `genes` field and one per exon.
 *
 * Output:
 *
 *     tsv\t<label>=<that file, escaped>
 *     out\t<label>=<the whole output, escaped>
 *     genes\t<label>=<the gene list beside it, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FuncotateSegmentsDump
 */

import org.broadinstitute.hellbender.tools.funcotator.FuncotateSegments;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class FuncotateSegmentsDump {

    static final int CONTIG_LENGTH = 199980;

    /** A data-source root: a manifest and one locatable XSV source under the given version. */
    static Path dataSources(final Path dir, final String name, final String manifestVersion,
                            final String refVersion, final boolean writeManifest) throws Exception {
        final Path root = dir.resolve(name);
        Files.createDirectories(root);
        if (writeManifest) {
            Files.writeString(root.resolve("MANIFEST.txt"),
                    "Version: " + manifestVersion + "\nSource: test\nAlternate Source: test\n",
                    StandardCharsets.UTF_8);
        }
        final Path source = root.resolve("regions").resolve(refVersion);
        Files.createDirectories(source);
        Files.writeString(source.resolve("regions.tsv"), String.join("\n",
                "CONTIG\tSTART\tEND\tLABEL\tSCORE",
                "chr1\t1000\t1999\tfirst\t10",
                "chr1\t1500\t2499\tsecond\t20",
                "chr1\t9000\t9999\tfar\t30",
                ""), StandardCharsets.UTF_8);
        Files.writeString(source.resolve("regions.config"), String.join("\n",
                "name = regions",
                "version = 1",
                "src_file = regions.tsv",
                "origin_location = test",
                "preprocessing_script =",
                "type = locatableXSV",
                // The three are column INDICES and not names, whatever the header row says.
                "contig_column = 0",
                "start_column = 1",
                "end_column = 2",
                "xsv_delimiter = \\t",
                ""), StandardCharsets.UTF_8);
        writeGencode(root, refVersion);
        return root;
    }

    /**
     * A Gencode source, which the tool REQUIRES: without one it refuses before any annotation.
     *
     * The GTF is the same shape as SVAnnotate's fixture, and the transcript FASTA beside it is
     * what `gencode_fasta_path` names.
     */
    static void writeGencode(final Path root, final String refVersion) throws Exception {
        final Path source = root.resolve("gencode").resolve(refVersion);
        Files.createDirectories(source);

        final String attributes = "gene_id \"ENSG00000000001.1\"; transcript_id "
                + "\"ENST00000000001.1\"; gene_type \"protein_coding\"; gene_name \"ALPHA\"; "
                + "transcript_type \"protein_coding\"; transcript_name \"ALPHA-201\"; "
                + "tag \"basic\"; transcript_status \"KNOWN\"; level 1;";
        final List<String> gtf = new ArrayList<>();
        // Exactly five header lines, in this order, or the codec refuses the file.
        gtf.add("##description: evidence-based annotation of the human genome, version 43 (test)");
        gtf.add("##provider: GENCODE");
        gtf.add("##contact: gencode@test");
        gtf.add("##format: gtf");
        gtf.add("##date: 2022-01-01");
        gtf.add("chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id \"ENSG00000000001.1\"; "
                + "gene_type \"protein_coding\"; gene_name \"ALPHA\"; level 1;");
        gtf.add("chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\t" + attributes);
        gtf.add("chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\t" + attributes + " exon_number 1; "
                + "exon_id \"ENSE00000000001.1\";");
        gtf.add("chr1\tHAVANA\tCDS\t1000\t1200\t.\t+\t0\t" + attributes + " exon_number 1; "
                + "exon_id \"ENSE00000000001.1\";");
        gtf.add("chr1\tHAVANA\tstart_codon\t1000\t1002\t.\t+\t0\t" + attributes
                + " exon_number 1; exon_id \"ENSE00000000001.1\";");
        gtf.add("chr1\tHAVANA\tstop_codon\t1198\t1200\t.\t+\t0\t" + attributes
                + " exon_number 1; exon_id \"ENSE00000000001.1\";");
        gtf.add("");
        final Path gtfPath = source.resolve("gencode.gtf");
        Files.writeString(gtfPath, String.join("\n", gtf), StandardCharsets.UTF_8);
        // The GTF is QUERIED by interval, so it needs an index beside it too.
        htsjdk.tribble.index.IndexFactory.createLinearIndex(gtfPath.toFile(),
                new org.broadinstitute.hellbender.utils.codecs.gtf.GencodeGtfCodec())
                .writeBasedOnFeatureFile(gtfPath.toFile());

        final StringBuilder transcripts = new StringBuilder(">ENST00000000001.1\n");
        for (int i = 0; i < 4; i++) {
            transcripts.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        final Path transcriptFasta = source.resolve("gencode.transcripts.fasta");
        Files.writeString(transcriptFasta, transcripts.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(transcriptFasta, true);
        // The transcript FASTA needs a dictionary beside it as well as an index.
        final htsjdk.samtools.SAMFileHeader transcriptHeader = new htsjdk.samtools.SAMFileHeader();
        transcriptHeader.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(List.of(
                new htsjdk.samtools.SAMSequenceRecord("ENST00000000001.1", 240))));
        try (final java.io.Writer writer = Files.newBufferedWriter(
                source.resolve("gencode.transcripts.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, transcriptHeader);
        }

        Files.writeString(source.resolve("gencode.config"), String.join("\n",
                // The field prefix is the datasource NAME, and the gene-list renderer looks for
                // `Gencode_<version>_genes` with a capital G: lowercase here finds nothing.
                "name = Gencode",
                "version = 1",
                "src_file = gencode.gtf",
                "origin_location = test",
                "preprocessing_script =",
                "type = gencode",
                "gencode_fasta_path = gencode.transcripts.fasta",
                "ncbi_build_version = " + refVersion,
                ""), StandardCharsets.UTF_8);
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("funcotate-segments-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FuncotateSegmentsDump: how a segment file is annotated from a folder "
                + "of data sources");

        final Path fasta = writeReference(dir);

        // The segments: one overlapping both source records, one overlapping neither, one
        // overlapping the far record alone.
        final String segments = String.join("\n",
                "CONTIG\tSTART\tEND\tCALL",
                "chr1\t1200\t1800\t+",
                "chr1\t5000\t5500\t0",
                // Over 150 bases, which is the minimum for a segment: see the short-segment run.
                "chr1\t9400\t9900\t-",
                "");
        final Path segmentsPath = write(dir, "segments.seg", segments);
        System.out.printf("tsv\tsegments=%s%n", ReferenceQueryDump.escape(segments));

        final Path good = dataSources(dir, "ds-good", "1.7.hg38.20220101", "hg38", true);
        final Path oldVersion = dataSources(dir, "ds-old", "1.2.hg38.20150101", "hg38", true);
        final Path noManifest = dataSources(dir, "ds-none", "1.7.hg38.20220101", "hg38", false);
        final Path hg19Only = dataSources(dir, "ds-hg19", "1.7.hg38.20220101", "hg19", true);

        run(dir, "annotated", segmentsPath, fasta, good, "hg38", List.of());
        // A version outside the accepted range, which is skipped in silence.
        run(dir, "old-version", segmentsPath, fasta, oldVersion, "hg38", List.of());
        // No manifest at all, which is also a skip.
        run(dir, "no-manifest", segmentsPath, fasta, noManifest, "hg38", List.of());
        // A source with no directory for the requested reference version.
        run(dir, "wrong-ref-version", segmentsPath, fasta, hg19Only, "hg38", List.of());
        // The same folder answered for the version it does have.
        run(dir, "hg19", segmentsPath, fasta, hg19Only, "hg19", List.of());

        // A config missing a required key for its type.
        final Path broken = dataSources(dir, "ds-broken", "1.7.hg38.20220101", "hg38", true);
        final Path config = broken.resolve("regions").resolve("hg38").resolve("regions.config");
        Files.writeString(config, Files.readString(config).replace("end_column = 2\n", ""),
                StandardCharsets.UTF_8);
        run(dir, "missing-config-key", segmentsPath, fasta, broken, "hg38", List.of());

        // A segment of 101 bases, under the 150-base minimum, which is REFUSED rather than
        // written through unannotated.
        final Path shortSegments = write(dir, "short.seg", String.join("\n",
                "CONTIG\tSTART\tEND\tCALL",
                "chr1\t1200\t1300\t+",
                ""));
        run(dir, "short-segment", shortSegments, fasta, good, "hg38", List.of());

    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    static void run(final Path dir, final String label, final Path segments, final Path fasta,
                    final Path dataSources, final String refVersion, final List<String> extra)
            throws Exception {
        final Path out = dir.resolve("out-" + label + ".seg");
        final List<String> argv = new ArrayList<>(List.of(
                "--segments", segments.toString(),
                "-O", out.toString(),
                "-R", fasta.toString(),
                "--data-sources-path", dataSources.toString(),
                "--ref-version", refVersion,
                "--output-file-format", "SEG"));
        argv.addAll(extra);
        try {
            new FuncotateSegments().instanceMain(argv.toArray(new String[0]));
        } catch (final Exception | AssertionError e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            System.out.printf("error\t%s\t%s:%s%n", label, cause.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(cause.getMessage()), dir)));
            return;
        }
        if (!Files.exists(out)) {
            System.out.printf("none\t%s=no output file%n", label);
            return;
        }
        System.out.printf("out\t%s=%s%n", label,
                ReferenceQueryDump.escape(masked(Files.readString(out), dir)));
        // Every run writes a SECOND file beside the first, one row per gene and exon.
        final Path geneList = dir.resolve(out.getFileName() + ".gene_list.txt");
        if (Files.exists(geneList)) {
            System.out.printf("genes\t%s=%s%n", label,
                    ReferenceQueryDump.escape(masked(Files.readString(geneList), dir)));
        } else {
            System.out.printf("none\t%s=no gene list%n", label);
        }
    }

    static Path writeReference(final Path dir) throws Exception {
        final Path fasta = dir.resolve("reference.fasta");
        final StringBuilder bases = new StringBuilder(">chr1\n");
        for (int i = 0; i < CONTIG_LENGTH / 60; i++) {
            bases.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n");
        }
        Files.writeString(fasta, bases.toString(), StandardCharsets.UTF_8);
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(fasta, true);
        final htsjdk.samtools.SAMFileHeader header = new htsjdk.samtools.SAMFileHeader();
        header.setSequenceDictionary(new htsjdk.samtools.SAMSequenceDictionary(List.of(
                new htsjdk.samtools.SAMSequenceRecord("chr1", CONTIG_LENGTH))));
        try (final java.io.Writer writer = Files.newBufferedWriter(dir.resolve("reference.dict"))) {
            new htsjdk.samtools.SAMTextHeaderCodec().encode(writer, header);
        }
        return fasta;
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>");
    }
}
