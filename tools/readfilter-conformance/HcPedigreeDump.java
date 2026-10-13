/*
 * The pedigree arguments of the callers' annotation plugin, end to end (G3.2k, #1343).
 *
 * --pedigree and --founder-id reach the pedigree annotations: ExcessHet and InbreedingCoeff count
 * only the founders' genotypes (--founder-id's, then the PED file's), and PossibleDeNovo and
 * TransmittedSingleton walk the PED file's trios. Either argument without a pedigree annotation
 * resolved is refused with the annotations named in the plugin's HashMap order.
 *
 * The reads and reference are HcDragenEndToEndDump's, with a third sample, s3, a copy of s2's
 * reads under other names: a trio whose child, s1, carries what its parents do not.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text the BAM was written from, escaped>
 *     ped\t<label>=<a PED file's text, escaped>
 *     out\t<label>=<the output's ##INFO and ##FORMAT lines and its body from #CHROM on, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: HcPedigreeDump
 */

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller;
import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class HcPedigreeDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("hc-pedigree-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# HcPedigreeDump: the pedigree arguments of the callers' annotations");

        // HcDragenEndToEndDump's contig, built the same way from the same seed.
        final Random random = new Random(20261013L);
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < HcDragenEndToEndDump.LENGTH; i++) {
            ref.append("ACGT".charAt(random.nextInt(4)));
        }
        HcDragenEndToEndDump.splice(ref, 401, "AAAAAAAAAA");
        HcDragenEndToEndDump.splice(ref, 701, "CACACACACACACA");
        HcDragenEndToEndDump.splice(ref, 1001, "TTTTTTTT");
        ref.setCharAt(399, 'G');
        ref.setCharAt(410, 'G');
        ref.setCharAt(699, 'T');
        ref.setCharAt(714, 'T');
        ref.setCharAt(999, 'G');
        ref.setCharAt(1008, 'G');
        final String reference = ref.toString();
        System.out.printf("fasta\t%s=%s%n", HcDragenEndToEndDump.CONTIG, reference);
        final Path fasta = HcDragenEndToEndDump.writeReference(dir, reference);

        // s3: s2's reads again, each right after its original so the coordinate order holds.
        final StringBuilder sam = new StringBuilder();
        for (final String line : HcDragenEndToEndDump.sam(reference).split("\n")) {
            sam.append(line).append("\n");
            if (line.equals("@RG\tID:rg2\tSM:s2")) {
                sam.append("@RG\tID:rg3\tSM:s3\n");
            } else if (line.contains("\tRG:Z:rg2\t")) {
                sam.append("t").append(line.replace("\tRG:Z:rg2\t", "\tRG:Z:rg3\t")).append("\n");
            }
        }
        final Path samPath = dir.resolve("input.sam");
        Files.writeString(samPath, sam.toString(), StandardCharsets.UTF_8);
        final Path bam = dir.resolve("input.bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(samPath);
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam.toString()));

        final Path trio = ped(dir, "trio", "fam1 s1 s3 s2 1 2\nfam1 s2 0 0 2 1\nfam1 s3 0 0 1 1\n");
        // Two families, a child with one parent known, a comment and a header-less phenotype.
        final Path mixed = ped(dir, "mixed",
                "# a comment\nfamB s2 0 0 2 -9\nfamA s1 s3 s2 2 0\nfamA s3 0 0 1 1\nfamB s9 0 s2 0 2\n");
        final Path malformed = ped(dir, "malformed", "fam1 s1 s3 s2 1\n");
        final Path unknownSex = ped(dir, "unknown-sex", "fam1 s1 s3 s2 7 2\n");

        final List<String> deNovo = List.of("-A", "PossibleDeNovo", "-A", "TransmittedSingleton");
        HcAllAnnotationsDump.run(dir, "hc-trio", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", trio.toString()), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-trio-thresholds", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", trio.toString(), "--denovo-parent-gq-threshold", "99",
                        "--denovo-depth-threshold", "40"), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-trio-default-annotations", HaplotypeCaller::new, bam, fasta,
                List.of("-ped", trio.toString()));
        HcAllAnnotationsDump.run(dir, "hc-founder-id", HaplotypeCaller::new, bam, fasta,
                List.of("--founder-id", "s2", "--founder-id", "s3"));
        HcAllAnnotationsDump.run(dir, "hc-founder-id-one", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--founder-id", "s1"), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-founder-id-and-pedigree", HaplotypeCaller::new, bam, fasta,
                List.of("--founder-id", "s1", "--pedigree", mixed.toString()));
        HcAllAnnotationsDump.run(dir, "hc-mixed", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", mixed.toString()), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-gvcf-trio", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", trio.toString(), "-ERC", "GVCF", "--sample-name", "s1"), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-malformed", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", malformed.toString()), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-unknown-sex", HaplotypeCaller::new, bam, fasta,
                concat(List.of("--pedigree", unknownSex.toString()), deNovo));
        HcAllAnnotationsDump.run(dir, "hc-no-pedigree-annotation", HaplotypeCaller::new, bam, fasta,
                List.of("--pedigree", trio.toString(), "--disable-tool-default-annotations"));
        HcAllAnnotationsDump.run(dir, "mutect2-pedigree", Mutect2::new, bam, fasta,
                List.of("--pedigree", trio.toString()));
        HcAllAnnotationsDump.run(dir, "mutect2-founder-id", Mutect2::new, bam, fasta,
                List.of("--founder-id", "s2"));
        HcAllAnnotationsDump.run(dir, "mutect2-trio", Mutect2::new, bam, fasta,
                concat(List.of("--pedigree", trio.toString(), "-A", "ExcessHet"), deNovo));
    }

    static Path ped(final Path dir, final String label, final String text) throws Exception {
        final Path path = dir.resolve(label + ".ped");
        Files.writeString(path, text, StandardCharsets.UTF_8);
        System.out.printf("ped\t%s=%s%n", label, ReferenceQueryDump.escape(text));
        return path;
    }

    static List<String> concat(final List<String> a, final List<String> b) {
        final List<String> out = new ArrayList<>(a);
        out.addAll(b);
        return out;
    }
}
