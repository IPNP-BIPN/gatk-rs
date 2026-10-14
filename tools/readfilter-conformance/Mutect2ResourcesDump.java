/*
 * Mutect2's germline resource and panel of normals, end to end (G3.2r, #1352).
 *
 * --germline-resource: each call's POPAF is minus log10 of the AF of the resource record whose
 * alleles match it once both are against the longer reference (GATKVariantContextUtils.alleleIndices),
 * the default allele frequency otherwise, a record whose AF is not one value per alternate skipped;
 * a reference of the same length but other bases is refused. In tumour-only calling, a locus with a
 * resource allele of AF at least --max-population-af is left inactive when that allele is an indel
 * and the likeliest tumour allele is one, or a substitution and the likeliest tumour base is the
 * resource REFERENCE's first base (what Mutect2Engine.isActive compares), unless
 * --genotype-germline-sites. --panel-of-normals: its loci are left inactive unless
 * --genotype-pon-sites, and a call starting on one is flagged PON.
 *
 * Read set: HcDragenEndToEndDump's, as Mutect2MitochondriaDump builds it. The resources are written
 * from a first plain run's calls.
 *
 * Output:
 *
 *     fasta\t<contig>=<the contig's bases>
 *     sam\tinput=<the SAM text of the one-sample BAM, escaped>
 *     sam\tboth=<the SAM text of the two-sample BAM, escaped>
 *     vcf\t<name>=<a resource VCF, escaped>
 *     out\t<label>=<the output's header lines but the command line, and its body, escaped>
 *     error\t<label>\t<exception class>:<message, the dump's directory as <dir>>
 *
 * Usage: Mutect2ResourcesDump
 */

import org.broadinstitute.hellbender.tools.walkers.mutect.Mutect2;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Mutect2ResourcesDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("mutect2-resources-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# Mutect2ResourcesDump: Mutect2's germline resource and panel of normals");

        // HcDragenEndToEndDump's contig, built as Mutect2MitochondriaDump builds it.
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

        final String both = HcDragenEndToEndDump.sam(reference);
        final String sam = java.util.Arrays.stream(both.split("\n"))
                .filter(line -> !line.contains("rg2"))
                .collect(java.util.stream.Collectors.joining("\n")) + "\n";
        final Path twoSamples = bam(dir, "both", both);
        final Path bam = bam(dir, "input", sam);
        System.out.printf("sam\tinput=%s%n", ReferenceQueryDump.escape(sam));
        System.out.printf("sam\tboth=%s%n", ReferenceQueryDump.escape(both));

        // The first plain run, whose calls the resources are written from.
        final Path plain = dir.resolve("plain.vcf");
        new Mutect2().instanceMain(new String[] {"-I", bam.toString(), "-O", plain.toString(), "-R", fasta.toString()});
        final List<String[]> calls = new ArrayList<>();
        for (final String line : Files.readAllLines(plain)) {
            if (!line.startsWith("#")) {
                calls.add(line.split("\t"));
            }
        }
        final String contig = HcDragenEndToEndDump.CONTIG;
        final StringBuilder germline = new StringBuilder();
        final StringBuilder pon = new StringBuilder();
        final StringBuilder badRef = new StringBuilder();
        final StringBuilder wrongLength = new StringBuilder();
        // Over the plain run's four calls (a SNV, a deletion, an insertion and a SNV), the resource
        // records come in file order and the last that matches a call wins.
        for (int i = 0; i < calls.size(); i++) {
            final String[] call = calls.get(i);
            final int pos = Integer.parseInt(call[1]);
            final String callRef = call[3];
            final String alt = call[4].split(",")[0];
            switch (i) {
                case 0 -> {
                    // An AF that is not one value per alternate, skipped (both values rare, so isActive
                    // never reads past the alternates); then a common substitution, which isActive does
                    // not exclude since it compares the resource's REF base.
                    germline.append(record(contig, pos, callRef, alt, "AF=0.001,0.002"));
                    germline.append(record(contig, pos, callRef, alt, "AF=0.5"));
                }
                case 1 -> {
                    germline.append(record(contig, pos, callRef, alt, "AF=."));
                    germline.append(record(contig, pos, callRef, alt, "AF=0.001"));
                }
                // A common indel: in tumour-only calling its locus is left inactive.
                case 2 -> germline.append(record(contig, pos, callRef, alt, "AF=0.5"));
                case 3 -> {
                    germline.append(record(contig, pos, callRef, alt + "," + otherThan(callRef, alt), "AF=0.2,0.003"));
                    // A longer representation of the same event: both alleles extended by the next base.
                    final String next = reference.substring(pos - 1 + callRef.length(), pos + callRef.length());
                    germline.append(record(contig, pos, callRef + next, alt + next, "AF=0.004"));
                }
                default -> germline.append(record(contig, pos, callRef, alt, "AF=0.05"));
            }
            if (i % 2 == 0) {
                pon.append(record(contig, pos, callRef, alt, "."));
            }
            if (i == 0) {
                // The same, common: isActive reads an alternate the record does not have.
                wrongLength.append(record(contig, pos, callRef, alt, "AF=0.1,0.2"));
                badRef.append(record(contig, pos, otherThan(callRef, alt) + callRef.substring(1), alt, "AF=0.3"));
            }
        }
        final Path germlineVcf = vcf(dir, "germline", germline.toString());
        final Path ponVcf = vcf(dir, "pon", pon.toString());
        final Path badRefVcf = vcf(dir, "bad-ref", badRef.toString());
        final Path wrongLengthVcf = vcf(dir, "wrong-length", wrongLength.toString());
        final String g = germlineVcf.toString();
        final String p = ponVcf.toString();

        run(dir, "plain", bam, fasta, List.of());
        run(dir, "germline", bam, fasta, List.of("--germline-resource", g));
        run(dir, "germline-max-af", bam, fasta, List.of("--germline-resource", g, "--max-population-af", "0.9"));
        run(dir, "germline-sites", bam, fasta, List.of("--germline-resource", g, "--genotype-germline-sites"));
        run(dir, "germline-default-af", bam, fasta, List.of("--germline-resource", g, "--af-of-alleles-not-in-resource", "1e-4"));
        run(dir, "germline-tumor-normal", twoSamples, fasta, List.of("--germline-resource", g, "--normal-sample", "s2"));
        run(dir, "germline-gvcf", bam, fasta, List.of("--germline-resource", g, "-ERC", "GVCF", "-L", contig + ":380-460"));
        run(dir, "pon", bam, fasta, List.of("--panel-of-normals", p));
        run(dir, "pon-genotyped", bam, fasta, List.of("--panel-of-normals", p, "--genotype-pon-sites"));
        run(dir, "both", bam, fasta, List.of("--germline-resource", g, "--panel-of-normals", p, "--genotype-pon-sites",
                "--max-population-af", "0.6"));
        run(dir, "bad-ref", bam, fasta, List.of("--germline-resource", badRefVcf.toString()));
        run(dir, "wrong-length", bam, fasta, List.of("--germline-resource", wrongLengthVcf.toString()));
        run(dir, "wrong-length-max-af", bam, fasta, List.of("--germline-resource", wrongLengthVcf.toString(),
                "--max-population-af", "0.5"));
        run(dir, "missing", bam, fasta, List.of("--germline-resource", dir.resolve("absent.vcf").toString()));
    }

    /** A base that is neither the reference's first nor the alternate. */
    static String otherThan(final String ref, final String alt) {
        for (final String base : List.of("A", "C", "G", "T")) {
            if (!ref.startsWith(base) && !alt.equals(base)) {
                return base;
            }
        }
        throw new IllegalStateException("no base left");
    }

    static String record(final String contig, final int pos, final String ref, final String alt, final String info) {
        return contig + "\t" + pos + "\t.\t" + ref + "\t" + alt + "\t.\t.\t" + info + "\n";
    }

    /** A resource VCF with its Tribble index, printed. */
    static Path vcf(final Path dir, final String name, final String records) throws Exception {
        final String text = "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n"
                + "##contig=<ID=" + HcDragenEndToEndDump.CONTIG + ",length=" + HcDragenEndToEndDump.LENGTH + ">\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + records;
        final Path path = dir.resolve(name + ".vcf");
        Files.writeString(path, text, StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.createDynamicIndex(path, new htsjdk.variant.vcf.VCFCodec(),
                        htsjdk.tribble.index.IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME)
                .write(Paths.get(path + ".idx"));
        System.out.printf("vcf\t%s=%s%n", name, ReferenceQueryDump.escape(text));
        return path;
    }

    static Path bam(final Path dir, final String name, final String sam) throws Exception {
        final Path samPath = Mutect2MitochondriaDump.writeSam(dir, name + ".sam", sam);
        final Path bam = dir.resolve(name + ".bam");
        try (SamReader reader = SamReaderFactory.makeDefault().open(samPath);
             SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final SAMRecord record : reader) {
                writer.addAlignment(record);
            }
        }
        return bam;
    }

    static void run(final Path dir, final String label, final Path bam, final Path fasta, final List<String> extra)
            throws Exception {
        Mutect2MitochondriaDump.run(dir, label, Mutect2::new, bam, fasta, extra);
    }
}
