/*
 * GenomicsDBImport through Main: its refusals, and the workspace it writes.
 *
 * The importer is native: htsjdk's VCF encoder feeds each VariantContext to a fork of htslib, whose
 * records `VCF2Binary` turns into cells that `LoaderArrayWriter` orders and TileDB writes as gzipped
 * tiles. Two runs of the reference differ in three things only, and this dump masks exactly those:
 *
 *   - A FRAGMENT'S NAME is `__<uuid><thread id>_<milliseconds>`, so it is listed as `fragment-<n>`,
 *     numbered by that timestamp.
 *   - A BATCH'S ROW BOUNDS are in `genomicsdb_meta_<uuid>.json`, listed as `genomicsdb_meta_<n>.json`,
 *     numbered by their content.
 *   - THE COMMAND LINE IN vcfheader.vcf records its `Date`.
 *
 * `__array_schema.tdb` also holds the array's absolute path, which moves with the run's directory,
 * and is listed with that length-prefixed path reduced to the array's own name.
 *
 * Everything else, every tile of every attribute included, is compared byte for byte.
 *
 * Output:
 *
 *     input\t<name>\t<base64 of a fixture file>
 *     argv\t<label>\t<the command line, tab-separated, escaped>
 *     returned\t<label>\t<what the tool returned, escaped>
 *     stdout\t<label>\t<what the tool printed, escaped>
 *     error\t<label>\t<exception class>:<message, escaped, cwd masked>
 *     dir\t<label>\t<a directory under the run's output, masked>
 *     file\t<label>\t<path under the run's output, masked>\t<base64 of its masked bytes>
 *
 * A refused run's files are not listed: the reference stops inside the native loader having written
 * part of a workspace, and which part depends on its buffers.
 *
 * Usage: GenomicsDBImportDump
 */

import org.broadinstitute.hellbender.Main;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GenomicsDBImportDump {

    static final Path DIR = Path.of("genomicsdb-dump");

    static final String HEADER = "##fileformat=VCFv4.2\n"
            + "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele\">\n"
            + "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##FILTER=<ID=FilterA,Description=\"a\">\n"
            + "##FILTER=<ID=FilterB,Description=\"b\">\n"
            + "##FILTER=<ID=zzz,Description=\"z\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">\n"
            + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"gq\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n"
            + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"ft\">\n"
            + "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"min dp\">\n"
            + "##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"pgt\">\n"
            + "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"pid\">\n"
            + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"pl\">\n"
            + "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"ps\">\n"
            + "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"sb\">\n"
            + "##FORMAT=<ID=XF,Number=R,Type=Float,Description=\"xf\">\n"
            + "##FORMAT=<ID=XI,Number=.,Type=Integer,Description=\"xi\">\n"
            + "##FORMAT=<ID=XS,Number=.,Type=String,Description=\"xs\">\n"
            + "##INFO=<ID=AS_QUALapprox,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=AS_VarDP,Number=1,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=F2,Number=2,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=FG,Number=G,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=FLG,Number=0,Type=Flag,Description=\"x\">\n"
            + "##INFO=<ID=FR,Number=R,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=MQ_DP,Number=1,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=NF,Number=.,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=NI,Number=.,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=QUALapprox,Number=1,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=RAW_GT_COUNT,Number=3,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description=\"x\">\n"
            + "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=TREE_SCORE,Number=1,Type=Float,Description=\"x\">\n"
            + "##INFO=<ID=TXT,Number=.,Type=String,Description=\"x\">\n"
            + "##INFO=<ID=VarDP,Number=1,Type=Integer,Description=\"x\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "##contig=<ID=chr2,length=50000>\n"
            + "##contig=<ID=chr3,length=2000>\n"
            + "##source=HaplotypeCaller\n";

    static final String BLOCK = "GT:DP:GQ:MIN_DP:PL";

    /** One line of a VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample. */
    static String line(final String contig, final int position, final String id, final String ref,
            final String alt, final String qual, final String filter, final String info,
            final String format, final String value) {
        return String.join("\t", contig, Integer.toString(position), id, ref, alt, qual, filter, info, format, value)
                + "\n";
    }

    static String block(final String contig, final int position, final int end, final String value) {
        return line(contig, position, ".", "A", "<NON_REF>", ".", ".", "END=" + end, BLOCK, value);
    }

    static String vcf(final String sample, final String header, final String... lines) {
        return header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
                + String.join("", lines);
    }

    static String rich() {
        return vcf("A", HEADER,
                block("chr1", 1, 99, "0/0:10:30:9:0,30,450"),
                line("chr1", 100, "rs1", "C", "T,<NON_REF>", "50.64", "PASS",
                        "AS_RAW_BaseQRankSum=|0.5,1,-1.2,2|;AS_RAW_MQ=100.00|200.00|0.00;AS_RAW_MQRankSum=|1.5,1|;"
                                + "AS_RAW_ReadPosRankSum=|-0.7,3,2.0,1|;AS_SB_TABLE=1,2|3,4|0,0;AS_QUALapprox=|50|0;"
                                + "AS_VarDP=5|7|0;BaseQRankSum=0.5;DP=12;MLEAC=1,0;MLEAF=0.500,0.00;RAW_MQandDP=43200,12;"
                                + "ReadPosRankSum=-0.3;QUALapprox=50;VarDP=12;RAW_GT_COUNT=0,1,0;TREE_SCORE=1.25",
                        "GT:AD:DP:GQ:PL:SB", "0/1:5,7,0:12:50:50,0,70,65,91,156:1,2,3,4"),
                block("chr1", 101, 200, "0/0:11:33:10:0,33,495"),
                // A deletion with a spanning allele, its span overlapped by the next two records.
                line("chr1", 201, ".", "ATT", "A,*,<NON_REF>", "99.9", "LowQual;FilterB;FilterA;zzz",
                        "DP=9;FLG;NI=1,.,3;NF=1.5,.,2.5;FR=1,2,3,4;F2=.,0.5;FG=1,2,3,4,5,6,7,8,9,10;TXT=one,two",
                        "GT:AD:DP:GQ:PGT:PID:PL:PS:SB",
                        "0|1:1,8,0,0:9:20:0|1:201_A_T:100,0,20,110,30,130,120,40,140,150:201:0,1,4,4"),
                line("chr1", 202, ".", "T", "*,<NON_REF>", "1e-3", "FilterA", "DP=4", "GT:AD:GQ:PL",
                        "1/1:0,4,0:12:120,12,0,120,12,120"),
                line("chr1", 203, ".", "T", "G,<NON_REF>", "0", ".",
                        "F2=1,2;NI=.;NF=.;TXT=.;BaseQRankSum=NaN;ReadPosRankSum=Infinity", "GT:AD:FT:XF:XI:XS",
                        "./.:1,2,3:bad:1.5,.,2.5:1,.,3:a,b"),
                line("chr1", 300, ".", "T", "<NON_REF>", ".", ".", "END=400", "GT:DP:GQ:FT", ".:1:.:PASS"),
                line("chr1", 401, ".", "T", "A,<NON_REF>", ".", ".", "MLEAC=.,1;FG=1;TXT=", "GT:XS:XI", "1:single:7"),
                line("chr1", 402, ".", "TAA", "T,TA,<NON_REF>", ".", ".", "AS_RAW_MQ=1|2|3|4;AS_SB_TABLE=1|2|3|4",
                        "GT", "0|1|2"),
                line("chr1", 410, ".", "A", "ATT,<NON_REF>", ".", ".", "AS_RAW_BaseQRankSum=NaN,1|,|;AS_QUALapprox=0x10|010|-3",
                        "GT", "0/1"),
                line("chr1", 420, ".", "AC", "A,<NON_REF>", ".", ".", ".", "GT", "1/1"),
                block("chr2", 1, 500, "0/0:8:24:7:0,24,360"),
                line("chr2", 600, ".", "A", "C,<NON_REF>", ".", ".", "DP=3;AS_VarDP=1|.|0", "GT:AD:DP", "0/1:1,2,0:3"));
    }

    static String blocks() {
        return vcf("B", HEADER,
                block("chr1", 1, 149, "0/0:10:30:9:0,30,450"),
                block("chr1", 150, 1000, "0/0:10:30:9:0,30,450"),
                line("chr1", 1001, ".", "A", "C,<NON_REF>", "12.5", ".", "DP=5", "GT:AD:DP:GQ:PL",
                        "0/1:2,3,0:5:9:9,0,9,20,30,40"),
                block("chr1", 1002, 2500, "0/0:11:33:10:0,33,495"),
                block("chr2", 400, 700, "0/0:8:24:7:0,24,360"));
    }

    /** Blocks just ahead of the partition chr1:150-1500, and one running past its end. */
    static String around() {
        return vcf("C", HEADER,
                block("chr1", 140, 145, "0/0:1:1:1:0,1,2"),
                block("chr1", 146, 149, "0/0:1:1:1:0,1,2"),
                block("chr1", 150, 150, "0/0:1:1:1:0,1,2"),
                block("chr1", 151, 151, "0/0:1:1:1:0,1,2"),
                line("chr1", 152, ".", "ATC", "A,<NON_REF>", ".", ".", "DP=2", "GT:AD", "0/1:1,1,0"),
                block("chr1", 1490, 1600, "0/0:1:1:1:0,1,2"));
    }

    /** More than two tiles' worth of cells, the blocks counting twice. */
    static String many() {
        final StringBuilder lines = new StringBuilder();
        long state = 7;
        int position = 1;
        final String[] bases = {"A", "C", "G", "T"};
        while (position < 9000) {
            state = (state * 1103515245L + 12345L) & 0x7fffffffL;
            final int a = (int) (state >> 8);
            final int depth = a % 97;
            final int likelihood = a % 2999;
            if (a % 2 == 0) {
                final int end = position + (a % 13);
                lines.append(block("chr1", position, end,
                        "0/0:" + depth + ":" + (a % 89) + ":" + (a % 83) + ":0," + likelihood + "," + (likelihood * 2)));
                position = end + 1;
            } else {
                final String ref = bases[a % 4];
                final String alt = bases[(a + 1 + a % 3) % 4];
                lines.append(line("chr1", position, ".", ref, alt + ",<NON_REF>", (a % 900) + "." + (a % 90), ".",
                        "AS_RAW_MQ=" + (a % 99991) + ".00|" + (a % 9973) + ".00|0.00;DP=" + depth + ";RAW_MQandDP="
                                + (a % 99991) + "," + depth,
                        "GT:AD:DP:GQ:PL:SB",
                        "0/1:" + (a % 50) + "," + (a % 40) + ",0:" + depth + ":" + (a % 89) + ":" + likelihood
                                + ",0," + (a % 1999) + "," + (a % 1499) + "," + (a % 999) + "," + (a % 499)
                                + ":" + (a % 7) + "," + (a % 5) + "," + (a % 3) + "," + (a % 11)));
                position += 1;
            }
        }
        return vcf("BIG", HEADER, lines.toString());
    }

    static final String ONE = line("chr1", 5, ".", "A", "C,<NON_REF>", ".", ".", ".", "GT", "0/1");

    public static void main(final String[] args) throws Exception {
        Files.createDirectories(DIR);
        input("a.g.vcf", rich());
        input("b.g.vcf", blocks());
        input("c.g.vcf", around());
        input("d.g.vcf", vcf("D", HEADER, block("chr1", 10, 5000, "0/0:1:1:1:0,1,2")));
        input("big.g.vcf", many());
        input("eq.g.vcf", vcf("S=1&<x>'q", HEADER, ONE));
        input("nopass.g.vcf", vcf("NP", HEADER.replace("##FILTER=<ID=PASS,Description=\"All filters passed\">\n", ""),
                line("chr1", 5, ".", "A", "C,<NON_REF>", ".", "PASS", ".", "GT", "0/1"),
                line("chr1", 6, ".", "A", "C,<NON_REF>", ".", "FilterA", ".", "GT", "0/1")));
        final StringBuilder quals = new StringBuilder();
        int at = 10;
        for (final String qual : List.of("0.001", "0.005", "0.0049", "1e-3", "2.5e1", "1000000.5", "0.0", "3.14159265358979",
                "1e10", "9.999", "7e-5", "12345.678")) {
            quals.append(line("chr1", at++, ".", "A", "C,<NON_REF>", qual, ".", ".", "GT", "0/1"));
        }
        input("qual.g.vcf", vcf("Q", HEADER, quals.toString()));
        input("values.g.vcf", vcf("V", HEADER,
                line("chr1", 10, ".", "A", "C,<NON_REF>", ".", ".", "TXT=a,.,b", "GT", "0/1"),
                line("chr1", 11, ".", "A", "C,<NON_REF>", ".", ".", "TXT=.,b;FLG", "GT:XS", "0/1:.,b"),
                line("chr1", 12, ".", "A", "C,<NON_REF>", ".", ".", "TXT=a,b", "GT:XS", "0/1:a,."),
                line("chr1", 13, ".", "A", "C,<NON_REF>", ".", ".", "FLG=1", "GT:GQ:DP:PL", "0/1:120:007:."),
                line("chr1", 14, ".", "A", "C,<NON_REF>", ".", ".", "NI=1;NF=2.5;MLEAC=.,.;MLEAF=.,.;FR=.,.,.", "GT:XI:XF", "0/1:.,.:.,.,."),
                line("chr1", 15, ".", "A", "C,<NON_REF>", ".", ".", "F2=.,.;RAW_MQandDP=.,.;DP=007;BaseQRankSum=1.50", "GT:XF:XI", "0/1:1.50,2e0,.:+3,04"),
                line("chr1", 16, ".", "A", "C,<NON_REF>", ".", ".", ".", "GT:DP:GQ", "0/1"),
                line("chr1", 17, ".", "A", "C,<NON_REF>", ".", ".", ".", "GT:AD:DP", "./.:.:.")));
        // A block-compressed input, read through its tabix index.
        final ByteArrayOutputStream compressed = new ByteArrayOutputStream();
        try (htsjdk.samtools.util.BlockCompressedOutputStream stream =
                new htsjdk.samtools.util.BlockCompressedOutputStream(compressed, (java.io.File) null)) {
            stream.write(blocks().getBytes(StandardCharsets.UTF_8));
        }
        inputBytes("bgz.vcf.gz", compressed.toByteArray());
        input("again.g.vcf", vcf("A", HEADER, ONE));
        input("two.g.vcf", HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"
                + "chr1\t5\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT\t0/1\t0/0\n");
        input("mnp.g.vcf", vcf("M", HEADER, line("chr1", 5, ".", "AC", "GT,<NON_REF>", ".", ".", ".", "GT", "0/1")));
        input("nocontig.g.vcf", vcf("NC", HEADER.replaceAll("##contig[^\n]*\n", "")));
        input("noindex.g.vcf", vcf("NI", HEADER, ONE));
        input("extra.g.vcf", vcf("EX", HEADER.replace("##source=HaplotypeCaller",
                "##INFO=<ID=EXTRA,Number=1,Type=Integer,Description=\"e\">\n"
                        + "##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"e\">\n##source=Other"),
                line("chr1", 5, ".", "A", "C,<NON_REF>", ".", ".", "EXTRA=3", "GT:DP2", "0/1:7")));
        input("map.txt", "zed\t" + at("b.g.vcf") + "\nalpha\t" + at("a.g.vcf") + "\t " + at("a.g.vcf.idx")
                + " \nmid dle\t" + at("c.g.vcf") + "\n");
        input("map_empty.txt", "");
        input("map_fields.txt", "a\tb\tc\td\n");
        input("map_space.txt", " a\t" + at("a.g.vcf") + "\n");
        input("map_index.txt", "a\t" + at("a.g.vcf") + "\t \n");
        input("map_twice.txt", "a\t" + at("a.g.vcf") + "\na\t" + at("b.g.vcf") + "\n");

        for (final String name : List.of("a", "b", "c", "d", "big", "eq", "again", "two", "mnp", "extra", "nopass", "qual", "values")) {
            final String[] index = {"IndexFeatureFile", "-I", at(name + ".g.vcf")};
            new Main().instanceMain(index);
        }
        new Main().instanceMain(new String[] {"IndexFeatureFile", "-I", at("bgz.vcf.gz")});

        final String a = at("a.g.vcf");
        final String b = at("b.g.vcf");
        final String c = at("c.g.vcf");
        final String d = at("d.g.vcf");

        // The workspace.
        ws("three-samples", "-V", a, "-V", b, "-V", c, "-L", "chr1", "-L", "chr2");
        ws("batches", "-V", a, "-V", b, "-V", c, "-V", d, "-L", "chr1", "-L", "chr2:1-3000", "--batch-size", "3");
        ws("one-per-batch", "-V", a, "-V", b, "-L", "chr3", "-L", "chr1:1-10", "--batch-size", "1");
        ws("across-the-begin", "-V", a, "-V", b, "-V", c, "-V", d, "-L", "chr1:150-1500", "--batch-size", "2");
        ws("single-column", "-V", a, "-V", b, "-V", c, "-V", d, "-L", "chr1:151-151", "-L", "chr1:202-205");
        ws("only-a-spanning-block", "-V", d, "-L", "chr1:2000-2100");
        ws("several-tiles", "-V", at("big.g.vcf"), "-V", a, "-L", "chr1");
        ws("small-segment", "-V", at("big.g.vcf"), "-L", "chr1", "--genomicsdb-segment-size", "4096");
        ws("small-buffer", "-V", a, "-L", "chr1", "--genomicsdb-vcf-buffer-size", "1024");
        ws("sample-map", "--sample-name-map", at("map.txt"), "-L", "chr1:1-2000");
        ws("sample-map-validated", "--sample-name-map", at("map.txt"), "-L", "chr1:1-2000", "--validate-sample-name-map");
        ws("merged-intervals", "-V", a, "-V", b, "-L", "chr1:10-20", "-L", "chr1:100-200", "-L", "chr2:5-6", "--merge-input-intervals");
        ws("padding", "-V", a, "-L", "chr1:100-110", "-ip", "5", "-L", "chr1:300-300");
        ws("excluded", "-V", a, "-L", "chr1", "-XL", "chr1:50-60");
        ws("shared-posix", "-V", a, "-V", b, "-L", "chr1", "--genomicsdb-shared-posixfs-optimizations");
        ws("reader-threads", "-V", a, "-V", b, "-L", "chr1", "--reader-threads", "2");
        ws("reader-threads-two-intervals", "-V", a, "-V", b, "-L", "chr1", "-L", "chr2", "--reader-threads", "2",
                "--max-num-intervals-to-import-in-parallel", "2");
        ws("merged-fields", "-V", a, "-V", at("extra.g.vcf"), "-L", "chr1:1-100");
        ws("escaped-name", "-V", at("eq.g.vcf"), "-L", "chr1");
        ws("no-pass-line", "-V", at("nopass.g.vcf"), "-L", "chr1");
        ws("quals", "-V", at("qual.g.vcf"), "-L", "chr1");
        ws("values", "-V", at("values.g.vcf"), "-L", "chr1");
        ws("compressed", "-V", at("bgz.vcf.gz"), "-L", "chr1", "-L", "chr2");
        run("workspace-there", "-V", b, "--genomicsdb-workspace-path", out("workspace-there"), "-L", "chr1:1-10");
        run("workspace-not-overwritten", "-V", b, "--genomicsdb-workspace-path", out("workspace-there"), "-L", "chr1:1-10");
        run("workspace-overwritten", "-V", b, "--genomicsdb-workspace-path", out("workspace-there"), "-L", "chr1:1-20",
                "--overwrite-existing-genomicsdb-workspace");
        run("workspace-inside-a-workspace", "-V", b, "--genomicsdb-workspace-path", out("workspace-there") + "/inner", "-L", "chr1:1-10");

        // The arguments.
        run("no-input", "--genomicsdb-workspace-path", out("no-input"), "-L", "chr1");
        run("no-intervals", "-V", a, "--genomicsdb-workspace-path", out("no-intervals"));
        run("only-excluded", "-V", a, "--genomicsdb-workspace-path", out("only-excluded"), "-XL", "chr1:50-60");
        run("no-workspace", "-V", a, "-L", "chr1");
        run("two-workspaces", "-V", a, "--genomicsdb-workspace-path", out("two-workspaces"),
                "--genomicsdb-update-workspace-path", out("two-workspaces"), "-L", "chr1");
        run("variant-and-map", "-V", a, "--sample-name-map", at("map.txt"), "--genomicsdb-workspace-path",
                out("variant-and-map"), "-L", "chr1");
        run("overwrite-an-update", "-V", a, "--genomicsdb-update-workspace-path", out("overwrite-an-update"),
                "--overwrite-existing-genomicsdb-workspace", "-L", "chr1");
        run("avoid-nio", "--sample-name-map", at("map.txt"), "--genomicsdb-workspace-path", out("avoid-nio"),
                "-L", "chr1", "--avoid-nio");
        run("avoid-nio-no-header", "--sample-name-map", at("map.txt"), "--genomicsdb-workspace-path",
                out("avoid-nio-no-header"), "-L", "chr1", "--avoid-nio", "--bypass-feature-reader");
        run("bypass-plain-vcf", "-V", a, "--genomicsdb-workspace-path", out("bypass-plain-vcf"), "-L", "chr1",
                "--bypass-feature-reader");
        run("bypass-with-indices", "--sample-name-map", at("map.txt"), "--genomicsdb-workspace-path",
                out("bypass-with-indices"), "-L", "chr1", "--bypass-feature-reader");

        // The inputs.
        run("repeated-sample", "-V", a, "-V", at("again.g.vcf"), "--genomicsdb-workspace-path",
                out("repeated-sample"), "-L", "chr1");
        run("two-samples", "-V", at("two.g.vcf"), "--genomicsdb-workspace-path", out("two-samples"), "-L", "chr1");
        run("no-dictionary", "-V", at("nocontig.g.vcf"), "--genomicsdb-workspace-path", out("no-dictionary"), "-L", "chr1");
        run("unknown-contig", "-V", a, "--genomicsdb-workspace-path", out("unknown-contig"), "-L", "chr77");
        run("no-such-input", "-V", at("none.g.vcf"), "--genomicsdb-workspace-path", out("no-such-input"), "-L", "chr1");
        run("no-index", "-V", at("noindex.g.vcf"), "--genomicsdb-workspace-path", out("no-index"), "-L", "chr1");
        run("an-mnp", "-V", at("mnp.g.vcf"), "--genomicsdb-workspace-path", out("an-mnp"), "-L", "chr1");
        ws("an-mnp-outside", "-V", at("mnp.g.vcf"), "-L", "chr2");
        run("map-missing", "--sample-name-map", at("none.txt"), "--genomicsdb-workspace-path", out("map-missing"), "-L", "chr1");
        run("map-empty", "--sample-name-map", at("map_empty.txt"), "--genomicsdb-workspace-path", out("map-empty"), "-L", "chr1");
        run("map-fields", "--sample-name-map", at("map_fields.txt"), "--genomicsdb-workspace-path", out("map-fields"), "-L", "chr1");
        run("map-space", "--sample-name-map", at("map_space.txt"), "--genomicsdb-workspace-path", out("map-space"), "-L", "chr1");
        run("map-index", "--sample-name-map", at("map_index.txt"), "--genomicsdb-workspace-path", out("map-index"), "-L", "chr1");
        run("map-twice", "--sample-name-map", at("map_twice.txt"), "--genomicsdb-workspace-path", out("map-twice"), "-L", "chr1");

        // The workspace that is already there.
        run("workspace-exists", "-V", a, "--genomicsdb-workspace-path", DIR.toString(), "-L", "chr1");
        run("workspace-is-a-file", "-V", a, "--genomicsdb-workspace-path", at("map.txt"), "-L", "chr1");
        run("workspace-no-parent", "-V", a, "--genomicsdb-workspace-path", out("workspace-no-parent") + "/no/such", "-L", "chr1");

        // The importer leaves threads behind that instanceMain, unlike main, does not exit past.
        System.out.flush();
        System.exit(0);
    }

    static String at(final String name) {
        return DIR.resolve(name).toString();
    }

    static String out(final String label) {
        return DIR.resolve("out").resolve(label).resolve("ws").toString();
    }

    static void input(final String name, final String content) throws Exception {
        inputBytes(name, content.getBytes(StandardCharsets.UTF_8));
    }

    static void inputBytes(final String name, final byte[] bytes) throws Exception {
        Files.write(DIR.resolve(name), bytes);
        System.out.printf("input\t%s\t%s%n", name, Base64.getEncoder().encodeToString(bytes));
    }

    /** A run with the workspace argument added, its label naming the directory. */
    static void ws(final String label, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(Arrays.asList(rest));
        argv.add("--genomicsdb-workspace-path");
        argv.add(out(label));
        run(label, argv.toArray(new String[0]));
    }

    static void run(final String label, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>();
        argv.add("GenomicsDBImport");
        argv.addAll(Arrays.asList(rest));
        System.out.printf("argv\t%s\t%s%n", label, ReferenceQueryDump.escape(String.join("\t", argv)));
        final Path outDir = DIR.resolve("out").resolve(label);
        Files.createDirectories(outDir);
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        String failure = null;
        Object returned = null;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            try {
                returned = new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception e) {
                Throwable cause = e;
                while (cause.getCause() != null && cause.getCause() != cause) {
                    cause = cause.getCause();
                }
                failure = e.getClass().getName() + ":" + String.valueOf(cause.getMessage());
            }
        } finally {
            System.out.flush();
            System.setOut(realOut);
        }
        System.out.printf("stdout\t%s\t%s%n", label,
                ReferenceQueryDump.escape(masked(outBytes.toString(StandardCharsets.UTF_8))));
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(masked(failure)));
            return;
        }
        System.out.printf("returned\t%s\t%s%n", label, ReferenceQueryDump.escape(String.valueOf(returned)));
        final Path workspace = outDir.resolve("ws");
        if (!Files.isDirectory(workspace)) {
            return;
        }
        final List<Path> produced;
        try (Stream<Path> walk = Files.walk(workspace)) {
            produced = walk.filter(path -> !path.equals(workspace)).sorted().collect(Collectors.toList());
        }
        final List<String> rows = new ArrayList<>();
        for (final Path path : produced) {
            final String name = stableName(workspace, path);
            if (Files.isDirectory(path)) {
                rows.add(String.format("dir\t%s\t%s", label, name));
            } else {
                rows.add(String.format("file\t%s\t%s\t%s", label, name,
                        Base64.getEncoder().encodeToString(stableBytes(path))));
            }
        }
        rows.sort(Comparator.naturalOrder());
        rows.forEach(System.out::println);
    }

    /** The path relative to the workspace with each run's own names numbered. */
    static String stableName(final Path workspace, final Path path) throws Exception {
        final List<String> parts = new ArrayList<>();
        Path walked = workspace;
        for (final Path part : workspace.relativize(path)) {
            walked = walked.resolve(part);
            final String name = part.toString();
            if (name.startsWith("__") && Files.isDirectory(walked)) {
                parts.add("fragment-" + indexAmong(walked, fragments(walked.getParent())));
            } else if (name.startsWith("genomicsdb_meta_") && !name.equals("genomicsdb_meta_dir")
                    && !name.equals("genomicsdb_meta_column_bounds.json")) {
                parts.add("genomicsdb_meta_" + indexAmong(walked, metas(walked.getParent())) + ".json");
            } else {
                parts.add(name);
            }
        }
        return String.join("/", parts);
    }

    static int indexAmong(final Path one, final List<Path> all) {
        return all.indexOf(one);
    }

    static List<Path> fragments(final Path array) throws Exception {
        try (Stream<Path> list = Files.list(array)) {
            return list.filter(path -> Files.isDirectory(path) && path.getFileName().toString().startsWith("__"))
                    .sorted(Comparator.comparingLong(path -> {
                        final String name = path.getFileName().toString();
                        return Long.parseLong(name.substring(name.lastIndexOf('_') + 1));
                    }))
                    .collect(Collectors.toList());
        }
    }

    static List<Path> metas(final Path directory) throws Exception {
        try (Stream<Path> list = Files.list(directory)) {
            return list.filter(path -> path.getFileName().toString().startsWith("genomicsdb_meta_"))
                    .sorted(Comparator.comparing(path -> {
                        try {
                            return new String(Files.readAllBytes(path), StandardCharsets.UTF_8);
                        } catch (final Exception e) {
                            throw new IllegalStateException(e);
                        }
                    }))
                    .collect(Collectors.toList());
        }
    }

    static byte[] stableBytes(final Path path) throws Exception {
        final byte[] raw = Files.readAllBytes(path);
        final String name = path.getFileName().toString();
        if (name.equals("vcfheader.vcf")) {
            return new String(raw, StandardCharsets.UTF_8).replaceAll("Date=\"[^\"]*\"", "Date=\"<masked>\"")
                    .getBytes(StandardCharsets.UTF_8);
        }
        if (name.equals("__array_schema.tdb")) {
            final int length = ByteBuffer.wrap(raw, 8, 4).order(ByteOrder.LITTLE_ENDIAN).getInt();
            final String array = new String(raw, 12, length, StandardCharsets.UTF_8);
            final byte[] base = array.substring(array.lastIndexOf('/') + 1).getBytes(StandardCharsets.UTF_8);
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            out.write(raw, 0, 8);
            out.write(base);
            out.write(raw, 12 + length, raw.length - 12 - length);
            return out.toByteArray();
        }
        return raw;
    }

    static String masked(final String text) {
        return text.replace(Path.of("").toAbsolutePath() + "/", "<cwd>/");
    }
}
