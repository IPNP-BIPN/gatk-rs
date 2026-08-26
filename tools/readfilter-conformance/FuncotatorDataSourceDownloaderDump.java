/*
 * FuncotatorDataSourceDownloader's copy, taken from the reference.
 *
 * A tar.gz fetched from a bucket, checked against a sha256 beside it, and optionally unpacked. The
 * bucket is never reached here: the tool's own testing override points it at local files, which is
 * what makes everything but the transport measurable.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE DESTINATION WITH NO -O IS THE WORKING DIRECTORY, not the source's directory:
 *     getOutputLocation falls back to `dataSourcesPath.getFileName()`, a RELATIVE path, so the
 *     copy lands wherever the process happens to be running;
 *   - THE CHECKSUM FILE IS CLEANED BEFORE COMPARISON: trimmed, lower-cased, then truncated at the
 *     first SPACE and then at the first TAB, in that order and on the already truncated string, so
 *     `<sum>  <name>` and `<sum>\t<name>` both work and an upper-case sum matches a lower-case one;
 *   - ONLY THE FIRST LINE OF THE CHECKSUM FILE IS READ, and a file with no lines at all is a
 *     UserException naming the path as a URI;
 *   - A MISMATCHED CHECKSUM NAMES BOTH SUMS, and the file is left on disk: validation happens
 *     AFTER the copy, so a corrupt download is reported rather than removed;
 *   - THE COPY REFUSES AN EXISTING DESTINATION unless --overwrite-output-file, and the refusal
 *     comes from the copier rather than from the tool;
 *   - --extract-after-download UNPACKS BESIDE THE ARCHIVE, into the destination's PARENT, and
 *     obeys the same overwrite flag;
 *   - THE TWO TESTING ARGUMENTS ARE ALL-OR-NOTHING: either alone is a UserException, and both are
 *     mutex with --somatic and --germline;
 *   - WITHOUT A DATA SOURCE OR A REFERENCE THE TOOL REFUSES IN onStartup, with two different
 *     messages in a fixed order;
 *   - AND THE BUCKET PATHS ARE BUILT FROM A VERSION STRING THAT ENCODES THE REFERENCE: the same
 *     major, minor and date, with 38 or 19 in the middle, which is what tells the two apart.
 *
 * Output:
 *
 *     constant\t<name>=<value>
 *     path\t<somatic|germline>\t<hg38|hg19>\t<data sources path>\t<checksum path>
 *     fixture\t<name>=<value>
 *     clean\t<label>\tin=<raw checksum line>\tout=<what the tool compared>
 *     run\t<label>\tcode=<the value doWork returned>
 *     file\t<label>\t<path>=<sha256 of the file that landed there>
 *     extracted\t<label>\t<path>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FuncotatorDataSourceDownloaderDump
 */

import org.broadinstitute.hellbender.tools.funcotator.FuncotatorDataSourceDownloader;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.DataSourceUtils;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

public class FuncotatorDataSourceDownloaderDump {

    public static void main(final String[] args) throws Exception {
        final Path dir = Path.of("funcotator-downloader-dump").toAbsolutePath();
        PrintReadsDump.emptyDirectory(dir);
        Files.createDirectories(dir);

        System.out.println("# FuncotatorDataSourceDownloaderDump: a data source copied and checked");

        // The constants the bucket paths are built from, and the four paths themselves. Only the
        // somatic ones are public; the germline pair is rebuilt here from the same public pieces,
        // which is also what shows that the reference number sits INSIDE the version string.
        for (final String[] constant : new String[][] {
                {"bucket", DataSourceUtils.DATA_SOURCES_BUCKET_PATH},
                {"prefix", DataSourceUtils.DATA_SOURCES_NAME_PREFIX},
                {"somatic-modifier", DataSourceUtils.DS_SOMATIC_NAME_MODIFIER},
                {"germline-modifier", DataSourceUtils.DS_GERMLINE_NAME_MODIFIER},
                {"extension", DataSourceUtils.DS_EXTENSION},
                {"checksum-extension", DataSourceUtils.DS_CHECKSUM_EXTENSION},
                {"max-version-38", DataSourceUtils.getDataSourceMaxVersionString(38)},
                {"max-version-19", DataSourceUtils.getDataSourceMaxVersionString(19)},
                {"min-version", DataSourceUtils.getDataSourceMinVersionString()},
                {"current-maximum", DataSourceUtils.CURRENT_MAXIMUM_DATA_SOURCE_VERSION},
                {"current-minimum", DataSourceUtils.CURRENT_MINIMUM_DATA_SOURCE_VERSION}}) {
            System.out.printf("constant\t%s=%s%n", constant[0], constant[1]);
        }
        for (final String kind : new String[] {"somatic", "germline"}) {
            for (final int reference : new int[] {38, 19}) {
                final String base = DataSourceUtils.DATA_SOURCES_BUCKET_PATH
                        + DataSourceUtils.DATA_SOURCES_NAME_PREFIX + "."
                        + DataSourceUtils.getDataSourceMaxVersionString(reference)
                        + (kind.equals("somatic")
                                ? DataSourceUtils.DS_SOMATIC_NAME_MODIFIER
                                : DataSourceUtils.DS_GERMLINE_NAME_MODIFIER);
                System.out.printf("path\t%s\thg%d\t%s\t%s%n", kind, reference,
                        base + DataSourceUtils.DS_EXTENSION,
                        base + DataSourceUtils.DS_CHECKSUM_EXTENSION);
            }
        }
        // The public constant, which is what the four paths above have to agree with.
        System.out.printf("constant\tsomatic-hg38-path=%s%n",
                FuncotatorDataSourceDownloader.HG38_SOMATIC_GCLOUD_DATASOURCES_PATH);

        // The archive that stands in for a data source bundle: one file, so extraction is
        // observable without being a fixture of its own.
        final Path archive = dir.resolve("funcotator_dataSources.tar.gz");
        writeArchive(archive, "MANIFEST.txt", "Version:          v1.8.hg38.20230908\n");
        final String sum = sha256(archive);
        System.out.printf("fixture\tarchive=%s%n", archive.getFileName());
        System.out.printf("fixture\tarchive-sha256=%s%n", sum);
        System.out.printf("fixture\tarchive-bytes=%d%n", Files.size(archive));

        // The checksum file, in four spellings the cleaner has to reduce to the same sum.
        final Path plain = write(dir, "plain.sha256", sum + "\n");
        final Path withName = write(dir, "with-name.sha256", sum + "  funcotator_dataSources.tar.gz\n");
        final Path tabbed = write(dir, "tabbed.sha256", sum + "\tfuncotator_dataSources.tar.gz\n");
        final Path upper = write(dir, "upper.sha256", "  " + sum.toUpperCase() + " \n");
        final Path wrong = write(dir, "wrong.sha256",
                "0000000000000000000000000000000000000000000000000000000000000000\n");
        final Path empty = write(dir, "empty.sha256", "");
        for (final Path path : List.of(plain, withName, tabbed, upper, wrong, empty)) {
            System.out.printf("clean\t%s\tin=%s%n", path.getFileName(),
                    ReferenceQueryDump.escape(Files.readString(path)));
        }

        // The copier does not create the directory it is pointed at, so a destination whose
        // parent is missing is a refusal rather than a mkdir. Measured once, then the directories
        // are made so the rest of the runs are about the tool and not about the path.
        run(dir, "missing-directory", archive, plain, List.of(
                "-O", dir.resolve("nowhere/copied.tar.gz").toString()));
        Files.createDirectories(dir.resolve("out"));
        Files.createDirectories(dir.resolve("unpack"));

        // A copy with an explicit destination, which is the simplest shape.
        run(dir, "explicit-output", archive, plain, List.of(
                "-O", dir.resolve("out/copied.tar.gz").toString()));

        // The same, with the integrity check, in each of the four spellings.
        run(dir, "validate-plain", archive, plain, List.of(
                "-O", dir.resolve("out/plain.tar.gz").toString(), "--validate-integrity", "true"));
        run(dir, "validate-with-name", archive, withName, List.of(
                "-O", dir.resolve("out/named.tar.gz").toString(), "--validate-integrity", "true"));
        run(dir, "validate-tabbed", archive, tabbed, List.of(
                "-O", dir.resolve("out/tabbed.tar.gz").toString(), "--validate-integrity", "true"));
        run(dir, "validate-upper", archive, upper, List.of(
                "-O", dir.resolve("out/upper.tar.gz").toString(), "--validate-integrity", "true"));

        // A checksum that does not match, and one that is not there at all.
        run(dir, "validate-wrong", archive, wrong, List.of(
                "-O", dir.resolve("out/wrong.tar.gz").toString(), "--validate-integrity", "true"));
        run(dir, "validate-empty", archive, empty, List.of(
                "-O", dir.resolve("out/empty.tar.gz").toString(), "--validate-integrity", "true"));

        // An existing destination, with and without the overwrite flag.
        final Path existing = dir.resolve("out/existing.tar.gz");
        Files.writeString(existing, "not an archive\n", StandardCharsets.UTF_8);
        run(dir, "existing-refused", archive, plain, List.of("-O", existing.toString()));
        run(dir, "existing-overwritten", archive, plain, List.of(
                "-O", existing.toString(), "--overwrite-output-file", "true"));

        // Extraction, which unpacks into the destination's parent.
        run(dir, "extracted", archive, plain, List.of(
                "-O", dir.resolve("unpack/copied.tar.gz").toString(),
                "--extract-after-download", "true"));

        // No destination at all, which lands in the working directory rather than beside the
        // source.
        run(dir, "default-output", archive, plain, List.of());
        // The destination is RELATIVE, so it lands in the working directory rather than beside
        // the source. Only the file itself is reported: its siblings are the harness's own.
        final Path landed = Path.of(archive.getFileName().toString()).toAbsolutePath();
        report(dir, "default-output", landed, false);
        Files.deleteIfExists(landed);

        // The refusals.
        raw(dir, "no-source", List.of("--hg38", "true"));
        raw(dir, "no-reference", List.of("--somatic", "true"));
        raw(dir, "testing-path-alone", List.of(
                "--testing-override-path-for-datasources", archive.toString()));
        raw(dir, "testing-sha-alone", List.of(
                "--testing-override-path-for-datasources-sha256", plain.toString()));
        raw(dir, "somatic-and-germline", List.of("--somatic", "true", "--germline", "true",
                "--hg38", "true"));
        raw(dir, "testing-and-somatic", List.of("--somatic", "true", "--hg38", "true",
                "--testing-override-path-for-datasources", archive.toString(),
                "--testing-override-path-for-datasources-sha256", plain.toString()));
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    /** A one-entry tar.gz, written by hand so its bytes are the fixture's own. */
    static void writeArchive(final Path path, final String name, final String content)
            throws Exception {
        final byte[] body = content.getBytes(StandardCharsets.UTF_8);
        final ByteArrayOutputStream tar = new ByteArrayOutputStream();
        final byte[] header = new byte[512];
        final byte[] nameBytes = name.getBytes(StandardCharsets.UTF_8);
        System.arraycopy(nameBytes, 0, header, 0, nameBytes.length);
        put(header, 100, "0000644");   // mode
        put(header, 108, "0000000");   // uid
        put(header, 116, "0000000");   // gid
        put(header, 124, String.format("%011o", body.length));
        put(header, 136, String.format("%011o", 0));  // mtime zero, so the bytes do not drift
        for (int i = 148; i < 156; i++) {
            header[i] = ' ';
        }
        header[156] = '0';
        put(header, 257, "ustar");
        header[263] = '0';
        header[264] = '0';
        int checksum = 0;
        for (final byte value : header) {
            checksum += value & 0xff;
        }
        put(header, 148, String.format("%06o", checksum));
        header[154] = 0;
        header[155] = ' ';
        tar.write(header);
        tar.write(body);
        tar.write(new byte[512 - (body.length % 512)]);
        tar.write(new byte[1024]);

        try (final java.util.zip.GZIPOutputStream out =
                     new java.util.zip.GZIPOutputStream(Files.newOutputStream(path))) {
            out.write(tar.toByteArray());
        }
    }

    static void put(final byte[] header, final int offset, final String text) {
        final byte[] bytes = text.getBytes(StandardCharsets.UTF_8);
        System.arraycopy(bytes, 0, header, offset, bytes.length);
    }

    static String sha256(final Path path) throws Exception {
        final MessageDigest digest = MessageDigest.getInstance("SHA-256");
        final byte[] hash = digest.digest(Files.readAllBytes(path));
        final StringBuilder text = new StringBuilder();
        for (final byte value : hash) {
            text.append(String.format("%02x", value));
        }
        return text.toString();
    }

    static void run(final Path dir, final String label, final Path archive, final Path checksum,
                    final List<String> extra) throws Exception {
        final List<String> argv = new ArrayList<>(List.of(
                "--testing-override-path-for-datasources", archive.toString(),
                "--testing-override-path-for-datasources-sha256", checksum.toString()));
        argv.addAll(extra);
        invoke(dir, label, argv);
        final int output = extra.indexOf("-O");
        if (output >= 0) {
            report(dir, label, Path.of(extra.get(output + 1)),
                    extra.contains("--extract-after-download"));
        }
    }

    static void raw(final Path dir, final String label, final List<String> argv) {
        invoke(dir, label, argv);
    }

    static void invoke(final Path dir, final String label, final List<String> argv) {
        try {
            final Object code = new FuncotatorDataSourceDownloader()
                    .instanceMain(argv.toArray(new String[0]));
            System.out.printf("run\t%s\tcode=%s%n", label, code);
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    ReferenceQueryDump.escape(masked(String.valueOf(e.getMessage()), dir)));
        }
    }

    /** What landed at the destination, and, when asked, what was unpacked beside it. */
    static void report(final Path dir, final String label, final Path destination,
                       final boolean listSiblings) throws Exception {
        if (Files.exists(destination)) {
            System.out.printf("file\t%s\t%s=%s%n", label, masked(destination.toString(), dir),
                    sha256(destination));
        }
        final Path parent = destination.getParent();
        if (listSiblings && parent != null && Files.isDirectory(parent)) {
            try (final Stream<Path> entries = Files.list(parent)) {
                final List<Path> sorted = entries.sorted(Comparator.comparing(Path::toString))
                        .toList();
                for (final Path entry : sorted) {
                    if (!entry.equals(destination)) {
                        System.out.printf("extracted\t%s\t%s%n", label,
                                masked(entry.toString(), dir));
                    }
                }
            }
        }
    }

    static String masked(final String text, final Path dir) {
        return text.replace(dir.toString(), "<dir>")
                .replace(dir.getParent().toString(), "<cwd>");
    }
}
