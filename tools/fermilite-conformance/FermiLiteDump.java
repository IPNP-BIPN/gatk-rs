import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.EOFException;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.channels.Channels;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;

/**
 * The real library's answers, in the format gatk-fermilite's {@code write_dump} writes.
 *
 * <pre>java -cp .:$GATK_JAR FermiLiteDump corpus.bin dump.txt [threads]</pre>
 *
 * An assert in libfml takes the JVM with it, which is an answer: this process relaunches itself as a
 * worker, and when a worker dies, writes the line for the assembly it died on and starts another past
 * it. The worker writes the record it starts before starting it. The worker's stderr goes to
 * {@code dump.txt.stderr}, where the assertion's message is.
 */
public final class FermiLiteDump {
    private static final class Read implements FermiLiteAssembler.BasesAndQuals {
        final byte[] b, q;
        Read(final byte[] b, final byte[] q) { this.b = b; this.q = q; }
        @Override public byte[] getBases() { return b; }
        @Override public byte[] getQuals() { return q; }
    }

    public static void main(final String[] a) throws Exception {
        if (a.length >= 1 && a[0].equals("--worker")) {
            try {
                worker(a[1], a[2], a[3], Long.parseLong(a[4]), Long.parseLong(a[5]), Integer.parseInt(a[6]));
            } catch (final Throwable t) {
                t.printStackTrace();
                System.exit(3);
            }
        } else {
            parent(a[0], a[1], a.length > 2 ? a[2] : "0");
        }
    }

    private static void parent(final String corpus, final String dump, final String threads) throws Exception {
        final String progress = dump + ".progress";
        new FileOutputStream(dump).close();
        long idx = 0, offset = 6;
        int deaths = 0;
        for (;;) {
            new File(progress).delete();
            final List<String> cmd = new ArrayList<>();
            cmd.add(System.getProperty("java.home") + "/bin/java");
            cmd.add("-XX:-CreateCoredumpOnCrash");
            cmd.add("-XX:TieredStopAtLevel=1");
            cmd.add("-Xshare:auto");
            cmd.add("-XX:ErrorFile=/dev/null");
            cmd.add("-cp");
            cmd.add(System.getProperty("java.class.path"));
            cmd.add("FermiLiteDump");
            cmd.add("--worker");
            cmd.add(corpus);
            cmd.add(dump);
            cmd.add(progress);
            cmd.add(Long.toString(idx));
            cmd.add(Long.toString(offset));
            cmd.add(threads);
            // HotSpot's own crash banner goes to stdout: keep it, because the difference between an
            // assert (SIGABRT, no banner) and anything else must stay visible in dump.txt.stderr
            final Process p = new ProcessBuilder(cmd)
                    .redirectErrorStream(true)
                    .redirectOutput(ProcessBuilder.Redirect.appendTo(new File(dump + ".stderr")))
                    .start();
            final int rc = p.waitFor();
            if (rc == 0) break;
            // 134 is SIGABRT from an assert or from malloc's checks. HotSpot's own crash handling, for
            // a SIGSEGV in native code, leaves with 1 when core dumps are off. Java exceptions are 3.
            if (rc != 134 && rc != 1) throw new IllegalStateException("worker exited with " + rc);
            final String[] f;
            try (BufferedReader r = new BufferedReader(new FileReader(progress))) { f = r.readLine().split(" "); }
            idx = Long.parseLong(f[0]) + 1;
            offset = Long.parseLong(f[1]);
            final String line = (rc == 134 ? "X " : "S ") + f[0] + " " + f[2] + " flag=" + f[3] + " reads=" + f[4]
                    + (rc == 134 ? " ABORT\n" : " SEGV\n");
            Files.write(Paths.get(dump), line.getBytes(StandardCharsets.ISO_8859_1), StandardOpenOption.APPEND);
            deaths++;
        }
        System.err.println("assemblies through " + idx + "; JVMs lost: " + deaths);
    }

    private static void worker(final String corpus, final String dump, final String progress,
                               long idx, final long offset, final int threads) throws Exception {
        final RandomAccessFile raf = new RandomAccessFile(corpus, "r");
        raf.seek(offset);
        final CountingIn cin = new CountingIn(new BufferedInputStream(Channels.newInputStream(raf.getChannel()), 1 << 20), offset);
        final DataInputStream in = new DataInputStream(cin);
        final OutputStream out = new BufferedOutputStream(new FileOutputStream(dump, true), 1 << 20);
        final RandomAccessFile prog = new RandomAccessFile(progress, "rw");
        for (;;) {
            final int nameLen;
            try { nameLen = in.readInt(); } catch (final EOFException e) { break; }
            final byte[] nb = new byte[nameLen];
            in.readFully(nb);
            final String name = new String(nb, StandardCharsets.US_ASCII);
            final int flag = in.readInt();
            final int n = in.readInt();
            final List<Read> reads = new ArrayList<>(n);
            for (int i = 0; i < n; i++) {
                final byte[] b = new byte[in.readInt()];
                in.readFully(b);
                final byte[] q = new byte[in.readInt()];
                in.readFully(q);
                reads.add(new Read(b, q));
            }
            final byte[] pl = (idx + " " + cin.position + " " + name + " " + flag + " " + n + "\n")
                    .getBytes(StandardCharsets.US_ASCII);
            prog.seek(0);
            prog.setLength(0);
            prog.write(pl);
            final ByteArrayOutputStream o = new ByteArrayOutputStream();
            try (FermiLiteAssembler assembler = new FermiLiteAssembler()) {
                assembler.setCleaningFlag(flag);
                if (threads > 0) assembler.setNThreads(threads);
                try {
                    write(o, idx, name, flag, n, assembler.createAssembly(reads));
                } catch (final java.nio.BufferOverflowException e) {
                    o.write(("J " + idx + " " + name + " flag=" + flag + " reads=" + n + " BufferOverflowException\n")
                            .getBytes(StandardCharsets.US_ASCII));
                }
            }
            out.write(o.toByteArray());
            out.flush();
            idx++;
        }
        out.close();
    }

    private static void write(final OutputStream o, final long idx, final String name, final int flag, final int n,
                              final FermiLiteAssembly assembly) throws Exception {
        final int m = assembly.getNContigs();
        final Map<FermiLiteAssembly.Contig, Integer> index = new IdentityHashMap<>();
        for (int i = 0; i < m; i++) index.put(assembly.getContig(i), i);
        o.write(("A " + idx + " " + name + " flag=" + flag + " reads=" + n + " contigs=" + m + "\n")
                .getBytes(StandardCharsets.US_ASCII));
        for (int i = 0; i < m; i++) {
            final FermiLiteAssembly.Contig c = assembly.getContig(i);
            o.write(("C " + i + " " + c.getSequence().length + " " + c.getNSupportingReads() + " ")
                    .getBytes(StandardCharsets.US_ASCII));
            o.write(c.getSequence());
            o.write(' ');
            o.write(c.getPerBaseCoverage());
            o.write('\n');
            for (final FermiLiteAssembly.Connection k : c.getConnections()) {
                o.write(("L " + index.get(k.getTarget()) + " " + k.getOverlapLen() + " " + (k.isRC() ? 1 : 0) + " "
                        + (k.isTargetRC() ? 1 : 0) + "\n").getBytes(StandardCharsets.US_ASCII));
            }
        }
    }

    private static final class CountingIn extends java.io.FilterInputStream {
        long position;
        CountingIn(final java.io.InputStream in, final long start) { super(in); position = start; }
        @Override public int read() throws java.io.IOException {
            final int r = super.read();
            if (r >= 0) position++;
            return r;
        }
        @Override public int read(final byte[] b, final int off, final int len) throws java.io.IOException {
            final int r = super.read(b, off, len);
            if (r > 0) position += r;
            return r;
        }
    }
}
