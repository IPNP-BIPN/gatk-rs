/*
 * F2J's `dgesdd('A', n, n, ...)` on symmetric matrices, the call Breeze makes for Spark's LocalLAPACK
 * `computeSVD`: what `crates/gatk-tools/src/f2j_gesdd.rs` ports, bit for bit.
 *
 * Breeze passes `dev.ludovic.netlib.lapack.F2jLAPACK`, which forwards to `org.netlib.lapack.Dgesdd`
 * (no native LAPACK loads in this image), with the workspace Breeze sizes:
 * `3 min^2 + max(max(m, n), 4 min^2 + 4 min)`. This calls it the same way. The matrices are Gramians
 * built from integers with only `+`, `*` and `/`, so the port rebuilds the same bits; past 25 the
 * bidiagonal SVD is `dbdsdc`'s divide and conquer, past 32 the left vectors go through the blocked
 * `dormqr`.
 *
 * Output, one row per matrix:
 *
 *     svd\t<kind>\t<n>\tinfo=<i>\trank=<r>\ts=<hash>\tu=<hash>\tfirst=<hex>\tlast=<hex>
 *
 * `rank` is Spark's: the square roots of the values kept while at least 1e-9 times the first. The
 * hashes fold `doubleToRawLongBits` of the values, then of U column major, as `h = (h ^ bits) * P`.
 *
 * Usage: DgesddDump. With `-Dfull=true` each row is followed by every value's bits, for a local
 * look at where a port parts from it; the golden never carries those.
 */

import org.netlib.lapack.Dgesdd;
import org.netlib.util.intW;

public class DgesddDump {

    static final long PRIME = 0x100000001b3L;

    /** `PanelRankDump.noise`. */
    static long noise(final int sample, final int interval) {
        final long i = interval + 7L;
        return 200 + ((sample + 1L) * i * i * 2654435761L + sample * 97L) % 1009;
    }

    /** `rows` by `n`, row major: what a panel's samples look like as columns. */
    static double[][] panel(final String kind, final int n) {
        final int rows = 2 * n + 3;
        final int families = Math.max(2, n / 5);
        final double[][] x = new double[rows][n];
        for (int r = 0; r < rows; r++) {
            for (int s = 0; s < n; s++) {
                switch (kind) {
                    case "families":
                        x[r][s] = noise(s % families, r) / 1009.0 * (1 + s / families);
                        break;
                    case "flat":
                        x[r][s] = s == 0 ? 0.0 : noise(s == n - 1 ? 1 : s, r) / 1009.0;
                        break;
                    default:
                        x[r][s] = noise(s, r) / 1009.0;
                }
            }
            if (kind.equals("centered") || kind.equals("families")) {
                double sum = 0.0;
                for (int s = 0; s < n; s++) {
                    sum += x[r][s];
                }
                final double mean = sum / n;
                for (int s = 0; s < n; s++) {
                    x[r][s] -= mean;
                }
            }
        }
        return x;
    }

    /** The Gramian, column major, each entry summed over the rows in order. */
    static double[] gramian(final double[][] x, final int n) {
        final double[] g = new double[n * n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (final double[] row : x) {
                    sum += row[i] * row[j];
                }
                g[i + j * n] = sum;
            }
        }
        return g;
    }

    static long fold(long h, final double[] values) {
        for (final double value : values) {
            h = (h ^ Double.doubleToRawLongBits(value)) * PRIME;
        }
        return h;
    }

    public static void main(final String[] args) {
        System.out.println("# DgesddDump: F2J dgesdd('A') on Gramians, as Breeze calls it");
        final int[] large = {62, 64, 65, 70, 75, 77, 80, 85, 90, 95, 99, 100, 104, 113, 127, 128};
        final int[] sizes = new int[59 + large.length];
        for (int n = 2; n <= 60; n++) {
            sizes[n - 2] = n;
        }
        System.arraycopy(large, 0, sizes, 59, large.length);
        for (final String kind : new String[]{"independent", "centered", "families", "flat"}) {
            for (final int n : sizes) {
                final double[] a = gramian(panel(kind, n), n);
                final double[] s = new double[n];
                final double[] u = new double[n * n];
                final double[] vt = new double[n * n];
                final int lwork = (int) (3L * n * n + Math.max(n, 4L * n * n + 4L * n));
                final double[] work = new double[lwork];
                final int[] iwork = new int[8 * n];
                final intW info = new intW(0);
                Dgesdd.dgesdd("A", n, n, a, 0, n, s, 0, u, 0, n, vt, 0, n, work, 0, lwork, iwork, 0, info);
                // RowMatrix.computeSVD: sigmas = sqrt(values), kept while >= rCond * sigmas(0).
                final double threshold = 1e-9 * Math.sqrt(s[0]);
                int rank = 0;
                while (rank < n && Math.sqrt(s[rank]) >= threshold) {
                    rank++;
                }
                System.out.printf("svd\t%s\t%d\tinfo=%d\trank=%d\ts=%016x\tu=%016x\tfirst=%016x\tlast=%016x%n",
                        kind, n, info.val, rank, fold(0, s), fold(0, u),
                        Double.doubleToRawLongBits(s[0]), Double.doubleToRawLongBits(s[n - 1]));
                if (Boolean.getBoolean("full")) {
                    final StringBuilder bits = new StringBuilder();
                    for (final double value : s) {
                        bits.append(String.format(" %016x", Double.doubleToRawLongBits(value)));
                    }
                    System.out.printf("values\t%s\t%d\t%s%n", kind, n, bits.toString().trim());
                }
            }
        }
    }
}
