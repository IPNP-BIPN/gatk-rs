/*
 * The eigenvalues Spark's ARPACK path gives a panel's Gramian, to the bit, and the rank it takes.
 *
 * `RowMatrix.computeSVD(k, true, 1e-9)` hands a request to `EigenValueDecomposition.symmetricEigs`
 * when `k < n / 3` under a hundred samples (LocalARPACK, the Gramian times a vector on the driver)
 * and when `k <= n / 2` from a hundred (DistARPACK, `multiplyGramianMatrixBy` over the rows). The
 * solver is F2J's ARPACK in the oracle image, so its bits are a function of its input and of the
 * order Spark sums in, and this dump calls the same code `computeSVD` calls on panels built here:
 * the Gramian from `computeGramianMatrix`, the local product from `computeSVD`'s own lambda, the
 * distributed one from `multiplyGramianMatrixBy`, over a hundred slices at `local[1]` as GATK's
 * `SparkConverter` makes them.
 *
 * `dgetv0`'s seed is a static field that carries over from call to call in one JVM. The tool makes
 * one call per JVM, so each case here resets it to where a fresh JVM starts.
 *
 * `local[1]` is the usual summation order, not a promise (docs/spark-parallel.md): the driver
 * deserialises task results on a pool of four threads and folds them as they land, so a result that
 * lands late changes the bits. The dump runs one result getter, which is that usual order.
 *
 * Output, one row per case:
 *
 *     eigen\t<label>\tn=<n>\tk=<k>\tmode=<local|distributed>\trank=<r>\t<raw bits of each value,
 *         descending, hex, comma separated>
 *     vectors\t<label>\t<for each vector: the raw bits of its first component and of its
 *         components summed in order, hex, ':' separated; ',' between vectors>
 *     error\t<label>\t<the exception's message>
 *
 * and then the tool itself on three panels, as `PanelRankDump` writes them.
 *
 * Usage: PanelArpackDump
 */

import breeze.linalg.DenseMatrix;
import breeze.linalg.DenseVector;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.linalg.EigenValueDecomposition$;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.linalg.Vectors;
import org.apache.spark.mllib.linalg.distributed.RowMatrix;
import scala.Function1;
import scala.Tuple2;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class PanelArpackDump {

    /**
     * Interval `i` of sample `s`: `families` profiles, centred, each sample a small multiple of its
     * family's, so the rank is `families` and every value past it is rounding noise. Exact in
     * binary: the port's test rebuilds the same doubles.
     */
    static double value(final int s, final int i, final int families) {
        return (PanelRankDump.noise(s % families, i) - 704) * (1 + (s / families) % 3) / 128.0;
    }

    static double[][] rows(final int n, final int intervals, final int families) {
        final double[][] rows = new double[intervals][n];
        for (int i = 0; i < intervals; i++) {
            for (int s = 0; s < n; s++) {
                rows[i][s] = value(s, i, families);
            }
        }
        return rows;
    }

    static String hex(final double value) {
        return Long.toHexString(Double.doubleToRawLongBits(value));
    }

    @SuppressWarnings("unchecked")
    static void run(final JavaSparkContext sc, final String label, final int n, final int intervals,
                    final int families, final int k) throws Exception {
        final List<Vector> list = new ArrayList<>();
        for (final double[] row : rows(n, intervals, families)) {
            list.add(Vectors.dense(row));
        }
        final RowMatrix matrix = new RowMatrix(sc.parallelize(list, 100).rdd());
        final boolean local = n < 100;
        final Function1<DenseVector<Object>, DenseVector<Object>> mul;
        if (local) {
            final DenseMatrix<Object> gramian = (DenseMatrix<Object>) matrix.computeGramianMatrix().asBreeze();
            // computeSVD's own `v => G * v`, synthetic, so reached by name.
            final java.lang.reflect.Method times = RowMatrix.class.getDeclaredMethod(
                    "$anonfun$computeSVD$4", DenseMatrix.class, DenseVector.class);
            mul = v -> {
                try {
                    return (DenseVector<Object>) times.invoke(null, gramian, v);
                } catch (final ReflectiveOperationException e) {
                    throw new IllegalStateException(e);
                }
            };
        } else {
            mul = matrix::multiplyGramianMatrixBy;
        }
        org.netlib.arpack.Dgetv0.inits = true;
        final Tuple2<DenseVector<Object>, DenseMatrix<Object>> eigen;
        try {
            eigen = EigenValueDecomposition$.MODULE$.symmetricEigs(mul, n, k, 1e-10, Math.max(300, 3 * k));
        } catch (final IllegalStateException e) {
            System.out.printf("error\t%s\t%s%n", label, e.getMessage());
            return;
        }
        final double[] values = (double[]) eigen._1().data();
        final double[] vectors = (double[]) eigen._2().data();
        // `computeSVD`'s rank: square roots, kept while at least rCond times the first.
        final double threshold = 1e-9 * Math.sqrt(values[0]);
        int rank = 0;
        while (rank < Math.min(k, values.length) && Math.sqrt(values[rank]) >= threshold) {
            rank++;
        }
        final StringBuilder bits = new StringBuilder();
        final StringBuilder prints = new StringBuilder();
        for (int c = 0; c < values.length; c++) {
            bits.append(c == 0 ? "" : ",").append(hex(values[c]));
            double sum = 0;
            for (int r = 0; r < n; r++) {
                sum += vectors[r + c * n];
            }
            prints.append(c == 0 ? "" : ",").append(hex(vectors[c * n])).append(':').append(hex(sum));
        }
        System.out.printf("eigen\t%s\tn=%d\tk=%d\tmode=%s\trank=%d\t%s%n", label, n, k,
                local ? "local" : "distributed", rank, bits);
        System.out.printf("vectors\t%s\t%s%n", label, prints);
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# PanelArpackDump: symmetricEigs on panels' Gramians, at local[1]");
        final JavaSparkContext sc = new JavaSparkContext(
                new SparkConf().setMaster("local[1]").setAppName("PanelArpackDump").set("spark.ui.enabled", "false")
                        // One result getter: the driver folds the task results in task order.
                        .set("spark.resultGetter.threads", "1"));
        try {
            // LocalARPACK: 20 < n / 3 from 63 samples.
            run(sc, "local-63-full", 63, 300, 63, 20);
            run(sc, "local-70-of-40", 70, 300, 40, 20);
            run(sc, "local-81-of-25", 81, 250, 25, 20);
            run(sc, "local-99-full-wide", 99, 1000, 99, 20);
            // Fewer intervals than slices: some slices are empty, and the rank is the intervals'.
            run(sc, "local-90-short", 90, 60, 90, 20);
            // Fewer families than the request: the values past them are noise ARPACK has to reach.
            run(sc, "local-70-of-12", 70, 300, 12, 20);
            // DistARPACK: from a hundred samples.
            run(sc, "distributed-100-full", 100, 400, 100, 20);
            run(sc, "distributed-120-of-60", 120, 400, 60, 20);
            run(sc, "distributed-150-of-30", 150, 333, 30, 20);
            run(sc, "distributed-110-short", 110, 90, 110, 20);
            run(sc, "distributed-120-of-15", 120, 300, 15, 20);
        } finally {
            sc.stop();
        }

        // The tool end to end at local[1], from counts through its own standardisation.
        final Path dir = Path.of("panel-arpack-dump").toAbsolutePath();
        Files.createDirectories(dir);
        PanelRankDump.run(dir, "tool-independent-70", PanelRankDump.independent(70, 120), PanelRankDump.NO_FILTERS);
        PanelRankDump.run(dir, "tool-families-70-of-30", PanelRankDump.families(70, 30, 120), PanelRankDump.NO_FILTERS);
        PanelRankDump.run(dir, "tool-independent-120", PanelRankDump.independent(120, 150), PanelRankDump.NO_FILTERS);
    }
}
