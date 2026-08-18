//! Ported from `org.apache.commons.math3.linear.SingularValueDecomposition` (commons-math3 3.5,
//! Apache 2.0, verified from the licence header of the sources jar in the resolved dependency).
//!
//! GATK pins commons-math **strictly** to 3.5 -- `build.gradle` says so and gives the reason
//! ("changing this breaks ModelSegmentsIntegrationTests, they're quite brittle") -- so 3.5 is the
//! version to transcribe, not the 3.6.1 that also sits in the dependency cache.
//!
//! `KernelSegmenter` is the only caller this port reaches, and through it
//! `CalculateContamination`. The decomposition is a JAMA-derived Golub-Reinsch: Householder
//! bidiagonalisation, then implicit-shift QR sweeps, all inside one constructor.
//!
//! # What byte identity means here, and where it does not
//!
//! `docs/what-the-kernel-segmenter-needs-from-the-decomposition.md` works this out. The short form:
//!
//!  * **the singular values must match to the bit.** The segmenter uses them as
//!    `1 / (sqrt(s) + 1e-10)`, which turns a rounding-level difference near zero into a factor of
//!    forty and moves the Gram matrix that decides the changepoints by ninety-six per cent;
//!  * **`U`'s null-space basis need not, and cannot.** Any orthonormal basis of the null space is
//!    as valid, and everything downstream sees `U` only through `Z Z^T`, where the freedom cancels
//!    exactly.
//!
//! This port is a transcription, so it reproduces the reference's basis as well, and the
//! conformance suite compares `U` entry by entry. The distinction above is about what a *different*
//! implementation could be held to, not about what this one delivers.
//!
//! # Why the loops are written the way they are
//!
//! Every arithmetic statement keeps the reference's order and its temporaries. The QR sweep
//! accumulates `f` and `g` across iterations and rotates `U` and `V` in place, so hoisting a
//! subexpression or reordering a rotation changes the last bits of a singular value and therefore
//! the amplification above. The array indexing is the reference's too, including the places where
//! JAMA writes past what a Rust reader would expect (`U[k][k] = 1 + U[k][k]` after the column has
//! been negated, `for i in 0..k-1` which is empty at `k = 0` rather than a wrap).
//!
//! # The one Java-specific decision
//!
//! The convergence test is written `if !(abs(e[k]) > threshold)` rather than
//! `if abs(e[k]) <= threshold`. The reference says why (MATH-947): with a NaN in the array the
//! second form never fires and the loop runs forever, while the first breaks out. Rust would behave
//! the same way, so the negation is kept rather than tidied.

/// `FastMath.getExponent`: the raw biased exponent, less the bias.
///
/// Raw bits, so NaN and infinity both come back as 1024, which is what the reference relies on.
fn get_exponent(d: f64) -> i32 {
    ((d.to_bits() >> 52) & 0x7ff) as i32 - 1023
}

/// `FastMath.scalb`, transcribed rather than delegated.
///
/// The fast path is a multiplication by `2^n` built from bits, and the slow paths handle the
/// subnormal rounding by hand. Rust's `f64` has no `scalb`, and `d * 2f64.powi(n)` is not the same
/// function: it rounds twice once the result is subnormal.
fn scalb(d: f64, n: i32) -> f64 {
    // The common case: `2^n` is representable as a normal number, so one multiplication does it.
    if n > -1023 && n < 1024 {
        return d * f64::from_bits(((n + 1023) as u64) << 52);
    }

    if d.is_nan() || d.is_infinite() || d == 0.0 {
        return d;
    }
    if n < -2098 {
        return if d > 0.0 { 0.0 } else { -0.0 };
    }
    if n > 2097 {
        return if d > 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
    }

    let bits = d.to_bits();
    let sign = bits & 0x8000_0000_0000_0000;
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut mantissa = bits & 0x000f_ffff_ffff_ffff;
    let mut scaled_exponent = exponent + n;

    if n < 0 {
        // Really the case `n <= -1023`.
        if scaled_exponent > 0 {
            f64::from_bits(sign | ((scaled_exponent as u64) << 52) | mantissa)
        } else if scaled_exponent > -53 {
            // Normal in, subnormal out: recover the hidden bit, shift, and round up on the most
            // significant bit that falls off the end.
            mantissa |= 1 << 52;
            let most_significant_lost_bit = mantissa & (1 << (-scaled_exponent));
            mantissa >>= 1 - scaled_exponent;
            if most_significant_lost_bit != 0 {
                mantissa += 1;
            }
            f64::from_bits(sign | mantissa)
        } else {
            if sign == 0 {
                0.0
            } else {
                -0.0
            }
        }
    } else {
        // Really the case `n >= 1024`.
        if exponent == 0 {
            // Subnormal in: normalise it first, which costs one exponent per shift.
            while (mantissa >> 52) != 1 {
                mantissa <<= 1;
                scaled_exponent -= 1;
            }
            scaled_exponent += 1;
            mantissa &= 0x000f_ffff_ffff_ffff;
        }
        if scaled_exponent < 2047 {
            f64::from_bits(sign | ((scaled_exponent as u64) << 52) | mantissa)
        } else if sign == 0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        }
    }
}

/// `FastMath.hypot`, which is **not** `f64::hypot`.
///
/// The reference scales both arguments to a middle exponent, squares and adds there, and scales
/// back. The host's `hypot` is a different algorithm with a different last bit, and this one is
/// reached inside the QR sweep, where a last bit becomes a singular value.
pub fn hypot(x: f64, y: f64) -> f64 {
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let exp_x = get_exponent(x);
    let exp_y = get_exponent(y);
    // Twenty-seven binary orders apart and the smaller one cannot reach the result's last bit.
    if exp_x > exp_y + 27 {
        return x.abs();
    }
    if exp_y > exp_x + 27 {
        return y.abs();
    }
    // An intermediate scale that avoids both overflow and underflow.
    let middle_exp = (exp_x + exp_y) / 2;
    let scaled_x = scalb(x, -middle_exp);
    let scaled_y = scalb(y, -middle_exp);
    let scaled_h = (scaled_x * scaled_x + scaled_y * scaled_y).sqrt();
    scalb(scaled_h, middle_exp)
}

/// `SingularValueDecomposition.EPS`: the relative threshold for a small singular value.
const EPS: f64 = 2.220446049250313e-16; // 0x1.0p-52

/// `SingularValueDecomposition.TINY`: the absolute threshold for a small singular value.
const TINY: f64 = 1.6033346880071782e-291; // 0x1.0p-966

/// The compact decomposition of a matrix: `A = U * S * V^T`.
///
/// Fields are the reference's `getSingularValues`, `getU` and `getV`, already untransposed. Rows
/// are the outer index, as in commons-math's `double[][]`.
#[derive(Debug, Clone)]
pub struct SingularValueDecomposition {
    /// The singular values, in non-increasing order.
    pub singular_values: Vec<f64>,
    /// `U`, `m` by `n`.
    pub u: Vec<Vec<f64>>,
    /// `V`, `n` by `n`.
    pub v: Vec<Vec<f64>>,
}

impl SingularValueDecomposition {
    /// The constructor, transcribed statement by statement.
    ///
    /// `matrix` is row-major. The reference transposes when there are fewer rows than columns so
    /// that `m` is always the larger dimension, and swaps `U` and `V` back at the end.
    ///
    /// The four allowances below are all the same decision: this function's shape is the
    /// reference's shape, and a reader checking the transcription reads them side by side.
    ///
    ///  * `needless_range_loop` -- the indices are the reference's, and several loops read one
    ///    array while writing another at the same index, or run over a sub-range chosen by the
    ///    algorithm rather than by the array's length;
    ///  * `assign_op_pattern` -- `U[k][k] = 1 + U[k][k]` is written that way in JAMA;
    ///  * `manual_swap` -- the ordering step swaps through a named temporary that the reference
    ///    reuses for the `U` and `V` swaps that follow;
    ///  * `neg_cmp_op_on_partial_ord` -- `!(abs(e[k]) > threshold)` is deliberate and documented at
    ///    the module level: MATH-947, the direct form never fires on a NaN and the loop never ends.
    #[allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::manual_swap,
        clippy::neg_cmp_op_on_partial_ord
    )]
    pub fn new(matrix: &[Vec<f64>]) -> Self {
        let rows = matrix.len();
        let columns = if rows == 0 { 0 } else { matrix[0].len() };

        // "m" is always the largest dimension.
        let (transposed, mut a, m, n) = if rows < columns {
            let mut transpose = vec![vec![0.0; rows]; columns];
            for (i, row) in matrix.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    transpose[j][i] = *value;
                }
            }
            (true, transpose, columns, rows)
        } else {
            (false, matrix.to_vec(), rows, columns)
        };

        let mut singular_values = vec![0.0; n];
        let mut u = vec![vec![0.0; n]; m];
        let mut v = vec![vec![0.0; n]; n];
        let mut e = vec![0.0; n];
        let mut work = vec![0.0; m];

        // Reduce A to bidiagonal form, the diagonal in `singular_values` and the super-diagonal
        // in `e`.
        let nct = (m as i64 - 1).min(n as i64);
        let nrt = 0.max(n as i64 - 2);
        let mut k: i64 = 0;
        while k < nct.max(nrt) {
            let ku = k as usize;
            if k < nct {
                // The transformation for the k-th column, whose 2-norm goes in `s[k]`. `hypot`
                // rather than a sum of squares, so that neither end of the range overflows.
                singular_values[ku] = 0.0;
                for i in ku..m {
                    singular_values[ku] = hypot(singular_values[ku], a[i][ku]);
                }
                if singular_values[ku] != 0.0 {
                    if a[ku][ku] < 0.0 {
                        singular_values[ku] = -singular_values[ku];
                    }
                    for i in ku..m {
                        a[i][ku] /= singular_values[ku];
                    }
                    a[ku][ku] += 1.0;
                }
                singular_values[ku] = -singular_values[ku];
            }
            for j in (ku + 1)..n {
                if k < nct && singular_values[ku] != 0.0 {
                    // Apply the transformation.
                    let mut t = 0.0;
                    for i in ku..m {
                        t += a[i][ku] * a[i][j];
                    }
                    t = -t / a[ku][ku];
                    for i in ku..m {
                        a[i][j] += t * a[i][ku];
                    }
                }
                // The k-th row of A, kept for the row transformation below.
                e[j] = a[ku][j];
            }
            if k < nct {
                // Keep the transformation in U for the back multiplication.
                for i in ku..m {
                    u[i][ku] = a[i][ku];
                }
            }
            if k < nrt {
                // The k-th row transformation, its 2-norm in `e[k]`.
                e[ku] = 0.0;
                for i in (ku + 1)..n {
                    e[ku] = hypot(e[ku], e[i]);
                }
                if e[ku] != 0.0 {
                    if e[ku + 1] < 0.0 {
                        e[ku] = -e[ku];
                    }
                    for i in (ku + 1)..n {
                        e[i] /= e[ku];
                    }
                    e[ku + 1] += 1.0;
                }
                e[ku] = -e[ku];
                if ku + 1 < m && e[ku] != 0.0 {
                    // Apply the transformation.
                    for item in work.iter_mut().take(m).skip(ku + 1) {
                        *item = 0.0;
                    }
                    for j in (ku + 1)..n {
                        for i in (ku + 1)..m {
                            work[i] += e[j] * a[i][j];
                        }
                    }
                    for j in (ku + 1)..n {
                        let t = -e[j] / e[ku + 1];
                        for i in (ku + 1)..m {
                            a[i][j] += t * work[i];
                        }
                    }
                }
                // Keep the transformation in V for the back multiplication.
                for i in (ku + 1)..n {
                    v[i][ku] = e[i];
                }
            }
            k += 1;
        }

        // The final bidiagonal matrix, of order p.
        let mut p = n;
        if (nct as usize) < n {
            singular_values[nct as usize] = a[nct as usize][nct as usize];
        }
        if m < p {
            singular_values[p - 1] = 0.0;
        }
        if (nrt as usize) + 1 < p {
            e[nrt as usize] = a[nrt as usize][p - 1];
        }
        e[p - 1] = 0.0;

        // Generate U.
        for j in (nct.max(0) as usize)..n {
            for row in u.iter_mut().take(m) {
                row[j] = 0.0;
            }
            u[j][j] = 1.0;
        }
        let mut k = nct - 1;
        while k >= 0 {
            let ku = k as usize;
            if singular_values[ku] != 0.0 {
                for j in (ku + 1)..n {
                    let mut t = 0.0;
                    for i in ku..m {
                        t += u[i][ku] * u[i][j];
                    }
                    t = -t / u[ku][ku];
                    for i in ku..m {
                        u[i][j] += t * u[i][ku];
                    }
                }
                for i in ku..m {
                    u[i][ku] = -u[i][ku];
                }
                u[ku][ku] = 1.0 + u[ku][ku];
                // JAMA stops one short of the diagonal here, so at `k = 0` this clears nothing.
                for i in 0..ku.saturating_sub(1) {
                    u[i][ku] = 0.0;
                }
            } else {
                for row in u.iter_mut().take(m) {
                    row[ku] = 0.0;
                }
                u[ku][ku] = 1.0;
            }
            k -= 1;
        }

        // Generate V.
        let mut k = n as i64 - 1;
        while k >= 0 {
            let ku = k as usize;
            if k < nrt && e[ku] != 0.0 {
                for j in (ku + 1)..n {
                    let mut t = 0.0;
                    for i in (ku + 1)..n {
                        t += v[i][ku] * v[i][j];
                    }
                    t = -t / v[ku + 1][ku];
                    for i in (ku + 1)..n {
                        v[i][j] += t * v[i][ku];
                    }
                }
            }
            for row in v.iter_mut().take(n) {
                row[ku] = 0.0;
            }
            v[ku][ku] = 1.0;
            k -= 1;
        }

        // The main iteration for the singular values.
        // `pp` is fixed before the loop and does not follow `p` down: the reference uses it for the
        // sign flip and the ordering swap, both of which address the original last column.
        let pp = p as i64 - 1;
        while p > 0 {
            let mut k: i64;
            let kase: i32;
            // The reference's case analysis, verbatim:
            //   1  s(p) and e[k-1] negligible and k < p
            //   2  s(k) negligible and k < p
            //   3  e[k-1] negligible, k < p, and s(k)..s(p) not negligible: a QR step
            //   4  e(p-1) negligible: convergence
            k = p as i64 - 2;
            while k >= 0 {
                let ku = k as usize;
                let threshold =
                    TINY + EPS * (singular_values[ku].abs() + singular_values[ku + 1].abs());
                // Written as a negation on purpose (MATH-947): with a NaN in `e` the direct form
                // never fires and the loop never ends.
                if !(e[ku].abs() > threshold) {
                    e[ku] = 0.0;
                    break;
                }
                k -= 1;
            }

            if k == p as i64 - 2 {
                kase = 4;
            } else {
                let mut ks = p as i64 - 1;
                while ks >= k {
                    if ks == k {
                        break;
                    }
                    let t = (if ks != p as i64 {
                        e[ks as usize].abs()
                    } else {
                        0.0
                    }) + (if ks != k + 1 {
                        e[(ks - 1) as usize].abs()
                    } else {
                        0.0
                    });
                    if singular_values[ks as usize].abs() <= TINY + EPS * t {
                        singular_values[ks as usize] = 0.0;
                        break;
                    }
                    ks -= 1;
                }
                if ks == k {
                    kase = 3;
                } else if ks == p as i64 - 1 {
                    kase = 1;
                } else {
                    kase = 2;
                    k = ks;
                }
            }
            k += 1;
            let ku = k as usize;

            match kase {
                // Deflate a negligible s(p).
                1 => {
                    let mut f = e[p - 2];
                    e[p - 2] = 0.0;
                    let mut j = p as i64 - 2;
                    while j >= k {
                        let ju = j as usize;
                        let mut t = hypot(singular_values[ju], f);
                        let cs = singular_values[ju] / t;
                        let sn = f / t;
                        singular_values[ju] = t;
                        if j != k {
                            f = -sn * e[ju - 1];
                            e[ju - 1] = cs * e[ju - 1];
                        }
                        for row in v.iter_mut().take(n) {
                            t = cs * row[ju] + sn * row[p - 1];
                            row[p - 1] = -sn * row[ju] + cs * row[p - 1];
                            row[ju] = t;
                        }
                        j -= 1;
                    }
                }
                // Split at a negligible s(k).
                2 => {
                    let mut f = e[ku - 1];
                    e[ku - 1] = 0.0;
                    for j in ku..p {
                        let mut t = hypot(singular_values[j], f);
                        let cs = singular_values[j] / t;
                        let sn = f / t;
                        singular_values[j] = t;
                        f = -sn * e[j];
                        e[j] = cs * e[j];
                        for row in u.iter_mut().take(m) {
                            t = cs * row[j] + sn * row[ku - 1];
                            row[ku - 1] = -sn * row[j] + cs * row[ku - 1];
                            row[j] = t;
                        }
                    }
                }
                // One QR step.
                3 => {
                    // The shift, computed on values scaled by the largest of the five so that the
                    // squares below cannot overflow.
                    let max_pm1_pm2 = singular_values[p - 1]
                        .abs()
                        .max(singular_values[p - 2].abs());
                    let scale = max_pm1_pm2
                        .max(e[p - 2].abs())
                        .max(singular_values[ku].abs())
                        .max(e[ku].abs());
                    let sp = singular_values[p - 1] / scale;
                    let spm1 = singular_values[p - 2] / scale;
                    let epm1 = e[p - 2] / scale;
                    let sk = singular_values[ku] / scale;
                    let ek = e[ku] / scale;
                    let b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                    let c = (sp * epm1) * (sp * epm1);
                    let mut shift = 0.0;
                    if b != 0.0 || c != 0.0 {
                        shift = (b * b + c).sqrt();
                        if b < 0.0 {
                            shift = -shift;
                        }
                        shift = c / (b + shift);
                    }
                    let mut f = (sk + sp) * (sk - sp) + shift;
                    let mut g = sk * ek;
                    // Chase the zeros down the bidiagonal.
                    for j in ku..(p - 1) {
                        let mut t = hypot(f, g);
                        let mut cs = f / t;
                        let mut sn = g / t;
                        if j as i64 != k {
                            e[j - 1] = t;
                        }
                        f = cs * singular_values[j] + sn * e[j];
                        e[j] = cs * e[j] - sn * singular_values[j];
                        g = sn * singular_values[j + 1];
                        singular_values[j + 1] = cs * singular_values[j + 1];
                        for row in v.iter_mut().take(n) {
                            t = cs * row[j] + sn * row[j + 1];
                            row[j + 1] = -sn * row[j] + cs * row[j + 1];
                            row[j] = t;
                        }
                        t = hypot(f, g);
                        cs = f / t;
                        sn = g / t;
                        singular_values[j] = t;
                        f = cs * e[j] + sn * singular_values[j + 1];
                        singular_values[j + 1] = -sn * e[j] + cs * singular_values[j + 1];
                        g = sn * e[j + 1];
                        e[j + 1] = cs * e[j + 1];
                        if j < m - 1 {
                            for row in u.iter_mut().take(m) {
                                t = cs * row[j] + sn * row[j + 1];
                                row[j + 1] = -sn * row[j] + cs * row[j + 1];
                                row[j] = t;
                            }
                        }
                    }
                    e[p - 2] = f;
                }
                // Convergence.
                _ => {
                    // Make the singular value positive. A negative zero becomes a positive zero
                    // here, and the column of V is only flipped when the value was truly negative.
                    if singular_values[ku] <= 0.0 {
                        singular_values[ku] = if singular_values[ku] < 0.0 {
                            -singular_values[ku]
                        } else {
                            0.0
                        };
                        for row in v.iter_mut().take((pp + 1).max(0) as usize) {
                            row[ku] = -row[ku];
                        }
                    }
                    // Order the singular values, one adjacent swap at a time.
                    let mut k = k;
                    while k < pp {
                        let ku = k as usize;
                        if singular_values[ku] >= singular_values[ku + 1] {
                            break;
                        }
                        let mut t = singular_values[ku];
                        singular_values[ku] = singular_values[ku + 1];
                        singular_values[ku + 1] = t;
                        if ku < n - 1 {
                            for row in v.iter_mut().take(n) {
                                t = row[ku + 1];
                                row[ku + 1] = row[ku];
                                row[ku] = t;
                            }
                        }
                        if ku < m - 1 {
                            for row in u.iter_mut().take(m) {
                                t = row[ku + 1];
                                row[ku + 1] = row[ku];
                                row[ku] = t;
                            }
                        }
                        k += 1;
                    }
                    p -= 1;
                }
            }
        }
        if transposed {
            SingularValueDecomposition {
                singular_values,
                u: v,
                v: u,
            }
        } else {
            SingularValueDecomposition {
                singular_values,
                u,
                v,
            }
        }
    }
}
