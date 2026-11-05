//! `dgesdd('A', n, n, ...)` as the F2J translation of LAPACK 3.1.1 runs it, for `n` up to 128.
//!
//! Spark's `RowMatrix.computeSVD` answers a panel of normals from `breeze.linalg.svd` over the
//! Gramian, and in the oracle image that is `dev.ludovic.netlib.lapack.F2jLAPACK`: no native LAPACK
//! loads, so the arithmetic is `org.netlib.lapack.Dgesdd` in a strict JVM, and its bits are a
//! function of its input. The Gramian of a standardised panel is rank deficient, so its trailing
//! values are rounding noise, and the noise is what `rCond` compares: only the same operations in
//! the same order reach the same rank.
//!
//! The path is the one a square matrix takes: `dgebd2` (`dgebrd` leaves its blocked code up to
//! `NX`, 128), then `dbdsdc` ([`crate::f2j_dbdsdc`]), then `dormqr` for the left vectors, which is
//! `dorm2r` up to `NB` (32) reflectors and blocks of 32 past it. The vectors on the right are not
//! accumulated: they feed nothing back into the values, and Spark throws them away.
//!
//! Ported from LAPACK 3.1.1 (`SRC/dgesdd.f`, `dgebd2.f`, `dlarfg.f`, `dlarf.f`, `dlartg.f`,
//! `dlas2.f`, `dlasv2.f`, `dlapy2.f`, `dormqr.f`, `dorm2r.f`, `dlarft.f`, `dlarfb.f`) and reference
//! BLAS (`dnrm2.f`, `dgemv.f`, `dger.f`, `dtrmv.f`, `dtrmm.f`, `dgemm.f`), BSD 3-Clause, with what
//! F2J does to them: `x**2` is `Math.pow(x, 2)`, which HotSpot answers as `x * x`, and `MAX` and
//! `MIN` are `Math.max` and `Math.min`.

/// The largest matrix `dgebrd` leaves to `dgebd2`: `NX`, `ilaenv(3, 'DGEBRD', ...)`.
pub const LARGEST: usize = 128;

/// `dlamch('E')`: F2J's is the rounding unit.
pub(crate) const EPS: f64 = f64::EPSILON / 2.0;
/// `dlamch('P')`, the epsilon times the base.
const PRECISION: f64 = f64::EPSILON;
/// `dlamch('S')`.
pub(crate) const SAFE_MINIMUM: f64 = f64::MIN_POSITIVE;
/// `Math.pow(EPS, -0.125)`, measured on the oracle: 2^(53/8), under `dbdsqr`'s cap of a hundred.
pub(crate) const TOLERANCE_MULTIPLIER: f64 = f64::from_bits(0x4058_ACE5_422A_A0DB);
/// `SAFMN2`, `2**INT(LOG(SAFMIN / EPS) / LOG(2) / 2)`: the logarithms give -484.5, and INT truncates.
const SAFE_MINIMUM_2: f64 = f64::from_bits((1023 - 484) << 52);

/// A square matrix, column major, read as Fortran reads it: one-based.
struct Square<'a> {
    n: usize,
    values: &'a mut [f64],
}

impl Square<'_> {
    fn at(&self, row: usize, column: usize) -> usize {
        (row - 1) + (column - 1) * self.n
    }
}

/// `SIGN(a, b)`. F2J's `Util.dsign` keeps the magnitude and takes `b >= 0`, so a negative zero
/// counts as positive.
pub(crate) fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// Fortran's `MAX`, which F2J writes as `Math.max`: a negative zero is under a positive one.
pub(crate) fn max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        a
    } else if a == 0.0 && b == 0.0 && a.is_sign_negative() {
        b
    } else if a >= b {
        a
    } else {
        b
    }
}

/// Fortran's `MIN`, `Math.min`.
pub(crate) fn min(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        a
    } else if a == 0.0 && b == 0.0 && b.is_sign_negative() {
        b
    } else if a <= b {
        a
    } else {
        b
    }
}

/// `MAX` of three, F2J's `Util.max`: `Math.max(a > b ? a : b, Math.max(b, c))`.
pub(crate) fn max3(a: f64, b: f64, c: f64) -> f64 {
    max(if a > b { a } else { b }, max(b, c))
}

/// `dnrm2`, the scaled sum of squares rather than a square root of a dot product.
pub(crate) fn dnrm2(x: &[f64]) -> f64 {
    match x.len() {
        0 => 0.0,
        1 => x[0].abs(),
        _ => {
            let mut scale = 0.0f64;
            let mut ssq = 1.0f64;
            for value in x {
                if *value != 0.0 {
                    let absxi = value.abs();
                    if scale < absxi {
                        let ratio = scale / absxi;
                        ssq = 1.0 + ssq * (ratio * ratio);
                        scale = absxi;
                    } else {
                        let ratio = absxi / scale;
                        ssq += ratio * ratio;
                    }
                }
            }
            scale * ssq.sqrt()
        }
    }
}

/// `dlapy2`.
pub(crate) fn dlapy2(x: f64, y: f64) -> f64 {
    let (xabs, yabs) = (x.abs(), y.abs());
    let w = max(xabs, yabs);
    let z = min(xabs, yabs);
    if z == 0.0 {
        w
    } else {
        let ratio = z / w;
        w * (1.0 + ratio * ratio).sqrt()
    }
}

/// `dlarfg`: `alpha` and `x` become beta and the reflector's tail. Returns tau.
fn dlarfg(alpha: &mut f64, x: &mut [f64]) -> f64 {
    let mut xnorm = dnrm2(x);
    if xnorm == 0.0 {
        return 0.0;
    }
    let mut beta = -sign(dlapy2(*alpha, xnorm), *alpha);
    let safmin = SAFE_MINIMUM / EPS;
    if beta.abs() < safmin {
        let rsafmn = 1.0 / safmin;
        let mut knt = 0;
        loop {
            knt += 1;
            for value in x.iter_mut() {
                *value *= rsafmn;
            }
            beta *= rsafmn;
            *alpha *= rsafmn;
            if beta.abs() >= safmin {
                break;
            }
        }
        xnorm = dnrm2(x);
        beta = -sign(dlapy2(*alpha, xnorm), *alpha);
        let tau = (beta - *alpha) / beta;
        let factor = 1.0 / (*alpha - beta);
        for value in x.iter_mut() {
            *value *= factor;
        }
        *alpha = beta;
        for _ in 0..knt {
            *alpha *= safmin;
        }
        tau
    } else {
        let tau = (beta - *alpha) / beta;
        let factor = 1.0 / (*alpha - beta);
        for value in x.iter_mut() {
            *value *= factor;
        }
        *alpha = beta;
        tau
    }
}

/// `dlarf('L', ...)` on the block of `c` from `(row, column)`: `dgemv('T')` then `dger`.
fn dlarf_left(v: &[f64], tau: f64, c: &mut Square, row: usize, column: usize) {
    if tau == 0.0 {
        return;
    }
    let (m, columns) = (v.len(), c.n + 1 - column);
    let mut work = vec![0.0; columns];
    for (j, slot) in work.iter_mut().enumerate() {
        let mut temp = 0.0;
        for (i, element) in v.iter().enumerate().take(m) {
            temp += c.values[c.at(row + i, column + j)] * element;
        }
        *slot += temp;
    }
    let alpha = -tau;
    for (j, slot) in work.iter().enumerate() {
        if *slot != 0.0 {
            let temp = alpha * slot;
            for (i, element) in v.iter().enumerate().take(m) {
                let at = c.at(row + i, column + j);
                c.values[at] += element * temp;
            }
        }
    }
}

/// `dlarf('R', ...)` on the block of `c` from `(row, column)`: `dgemv('N')` then `dger`.
fn dlarf_right(v: &[f64], tau: f64, c: &mut Square, row: usize, column: usize) {
    if tau == 0.0 || row > c.n {
        return;
    }
    let rows = c.n + 1 - row;
    let mut work = vec![0.0; rows];
    for (j, element) in v.iter().enumerate() {
        if *element != 0.0 {
            let temp = *element;
            for (i, slot) in work.iter_mut().enumerate() {
                *slot += temp * c.values[c.at(row + i, column + j)];
            }
        }
    }
    let alpha = -tau;
    for (j, element) in v.iter().enumerate() {
        if *element != 0.0 {
            let temp = alpha * element;
            for (i, slot) in work.iter().enumerate() {
                let at = c.at(row + i, column + j);
                c.values[at] += slot * temp;
            }
        }
    }
}

/// `dgebd2`, for a square matrix: the diagonal, the superdiagonal and the left taus. The left
/// reflectors stay below the diagonal of `a`.
fn dgebd2(a: &mut Square) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = a.n;
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n.saturating_sub(1)];
    let mut tauq = vec![0.0; n];
    for i in 1..=n {
        let mut alpha = a.values[a.at(i, i)];
        let mut tail: Vec<f64> = (i + 1..=n).map(|r| a.values[a.at(r, i)]).collect();
        tauq[i - 1] = dlarfg(&mut alpha, &mut tail);
        d[i - 1] = alpha;
        for (offset, value) in tail.iter().enumerate() {
            let at = a.at(i + 1 + offset, i);
            a.values[at] = *value;
        }
        if i < n {
            let mut v = vec![1.0];
            v.extend_from_slice(&tail);
            dlarf_left(&v, tauq[i - 1], a, i, i + 1);
            let at = a.at(i, i);
            a.values[at] = d[i - 1];

            let mut alpha = a.values[a.at(i, i + 1)];
            let mut tail: Vec<f64> = (i + 2..=n).map(|c| a.values[a.at(i, c)]).collect();
            let taup = dlarfg(&mut alpha, &mut tail);
            e[i - 1] = alpha;
            for (offset, value) in tail.iter().enumerate() {
                let at = a.at(i, i + 2 + offset);
                a.values[at] = *value;
            }
            let mut v = vec![1.0];
            v.extend_from_slice(&tail);
            dlarf_right(&v, taup, a, i + 1, i + 1);
            let at = a.at(i, i + 1);
            a.values[at] = e[i - 1];
        } else {
            let at = a.at(i, i);
            a.values[at] = d[i - 1];
        }
    }
    (d, e, tauq)
}

/// `dlartg`: `(cs, sn, r)`.
pub(crate) fn dlartg(f: f64, g: f64) -> (f64, f64, f64) {
    if g == 0.0 {
        return (1.0, 0.0, f);
    }
    if f == 0.0 {
        return (0.0, 1.0, g);
    }
    let safmn2 = SAFE_MINIMUM_2;
    let safmx2 = 1.0 / safmn2;
    let (mut f1, mut g1) = (f, g);
    let mut scale = max(f1.abs(), g1.abs());
    let (mut cs, mut sn, mut r);
    if scale >= safmx2 {
        let mut count = 0;
        loop {
            count += 1;
            f1 *= safmn2;
            g1 *= safmn2;
            scale = max(f1.abs(), g1.abs());
            if scale < safmx2 {
                break;
            }
        }
        r = (f1 * f1 + g1 * g1).sqrt();
        cs = f1 / r;
        sn = g1 / r;
        for _ in 0..count {
            r *= safmx2;
        }
    } else if scale <= safmn2 {
        let mut count = 0;
        loop {
            count += 1;
            f1 *= safmx2;
            g1 *= safmx2;
            scale = max(f1.abs(), g1.abs());
            if scale > safmn2 {
                break;
            }
        }
        r = (f1 * f1 + g1 * g1).sqrt();
        cs = f1 / r;
        sn = g1 / r;
        for _ in 0..count {
            r *= safmn2;
        }
    } else {
        r = (f1 * f1 + g1 * g1).sqrt();
        cs = f1 / r;
        sn = g1 / r;
    }
    if f.abs() > g.abs() && cs < 0.0 {
        cs = -cs;
        sn = -sn;
        r = -r;
    }
    (cs, sn, r)
}

/// `dlas2`: `(ssmin, ssmax)`.
pub(crate) fn dlas2(f: f64, g: f64, h: f64) -> (f64, f64) {
    let (fa, ga, ha) = (f.abs(), g.abs(), h.abs());
    let fhmn = min(fa, ha);
    let fhmx = max(fa, ha);
    if fhmn == 0.0 {
        let ssmax = if fhmx == 0.0 {
            ga
        } else {
            let ratio = min(fhmx, ga) / max(fhmx, ga);
            max(fhmx, ga) * (1.0 + ratio * ratio).sqrt()
        };
        return (0.0, ssmax);
    }
    if ga < fhmx {
        let r#as = 1.0 + fhmn / fhmx;
        let at = (fhmx - fhmn) / fhmx;
        let ratio = ga / fhmx;
        let au = ratio * ratio;
        let c = 2.0 / ((r#as * r#as + au).sqrt() + (at * at + au).sqrt());
        (fhmn * c, fhmx / c)
    } else {
        let au = fhmx / ga;
        if au == 0.0 {
            ((fhmn * fhmx) / ga, ga)
        } else {
            let r#as = 1.0 + fhmn / fhmx;
            let at = (fhmx - fhmn) / fhmx;
            let (x, y) = (r#as * au, at * au);
            let c = 1.0 / ((1.0 + x * x).sqrt() + (1.0 + y * y).sqrt());
            let ssmin = (fhmn * c) * au;
            (ssmin + ssmin, ga / (c + c))
        }
    }
}

/// `dlasv2`: `(ssmin, ssmax, snr, csr, snl, csl)`.
pub(crate) fn dlasv2(f: f64, g: f64, h: f64) -> (f64, f64, f64, f64, f64, f64) {
    let mut ft = f;
    let mut fa = ft.abs();
    let mut ht = h;
    let mut ha = h.abs();
    let mut pmax = 1;
    let swap = ha > fa;
    if swap {
        pmax = 3;
        std::mem::swap(&mut ft, &mut ht);
        std::mem::swap(&mut fa, &mut ha);
    }
    let gt = g;
    let ga = gt.abs();
    let (mut ssmin, mut ssmax);
    let (mut clt, mut crt, mut slt, mut srt);
    if ga == 0.0 {
        ssmin = ha;
        ssmax = fa;
        clt = 1.0;
        crt = 1.0;
        slt = 0.0;
        srt = 0.0;
    } else {
        let mut gasmal = true;
        ssmin = 0.0;
        ssmax = 0.0;
        clt = 0.0;
        crt = 0.0;
        slt = 0.0;
        srt = 0.0;
        if ga > fa {
            pmax = 2;
            if fa / ga < EPS {
                gasmal = false;
                ssmax = ga;
                ssmin = if ha > 1.0 {
                    fa / (ga / ha)
                } else {
                    (fa / ga) * ha
                };
                clt = 1.0;
                slt = ht / gt;
                srt = 1.0;
                crt = ft / gt;
            }
        }
        if gasmal {
            let d = fa - ha;
            let mut l = if d == fa { 1.0 } else { d / fa };
            let m = gt / ft;
            let mut t = 2.0 - l;
            let mm = m * m;
            let tt = t * t;
            let s = (tt + mm).sqrt();
            let r = if l == 0.0 {
                m.abs()
            } else {
                (l * l + mm).sqrt()
            };
            let a = 0.5 * (s + r);
            ssmin = ha / a;
            ssmax = fa * a;
            if mm == 0.0 {
                if l == 0.0 {
                    t = sign(2.0, ft) * sign(1.0, gt);
                } else {
                    t = gt / sign(d, ft) + m / t;
                }
            } else {
                t = (m / (s + t) + m / (r + l)) * (1.0 + a);
            }
            l = (t * t + 4.0).sqrt();
            crt = 2.0 / l;
            srt = t / l;
            clt = (crt + srt * m) / a;
            slt = (ht / ft) * srt / a;
        }
    }
    let (csl, snl, csr, snr) = if swap {
        (srt, crt, slt, clt)
    } else {
        (clt, slt, crt, srt)
    };
    let tsign = match pmax {
        1 => sign(1.0, csr) * sign(1.0, csl) * sign(1.0, f),
        2 => sign(1.0, snr) * sign(1.0, csl) * sign(1.0, g),
        _ => sign(1.0, snr) * sign(1.0, snl) * sign(1.0, h),
    };
    ssmax = sign(ssmax, tsign);
    ssmin = sign(ssmin, tsign * sign(1.0, f) * sign(1.0, h));
    (ssmin, ssmax, snr, csr, snl, csl)
}

/// `ilaenv(1, 'DORMQR', ...)`: past this many reflectors `dormqr` applies them in blocks.
const DORMQR_BLOCK: usize = 32;

/// `dlarft('F', 'C', rows, k, V, tau, T)`, V being `a` from `(first, first)`: the triangular
/// factor of k reflectors, k by k.
fn dlarft(a: &mut [f64], n: usize, first: usize, rows: usize, k: usize, tau: &[f64]) -> Vec<f64> {
    let v = |r: usize, s: usize| (first + r - 2) + (first + s - 2) * n;
    let at = |r: usize, s: usize| (r - 1) + (s - 1) * k;
    let mut t = vec![0.0; k * k];
    for i in 1..=k {
        if tau[i - 1] == 0.0 {
            continue;
        }
        let vii = a[v(i, i)];
        a[v(i, i)] = 1.0;
        // `dgemv('T', rows - i + 1, i - 1, -tau, V(i, 1), V(i, i), 0, T(1, i))`.
        let alpha = -tau[i - 1];
        for j in 1..i {
            let mut temp = 0.0;
            for r in i..=rows {
                temp += a[v(r, j)] * a[v(r, i)];
            }
            t[at(j, i)] = 0.0;
            t[at(j, i)] += alpha * temp;
        }
        a[v(i, i)] = vii;
        // `dtrmv('U', 'N', 'N', i - 1, T, T(1, i))`.
        for j in 1..i {
            let temp = t[at(j, i)];
            if temp != 0.0 {
                for r in 1..j {
                    t[at(r, i)] += temp * t[at(r, j)];
                }
                t[at(j, i)] *= t[at(j, j)];
            }
        }
        t[at(i, i)] = tau[i - 1];
    }
    t
}

/// `dlarfb('L', 'N', 'F', 'C', m, n, k, V, T, C)`, V being `a` from `(first, first)` and C the rows
/// of `c` from `first`: `dtrmm`, `dgemm`, `dtrmm`, `dgemm`, `dtrmm` through the workspace.
fn dlarfb(a: &[f64], n: usize, first: usize, m: usize, k: usize, t: &[f64], c: &mut [f64]) {
    let v = |r: usize, s: usize| (first + r - 2) + (first + s - 2) * n;
    let cc = |r: usize, s: usize| (first + r - 2) + (s - 1) * n;
    let tt = |r: usize, s: usize| (r - 1) + (s - 1) * k;
    let wa = |r: usize, s: usize| (r - 1) + (s - 1) * n;
    let mut w = vec![0.0; n * k];
    for j in 1..=k {
        for r in 1..=n {
            w[wa(r, j)] = c[cc(j, r)];
        }
    }
    // `dtrmm('R', 'L', 'N', 'U')`: W V1.
    for j in 1..=k {
        for l in j + 1..=k {
            let x = a[v(l, j)];
            if x != 0.0 {
                for r in 1..=n {
                    w[wa(r, j)] += x * w[wa(r, l)];
                }
            }
        }
    }
    if m > k {
        // `dgemm('T', 'N', n, k, m - k, 1, C2, V2, 1, W)`.
        for j in 1..=k {
            for r in 1..=n {
                let mut temp = 0.0;
                for l in k + 1..=m {
                    temp += c[cc(l, r)] * a[v(l, j)];
                }
                w[wa(r, j)] = temp + w[wa(r, j)];
            }
        }
    }
    // `dtrmm('R', 'U', 'T', 'N')`: W T'.
    for l in 1..=k {
        for j in 1..l {
            let x = t[tt(j, l)];
            if x != 0.0 {
                for r in 1..=n {
                    w[wa(r, j)] += x * w[wa(r, l)];
                }
            }
        }
        let x = t[tt(l, l)];
        if x != 1.0 {
            for r in 1..=n {
                w[wa(r, l)] *= x;
            }
        }
    }
    if m > k {
        // `dgemm('N', 'T', m - k, n, k, -1, V2, W, 1, C2)`.
        for j in 1..=n {
            for l in 1..=k {
                let x = w[wa(j, l)];
                if x != 0.0 {
                    let temp = -x;
                    for r in k + 1..=m {
                        c[cc(r, j)] += temp * a[v(r, l)];
                    }
                }
            }
        }
    }
    // `dtrmm('R', 'L', 'T', 'U')`: W V1'.
    for l in (1..=k).rev() {
        for j in l + 1..=k {
            let x = a[v(j, l)];
            if x != 0.0 {
                for r in 1..=n {
                    w[wa(r, j)] += x * w[wa(r, l)];
                }
            }
        }
    }
    for j in 1..=k {
        for r in 1..=n {
            c[cc(j, r)] -= w[wa(r, j)];
        }
    }
}

/// `dormqr('L', 'N', n, n, n, A, tau, C)`: the reflectors `dgebd2` left below the diagonal of `a`,
/// last first, applied to `c`. Up to `DORMQR_BLOCK` it is `dorm2r`, one reflector at a time; past
/// it the blocks of 32 from the last, each through `dlarft` and `dlarfb`.
fn dormqr(n: usize, a: &mut [f64], tau: &[f64], c: &mut [f64]) {
    if n <= DORMQR_BLOCK {
        let matrix = Square { n, values: a };
        let mut c = Square { n, values: c };
        for i in (1..=n).rev() {
            let mut v = vec![1.0];
            v.extend((i + 1..=n).map(|r| matrix.values[matrix.at(r, i)]));
            dlarf_left(&v, tau[i - 1], &mut c, i, 1);
        }
        return;
    }
    let mut i = ((n - 1) / DORMQR_BLOCK) * DORMQR_BLOCK + 1;
    loop {
        let ib = DORMQR_BLOCK.min(n - i + 1);
        let t = dlarft(a, n, i, n - i + 1, ib, &tau[i - 1..]);
        dlarfb(a, n, i, n - i + 1, ib, &t, c);
        if i == 1 {
            break;
        }
        i -= DORMQR_BLOCK;
    }
}

/// Whether [`dgesdd`] takes this matrix down the path the reference takes: small enough, finite,
/// and inside the range `dgesdd` does not scale. `SMLNUM`, `sqrt(safe minimum) / precision`.
pub fn reachable(n: usize, values: &[f64]) -> bool {
    let smlnum = SAFE_MINIMUM.sqrt() / PRECISION;
    let bignum = 1.0 / smlnum;
    let mut anrm = 0.0f64;
    for value in values {
        if !value.is_finite() {
            return false;
        }
        anrm = max(anrm, value.abs());
    }
    (2..=LARGEST).contains(&n) && (anrm == 0.0 || (anrm >= smlnum && anrm <= bignum))
}

/// The singular values, decreasing, and the left singular vectors, column major. `None` where
/// `dbdsdc` did not converge, which is `info > 0` in the reference and a `NotConvergedException` in
/// Breeze. `a` is column major, and is overwritten.
pub fn dgesdd(n: usize, a: &mut [f64]) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut matrix = Square { n, values: a };
    let (mut d, mut e, tauq) = dgebd2(&mut matrix);
    let mut u = vec![0.0; n * n];
    let mut vt = vec![0.0; n * n];
    if crate::f2j_dbdsdc::dbdsdc(n, &mut d, &mut e, &mut u, &mut vt) != 0 {
        return None;
    }
    dormqr(n, matrix.values, &tauq, &mut u);
    Some((d, u))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn product(n: usize, a: &[f64], b: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    out[i + j * n] += a[i + k * n] * b[k + j * n];
                }
            }
        }
        out
    }

    /// The values are those of the matrix: `G = U S U'` for a symmetric positive semi-definite G.
    #[test]
    fn a_gramian_comes_back_from_its_values_and_vectors() {
        let n = 5;
        let a: Vec<f64> = (0..n * n)
            .map(|i| ((i * 37 % 11) as f64 - 5.0) / 3.0)
            .collect();
        let mut transpose = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                transpose[i + j * n] = a[j + i * n];
            }
        }
        let g = product(n, &transpose, &a);
        let mut work = g.clone();
        let (s, u) = dgesdd(n, &mut work).expect("it converges");
        assert!(s.windows(2).all(|w| w[0] >= w[1]));
        let mut us = u.clone();
        for j in 0..n {
            for i in 0..n {
                us[i + j * n] *= s[j];
            }
        }
        let mut ut = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                ut[i + j * n] = u[j + i * n];
            }
        }
        let back = product(n, &us, &ut);
        for (x, y) in back.iter().zip(&g) {
            assert!((x - y).abs() < 1e-12, "{x} vs {y}");
        }
    }

    /// A zero column is zero through every reflector, so its value is EXACTLY zero, which is what
    /// lets a flat sample lower the rank whatever order the Gramian was summed in.
    #[test]
    fn a_zero_row_and_column_give_an_exact_zero() {
        let n = 3;
        let mut g = vec![2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 3.0];
        let (s, _) = dgesdd(n, &mut g).expect("it converges");
        assert_eq!(s[2], 0.0);
        assert!(s[1] > 0.0);
    }

    #[test]
    fn the_range_it_does_not_scale_is_the_one_it_takes() {
        assert!(reachable(2, &[1.0, 0.0, 0.0, 1.0]));
        assert!(!reachable(129, &[1.0]));
        assert!(!reachable(2, &[f64::NAN, 0.0, 0.0, 1.0]));
        assert!(!reachable(2, &[1e-300, 0.0, 0.0, 1e-300]));
    }
}
