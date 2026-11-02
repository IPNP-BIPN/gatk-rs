//! `dgesdd('A', n, n, ...)` as the F2J translation of LAPACK runs it, for `n` up to 25.
//!
//! Spark's `RowMatrix.computeSVD` answers a panel of normals from `breeze.linalg.svd` over the
//! Gramian, and in the oracle image that is `dev.ludovic.netlib.lapack.F2jLAPACK`: no native LAPACK
//! loads, so the arithmetic is `org.netlib.lapack.Dgesdd` in a strict JVM, and its bits are a
//! function of its input. The Gramian of a standardised panel is rank deficient, so its trailing
//! values are rounding noise, and the noise is what `rCond` compares: only the same operations in
//! the same order reach the same rank.
//!
//! The path is the one a square matrix takes up to `SMLSIZ` (25) in `dbdsdc`: `dgebd2` (`dgebrd`
//! leaves its blocked code above `NX`, 128), then `dbdsqr`, then `dorm2r` for the left vectors
//! (`dormqr` leaves its blocked code below `NB`, 32). A larger matrix is divide and conquer
//! (`dlasd0` and the rest), which is not ported; the caller asks [`reachable`] first.
//!
//! The vectors on the right are not accumulated. They feed nothing back into the values, and Spark
//! throws them away.
//!
//! Ported from LAPACK 3.1.1 (`SRC/dgesdd.f`, `dgebd2.f`, `dlarfg.f`, `dlarf.f`, `dbdsdc.f`,
//! `dbdsqr.f`, `dlasr.f`, `dlas2.f`, `dlasv2.f`, `dlartg.f`, `dlapy2.f`, `dorm2r.f`) and reference
//! BLAS (`dnrm2.f`, `dgemv.f`, `dger.f`), BSD 3-Clause, with what F2J does to them: `x**2` is
//! `Math.pow(x, 2)`, which HotSpot answers as `x * x`.

/// `SMLSIZ`, `ilaenv(9, ...)`: the largest matrix `dbdsdc` hands to `dbdsqr`.
pub const LARGEST: usize = 25;

/// `dlamch('E')`: F2J's is the rounding unit.
const EPS: f64 = f64::EPSILON / 2.0;
/// `dlamch('P')`, the epsilon times the base.
const PRECISION: f64 = f64::EPSILON;
/// `dlamch('S')`.
const SAFE_MINIMUM: f64 = f64::MIN_POSITIVE;
/// `Math.pow(EPS, -0.125)`, measured on the oracle: 2^(53/8), under `dbdsqr`'s cap of a hundred.
const TOLERANCE_MULTIPLIER: f64 = f64::from_bits(0x4058_ACE5_422A_A0DB);
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
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// Fortran's `MAX`, as F2J writes it: `a >= b ? a : b` over values that are not NaN here.
fn max(a: f64, b: f64) -> f64 {
    if a >= b {
        a
    } else {
        b
    }
}

/// Fortran's `MIN`.
fn min(a: f64, b: f64) -> f64 {
    if a <= b {
        a
    } else {
        b
    }
}

/// `dnrm2`, the scaled sum of squares rather than a square root of a dot product.
fn dnrm2(x: &[f64]) -> f64 {
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
fn dlapy2(x: f64, y: f64) -> f64 {
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
fn dlartg(f: f64, g: f64) -> (f64, f64, f64) {
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
fn dlas2(f: f64, g: f64, h: f64) -> (f64, f64) {
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
fn dlasv2(f: f64, g: f64, h: f64) -> (f64, f64, f64, f64, f64, f64) {
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

/// `dlasr('R', 'V', direct, n, count, c, s, u(1, first))`: rotate the columns from `first`.
fn rotate_columns(u: &mut Square, first: usize, c: &[f64], s: &[f64], forward: bool) {
    let n = u.n;
    let apply = |u: &mut Square, j: usize| {
        let (ctemp, stemp) = (c[j], s[j]);
        if ctemp != 1.0 || stemp != 0.0 {
            for i in 1..=n {
                let left = u.at(i, first + j);
                let right = u.at(i, first + j + 1);
                let temp = u.values[right];
                u.values[right] = ctemp * temp - stemp * u.values[left];
                u.values[left] = stemp * temp + ctemp * u.values[left];
            }
        }
    };
    if forward {
        for j in 0..c.len() {
            apply(u, j);
        }
    } else {
        for j in (0..c.len()).rev() {
            apply(u, j);
        }
    }
}

/// `dbdsqr('U', n, ncvt > 0, nru = n, 0, ...)`, without accumulating the right vectors: the
/// singular values in `d`, decreasing, and `u` rotated. `false` when it did not converge.
fn dbdsqr(d: &mut [f64], e: &mut [f64], u: &mut Square) -> bool {
    let n = d.len();
    // One-based, as the Fortran's two pointers are.
    macro_rules! d {
        ($i:expr) => {
            d[$i - 1]
        };
    }
    macro_rules! e {
        ($i:expr) => {
            e[$i - 1]
        };
    }
    const MAXITR: usize = 6;
    if n > 1 {
        let unfl = SAFE_MINIMUM;
        let tolmul = max(10.0, min(100.0, TOLERANCE_MULTIPLIER));
        let tol = tolmul * EPS;
        let mut smax = 0.0;
        for i in 1..=n {
            smax = max(smax, d!(i).abs());
        }
        for i in 1..n {
            smax = max(smax, e!(i).abs());
        }
        let mut sminl;
        let mut sminoa = d!(1).abs();
        if sminoa != 0.0 {
            let mut mu = sminoa;
            for i in 2..=n {
                mu = d!(i).abs() * (mu / (mu + e!(i - 1).abs()));
                sminoa = min(sminoa, mu);
                if sminoa == 0.0 {
                    break;
                }
            }
        }
        sminoa /= (n as f64).sqrt();
        let thresh = max(tol * sminoa, (MAXITR * n * n) as f64 * unfl);
        let maxit = MAXITR * n * n;
        let mut iter = 0;
        let mut oldll: isize = -1;
        let mut oldm: isize = -1;
        let mut idir = 0;
        let mut m = n;
        'main: loop {
            if m <= 1 {
                break;
            }
            if iter > maxit {
                return false;
            }
            smax = d!(m).abs();
            let mut ll = 0;
            let mut split = false;
            for lll in 1..m {
                let candidate = m - lll;
                let abss = d!(candidate).abs();
                let abse = e!(candidate).abs();
                if abse <= thresh {
                    ll = candidate;
                    split = true;
                    break;
                }
                smax = max(max(smax, abss), abse);
            }
            if split {
                e!(ll) = 0.0;
                if ll == m - 1 {
                    m -= 1;
                    continue 'main;
                }
            }
            ll += 1;
            if ll == m - 1 {
                let (sigmn, sigmx, _sinr, _cosr, sinl, cosl) = dlasv2(d!(m - 1), e!(m - 1), d!(m));
                d!(m - 1) = sigmx;
                e!(m - 1) = 0.0;
                d!(m) = sigmn;
                // `drot` over the two columns of U.
                for i in 1..=n {
                    let (x, y) = (u.at(i, m - 1), u.at(i, m));
                    let temp = cosl * u.values[x] + sinl * u.values[y];
                    u.values[y] = cosl * u.values[y] - sinl * u.values[x];
                    u.values[x] = temp;
                }
                m -= 2;
                continue 'main;
            }
            if (ll as isize) > oldm || (m as isize) < oldll {
                idir = if d!(ll).abs() >= d!(m).abs() { 1 } else { 2 };
            }
            if idir == 1 {
                if e!(m - 1).abs() <= tol.abs() * d!(m).abs() {
                    e!(m - 1) = 0.0;
                    continue 'main;
                }
                let mut mu = d!(ll).abs();
                sminl = mu;
                for lll in ll..m {
                    if e!(lll).abs() <= tol * mu {
                        e!(lll) = 0.0;
                        continue 'main;
                    }
                    mu = d!(lll + 1).abs() * (mu / (mu + e!(lll).abs()));
                    sminl = min(sminl, mu);
                }
            } else {
                if e!(ll).abs() <= tol.abs() * d!(ll).abs() {
                    e!(ll) = 0.0;
                    continue 'main;
                }
                let mut mu = d!(m).abs();
                sminl = mu;
                for lll in (ll..m).rev() {
                    if e!(lll).abs() <= tol * mu {
                        e!(lll) = 0.0;
                        continue 'main;
                    }
                    mu = d!(lll).abs() * (mu / (mu + e!(lll).abs()));
                    sminl = min(sminl, mu);
                }
            }
            oldll = ll as isize;
            oldm = m as isize;

            let mut shift;
            if (n as f64) * tol * (sminl / smax) <= max(EPS, 0.01 * tol) {
                shift = 0.0;
            } else {
                let sll;
                if idir == 1 {
                    sll = d!(ll).abs();
                    shift = dlas2(d!(m - 1), e!(m - 1), d!(m)).0;
                } else {
                    sll = d!(m).abs();
                    shift = dlas2(d!(ll), e!(ll), d!(ll + 1)).0;
                }
                if sll > 0.0 {
                    let ratio = shift / sll;
                    if ratio * ratio < EPS {
                        shift = 0.0;
                    }
                }
            }
            iter += m - ll;
            let count = m - ll;
            // WORK(1..) and WORK(N..), then WORK(NM12+1..) and WORK(NM13+1..).
            let mut first_c = vec![0.0; count];
            let mut first_s = vec![0.0; count];
            let mut second_c = vec![0.0; count];
            let mut second_s = vec![0.0; count];
            if shift == 0.0 {
                if idir == 1 {
                    let (mut cs, mut oldcs) = (1.0, 1.0);
                    let mut oldsn = 0.0;
                    for i in ll..m {
                        let (c, sn, r) = dlartg(d!(i) * cs, e!(i));
                        cs = c;
                        if i > ll {
                            e!(i - 1) = oldsn * r;
                        }
                        let (c, s, r2) = dlartg(oldcs * r, d!(i + 1) * sn);
                        oldcs = c;
                        oldsn = s;
                        d!(i) = r2;
                        first_c[i - ll] = cs;
                        first_s[i - ll] = sn;
                        second_c[i - ll] = oldcs;
                        second_s[i - ll] = oldsn;
                    }
                    let h = d!(m) * cs;
                    d!(m) = h * oldcs;
                    e!(m - 1) = h * oldsn;
                    rotate_columns(u, ll, &second_c, &second_s, true);
                    if e!(m - 1).abs() <= thresh {
                        e!(m - 1) = 0.0;
                    }
                } else {
                    let (mut cs, mut oldcs) = (1.0, 1.0);
                    let mut oldsn = 0.0;
                    for i in (ll + 1..=m).rev() {
                        let (c, sn, r) = dlartg(d!(i) * cs, e!(i - 1));
                        cs = c;
                        if i < m {
                            e!(i) = oldsn * r;
                        }
                        let (c, s, r2) = dlartg(oldcs * r, d!(i - 1) * sn);
                        oldcs = c;
                        oldsn = s;
                        d!(i) = r2;
                        first_c[i - ll - 1] = cs;
                        first_s[i - ll - 1] = -sn;
                        second_c[i - ll - 1] = oldcs;
                        second_s[i - ll - 1] = -oldsn;
                    }
                    let h = d!(ll) * cs;
                    d!(ll) = h * oldcs;
                    e!(ll) = h * oldsn;
                    rotate_columns(u, ll, &first_c, &first_s, false);
                    if e!(ll).abs() <= thresh {
                        e!(ll) = 0.0;
                    }
                }
            } else if idir == 1 {
                let mut f = (d!(ll).abs() - shift) * (sign(1.0, d!(ll)) + shift / d!(ll));
                let mut g = e!(ll);
                for i in ll..m {
                    let (cosr, sinr, r) = dlartg(f, g);
                    if i > ll {
                        e!(i - 1) = r;
                    }
                    f = cosr * d!(i) + sinr * e!(i);
                    e!(i) = cosr * e!(i) - sinr * d!(i);
                    g = sinr * d!(i + 1);
                    d!(i + 1) *= cosr;
                    let (cosl, sinl, r) = dlartg(f, g);
                    d!(i) = r;
                    f = cosl * e!(i) + sinl * d!(i + 1);
                    d!(i + 1) = cosl * d!(i + 1) - sinl * e!(i);
                    if i < m - 1 {
                        g = sinl * e!(i + 1);
                        e!(i + 1) *= cosl;
                    }
                    second_c[i - ll] = cosl;
                    second_s[i - ll] = sinl;
                }
                e!(m - 1) = f;
                rotate_columns(u, ll, &second_c, &second_s, true);
                if e!(m - 1).abs() <= thresh {
                    e!(m - 1) = 0.0;
                }
            } else {
                let mut f = (d!(m).abs() - shift) * (sign(1.0, d!(m)) + shift / d!(m));
                let mut g = e!(m - 1);
                for i in (ll + 1..=m).rev() {
                    let (cosr, sinr, r) = dlartg(f, g);
                    if i < m {
                        e!(i) = r;
                    }
                    f = cosr * d!(i) + sinr * e!(i - 1);
                    e!(i - 1) = cosr * e!(i - 1) - sinr * d!(i);
                    g = sinr * d!(i - 1);
                    d!(i - 1) *= cosr;
                    let (cosl, sinl, r) = dlartg(f, g);
                    d!(i) = r;
                    f = cosl * e!(i - 1) + sinl * d!(i - 1);
                    d!(i - 1) = cosl * d!(i - 1) - sinl * e!(i - 1);
                    if i > ll + 1 {
                        g = sinl * e!(i - 2);
                        e!(i - 2) *= cosl;
                    }
                    first_c[i - ll - 1] = cosr;
                    first_s[i - ll - 1] = -sinr;
                }
                e!(ll) = f;
                if e!(ll).abs() <= thresh {
                    e!(ll) = 0.0;
                }
                rotate_columns(u, ll, &first_c, &first_s, false);
            }
        }
    }
    for value in d.iter_mut() {
        if *value < 0.0 {
            *value = -*value;
        }
    }
    // Smallest to the end, one transposition per vector.
    for i in 1..n {
        let last = n + 1 - i;
        let mut isub = 1;
        let mut smin = d[0];
        for j in 2..=last {
            if d[j - 1] <= smin {
                isub = j;
                smin = d[j - 1];
            }
        }
        if isub != last {
            d[isub - 1] = d[last - 1];
            d[last - 1] = smin;
            for row in 1..=n {
                let (x, y) = (u.at(row, isub), u.at(row, last));
                u.values.swap(x, y);
            }
        }
    }
    true
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
/// `dbdsqr` did not converge, which is `info > 0` in the reference and a `NotConvergedException` in
/// Breeze. `a` is column major, and is overwritten.
pub fn dgesdd(n: usize, a: &mut [f64]) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut matrix = Square { n, values: a };
    let (mut d, mut e, tauq) = dgebd2(&mut matrix);
    let mut u_values = vec![0.0; n * n];
    for i in 0..n {
        u_values[i + i * n] = 1.0;
    }
    let mut u = Square {
        n,
        values: &mut u_values,
    };
    if !dbdsqr(&mut d, &mut e, &mut u) {
        return None;
    }
    let swap_columns = |u: &mut Square, a: usize, b: usize| {
        for row in 1..=n {
            let (x, y) = (u.at(row, a), u.at(row, b));
            u.values.swap(x, y);
        }
    };
    // `dlasdq` sorts what `dbdsqr` sorted into ascending order, and `dbdsdc` back. The values end
    // where they were; a tie's vectors do not.
    for i in 1..=n {
        let mut isub = i;
        let mut smin = d[i - 1];
        for j in i + 1..=n {
            if d[j - 1] < smin {
                isub = j;
                smin = d[j - 1];
            }
        }
        if isub != i {
            d[isub - 1] = d[i - 1];
            d[i - 1] = smin;
            swap_columns(&mut u, isub, i);
        }
    }
    for i in 1..n {
        let mut kk = i;
        let mut p = d[i - 1];
        for j in i + 1..=n {
            if d[j - 1] > p {
                kk = j;
                p = d[j - 1];
            }
        }
        if kk != i {
            d[kk - 1] = d[i - 1];
            d[i - 1] = p;
            swap_columns(&mut u, i, kk);
        }
    }
    // `dormbr('Q', 'L', 'N')` is `dorm2r` backwards, H(1) ... H(k) applied to U.
    for i in (1..=n).rev() {
        let mut v = vec![1.0];
        v.extend((i + 1..=n).map(|r| matrix.values[matrix.at(r, i)]));
        dlarf_left(&v, tauq[i - 1], &mut u, i, 1);
    }
    Some((d, u_values))
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
        assert!(!reachable(26, &[1.0]));
        assert!(!reachable(2, &[f64::NAN, 0.0, 0.0, 1.0]));
        assert!(!reachable(2, &[1e-300, 0.0, 0.0, 1e-300]));
    }
}
