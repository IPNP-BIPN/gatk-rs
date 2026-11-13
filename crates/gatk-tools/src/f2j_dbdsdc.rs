//! `dbdsdc('U', 'I', n, ...)` as the F2J translation of LAPACK 3.1.1 runs it: the singular values of
//! an upper bidiagonal matrix and both sets of vectors, which [`crate::f2j_gesdd`] calls.
//!
//! Up to `SMLSIZ` (25) it is `dlasdq`, that is `dbdsqr`. Past it the matrix is scaled to a largest
//! entry of one, every diagonal entry under the epsilon is raised to it, and `dlasd0` divides: the
//! tree of `dlasdt`, `dlasdq` on its leaves, and `dlasd1` merging two halves through a row, which
//! deflates (`dlasd2`), solves the secular equation (`dlasd3`, `dlasd4`, `dlasd5`, `dlaed6`) and
//! multiplies the vectors back. The vectors on the right are not optional here, as they are for
//! `dbdsqr` alone: each merge reads its `z` from them.
//!
//! The raising is why a Gramian of more than 25 samples keeps every value under Spark's `rCond`:
//! its rounding noise becomes the epsilon times the largest entry, which is 1e-8 of the first
//! singular value once Spark takes the square root, and the threshold is 1e-9.
//!
//! Ported from LAPACK 3.1.1 (`SRC/dbdsdc.f`, `dbdsqr.f`, `dlasdq.f`, `dlasd0.f` to `dlasd5.f`,
//! `dlasdt.f`, `dlaed6.f`, `dlamrg.f`, `dlascl.f`, `dlanst.f`, `dlasr.f`) and reference BLAS
//! (`drot.f`, `dswap.f`, `dscal.f`, `dcopy.f`, `dgemm.f`), BSD 3-Clause, with what F2J does to them:
//! `MAX` and `MIN` are `Math.max` and `Math.min`, `SIGN` is `Util.dsign`, `x**2` is `Math.pow(x, 2)`,
//! which HotSpot answers as `x * x`, and `DLAMC3(a, b)` is `a + b`.

use crate::f2j_gesdd::{
    dlapy2, dlartg, dlas2, dlasv2, dnrm2, max, max3, min, sign, EPS, SAFE_MINIMUM,
    TOLERANCE_MULTIPLIER,
};

/// `ilaenv(9, ...)`: the largest problem `dbdsdc` and `dlasd0` hand to `dlasdq`.
pub const SMLSIZ: usize = 25;

/// Where `(1, 1)` of a column-major block sits in its buffer, and the buffer's leading dimension.
#[derive(Clone, Copy)]
pub(crate) struct View {
    pub at: usize,
    pub ld: usize,
}

impl View {
    pub fn new(ld: usize) -> View {
        View { at: 0, ld }
    }

    /// `A(i, j)`, one-based.
    pub fn ix(self, i: usize, j: usize) -> usize {
        self.at + (i - 1) + (j - 1) * self.ld
    }

    /// The block from `A(i, j)`, as Fortran passes `A(I, J)`.
    pub fn sub(self, i: usize, j: usize) -> View {
        View {
            at: self.ix(i, j),
            ld: self.ld,
        }
    }
}

/// `drot` over two rows of `a`, `count` columns long.
fn rot_rows(a: &mut [f64], v: View, x: usize, y: usize, count: usize, c: f64, s: f64) {
    for j in 1..=count {
        let (ix, iy) = (v.ix(x, j), v.ix(y, j));
        let temp = c * a[ix] + s * a[iy];
        a[iy] = c * a[iy] - s * a[ix];
        a[ix] = temp;
    }
}

/// `drot` over two columns of `a`, `count` rows long.
fn rot_columns(a: &mut [f64], v: View, x: usize, y: usize, count: usize, c: f64, s: f64) {
    for i in 1..=count {
        let (ix, iy) = (v.ix(i, x), v.ix(i, y));
        let temp = c * a[ix] + s * a[iy];
        a[iy] = c * a[iy] - s * a[ix];
        a[ix] = temp;
    }
}

fn swap_rows(a: &mut [f64], v: View, x: usize, y: usize, count: usize) {
    for j in 1..=count {
        a.swap(v.ix(x, j), v.ix(y, j));
    }
}

fn swap_columns(a: &mut [f64], v: View, x: usize, y: usize, count: usize) {
    for i in 1..=count {
        a.swap(v.ix(i, x), v.ix(i, y));
    }
}

/// `dlasr(side, 'V', direct, m, n, c, s, a)`: rotations between neighbouring rows (`left`) or
/// columns, the first pair first (`forward`) or last.
#[allow(clippy::too_many_arguments)]
fn dlasr(
    left: bool,
    forward: bool,
    m: usize,
    n: usize,
    c: &[f64],
    s: &[f64],
    a: &mut [f64],
    v: View,
) {
    if m == 0 || n == 0 {
        return;
    }
    let count = if left { m - 1 } else { n - 1 };
    for step in 0..count {
        let j = if forward { step + 1 } else { count - step };
        let (ctemp, stemp) = (c[j - 1], s[j - 1]);
        if ctemp != 1.0 || stemp != 0.0 {
            if left {
                for i in 1..=n {
                    let (lo, hi) = (v.ix(j, i), v.ix(j + 1, i));
                    let temp = a[hi];
                    a[hi] = ctemp * temp - stemp * a[lo];
                    a[lo] = stemp * temp + ctemp * a[lo];
                }
            } else {
                for i in 1..=m {
                    let (lo, hi) = (v.ix(i, j), v.ix(i, j + 1));
                    let temp = a[hi];
                    a[hi] = ctemp * temp - stemp * a[lo];
                    a[lo] = stemp * temp + ctemp * a[lo];
                }
            }
        }
    }
}

/// `dlascl('G', 0, 0, cfrom, cto, ...)`: one multiplication by `cto / cfrom` unless that would
/// leave the range, then in steps.
fn dlascl(cfrom: f64, cto: f64, values: &mut [f64]) {
    let smlnum = SAFE_MINIMUM;
    let bignum = 1.0 / smlnum;
    let (mut cfromc, mut ctoc) = (cfrom, cto);
    loop {
        let cfrom1 = cfromc * smlnum;
        let cto1 = ctoc / bignum;
        let (mul, done);
        if cfrom1.abs() > ctoc.abs() && ctoc != 0.0 {
            mul = smlnum;
            done = false;
            cfromc = cfrom1;
        } else if cto1.abs() > cfromc.abs() {
            mul = bignum;
            done = false;
            ctoc = cto1;
        } else {
            mul = ctoc / cfromc;
            done = true;
        }
        for value in values.iter_mut() {
            *value *= mul;
        }
        if done {
            return;
        }
    }
}

/// `dbdsqr('U', n, ncvt, nru, 0, ...)` with vectors: the values in `d`, decreasing, rows of `vt`
/// and columns of `u` rotated with them. Returns `info`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn dbdsqr(
    n: usize,
    ncvt: usize,
    nru: usize,
    d: &mut [f64],
    e: &mut [f64],
    vt: &mut [f64],
    vtv: View,
    u: &mut [f64],
    uv: View,
) -> usize {
    const MAXITR: usize = 6;
    if n == 0 {
        return 0;
    }
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
    if n > 1 {
        // `ncvt` and `nru` are never both zero here, so the qd algorithm (`dlasq1`) is not taken.
        let nm1 = n - 1;
        let mut work = vec![0.0; 4 * nm1];
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
                return (1..n).filter(|&i| e!(i) != 0.0).count();
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
                smax = max3(smax, abss, abse);
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
                let (sigmn, sigmx, sinr, cosr, sinl, cosl) = dlasv2(d!(m - 1), e!(m - 1), d!(m));
                d!(m - 1) = sigmx;
                e!(m - 1) = 0.0;
                d!(m) = sigmn;
                if ncvt > 0 {
                    rot_rows(vt, vtv, m - 1, m, ncvt, cosr, sinr);
                }
                if nru > 0 {
                    rot_columns(u, uv, m - 1, m, nru, cosl, sinl);
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
            let rows = m - ll + 1;
            // WORK(1), WORK(N), WORK(NM12 + 1) and WORK(NM13 + 1).
            let (first_c, rest) = work.split_at_mut(nm1);
            let (first_s, rest) = rest.split_at_mut(nm1);
            let (second_c, second_s) = rest.split_at_mut(nm1);
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
                    if ncvt > 0 {
                        dlasr(true, true, rows, ncvt, first_c, first_s, vt, vtv.sub(ll, 1));
                    }
                    if nru > 0 {
                        dlasr(false, true, nru, rows, second_c, second_s, u, uv.sub(1, ll));
                    }
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
                    if ncvt > 0 {
                        dlasr(
                            true,
                            false,
                            rows,
                            ncvt,
                            second_c,
                            second_s,
                            vt,
                            vtv.sub(ll, 1),
                        );
                    }
                    if nru > 0 {
                        dlasr(false, false, nru, rows, first_c, first_s, u, uv.sub(1, ll));
                    }
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
                    first_c[i - ll] = cosr;
                    first_s[i - ll] = sinr;
                    second_c[i - ll] = cosl;
                    second_s[i - ll] = sinl;
                }
                e!(m - 1) = f;
                if ncvt > 0 {
                    dlasr(true, true, rows, ncvt, first_c, first_s, vt, vtv.sub(ll, 1));
                }
                if nru > 0 {
                    dlasr(false, true, nru, rows, second_c, second_s, u, uv.sub(1, ll));
                }
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
                    second_c[i - ll - 1] = cosl;
                    second_s[i - ll - 1] = -sinl;
                }
                e!(ll) = f;
                if e!(ll).abs() <= thresh {
                    e!(ll) = 0.0;
                }
                if ncvt > 0 {
                    dlasr(
                        true,
                        false,
                        rows,
                        ncvt,
                        second_c,
                        second_s,
                        vt,
                        vtv.sub(ll, 1),
                    );
                }
                if nru > 0 {
                    dlasr(false, false, nru, rows, first_c, first_s, u, uv.sub(1, ll));
                }
            }
        }
    }
    for i in 1..=n {
        if d!(i) < 0.0 {
            d!(i) = -d!(i);
            // `dscal` by -1.
            for j in 1..=ncvt {
                let at = vtv.ix(i, j);
                vt[at] *= -1.0;
            }
        }
    }
    // Smallest to the end, one transposition per vector.
    for i in 1..n {
        let last = n + 1 - i;
        let mut isub = 1;
        let mut smin = d!(1);
        for j in 2..=last {
            if d!(j) <= smin {
                isub = j;
                smin = d!(j);
            }
        }
        if isub != last {
            d!(isub) = d!(last);
            d!(last) = smin;
            swap_rows(vt, vtv, isub, last, ncvt);
            swap_columns(u, uv, isub, last, nru);
        }
    }
    0
}

/// `dlasdq('U', sqre, n, ncvt, nru, 0, ...)`: an upper bidiagonal matrix with `sqre` more columns
/// than rows made square, `dbdsqr`, then the values in INCREASING order. Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dlasdq(
    sqre: usize,
    n: usize,
    ncvt: usize,
    nru: usize,
    d: &mut [f64],
    e: &mut [f64],
    vt: &mut [f64],
    vtv: View,
    u: &mut [f64],
    uv: View,
) -> usize {
    if n == 0 {
        return 0;
    }
    let np1 = n + 1;
    // WORK(1..N) and WORK(NP1..).
    let mut cs = vec![0.0; n];
    let mut sn = vec![0.0; n];
    if sqre == 1 {
        // Non-square upper bidiagonal, rotated to square lower bidiagonal from the right.
        for i in 1..n {
            let (c, s, r) = dlartg(d[i - 1], e[i - 1]);
            d[i - 1] = r;
            e[i - 1] = s * d[i];
            d[i] *= c;
            cs[i - 1] = c;
            sn[i - 1] = s;
        }
        let (c, s, r) = dlartg(d[n - 1], e[n - 1]);
        d[n - 1] = r;
        e[n - 1] = 0.0;
        cs[n - 1] = c;
        sn[n - 1] = s;
        if ncvt > 0 {
            dlasr(true, true, np1, ncvt, &cs, &sn, vt, vtv);
        }
        // Lower bidiagonal, rotated back to upper from the left.
        for i in 1..n {
            let (c, s, r) = dlartg(d[i - 1], e[i - 1]);
            d[i - 1] = r;
            e[i - 1] = s * d[i];
            d[i] *= c;
            cs[i - 1] = c;
            sn[i - 1] = s;
        }
        if nru > 0 {
            dlasr(false, true, nru, n, &cs, &sn, u, uv);
        }
    }
    let info = dbdsqr(n, ncvt, nru, d, e, vt, vtv, u, uv);
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
            swap_rows(vt, vtv, isub, i, ncvt);
            swap_columns(u, uv, isub, i, nru);
        }
    }
    info
}

/// `dlamrg`: the permutation that merges two sorted runs of `a` into increasing order, one-based.
/// A run read backwards has a stride of -1.
fn dlamrg(n1: usize, n2: usize, a: &[f64], forward1: bool, forward2: bool, index: &mut [usize]) {
    let (mut n1sv, mut n2sv) = (n1, n2);
    let (step1, step2): (isize, isize) =
        (if forward1 { 1 } else { -1 }, if forward2 { 1 } else { -1 });
    let mut ind1 = if forward1 { 1 } else { n1 as isize };
    let mut ind2 = if forward2 {
        1 + n1 as isize
    } else {
        (n1 + n2) as isize
    };
    let mut i = 0;
    while n1sv > 0 && n2sv > 0 {
        if a[ind1 as usize - 1] <= a[ind2 as usize - 1] {
            index[i] = ind1 as usize;
            ind1 += step1;
            n1sv -= 1;
        } else {
            index[i] = ind2 as usize;
            ind2 += step2;
            n2sv -= 1;
        }
        i += 1;
    }
    if n1sv == 0 {
        for _ in 0..n2sv {
            index[i] = ind2 as usize;
            ind2 += step2;
            i += 1;
        }
    } else {
        for _ in 0..n1sv {
            index[i] = ind1 as usize;
            ind1 += step1;
            i += 1;
        }
    }
}

/// `dlasdt`: the tree. One-based arrays of each node's centre row and the sizes on its left and
/// right, and the number of levels.
///
/// `LOG` is `Math.log`, here `jmath::math::log`.
fn dlasdt(n: usize, msub: usize) -> (usize, usize, Vec<usize>, Vec<usize>, Vec<usize>) {
    let maxn = n.max(1);
    let temp = jmath::math::log(maxn as f64 / (msub + 1) as f64) / jmath::math::log(2.0);
    let lvl = temp as usize + 1;
    let mut inode = vec![0; n + 1];
    let mut ndiml = vec![0; n + 1];
    let mut ndimr = vec![0; n + 1];
    let i = n / 2;
    inode[1] = i + 1;
    ndiml[1] = i;
    ndimr[1] = n - i - 1;
    let (mut il, mut ir) = (0, 1);
    let mut llst = 1;
    for _ in 1..lvl {
        for i in 0..llst {
            il += 2;
            ir += 2;
            let ncrnt = llst + i;
            ndiml[il] = ndiml[ncrnt] / 2;
            ndimr[il] = ndiml[ncrnt] - ndiml[il] - 1;
            inode[il] = inode[ncrnt] - ndimr[il] - 1;
            ndiml[ir] = ndimr[ncrnt] / 2;
            ndimr[ir] = ndimr[ncrnt] - ndiml[ir] - 1;
            inode[ir] = inode[ncrnt] + ndiml[ir] + 1;
        }
        llst *= 2;
    }
    (lvl, llst * 2 - 1, inode, ndiml, ndimr)
}

/// `dlasd5`: the `i`th root of the secular equation of order two. Returns sigma.
fn dlasd5(i: usize, d: &[f64], z: &[f64], delta: &mut [f64], rho: f64, work: &mut [f64]) -> f64 {
    let del = d[1] - d[0];
    let delsq = del * (d[1] + d[0]);
    let tau = if i == 1 {
        let w = 1.0
            + 4.0 * rho * (z[1] * z[1] / (d[0] + 3.0 * d[1]) - z[0] * z[0] / (3.0 * d[0] + d[1]))
                / del;
        if w > 0.0 {
            let b = delsq + rho * (z[0] * z[0] + z[1] * z[1]);
            let c = rho * z[0] * z[0] * delsq;
            let t = 2.0 * c / (b + (b * b - 4.0 * c).abs().sqrt());
            let tau = t / (d[0] + (d[0] * d[0] + t).sqrt());
            delta[0] = -tau;
            delta[1] = del - tau;
            work[0] = 2.0 * d[0] + tau;
            work[1] = (d[0] + tau) + d[1];
            return d[0] + tau;
        }
        let b = -delsq + rho * (z[0] * z[0] + z[1] * z[1]);
        let c = rho * z[1] * z[1] * delsq;
        let t = if b > 0.0 {
            -2.0 * c / (b + (b * b + 4.0 * c).sqrt())
        } else {
            (b - (b * b + 4.0 * c).sqrt()) / 2.0
        };
        t / (d[1] + (d[1] * d[1] + t).abs().sqrt())
    } else {
        let b = -delsq + rho * (z[0] * z[0] + z[1] * z[1]);
        let c = rho * z[1] * z[1] * delsq;
        let t = if b > 0.0 {
            (b + (b * b + 4.0 * c).sqrt()) / 2.0
        } else {
            2.0 * c / (-b + (b * b + 4.0 * c).sqrt())
        };
        t / (d[1] + (d[1] * d[1] + t).sqrt())
    };
    delta[0] = -(del + tau);
    delta[1] = -tau;
    work[0] = d[0] + tau + d[1];
    work[1] = 2.0 * d[1] + tau;
    d[1] + tau
}

/// `dlaed6`: one step of the root of a secular equation with three poles. Returns tau and `info`.
fn dlaed6(
    kniter: usize,
    orgati: bool,
    rho: f64,
    d: &[f64; 3],
    z: &[f64; 3],
    finit: f64,
) -> (f64, usize) {
    const MAXIT: usize = 40;
    // `BASE**INT(LOG(SAFMIN) / LOG(BASE) / 3)`: -1022 / 3 truncates to -340.
    const SMALL1: f64 = f64::from_bits((1023 - 340) << 52);
    const SMINV1: f64 = f64::from_bits((1023 + 340) << 52);
    let (mut lbd, mut ubd) = if orgati { (d[1], d[2]) } else { (d[0], d[1]) };
    if finit < 0.0 {
        lbd = 0.0;
    } else {
        ubd = 0.0;
    }
    let mut tau = 0.0;
    if kniter == 2 {
        let (mut a, mut b, mut c);
        if orgati {
            let temp = (d[2] - d[1]) / 2.0;
            c = rho + z[0] / ((d[0] - d[1]) - temp);
            a = c * (d[1] + d[2]) + z[1] + z[2];
            b = c * d[1] * d[2] + z[1] * d[2] + z[2] * d[1];
        } else {
            let temp = (d[0] - d[1]) / 2.0;
            c = rho + z[2] / ((d[2] - d[1]) - temp);
            a = c * (d[0] + d[1]) + z[0] + z[1];
            b = c * d[0] * d[1] + z[0] * d[1] + z[1] * d[0];
        }
        let temp = max3(a.abs(), b.abs(), c.abs());
        a /= temp;
        b /= temp;
        c /= temp;
        tau = if c == 0.0 {
            b / a
        } else if a <= 0.0 {
            (a - (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c)
        } else {
            2.0 * b / (a + (a * a - 4.0 * b * c).abs().sqrt())
        };
        if tau < lbd || tau > ubd {
            tau = (lbd + ubd) / 2.0;
        }
        if d[0] == tau || d[1] == tau || d[2] == tau {
            tau = 0.0;
        } else {
            let temp = finit
                + tau * z[0] / (d[0] * (d[0] - tau))
                + tau * z[1] / (d[1] * (d[1] - tau))
                + tau * z[2] / (d[2] * (d[2] - tau));
            if temp <= 0.0 {
                lbd = tau;
            } else {
                ubd = tau;
            }
            if finit.abs() <= temp.abs() {
                tau = 0.0;
            }
        }
    }
    let small2 = SMALL1 * SMALL1;
    let sminv2 = SMINV1 * SMINV1;
    let temp = if orgati {
        min((d[1] - tau).abs(), (d[2] - tau).abs())
    } else {
        min((d[0] - tau).abs(), (d[1] - tau).abs())
    };
    let mut dscale = *d;
    let mut zscale = *z;
    let mut sclinv = 1.0;
    let scale = temp <= SMALL1;
    if scale {
        let sclfac;
        if temp <= small2 {
            sclfac = sminv2;
            sclinv = small2;
        } else {
            sclfac = SMINV1;
            sclinv = SMALL1;
        }
        for i in 0..3 {
            dscale[i] = d[i] * sclfac;
            zscale[i] = z[i] * sclfac;
        }
        tau *= sclfac;
        lbd *= sclfac;
        ubd *= sclfac;
    }
    let (mut fc, mut df, mut ddf) = (0.0, 0.0, 0.0);
    for i in 0..3 {
        let temp = 1.0 / (dscale[i] - tau);
        let temp1 = zscale[i] * temp;
        let temp2 = temp1 * temp;
        let temp3 = temp2 * temp;
        fc += temp1 / dscale[i];
        df += temp2;
        ddf += temp3;
    }
    let mut f = finit + tau * fc;
    let mut info = 0;
    if f.abs() > 0.0 {
        if f <= 0.0 {
            lbd = tau;
        } else {
            ubd = tau;
        }
        info = 1;
        for _ in 2..=MAXIT {
            let (temp1, temp2) = if orgati {
                (dscale[1] - tau, dscale[2] - tau)
            } else {
                (dscale[0] - tau, dscale[1] - tau)
            };
            let mut a = (temp1 + temp2) * f - temp1 * temp2 * df;
            let mut b = temp1 * temp2 * f;
            let mut c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
            let temp = max3(a.abs(), b.abs(), c.abs());
            a /= temp;
            b /= temp;
            c /= temp;
            let mut eta = if c == 0.0 {
                b / a
            } else if a <= 0.0 {
                (a - (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c)
            } else {
                2.0 * b / (a + (a * a - 4.0 * b * c).abs().sqrt())
            };
            if f * eta >= 0.0 {
                eta = -f / df;
            }
            tau += eta;
            if tau < lbd || tau > ubd {
                tau = (lbd + ubd) / 2.0;
            }
            fc = 0.0;
            let mut erretm = 0.0;
            df = 0.0;
            ddf = 0.0;
            for i in 0..3 {
                let temp = 1.0 / (dscale[i] - tau);
                let temp1 = zscale[i] * temp;
                let temp2 = temp1 * temp;
                let temp3 = temp2 * temp;
                let temp4 = temp1 / dscale[i];
                fc += temp4;
                erretm += temp4.abs();
                df += temp2;
                ddf += temp3;
            }
            f = finit + tau * fc;
            erretm = 8.0 * (finit.abs() + tau.abs() * erretm) + tau.abs() * df;
            if f.abs() <= EPS * erretm {
                info = 0;
                break;
            }
            if f <= 0.0 {
                lbd = tau;
            } else {
                ubd = tau;
            }
        }
    }
    if scale {
        tau *= sclinv;
    }
    (tau, info)
}

/// The sums `dlasd4` takes over the poles left of `ii` (psi) and right of it (phi), with their
/// derivatives and the error bound's running part: `(psi, dpsi, phi, dphi, erretm)`.
fn psi_phi(
    n: usize,
    ii: usize,
    z: &[f64],
    delta: &[f64],
    work: &[f64],
) -> (f64, f64, f64, f64, f64) {
    let (mut psi, mut dpsi, mut erretm) = (0.0, 0.0, 0.0);
    for j in 1..ii {
        let temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
        psi += z[j - 1] * temp;
        dpsi += temp * temp;
        erretm += psi;
    }
    erretm = erretm.abs();
    let (mut phi, mut dphi) = (0.0, 0.0);
    for j in (ii + 1..=n).rev() {
        let temp = z[j - 1] / (work[j - 1] * delta[j - 1]);
        phi += z[j - 1] * temp;
        dphi += temp * temp;
        erretm += phi;
    }
    (psi, dpsi, phi, dphi, erretm)
}

/// `dlasd4`: the `i`th root, increasing, of `1 + rho * sum(z_j^2 / ((d_j - s)(d_j + s)))`, with
/// `delta_j = d_j - sigma` and `work_j = d_j + sigma`. Returns sigma and `info`.
#[allow(clippy::too_many_arguments)]
fn dlasd4(
    n: usize,
    i: usize,
    d: &[f64],
    z: &[f64],
    delta: &mut [f64],
    rho: f64,
    work: &mut [f64],
) -> (f64, usize) {
    const MAXIT: usize = 20;
    macro_rules! d {
        ($i:expr) => {
            d[$i - 1]
        };
    }
    macro_rules! z {
        ($i:expr) => {
            z[$i - 1]
        };
    }
    macro_rules! delta {
        ($i:expr) => {
            delta[$i - 1]
        };
    }
    macro_rules! work {
        ($i:expr) => {
            work[$i - 1]
        };
    }
    if n == 1 {
        delta!(1) = 1.0;
        work!(1) = 1.0;
        return ((d!(1) * d!(1) + rho * z!(1) * z!(1)).sqrt(), 0);
    }
    if n == 2 {
        return (dlasd5(i, d, z, delta, rho, work), 0);
    }
    let eps = EPS;
    let rhoinv = 1.0 / rho;
    let mut sigma;
    if i == n {
        let ii = n - 1;
        let mut niter = 1;
        let temp = rho / 2.0;
        let temp1 = temp / (d!(n) + (d!(n) * d!(n) + temp).sqrt());
        for j in 1..=n {
            work!(j) = d!(j) + d!(n) + temp1;
            delta!(j) = (d!(j) - d!(n)) - temp1;
        }
        let mut psi = 0.0;
        for j in 1..=n - 2 {
            psi += z!(j) * z!(j) / (delta!(j) * work!(j));
        }
        let mut c = rhoinv + psi;
        let mut w =
            c + z!(ii) * z!(ii) / (delta!(ii) * work!(ii)) + z!(n) * z!(n) / (delta!(n) * work!(n));
        let mut tau;
        let solve = |c: f64| {
            let delsq = (d!(n) - d!(n - 1)) * (d!(n) + d!(n - 1));
            let a = -c * delsq + z!(n - 1) * z!(n - 1) + z!(n) * z!(n);
            let b = z!(n) * z!(n) * delsq;
            if a < 0.0 {
                2.0 * b / ((a * a + 4.0 * b * c).sqrt() - a)
            } else {
                (a + (a * a + 4.0 * b * c).sqrt()) / (2.0 * c)
            }
        };
        if w <= 0.0 {
            let temp1 = (d!(n) * d!(n) + rho).sqrt();
            let temp = z!(n - 1) * z!(n - 1)
                / ((d!(n - 1) + temp1) * (d!(n) - d!(n - 1) + rho / (d!(n) + temp1)))
                + z!(n) * z!(n) / rho;
            tau = if c <= temp { rho } else { solve(c) };
        } else {
            tau = solve(c);
        }
        let mut eta = tau / (d!(n) + (d!(n) * d!(n) + tau).sqrt());
        sigma = d!(n) + eta;
        for j in 1..=n {
            delta!(j) = (d!(j) - d!(i)) - eta;
            work!(j) = d!(j) + d!(i) + eta;
        }
        let sums = |delta: &[f64], work: &[f64], tau: f64| {
            let (psi, dpsi, _, _, erretm) = psi_phi(n, n, z, delta, work);
            let temp = z!(n) / (delta[n - 1] * work[n - 1]);
            let phi = z!(n) * temp;
            let dphi = temp * temp;
            let erretm = 8.0 * (-phi - psi) + erretm - phi + rhoinv + tau.abs() * (dpsi + dphi);
            (psi, dpsi, phi, dphi, erretm)
        };
        // The first step's loop runs to II, the others' too: `psi_phi` over all but the last pole.
        let (mut psi, mut dpsi, mut phi, mut dphi, mut erretm) = sums(delta, work, tau);
        w = rhoinv + phi + psi;
        if w.abs() <= eps * erretm {
            return (sigma, 0);
        }
        niter += 1;
        let dtnsq1 = work!(n - 1) * delta!(n - 1);
        let dtnsq = work!(n) * delta!(n);
        c = w - dtnsq1 * dpsi - dtnsq * dphi;
        let a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
        let b = dtnsq * dtnsq1 * w;
        if c < 0.0 {
            c = c.abs();
        }
        if c == 0.0 {
            eta = rho - sigma * sigma;
        } else if a >= 0.0 {
            eta = (a + (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c);
        } else {
            eta = 2.0 * b / (a - (a * a - 4.0 * b * c).abs().sqrt());
        }
        if w * eta > 0.0 {
            eta = -w / (dpsi + dphi);
        }
        let temp = eta - dtnsq;
        if temp > rho {
            eta = rho + dtnsq;
        }
        tau += eta;
        eta /= sigma + (eta + sigma * sigma).sqrt();
        for j in 1..=n {
            delta!(j) -= eta;
            work!(j) += eta;
        }
        sigma += eta;
        (psi, dpsi, phi, dphi, erretm) = sums(delta, work, tau);
        w = rhoinv + phi + psi;
        for _ in niter + 1..=MAXIT {
            if w.abs() <= eps * erretm {
                return (sigma, 0);
            }
            let dtnsq1 = work!(n - 1) * delta!(n - 1);
            let dtnsq = work!(n) * delta!(n);
            let c = w - dtnsq1 * dpsi - dtnsq * dphi;
            let a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
            let b = dtnsq1 * dtnsq * w;
            let mut eta = if a >= 0.0 {
                (a + (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c)
            } else {
                2.0 * b / (a - (a * a - 4.0 * b * c).abs().sqrt())
            };
            if w * eta > 0.0 {
                eta = -w / (dpsi + dphi);
            }
            let temp = eta - dtnsq;
            if temp <= 0.0 {
                eta /= 2.0;
            }
            tau += eta;
            eta /= sigma + (eta + sigma * sigma).sqrt();
            for j in 1..=n {
                delta!(j) -= eta;
                work!(j) += eta;
            }
            sigma += eta;
            (psi, dpsi, phi, dphi, erretm) = sums(delta, work, tau);
            w = rhoinv + phi + psi;
        }
        return (sigma, 1);
    }

    // The root between d(i) and d(i + 1).
    let mut niter = 1;
    let ip1 = i + 1;
    let delsq = (d!(ip1) - d!(i)) * (d!(ip1) + d!(i));
    let delsq2 = delsq / 2.0;
    let temp = delsq2 / (d!(i) + (d!(i) * d!(i) + delsq2).sqrt());
    for j in 1..=n {
        work!(j) = d!(j) + d!(i) + temp;
        delta!(j) = (d!(j) - d!(i)) - temp;
    }
    let mut psi = 0.0;
    for j in 1..i {
        psi += z!(j) * z!(j) / (work!(j) * delta!(j));
    }
    let mut phi = 0.0;
    for j in (i + 2..=n).rev() {
        phi += z!(j) * z!(j) / (work!(j) * delta!(j));
    }
    let c = rhoinv + psi + phi;
    let mut w =
        c + z!(i) * z!(i) / (work!(i) * delta!(i)) + z!(ip1) * z!(ip1) / (work!(ip1) * delta!(ip1));
    let orgati;
    let (mut sg2lb, mut sg2ub);
    let mut tau;
    let mut eta;
    if w > 0.0 {
        orgati = true;
        sg2lb = 0.0;
        sg2ub = delsq2;
        let a = c * delsq + z!(i) * z!(i) + z!(ip1) * z!(ip1);
        let b = z!(i) * z!(i) * delsq;
        tau = if a > 0.0 {
            2.0 * b / (a + (a * a - 4.0 * b * c).abs().sqrt())
        } else {
            (a - (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c)
        };
        eta = tau / (d!(i) + (d!(i) * d!(i) + tau).sqrt());
    } else {
        orgati = false;
        sg2lb = -delsq2;
        sg2ub = 0.0;
        let a = c * delsq - z!(i) * z!(i) - z!(ip1) * z!(ip1);
        let b = z!(ip1) * z!(ip1) * delsq;
        tau = if a < 0.0 {
            2.0 * b / (a - (a * a + 4.0 * b * c).abs().sqrt())
        } else {
            -(a + (a * a + 4.0 * b * c).abs().sqrt()) / (2.0 * c)
        };
        eta = tau / (d!(ip1) + (d!(ip1) * d!(ip1) + tau).abs().sqrt());
    }
    let ii;
    if orgati {
        ii = i;
        sigma = d!(i) + eta;
        for j in 1..=n {
            work!(j) = d!(j) + d!(i) + eta;
            delta!(j) = (d!(j) - d!(i)) - eta;
        }
    } else {
        ii = i + 1;
        sigma = d!(ip1) + eta;
        for j in 1..=n {
            work!(j) = d!(j) + d!(ip1) + eta;
            delta!(j) = (d!(j) - d!(ip1)) - eta;
        }
    }
    let iim1 = ii - 1;
    let iip1 = ii + 1;
    let (mut psi, mut dpsi, mut phi, mut dphi, mut erretm) = psi_phi(n, ii, z, delta, work);
    w = rhoinv + phi + psi;
    let mut swtch3 = if orgati { w < 0.0 } else { w > 0.0 };
    if ii == 1 || ii == n {
        swtch3 = false;
    }
    let temp = z!(ii) / (work!(ii) * delta!(ii));
    let mut dw = dpsi + dphi + temp * temp;
    let temp = z!(ii) * temp;
    w += temp;
    erretm = 8.0 * (phi - psi) + erretm + 2.0 * rhoinv + 3.0 * temp.abs() + tau.abs() * dw;
    if w.abs() <= eps * erretm {
        return (sigma, 0);
    }
    if w <= 0.0 {
        sg2lb = max(sg2lb, tau);
    } else {
        sg2ub = min(sg2ub, tau);
    }
    niter += 1;
    let mut dd = [0.0; 3];
    let mut zz = [0.0; 3];
    if !swtch3 {
        let dtipsq = work!(ip1) * delta!(ip1);
        let dtisq = work!(i) * delta!(i);
        let c = if orgati {
            let ratio = z!(i) / dtisq;
            w - dtipsq * dw + delsq * (ratio * ratio)
        } else {
            let ratio = z!(ip1) / dtipsq;
            w - dtisq * dw - delsq * (ratio * ratio)
        };
        let mut a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
        let b = dtipsq * dtisq * w;
        if c == 0.0 {
            if a == 0.0 {
                a = if orgati {
                    z!(i) * z!(i) + dtipsq * dtipsq * (dpsi + dphi)
                } else {
                    z!(ip1) * z!(ip1) + dtisq * dtisq * (dpsi + dphi)
                };
            }
            eta = b / a;
        } else if a <= 0.0 {
            eta = (a - (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c);
        } else {
            eta = 2.0 * b / (a + (a * a - 4.0 * b * c).abs().sqrt());
        }
    } else {
        let dtiim = work!(iim1) * delta!(iim1);
        let dtiip = work!(iip1) * delta!(iip1);
        let temp = rhoinv + psi + phi;
        let c;
        if orgati {
            let mut temp1 = z!(iim1) / dtiim;
            temp1 *= temp1;
            c = (temp - dtiip * (dpsi + dphi))
                - (d!(iim1) - d!(iip1)) * (d!(iim1) + d!(iip1)) * temp1;
            zz[0] = z!(iim1) * z!(iim1);
            zz[2] = if dpsi < temp1 {
                dtiip * dtiip * dphi
            } else {
                dtiip * dtiip * ((dpsi - temp1) + dphi)
            };
        } else {
            let mut temp1 = z!(iip1) / dtiip;
            temp1 *= temp1;
            c = (temp - dtiim * (dpsi + dphi))
                - (d!(iip1) - d!(iim1)) * (d!(iim1) + d!(iip1)) * temp1;
            zz[0] = if dphi < temp1 {
                dtiim * dtiim * dpsi
            } else {
                dtiim * dtiim * (dpsi + (dphi - temp1))
            };
            zz[2] = z!(iip1) * z!(iip1);
        }
        zz[1] = z!(ii) * z!(ii);
        dd[0] = dtiim;
        dd[1] = delta!(ii) * work!(ii);
        dd[2] = dtiip;
        let info;
        (eta, info) = dlaed6(niter, orgati, c, &dd, &zz, w);
        if info != 0 {
            return (sigma, info);
        }
    }
    if w * eta >= 0.0 {
        eta = -w / dw;
    }
    let temp = if orgati {
        eta - work!(i) * delta!(i)
    } else {
        eta - work!(ip1) * delta!(ip1)
    };
    if temp > sg2ub || temp < sg2lb {
        eta = if w < 0.0 {
            (sg2ub - tau) / 2.0
        } else {
            (sg2lb - tau) / 2.0
        };
    }
    tau += eta;
    eta /= sigma + (sigma * sigma + eta).sqrt();
    let mut prew = w;
    sigma += eta;
    for j in 1..=n {
        work!(j) += eta;
        delta!(j) -= eta;
    }
    (psi, dpsi, phi, dphi, erretm) = psi_phi(n, ii, z, delta, work);
    let temp = z!(ii) / (work!(ii) * delta!(ii));
    dw = dpsi + dphi + temp * temp;
    let temp = z!(ii) * temp;
    w = rhoinv + phi + psi + temp;
    erretm = 8.0 * (phi - psi) + erretm + 2.0 * rhoinv + 3.0 * temp.abs() + tau.abs() * dw;
    if w <= 0.0 {
        sg2lb = max(sg2lb, tau);
    } else {
        sg2ub = min(sg2ub, tau);
    }
    let mut swtch = if orgati {
        -w > prew.abs() / 10.0
    } else {
        w > prew.abs() / 10.0
    };
    for niter in niter + 1..=MAXIT {
        if w.abs() <= eps * erretm {
            return (sigma, 0);
        }
        if !swtch3 {
            let dtipsq = work!(ip1) * delta!(ip1);
            let dtisq = work!(i) * delta!(i);
            let c = if !swtch {
                if orgati {
                    let ratio = z!(i) / dtisq;
                    w - dtipsq * dw + delsq * (ratio * ratio)
                } else {
                    let ratio = z!(ip1) / dtipsq;
                    w - dtisq * dw - delsq * (ratio * ratio)
                }
            } else {
                let temp = z!(ii) / (work!(ii) * delta!(ii));
                if orgati {
                    dpsi += temp * temp;
                } else {
                    dphi += temp * temp;
                }
                w - dtisq * dpsi - dtipsq * dphi
            };
            let mut a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
            let b = dtipsq * dtisq * w;
            if c == 0.0 {
                if a == 0.0 {
                    a = if !swtch {
                        if orgati {
                            z!(i) * z!(i) + dtipsq * dtipsq * (dpsi + dphi)
                        } else {
                            z!(ip1) * z!(ip1) + dtisq * dtisq * (dpsi + dphi)
                        }
                    } else {
                        dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi
                    };
                }
                eta = b / a;
            } else if a <= 0.0 {
                eta = (a - (a * a - 4.0 * b * c).abs().sqrt()) / (2.0 * c);
            } else {
                eta = 2.0 * b / (a + (a * a - 4.0 * b * c).abs().sqrt());
            }
        } else {
            let dtiim = work!(iim1) * delta!(iim1);
            let dtiip = work!(iip1) * delta!(iip1);
            let temp = rhoinv + psi + phi;
            let c;
            if swtch {
                c = temp - dtiim * dpsi - dtiip * dphi;
                zz[0] = dtiim * dtiim * dpsi;
                zz[2] = dtiip * dtiip * dphi;
            } else if orgati {
                let mut temp1 = z!(iim1) / dtiim;
                temp1 *= temp1;
                let temp2 = (d!(iim1) - d!(iip1)) * (d!(iim1) + d!(iip1)) * temp1;
                c = temp - dtiip * (dpsi + dphi) - temp2;
                zz[0] = z!(iim1) * z!(iim1);
                zz[2] = if dpsi < temp1 {
                    dtiip * dtiip * dphi
                } else {
                    dtiip * dtiip * ((dpsi - temp1) + dphi)
                };
            } else {
                let mut temp1 = z!(iip1) / dtiip;
                temp1 *= temp1;
                let temp2 = (d!(iip1) - d!(iim1)) * (d!(iim1) + d!(iip1)) * temp1;
                c = temp - dtiim * (dpsi + dphi) - temp2;
                zz[0] = if dphi < temp1 {
                    dtiim * dtiim * dpsi
                } else {
                    dtiim * dtiim * (dpsi + (dphi - temp1))
                };
                zz[2] = z!(iip1) * z!(iip1);
            }
            dd[0] = dtiim;
            dd[1] = delta!(ii) * work!(ii);
            dd[2] = dtiip;
            let info;
            (eta, info) = dlaed6(niter, orgati, c, &dd, &zz, w);
            if info != 0 {
                return (sigma, info);
            }
        }
        if w * eta >= 0.0 {
            eta = -w / dw;
        }
        let temp = if orgati {
            eta - work!(i) * delta!(i)
        } else {
            eta - work!(ip1) * delta!(ip1)
        };
        if temp > sg2ub || temp < sg2lb {
            eta = if w < 0.0 {
                (sg2ub - tau) / 2.0
            } else {
                (sg2lb - tau) / 2.0
            };
        }
        tau += eta;
        eta /= sigma + (sigma * sigma + eta).sqrt();
        sigma += eta;
        for j in 1..=n {
            work!(j) += eta;
            delta!(j) -= eta;
        }
        prew = w;
        (psi, dpsi, phi, dphi, erretm) = psi_phi(n, ii, z, delta, work);
        let temp = z!(ii) / (work!(ii) * delta!(ii));
        dw = dpsi + dphi + temp * temp;
        let temp = z!(ii) * temp;
        w = rhoinv + phi + psi + temp;
        erretm = 8.0 * (phi - psi) + erretm + 2.0 * rhoinv + 3.0 * temp.abs() + tau.abs() * dw;
        if w * prew > 0.0 && w.abs() > prew.abs() / 10.0 {
            swtch = !swtch;
        }
        if w <= 0.0 {
            sg2lb = max(sg2lb, tau);
        } else {
            sg2ub = min(sg2ub, tau);
        }
    }
    (sigma, 1)
}

/// `dgemm('N', 'N', m, n, k, alpha, a, b, beta, c)`.
#[allow(clippy::too_many_arguments)]
fn dgemm(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    av: View,
    b: &[f64],
    bv: View,
    beta: f64,
    c: &mut [f64],
    cv: View,
) {
    if m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0) {
        return;
    }
    if alpha == 0.0 {
        for j in 1..=n {
            for i in 1..=m {
                let at = cv.ix(i, j);
                c[at] = if beta == 0.0 { 0.0 } else { beta * c[at] };
            }
        }
        return;
    }
    for j in 1..=n {
        if beta == 0.0 {
            for i in 1..=m {
                c[cv.ix(i, j)] = 0.0;
            }
        } else if beta != 1.0 {
            for i in 1..=m {
                let at = cv.ix(i, j);
                c[at] *= beta;
            }
        }
        for l in 1..=k {
            let x = b[bv.ix(l, j)];
            if x != 0.0 {
                let temp = alpha * x;
                for i in 1..=m {
                    let at = cv.ix(i, j);
                    c[at] += temp * a[av.ix(i, l)];
                }
            }
        }
    }
}

/// What `dlasd2` leaves for `dlasd3`, one-based as the workspace is.
struct Deflated {
    k: usize,
    z: Vec<f64>,
    dsigma: Vec<f64>,
    u2: Vec<f64>,
    vt2: Vec<f64>,
    idxc: Vec<usize>,
    ctot: [usize; 5],
}

/// `dlasd2`: the merged problem sorted, its `z` read from the vectors on the right, and the
/// entries with a negligible `z` or a near-equal neighbour deflated to the end.
#[allow(clippy::too_many_arguments)]
fn dlasd2(
    nl: usize,
    nr: usize,
    sqre: usize,
    d: &mut [f64],
    alpha: f64,
    beta: f64,
    u: &mut [f64],
    uv: View,
    vt: &mut [f64],
    vtv: View,
    idxq: &mut [usize],
) -> Deflated {
    let n = nl + nr + 1;
    let m = n + sqre;
    let (nlp1, nlp2) = (nl + 1, nl + 2);
    let u2v = View::new(n);
    let vt2v = View::new(m);
    let mut z = vec![0.0; m + 1];
    let mut dsigma = vec![0.0; n + 1];
    let mut u2 = vec![0.0; n * n];
    let mut vt2 = vec![0.0; m * m];
    let mut idxp = vec![0usize; n + 1];
    let mut idx = vec![0usize; n + 1];
    let mut idxc = vec![0usize; n + 1];
    let mut coltyp = vec![0usize; n + 1];

    let z1 = alpha * vt[vtv.ix(nlp1, nlp1)];
    z[1] = z1;
    for i in (1..=nl).rev() {
        z[i + 1] = alpha * vt[vtv.ix(i, nlp1)];
        d[i] = d[i - 1];
        idxq[i] = idxq[i - 1] + 1;
    }
    for i in nlp2..=m {
        z[i] = beta * vt[vtv.ix(i, nlp2)];
    }
    for value in &mut coltyp[2..=nlp1] {
        *value = 1;
    }
    for value in &mut coltyp[nlp2..=n] {
        *value = 2;
    }
    for value in &mut idxq[nlp2 - 1..n] {
        *value += nlp1;
    }
    for i in 2..=n {
        let q = idxq[i - 1];
        dsigma[i] = d[q - 1];
        u2[u2v.ix(i, 1)] = z[q];
        idxc[i] = coltyp[q];
    }
    dlamrg(nl, nr, &dsigma[2..], true, true, &mut idx[2..]);
    for i in 2..=n {
        let idxi = 1 + idx[i];
        d[i - 1] = dsigma[idxi];
        z[i] = u2[u2v.ix(idxi, 1)];
        coltyp[i] = idxc[idxi];
    }

    let tol = max(alpha.abs(), beta.abs());
    let tol = 8.0 * EPS * max(d[n - 1].abs(), tol);
    let mut k = 1;
    let mut k2 = n + 1;
    let mut jprev = 0;
    let mut all_deflated = false;
    for j in 2..=n {
        if z[j].abs() <= tol {
            k2 -= 1;
            idxp[k2] = j;
            coltyp[j] = 4;
            if j == n {
                all_deflated = true;
            }
        } else {
            jprev = j;
            break;
        }
    }
    if !all_deflated {
        let mut j = jprev;
        loop {
            j += 1;
            if j > n {
                break;
            }
            if z[j].abs() <= tol {
                k2 -= 1;
                idxp[k2] = j;
                coltyp[j] = 4;
            } else if (d[j - 1] - d[jprev - 1]).abs() <= tol {
                let s = z[jprev];
                let c = z[j];
                let tau = dlapy2(c, s);
                let c = c / tau;
                let s = -s / tau;
                z[j] = tau;
                z[jprev] = 0.0;
                let mut idxjp = idxq[idx[jprev]];
                let mut idxj = idxq[idx[j]];
                if idxjp <= nlp1 {
                    idxjp -= 1;
                }
                if idxj <= nlp1 {
                    idxj -= 1;
                }
                rot_columns(u, uv, idxjp, idxj, n, c, s);
                rot_rows(vt, vtv, idxjp, idxj, m, c, s);
                if coltyp[j] != coltyp[jprev] {
                    coltyp[j] = 3;
                }
                coltyp[jprev] = 4;
                k2 -= 1;
                idxp[k2] = jprev;
                jprev = j;
            } else {
                k += 1;
                u2[u2v.ix(k, 1)] = z[jprev];
                dsigma[k] = d[jprev - 1];
                idxp[k] = jprev;
                jprev = j;
            }
        }
        k += 1;
        u2[u2v.ix(k, 1)] = z[jprev];
        dsigma[k] = d[jprev - 1];
        idxp[k] = jprev;
    }

    let mut ctot = [0usize; 5];
    for &ct in &coltyp[2..=n] {
        ctot[ct] += 1;
    }
    let mut psm = [0usize; 5];
    psm[1] = 2;
    psm[2] = 2 + ctot[1];
    psm[3] = psm[2] + ctot[2];
    psm[4] = psm[3] + ctot[3];
    for j in 2..=n {
        let ct = coltyp[idxp[j]];
        idxc[psm[ct]] = j;
        psm[ct] += 1;
    }
    for j in 2..=n {
        dsigma[j] = d[idxp[j] - 1];
        let mut idxj = idxq[idx[idxp[idxc[j]]]];
        if idxj <= nlp1 {
            idxj -= 1;
        }
        for i in 1..=n {
            u2[u2v.ix(i, j)] = u[uv.ix(i, idxj)];
        }
        for i in 1..=m {
            vt2[vt2v.ix(j, i)] = vt[vtv.ix(idxj, i)];
        }
    }
    dsigma[1] = 0.0;
    let hlftol = tol / 2.0;
    if dsigma[2].abs() <= hlftol {
        dsigma[2] = hlftol;
    }
    let (mut c, mut s) = (1.0, 0.0);
    if m > n {
        z[1] = dlapy2(z1, z[m]);
        if z[1] <= tol {
            z[1] = tol;
        } else {
            c = z1 / z[1];
            s = z[m] / z[1];
        }
    } else if z1.abs() <= tol {
        z[1] = tol;
    } else {
        z[1] = z1;
    }
    for i in 2..=k {
        z[i] = u2[u2v.ix(i, 1)];
    }
    for i in 1..=n {
        u2[u2v.ix(i, 1)] = 0.0;
    }
    u2[u2v.ix(nlp1, 1)] = 1.0;
    if m > n {
        for i in 1..=nlp1 {
            let x = vt[vtv.ix(nlp1, i)];
            vt[vtv.ix(m, i)] = -s * x;
            vt2[vt2v.ix(1, i)] = c * x;
        }
        for i in nlp2..=m {
            let x = vt[vtv.ix(m, i)];
            vt2[vt2v.ix(1, i)] = s * x;
            vt[vtv.ix(m, i)] = c * x;
        }
        for i in 1..=m {
            vt2[vt2v.ix(m, i)] = vt[vtv.ix(m, i)];
        }
    } else {
        for i in 1..=m {
            vt2[vt2v.ix(1, i)] = vt[vtv.ix(nlp1, i)];
        }
    }
    if n > k {
        d[k..n].copy_from_slice(&dsigma[k + 1..=n]);
        for j in k + 1..=n {
            for i in 1..=n {
                u[uv.ix(i, j)] = u2[u2v.ix(i, j)];
            }
        }
        for i in 1..=m {
            for j in k + 1..=n {
                vt[vtv.ix(j, i)] = vt2[vt2v.ix(j, i)];
            }
        }
    }
    Deflated {
        k,
        z,
        dsigma,
        u2,
        vt2,
        idxc,
        ctot,
    }
}

/// `dlasd3`: the roots of the secular equation, and the vectors multiplied back. Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dlasd3(
    nl: usize,
    nr: usize,
    sqre: usize,
    d: &mut [f64],
    u: &mut [f64],
    uv: View,
    vt: &mut [f64],
    vtv: View,
    deflated: Deflated,
) -> usize {
    let Deflated {
        k,
        mut z,
        mut dsigma,
        u2,
        mut vt2,
        idxc,
        ctot,
    } = deflated;
    let n = nl + nr + 1;
    let m = n + sqre;
    let (nlp1, nlp2) = (nl + 1, nl + 2);
    let u2v = View::new(n);
    let vt2v = View::new(m);
    if k == 1 {
        d[0] = z[1].abs();
        for i in 1..=m {
            vt[vtv.ix(1, i)] = vt2[vt2v.ix(1, i)];
        }
        for i in 1..=n {
            u[uv.ix(i, 1)] = if z[1] > 0.0 {
                u2[u2v.ix(i, 1)]
            } else {
                -u2[u2v.ix(i, 1)]
            };
        }
        return 0;
    }
    for value in &mut dsigma[1..=k] {
        // `DLAMC3(DSIGMA(I), DSIGMA(I)) - DSIGMA(I)`.
        *value = (*value + *value) - *value;
    }
    let qv = View::new(k);
    let mut q = vec![0.0; k * k];
    q[..k].copy_from_slice(&z[1..=k]);
    let rho = dnrm2(&z[1..=k]);
    dlascl(rho, 1.0, &mut z[1..=k]);
    let rho = rho * rho;
    for j in 1..=k {
        let (column, row) = (uv.ix(1, j), vtv.ix(1, j));
        let (sigma, info) = dlasd4(
            k,
            j,
            &dsigma[1..=k],
            &z[1..=k],
            &mut u[column..column + k],
            rho,
            &mut vt[row..row + k],
        );
        d[j - 1] = sigma;
        if info != 0 {
            return info;
        }
    }
    for i in 1..=k {
        let mut value = u[uv.ix(i, k)] * vt[vtv.ix(i, k)];
        for j in 1..i {
            value *= u[uv.ix(i, j)] * vt[vtv.ix(i, j)]
                / (dsigma[i] - dsigma[j])
                / (dsigma[i] + dsigma[j]);
        }
        for j in i..k {
            value *= u[uv.ix(i, j)] * vt[vtv.ix(i, j)]
                / (dsigma[i] - dsigma[j + 1])
                / (dsigma[i] + dsigma[j + 1]);
        }
        z[i] = sign(value.abs().sqrt(), q[qv.ix(i, 1)]);
    }
    for i in 1..=k {
        let at = vtv.ix(1, i);
        vt[at] = z[1] / u[uv.ix(1, i)] / vt[at];
        u[uv.ix(1, i)] = -1.0;
        for j in 2..=k {
            let at = vtv.ix(j, i);
            vt[at] = z[j] / u[uv.ix(j, i)] / vt[at];
            u[uv.ix(j, i)] = dsigma[j] * vt[at];
        }
        let column = uv.ix(1, i);
        let temp = dnrm2(&u[column..column + k]);
        q[qv.ix(1, i)] = u[uv.ix(1, i)] / temp;
        for j in 2..=k {
            q[qv.ix(j, i)] = u[uv.ix(idxc[j], i)] / temp;
        }
    }
    if k == 2 {
        dgemm(n, k, k, 1.0, &u2, u2v, &q, qv, 0.0, u, uv);
    } else {
        if ctot[1] > 0 {
            dgemm(
                nl,
                k,
                ctot[1],
                1.0,
                &u2,
                u2v.sub(1, 2),
                &q,
                qv.sub(2, 1),
                0.0,
                u,
                uv,
            );
            if ctot[3] > 0 {
                let ktemp = 2 + ctot[1] + ctot[2];
                dgemm(
                    nl,
                    k,
                    ctot[3],
                    1.0,
                    &u2,
                    u2v.sub(1, ktemp),
                    &q,
                    qv.sub(ktemp, 1),
                    1.0,
                    u,
                    uv,
                );
            }
        } else if ctot[3] > 0 {
            let ktemp = 2 + ctot[1] + ctot[2];
            dgemm(
                nl,
                k,
                ctot[3],
                1.0,
                &u2,
                u2v.sub(1, ktemp),
                &q,
                qv.sub(ktemp, 1),
                0.0,
                u,
                uv,
            );
        } else {
            for j in 1..=k {
                for i in 1..=nl {
                    u[uv.ix(i, j)] = u2[u2v.ix(i, j)];
                }
            }
        }
        for j in 1..=k {
            u[uv.ix(nlp1, j)] = q[qv.ix(1, j)];
        }
        let ktemp = 2 + ctot[1];
        let ctemp = ctot[2] + ctot[3];
        dgemm(
            nr,
            k,
            ctemp,
            1.0,
            &u2,
            u2v.sub(nlp2, ktemp),
            &q,
            qv.sub(ktemp, 1),
            0.0,
            u,
            uv.sub(nlp2, 1),
        );
    }
    for i in 1..=k {
        let column = vtv.ix(1, i);
        let temp = dnrm2(&vt[column..column + k]);
        q[qv.ix(i, 1)] = vt[vtv.ix(1, i)] / temp;
        for j in 2..=k {
            q[qv.ix(i, j)] = vt[vtv.ix(idxc[j], i)] / temp;
        }
    }
    if k == 2 {
        dgemm(k, m, k, 1.0, &q, qv, &vt2, vt2v, 0.0, vt, vtv);
        return 0;
    }
    let ktemp = 1 + ctot[1];
    dgemm(k, nlp1, ktemp, 1.0, &q, qv, &vt2, vt2v, 0.0, vt, vtv);
    let ktemp = 2 + ctot[1] + ctot[2];
    if ktemp <= m {
        dgemm(
            k,
            nlp1,
            ctot[3],
            1.0,
            &q,
            qv.sub(1, ktemp),
            &vt2,
            vt2v.sub(ktemp, 1),
            1.0,
            vt,
            vtv,
        );
    }
    let ktemp = ctot[1] + 1;
    let nrp1 = nr + sqre;
    if ktemp > 1 {
        for i in 1..=k {
            q[qv.ix(i, ktemp)] = q[qv.ix(i, 1)];
        }
        for i in nlp2..=m {
            vt2[vt2v.ix(ktemp, i)] = vt2[vt2v.ix(1, i)];
        }
    }
    let ctemp = 1 + ctot[2] + ctot[3];
    dgemm(
        k,
        nrp1,
        ctemp,
        1.0,
        &q,
        qv.sub(1, ktemp),
        &vt2,
        vt2v.sub(ktemp, nlp2),
        0.0,
        vt,
        vtv.sub(1, nlp2),
    );
    0
}

/// `dlasd1`: two solved halves and the row between them, merged. `idxq` holds, on entry, each
/// half's sorting permutation and, on exit, the merged one's. Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dlasd1(
    nl: usize,
    nr: usize,
    sqre: usize,
    d: &mut [f64],
    alpha: f64,
    beta: f64,
    u: &mut [f64],
    uv: View,
    vt: &mut [f64],
    vtv: View,
    idxq: &mut [usize],
) -> usize {
    let n = nl + nr + 1;
    let mut orgnrm = max(alpha.abs(), beta.abs());
    d[nl] = 0.0;
    for value in &d[..n] {
        if value.abs() > orgnrm {
            orgnrm = value.abs();
        }
    }
    dlascl(orgnrm, 1.0, &mut d[..n]);
    let (alpha, beta) = (alpha / orgnrm, beta / orgnrm);
    let deflated = dlasd2(nl, nr, sqre, d, alpha, beta, u, uv, vt, vtv, idxq);
    let k = deflated.k;
    let info = dlasd3(nl, nr, sqre, d, u, uv, vt, vtv, deflated);
    if info != 0 {
        return info;
    }
    dlascl(1.0, orgnrm, &mut d[..n]);
    dlamrg(k, n - k, &d[..n], true, false, idxq);
    0
}

/// `dlasd0`: divide and conquer over a problem of `n` rows and `n + sqre` columns. Returns `info`.
#[allow(clippy::too_many_arguments)]
fn dlasd0(
    n: usize,
    sqre: usize,
    d: &mut [f64],
    e: &mut [f64],
    u: &mut [f64],
    uv: View,
    vt: &mut [f64],
    vtv: View,
) -> usize {
    let m = n + sqre;
    if n <= SMLSIZ {
        return dlasdq(sqre, n, m, n, d, e, vt, vtv, u, uv);
    }
    let (nlvl, nd, inode, ndiml, ndimr) = dlasdt(n, SMLSIZ);
    // IWORK(IDXQ), one-based.
    let mut idxq = vec![0usize; n + 1];
    for i in nd.div_ceil(2)..=nd {
        let ic = inode[i];
        let nl = ndiml[i];
        let nr = ndimr[i];
        let (nlf, nrf) = (ic - nl, ic + 1);
        let info = dlasdq(
            1,
            nl,
            nl + 1,
            nl,
            &mut d[nlf - 1..],
            &mut e[nlf - 1..],
            vt,
            vtv.sub(nlf, nlf),
            u,
            uv.sub(nlf, nlf),
        );
        if info != 0 {
            return info;
        }
        for j in 1..=nl {
            idxq[nlf - 1 + j] = j;
        }
        let sqrei = if i == nd { sqre } else { 1 };
        let info = dlasdq(
            sqrei,
            nr,
            nr + sqrei,
            nr,
            &mut d[nrf - 1..],
            &mut e[nrf - 1..],
            vt,
            vtv.sub(nrf, nrf),
            u,
            uv.sub(nrf, nrf),
        );
        if info != 0 {
            return info;
        }
        for j in 1..=nr {
            idxq[ic + j] = j;
        }
    }
    for lvl in (1..=nlvl).rev() {
        let (lf, ll) = if lvl == 1 {
            (1, 1)
        } else {
            let lf = 1 << (lvl - 1);
            (lf, 2 * lf - 1)
        };
        for i in lf..=ll {
            let ic = inode[i];
            let nl = ndiml[i];
            let nr = ndimr[i];
            let nlf = ic - nl;
            let sqrei = if sqre == 0 && i == ll { sqre } else { 1 };
            let (alpha, beta) = (d[ic - 1], e[ic - 1]);
            let info = dlasd1(
                nl,
                nr,
                sqrei,
                &mut d[nlf - 1..],
                alpha,
                beta,
                u,
                uv.sub(nlf, nlf),
                vt,
                vtv.sub(nlf, nlf),
                &mut idxq[nlf..],
            );
            if info != 0 {
                return info;
            }
        }
    }
    0
}

/// `dbdsdc('U', 'I', n, d, e, u, n, vt, n, ...)`: the values in `d`, decreasing, `u` and `vt`
/// n by n and column major. Returns `info`.
pub(crate) fn dbdsdc(
    n: usize,
    d: &mut [f64],
    e: &mut [f64],
    u: &mut [f64],
    vt: &mut [f64],
) -> usize {
    let full = View::new(n);
    if n == 0 {
        return 0;
    }
    if n == 1 {
        u[0] = sign(1.0, d[0]);
        vt[0] = 1.0;
        d[0] = d[0].abs();
        return 0;
    }
    for matrix in [&mut *u, &mut *vt] {
        matrix.fill(0.0);
        for i in 0..n {
            matrix[i + i * n] = 1.0;
        }
    }
    if n <= SMLSIZ {
        let info = dlasdq(0, n, n, n, d, e, vt, full, u, full);
        if info != 0 {
            return info;
        }
    } else {
        // `dlanst('M')`.
        let mut orgnrm = d[n - 1].abs();
        for i in 0..n - 1 {
            orgnrm = max(orgnrm, d[i].abs());
            orgnrm = max(orgnrm, e[i].abs());
        }
        if orgnrm == 0.0 {
            return 0;
        }
        dlascl(orgnrm, 1.0, &mut d[..n]);
        dlascl(orgnrm, 1.0, &mut e[..n - 1]);
        for value in d[..n].iter_mut() {
            if value.abs() < EPS {
                *value = sign(EPS, *value);
            }
        }
        let nm1 = n - 1;
        let mut start = 1;
        for i in 1..=nm1 {
            if e[i - 1].abs() < EPS || i == nm1 {
                let nsize;
                if i < nm1 {
                    nsize = i - start + 1;
                } else if e[i - 1].abs() >= EPS {
                    nsize = n - start + 1;
                } else {
                    nsize = i - start + 1;
                    u[full.ix(n, n)] = sign(1.0, d[n - 1]);
                    vt[full.ix(n, n)] = 1.0;
                    d[n - 1] = d[n - 1].abs();
                }
                let info = dlasd0(
                    nsize,
                    0,
                    &mut d[start - 1..],
                    &mut e[start - 1..],
                    u,
                    full.sub(start, start),
                    vt,
                    full.sub(start, start),
                );
                if info != 0 {
                    return info;
                }
                start = i + 1;
            }
        }
        dlascl(1.0, orgnrm, &mut d[..n]);
    }
    // Largest first, one transposition per vector.
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
            swap_columns(u, full, i, kk, n);
            swap_rows(vt, full, i, kk, n);
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `B = U S VT` for an upper bidiagonal B past `SMLSIZ`, through the divide and conquer.
    #[test]
    fn a_bidiagonal_comes_back_from_its_values_and_vectors() {
        for n in [26, 40, 77] {
            let d0: Vec<f64> = (0..n).map(|i| 1.0 + ((i * 37) % 11) as f64 / 7.0).collect();
            let e0: Vec<f64> = (0..n - 1)
                .map(|i| 0.5 - ((i * 13) % 5) as f64 / 9.0)
                .collect();
            let (mut d, mut e) = (d0.clone(), e0.clone());
            let (mut u, mut vt) = (vec![0.0; n * n], vec![0.0; n * n]);
            assert_eq!(dbdsdc(n, &mut d, &mut e, &mut u, &mut vt), 0);
            assert!(d.windows(2).all(|w| w[0] >= w[1]), "decreasing");
            for i in 0..n {
                for j in 0..n {
                    let back: f64 = (0..n).map(|l| u[i + l * n] * d[l] * vt[l + j * n]).sum();
                    let expected = if i == j {
                        d0[i]
                    } else if j == i + 1 {
                        e0[i]
                    } else {
                        0.0
                    };
                    assert!(
                        (back - expected).abs() < 1e-12,
                        "{n} ({i}, {j}): {back} {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_tree_splits_as_dlasdt_does() {
        let (lvl, nd, inode, ndiml, ndimr) = dlasdt(60, SMLSIZ);
        assert_eq!((lvl, nd), (2, 3));
        assert_eq!((inode[1], ndiml[1], ndimr[1]), (31, 30, 29));
        assert_eq!((inode[2], ndiml[2], ndimr[2]), (16, 15, 14));
        assert_eq!((inode[3], ndiml[3], ndimr[3]), (46, 14, 14));
    }
}
