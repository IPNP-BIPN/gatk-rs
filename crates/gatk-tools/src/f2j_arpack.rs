//! `EigenValueDecomposition.symmetricEigs` as Spark 3.5.0 runs it in the oracle image.
//!
//! `RowMatrix.computeSVD` answers a request of `k` eigensamples from ARPACK whenever `k < n / 3`
//! under a hundred samples (LocalARPACK) or `k <= n / 2` from a hundred (DistARPACK). Both drive
//! `dsaupd` and `dseupd` by reverse communication, with `which = "LM"`, `ncv = min(2k, n)`,
//! `tol = 1e-10` and `max(300, 3k)` iterations, and in the image the solver is
//! `dev.ludovic.netlib.arpack.F2jARPACK`: no native library loads, so the arithmetic is
//! `org.netlib.arpack` and the F2J LAPACK and BLAS it calls, in a strict JVM. The trailing
//! eigenvalues of a rank-deficient Gramian are rounding noise, and `rCond` compares that noise, so
//! only the same operations in the same order reach the same rank.
//!
//! Reverse communication is a call here: where ARPACK would return `ido = 1` the operator is
//! applied to the same slice of `workd`, which is all Spark's loop does. What `dsaupd` reaches with
//! `bmat = "I"`, `mode = 1` and exact shifts is ported, nothing else: no `B` products, no user
//! shifts, no `BE`. `dseupd` keeps what the values and the Ritz vectors need and drops the Ritz
//! estimates it computes after them, which Spark never reads.
//!
//! `dgetv0` draws its starting vectors from `dlarnv`, whose seed F2J keeps in a static field: it
//! starts at `1, 3, 5, 7` once per JVM and carries over from call to call. The tool decomposes one
//! panel per JVM, so every call here starts from that seed.
//!
//! Ported from ARPACK (`SRC/dsaupd.f`, `dsaup2.f`, `dsaitr.f`, `dgetv0.f`, `dsapps.f`, `dseigt.f`,
//! `dstqrb.f`, `dsconv.f`, `dsgets.f`, `dsortr.f`, `dseupd.f`, `dsesrt.f`, the initial import of
//! arpack-ng, BSD 3-Clause), LAPACK 3.1.1 (`SRC/dsteqr.f`, `dlaev2.f`, `dlasr.f`, `dlartg.f`,
//! `dlascl.f`, `dlanst.f`, `dlapy2.f`, `dlarnv.f`, `dlaruv.f`, `dgeqr2.f`, `dorm2r.f`, `dlarf.f`,
//! `dlarfg.f`) and reference BLAS (`dnrm2.f`, `dgemv.f`, `dger.f`), BSD 3-Clause, with what F2J does
//! to them: `x**2` is `Math.pow(x, 2)`, which HotSpot answers as `x * x`, and ARPACK's REAL
//! constant `0.717` stays a float.

/// `dlamch('E')`: F2J's is the rounding unit.
const EPS: f64 = f64::EPSILON / 2.0;
/// `dlamch('S')`.
const SAFE_MINIMUM: f64 = f64::MIN_POSITIVE;
/// `Math.pow(dlamch('E'), 2.0 / 3.0)`, measured on the oracle.
const EPS23: f64 = f64::from_bits(0x3DB9_65FE_A53D_6E46);
/// `SAFMN2`, as `f2j_gesdd` derives it.
const SAFE_MINIMUM_2: f64 = f64::from_bits((1023 - 484) << 52);
/// ARPACK's re-orthogonalisation test, `0.717`, a REAL literal that F2J keeps as `0.717f`.
const ORTHOGONAL: f64 = 0.717f32 as f64;
/// `computeSVD`'s `tol`.
const TOLERANCE: f64 = 1e-10;

/// The operator Spark hands ARPACK: `y = A x`.
pub type Operator<'a> = dyn FnMut(&[f64], &mut [f64]) + 'a;

/// `symmetricEigs(mul, n, k, 1e-10, max(300, 3k))`: the eigenvalues in Spark's descending order
/// and their vectors, column major, `n` rows each. `None` where Spark throws: ARPACK's `info` is not
/// zero after `dsaupd`, or nothing converged, which leaves `computeSVD` no first value to read.
pub fn symmetric_eigs(n: usize, k: usize, op: &mut Operator) -> Option<(Vec<f64>, Vec<f64>)> {
    assert!(k > 0 && k < n, "symmetricEigs requires 0 < k < n");
    let ncv = (2 * k).min(n);
    let mxiter = 300.max(3 * k);
    let mut arpack = Lanczos {
        n,
        ncv,
        v: vec![0.0; n * ncv],
        h: vec![0.0; 2 * ncv],
        resid: vec![0.0; n],
        workd: vec![0.0; 3 * n],
        rnorm: 0.0,
        seed: [1, 3, 5, 7],
        op,
    };
    let mut ritz = vec![0.0; ncv];
    let mut bounds = vec![0.0; ncv];
    let mut q = vec![0.0; ncv * ncv];
    let mut workl = vec![0.0; 3 * ncv];
    let nconv = arpack.saup2(k, ncv - k, mxiter, &mut ritz, &mut bounds, &mut q, &mut workl)?;
    if nconv == 0 {
        return None;
    }
    let (d, z) = arpack.seupd(k, nconv, &workl);
    // `sortBy(-_._1)`: stable, under `java.lang.Double.compare`, which has one NaN.
    let key = |value: f64| {
        let negated = -value;
        if negated.is_nan() {
            f64::NAN
        } else {
            negated
        }
    };
    let mut order: Vec<usize> = (0..nconv).collect();
    order.sort_by(|a, b| key(d[*a]).total_cmp(&key(d[*b])));
    let values = order.iter().map(|i| d[*i]).collect();
    let vectors = order
        .iter()
        .flat_map(|i| z[i * n..(i + 1) * n].iter().copied())
        .collect();
    Some((values, vectors))
}

/// `dsaupd`'s state between calls, the parts `bmat = "I"` reads.
struct Lanczos<'a, 'b> {
    n: usize,
    ncv: usize,
    /// The Lanczos basis, `n` by `ncv`, column major.
    v: Vec<f64>,
    /// The tridiagonal `H`, `ldh = ncv`: the subdiagonal in the first column (`h(1,1)` unused),
    /// the diagonal in the second.
    h: Vec<f64>,
    resid: Vec<f64>,
    workd: Vec<f64>,
    rnorm: f64,
    /// `dgetv0`'s `iseed`.
    seed: [i32; 4],
    op: &'a mut Operator<'b>,
}

impl Lanczos<'_, '_> {
    fn sub(&self, j: usize) -> f64 {
        self.h[j - 1]
    }

    fn diag(&self, j: usize) -> f64 {
        self.h[self.ncv + j - 1]
    }

    fn column(&self, j: usize) -> std::ops::Range<usize> {
        (j - 1) * self.n..j * self.n
    }

    /// `dgetv0` with `bmat = "I"` and a random start: the B-norm of the vector, and `ierr`.
    fn getv0(&mut self, j: usize) -> (f64, i32) {
        let n = self.n;
        dlarnv_uniform(&mut self.seed, &mut self.resid);
        self.workd[..n].copy_from_slice(&self.resid);
        let mut rnorm0 = dnrm2(&self.resid);
        if j == 1 {
            return (rnorm0, 0);
        }
        let mut iter = 0;
        loop {
            let (x, y) = self.workd.split_at_mut(n);
            gemv_t(n, j - 1, &self.v, x, &mut y[..j - 1]);
            gemv_n_minus(n, j - 1, &self.v, &y[..j - 1], &mut self.resid);
            self.workd[..n].copy_from_slice(&self.resid);
            let rnorm = dnrm2(&self.resid);
            if rnorm > ORTHOGONAL * rnorm0 {
                return (rnorm, 0);
            }
            iter += 1;
            if iter <= 5 {
                rnorm0 = rnorm;
            } else {
                self.resid.fill(0.0);
                return (0.0, -1);
            }
        }
    }

    /// `dsaitr`: extends the factorization from step `k` by `np` steps. Returns `info`.
    fn saitr(&mut self, k: usize, np: usize) -> usize {
        let n = self.n;
        let (ipj, irj, ivj) = (0, n, 2 * n);
        let mut j = k + 1;
        let mut rstart = false;
        loop {
            if self.rnorm <= 0.0 {
                let mut itry = 1;
                loop {
                    rstart = true;
                    let (rnorm, ierr) = self.getv0(j);
                    self.rnorm = rnorm;
                    if ierr < 0 {
                        itry += 1;
                        if itry <= 3 {
                            continue;
                        }
                        return j - 1;
                    }
                    break;
                }
            }
            let column = self.column(j);
            self.v[column.clone()].copy_from_slice(&self.resid);
            if self.rnorm >= SAFE_MINIMUM {
                let temp1 = 1.0 / self.rnorm;
                for value in &mut self.v[column.clone()] {
                    *value *= temp1;
                }
                for value in &mut self.workd[ipj..ipj + n] {
                    *value *= temp1;
                }
            } else {
                dlascl(self.rnorm, 1.0, &mut self.v[column.clone()]);
                dlascl(self.rnorm, 1.0, &mut self.workd[ipj..ipj + n]);
            }
            self.workd[ivj..ivj + n].copy_from_slice(&self.v[column]);
            {
                let (low, high) = self.workd.split_at_mut(ivj);
                (self.op)(&high[..n], &mut low[irj..irj + n]);
            }
            self.resid.copy_from_slice(&self.workd[irj..irj + n]);
            self.workd[ipj..ipj + n].copy_from_slice(&self.resid);
            let wnorm = dnrm2(&self.resid);
            self.orthogonalize(j);
            self.h[self.ncv + j - 1] = self.workd[irj + j - 1];
            self.h[j - 1] = if j == 1 || rstart { 0.0 } else { self.rnorm };
            self.workd[ipj..ipj + n].copy_from_slice(&self.resid);
            self.rnorm = dnrm2(&self.resid);
            if self.rnorm <= ORTHOGONAL * wnorm {
                let mut iter = 0;
                loop {
                    self.orthogonalize(j);
                    if j == 1 || rstart {
                        self.h[j - 1] = 0.0;
                    }
                    self.h[self.ncv + j - 1] += self.workd[irj + j - 1];
                    self.workd[ipj..ipj + n].copy_from_slice(&self.resid);
                    let rnorm1 = dnrm2(&self.resid);
                    if rnorm1 > ORTHOGONAL * self.rnorm {
                        self.rnorm = rnorm1;
                        break;
                    }
                    self.rnorm = rnorm1;
                    iter += 1;
                    if iter <= 1 {
                        continue;
                    }
                    self.resid.fill(0.0);
                    self.rnorm = 0.0;
                    break;
                }
            }
            rstart = false;
            if self.sub(j) < 0.0 {
                self.h[j - 1] = -self.h[j - 1];
                let target = if j < k + np {
                    &mut self.v[j * n..(j + 1) * n]
                } else {
                    &mut self.resid[..]
                };
                for value in target {
                    *value *= -1.0;
                }
            }
            j += 1;
            if j > k + np {
                return 0;
            }
        }
    }

    /// One classical Gram-Schmidt pass against the first `j` vectors: `workd(irj)` gets the
    /// coefficients of `workd(ipj)`, and `resid` loses its projection.
    fn orthogonalize(&mut self, j: usize) {
        let n = self.n;
        let (x, y) = self.workd.split_at_mut(n);
        gemv_t(n, j, &self.v, x, &mut y[..j]);
        gemv_n_minus(n, j, &self.v, &y[..j], &mut self.resid);
    }

    /// `dsaup2` with exact shifts and `which = "LM"`: the number of converged values, or `None`
    /// where `dsaupd` leaves a non-zero `info` behind.
    #[allow(clippy::too_many_arguments)]
    fn saup2(
        &mut self,
        nev0: usize,
        np0: usize,
        mxiter: usize,
        ritz: &mut [f64],
        bounds: &mut [f64],
        q: &mut [f64],
        workl: &mut [f64],
    ) -> Option<usize> {
        let n = self.n;
        let kplusp = nev0 + np0;
        let (rnorm, _) = self.getv0(1);
        self.rnorm = rnorm;
        if rnorm == 0.0 {
            return None;
        }
        if self.saitr(0, nev0) > 0 {
            return None;
        }
        let (mut nev, mut np) = (nev0, np0);
        let mut iter = 0;
        loop {
            iter += 1;
            if self.saitr(nev, np) > 0 {
                return None;
            }
            if self.seigt(kplusp, ritz, bounds, workl) != 0 {
                return None;
            }
            workl[kplusp..2 * kplusp].copy_from_slice(&ritz[..kplusp]);
            workl[2 * kplusp..3 * kplusp].copy_from_slice(&bounds[..kplusp]);
            nev = nev0;
            np = np0;
            dsgets(true, nev, np, ritz, bounds, workl);
            workl[np..np + nev].copy_from_slice(&bounds[np..np + nev]);
            let nconv = dsconv(&ritz[np..np + nev], &workl[np..np + nev]);
            for j in 0..np {
                if bounds[j] == 0.0 {
                    np -= 1;
                    nev += 1;
                }
            }
            if nconv >= nev0 || iter > mxiter || np == 0 {
                dsortr(Which::SmallestMagnitude, kplusp, ritz, Some(bounds));
                for j in 0..nev0 {
                    bounds[j] /= max(EPS23, ritz[j].abs());
                }
                dsortr(Which::LargestAlgebraic, nev0, bounds, Some(ritz));
                for j in 0..nev0 {
                    bounds[j] *= max(EPS23, ritz[j].abs());
                }
                dsortr(Which::LargestMagnitude, nconv, ritz, Some(bounds));
                self.h[0] = self.rnorm;
                if (iter > mxiter && nconv < nev) || (np == 0 && nconv < nev0) {
                    // `info` 1 or 3, which Spark throws on.
                    return None;
                }
                return Some(nconv);
            } else if nconv < nev {
                let nevbef = nev;
                nev += nconv.min(np / 2);
                if nev == 1 && kplusp >= 6 {
                    nev = kplusp / 2;
                } else if nev == 1 && kplusp > 2 {
                    nev = 2;
                }
                np = kplusp - nev;
                if nevbef < nev {
                    dsgets(true, nev, np, ritz, bounds, workl);
                }
            }
            self.sapps(nev, np, ritz, q);
            self.workd[..n].copy_from_slice(&self.resid);
            self.rnorm = dnrm2(&self.resid);
        }
    }

    /// `dseigt`: the Ritz values of `H` and their error bounds. Returns `ierr`.
    fn seigt(&self, size: usize, eig: &mut [f64], bounds: &mut [f64], workl: &mut [f64]) -> usize {
        eig[..size].copy_from_slice(&self.h[self.ncv..self.ncv + size]);
        workl[..size - 1].copy_from_slice(&self.h[1..size]);
        let ierr = steqr(&mut eig[..size], &mut workl[..size - 1], &mut bounds[..size], 1);
        for bound in &mut bounds[..size] {
            *bound = self.rnorm * bound.abs();
        }
        ierr
    }

    /// `dsapps`: `np` implicit shifts applied to the factorization of length `kev + np`.
    fn sapps(&mut self, kev: usize, np: usize, shift: &[f64], q: &mut [f64]) {
        let n = self.n;
        let ldq = self.ncv;
        let kplusp = kev + np;
        let ncv = self.ncv;
        // `h(i,1)` and `h(i,2)`, one-based.
        let s1 = |i: usize| i - 1;
        let s2 = |i: usize| ncv + i - 1;
        let qa = |r: usize, c: usize| (r - 1) + (c - 1) * ldq;
        for c in 1..=kplusp {
            for r in 1..=kplusp {
                q[qa(r, c)] = if r == c { 1.0 } else { 0.0 };
            }
        }
        if np == 0 {
            return;
        }
        let h = &mut self.h;
        let mut itop = 1;
        'shifts: for jj in 1..=np {
            let mut istart = itop;
            loop {
                let mut iend = kplusp;
                for i in istart..kplusp {
                    let big = h[s2(i)].abs() + h[s2(i + 1)].abs();
                    if h[s1(i + 1)] <= EPS * big {
                        h[s1(i + 1)] = 0.0;
                        iend = i;
                        break;
                    }
                }
                if istart < iend {
                    let f = h[s2(istart)] - shift[jj - 1];
                    let g = h[s1(istart + 1)];
                    let (mut c, mut s, _) = dlartg(f, g);
                    let a1 = c * h[s2(istart)] + s * h[s1(istart + 1)];
                    let a2 = c * h[s1(istart + 1)] + s * h[s2(istart + 1)];
                    let a4 = c * h[s2(istart + 1)] - s * h[s1(istart + 1)];
                    let a3 = c * h[s1(istart + 1)] - s * h[s2(istart)];
                    h[s2(istart)] = c * a1 + s * a2;
                    h[s2(istart + 1)] = c * a4 - s * a3;
                    h[s1(istart + 1)] = c * a3 + s * a4;
                    for j in 1..=(istart + jj).min(kplusp) {
                        let a1 = c * q[qa(j, istart)] + s * q[qa(j, istart + 1)];
                        q[qa(j, istart + 1)] = -s * q[qa(j, istart)] + c * q[qa(j, istart + 1)];
                        q[qa(j, istart)] = a1;
                    }
                    for i in istart + 1..iend {
                        let f = h[s1(i)];
                        let g = s * h[s1(i + 1)];
                        h[s1(i + 1)] *= c;
                        let (c1, s1_, r) = dlartg(f, g);
                        let (r, c1, s1_) = if r < 0.0 { (-r, -c1, -s1_) } else { (r, c1, s1_) };
                        c = c1;
                        s = s1_;
                        h[s1(i)] = r;
                        let a1 = c * h[s2(i)] + s * h[s1(i + 1)];
                        let a2 = c * h[s1(i + 1)] + s * h[s2(i + 1)];
                        let a3 = c * h[s1(i + 1)] - s * h[s2(i)];
                        let a4 = c * h[s2(i + 1)] - s * h[s1(i + 1)];
                        h[s2(i)] = c * a1 + s * a2;
                        h[s2(i + 1)] = c * a4 - s * a3;
                        h[s1(i + 1)] = c * a3 + s * a4;
                        for j in 1..=(i + jj).min(kplusp) {
                            let a1 = c * q[qa(j, i)] + s * q[qa(j, i + 1)];
                            q[qa(j, i + 1)] = -s * q[qa(j, i)] + c * q[qa(j, i + 1)];
                            q[qa(j, i)] = a1;
                        }
                    }
                }
                istart = iend + 1;
                if h[s1(iend)] < 0.0 {
                    h[s1(iend)] = -h[s1(iend)];
                    for r in 1..=kplusp {
                        q[qa(r, iend)] *= -1.0;
                    }
                }
                if iend >= kplusp {
                    break;
                }
            }
            for i in itop..kplusp {
                if h[s1(i + 1)] > 0.0 {
                    continue 'shifts;
                }
                itop += 1;
            }
        }
        for i in itop..kplusp {
            let big = h[s2(i)].abs() + h[s2(i + 1)].abs();
            if h[s1(i + 1)] <= EPS * big {
                h[s1(i + 1)] = 0.0;
            }
        }
        let beta = h[s1(kev + 1)];
        if beta > 0.0 {
            let x = &q[qa(1, kev + 1)..qa(1, kev + 1) + kplusp];
            gemv_n(n, kplusp, &self.v, x, &mut self.workd[n..2 * n]);
        }
        for i in 1..=kev {
            let x = &q[qa(1, kev - i + 1)..qa(1, kev - i + 1) + kplusp - i + 1];
            gemv_n(n, kplusp - i + 1, &self.v, x, &mut self.workd[..n]);
            let target = (kplusp - i) * n;
            self.v[target..target + n].copy_from_slice(&self.workd[..n]);
        }
        for j in 0..kev {
            for i in 0..n {
                self.v[i + j * n] = self.v[i + (np + j) * n];
            }
        }
        if beta > 0.0 {
            self.v[kev * n..(kev + 1) * n].copy_from_slice(&self.workd[n..2 * n]);
        }
        let sigma = q[qa(kplusp, kev)];
        for value in &mut self.resid {
            *value *= sigma;
        }
        if beta > 0.0 {
            daxpy(beta, &self.v[kev * n..(kev + 1) * n], &mut self.resid);
        }
    }

    /// `dseupd(rvec = true, howmny = "A")`: the converged values in ascending order and their
    /// Ritz vectors. `workl` is what `dsaup2` left in `workl(iw)`: the last unsorted Ritz values
    /// and their bounds at `ncv` and `2 ncv`.
    fn seupd(&mut self, nev: usize, nconv: usize, workl: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let (n, ncv) = (self.n, self.ncv);
        let mut d = vec![0.0; nconv];
        let mut irz = workl[ncv..2 * ncv].to_vec();
        let ibd = &workl[2 * ncv..3 * ncv];
        let mut indices: Vec<f64> = (1..=ncv).map(|j| j as f64).collect();
        let mut select = vec![false; ncv];
        dsgets(false, nev, ncv - nev, &mut irz, &mut indices, &mut []);
        let mut numcnv = 0;
        let mut reord = false;
        for j in 1..=ncv {
            let temp1 = max(EPS23, irz[ncv - j].abs());
            let jj = indices[ncv - j] as usize;
            if numcnv < nconv && ibd[jj - 1] <= TOLERANCE * temp1 {
                select[jj - 1] = true;
                numcnv += 1;
                if jj > nev {
                    reord = true;
                }
            }
        }
        // Spark reads neither `info` here nor anything but `d` and `z`: on these refusals they
        // stay what it allocated, zeros and the first columns of `v`.
        let untouched = |v: &[f64]| (vec![0.0; nconv], v[..nconv * n].to_vec());
        if numcnv != nconv {
            return untouched(&self.v);
        }
        let mut ihb = self.h[1..ncv].to_vec();
        let mut ihd = self.h[ncv..2 * ncv].to_vec();
        let mut iq = vec![0.0; ncv * ncv];
        if steqr(&mut ihd, &mut ihb, &mut iq, ncv) != 0 {
            return untouched(&self.v);
        }
        if reord && ncv > 1 {
            let (mut left, mut right) = (1, ncv);
            loop {
                if select[left - 1] {
                    left += 1;
                } else if !select[right - 1] {
                    right -= 1;
                } else {
                    ihd.swap(left - 1, right - 1);
                    for r in 0..ncv {
                        iq.swap(r + ncv * (left - 1), r + ncv * (right - 1));
                    }
                    left += 1;
                    right -= 1;
                }
                if left >= right {
                    break;
                }
            }
        }
        d.copy_from_slice(&ihd[..nconv]);
        dsesrt(nconv, &mut d, ncv, &mut iq);
        let mut tau = vec![0.0; nconv];
        dgeqr2(ncv, nconv, &mut iq, ncv, &mut tau);
        dorm2r_right(n, ncv, nconv, &mut iq, ncv, &tau, &mut self.v);
        let z = self.v[..nconv * n].to_vec();
        (d, z)
    }
}

/// `SIGN(a, b)`, as `f2j_gesdd` has it: a negative zero counts as positive.
fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// Fortran's `MAX`, as F2J writes it.
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

/// `dnrm2`, the scaled sum of squares.
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

/// `daxpy`: `y = y + a x`.
fn daxpy(a: f64, x: &[f64], y: &mut [f64]) {
    if a == 0.0 {
        return;
    }
    for (y, x) in y.iter_mut().zip(x) {
        *y += a * x;
    }
}

/// `dgemv('T', m, columns, 1, a, m, x, 1, 0, y, 1)`.
fn gemv_t(m: usize, columns: usize, a: &[f64], x: &[f64], y: &mut [f64]) {
    for j in 0..columns {
        let mut temp = 0.0;
        for i in 0..m {
            temp += a[i + j * m] * x[i];
        }
        y[j] = 0.0 + temp;
    }
}

/// `dgemv('N', m, columns, -1, a, m, x, 1, 1, y, 1)`.
fn gemv_n_minus(m: usize, columns: usize, a: &[f64], x: &[f64], y: &mut [f64]) {
    for j in 0..columns {
        if x[j] != 0.0 {
            let temp = -x[j];
            for i in 0..m {
                y[i] += temp * a[i + j * m];
            }
        }
    }
}

/// `dgemv('N', m, columns, 1, a, m, x, 1, 0, y, 1)`.
fn gemv_n(m: usize, columns: usize, a: &[f64], x: &[f64], y: &mut [f64]) {
    y[..m].fill(0.0);
    for j in 0..columns {
        if x[j] != 0.0 {
            let temp = x[j];
            for i in 0..m {
                y[i] += temp * a[i + j * m];
            }
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

/// `dlascl('G', 0, 0, cfrom, cto, m, 1, a, ...)` over one column.
fn dlascl(cfrom: f64, cto: f64, a: &mut [f64]) {
    if a.is_empty() || cfrom == 0.0 {
        return;
    }
    let smlnum = SAFE_MINIMUM;
    let bignum = 1.0 / smlnum;
    let (mut cfromc, mut ctoc) = (cfrom, cto);
    loop {
        let cfrom1 = cfromc * smlnum;
        let cto1 = ctoc / bignum;
        let (mul, done) = if cfrom1.abs() > ctoc.abs() && ctoc != 0.0 {
            cfromc = cfrom1;
            (smlnum, false)
        } else if cto1.abs() > cfromc.abs() {
            ctoc = cto1;
            (bignum, false)
        } else {
            (ctoc / cfromc, true)
        };
        for value in a.iter_mut() {
            *value *= mul;
        }
        if done {
            return;
        }
    }
}

/// `dlanst('I', ...)` of a symmetric tridiagonal matrix.
fn dlanst_i(d: &[f64], e: &[f64]) -> f64 {
    let n = d.len();
    if n == 1 {
        return d[0].abs();
    }
    let mut anorm = max(d[0].abs() + e[0].abs(), e[n - 2].abs() + d[n - 1].abs());
    for i in 1..n - 1 {
        anorm = max(anorm, d[i].abs() + e[i].abs() + e[i - 1].abs());
    }
    anorm
}

/// `dlaev2`: `(rt1, rt2, cs1, sn1)`.
fn dlaev2(a: f64, b: f64, c: f64) -> (f64, f64, f64, f64) {
    let sm = a + c;
    let df = a - c;
    let adf = df.abs();
    let tb = b + b;
    let ab = tb.abs();
    let (acmx, acmn) = if a.abs() > c.abs() { (a, c) } else { (c, a) };
    let rt = if adf > ab {
        let ratio = ab / adf;
        adf * (1.0 + ratio * ratio).sqrt()
    } else if adf < ab {
        let ratio = adf / ab;
        ab * (1.0 + ratio * ratio).sqrt()
    } else {
        ab * 2.0f64.sqrt()
    };
    let (rt1, rt2, sgn1) = if sm < 0.0 {
        let rt1 = 0.5 * (sm - rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b, -1)
    } else if sm > 0.0 {
        let rt1 = 0.5 * (sm + rt);
        (rt1, (acmx / rt1) * acmn - (b / rt1) * b, 1)
    } else {
        (0.5 * rt, -0.5 * rt, 1)
    };
    let (cs, sgn2) = if df >= 0.0 { (df + rt, 1) } else { (df - rt, -1) };
    let acs = cs.abs();
    let (mut cs1, mut sn1);
    if acs > ab {
        let ct = -tb / cs;
        sn1 = 1.0 / (1.0 + ct * ct).sqrt();
        cs1 = ct * sn1;
    } else if ab == 0.0 {
        cs1 = 1.0;
        sn1 = 0.0;
    } else {
        let tn = -cs / tb;
        cs1 = 1.0 / (1.0 + tn * tn).sqrt();
        sn1 = tn * cs1;
    }
    if sgn1 == sgn2 {
        let tn = cs1;
        cs1 = -sn1;
        sn1 = tn;
    }
    (rt1, rt2, cs1, sn1)
}

/// `dlasr('R', 'V', direct, rows, columns, c, s, a, rows)`: the plane rotations applied to the
/// columns of `a`, backward or forward.
fn dlasr(backward: bool, rows: usize, columns: usize, c: &[f64], s: &[f64], a: &mut [f64]) {
    let rotate = |j: usize, a: &mut [f64]| {
        let (ctemp, stemp) = (c[j], s[j]);
        if ctemp != 1.0 || stemp != 0.0 {
            for i in 0..rows {
                let temp = a[i + (j + 1) * rows];
                a[i + (j + 1) * rows] = ctemp * temp - stemp * a[i + j * rows];
                a[i + j * rows] = stemp * temp + ctemp * a[i + j * rows];
            }
        }
    };
    if backward {
        for j in (0..columns.saturating_sub(1)).rev() {
            rotate(j, a);
        }
    } else {
        for j in 0..columns.saturating_sub(1) {
            rotate(j, a);
        }
    }
}

/// `dsteqr('I', ...)` when `rows` is `n`, and ARPACK's `dstqrb` when `rows` is 1: the
/// eigenvalues of the tridiagonal `(d, e)` in ascending order, and in `z`, column major and `rows`
/// by `n`, the eigenvectors or only their last components. Returns `info`.
fn steqr(d: &mut [f64], e: &mut [f64], z: &mut [f64], rows: usize) -> usize {
    let n = d.len();
    let arpack = rows == 1;
    if n == 0 {
        return 0;
    }
    if n == 1 {
        z[0] = 1.0;
        return 0;
    }
    let eps = EPS;
    let eps2 = eps * eps;
    let safmin = SAFE_MINIMUM;
    let safmax = 1.0 / safmin;
    let ssfmax = safmax.sqrt() / 3.0;
    let ssfmin = safmin.sqrt() / eps2;
    for c in 0..n {
        for r in 0..rows {
            z[r + c * rows] = if arpack { (c == n - 1) as u8 as f64 } else { (r == c) as u8 as f64 };
        }
    }
    // One-based views: d(i), e(i), work(i).
    let mut work = vec![0.0; 2 * (n - 1)];
    let nmaxit = n * 30;
    let mut jtot = 0;
    let mut l1 = 1;
    let nm1 = n - 1;
    let mut info = 0;
    // A rotation of columns `j` and `j + 1` (one-based) by one pair, `dlasr` over two columns.
    let pair = |z: &mut [f64], j: usize, c: f64, s: f64| {
        if arpack || c != 1.0 || s != 0.0 {
            for i in 0..rows {
                let temp = z[i + j * rows];
                z[i + j * rows] = c * temp - s * z[i + (j - 1) * rows];
                z[i + (j - 1) * rows] = s * temp + c * z[i + (j - 1) * rows];
            }
        }
    };
    loop {
        if l1 > n {
            break;
        }
        if l1 > 1 {
            e[l1 - 2] = 0.0;
        }
        let mut m = n;
        if l1 <= nm1 {
            for mm in l1..=nm1 {
                let tst = e[mm - 1].abs();
                if tst == 0.0 {
                    m = mm;
                    break;
                }
                if tst <= (d[mm - 1].abs().sqrt() * d[mm].abs().sqrt()) * eps {
                    e[mm - 1] = 0.0;
                    m = mm;
                    break;
                }
            }
        }
        let mut l = l1;
        let lsv = l;
        let mut lend = m;
        let lendsv = lend;
        l1 = m + 1;
        if lend == l {
            continue;
        }
        let anorm = dlanst_i(&d[l - 1..lend], &e[l - 1..lend - 1]);
        let mut iscale = 0;
        if anorm == 0.0 {
            continue;
        }
        if anorm > ssfmax {
            iscale = 1;
            dlascl(anorm, ssfmax, &mut d[l - 1..lend]);
            dlascl(anorm, ssfmax, &mut e[l - 1..lend - 1]);
        } else if anorm < ssfmin {
            iscale = 2;
            dlascl(anorm, ssfmin, &mut d[l - 1..lend]);
            dlascl(anorm, ssfmin, &mut e[l - 1..lend - 1]);
        }
        if d[lend - 1].abs() < d[l - 1].abs() {
            lend = lsv;
            l = lendsv;
        }
        if lend > l {
            // QL iteration.
            loop {
                let mut m = lend;
                if l != lend {
                    for mm in l..lend {
                        let tst = e[mm - 1].abs() * e[mm - 1].abs();
                        if tst <= (eps2 * d[mm - 1].abs()) * d[mm].abs() + safmin {
                            m = mm;
                            break;
                        }
                    }
                }
                if m < lend {
                    e[m - 1] = 0.0;
                }
                let mut p = d[l - 1];
                if m == l {
                    d[l - 1] = p;
                    l += 1;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if m == l + 1 {
                    let (rt1, rt2, c, s) = dlaev2(d[l - 1], e[l - 1], d[l]);
                    pair(z, l, c, s);
                    d[l - 1] = rt1;
                    d[l] = rt2;
                    e[l - 1] = 0.0;
                    l += 2;
                    if l <= lend {
                        continue;
                    }
                    break;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let mut g = (d[l] - p) / (2.0 * e[l - 1]);
                let mut r = dlapy2(g, 1.0);
                g = d[m - 1] - p + (e[l - 1] / (g + sign(r, g)));
                let (mut s, mut c) = (1.0, 1.0);
                p = 0.0;
                for i in (l..m).rev() {
                    let f = s * e[i - 1];
                    let b = c * e[i - 1];
                    (c, s, r) = dlartg(g, f);
                    if i != m - 1 {
                        e[i] = r;
                    }
                    g = d[i] - p;
                    r = (d[i - 1] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i] = g + p;
                    g = c * r - b;
                    work[i - 1] = c;
                    work[n - 2 + i] = -s;
                }
                let mm = m - l + 1;
                dlasr(
                    true,
                    rows,
                    mm,
                    &work[l - 1..],
                    &work[n - 2 + l..],
                    &mut z[(l - 1) * rows..],
                );
                d[l - 1] -= p;
                e[l - 1] = g;
            }
        } else {
            // QR iteration.
            loop {
                let mut m = lend;
                if l != lend {
                    for mm in (lend + 1..=l).rev() {
                        let tst = e[mm - 2].abs() * e[mm - 2].abs();
                        if tst <= (eps2 * d[mm - 1].abs()) * d[mm - 2].abs() + safmin {
                            m = mm;
                            break;
                        }
                    }
                }
                if m > lend {
                    e[m - 2] = 0.0;
                }
                let mut p = d[l - 1];
                if m == l {
                    d[l - 1] = p;
                    l -= 1;
                    if l >= lend {
                        continue;
                    }
                    break;
                }
                if m + 1 == l {
                    let (rt1, rt2, c, s) = dlaev2(d[l - 2], e[l - 2], d[l - 1]);
                    pair(z, l - 1, c, s);
                    d[l - 2] = rt1;
                    d[l - 1] = rt2;
                    e[l - 2] = 0.0;
                    if l < lend + 2 {
                        break;
                    }
                    l -= 2;
                    continue;
                }
                if jtot == nmaxit {
                    break;
                }
                jtot += 1;
                let mut g = (d[l - 2] - p) / (2.0 * e[l - 2]);
                let mut r = dlapy2(g, 1.0);
                g = d[m - 1] - p + (e[l - 2] / (g + sign(r, g)));
                let (mut s, mut c) = (1.0, 1.0);
                p = 0.0;
                for i in m..l {
                    let f = s * e[i - 1];
                    let b = c * e[i - 1];
                    (c, s, r) = dlartg(g, f);
                    if i != m {
                        e[i - 2] = r;
                    }
                    g = d[i - 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i - 1] = g + p;
                    g = c * r - b;
                    work[i - 1] = c;
                    work[n - 2 + i] = s;
                }
                let mm = l - m + 1;
                dlasr(
                    false,
                    rows,
                    mm,
                    &work[m - 1..],
                    &work[n - 2 + m..],
                    &mut z[(m - 1) * rows..],
                );
                d[l - 1] -= p;
                e[l - 2] = g;
            }
        }
        if iscale == 1 {
            dlascl(ssfmax, anorm, &mut d[lsv - 1..lendsv]);
            dlascl(ssfmax, anorm, &mut e[lsv - 1..lendsv - 1]);
        } else if iscale == 2 {
            dlascl(ssfmin, anorm, &mut d[lsv - 1..lendsv]);
            dlascl(ssfmin, anorm, &mut e[lsv - 1..lendsv - 1]);
        }
        if jtot < nmaxit {
            continue;
        }
        for value in e.iter() {
            if *value != 0.0 {
                info += 1;
            }
        }
        return info;
    }
    for ii in 2..=n {
        let i = ii - 1;
        let mut k = i;
        let mut p = d[i - 1];
        for j in ii..=n {
            if d[j - 1] < p {
                k = j;
                p = d[j - 1];
            }
        }
        if k != i {
            d[k - 1] = d[i - 1];
            d[i - 1] = p;
            for r in 0..rows {
                z.swap(r + (i - 1) * rows, r + (k - 1) * rows);
            }
        }
    }
    0
}

/// `dsortr`'s and `dsesrt`'s orders: the one each sorts TOWARD, the wanted values last.
#[derive(Clone, Copy)]
enum Which {
    SmallestMagnitude,
    LargestAlgebraic,
    LargestMagnitude,
}

impl Which {
    /// Whether `a` at the lower index and `b` above it are out of order.
    fn swaps(self, a: f64, b: f64) -> bool {
        match self {
            Which::SmallestMagnitude => a.abs() < b.abs(),
            Which::LargestAlgebraic => a > b,
            Which::LargestMagnitude => a.abs() > b.abs(),
        }
    }
}

/// `dsortr`: a shell sort of `x1[..n]`, `x2` following.
fn dsortr(which: Which, n: usize, x1: &mut [f64], mut x2: Option<&mut [f64]>) {
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let ju = j as usize;
                if !which.swaps(x1[ju], x1[ju + igap]) {
                    break;
                }
                x1.swap(ju, ju + igap);
                if let Some(x2) = x2.as_deref_mut() {
                    x2.swap(ju, ju + igap);
                }
                j -= igap as isize;
            }
        }
        igap /= 2;
    }
}

/// `dsesrt('LA', .true., n, x, na, a, lda)`: `x` ascending, the columns of `a` following.
fn dsesrt(n: usize, x: &mut [f64], lda: usize, a: &mut [f64]) {
    let mut igap = n / 2;
    while igap != 0 {
        for i in igap..n {
            let mut j = i as isize - igap as isize;
            while j >= 0 {
                let ju = j as usize;
                if !Which::LargestAlgebraic.swaps(x[ju], x[ju + igap]) {
                    break;
                }
                x.swap(ju, ju + igap);
                for r in 0..lda {
                    a.swap(r + ju * lda, r + (ju + igap) * lda);
                }
                j -= igap as isize;
            }
        }
        igap /= 2;
    }
}

/// `dsgets` with `which = "LM"`: the unwanted values first, and with exact shifts the shifts
/// ordered by their bounds.
fn dsgets(
    exact_shifts: bool,
    kev: usize,
    np: usize,
    ritz: &mut [f64],
    bounds: &mut [f64],
    shifts: &mut [f64],
) {
    dsortr(Which::LargestMagnitude, kev + np, ritz, Some(bounds));
    if exact_shifts && np > 0 {
        dsortr(Which::SmallestMagnitude, np, bounds, Some(ritz));
        shifts[..np].copy_from_slice(&ritz[..np]);
    }
}

/// `dsconv`: how many Ritz values have converged.
fn dsconv(ritz: &[f64], bounds: &[f64]) -> usize {
    ritz.iter()
        .zip(bounds)
        .filter(|(value, bound)| **bound <= TOLERANCE * max(EPS23, value.abs()))
        .count()
}

/// `dlarnv(2, iseed, n, x)`: uniform on (-1, 1), 64 at a time.
fn dlarnv_uniform(seed: &mut [i32; 4], x: &mut [f64]) {
    for chunk in x.chunks_mut(64) {
        let u = dlaruv(seed, chunk.len());
        for (x, u) in chunk.iter_mut().zip(u) {
            *x = 2.0 * u - 1.0;
        }
    }
}

/// The first 64 rows of `dlaruv`'s multipliers, all `dlarnv` asks of it.
const MULTIPLIERS: [[i32; 4]; 64] = [
    [494, 322, 2508, 2549], [2637, 789, 3754, 1145], [255, 1440, 1766, 2253], [2008, 752, 3572, 305],
    [1253, 2859, 2893, 3301], [3344, 123, 307, 1065], [4084, 1848, 1297, 3133], [1739, 643, 3966, 2913],
    [3143, 2405, 758, 3285], [3468, 2638, 2598, 1241], [688, 2344, 3406, 1197], [1657, 46, 2922, 3729],
    [1238, 3814, 1038, 2501], [3166, 913, 2934, 1673], [1292, 3649, 2091, 541], [3422, 339, 2451, 2753],
    [1270, 3808, 1580, 949], [2016, 822, 1958, 2361], [154, 2832, 2055, 1165], [2862, 3078, 1507, 4081],
    [697, 3633, 1078, 2725], [1706, 2970, 3273, 3305], [491, 637, 17, 3069], [931, 2249, 854, 3617],
    [1444, 2081, 2916, 3733], [444, 4019, 3971, 409], [3577, 1478, 2889, 2157], [3944, 242, 3831, 1361],
    [2184, 481, 2621, 3973], [1661, 2075, 1541, 1865], [3482, 4058, 893, 2525], [657, 622, 736, 1409],
    [3023, 3376, 3992, 3445], [3618, 812, 787, 3577], [1267, 234, 2125, 77], [1828, 641, 2364, 3761],
    [164, 4005, 2460, 2149], [3798, 1122, 257, 1449], [3087, 3135, 1574, 3005], [2400, 2640, 3912, 225],
    [2870, 2302, 1216, 85], [3876, 40, 3248, 3673], [1905, 1832, 3401, 3117], [1593, 2247, 2124, 3089],
    [1797, 2034, 2762, 1349], [1234, 2637, 149, 2057], [3460, 1287, 2245, 413], [328, 1691, 166, 65],
    [2861, 496, 466, 1845], [1950, 1597, 4018, 697], [617, 2394, 1399, 3085], [2070, 2584, 190, 3441],
    [3331, 1843, 2879, 1573], [769, 336, 153, 3689], [1558, 1472, 2320, 2941], [2412, 2407, 18, 929],
    [2800, 433, 712, 533], [189, 2096, 2159, 2841], [287, 1761, 2318, 4077], [2045, 2810, 2091, 721],
    [1227, 566, 3443, 2821], [2838, 442, 1510, 2249], [209, 41, 449, 2397], [2770, 1238, 1956, 2817],
];

/// `dlaruv`: `n` uniform values on (0, 1), at most 64, and the seed moved on.
fn dlaruv(seed: &mut [i32; 4], n: usize) -> Vec<f64> {
    const IPW2: i32 = 4096;
    let r = 1.0 / f64::from(IPW2);
    let [mut i1, mut i2, mut i3, mut i4] = *seed;
    let (mut it1, mut it2, mut it3, mut it4) = (0, 0, 0, 0);
    let mut x = vec![0.0; n];
    for (i, value) in x.iter_mut().enumerate() {
        let mm = MULTIPLIERS[i];
        loop {
            it4 = i4 * mm[3];
            it3 = it4 / IPW2;
            it4 -= IPW2 * it3;
            it3 += i3 * mm[3] + i4 * mm[2];
            it2 = it3 / IPW2;
            it3 -= IPW2 * it2;
            it2 += i2 * mm[3] + i3 * mm[2] + i4 * mm[1];
            it1 = it2 / IPW2;
            it2 -= IPW2 * it1;
            it1 += i1 * mm[3] + i2 * mm[2] + i3 * mm[1] + i4 * mm[0];
            it1 %= IPW2;
            *value = r
                * (f64::from(it1)
                    + r * (f64::from(it2) + r * (f64::from(it3) + r * f64::from(it4))));
            if *value != 1.0 {
                break;
            }
            i1 += 2;
            i2 += 2;
            i3 += 2;
            i4 += 2;
        }
    }
    *seed = [it1, it2, it3, it4];
    x
}

/// `dlarfg`: `alpha` and `x` become beta and the reflector's tail. Returns tau.
fn dlarfg(alpha: &mut f64, x: &mut [f64]) -> f64 {
    let mut xnorm = dnrm2(x);
    if xnorm == 0.0 {
        return 0.0;
    }
    let mut beta = -sign(dlapy2(*alpha, xnorm), *alpha);
    let safmin = SAFE_MINIMUM / EPS;
    let mut knt = 0;
    if beta.abs() < safmin {
        let rsafmn = 1.0 / safmin;
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
    }
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
}

/// `dgeqr2(m, n, a, lda, tau, work, info)`.
fn dgeqr2(m: usize, n: usize, a: &mut [f64], lda: usize, tau: &mut [f64]) {
    for i in 1..=m.min(n) {
        let ii = (i - 1) + (i - 1) * lda;
        let (head, tail) = a.split_at_mut(ii + 1);
        tau[i - 1] = dlarfg(&mut head[ii], &mut tail[..m - i]);
        if i < n {
            let aii = a[ii];
            a[ii] = 1.0;
            let v = a[ii..ii + m - i + 1].to_vec();
            let tau = tau[i - 1];
            if tau != 0.0 {
                // `dlarf('Left', ...)`: `work = C' v`, then `C = C - tau v work'`.
                let c = |r: usize, col: usize| (i - 1 + r) + (i + col) * lda;
                let mut work = vec![0.0; n - i];
                for (col, w) in work.iter_mut().enumerate() {
                    let mut temp = 0.0;
                    for (r, vr) in v.iter().enumerate() {
                        temp += a[c(r, col)] * vr;
                    }
                    *w = 0.0 + temp;
                }
                for (col, w) in work.iter().enumerate() {
                    if *w != 0.0 {
                        let temp = -tau * w;
                        for (r, vr) in v.iter().enumerate() {
                            a[c(r, col)] += vr * temp;
                        }
                    }
                }
            }
            a[ii] = aii;
        }
    }
}

/// `dorm2r('Right', 'Notranspose', m, n, k, a, lda, tau, c, m, work, info)`: `c` times the `k`
/// reflectors `dgeqr2` left in `a`.
fn dorm2r_right(m: usize, n: usize, k: usize, a: &mut [f64], lda: usize, tau: &[f64], c: &mut [f64]) {
    let mut work = vec![0.0; m];
    for i in 1..=k {
        let ni = n - i + 1;
        let ii = (i - 1) + (i - 1) * lda;
        let aii = a[ii];
        a[ii] = 1.0;
        let v = &a[ii..ii + ni];
        let tau = tau[i - 1];
        if tau != 0.0 {
            // `dlarf('Right', ...)`: `work = C v`, then `C = C - tau work v'`.
            let first = (i - 1) * m;
            gemv_n(m, ni, &c[first..], v, &mut work);
            for (j, vj) in v.iter().enumerate() {
                if *vj != 0.0 {
                    let temp = -tau * vj;
                    for r in 0..m {
                        c[first + r + j * m] += work[r] * temp;
                    }
                }
            }
        }
        a[ii] = aii;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A diagonal matrix's top eigenvalues are its largest entries.
    #[test]
    fn a_diagonal_matrix_gives_its_largest_entries() {
        let n = 40;
        let diagonal: Vec<f64> = (0..n).map(|i| (i * 7 % n) as f64 + 1.0).collect();
        let (values, vectors) = symmetric_eigs(n, 5, &mut |x: &[f64], y: &mut [f64]| {
            for i in 0..n {
                y[i] = diagonal[i] * x[i];
            }
        })
        .expect("converges");
        let expected = [40.0, 39.0, 38.0, 37.0, 36.0];
        for (value, expected) in values.iter().zip(expected) {
            assert!((value - expected).abs() < 1e-8, "{values:?}");
        }
        // Each vector is a unit vector on the entry holding its value.
        for (column, value) in values.iter().enumerate() {
            let at = diagonal.iter().position(|d| d == &value.round()).unwrap();
            assert!((vectors[at + column * n].abs() - 1.0).abs() < 1e-8);
        }
    }

    /// `dlaruv` from LAPACK's own seed: the first values F2J draws.
    #[test]
    fn the_generator_moves_its_seed() {
        let mut seed = [1, 3, 5, 7];
        let first = dlaruv(&mut seed, 2);
        assert!(first.iter().all(|u| *u > 0.0 && *u < 1.0));
        assert_ne!(seed, [1, 3, 5, 7]);
    }
}

