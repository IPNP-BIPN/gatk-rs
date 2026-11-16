//! `segment_gcnv_calls.py`, the Python half of `PostprocessGermlineCNVCalls`: gcnvkernel's
//! `ViterbiSegmentationEngine` over one sample of the shards, as the frozen reference computes it
//! (`docs/gcnv-postprocess.md`).
//!
//! Four compiled pytensor functions do the arithmetic, and what is reproduced is the graph
//! `fast_run` leaves and the C pytensor generates for it, not the model's formulas:
//!
//! * **every `logsumexp` is rewritten twice over**: gcnvkernel's max-shifted form, then
//!   `local_log_sum_exp` around it, so it is `(m2 + log(sum(e)) + m1)` with `m1 = max(x)`,
//!   `m2 = max(x - m1)` and `e = isinf(m2) ? exp(m2) : exp(x - m1 - m2)`, added left to right as
//!   the composite's C does under `-ffp-contract=off`;
//! * **reductions are sequential**: pytensor's `CAReduce` runs along the axis in order;
//! * **the one BLAS call** (`CGemv`, the first interval's prior) has an inner dimension of two
//!   classes, which no order of addition can change;
//! * **`exp` and `log` are glibc's SSE2 builds** ([`crate::glibc_libm`]).
//!
//! The Python around them is numpy and scipy: `scipy.special.logsumexp` with numpy's pairwise sum,
//! `np.logaddexp`, and pandas' `precise_xstrtod` reading the shards' tables.

use crate::glibc_libm::{exp, expm1, log, log1p};

/// pandas 2.1.3's `precise_xstrtod` (`tokenizer.c`), what `read_csv` parses a float with by
/// default: at most 17 digits, leading zeros included, accumulated in a double, then one multiplication or
/// division by a power of ten. Not correctly rounded: on the measured shards it differs from the
/// nearest double on 41 % of `log_q_tau_tk`'s entries. `None` where pandas would not see a number.
pub fn precise_xstrtod(text: &str) -> Option<f64> {
    let bytes = text.as_bytes();
    let mut p = 0;
    while p < bytes.len() && bytes[p].is_ascii_whitespace() {
        p += 1;
    }
    let mut negative = false;
    if p < bytes.len() && (bytes[p] == b'-' || bytes[p] == b'+') {
        negative = bytes[p] == b'-';
        p += 1;
    }
    let mut number = 0.0f64;
    let mut exponent: i32 = 0;
    let mut digits = 0;
    let mut decimals = 0;
    while p < bytes.len() && bytes[p].is_ascii_digit() {
        if digits < 17 {
            number = number * 10.0 + f64::from(bytes[p] - b'0');
            digits += 1;
        } else {
            exponent += 1;
        }
        p += 1;
    }
    if p < bytes.len() && bytes[p] == b'.' {
        p += 1;
        while digits < 17 && p < bytes.len() && bytes[p].is_ascii_digit() {
            number = number * 10.0 + f64::from(bytes[p] - b'0');
            p += 1;
            digits += 1;
            decimals += 1;
        }
        while p < bytes.len() && bytes[p].is_ascii_digit() {
            p += 1;
        }
        exponent -= decimals;
    }
    if digits == 0 {
        return None;
    }
    if negative {
        number = -number;
    }
    if p < bytes.len() && (bytes[p] == b'e' || bytes[p] == b'E') {
        let mut q = p + 1;
        let mut exponent_negative = false;
        if q < bytes.len() && (bytes[q] == b'-' || bytes[q] == b'+') {
            exponent_negative = bytes[q] == b'-';
            q += 1;
        }
        let mut n: i32 = 0;
        let mut exponent_digits = 0;
        while exponent_digits < 9 && q < bytes.len() && bytes[q].is_ascii_digit() {
            n = n * 10 + i32::from(bytes[q] - b'0');
            exponent_digits += 1;
            q += 1;
        }
        if exponent_digits > 0 {
            exponent += if exponent_negative { -n } else { n };
            p = q;
        }
    }
    while p < bytes.len() && bytes[p].is_ascii_whitespace() {
        p += 1;
    }
    if p != bytes.len() {
        return None;
    }
    let power = |e: i32| -> f64 { format!("1e{e}").parse().expect("a power of ten") };
    Some(if exponent > 308 {
        f64::INFINITY.copysign(number)
    } else if exponent > 0 {
        number * power(exponent)
    } else if exponent < -308 {
        if exponent < -616 {
            0.0
        } else {
            number / power(-308 - exponent) / power(308)
        }
    } else {
        number / power(-exponent)
    })
}

/// pytensor's `scalar_maximum` reduced in order, NaN winning.
fn pt_max(values: impl Iterator<Item = f64>) -> f64 {
    let mut best = f64::NEG_INFINITY;
    let mut first = true;
    for value in values {
        if first || value > best || value.is_nan() {
            best = value;
            first = false;
        }
        if best.is_nan() {
            break;
        }
    }
    best
}

/// The rewritten `logsumexp` of `x`, as its three parts: `(m1, m2, sum)`.
fn pt_logsumexp_parts(x: &[f64]) -> (f64, f64, f64) {
    let m1 = pt_max(x.iter().copied());
    let m2 = pt_max(x.iter().map(|v| v - m1));
    let mut sum = 0.0;
    for v in x {
        sum += if m2.is_infinite() {
            exp(m2)
        } else {
            exp((v - m1) - m2)
        };
    }
    (m1, m2, sum)
}

/// `Composite{(i1 + log(i0) + i2)}` over the rewritten `logsumexp`.
fn pt_logsumexp(x: &[f64]) -> f64 {
    let (m1, m2, sum) = pt_logsumexp_parts(x);
    m2 + log(sum) + m1
}

/// One step of a forward recursion: `Composite{(i1 + i2 + log(i0) + i3)}` with the emission first.
fn pt_step(emission: f64, x: &[f64]) -> f64 {
    let (m1, m2, sum) = pt_logsumexp_parts(x);
    emission + m2 + log(sum) + m1
}

/// The HMM's shape: `t` sites, `c` hidden states, matrices in row-major order.
#[derive(Debug, Clone)]
pub struct Hmm {
    pub states: usize,
    /// `log_prior_c`.
    pub log_prior: Vec<f64>,
    /// `log_trans_tcc`, `(t - 1) x c x c`, departure then destination.
    pub log_trans: Vec<f64>,
    /// `log_emission_tc`, `t x c`.
    pub log_emission: Vec<f64>,
}

impl Hmm {
    fn sites(&self) -> usize {
        self.log_emission.len() / self.states
    }
    fn trans(&self, t: usize, a: usize, b: usize) -> f64 {
        self.log_trans[(t * self.states + a) * self.states + b]
    }
    fn emission(&self, t: usize, c: usize) -> f64 {
        self.log_emission[t * self.states + c]
    }
}

/// `get_compiled_copy_number_hmm_specs_pytensor_func` for one contig: the first interval's log
/// prior and the class-averaged, normalised log transition matrices. `pi_kc` is `classes x
/// states`, `log_q_tau_tk` is `t x classes`, `stay` has `t - 1` entries.
pub fn hmm_specs(
    pi_kc: &[Vec<f64>],
    stay: &[f64],
    log_q_tau_tk: &[Vec<f64>],
) -> (Vec<f64>, Vec<f64>) {
    let states = pi_kc[0].len();
    // `CGemv`: log(pi.T) . exp(q_0), two terms.
    let q0: Vec<f64> = log_q_tau_tk[0].iter().map(|&v| exp(v)).collect();
    let mut prior: Vec<f64> = (0..states)
        .map(|c| {
            let mut y = 0.0;
            for (k, pi) in pi_kc.iter().enumerate() {
                y += q0[k] * log(pi[c]);
            }
            y
        })
        .collect();
    let normaliser = pt_logsumexp(&prior);
    for value in &mut prior {
        *value -= normaliser;
    }

    let mut trans = Vec::with_capacity(stay.len() * states * states);
    for (t, &s) in stay.iter().enumerate() {
        let not_stay = 1.0 - s;
        let q: Vec<f64> = log_q_tau_tk[t + 1].iter().map(|&v| exp(v)).collect();
        for a in 0..states {
            let row: Vec<f64> = (0..states)
                .map(|b| {
                    let delta = if a == b { 1.0 } else { 0.0 };
                    let mut sum = 0.0;
                    for (k, pi) in pi_kc.iter().enumerate() {
                        sum += q[k] * log(not_stay * pi[b] + s * delta);
                    }
                    sum
                })
                .collect();
            let normaliser = pt_logsumexp(&row);
            trans.extend(row.iter().map(|v| v - normaliser));
        }
    }
    (prior, trans)
}

/// What the forward-backward function returns.
#[derive(Debug, Clone)]
pub struct ForwardBackward {
    pub log_posterior: Vec<f64>,
    pub log_data_likelihood: f64,
    pub alpha: Vec<f64>,
    pub beta: Vec<f64>,
}

/// `PytensorForwardBackward` with `include_alpha_beta_output` and nothing else.
pub fn forward_backward(hmm: &Hmm) -> ForwardBackward {
    let (n, c) = (hmm.sites(), hmm.states);
    let mut alpha = vec![0.0; n * c];
    for (s, value) in alpha[..c].iter_mut().enumerate() {
        *value = hmm.log_prior[s] + hmm.emission(0, s);
    }
    let mut work = vec![0.0; c];
    for t in 1..n {
        for b in 0..c {
            for a in 0..c {
                work[a] = alpha[(t - 1) * c + a] + hmm.trans(t - 1, a, b);
            }
            alpha[t * c + b] = pt_step(hmm.emission(t, b), &work);
        }
    }
    let mut beta = vec![0.0; n * c];
    for t in (1..n).rev() {
        for a in 0..c {
            for b in 0..c {
                work[b] = (beta[t * c + b] + hmm.emission(t, b)) + hmm.trans(t - 1, a, b);
            }
            beta[(t - 1) * c + a] = pt_logsumexp(&work);
        }
    }
    let mut log_posterior = vec![0.0; n * c];
    let mut likelihood = 0.0;
    for t in 0..n {
        let row: Vec<f64> = (0..c).map(|s| alpha[t * c + s] + beta[t * c + s]).collect();
        likelihood = pt_logsumexp(&row);
        for s in 0..c {
            log_posterior[t * c + s] = row[s] - likelihood;
        }
    }
    ForwardBackward {
        log_posterior,
        log_data_likelihood: likelihood,
        alpha,
        beta,
    }
}

/// numpy's `argmax`: the first maximum, or the first NaN.
fn first_argmax(values: impl Iterator<Item = f64>) -> usize {
    let mut best = (0, f64::NEG_INFINITY);
    for (i, value) in values.enumerate() {
        if value.is_nan() {
            return i;
        }
        if i == 0 || value > best.1 {
            best = (i, value);
        }
    }
    best.0
}

/// `PytensorViterbi`: `MaxAndArgmax` over departures, then the backtrack.
pub fn viterbi(hmm: &Hmm) -> Vec<usize> {
    let (n, c) = (hmm.sites(), hmm.states);
    let mut omega: Vec<f64> = (0..c)
        .map(|s| hmm.emission(0, s) + hmm.log_prior[s])
        .collect();
    let mut psi = vec![0usize; n * c];
    for t in 1..n {
        let mut next = vec![0.0; c];
        for b in 0..c {
            let tau = |a: usize| hmm.trans(t - 1, a, b) + omega[a];
            let best = first_argmax((0..c).map(tau));
            psi[t * c + b] = best;
            next[b] = hmm.emission(t, b) + tau(best);
        }
        omega = next;
    }
    let mut path = vec![first_argmax(omega.iter().copied()); n];
    for t in (1..n).rev() {
        path[t - 1] = psi[t * c + path[t]];
    }
    path
}

/// `numpy.add.reduce` over a contiguous float64 vector: `pairwise_sum` from the identity.
pub fn np_sum(a: &[f64]) -> f64 {
    fn pairwise(a: &[f64]) -> f64 {
        let n = a.len();
        if n < 8 {
            let mut res = 0.0;
            for &x in a {
                res += x;
            }
            res
        } else if n <= 128 {
            let mut r = [0.0; 8];
            r.copy_from_slice(&a[..8]);
            let mut i = 8;
            while i < n - (n % 8) {
                for j in 0..8 {
                    r[j] += a[i + j];
                }
                i += 8;
            }
            let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
            while i < n {
                res += a[i];
                i += 1;
            }
            res
        } else {
            let mut half = n / 2;
            half -= half % 8;
            pairwise(&a[..half]) + pairwise(&a[half..])
        }
    }
    0.0 + pairwise(a)
}

/// `scipy.special.logsumexp` (1.11.4) of a vector.
pub fn scipy_logsumexp(a: &[f64]) -> f64 {
    let mut a_max = f64::NEG_INFINITY;
    for &x in a {
        if x > a_max || x.is_nan() {
            a_max = x;
        }
        if a_max.is_nan() {
            break;
        }
    }
    if !a_max.is_finite() {
        a_max = 0.0;
    }
    let tmp: Vec<f64> = a.iter().map(|x| exp(x - a_max)).collect();
    log(np_sum(&tmp)) + a_max
}

/// `npy_logaddexp`.
pub fn np_logaddexp(x: f64, y: f64) -> f64 {
    if x == y {
        return x + std::f64::consts::LN_2;
    }
    let tmp = x - y;
    if tmp > 0.0 {
        x + log1p(exp(-tmp))
    } else if tmp <= 0.0 {
        y + log1p(exp(tmp))
    } else {
        tmp
    }
}

const INV_LN_10_TIMES_10: f64 = 4.342944819032518;
const LN_HALF: f64 = -std::f64::consts::LN_2;

/// Python's `min(a, b)`: `a` unless `b < a`.
fn py_min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}

/// `gcnvkernel.utils.math.logp_complement`.
fn logp_complement(logp: f64) -> f64 {
    let capped = py_min(0.0, logp);
    if capped >= LN_HALF {
        log(-expm1(capped))
    } else {
        log1p(-exp(capped))
    }
}

/// `gcnvkernel.utils.math.logp_to_phred`.
pub fn logp_to_phred(logp: f64, complement: bool) -> f64 {
    let max_phred = -INV_LN_10_TIMES_10 * log(f64::MIN_POSITIVE);
    let capped = py_min(0.0, logp);
    let value = -INV_LN_10_TIMES_10
        * if complement {
            logp_complement(capped)
        } else {
            capped
        };
    py_min(value, max_phred)
}

/// `gcnvkernel.utils.math.logsumexp_double_complement` with its default tolerance.
pub fn logsumexp_double_complement(a: &[f64]) -> f64 {
    let rel_tol = 1e-3;
    let mut a: Vec<f64> = a.iter().map(|&x| if x > 0.0 { 0.0 } else { x }).collect();
    if a.len() == 1 {
        return a[0];
    }
    // `np.sort` then reversed; NaN sorts last, so first once reversed.
    a.sort_by(|x, y| {
        x.partial_cmp(y)
            .unwrap_or_else(|| x.is_nan().cmp(&y.is_nan()))
    });
    a.reverse();
    let size = a.len();
    let mut x = a[0];
    for (i, &entry) in a[1..].iter().enumerate() {
        let x_new = np_logaddexp(x, entry + logp_complement(x));
        if (x_new - x).abs() * ((size - i - 1) as f64) < rel_tol * x.abs() {
            return x_new;
        }
        x = x_new;
    }
    x
}

/// `_constrained_path_logp_pytensor_func`: `alpha_first` and `beta_last` restricted to the
/// allowed states, `m` of them; `emission` is `steps x m` and `trans` `steps x m x m`.
pub fn constrained_path_logp(
    alpha_first: &[f64],
    beta_last: &[f64],
    emission: &[f64],
    trans: &[f64],
    log_data_likelihood: f64,
) -> f64 {
    let m = alpha_first.len();
    let steps = emission.len() / m;
    let mut alpha = alpha_first.to_vec();
    let mut work = vec![0.0; m];
    for j in 0..steps {
        let next: Vec<f64> = (0..m)
            .map(|b| {
                for a in 0..m {
                    work[a] = alpha[a] + trans[(j * m + a) * m + b];
                }
                pt_step(emission[j * m + b], &work)
            })
            .collect();
        alpha = next;
    }
    let row: Vec<f64> = (0..m).map(|s| alpha[s] + beta_last[s]).collect();
    pt_logsumexp(&row) - log_data_likelihood
}

/// `HMMSegmentationQualityCalculator` over one contig's tables.
pub struct Qualities<'a> {
    pub hmm: &'a Hmm,
    pub fb: &'a ForwardBackward,
}

impl Qualities<'_> {
    fn others(&self, call: usize) -> Vec<usize> {
        (0..self.hmm.states).filter(|&s| s != call).collect()
    }

    fn at(&self, table: &[f64], t: usize, states: &[usize]) -> Vec<f64> {
        states
            .iter()
            .map(|&s| table[t * self.hmm.states + s])
            .collect()
    }

    /// `get_log_constrained_posterior_prob`.
    fn constrained(&self, start: usize, end: usize, allowed: &[usize]) -> f64 {
        let alpha_first = self.at(&self.fb.alpha, start, allowed);
        let beta_last = self.at(&self.fb.beta, end, allowed);
        if start == end {
            let sum: Vec<f64> = alpha_first
                .iter()
                .zip(&beta_last)
                .map(|(a, b)| a + b)
                .collect();
            return scipy_logsumexp(&sum) - self.fb.log_data_likelihood;
        }
        let mut emission = Vec::new();
        for t in start + 1..=end {
            emission.extend(self.at(&self.hmm.log_emission, t, allowed));
        }
        let mut trans = Vec::new();
        for t in start..end {
            for &a in allowed {
                for &b in allowed {
                    trans.push(self.hmm.trans(t, a, b));
                }
            }
        }
        constrained_path_logp(
            &alpha_first,
            &beta_last,
            &emission,
            &trans,
            self.fb.log_data_likelihood,
        )
    }

    pub fn some_called(&self, start: usize, end: usize, call: usize) -> f64 {
        logp_to_phred(self.constrained(start, end, &self.others(call)), false)
    }

    pub fn all_called(&self, start: usize, end: usize, call: usize) -> f64 {
        let others = self.others(call);
        let uncorrelated: Vec<f64> = (start..=end)
            .map(|t| scipy_logsumexp(&self.at(&self.fb.log_posterior, t, &others)))
            .collect();
        let uncorrelated_quality = logp_to_phred(logsumexp_double_complement(&uncorrelated), false);
        if uncorrelated_quality < 60.0 {
            let exact = logp_to_phred(self.constrained(start, end, &[call]), true);
            if !exact.is_nan() {
                return exact;
            }
        }
        uncorrelated_quality
    }

    pub fn start(&self, start: usize, call: usize) -> f64 {
        let others = self.others(call);
        let logp = if start == 0 {
            scipy_logsumexp(&self.at(&self.fb.log_posterior, 0, &others))
        } else {
            let path = |prev: usize, next: usize| {
                self.fb.alpha[(start - 1) * self.hmm.states + prev]
                    + self.hmm.trans(start - 1, prev, next)
                    + self.hmm.emission(start, next)
                    + self.fb.beta[start * self.hmm.states + next]
            };
            let mut paths = Vec::new();
            for prev in 0..self.hmm.states {
                for &next in &others {
                    paths.push(path(prev, next));
                }
            }
            paths.push(path(call, call));
            scipy_logsumexp(&paths) - self.fb.log_data_likelihood
        };
        logp_to_phred(logp, false)
    }

    pub fn end(&self, end: usize, call: usize) -> f64 {
        let others = self.others(call);
        let sites = self.hmm.sites();
        let logp = if end == sites - 1 {
            scipy_logsumexp(&self.at(&self.fb.log_posterior, sites - 1, &others))
        } else {
            // `beta_tc[end_index + 1, end_state]`: the departure's state, as gcnvkernel indexes it.
            let path = |from: usize, next: usize| {
                self.fb.alpha[end * self.hmm.states + from]
                    + self.hmm.trans(end, from, next)
                    + self.hmm.emission(end + 1, next)
                    + self.fb.beta[(end + 1) * self.hmm.states + from]
            };
            let mut paths = Vec::new();
            for &from in &others {
                for next in 0..self.hmm.states {
                    paths.push(path(from, next));
                }
            }
            paths.push(path(call, call));
            scipy_logsumexp(&paths) - self.fb.log_data_likelihood
        };
        logp_to_phred(logp, false)
    }
}

/// A row of `copy_number_segments.tsv`, before its qualities are printed.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub num_points: usize,
    pub call: usize,
    pub baseline: usize,
    /// QS, QA, QSS, QSE.
    pub qualities: [f64; 4],
}

/// One sample's inputs as the engine holds them, in the shards' genomic order.
pub struct SampleInputs {
    pub intervals: Vec<(String, i64, i64)>,
    pub log_q_tau_tk: Vec<Vec<f64>>,
    /// `log_c_emission_tc`, one row per interval.
    pub log_emission: Vec<Vec<f64>>,
    /// The sample's ploidy on each contig.
    pub ploidy: Vec<(String, usize)>,
    pub num_states: usize,
    pub p_alt: f64,
    pub cnv_coherence_length: f64,
}

/// `ViterbiSegmentationEngine._viterbi_segments_generator`, without `--clustered-breakpoints`.
pub fn segments(inputs: &SampleInputs) -> Result<Vec<Segment>, String> {
    let mut contigs: Vec<&str> = Vec::new();
    for (contig, _, _) in &inputs.intervals {
        if !contigs.contains(&contig.as_str()) {
            contigs.push(contig);
        }
    }
    let states = inputs.num_states;
    let mut out = Vec::new();
    for contig in contigs {
        let indices: Vec<usize> = (0..inputs.intervals.len())
            .filter(|&i| inputs.intervals[i].0 == contig)
            .collect();
        let baseline = inputs
            .ploidy
            .iter()
            .find(|(name, _)| name == contig)
            .map(|(_, p)| *p)
            .ok_or_else(|| format!("no ploidy for contig {contig}"))?;
        // `get_copy_number_prior_for_sample_jkc`.
        let p_baseline = 1.0 - ((states - 1) as f64) * inputs.p_alt;
        let mut silent = vec![inputs.p_alt; states];
        *silent
            .get_mut(baseline)
            .ok_or_else(|| format!("ploidy {baseline} past the copy-number states"))? = p_baseline;
        let pi_kc = vec![silent, vec![1.0 / states as f64; states]];
        let midpoint = |i: usize| 0.5 * ((inputs.intervals[i].1 + inputs.intervals[i].2) as f64);
        let stay: Vec<f64> = indices
            .windows(2)
            .map(|w| exp(-(midpoint(w[1]) - midpoint(w[0])).abs() / inputs.cnv_coherence_length))
            .collect();
        let log_q_tau: Vec<Vec<f64>> = indices
            .iter()
            .map(|&i| inputs.log_q_tau_tk[i].clone())
            .collect();
        let (log_prior, log_trans) = hmm_specs(&pi_kc, &stay, &log_q_tau);
        let hmm = Hmm {
            states,
            log_prior,
            log_trans,
            log_emission: indices
                .iter()
                .flat_map(|&i| inputs.log_emission[i].iter().copied())
                .collect(),
        };
        let fb = forward_backward(&hmm);
        let path = viterbi(&hmm);
        let quality = Qualities { hmm: &hmm, fb: &fb };
        let mut first = 0;
        while first < path.len() {
            let call = path[first];
            let mut last = first;
            while last + 1 < path.len() && path[last + 1] == call {
                last += 1;
            }
            let qs = quality.some_called(first, last, call);
            let qualities = if last > first {
                [
                    qs,
                    quality.all_called(first, last, call),
                    quality.start(first, call),
                    quality.end(last, call),
                ]
            } else {
                [qs; 4]
            };
            out.push(Segment {
                contig: contig.to_string(),
                start: inputs.intervals[indices[first]].1,
                end: inputs.intervals[indices[last]].2,
                num_points: last - first + 1,
                call,
                baseline,
                qualities,
            });
            first = last + 1;
        }
    }
    Ok(out)
}

/// Python's `'{0:.2f}'.format(q)`: the exact binary value rounded half to even.
pub fn python_fixed2(value: f64) -> String {
    if !value.is_finite() {
        return if value.is_nan() {
            "nan".to_string()
        } else if value > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        };
    }
    // 1100 places hold every double's expansion exactly, so the tie test sees the true digits.
    let exact = format!("{:.1100}", value.abs());
    let (whole, fraction) = exact.split_once('.').expect("a point");
    let mut digits: Vec<u8> = format!("{whole}{}", &fraction[..2]).into_bytes();
    let rest = &fraction[2..];
    let round_up = match rest.as_bytes()[0] {
        b'6'..=b'9' => true,
        b'5' => rest[1..].bytes().any(|d| d != b'0') || (digits[digits.len() - 1] - b'0') % 2 == 1,
        _ => false,
    };
    if round_up {
        let mut i = digits.len();
        loop {
            if i == 0 {
                digits.insert(0, b'1');
                break;
            }
            i -= 1;
            if digits[i] == b'9' {
                digits[i] = b'0';
            } else {
                digits[i] += 1;
                break;
            }
        }
    }
    let text = String::from_utf8(digits).expect("ascii");
    let (int_part, frac_part) = text.split_at(text.len() - 2);
    let sign = if value.is_sign_negative() { "-" } else { "" };
    format!(
        "{sign}{}.{frac_part}",
        if int_part.is_empty() { "0" } else { int_part }
    )
}

/// `io_commons.read_csv`: past the leading `@` lines, a header and tab-separated rows.
fn read_table(path: &str) -> Result<(Vec<String>, Vec<Vec<String>>), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut lines = text.lines().skip_while(|line| line.starts_with('@'));
    let header: Vec<String> = lines
        .next()
        .ok_or_else(|| format!("{path}: no header"))?
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .filter(|line| !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    Ok((header, rows))
}

fn column(header: &[String], name: &str, path: &str) -> Result<usize, String> {
    header
        .iter()
        .position(|h| h == name)
        .ok_or_else(|| format!("{path}: no column {name}"))
}

fn float_rows(path: &str) -> Result<Vec<Vec<f64>>, String> {
    let (_, rows) = read_table(path)?;
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|v| precise_xstrtod(v).ok_or_else(|| format!("{path}: {v} is not a number")))
                .collect()
        })
        .collect()
}

/// The `@RG ... SM:` of a gcnvkernel table.
fn header_sample_name(path: &str) -> Result<String, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    text.lines()
        .take_while(|line| line.starts_with('@'))
        .filter(|line| line.starts_with("@RG"))
        .flat_map(|line| line.split('\t'))
        .find_map(|field| field.strip_prefix("SM:").map(str::to_string))
        .ok_or_else(|| format!("{path}: no sample name"))
}

/// A number of a gcnvkernel config, read as Python's `json` reads it: correctly rounded.
fn json_number(path: &str, key: &str) -> Result<f64, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let after = text
        .split_once(&format!("\"{key}\""))
        .and_then(|(_, rest)| rest.trim_start().strip_prefix(':'))
        .ok_or_else(|| format!("{path}: no {key}"))?
        .trim_start();
    let end = after
        .find(|c: char| !(c.is_ascii_digit() || "+-.eE".contains(c)))
        .unwrap_or(after.len());
    after[..end]
        .parse()
        .map_err(|_| format!("{path}: {key} is not a number"))
}

/// What `ViterbiSegmentationEngine.__init__` and the generator read, for `sample_index`, from
/// shards already in genomic order.
pub fn read_sample_inputs(
    model_shards: &[String],
    calls_shards: &[String],
    ploidy_calls: &str,
    sample_index: usize,
) -> Result<SampleInputs, String> {
    let mut intervals = Vec::new();
    let mut log_q_tau_tk = Vec::new();
    for shard in model_shards {
        let path = format!("{shard}/interval_list.tsv");
        let (header, rows) = read_table(&path)?;
        let (c, s, e) = (
            column(&header, "CONTIG", &path)?,
            column(&header, "START", &path)?,
            column(&header, "END", &path)?,
        );
        for row in rows {
            let number = |i: usize| -> Result<i64, String> {
                row[i]
                    .parse()
                    .map_err(|_| format!("{path}: {} is not an integer", row[i]))
            };
            intervals.push((row[c].clone(), number(s)?, number(e)?));
        }
        log_q_tau_tk.extend(float_rows(&format!("{shard}/log_q_tau_tk.tsv"))?);
    }
    let mut log_emission = Vec::new();
    for shard in calls_shards {
        log_emission.extend(float_rows(&format!(
            "{shard}/SAMPLE_{sample_index}/log_c_emission_tc.tsv"
        ))?);
    }
    let name_path = format!("{}/SAMPLE_{sample_index}/sample_name.txt", calls_shards[0]);
    let sample_name = std::fs::read_to_string(&name_path)
        .map_err(|e| format!("{name_path}: {e}"))?
        .lines()
        .next()
        .unwrap_or_default()
        .trim()
        .to_string();
    let mut ploidy = None;
    let mut entries: Vec<_> = std::fs::read_dir(ploidy_calls)
        .map_err(|e| format!("{ploidy_calls}: {e}"))?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.join("contig_ploidy.tsv").exists())
        .collect();
    entries.sort();
    for dir in entries {
        let path = dir.join("contig_ploidy.tsv").to_string_lossy().into_owned();
        if header_sample_name(&path)? != sample_name {
            continue;
        }
        let (header, rows) = read_table(&path)?;
        let (c, p) = (
            column(&header, "CONTIG", &path)?,
            column(&header, "PLOIDY", &path)?,
        );
        ploidy = Some(
            rows.iter()
                .map(|row| {
                    Ok((
                        row[c].clone(),
                        row[p].parse().map_err(|_| format!("{path}: ploidy"))?,
                    ))
                })
                .collect::<Result<Vec<_>, String>>()?,
        );
    }
    let config = format!("{}/calling_config.json", model_shards[0]);
    Ok(SampleInputs {
        intervals,
        log_q_tau_tk,
        log_emission,
        ploidy: ploidy.ok_or_else(|| format!("no ploidy calls for {sample_name}"))?,
        num_states: json_number(&config, "num_copy_number_states")? as usize,
        p_alt: json_number(&config, "p_alt")?,
        cnv_coherence_length: json_number(&config, "cnv_coherence_length")?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pandas_parses_seventeen_digits_and_scales_once() {
        assert_eq!(
            precise_xstrtod("-1.882712420045912e-06"),
            Some(-1.882712420045912e-06)
        );
        assert_eq!(precise_xstrtod("2"), Some(2.0));
        assert_eq!(precise_xstrtod("x"), None);
        // Leading zeros count among the 17 digits, so this one keeps 13 and lands 1190 ULPs off.
        assert_eq!(
            precise_xstrtod("-0.00047008861787556453").map(f64::to_bits),
            Some((-0.00047008861787556453f64).to_bits() - 1190)
        );
    }

    #[test]
    fn python_rounds_the_exact_value() {
        assert_eq!(python_fixed2(0.125), "0.12");
        assert_eq!(python_fixed2(0.375), "0.38");
        assert_eq!(python_fixed2(73.575), "73.58");
        assert_eq!(python_fixed2(99.999), "100.00");
        assert_eq!(python_fixed2(3076.0), "3076.00");
        assert_eq!(python_fixed2(-0.001), "-0.00");
        assert_eq!(python_fixed2(2.675), "2.67");
    }

    #[test]
    fn numpy_sums_in_eights() {
        let a: Vec<f64> = (0..31).map(|i| 1.0 / f64::from(i + 1)).collect();
        let mut r = [0.0; 8];
        r.copy_from_slice(&a[..8]);
        for i in (8..24).step_by(8) {
            for j in 0..8 {
                r[j] += a[i + j];
            }
        }
        let mut expected = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for x in &a[24..] {
            expected += x;
        }
        assert_eq!(np_sum(&a), expected);
    }
}
