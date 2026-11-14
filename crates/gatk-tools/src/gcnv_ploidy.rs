//! `cohort_determine_ploidy_and_depth.py` (gcnvkernel 0.9, pymc 5.10.1, pytensor 2.18) under the
//! frozen reference of `docs/gcnv-ploidy.md`: the contig-ploidy model's mean-field ADVI with
//! gcnvkernel's Adamax and linear annealing, the log-emission sampler, the caller, and the writers.
//!
//! Each step below is the C pytensor generated for the reference's compiled graphs (dumped by
//! `tools/python-tools/ploidy_trace.sh`), operation for operation: the same folded constants, the
//! same association, sums in index order from zero, glibc's SSE2 `exp`/`log`/`log1p`/`expm1` and
//! `lgamma`, pytensor's own digamma.
//!
//! The C's own spelling is kept where clippy would shorten it (`x * -1.0`, AS 103's truncated
//! Euler constant, fdlibm's long literals, index loops over several arrays).
#![allow(
    clippy::approx_constant,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::neg_multiply
)]

use crate::gcnv_ploidy_io::{json_floats, ndarray_tsv, python_repr};
use crate::gcnv_segments::{np_sum, precise_xstrtod};
use crate::glibc_libm::{exp, expm1, log, log1p};
use crate::numpy_random::{reseed, Pcg64, RandomState};

/// The script's arguments, as Java hands them over (`%e`, read back by Python's `float`).
#[derive(Clone, Debug)]
pub struct Settings {
    pub mapping_error_rate: f64,
    pub psi_s_scale: f64,
    pub mean_bias_sd: f64,
    pub psi_j_scale: f64,
    pub learning_rate: f64,
    pub adamax_beta1: f64,
    pub adamax_beta2: f64,
    pub log_emission_samples_per_round: i64,
    pub log_emission_sampling_rounds: i64,
    pub log_emission_sampling_median_rel_error: f64,
    pub max_advi_iter_first_epoch: i64,
    pub max_advi_iter_subsequent_epochs: i64,
    pub min_training_epochs: i64,
    pub max_training_epochs: i64,
    pub initial_temperature: f64,
    pub num_thermal_advi_iters: i64,
    pub convergence_snr_averaging_window: i64,
    pub convergence_snr_trigger_threshold: f64,
    pub convergence_snr_countdown_window: i64,
    pub max_calling_iters: i64,
    pub caller_update_convergence_threshold: f64,
    pub caller_external_admixing_rate: f64,
    pub disable_annealing: bool,
    pub disable_caller: bool,
    pub disable_sampler: bool,
}

/// What Java hands the script: contigs in interval order, intervals per contig, each sample's
/// counts per contig, and the priors table's text.
pub struct Cohort {
    pub contigs: Vec<String>,
    pub intervals_per_contig: Vec<u64>,
    pub samples: Vec<(String, Vec<i64>)>,
    pub priors_text: String,
}

/// The files the script writes besides the two it copies: `(name, text)` under the model
/// directory, and under the calls directory.
pub struct Files {
    pub model: Vec<(String, String)>,
    pub calls: Vec<(String, String)>,
}

// ---------------------------------------------------------------------------------------------
// Special functions.

/// glibc 2.35's `__ieee754_lgamma_r` (fdlibm's `e_lgamma_r.c`) for positive finite `x`, the only
/// arguments the model gives it, calling the SSE2 `log` the ifunc picks.
pub fn lgamma(x: f64) -> f64 {
    const A: [f64; 12] = [
        7.72156649015328655494e-02,
        3.22467033424113591611e-01,
        6.73523010531292681824e-02,
        2.05808084325167332806e-02,
        7.38555086081402883957e-03,
        2.89051383673415629091e-03,
        1.19270763183362067845e-03,
        5.10069792153511336608e-04,
        2.20862790713908385557e-04,
        1.08011567247583939954e-04,
        2.52144565451257326939e-05,
        4.48640949618915160150e-05,
    ];
    const TC: f64 = 1.46163214496836224576e+00;
    const TF: f64 = -1.21486290535849611461e-01;
    const TT: f64 = -3.63867699703950536541e-18;
    const T: [f64; 15] = [
        4.83836122723810047042e-01,
        -1.47587722994593911752e-01,
        6.46249402391333854778e-02,
        -3.27885410759859649565e-02,
        1.79706750811820387126e-02,
        -1.03142241298341437450e-02,
        6.10053870246291332635e-03,
        -3.68452016781138256760e-03,
        2.25964780900612472250e-03,
        -1.40346469989232843813e-03,
        8.81081882437654011382e-04,
        -5.38595305356740546715e-04,
        3.15632070903625950361e-04,
        -3.12754168375120860518e-04,
        3.35529192635519073543e-04,
    ];
    const U: [f64; 6] = [
        -7.72156649015328655494e-02,
        6.32827064025093366517e-01,
        1.45492250137234768737e+00,
        9.77717527963372745603e-01,
        2.28963728064692451092e-01,
        1.33810918536787660377e-02,
    ];
    const V: [f64; 6] = [
        0.0,
        2.45597793713041134822e+00,
        2.12848976379893395361e+00,
        7.69285150456672783825e-01,
        1.04222645593369134254e-01,
        3.21709242282423911810e-03,
    ];
    const S: [f64; 7] = [
        -7.72156649015328655494e-02,
        2.14982415960608852501e-01,
        3.25778796408930981787e-01,
        1.46350472652464452805e-01,
        2.66422703033638609560e-02,
        1.84028451407337715652e-03,
        3.19475326584100867617e-05,
    ];
    const R: [f64; 7] = [
        0.0,
        1.39200533467621045958e+00,
        7.21935547567138069525e-01,
        1.71933865632803078993e-01,
        1.86459191715652901344e-02,
        7.77942496381893596434e-04,
        7.32668430744625636189e-06,
    ];
    const W: [f64; 7] = [
        4.18938533204672725052e-01,
        8.33333333333329678849e-02,
        -2.77777777728775536470e-03,
        7.93650558643019558500e-04,
        -5.95187557450339963135e-04,
        8.36339918996282139126e-04,
        -1.63092934096575273989e-03,
    ];
    let bits = x.to_bits();
    let hx = (bits >> 32) as u32 as i32;
    let lx = bits as u32;
    let ix = hx & 0x7fff_ffff;
    if ix >= 0x7ff0_0000 {
        return x * x;
    }
    if (ix as u32 | lx) == 0 {
        return 1.0 / x.abs();
    }
    assert!(hx >= 0, "lgamma of a negative argument is not ported");
    if ix < 0x3b90_0000 {
        return -log(x);
    }
    if (((ix - 0x3ff0_0000) as u32 | lx) == 0) || (((ix - 0x4000_0000) as u32 | lx) == 0) {
        return 0.0;
    }
    if ix < 0x4000_0000 {
        let (mut r, y, i);
        if ix <= 0x3fec_cccc {
            r = -log(x);
            if ix >= 0x3FE7_6944 {
                y = 1.0 - x;
                i = 0;
            } else if ix >= 0x3FCD_A661 {
                y = x - (TC - 1.0);
                i = 1;
            } else {
                y = x;
                i = 2;
            }
        } else {
            r = 0.0;
            if ix >= 0x3FFB_B4C3 {
                y = 2.0 - x;
                i = 0;
            } else if ix >= 0x3FF3_B4C4 {
                y = x - TC;
                i = 1;
            } else {
                y = x - 1.0;
                i = 2;
            }
        }
        match i {
            0 => {
                let z = y * y;
                let p1 = A[0] + z * (A[2] + z * (A[4] + z * (A[6] + z * (A[8] + z * A[10]))));
                let p2 = z * (A[1] + z * (A[3] + z * (A[5] + z * (A[7] + z * (A[9] + z * A[11])))));
                let p = y * p1 + p2;
                r += p - 0.5 * y;
            }
            1 => {
                let z = y * y;
                let w = z * y;
                let p1 = T[0] + w * (T[3] + w * (T[6] + w * (T[9] + w * T[12])));
                let p2 = T[1] + w * (T[4] + w * (T[7] + w * (T[10] + w * T[13])));
                let p3 = T[2] + w * (T[5] + w * (T[8] + w * (T[11] + w * T[14])));
                let p = z * p1 - (TT - w * (p2 + y * p3));
                r += TF + p;
            }
            _ => {
                let p1 = y * (U[0] + y * (U[1] + y * (U[2] + y * (U[3] + y * (U[4] + y * U[5])))));
                let p2 = 1.0 + y * (V[1] + y * (V[2] + y * (V[3] + y * (V[4] + y * V[5]))));
                r += -0.5 * y + p1 / p2;
            }
        }
        r
    } else if ix < 0x4020_0000 {
        let i = x as i32;
        let y = x - f64::from(i);
        let p =
            y * (S[0] + y * (S[1] + y * (S[2] + y * (S[3] + y * (S[4] + y * (S[5] + y * S[6]))))));
        let q = 1.0 + y * (R[1] + y * (R[2] + y * (R[3] + y * (R[4] + y * (R[5] + y * R[6])))));
        let mut r = 0.5 * y + p / q;
        let mut z = 1.0;
        for k in (3..=i).rev() {
            z *= y + f64::from(k - 1);
        }
        if i >= 3 {
            r += log(z);
        }
        r
    } else if ix < 0x4390_0000 {
        let t = log(x);
        let z = 1.0 / x;
        let y = z * z;
        let w = W[0] + z * (W[1] + y * (W[2] + y * (W[3] + y * (W[4] + y * (W[5] + y * W[6])))));
        (x - 0.5) * (t - 1.0) + w
    } else {
        x * (log(x) - 1.0)
    }
}

/// pytensor's `_psi` (Bernardo's AS 103), in its C.
pub fn psi(x: f64) -> f64 {
    let mut y = x;
    let mut psi = 0.0;
    if y <= 0.0 {
        return psi;
    }
    if y <= 1.0e-5 {
        return -0.5772156649 - 1.0 / y;
    }
    while y < 8.5 {
        psi -= 1.0 / y;
        y += 1.0;
    }
    let mut r = 1.0 / y;
    psi = psi + log(y) - 0.5 * r;
    r *= r;
    psi - r * (8.333333333e-2 - r * (8.333333333e-3 - r * 3.968253968e-3))
}

/// pytensor's `softplus` C.
fn softplus(x: f64) -> f64 {
    if x < -37.0 {
        exp(x)
    } else if x < 18.0 {
        log1p(exp(x))
    } else if x < 33.3 {
        x + exp(-x)
    } else {
        x
    }
}

/// pytensor's `maximum` C: `a > b ? a : (b >= a ? b : nan)`.
fn maximum(a: f64, b: f64) -> f64 {
    if a > b {
        a
    } else if b >= a {
        b
    } else {
        f64::NAN
    }
}

fn seq_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0, |acc, x| acc + x)
}

/// `np.median` of a vector without NaN: the middle order statistic, or the mean of the two.
fn np_median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let n = values.len();
    if n % 2 == 1 {
        np_sum(&values[n / 2..n / 2 + 1]) / 1.0
    } else {
        np_sum(&values[n / 2 - 1..n / 2 + 1]) / 2.0
    }
}

// ---------------------------------------------------------------------------------------------
// The model.

/// The workspace's shared values and the constants pytensor folded into the graphs.
struct Model {
    s: usize,
    j: usize,
    k: usize,
    n_sj: Vec<f64>,
    n_s: Vec<f64>,
    t_j: Vec<f64>,
    /// `contig_exclusion_mask_jj` as the float64 `dot` reads it.
    mask: Vec<f64>,
    eps_j: Vec<f64>,
    one_minus_error: f64,
    lam_s: f64,
    log_lam_s: f64,
    lam_j: f64,
    log_lam_j: f64,
    log_p_jk: Vec<f64>,
    lgamma_1n: Vec<f64>,
    switch_n: Vec<f64>,
}

/// The approximation, the optimizer's slots, and the discrete sector's shared values.
struct State {
    mu: Vec<f64>,
    rho: Vec<f64>,
    m_mu: Vec<f64>,
    u_mu: Vec<f64>,
    m_rho: Vec<f64>,
    u_rho: Vec<f64>,
    beta1_power: f64,
    temperature: f64,
    log_q: Vec<f64>,
    log_emission: Vec<f64>,
}

impl Model {
    fn d(&self) -> usize {
        2 * self.j + self.s
    }

    /// `mean_ploidy_sj`: `sum_k exp(log_q) * k`.
    fn mean_ploidy(&self, log_q: &[f64]) -> Vec<f64> {
        (0..self.s * self.j)
            .map(|sj| seq_sum((0..self.k).map(|k| exp(log_q[sj * self.k + k]) * k as f64)))
            .collect()
    }

    /// `dot(gamma_sj, mask_jj)` (or its transpose), the inner index ascending from zero.
    fn dot_mask(&self, a: &[f64], transpose: bool) -> Vec<f64> {
        let (s, j) = (self.s, self.j);
        let mut out = vec![0.0; s * j];
        for si in 0..s {
            for l in 0..j {
                let mut acc = 0.0;
                for jj in 0..j {
                    let m = if transpose {
                        self.mask[l * j + jj]
                    } else {
                        self.mask[jj * j + l]
                    };
                    acc += a[si * j + jj] * m;
                }
                out[si * j + l] = acc;
            }
        }
        out
    }

    /// The step function (graph 2): one draw, the loss before normalization, and Adamax's update
    /// of `mu` and `rho` (gcnvkernel's `FancyAdamax`), with the temperature's linear drop.
    fn step(&self, st: &mut State, rng: &mut Pcg64, settings: &Settings, drop: Option<f64>) -> f64 {
        let (s, j, k, d) = (self.s, self.j, self.k, self.d());
        let temperature = st.temperature;
        let e: Vec<f64> = (0..d).map(|_| 0.0 + 1.0 * rng.standard_normal()).collect();
        let b1 = settings.adamax_beta1;
        let b2 = settings.adamax_beta2;
        let beta1_power = b1 * st.beta1_power;
        let a = settings.learning_rate / (1.0 - beta1_power);
        let q: Vec<f64> = st.log_q.iter().map(|&x| exp(x)).collect();
        let mean_ploidy = self.mean_ploidy(&st.log_q);

        // n50 to n54: sigma, the entropy's gradient, log q(z) and z.
        let sigma: Vec<f64> = st.rho.iter().map(|&r| softplus(r)).collect();
        let sigmoid: Vec<f64> = st.rho.iter().map(|&r| 1.0 / (1.0 + exp(-r))).collect();
        let grad_entropy: Vec<f64> = (0..d)
            .map(|i| (-temperature / sigma[i]) * sigmoid[i])
            .collect();
        let log_q_z =
            seq_sum((0..d).map(|i| ((-0.5 * (e[i] * e[i])) - 0.9189385332046727) - log(sigma[i])));
        let z: Vec<f64> = (0..d).map(|i| sigma[i] * e[i] + st.mu[i]).collect();

        // n64 to n66: the priors' logp and gradients.
        let exponential = |x: f64, lam: f64, log_lam: f64| {
            let ex = exp(x);
            let positive = ex >= 0.0;
            let grad = if positive { (-lam * ex) * -1.0 } else { 0.0 };
            let logp = (if positive {
                log_lam - lam * ex
            } else {
                f64::NEG_INFINITY
            }) + x;
            (ex, grad, logp)
        };
        let psi_s: Vec<(f64, f64, f64)> = (0..s)
            .map(|i| exponential(z[2 * j + i], self.lam_s, self.log_lam_s))
            .collect();
        let psi_j: Vec<(f64, f64, f64)> = (0..j)
            .map(|i| exponential(z[j + i], self.lam_j, self.log_lam_j))
            .collect();
        let bias: Vec<(f64, f64, f64)> = (0..j)
            .map(|i| {
                let x = z[i];
                let ex = exp(x);
                let em = expm1(x);
                let negative = 0.0 > ex;
                let grad = if negative {
                    0.0
                } else {
                    ((-1.0 * em) * ex) * -1.0
                };
                let logp = (if negative {
                    f64::NEG_INFINITY
                } else {
                    (-0.5 * (em * em)) - 0.7461847591032866
                }) + x;
                (ex, grad, logp)
            })
            .collect();
        let sum_psi_s = seq_sum(psi_s.iter().map(|p| p.2));
        let sum_psi_j = seq_sum(psi_j.iter().map(|p| p.2));
        let sum_bias = seq_sum(bias.iter().map(|p| p.2));

        // n75: over-dispersion.
        let mut psi_sum = vec![0.0; s * j];
        let mut psi_expm1 = vec![0.0; s * j];
        let mut psi_inv = vec![0.0; s * j];
        let mut alpha = vec![0.0; s * j];
        for si in 0..s {
            for jj in 0..j {
                let i = si * j + jj;
                psi_sum[i] = psi_j[jj].0 + psi_s[si].0;
                psi_expm1[i] = expm1(psi_sum[i]);
                psi_inv[i] = 1.0 / psi_expm1[i];
                alpha[i] = maximum(1e-10, psi_inv[i]);
            }
        }
        // n77, n78, n81: the NB means' pieces.
        let mu_num: Vec<f64> = (0..j * k)
            .map(|i| (self.t_j[i / k] * bias[i / k].0) * (i % k) as f64)
            .collect();
        let gamma: Vec<f64> = (0..s * j)
            .map(|i| (mean_ploidy[i] * self.t_j[i % j]) * bias[i % j].0)
            .collect();
        let rest = self.dot_mask(&gamma, false);
        // n82, n83.
        let mut binom = vec![0.0; s * j];
        let mut n_alpha = vec![0.0; s * j];
        for i in 0..s * j {
            let n = self.n_sj[i];
            n_alpha[i] = n + alpha[i];
            binom[i] = (lgamma(n_alpha[i]) - self.lgamma_1n[i]) - lgamma(alpha[i]);
        }

        // n88: the NB logp and its gradients, then their sums over k (n89 to n96).
        let sjk = s * j * k;
        let (mut o0, mut o2, mut o3, mut o4, mut o6, mut o7, mut o8) = (
            vec![0.0; sjk],
            vec![0.0; sjk],
            vec![0.0; sjk],
            vec![0.0; sjk],
            vec![0.0; sjk],
            vec![0.0; sjk],
            vec![0.0; sjk],
        );
        let mut all_mean_positive = true;
        for si in 0..s {
            for jj in 0..j {
                let sj = si * j + jj;
                let al = alpha[sj];
                let n = self.n_sj[sj];
                let switch_alpha = if al == 0.0 { 0.0 } else { f64::NEG_INFINITY };
                for kk in 0..k {
                    let x = sj * k + kk;
                    let num = mu_num[jj * k + kk];
                    let tmp1 = rest[sj] + num;
                    let tmp2 = num / tmp1;
                    let tmp3 = self.one_minus_error * tmp2;
                    let tmp4 = tmp1 == 0.0;
                    let tmp5 = if tmp4 { 0.0 } else { tmp3 };
                    let tmp6 = tmp5 + self.eps_j[jj];
                    let tmp7 = tmp6 * self.n_s[si];
                    let tmp8 = tmp7 + al;
                    let qx = q[x];
                    let tmp9 = (al * -1.0) * qx;
                    let tmp10 = al / tmp8;
                    let tmp11 = tmp10 == 0.0;
                    let tmp12 = if tmp11 { 0.0 } else { tmp9 };
                    o0[x] = -tmp12 / tmp8;
                    let tmp14 = tmp1 * tmp1;
                    let tmp15 = (n * -1.0) * qx;
                    let tmp16 = tmp7 / tmp8;
                    let tmp17 = tmp16 == 0.0;
                    let tmp18 = if tmp17 { 0.0 } else { tmp15 };
                    o4[x] = -tmp18 / tmp8;
                    let tmp20 = tmp15 / tmp7;
                    let tmp21 = if tmp17 { 0.0 } else { tmp20 };
                    let tmp22 = (tmp21 + o4[x]) + o0[x];
                    let tmp23 = ((((-self.one_minus_error * self.t_j[jj]) * bias[jj].0)
                        * kk as f64)
                        * tmp22)
                        * self.n_s[si];
                    o2[x] = if tmp4 { 0.0 } else { tmp23 / tmp14 };
                    let tmp25 = (self.one_minus_error * tmp22) * self.n_s[si];
                    o3[x] = if tmp4 { 0.0 } else { tmp25 / tmp1 };
                    all_mean_positive &= tmp7 > 0.0;
                    let tmp27 = log(tmp10);
                    o6[x] = if tmp11 { 0.0 } else { (tmp27 * -1.0) * qx };
                    o7[x] = if tmp11 { 0.0 } else { -1.0 * qx };
                    let tmp31 = if tmp11 { switch_alpha } else { al * tmp27 };
                    let tmp33 = n * log(tmp16);
                    let tmp34 = if tmp17 { self.switch_n[sj] } else { tmp33 };
                    o8[x] = (binom[sj] + tmp34) + tmp31;
                }
            }
        }
        let over_k = |v: &[f64]| -> Vec<f64> {
            (0..s * j)
                .map(|sj| seq_sum((0..k).map(|kk| v[sj * k + kk])))
                .collect()
        };
        let over_s = |v: &[f64]| -> Vec<f64> {
            (0..j * k)
                .map(|jk| seq_sum((0..s).map(|si| v[si * j * k + jk])))
                .collect()
        };
        let n89 = over_k(&q);
        let n90 = over_k(&o6);
        let n91 = over_k(&o0);
        let n92 = over_k(&o7);
        let n93 = over_k(&o4);
        let n94 = over_s(&o2);
        let n95 = over_k(&o2);
        let n96 = over_s(&o3);
        let valid = all_mean_positive
            && self.n_sj.iter().all(|&n| n >= 0.0)
            && alpha.iter().all(|&a| a > 0.0);

        // n98: the gradient through alpha.
        let n98: Vec<f64> = (0..s * j)
            .map(|i| {
                let tmp1 = exp(psi_sum[i]);
                let tmp2 = psi_expm1[i] * psi_expm1[i];
                let tmp4 = -1.0 * n89[i];
                let tmp5 = (-1.0 * tmp4) * psi(alpha[i]);
                let tmp7 = tmp4 * psi(n_alpha[i]);
                let tmp8 = ((((tmp7 + tmp5) + n93[i]) + n92[i]) + n91[i]) + n90[i];
                let tmp9 = if alpha[i] == psi_inv[i] { 1.0 } else { 0.0 };
                let tmp11 = -(tmp9 * tmp8);
                (tmp11 / tmp2) * tmp1
            })
            .collect();
        let n100 = self.dot_mask(&n95, true);
        let n103: Vec<f64> = (0..s)
            .map(|si| seq_sum((0..j).map(|jj| n98[si * j + jj])))
            .collect();
        let n104: Vec<f64> = (0..j)
            .map(|jj| seq_sum((0..s).map(|si| n98[si * j + jj])))
            .collect();
        let n105: Vec<f64> = (0..s * j).map(|i| n100[i] * mean_ploidy[i]).collect();
        let n106: Vec<f64> = (0..j * k)
            .map(|i| (n94[i] + n96[i]) * (i % k) as f64)
            .collect();
        let n110: Vec<f64> = (0..j)
            .map(|jj| seq_sum((0..s).map(|si| n105[si * j + jj])))
            .collect();
        let n111: Vec<f64> = (0..j)
            .map(|jj| seq_sum((0..k).map(|kk| n106[jj * k + kk])))
            .collect();
        let nb_logp =
            seq_sum((0..sjk).map(|x| q[x] * if valid { o8[x] } else { f64::NEG_INFINITY }));

        // n108, n109, n115: d(-logp)/dz.
        let mut grad = vec![0.0; d];
        for jj in 0..j {
            let t = self.t_j[jj];
            let v = (bias[jj].1 + -1.0) + (((t * n110[jj]) + (n111[jj] * t)) * bias[jj].0);
            grad[jj] = 0.0 + (0.0 + v);
            let v = (psi_j[jj].1 + -1.0) + (n104[jj] * psi_j[jj].0);
            grad[j + jj] = 0.0 + (0.0 + v);
        }
        for si in 0..s {
            let v = (psi_s[si].1 + -1.0) + (n103[si] * psi_s[si].0);
            grad[2 * j + si] = 0.0 + (0.0 + v);
        }
        let logp = ((sum_bias + sum_psi_j) + sum_psi_s) + nb_logp;
        let loss = (temperature * log_q_z) - logp;

        // n127, n128: Adamax.
        for i in 0..d {
            let g_rho = grad_entropy[i] + (grad[i] * (0.0 + e[i])) * sigmoid[i];
            st.u_rho[i] = maximum(g_rho.abs(), b2 * st.u_rho[i]);
            st.m_rho[i] = (b1 * st.m_rho[i]) + ((1.0 - b1) * g_rho);
            st.rho[i] -= (a * st.m_rho[i]) / (1e-08 + st.u_rho[i]);
            st.u_mu[i] = maximum(grad[i].abs(), b2 * st.u_mu[i]);
            st.m_mu[i] = (b1 * st.m_mu[i]) + ((1.0 - b1) * grad[i]);
            st.mu[i] -= (a * st.m_mu[i]) / (1e-08 + st.u_mu[i]);
        }
        st.beta1_power = beta1_power;
        if let Some(drop) = drop {
            st.temperature = maximum(temperature - drop, 1.0);
        }
        loss
    }

    /// The log-emission sampler (graph 3): the mean of `logp_sjk` over `samples` draws.
    fn sample(&self, st: &State, rng: &mut Pcg64, samples: usize) -> Result<Vec<f64>, String> {
        let (s, j, k, d) = (self.s, self.j, self.k, self.d());
        let sigma: Vec<f64> = st.rho.iter().map(|&r| softplus(r)).collect();
        let mean_ploidy_t: Vec<f64> = self
            .mean_ploidy(&st.log_q)
            .iter()
            .enumerate()
            .map(|(i, &m)| m * self.t_j[i % j])
            .collect();
        let mut acc = vec![0.0; s * j * k];
        let mut e = vec![0.0; d];
        for _ in 0..samples {
            for x in e.iter_mut() {
                *x = 0.0 + 1.0 * rng.standard_normal();
            }
            let z: Vec<f64> = (0..d).map(|i| sigma[i] * e[i] + st.mu[i]).collect();
            let e_s: Vec<f64> = (0..s).map(|i| exp(z[2 * j + i])).collect();
            let e_j: Vec<f64> = (0..j).map(|i| exp(z[j + i])).collect();
            let bias: Vec<f64> = (0..j).map(|i| exp(z[i])).collect();
            let mut valid = true;
            let mut alpha = vec![0.0; s * j];
            let mut binom = vec![0.0; s * j];
            for si in 0..s {
                for jj in 0..j {
                    let i = si * j + jj;
                    let a = maximum(1e-10, 1.0 / expm1(e_j[jj] + e_s[si]));
                    valid &= a > 0.0;
                    alpha[i] = a;
                    binom[i] = (lgamma(self.n_sj[i] + a) - self.lgamma_1n[i]) - lgamma(a);
                }
            }
            let gamma: Vec<f64> = (0..s * j).map(|i| mean_ploidy_t[i] * bias[i % j]).collect();
            let rest = self.dot_mask(&gamma, false);
            let mut logp = vec![0.0; s * j * k];
            for si in 0..s {
                for jj in 0..j {
                    let sj = si * j + jj;
                    let al = alpha[sj];
                    let switch_alpha = if al == 0.0 { 0.0 } else { f64::NEG_INFINITY };
                    for kk in 0..k {
                        let num = (self.t_j[jj] * bias[jj]) * kk as f64;
                        let tmp1 = rest[sj] + num;
                        let tmp3 = self.one_minus_error * (num / tmp1);
                        let tmp5 = if tmp1 == 0.0 { 0.0 } else { tmp3 };
                        let tmp7 = (tmp5 + self.eps_j[jj]) * self.n_s[si];
                        valid &= tmp7 > 0.0;
                        let tmp8 = tmp7 + al;
                        let tmp9 = al / tmp8;
                        let tmp13 = if tmp9 == 0.0 {
                            switch_alpha
                        } else {
                            al * log(tmp9)
                        };
                        let tmp14 = tmp7 / tmp8;
                        let tmp18 = if tmp14 == 0.0 {
                            self.switch_n[sj]
                        } else {
                            self.n_sj[sj] * log(tmp14)
                        };
                        logp[sj * k + kk] = (binom[sj] + tmp18) + tmp13;
                    }
                }
            }
            if !valid || self.n_sj.iter().any(|&n| n.is_nan() || n < 0.0) {
                return Err("the log-emission sampler's parameter check failed".into());
            }
            for (a, x) in acc.iter_mut().zip(&logp) {
                *a += x;
            }
        }
        let scale = 1.0 / samples as f64;
        Ok(acc.iter().map(|&x| scale * x).collect())
    }

    /// The caller (graph 4): admixes the new ploidy posterior into `log_q` and returns the
    /// Hellinger distance per sample and contig.
    fn call(&self, st: &mut State, admixing: f64) -> Vec<f64> {
        let (s, j, k) = (self.s, self.j, self.k);
        let log_r = log(admixing);
        let log_1mr = log(1.0 - admixing);
        let diff_constant = log_1mr - log_r;
        let mut hellinger = vec![0.0; s * j];
        for sj in 0..s * j {
            let jj = sj % j;
            let x: Vec<f64> = (0..k)
                .map(|kk| self.log_p_jk[jj * k + kk] + st.log_emission[sj * k + kk])
                .collect();
            let m1 = x
                .iter()
                .fold(f64::NEG_INFINITY, |m, &v| if v > m { v } else { m });
            let shifted: Vec<f64> = x.iter().map(|&v| v - m1).collect();
            let m2 = shifted
                .iter()
                .fold(f64::NEG_INFINITY, |m, &v| if v > m { v } else { m });
            let e_m2 = exp(m2);
            let total =
                seq_sum(
                    shifted
                        .iter()
                        .map(|&v| if m2.is_infinite() { e_m2 } else { exp(v - m2) }),
                );
            let lse_low = m2 + log(total);
            let lse = lse_low + m1;
            let mut squares = vec![0.0; k];
            for kk in 0..k {
                let i = sj * k + kk;
                let old = st.log_q[i];
                let tmp2 = (((diff_constant + old) + lse_low) + m1) - x[kk];
                let nan = tmp2.is_nan();
                let tmp4 = if nan { 0.0 } else { tmp2 };
                let tmp7 = (((log_r + self.log_p_jk[jj * k + kk]) + st.log_emission[i])
                    + softplus(tmp4))
                    - lse;
                let tmp9 = if nan { 0.0 } else { -tmp2 };
                let tmp11 = (log_1mr + old) + softplus(tmp9);
                let new = if tmp4 >= 0.0 { tmp11 } else { tmp7 };
                let delta = exp(new).sqrt() - exp(old).sqrt();
                squares[kk] = delta * delta;
                st.log_q[i] = new;
            }
            hellinger[sj] = 0.7071067932881648 * seq_sum(squares).sqrt();
        }
        hellinger
    }
}

// ---------------------------------------------------------------------------------------------
// The script.

/// `io_ploidy.get_contig_ploidy_prior_map_from_tsv_file`, then `PloidyModelConfig`'s validation.
/// Each contig's prior, and the number of ploidy states.
type Priors = (Vec<(String, Vec<f64>)>, usize);

fn priors(text: &str) -> Result<Priors, String> {
    let mut lines = text
        .lines()
        .map(|line| line.split('@').next().unwrap_or(""))
        .filter(|line| !line.is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or("an empty priors table")?
        .split('\t')
        .collect();
    if header.len() < 2 || header[0] != "CONTIG_NAME" {
        return Err("the priors table's first column must be CONTIG_NAME".into());
    }
    let mut ploidies = Vec::new();
    for column in &header[1..] {
        let k = column
            .strip_prefix("PLOIDY_PRIOR_")
            .and_then(|k| k.parse::<usize>().ok())
            .ok_or_else(|| format!("bad priors column {column}"))?;
        ploidies.push(k);
    }
    let states = ploidies.iter().max().copied().unwrap_or(0) + 1;
    let mut map = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let mut prior = vec![0.0; states];
        for (column, &k) in ploidies.iter().enumerate() {
            let text = fields.get(column + 1).ok_or("a short priors row")?;
            prior[k] = precise_xstrtod(text).ok_or_else(|| format!("bad prior {text}"))?;
        }
        let padded: Vec<f64> = prior.iter().map(|&p| 0.0 + p).collect();
        if padded.iter().any(|&p| p.is_nan() || p < 0.0) {
            return Err("Probabilities must be non-negative".into());
        }
        let total = np_sum(&padded);
        let normalized = if (total - 1.0).abs() < 1e-10 {
            padded
        } else {
            padded.iter().map(|&p| p / total).collect()
        };
        map.push((fields[0].to_string(), normalized));
    }
    Ok((map, states))
}

/// The configurations the port follows; the rest it refuses rather than guess.
fn refusals(settings: &Settings) -> Result<(), String> {
    if settings.disable_annealing || settings.disable_caller || settings.disable_sampler {
        return Err("Disabling the annealing, the caller or the sampler".into());
    }
    // The truncated normal's normalizing constant is folded for this deviation only.
    if settings.mean_bias_sd != 1.0 {
        return Err("A --mean-bias-standard-deviation other than 1".into());
    }
    // shortcut: the ELBO convergence tracker (a recursive least squares through MKL's matrix
    // inverse) is not ported, so runs where it could stop ADVI are refused; port it to lift this.
    let total = settings.max_advi_iter_first_epoch
        + (settings.max_training_epochs - 1).max(0) * settings.max_advi_iter_subsequent_epochs;
    if total + 1 > settings.convergence_snr_averaging_window {
        return Err(format!(
            "More ADVI iterations in all ({total}) than --convergence-snr-averaging-window ({})",
            settings.convergence_snr_averaging_window
        ));
    }
    if settings.log_emission_sampling_rounds < 1 || settings.log_emission_samples_per_round < 1 {
        return Err("No log-emission sampling round".into());
    }
    Ok(())
}

/// What the replay test watches: each step's loss, `mu` and `rho`, each draw, each call.
pub enum Event<'a> {
    Step {
        loss: f64,
        mu: &'a [f64],
        rho: &'a [f64],
    },
    Draw(&'a [f64]),
    Call(&'a [f64]),
}

/// The script, from the arguments to the files it writes.
pub fn run(cohort: &Cohort, settings: &Settings) -> Result<Files, String> {
    run_watched(cohort, settings, &mut |_| {})
}

pub fn run_watched(
    cohort: &Cohort,
    settings: &Settings,
    watch: &mut dyn FnMut(Event),
) -> Result<Files, String> {
    refusals(settings)?;
    let (prior_map, k) = priors(&cohort.priors_text)?;
    if k < 2 {
        return Err("At least two ploidy states are required for genotyping".into());
    }
    let s = cohort.samples.len();
    let j = cohort.contigs.len();
    let mut log_p_jk = Vec::with_capacity(j * k);
    for contig in &cohort.contigs {
        let (_, prior) = prior_map
            .iter()
            .find(|(c, _)| c == contig)
            .ok_or("Some contigs do not have ploidy priors")?;
        log_p_jk.extend(
            prior
                .iter()
                .map(|&p| log(if p > 1e-10 { p } else { 1e-10 })),
        );
    }
    let t_j: Vec<f64> = cohort
        .intervals_per_contig
        .iter()
        .map(|&t| t as f64)
        .collect();
    let n_sj: Vec<f64> = cohort
        .samples
        .iter()
        .flat_map(|(_, n)| n.iter().map(|&c| c as f64))
        .collect();
    let n_total: Vec<i64> = cohort.samples.iter().map(|(_, n)| n.iter().sum()).collect();
    let t_sum = seq_sum(t_j.iter().copied());
    let error = settings.mapping_error_rate;
    // pymc keeps each Exponential's scale in float32.
    let scale_s = settings.psi_s_scale as f32 as f64;
    let scale_j = settings.psi_j_scale as f32 as f64;
    let model = Model {
        s,
        j,
        k,
        n_s: n_total.iter().map(|&n| n as f64).collect(),
        mask: (0..j * j)
            .map(|i| if i / j == i % j { 0.0 } else { 1.0 })
            .collect(),
        eps_j: t_j.iter().map(|&t| (error * t) / t_sum).collect(),
        one_minus_error: 1.0 - error,
        lam_s: 1.0 / scale_s,
        log_lam_s: -log(scale_s),
        lam_j: 1.0 / scale_j,
        log_lam_j: -log(scale_j),
        lgamma_1n: n_sj.iter().map(|&n| lgamma(1.0 + n)).collect(),
        switch_n: n_sj
            .iter()
            .map(|&n| if n == 0.0 { 0.0 } else { f64::NEG_INFINITY })
            .collect(),
        t_j,
        n_sj,
        log_p_jk,
    };
    let d = model.d();
    let mut mu = vec![0.0; d];
    for i in 0..j {
        mu[j + i] = log(scale_j);
    }
    for i in 0..s {
        mu[2 * j + i] = log(scale_s);
    }
    let mut st = State {
        mu,
        rho: vec![0.0; d],
        m_mu: vec![0.0; d],
        u_mu: vec![0.0; d],
        m_rho: vec![0.0; d],
        u_rho: vec![0.0; d],
        beta1_power: settings.adamax_beta1,
        temperature: settings.initial_temperature,
        log_q: (0..s)
            .flat_map(|_| model.log_p_jk.iter().copied())
            .collect(),
        log_emission: vec![0.0; s * j * k],
    };

    // `approx.rng`: the initial point's seed, the step function's, the sampler's.
    let mut approx_rng = RandomState::new(1984);
    approx_rng.randint(0, 1 << 30);
    let mut step_rng = reseed(approx_rng.randint(0, 1 << 30) as u32, 1).remove(0);
    let mut sampler_rng = reseed(approx_rng.randint(0, 1 << 30) as u32, 4).remove(0);

    let drop = ((settings.initial_temperature - 1.0).abs() >= 1e-6)
        .then(|| (settings.initial_temperature - 1.0) / settings.num_thermal_advi_iters as f64);
    let normalization = (s * j) as f64;
    for epoch in 1..=settings.max_training_epochs {
        let iterations = if epoch > 1 {
            settings.max_advi_iter_subsequent_epochs
        } else {
            settings.max_advi_iter_first_epoch
        };
        for _ in 0..iterations {
            let loss = model.step(&mut st, &mut step_rng, settings, drop);
            watch(Event::Step {
                loss,
                mu: &st.mu,
                rho: &st.rho,
            });
            if (loss / normalization).is_nan() {
                return Err("ADVI's loss is NaN (gcnvkernel's ConvergenceError)".into());
            }
        }
        if epoch == 1 {
            st.log_emission = vec![0.0; s * j * k];
        }
        for round in 0..settings.log_emission_sampling_rounds {
            let draw = model.sample(
                &st,
                &mut sampler_rng,
                settings.log_emission_samples_per_round as usize,
            )?;
            watch(Event::Draw(&draw));
            let divisor = (round + 1) as f64;
            let mut relative = Vec::with_capacity(draw.len());
            for (latest, new) in st.log_emission.iter_mut().zip(&draw) {
                let update = (new - *latest) / divisor;
                *latest += update;
                let r = (update / *latest).abs();
                relative.push(if r.is_nan() {
                    0.0
                } else if r.is_infinite() {
                    f64::MAX
                } else {
                    r
                });
            }
            if np_median(&mut relative) < settings.log_emission_sampling_median_rel_error {
                break;
            }
        }
        for _ in 0..settings.max_calling_iters {
            let hellinger = model.call(&mut st, settings.caller_external_admixing_rate);
            watch(Event::Call(&hellinger));
            let mean = np_sum(&hellinger) / (hellinger.len() as f64);
            if mean < settings.caller_update_convergence_threshold {
                break;
            }
        }
    }

    // The writers.
    let std: Vec<f64> = st.rho.iter().map(|&r| softplus(r)).collect();
    let mut model_files = vec![
        (
            "gcnvkernel_version.json".to_string(),
            "{\n \"version\": \"0.9\"\n}\n".to_string(),
        ),
        (
            "ploidy_config.json".to_string(),
            json_floats(&[
                ("mapping_error_rate", settings.mapping_error_rate),
                ("mean_bias_sd", settings.mean_bias_sd),
                ("psi_j_scale", settings.psi_j_scale),
                ("psi_s_scale", settings.psi_s_scale),
            ]),
        ),
    ];
    for (name, range) in [("mean_bias_j_interval__", 0..j), ("psi_j_log__", j..2 * j)] {
        model_files.push((
            format!("mu_{name}.tsv"),
            ndarray_tsv(&st.mu[range.clone()], &[]),
        ));
        model_files.push((format!("std_{name}.tsv"), ndarray_tsv(&std[range], &[])));
    }
    let mut calls = Vec::new();
    let total_intervals: u64 = cohort.intervals_per_contig.iter().sum();
    for (si, (name, _)) in cohort.samples.iter().enumerate() {
        let dir = format!("SAMPLE_{si}");
        let comment = format!("RG\tID:GATKCopyNumber\tSM:{name}");
        let mut ploidy_table = format!("@{comment}\nCONTIG\tPLOIDY\tPLOIDY_GQ\n");
        let mut copies: u64 = 0;
        for (jj, contig) in cohort.contigs.iter().enumerate() {
            let log_q = &st.log_q[(si * j + jj) * k..(si * j + jj + 1) * k];
            let mut order: Vec<usize> = (0..k).collect();
            order.sort_by(|&a, &b| (-log_q[a]).partial_cmp(&-log_q[b]).expect("no NaN"));
            let gq = 4.342944819032518 * (log_q[order[0]] - log_q[order[1]]);
            copies += cohort.intervals_per_contig[jj] * order[0] as u64;
            ploidy_table.push_str(&format!("{contig}\t{}\t{}\n", order[0], python_repr(gq)));
        }
        let depth = n_total[si] as f64 / copies as f64;
        let average = copies as f64 / total_intervals as f64;
        calls.push((format!("{dir}/contig_ploidy.tsv"), ploidy_table));
        calls.push((
            format!("{dir}/global_read_depth.tsv"),
            format!(
                "@{comment}\nGLOBAL_READ_DEPTH\tAVERAGE_PLOIDY\n{}\t{}\n",
                python_repr(depth),
                python_repr(average)
            ),
        ));
        calls.push((format!("{dir}/sample_name.txt"), format!("{name}\n")));
        let i = 2 * j + si;
        calls.push((
            format!("{dir}/mu_psi_s_log__.tsv"),
            ndarray_tsv(&st.mu[i..i + 1], std::slice::from_ref(&comment)),
        ));
        calls.push((
            format!("{dir}/std_psi_s_log__.tsv"),
            ndarray_tsv(&std[i..i + 1], std::slice::from_ref(&comment)),
        ));
    }
    Ok(Files {
        model: model_files,
        calls,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_folded_constants_are_the_references() {
        let s = 1e-4f32 as f64;
        let j = 1e-3f32 as f64;
        assert_eq!(1.0 / s, 10000.000252621256);
        assert_eq!(-log(s), 9.210340397238308);
        assert_eq!(1.0 / j, 999.999952502551);
        assert_eq!(-log(j), 6.907755231484687);
        assert_eq!(log(j).to_bits(), 0xc01ba18a965fffa2);
        assert_eq!(log(s).to_bits(), 0xc0226bb1bc8e5516);
        assert_eq!(log(0.25) - log(0.75), -1.0986122886681096);
        assert_eq!(log(0.75), -0.2876820724517809);
        assert_eq!(log(0.25), -1.3862943611198906);
    }

    #[test]
    fn lgamma_and_psi_agree_with_known_values() {
        assert_eq!(lgamma(1.0), 0.0);
        assert_eq!(lgamma(2.0), 0.0);
        // glibc 2.35's `lgamma` in the image, under the frozen reference's tunables.
        for (x, want) in [
            (0.5, 0x3fe250d048e7a1bd),
            (0.9, 0x3fb0fe0882e465b2),
            (1.3, 0xbfbbb158251658ad),
            (1.8, 0xbfb2328d7d70b422),
            (2.5, 0x3fd2383e809a67e8),
            (5.5, 0x400fa99a5e94985a),
            (7.9, 0x4020a6062d049d21),
            (6365.0, 0x40e81c7ab5142463),
            (1e5, 0x41300a97b57f4c2e),
            (1e-10, 0x4037069e2aa26ae4),
            (1e-25, 0x404cc845b54b54f2),
            (0.2, 0x3ff86290bf25b626),
            (1e20, 0x446e873977d7ba97u64),
        ] {
            assert_eq!(lgamma(x).to_bits(), want, "{x}");
        }
        assert!((psi(1.0) + 0.5772156649015329).abs() < 1e-9);
        assert!((psi(1000.0) - 6.907255195648812).abs() < 1e-12);
    }
}
