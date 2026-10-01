//! The Markov-chain machinery `ModelSegments` fits its two models with.
//!
//! Ported from `org.broadinstitute.hellbender.utils.mcmc` (`AbstractSliceSampler`,
//! `MinibatchSliceSampler`, `DecileCollection`, the generator `GibbsSampler` holds) and from the
//! commons-math3 3.5 pieces those reach: `ExponentialDistribution.sample`,
//! `BetaDistribution.sample` and `logDensity`, `TDistribution.cumulativeProbability`,
//! `NormalDistribution.density`, `Primes.nextPrime`, `Mean`, `Variance`, `FastMath.pow(double, int)`
//! and the `BrentSolver` the inverse cumulative probability runs. GATK pins commons-math3 to 3.5
//! strictly, and 3.5 is what this follows: in 3.6 `BetaDistribution.sample` became Cheng's
//! algorithm, which draws a different number of values from the generator.
//!
//! # One generator, restarted at every chain
//!
//! `GibbsSampler` holds a STATIC `RandomGeneratorFactory` adapter over `new Random(42)` and calls
//! `setSeed(42)` at the top of every `runMCMC`, so each fit starts the same stream whatever ran
//! before it. The adapter forwards `nextDouble`, `nextBoolean` and `nextInt(n)` to the wrapped
//! `java.util.Random`, which is [`crate::java_random::JavaRandom`]. The distributions a sampler
//! builds are handed the same generator, so every draw of every sampler interleaves in one stream
//! and the order of the draws is the whole comparison.
//!
//! # `Math.pow(x, 2)` is the square
//!
//! `MinibatchSliceSampler` squares its running mean with `Math.pow(mean, 2)`. The correctly rounded
//! square of a double is the one IEEE multiplication `x * x`, and the platform's `pow` returns that
//! at an exponent of two, so the port multiplies.
//!
//! # Caches are values here
//!
//! The reference caches the prior and the likelihoods at the current sample, in a `HashMap` keyed
//! by the data point. The functions cached are pure and the points distinct, so a cached value is
//! the value recomputed; the port keeps the cache by index for speed and nothing else.

use crate::java_random::JavaRandom;

/// `IllegalArgumentException` from a `Utils.validateArg` or a `ParamUtils` check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IllegalArgument(pub String);

/// `GibbsSampler.RANDOM_SEED`.
pub const GIBBS_RANDOM_SEED: i64 = 42;

/// `FastMath.pow(double, int)`, the Veltkamp split product of commons-math3 3.5.
pub fn fast_math_pow_int(d: f64, e: i32) -> f64 {
    if e == 0 {
        return 1.0;
    }
    let (mut d, mut e) = (d, e);
    if e < 0 {
        e = -e;
        d = 1.0 / d;
    }
    let split_factor = 134_217_729.0; // 0x8000001
    let cd = split_factor * d;
    let d1_high = cd - (cd - d);
    let d1_low = d - d1_high;

    let mut result_high = 1.0;
    let mut result_low = 0.0;

    let mut d2p = d;
    let mut d2p_high = d1_high;
    let mut d2p_low = d1_low;

    while e != 0 {
        if (e & 0x1) != 0 {
            let tmp_high = result_high * d2p;
            let c_rh = split_factor * result_high;
            let r_hh = c_rh - (c_rh - result_high);
            let r_hl = result_high - r_hh;
            let tmp_low =
                r_hl * d2p_low - (((tmp_high - r_hh * d2p_high) - r_hl * d2p_high) - r_hh * d2p_low);
            result_high = tmp_high;
            result_low = result_low * d2p + tmp_low;
        }
        let tmp_high = d2p_high * d2p;
        let c_d2p_h = split_factor * d2p_high;
        let d2p_hh = c_d2p_h - (c_d2p_h - d2p_high);
        let d2p_hl = d2p_high - d2p_hh;
        let tmp_low =
            d2p_hl * d2p_low - (((tmp_high - d2p_hh * d2p_high) - d2p_hl * d2p_high) - d2p_hh * d2p_low);
        let c_tmp_h = split_factor * tmp_high;
        d2p_high = c_tmp_h - (c_tmp_h - tmp_high);
        d2p_low = d2p_low * d2p + tmp_low + (tmp_high - d2p_high);
        d2p = d2p_high + d2p_low;
        e >>= 1;
    }
    result_high + result_low
}

/// `ExponentialDistribution.EXPONENTIAL_SA_QI`, the table the static initializer fills:
/// `qi += FastMath.pow(LN2, i) / CombinatoricsUtils.factorial(i)` until `qi` reaches one.
fn exponential_sa_qi() -> &'static [f64] {
    use std::sync::OnceLock;
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let ln2 = jmath::fast_math::log(2.0);
        let mut table = Vec::new();
        let mut qi = 0.0;
        let mut i = 1;
        let mut factorial: i64 = 1;
        while qi < 1.0 {
            factorial *= i64::from(i);
            qi += fast_math_pow_int(ln2, i) / factorial as f64;
            table.push(qi);
            i += 1;
        }
        table
    })
}

/// `ExponentialDistribution.sample()` in 3.5: Ahrens and Dieter's algorithm SA, which draws one
/// value and then as many more as the table walk needs.
pub fn exponential_sample(rng: &mut JavaRandom, mean: f64) -> f64 {
    let qi = exponential_sa_qi();
    let mut a = 0.0;
    let mut u = rng.next_double();
    while u < 0.5 {
        a += qi[0];
        u *= 2.0;
    }
    u += u - 1.0;
    if u <= qi[0] {
        return mean * (a + u);
    }
    let mut i = 0;
    let mut u2 = rng.next_double();
    let mut umin = u2;
    loop {
        i += 1;
        u2 = rng.next_double();
        if u2 < umin {
            umin = u2;
        }
        if u <= qi[i] {
            break;
        }
    }
    mean * (a + umin * qi[0])
}

/// `Beta.regularizedBeta(x, a, b)`, whose only error is an iteration budget of
/// `Integer.MAX_VALUE`, which no argument here reaches.
fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    jmath::beta::regularized_beta(x, a, b).unwrap_or(f64::NAN)
}

/// `BetaDistribution.cumulativeProbability`.
fn beta_cumulative(x: f64, alpha: f64, beta: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if x >= 1.0 {
        1.0
    } else {
        regularized_beta(x, alpha, beta)
    }
}

/// `Precision.equals(x, 0)`: within one ulp of either zero, and not NaN.
fn within_one_ulp_of_zero(x: f64) -> bool {
    if x.is_nan() {
        return false;
    }
    let bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    bits <= 1
}

/// `BrentSolver.solve(Integer.MAX_VALUE, f, min, max)` with the given absolute accuracy and the
/// defaults for the other two (relative `1e-14`, function value `1e-15`), started at the midpoint.
///
/// The one refusal, `NoBracketingException`, is `None`: a distribution function minus a
/// probability in `(0, 1)` always brackets on `[0, 1]`, so the inverse this serves never reaches it.
pub fn brent_solve(f: impl Fn(f64) -> f64, min: f64, max: f64, absolute: f64) -> Option<f64> {
    let relative = 1e-14;
    let function_value_accuracy = 1e-15;
    let initial = min + 0.5 * (max - min);
    let y_initial = f(initial);
    if y_initial.abs() <= function_value_accuracy {
        return Some(initial);
    }
    let y_min = f(min);
    if y_min.abs() <= function_value_accuracy {
        return Some(min);
    }
    if y_initial * y_min < 0.0 {
        return Some(brent(&f, min, initial, y_min, y_initial, absolute, relative));
    }
    let y_max = f(max);
    if y_max.abs() <= function_value_accuracy {
        return Some(max);
    }
    if y_initial * y_max < 0.0 {
        return Some(brent(&f, initial, max, y_initial, y_max, absolute, relative));
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn brent(
    f: &impl Fn(f64) -> f64,
    lo: f64,
    hi: f64,
    f_lo: f64,
    f_hi: f64,
    t: f64,
    eps: f64,
) -> f64 {
    let mut a = lo;
    let mut fa = f_lo;
    let mut b = hi;
    let mut fb = f_hi;
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;
    loop {
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol = 2.0 * eps * b.abs() + t;
        let m = 0.5 * (c - b);
        if m.abs() <= tol || within_one_ulp_of_zero(fb) {
            return b;
        }
        if e.abs() < tol || fa.abs() <= fb.abs() {
            d = m;
            e = d;
        } else {
            let mut s = fb / fa;
            let mut p;
            let mut q;
            if a == c {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            s = e;
            e = d;
            if p >= 1.5 * m * q - (tol * q).abs() || p >= (0.5 * s * q).abs() {
                d = m;
                e = d;
            } else {
                d = p / q;
            }
        }
        a = b;
        fa = fb;
        if d.abs() > tol {
            b += d;
        } else if m > 0.0 {
            b += tol;
        } else {
            b -= tol;
        }
        fb = f(b);
        if (fb > 0.0 && fc > 0.0) || (fb <= 0.0 && fc <= 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
    }
}

/// `new BetaDistribution(rng, alpha, beta).sample()` in 3.5, which is
/// `inverseCumulativeProbability(random.nextDouble())` solved by Brent at `1e-9` on the support.
pub fn beta_sample(rng: &mut JavaRandom, alpha: f64, beta: f64) -> f64 {
    let p = rng.next_double();
    if p == 0.0 {
        return 0.0;
    }
    brent_solve(|x| beta_cumulative(x, alpha, beta) - p, 0.0, 1.0, 1e-9).unwrap_or(f64::NAN)
}

/// `new BetaDistribution(null, alpha, beta).logDensity(x)`.
///
/// The two refusals at the ends of the support are the `NumberIsTooSmallException`s the
/// reference raises for a shape below one there.
pub fn beta_log_density(x: f64, alpha: f64, beta: f64) -> Result<f64, IllegalArgument> {
    let z = jmath::gamma::log_gamma(alpha) + jmath::gamma::log_gamma(beta)
        - jmath::gamma::log_gamma(alpha + beta);
    if x < 0.0 || x > 1.0 {
        Ok(f64::NEG_INFINITY)
    } else if x == 0.0 {
        if alpha < 1.0 {
            return Err(IllegalArgument(format!(
                "Cannot compute beta density at 0 when alpha = {}",
                crate::tsv_table::java_double_to_string(alpha)
            )));
        }
        Ok(f64::NEG_INFINITY)
    } else if x == 1.0 {
        if beta < 1.0 {
            return Err(IllegalArgument(format!(
                "Cannot compute beta density at 1 when beta = {}",
                crate::tsv_table::java_double_to_string(beta)
            )));
        }
        Ok(f64::NEG_INFINITY)
    } else {
        let log_x = jmath::fast_math::log(x);
        let log1m_x = jmath::fast_math::log1p(-x);
        Ok((alpha - 1.0) * log_x + (beta - 1.0) * log1m_x - z)
    }
}

/// `new TDistribution(null, degreesOfFreedom).cumulativeProbability(x)`.
pub fn t_cumulative(x: f64, degrees_of_freedom: f64) -> f64 {
    if x == 0.0 {
        return 0.5;
    }
    let t = regularized_beta(
        degrees_of_freedom / (degrees_of_freedom + (x * x)),
        0.5 * degrees_of_freedom,
        0.5,
    );
    if x < 0.0 {
        0.5 * t
    } else {
        1.0 - 0.5 * t
    }
}

/// `new NormalDistribution(null, mean, sd).density(x)` in 3.5: `FastMath.exp` of the log density,
/// whose constant is `FastMath.log(sd) + 0.5 * FastMath.log(2 * PI)`.
pub fn normal_density(mean: f64, standard_deviation: f64, x: f64) -> f64 {
    let log_sd_plus_half_log_2pi = jmath::fast_math::log(standard_deviation)
        + 0.5 * jmath::fast_math::log(2.0 * std::f64::consts::PI);
    let x0 = x - mean;
    let x1 = x0 / standard_deviation;
    jmath::fast_math::exp(-0.5 * x1 * x1 - log_sd_plus_half_log_2pi)
}

/// `Primes.nextPrime(n)`: the smallest prime at least `n`, and two for zero and one.
pub fn next_prime(n: i32) -> i32 {
    let mut candidate = n.max(2);
    while !is_prime(candidate) {
        candidate += 1;
    }
    candidate
}

fn is_prime(n: i32) -> bool {
    if n < 2 {
        return false;
    }
    let n = i64::from(n);
    let mut divisor = 2i64;
    while divisor * divisor <= n {
        if n % divisor == 0 {
            return false;
        }
        divisor += 1;
    }
    true
}

/// `new Mean().evaluate(values)`: the definitional mean, then the second-pass correction.
pub fn commons_mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let sample_size = values.len() as f64;
    let mut sum = 0.0;
    for value in values {
        sum += value;
    }
    let xbar = sum / sample_size;
    let mut correction = 0.0;
    for value in values {
        correction += value - xbar;
    }
    xbar + (correction / sample_size)
}

/// `new Variance().evaluate(values)`, bias-corrected, zero for one value and NaN for none.
pub fn commons_variance(values: &[f64]) -> f64 {
    match values.len() {
        0 => f64::NAN,
        1 => 0.0,
        length => {
            let mean = commons_mean(values);
            let mut accum = 0.0;
            let mut accum2 = 0.0;
            for value in values {
                let dev = value - mean;
                accum += dev * dev;
                accum2 += dev;
            }
            let len = length as f64;
            (accum - (accum2 * accum2 / len)) / (len - 1.0)
        }
    }
}

/// `DoubleStream.sum()` as the oracle's JDK 17.0.19 computes it: Kahan summation whose final
/// step SUBTRACTS the stored compensation.
///
/// The accumulator keeps the rounding error with the opposite sign to the running sum, and the
/// final sum is `sum - compensation`. [`crate::allele_fraction_cluster::double_stream_sum`] adds
/// it instead, which is a different double whenever the compensation is not zero: two
/// allele-fraction likelihoods of one segment, `-20.337...` and `-32.733...`, sum to
/// `0xc04a890dc9de42ed` on the oracle and `...42ef` when the compensation is added. Measured by
/// tracing `AlleleFractionLikelihoods.segmentLogLikelihood` inside the oracle container.
pub fn double_stream_sum(values: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut compensation = 0.0;
    let mut simple = 0.0;
    for value in values {
        let tmp = value - compensation;
        let velvel = sum + tmp;
        compensation = (velvel - sum) - tmp;
        sum = velvel;
        simple += value;
    }
    let total = sum - compensation;
    if total.is_nan() && simple.is_infinite() {
        simple
    } else {
        total
    }
}

/// `DoubleStream.average()`: the compensated sum over the count, or `None` for no values.
pub fn double_stream_average(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        None
    } else {
        Some(double_stream_sum(values) / values.len() as f64)
    }
}

/// `DecileCollection`: the nine deciles of the samples under commons-math3's default
/// `Percentile`, first to ninth.
pub fn deciles(samples: &[f64]) -> [f64; 9] {
    let mut out = [0.0; 9];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = jmath::percentile::evaluate(
            samples,
            10.0 * (i + 1) as f64,
            jmath::percentile::EstimationType::Legacy,
        );
    }
    out
}

/// `AbstractSliceSampler.MAXIMUM_NUMBER_OF_DOUBLINGS`.
const MAXIMUM_NUMBER_OF_DOUBLINGS: i32 = 16;
/// `AbstractSliceSampler.MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS`.
const MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS: i32 = 100;
/// `AbstractSliceSampler.EPSILON`.
const SLICE_EPSILON: f64 = 1e-10;

/// `MinibatchSliceSampler`: slice sampling of a posterior given as a prior and a per-point
/// likelihood, with Dubois et al.'s minibatch test deciding when enough points have been seen.
pub struct MinibatchSliceSampler<'a, D, P, L>
where
    P: Fn(f64) -> Result<f64, IllegalArgument>,
    L: Fn(&D, f64) -> f64,
{
    data: &'a [D],
    log_prior: P,
    log_likelihood: L,
    x_min: f64,
    x_max: f64,
    width: f64,
    minibatch_size: usize,
    approx_threshold: f64,
    x_sample_cache: Option<f64>,
    log_prior_cache: f64,
    log_likelihoods_cache: Vec<Option<f64>>,
}

impl<'a, D, P, L> MinibatchSliceSampler<'a, D, P, L>
where
    P: Fn(f64) -> Result<f64, IllegalArgument>,
    L: Fn(&D, f64) -> f64,
{
    /// The constructor's checks, in the reference's order: the superclass's bounds and width,
    /// then the minibatch size and the threshold.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        data: &'a [D],
        log_prior: P,
        log_likelihood: L,
        x_min: f64,
        x_max: f64,
        width: f64,
        minibatch_size: usize,
        approx_threshold: f64,
    ) -> Result<Self, IllegalArgument> {
        if !(x_min < x_max) {
            return Err(IllegalArgument(
                "Maximum bound must be greater than minimum bound.".to_string(),
            ));
        }
        if !(width > 0.0) {
            return Err(IllegalArgument(
                "Slice-sampling width must be positive.".to_string(),
            ));
        }
        if minibatch_size <= 1 {
            return Err(IllegalArgument(
                "Minibatch size must be greater than 1.".to_string(),
            ));
        }
        if !(approx_threshold >= 0.0) {
            return Err(IllegalArgument(
                "Minibatch approximation threshold must be non-negative.".to_string(),
            ));
        }
        Ok(Self {
            data,
            log_prior,
            log_likelihood,
            x_min,
            x_max,
            width,
            minibatch_size,
            approx_threshold,
            x_sample_cache: None,
            log_prior_cache: 0.0,
            log_likelihoods_cache: Vec::new(),
        })
    }

    /// `AbstractSliceSampler.sample(xInitial)`: one draw, Neal's doubling then shrinking.
    pub fn sample(&mut self, rng: &mut JavaRandom, x_initial: f64) -> Result<f64, IllegalArgument> {
        if !(self.x_min <= x_initial && x_initial <= self.x_max) {
            return Err(IllegalArgument(
                "Initial point in slice sampler is not within specified range.".to_string(),
            ));
        }
        let x_sample = java_min(
            java_max(x_initial, self.x_min + SLICE_EPSILON),
            self.x_max - SLICE_EPSILON,
        );
        let z = exponential_sample(rng, 1.0);
        let mut x_left = x_sample - self.width * rng.next_double();
        let mut x_right = x_left + self.width;

        let mut k = MAXIMUM_NUMBER_OF_DOUBLINGS;
        while k > 0 && {
            let left = self.is_greater_than_slice_height(rng, x_left, x_sample, z)?;
            left || self.is_greater_than_slice_height(rng, x_right, x_sample, z)?
        } {
            if rng.next_boolean() {
                x_left -= x_right - x_left;
            } else {
                x_right += x_right - x_left;
            }
            k -= 1;
        }

        let mut num_iterations = 1;
        let mut x_proposed = rng.next_double() * (x_right - x_left) + x_left;
        while num_iterations <= MAXIMUM_NUMBER_OF_SLICE_SAMPLINGS {
            if self.is_greater_than_slice_height(rng, x_proposed, x_sample, z)? {
                break;
            }
            if x_proposed < x_sample {
                x_left = x_proposed;
            } else {
                x_right = x_proposed;
            }
            x_proposed = rng.next_double() * (x_right - x_left) + x_left;
            num_iterations += 1;
        }
        Ok(java_min(
            java_max(x_proposed, self.x_min + SLICE_EPSILON),
            self.x_max - SLICE_EPSILON,
        ))
    }

    /// `isGreaterThanSliceHeight`, the OnSlice procedure of Dubois et al.
    fn is_greater_than_slice_height(
        &mut self,
        rng: &mut JavaRandom,
        x_proposed: f64,
        x_sample: f64,
        z: f64,
    ) -> Result<bool, IllegalArgument> {
        if x_proposed < self.x_min || self.x_max < x_proposed {
            return Ok(false);
        }
        // `xSampleCache != xSample` unboxes, so a NaN sample never matches its cache.
        if self.x_sample_cache.is_none_or(|cached| cached != x_sample) {
            self.x_sample_cache = Some(x_sample);
            self.log_prior_cache = (self.log_prior)(x_sample)?;
            self.log_likelihoods_cache = vec![None; self.data.len()];
        }
        let num_data_points = self.data.len();
        if num_data_points == 0 {
            return Ok((self.log_prior)(x_proposed)? > self.log_prior_cache - z);
        }
        let mu0 = (self.log_prior_cache - (self.log_prior)(x_proposed)? - z) / num_data_points as f64;

        let num_minibatches = (num_data_points / self.minibatch_size).max(1);
        let mut order = if num_minibatches > 1 {
            DataOrder::shuffled(rng, num_data_points)
        } else {
            DataOrder::Sequential(0)
        };

        let mut num_data_indices_seen: usize = 0;
        let mut mean = 0.0;
        let mut squared_mean = 0.0;
        for minibatch_index in 0..num_minibatches {
            let start = minibatch_index * self.minibatch_size;
            let end = ((minibatch_index + 1) * self.minibatch_size).min(num_data_points);
            let actual = end - start;
            let batch: Vec<usize> = (0..actual).map(|_| order.next()).collect();

            let mut sum = 0.0;
            let mut squared_sum = 0.0;
            for index in batch {
                let at_sample = match self.log_likelihoods_cache[index] {
                    Some(value) => value,
                    None => {
                        let value = (self.log_likelihood)(&self.data[index], x_sample);
                        self.log_likelihoods_cache[index] = Some(value);
                        value
                    }
                };
                let at_proposed = (self.log_likelihood)(&self.data[index], x_proposed);
                let difference = at_proposed - at_sample;
                sum += difference;
                squared_sum += difference * difference;
            }

            let seen = num_data_indices_seen as f64;
            mean = (seen * mean + sum) / (seen + actual as f64);
            squared_mean = (seen * squared_mean + squared_sum) / (seen + actual as f64);
            num_data_indices_seen += actual;

            if num_minibatches == 1 {
                break;
            }
            let seen = num_data_indices_seen as f64;
            let s = (1.0 - seen / num_data_points as f64).sqrt()
                * ((squared_mean - mean * mean) / (seen - 1.0)).sqrt();
            let delta = 1.0 - t_cumulative(((mean - mu0) / s).abs(), seen - 1.0);
            if delta < self.approx_threshold {
                break;
            }
        }
        Ok(mean > mu0)
    }
}

/// The order `isGreaterThanSliceHeight` visits the data in: straight through for one minibatch,
/// else `lazyShuffleIterator`, a walk by a random stride modulo the next prime.
enum DataOrder {
    Sequential(usize),
    Shuffled {
        index: i64,
        increment: i64,
        prime: i64,
        size: i64,
    },
}

impl DataOrder {
    fn shuffled(rng: &mut JavaRandom, size: usize) -> Self {
        let prime = i64::from(next_prime(size as i32));
        let index = i64::from(rng.next_int_bound(size as i32)) + 1;
        DataOrder::Shuffled {
            index,
            increment: index,
            prime,
            size: size as i64,
        }
    }

    fn next(&mut self) -> usize {
        match self {
            DataOrder::Sequential(position) => {
                let current = *position;
                *position += 1;
                current
            }
            DataOrder::Shuffled {
                index,
                increment,
                prime,
                size,
            } => loop {
                *index = (*index + *increment) % *prime;
                if *index < *size {
                    return *index as usize;
                }
            },
        }
    }
}

/// `Math.min(double, double)`, which propagates NaN.
pub fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a <= b {
        if a == 0.0 && b == 0.0 && (a.is_sign_negative() || b.is_sign_negative()) {
            -0.0
        } else {
            a
        }
    } else {
        b
    }
}

/// `Math.max(double, double)`, which propagates NaN.
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a >= b {
        if a == 0.0 && b == 0.0 && (a.is_sign_positive() || b.is_sign_positive()) {
            0.0
        } else {
            a
        }
    } else {
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_starts_at_ln2_and_ends_at_one() {
        let table = exponential_sa_qi();
        assert_eq!(table[0], jmath::fast_math::log(2.0));
        assert!(*table.last().unwrap() >= 1.0);
    }

    #[test]
    fn next_prime_is_at_least_its_argument() {
        assert_eq!(next_prime(0), 2);
        assert_eq!(next_prime(1), 2);
        assert_eq!(next_prime(2), 2);
        assert_eq!(next_prime(24), 29);
        assert_eq!(next_prime(29), 29);
    }

    #[test]
    fn the_integer_power_squares_exactly() {
        assert_eq!(fast_math_pow_int(3.0, 2), 9.0);
        assert_eq!(fast_math_pow_int(2.0, -2), 0.25);
    }
}
