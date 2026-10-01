//! `VariantRecalibrator`'s model: the annotations normalised, the two Gaussian mixtures trained by
//! variational Bayes, and every variant scored against both.
//!
//! [`crate::variant_recalibrator`] holds the tranche arithmetic that reads the scores; this holds
//! what produces them, transcribed so that the ORDER of the floating-point operations and of the
//! draws from GATK's shared `java.util.Random` is the reference's:
//!
//! * **one random stream serves the whole run**: the replicate draws at startup, the jitter each
//!   annotation takes as it is decoded, the noise a missing annotation is normalised to, the k-means
//!   seeds, the covariance seeds, the marginalisation of a missing dimension and the score of a
//!   variant the positive model put at minus infinity. Every function below that draws takes the
//!   stream as an argument, so the order of the calls is the order of the draws;
//! * **the matrices are Jama's**: an inverse is an LU solve against the identity and a determinant
//!   the product of the LU diagonal, both with Jama's partial pivoting, because a different
//!   factorisation rounds differently;
//! * **a stream's `sum()` is compensated**: `DoubleStream.sum` is Kahan summation, where
//!   `MathUtils.sum` is a plain loop, and the mixture uses both;
//! * **`Math.pow` is the host's**, as everywhere in this port ([`gatk_engine::math_utils::pow10`]
//!   says why); `Math.log` and `Math.log10` are correctly rounded and are `jmath`'s.
//!
//! What this cannot promise is the last bit of a number the model report prints at `%.16E`: the
//! host `pow` is measured at 99.94% of `Math.pow`, and the E step takes it millions of times. The
//! VQSLOD a recal file carries is printed at four decimals, and that is what the scores are
//! measured by.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.vqsr.VariantDataManager`,
//! `VariantRecalibratorEngine`, `GaussianMixtureModel`, `MultivariateGaussian`, `VariantDatum` and
//! `TrainingSet`, and from `Jama.Matrix` and `Jama.LUDecomposition`, in GATK 4.6.2.0.

// The loops keep Java's indices and Java's `-1.0 * x` because the order of the operations is what
// is being reproduced; an iterator or a negation would be the same numbers today and a different
// transcription to check against the reference tomorrow.
#![allow(
    clippy::needless_range_loop,
    clippy::neg_multiply,
    clippy::assign_op_pattern
)]

use gatk_engine::java_random::JavaRandom;

/// `VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE`.
pub const MIN_ACCEPTABLE_LOD_SCORE: f64 = -20000.0;
/// `VariantRecalibratorEngine.MIN_PROB_CONVERGENCE`.
const MIN_PROB_CONVERGENCE: f64 = 2e-3;
/// `MultivariateGaussian.EPSILON`, the smallest determinant a Gaussian may be evaluated with.
const EPSILON: f64 = 1e-200;
/// `MultivariateGaussian.COVARIANCE_REGULARIZATION_EPSILON`.
const COVARIANCE_REGULARIZATION_EPSILON: f64 = 1e-6;

pub const USER_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.UserException";
pub const BAD_INPUT: &str = "org.broadinstitute.hellbender.exceptions.UserException$BadInput";
pub const POSITIVE_FAILURE: &str =
    "org.broadinstitute.hellbender.exceptions.UserException$VQSRPositiveModelFailure";
pub const NEGATIVE_FAILURE: &str =
    "org.broadinstitute.hellbender.exceptions.UserException$VQSRNegativeModelFailure";
pub const GATK_EXCEPTION: &str = "org.broadinstitute.hellbender.exceptions.GATKException";
pub const ILLEGAL_ARGUMENT: &str = "java.lang.IllegalArgumentException";
pub const NULL_POINTER: &str = "java.lang.NullPointerException";

/// An exception the model threw: its binary class name and its message.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelError {
    pub class: &'static str,
    pub message: String,
}

impl ModelError {
    fn new(class: &'static str, message: impl Into<String>) -> Self {
        Self {
            class,
            message: message.into(),
        }
    }

    /// Whether `Main` reports it as a user error rather than as a stack trace.
    pub fn is_user(&self) -> bool {
        self.class.contains("UserException")
    }
}

// ================================================================================================
// The numeric functions the reference names.
// ================================================================================================

/// `Math.pow`, the host's and never folded at compile time.
fn pow(x: f64, y: f64) -> f64 {
    std::hint::black_box(x).powf(y)
}

fn log10(x: f64) -> f64 {
    jmath::math::log10(x)
}

fn ln(x: f64) -> f64 {
    jmath::math::log(x)
}

/// `Gamma.digamma`, which answers NaN where the port's own refuses a non-terminating input.
fn digamma(x: f64) -> f64 {
    jmath::gamma::digamma(x).unwrap_or(f64::NAN)
}

/// `Math.max(a, b)`: NaN wins, and 0.0 is above -0.0.
pub fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && a.is_sign_negative() {
        return b;
    }
    if a >= b {
        a
    } else {
        b
    }
}

/// `Math.min(a, b)`: NaN wins, and -0.0 is below 0.0.
pub fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() {
        return a;
    }
    if a == 0.0 && b == 0.0 && b.is_sign_negative() {
        return b;
    }
    if a <= b {
        a
    } else {
        b
    }
}

/// `Double.compare`: numeric, with -0.0 below 0.0 and NaN above everything.
pub fn double_compare(a: f64, b: f64) -> std::cmp::Ordering {
    if a < b {
        return std::cmp::Ordering::Less;
    }
    if a > b {
        return std::cmp::Ordering::Greater;
    }
    let bits = |value: f64| {
        if value.is_nan() {
            f64::NAN.to_bits() as i64
        } else {
            value.to_bits() as i64
        }
    };
    bits(a).cmp(&bits(b))
}

/// `DoubleStream.sum()`: Kahan summation, with the plain sum kept beside it for the one case the
/// compensated one answers NaN and the plain one an infinity.
pub fn stream_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut compensation, mut simple) = (0.0f64, 0.0f64, 0.0f64);
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

/// `MathUtils.maxElementIndex`: the FIRST of equal maxima.
fn max_element_index(array: &[f64]) -> usize {
    let mut max = 0;
    for i in 1..array.len() {
        if array[i] > array[max] {
            max = i;
        }
    }
    max
}

/// `MathUtils.log10sumLog10(double[])`.
pub fn log10_sum_log10(values: &[f64]) -> Result<f64, ModelError> {
    if values.len() < 2 {
        return Ok(values.first().copied().unwrap_or(f64::NEG_INFINITY));
    }
    let max_index = max_element_index(values);
    let max = values[max_index];
    if max == f64::NEG_INFINITY {
        return Ok(max);
    }
    // `1.0 + IndexRange.sum(...)`, the range's own sum starting at zero.
    let mut range = 0.0;
    for (i, value) in values.iter().enumerate() {
        range += if i == max_index {
            0.0
        } else {
            pow(10.0, value - max)
        };
    }
    let sum = 1.0 + range;
    if sum.is_nan() || sum == f64::INFINITY {
        return Err(ModelError::new(
            ILLEGAL_ARGUMENT,
            "log10p values must be non-infinite and non-NAN",
        ));
    }
    Ok(max + log10(sum))
}

/// `GaussianMixtureModel.nanTolerantLog10SumLog10`.
fn nan_tolerant_log10_sum_log10(values: &[f64]) -> Result<f64, ModelError> {
    if values.iter().any(|value| value.is_nan()) {
        return Ok(f64::NAN);
    }
    log10_sum_log10(values)
}

/// `MathUtils.normalizeLog10DeleteMePlease(array, false)`: the linear probabilities.
fn normalize_to_linear(array: &[f64]) -> Vec<f64> {
    let max = array[max_element_index(array)];
    let mut normalized: Vec<f64> = array.iter().map(|x| pow(10.0, x - max)).collect();
    let mut sum = 0.0;
    for value in &normalized {
        sum += value;
    }
    for value in &mut normalized {
        *value /= sum;
    }
    normalized
}

/// `MathUtils.normalizeLog10DeleteMePlease(array, true)`, which rewrites `array` in place.
fn normalize_log10_in_place(array: &mut [f64]) {
    let max = array[max_element_index(array)];
    let mut sum = 0.0;
    for x in array.iter() {
        sum += pow(10.0, x - max);
    }
    let log10_sum = log10(sum);
    for x in array.iter_mut() {
        *x = *x - max - log10_sum;
    }
}

/// `MathUtils.normalDistributionLog10(mean, sd, x)`.
fn normal_distribution_log10(mean: f64, sd: f64, x: f64) -> Result<f64, ModelError> {
    if sd < 0.0 || sd.is_nan() {
        return Err(ModelError::new(
            ILLEGAL_ARGUMENT,
            "sd: Standard deviation of normal must be > 0",
        ));
    }
    if !mean.is_finite() || !sd.is_finite() || !x.is_finite() {
        return Err(ModelError::new(
            ILLEGAL_ARGUMENT,
            "mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)",
        ));
    }
    let root_two_pi = (2.0 * std::f64::consts::PI).sqrt();
    let a = -1.0 * log10(sd * root_two_pi);
    let b = -1.0 * ((x - mean) * (x - mean) / (2.0 * (sd * sd))) / ln(10.0);
    Ok(a + b)
}

/// `MathUtils.distanceSquared`, an `IndexRange` sum from zero.
fn distance_squared(x: &[f64], y: &[f64]) -> f64 {
    let mut result = 0.0;
    for n in 0..x.len() {
        let d = x[n] - y[n];
        result += d * d;
    }
    result
}

/// `MathUtils.compareDoubles(a, b, epsilon) == 0`.
fn close(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

/// `Collections.shuffle(list, random)` over a random-access list.
pub fn shuffle<T>(list: &mut [T], random: &mut JavaRandom) {
    let mut i = list.len();
    while i > 1 {
        let j = random.next_int_bound(i as i32) as usize;
        list.swap(i - 1, j);
        i -= 1;
    }
}

// ================================================================================================
// Jama.
// ================================================================================================

/// `Jama.Matrix`, square or not, row-major as Jama keeps it.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub a: Vec<Vec<f64>>,
}

impl Matrix {
    pub fn zeros(rows: usize, columns: usize) -> Self {
        Self {
            a: vec![vec![0.0; columns]; rows],
        }
    }

    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.a[i][i] = 1.0;
        }
        m
    }

    pub fn rows(&self) -> usize {
        self.a.len()
    }

    pub fn columns(&self) -> usize {
        self.a.first().map_or(0, Vec::len)
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.a[i][j]
    }

    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.a[i][j] = value;
    }

    /// `times(double)`: a new matrix.
    pub fn times_scalar(&self, s: f64) -> Matrix {
        Matrix {
            a: self
                .a
                .iter()
                .map(|row| row.iter().map(|x| s * x).collect())
                .collect(),
        }
    }

    /// `timesEquals(double)`.
    pub fn times_equals(&mut self, s: f64) {
        for row in &mut self.a {
            for x in row.iter_mut() {
                *x = s * *x;
            }
        }
    }

    /// `plusEquals(Matrix)`.
    pub fn plus_equals(&mut self, other: &Matrix) {
        for (row, other_row) in self.a.iter_mut().zip(&other.a) {
            for (x, y) in row.iter_mut().zip(other_row) {
                *x += y;
            }
        }
    }

    pub fn transpose(&self) -> Matrix {
        let mut t = Matrix::zeros(self.columns(), self.rows());
        for i in 0..self.rows() {
            for j in 0..self.columns() {
                t.a[j][i] = self.a[i][j];
            }
        }
        t
    }

    /// `times(Matrix)`: each column of `b` copied out, then dotted with each row.
    pub fn times(&self, b: &Matrix) -> Matrix {
        let (m, n, p) = (self.rows(), self.columns(), b.columns());
        let mut c = Matrix::zeros(m, p);
        let mut column = vec![0.0; n];
        for j in 0..p {
            for k in 0..n {
                column[k] = b.a[k][j];
            }
            for i in 0..m {
                let mut s = 0.0;
                for k in 0..n {
                    s += self.a[i][k] * column[k];
                }
                c.a[i][j] = s;
            }
        }
        c
    }

    /// `inverse()`, which is `solve(identity)` through an LU decomposition for a square matrix.
    pub fn inverse(&self) -> Result<Matrix, String> {
        Lu::new(self).solve(&Matrix::identity(self.rows()))
    }

    /// `det()`, the LU diagonal's product signed by the pivots.
    pub fn det(&self) -> f64 {
        Lu::new(self).det()
    }
}

/// `Jama.LUDecomposition`: a left-looking, dot-product Crout/Doolittle with partial pivoting.
struct Lu {
    lu: Vec<Vec<f64>>,
    pivots: Vec<usize>,
    sign: f64,
}

impl Lu {
    fn new(matrix: &Matrix) -> Lu {
        let mut lu = matrix.a.clone();
        let m = matrix.rows();
        let n = matrix.columns();
        let mut pivots: Vec<usize> = (0..m).collect();
        let mut sign = 1.0;
        let mut column = vec![0.0; m];
        for j in 0..n {
            for i in 0..m {
                column[i] = lu[i][j];
            }
            for i in 0..m {
                let kmax = i.min(j);
                let mut s = 0.0;
                for k in 0..kmax {
                    s += lu[i][k] * column[k];
                }
                column[i] -= s;
                lu[i][j] = column[i];
            }
            let mut p = j;
            for i in (j + 1)..m {
                if column[i].abs() > column[p].abs() {
                    p = i;
                }
            }
            if p != j {
                lu.swap(p, j);
                pivots.swap(p, j);
                sign = -sign;
            }
            if j < m && lu[j][j] != 0.0 {
                for i in (j + 1)..m {
                    lu[i][j] /= lu[j][j];
                }
            }
        }
        Lu { lu, pivots, sign }
    }

    fn det(&self) -> f64 {
        let mut d = self.sign;
        for j in 0..self.lu.len() {
            d *= self.lu[j][j];
        }
        d
    }

    fn solve(&self, b: &Matrix) -> Result<Matrix, String> {
        let n = self.lu.len();
        if (0..n).any(|j| self.lu[j][j] == 0.0) {
            return Err("Matrix is singular.".to_string());
        }
        let nx = b.columns();
        let mut x: Vec<Vec<f64>> = self.pivots.iter().map(|p| b.a[*p].clone()).collect();
        for k in 0..n {
            for i in (k + 1)..n {
                for j in 0..nx {
                    x[i][j] -= x[k][j] * self.lu[i][k];
                }
            }
        }
        for k in (0..n).rev() {
            for j in 0..nx {
                x[k][j] /= self.lu[k][k];
            }
            for i in 0..k {
                for j in 0..nx {
                    x[i][j] -= x[k][j] * self.lu[i][k];
                }
            }
        }
        Ok(Matrix { a: x })
    }
}

// ================================================================================================
// The data.
// ================================================================================================

/// `VariantDatum`.
#[derive(Debug, Clone, PartialEq)]
pub struct Datum {
    pub annotations: Vec<f64>,
    pub is_null: Vec<bool>,
    pub is_known: bool,
    pub lod: f64,
    pub at_truth_site: bool,
    pub at_training_site: bool,
    pub at_anti_training_site: bool,
    pub is_transition: bool,
    pub is_snp: bool,
    pub failing_std_threshold: bool,
    pub prior: f64,
    /// `loc`: contig, start and end. `None` for an aggregate datum.
    pub loc: Option<(String, i64, i64)>,
    pub worst_annotation: i32,
    pub worst_value: f64,
    pub is_aggregate: bool,
    /// The alleles an allele-specific datum stands for, as display strings.
    pub reference_allele: Option<String>,
    pub alternate_allele: Option<String>,
}

impl Datum {
    pub fn new(annotations: Vec<f64>, is_null: Vec<bool>) -> Self {
        Self {
            annotations,
            is_null,
            is_known: false,
            lod: 0.0,
            at_truth_site: false,
            at_training_site: false,
            at_anti_training_site: false,
            is_transition: false,
            is_snp: false,
            failing_std_threshold: false,
            prior: 0.0,
            loc: None,
            worst_annotation: 0,
            worst_value: 0.0,
            is_aggregate: false,
            reference_allele: None,
            alternate_allele: None,
        }
    }
}

/// `VariantRecalibratorArgumentCollection`, the model's half of the arguments.
#[derive(Debug, Clone, PartialEq)]
pub struct ModelArguments {
    pub use_as_annotations: bool,
    pub max_gaussians: i32,
    pub max_negative_gaussians: i32,
    pub max_iterations: i32,
    pub k_means_iterations: i32,
    pub std_threshold: f64,
    pub shrinkage: f64,
    pub dirichlet: f64,
    pub prior_counts: f64,
    pub max_training: i32,
    pub min_bad_variants: i32,
    pub bad_lod_cutoff: f64,
    pub mq_cap: i32,
    pub mq_jitter: f64,
}

/// `TrainingSet`, read off a resource's tags.
#[derive(Debug, Clone, PartialEq)]
pub struct TrainingSet {
    pub is_known: bool,
    pub is_training: bool,
    pub is_anti_training: bool,
    pub is_truth: bool,
    pub prior: f64,
}

impl TrainingSet {
    /// The tags a `--resource:name,key=value` carries. A boolean is set only by the exact string
    /// `true`, and a prior that is not a number is refused.
    pub fn from_attributes(attributes: &[(String, String)]) -> Result<TrainingSet, ModelError> {
        let get = |key: &str| {
            attributes
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value.as_str())
        };
        let flag = |key: &str| get(key) == Some("true");
        let prior = match get("prior") {
            None => 0.0,
            Some(text) => java_parse_double(text).ok_or_else(|| {
                ModelError::new(
                    "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile",
                    "Unknown file is malformed: Malformed floating point valueprior",
                )
            })?,
        };
        Ok(TrainingSet {
            is_known: flag("known"),
            is_training: flag("training"),
            is_anti_training: flag("bad"),
            is_truth: flag("truth"),
            prior,
        })
    }
}

/// `Double.valueOf(String)` for the forms a VCF or a tag writes: surrounding blanks are trimmed,
/// and a trailing `d` or `f` is Java's own suffix.
pub fn java_parse_double(text: &str) -> Option<f64> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    let body = trimmed
        .strip_suffix(['d', 'D', 'f', 'F'])
        .unwrap_or(trimmed);
    match body {
        "NaN" | "+NaN" | "-NaN" => return Some(f64::NAN),
        "Infinity" | "+Infinity" => return Some(f64::INFINITY),
        "-Infinity" => return Some(f64::NEG_INFINITY),
        _ => {}
    }
    if body.is_empty()
        || body.contains(|c: char| c.is_ascii_alphabetic() && c != 'e' && c != 'E')
        || body.eq_ignore_ascii_case("inf")
    {
        return None;
    }
    body.parse::<f64>().ok()
}

/// What `decodeAnnotation` reads off a record for one key.
pub enum AttributeValue<'a> {
    Absent,
    /// A single value, as the VCF wrote it.
    One(&'a str),
    /// Several comma-separated values.
    Many(Vec<&'a str>),
}

/// The one-in-a-hundred jitter keys `decodeAnnotation` names, compared ignoring case.
const ZERO_JITTER_KEYS: [&str; 4] = ["HaplotypeScore", "FS", "AS_FilterStatus", "InbreedingCoeff"];

/// `VariantDataManager.logitTransform`.
fn logit_transform(x: f64, xmin: f64, xmax: f64) -> f64 {
    ln((x - xmin) / (xmax - x))
}

/// `VariantDataManager.decodeAnnotation`: the value, NaN when it is missing or not a number, and
/// then the jitter the key asks for, drawn from `random`.
///
/// `allele_index` is the alternate an allele-specific datum stands for, which is read out of an
/// `AS_` annotation's list. A list read as one value is `getAttributeAsDouble`'s
/// `ClassCastException`, which escapes the tool's own `NumberFormatException` handler.
pub fn decode_annotation(
    key: &str,
    value: AttributeValue,
    allele_index: Option<usize>,
    arguments: &ModelArguments,
    random: &mut JavaRandom,
) -> Result<f64, ModelError> {
    // The reference's own seven-digit literal, which is what the comparison is against.
    #[allow(clippy::approx_constant)]
    const LOG_OF_TWO: f64 = 0.6931472;
    const SAFETY_OFFSET: f64 = 0.01;
    const PRECISION: f64 = 0.01;
    let parse = |text: &str| java_parse_double(text).unwrap_or(f64::NAN);
    let mut value = if arguments.use_as_annotations && key.starts_with("AS_") {
        let list: Vec<&str> = match value {
            AttributeValue::Absent => Vec::new(),
            AttributeValue::One(text) => vec![text],
            AttributeValue::Many(list) => list,
        };
        let index = allele_index.unwrap_or(0);
        match list.get(index) {
            Some(text) => parse(text),
            None => {
                return Err(ModelError::new(
                    "java.lang.IndexOutOfBoundsException",
                    format!("Index {index} out of bounds for length {}", list.len()),
                ))
            }
        }
    } else {
        match value {
            AttributeValue::Absent => f64::NAN,
            AttributeValue::One(text) => parse(text),
            AttributeValue::Many(list) if list.len() == 1 => parse(list[0]),
            AttributeValue::Many(_) => {
                return Err(ModelError::new(
                    "java.lang.ClassCastException",
                    "class java.util.ArrayList cannot be cast to class java.lang.String \
                     (java.util.ArrayList and java.lang.String are in module java.base of loader \
                     'bootstrap')",
                ))
            }
        }
    };
    if value.is_infinite() {
        value = f64::NAN;
    }
    let is = |name: &str| key.eq_ignore_ascii_case(name);
    if is(ZERO_JITTER_KEYS[0]) && close(value, 0.0, PRECISION) {
        value += 0.01 * random.next_gaussian();
    }
    if (is(ZERO_JITTER_KEYS[1]) || is(ZERO_JITTER_KEYS[2])) && close(value, 0.0, PRECISION) {
        value += 0.01 * random.next_gaussian();
    }
    if is(ZERO_JITTER_KEYS[3]) && close(value, 0.0, PRECISION) {
        value += 0.01 * random.next_gaussian();
    }
    if (is("SOR") || is("AS_SOR")) && close(value, LOG_OF_TWO, PRECISION) {
        value += 0.01 * random.next_gaussian();
    }
    if is("MQ") {
        if arguments.mq_cap > 0 {
            let cap = arguments.mq_cap as f64;
            value = logit_transform(value, -SAFETY_OFFSET, cap + SAFETY_OFFSET);
            if close(
                value,
                logit_transform(cap, -SAFETY_OFFSET, cap + SAFETY_OFFSET),
                PRECISION,
            ) {
                value += arguments.mq_jitter * random.next_gaussian();
            }
        } else if close(value, arguments.mq_cap as f64, PRECISION) {
            value += arguments.mq_jitter * random.next_gaussian();
        }
    }
    if is("AS_MQ") {
        value += arguments.mq_jitter * random.next_gaussian();
    }
    Ok(value)
}

/// `QualityUtils.qualToProb`, and the log-odds prior `addDatum` turns it into.
pub fn prior_log_odds(prior: f64) -> f64 {
    let factor = 1.0 - pow(10.0, prior / -10.0);
    log10(factor) - log10(1.0 - factor)
}

/// `VariantDataManager`: the data, the annotation names in their current order, and the
/// normalisation.
///
/// Every datum ever built stays in `data`, and `active` is the reference's `data` LIST: the
/// indices it holds, in its order. A datum dropped from the list (an aggregate one, once the models
/// are trained) is still reachable from another list that held it, as a Java reference is.
#[derive(Debug, Clone)]
pub struct DataManager {
    pub data: Vec<Datum>,
    pub active: Vec<usize>,
    pub annotation_keys: Vec<String>,
    pub mean_vector: Vec<f64>,
    pub variance_vector: Vec<f64>,
}

impl DataManager {
    /// Duplicates are dropped, the first occurrence keeping its place.
    pub fn new(keys: &[String]) -> Self {
        let mut unique: Vec<String> = Vec::new();
        for key in keys {
            if !unique.contains(key) {
                unique.push(key.clone());
            }
        }
        let n = unique.len();
        Self {
            data: Vec::new(),
            active: Vec::new(),
            annotation_keys: unique,
            mean_vector: vec![0.0; n],
            variance_vector: vec![0.0; n],
        }
    }

    /// `setNormalization`, from a model report's two tables, keyed by annotation name.
    pub fn set_normalization(&mut self, means: &[(String, f64)], stdevs: &[(String, f64)]) {
        let lookup = |table: &[(String, f64)], key: &str| {
            table
                .iter()
                .rev()
                .find(|(name, _)| name == key)
                .map_or(f64::NAN, |(_, value)| *value)
        };
        for i in 0..self.annotation_keys.len() {
            self.mean_vector[i] = lookup(means, &self.annotation_keys[i]);
            self.variance_vector[i] = lookup(stdevs, &self.annotation_keys[i]);
        }
    }

    fn mean(&self, index: usize, training: bool) -> f64 {
        let (mut sum, mut n) = (0.0, 0i32);
        for datum in self.active.iter().map(|i| &self.data[*i]) {
            if training == datum.at_training_site && !datum.is_null[index] {
                sum += datum.annotations[index];
                n += 1;
            }
        }
        sum / n as f64
    }

    fn standard_deviation(&self, mean: f64, index: usize, training: bool) -> f64 {
        let (mut sum, mut n) = (0.0, 0i32);
        for datum in self.active.iter().map(|i| &self.data[*i]) {
            if training == datum.at_training_site && !datum.is_null[index] {
                let d = datum.annotations[index] - mean;
                sum += d * d;
                n += 1;
            }
        }
        (sum / n as f64).sqrt()
    }

    /// `normalizeData`: every annotation to mean zero and deviation one over the training data, a
    /// missing one to a tenth of a Gaussian draw, the outliers marked, and the annotations put in
    /// `order` or, without one, in decreasing distance of the training mean from the rest.
    pub fn normalize(
        &mut self,
        calculate_means: bool,
        order: Option<&[usize]>,
        std_threshold: f64,
        random: &mut JavaRandom,
    ) -> Result<(), ModelError> {
        let mut zero_variance = false;
        for i in 0..self.mean_vector.len() {
            let (mean, std) = if calculate_means {
                let mean = self.mean(i, true);
                let std = self.standard_deviation(mean, i, true);
                if mean.is_nan() {
                    return Err(ModelError::new(
                        BAD_INPUT,
                        format!(
                            "Bad input: Values for {} annotation not detected for ANY training \
                             variant in the input callset. VariantAnnotator may be used to add \
                             these annotations.",
                            self.annotation_keys[i]
                        ),
                    ));
                }
                zero_variance = zero_variance || std < 1e-5;
                self.mean_vector[i] = mean;
                self.variance_vector[i] = std;
                (mean, std)
            } else {
                (self.mean_vector[i], self.variance_vector[i])
            };
            for &k in &self.active {
                let datum = &mut self.data[k];
                datum.annotations[i] = if datum.is_null[i] {
                    0.1 * random.next_gaussian()
                } else {
                    (datum.annotations[i] - mean) / std
                };
            }
        }
        if zero_variance {
            return Err(ModelError::new(
                BAD_INPUT,
                "Bad input: Found annotations with zero variance. They must be excluded before \
                 proceeding.",
            ));
        }
        for &k in &self.active {
            let datum = &mut self.data[k];
            let mut remove = false;
            for value in &datum.annotations {
                remove = remove || value.abs() > std_threshold;
            }
            datum.failing_std_threshold = remove;
        }
        let order: Vec<usize> = match order {
            Some(order) => order.to_vec(),
            None => {
                let mut keyed: Vec<(f64, usize)> = (0..self.mean_vector.len())
                    .map(|i| (-1.0 * (self.mean_vector[i] - self.mean(i, false)).abs(), i))
                    .collect();
                keyed.sort_by(|a, b| double_compare(a.0, b.0));
                keyed.into_iter().map(|(_, i)| i).collect()
            }
        };
        let reorder = |values: &[f64]| order.iter().map(|i| values[*i]).collect::<Vec<f64>>();
        self.annotation_keys = order
            .iter()
            .map(|i| self.annotation_keys[*i].clone())
            .collect();
        self.variance_vector = reorder(&self.variance_vector);
        self.mean_vector = reorder(&self.mean_vector);
        for &k in &self.active {
            let datum = &mut self.data[k];
            datum.annotations = order.iter().map(|i| datum.annotations[*i]).collect();
            datum.is_null = order.iter().map(|i| datum.is_null[*i]).collect();
        }
        Ok(())
    }

    /// `getTrainingData`: the training sites inside the deviation threshold, shuffled and cut when
    /// there are more than the maximum and at least the minimum of bad variants.
    pub fn training_data(&self, arguments: &ModelArguments, random: &mut JavaRandom) -> Vec<usize> {
        let mut training: Vec<usize> = self
            .active
            .iter()
            .copied()
            .filter(|i| self.data[*i].at_training_site && !self.data[*i].failing_std_threshold)
            .collect();
        if (training.len() as i64) < arguments.min_bad_variants as i64 {
            return training;
        }
        if training.len() as i64 > arguments.max_training as i64 {
            shuffle(&mut training, random);
            training.truncate(arguments.max_training.max(0) as usize);
        }
        training
    }

    /// `selectWorstVariants`, which marks what it selects.
    pub fn select_worst(&mut self, cutoff: f64) -> Vec<usize> {
        let mut worst = Vec::new();
        for &i in &self.active {
            let datum = &mut self.data[i];
            if !datum.failing_std_threshold && !datum.lod.is_infinite() && datum.lod < cutoff {
                datum.at_anti_training_site = true;
                worst.push(i);
            }
        }
        worst
    }

    /// `getEvaluationData`.
    pub fn evaluation_data(&self) -> Vec<usize> {
        self.active
            .iter()
            .copied()
            .filter(|i| {
                let d = &self.data[*i];
                !d.failing_std_threshold && !d.at_training_site && !d.at_anti_training_site
            })
            .collect()
    }

    /// `denormalizeDatum`.
    pub fn denormalize(&self, value: f64, index: usize) -> f64 {
        value * self.variance_vector[index] + self.mean_vector[index]
    }
}

/// `getRandomDataForPlotting`: each list shuffled, the first `n` of each, and those shuffled.
pub fn random_data_for_plotting(
    n: usize,
    training: &mut [usize],
    anti_training: &mut [usize],
    evaluation: &mut [usize],
    random: &mut JavaRandom,
) -> Vec<usize> {
    shuffle(training, random);
    shuffle(anti_training, random);
    shuffle(evaluation, random);
    let mut out: Vec<usize> = Vec::new();
    out.extend_from_slice(&training[..n.min(training.len())]);
    out.extend_from_slice(&anti_training[..n.min(anti_training.len())]);
    out.extend_from_slice(&evaluation[..n.min(evaluation.len())]);
    shuffle(&mut out, random);
    out
}

// ================================================================================================
// The mixture.
// ================================================================================================

/// `MultivariateGaussian`.
#[derive(Debug, Clone, PartialEq)]
pub struct Gaussian {
    pub p_mixture_log10: f64,
    pub sum_prob: f64,
    pub mu: Vec<f64>,
    pub sigma: Matrix,
    pub hyper_a: f64,
    pub hyper_b: f64,
    pub hyper_lambda: f64,
    cached_denom_log10: f64,
    cached_sigma_inverse: Option<Matrix>,
    p_var: Vec<f64>,
    p_var_index: usize,
}

impl Gaussian {
    pub fn new(num_variants: usize, num_annotations: usize) -> Self {
        Self {
            p_mixture_log10: 0.0,
            sum_prob: 0.0,
            mu: vec![0.0; num_annotations],
            sigma: Matrix::zeros(num_annotations, num_annotations),
            hyper_a: 0.0,
            hyper_b: 0.0,
            hyper_lambda: 0.0,
            cached_denom_log10: 0.0,
            cached_sigma_inverse: None,
            p_var: vec![0.0; num_variants],
            p_var_index: 0,
        }
    }

    fn initialize_random_mu(&mut self, random: &mut JavaRandom) {
        for value in &mut self.mu {
            *value = -4.0 + 8.0 * random.next_double();
        }
    }

    fn initialize_random_sigma(&mut self, random: &mut JavaRandom) {
        let n = self.mu.len();
        let mut rand = Matrix::zeros(n, n);
        for i in 0..n {
            for j in i..n {
                rand.a[j][i] = 0.55 + 1.25 * random.next_double();
                if random.next_boolean() {
                    rand.a[j][i] *= -1.0;
                }
                if i != j {
                    rand.a[i][j] = 0.0;
                }
            }
        }
        self.sigma = rand.times(&rand.transpose());
    }

    fn precompute_inverse(&mut self) -> Result<(), ModelError> {
        match self.sigma.inverse() {
            Ok(inverse) => {
                self.cached_sigma_inverse = Some(inverse);
                Ok(())
            }
            Err(_) => Err(ModelError::new(
                USER_EXCEPTION,
                "Error during clustering. Most likely there are too few variants used during \
                 Gaussian mixture modeling. Please consider raising the number of variants used \
                 to train the negative model (via --percentBadVariants 0.05, for example) or \
                 lowering the maximum number of Gaussians to use in the model (via --maxGaussians \
                 4, for example).",
            )),
        }
    }

    fn precompute_denominator_for_evaluation(&mut self) -> Result<(), ModelError> {
        if self.p_mixture_log10 == f64::NEG_INFINITY {
            return Ok(());
        }
        self.precompute_inverse()?;
        let d = self.mu.len() as f64;
        self.cached_denom_log10 = log10(pow(2.0 * std::f64::consts::PI, -1.0 * d / 2.0))
            + log10(pow(self.sigma.det(), -0.5));
        let det = self.sigma.det();
        if self.cached_denom_log10.is_nan() || det < EPSILON {
            return Err(ModelError::new(
                GATK_EXCEPTION,
                format!(
                    "Denominator for gaussian evaluation cannot be computed. Covariance \
                     determinant is {}. One or more annotations (usually MQ) may have \
                     insufficient variance.",
                    java_double_string(det)
                ),
            ));
        }
        Ok(())
    }

    fn precompute_denominator_for_variational_bayes(
        &mut self,
        sum_lambda: f64,
    ) -> Result<(), ModelError> {
        self.precompute_inverse()?;
        let a = self.hyper_a;
        if let Some(inverse) = &mut self.cached_sigma_inverse {
            inverse.times_equals(a);
        }
        let n = self.mu.len();
        let mut sum = 0.0;
        for j in 1..=n {
            sum += digamma((a + 1.0 - j as f64) / 2.0);
        }
        sum -= ln(self.sigma.det());
        sum += ln(2.0) * n as f64;
        let lambda = 0.5 * sum;
        let pi = digamma(self.hyper_lambda) - digamma(sum_lambda);
        let beta = (-1.0 * n as f64) / (2.0 * self.hyper_b);
        let ln10 = ln(10.0);
        self.cached_denom_log10 = (pi / ln10) + (lambda / ln10) + (beta / ln10);
        Ok(())
    }

    fn evaluate_datum_log10(&self, annotations: &[f64]) -> f64 {
        if self.p_mixture_log10 == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        let n = self.mu.len();
        let inverse = self
            .cached_sigma_inverse
            .as_ref()
            .expect("a precomputed inverse");
        let mut cross = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                cross[i] += (annotations[j] - self.mu[j]) * inverse.a[j][i];
            }
        }
        let mut kernel = 0.0;
        for i in 0..n {
            kernel += cross[i] * (annotations[i] - self.mu[i]);
        }
        ((-0.5 * kernel) / ln(10.0)) + self.cached_denom_log10
    }

    fn maximize(
        &mut self,
        data: &[&Datum],
        empirical_mu: &[f64],
        empirical_sigma: &Matrix,
        shrinkage: f64,
        dirichlet: f64,
        degrees_of_freedom: f64,
    ) {
        let n = self.mu.len();
        self.sum_prob = 1e-10;
        let mut wishart = Matrix::zeros(n, n);
        self.mu.iter_mut().for_each(|x| *x = 0.0);
        self.sigma = Matrix::zeros(n, n);
        for (index, datum) in data.iter().enumerate() {
            let prob = self.p_var[index];
            self.sum_prob += prob;
            for j in 0..n {
                self.mu[j] += prob * datum.annotations[j];
            }
        }
        for j in 0..n {
            self.mu[j] /= self.sum_prob;
        }
        let factor = (shrinkage * self.sum_prob) / (shrinkage + self.sum_prob);
        for i in 0..n {
            let delta = factor * (self.mu[i] - empirical_mu[i]);
            for j in 0..n {
                wishart.a[i][j] = delta * (self.mu[j] - empirical_mu[j]);
            }
        }
        self.accumulate_sigma(data);
        self.sigma.plus_equals(empirical_sigma);
        self.sigma.plus_equals(&wishart);
        for i in 0..n {
            self.mu[i] = (self.sum_prob * self.mu[i] + shrinkage * empirical_mu[i])
                / (self.sum_prob + shrinkage);
        }
        self.hyper_a = self.sum_prob + degrees_of_freedom;
        self.hyper_b = self.sum_prob + shrinkage;
        self.hyper_lambda = self.sum_prob + dirichlet;
        self.reset_p_var();
    }

    /// The weighted scatter matrix, one datum's outer product at a time, each regularised.
    fn accumulate_sigma(&mut self, data: &[&Datum]) {
        let n = self.mu.len();
        let mut term = Matrix::zeros(n, n);
        for (index, datum) in data.iter().enumerate() {
            let prob = self.p_var[index];
            for i in 0..n {
                for j in 0..n {
                    let reg = if i == j {
                        COVARIANCE_REGULARIZATION_EPSILON
                    } else {
                        0.0
                    };
                    term.a[i][j] = prob
                        * (datum.annotations[i] - self.mu[i])
                        * (datum.annotations[j] - self.mu[j])
                        + reg;
                }
            }
            self.sigma.plus_equals(&term);
        }
    }

    fn evaluate_final_parameters(&mut self, data: &[&Datum]) {
        let n = self.mu.len();
        self.sum_prob = 0.0;
        self.mu.iter_mut().for_each(|x| *x = 0.0);
        self.sigma = Matrix::zeros(n, n);
        for (index, datum) in data.iter().enumerate() {
            let prob = self.p_var[index];
            self.sum_prob += prob;
            for j in 0..n {
                self.mu[j] += prob * datum.annotations[j];
            }
        }
        for j in 0..n {
            self.mu[j] /= self.sum_prob;
        }
        self.accumulate_sigma(data);
        self.sigma.times_equals(1.0 / self.sum_prob);
        self.reset_p_var();
    }

    fn reset_p_var(&mut self) {
        self.p_var.iter_mut().for_each(|x| *x = 0.0);
        self.p_var_index = 0;
    }
}

/// `Double.toString`, for the one message that prints a determinant.
pub fn java_double_string(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// `GaussianMixtureModel`.
#[derive(Debug, Clone, PartialEq)]
pub struct Mixture {
    pub gaussians: Vec<Gaussian>,
    shrinkage: f64,
    dirichlet: f64,
    prior_counts: f64,
    empirical_mu: Vec<f64>,
    empirical_sigma: Matrix,
    pub ready_for_evaluation: bool,
    pub failed_to_converge: bool,
}

impl Mixture {
    pub fn new(
        num_gaussians: usize,
        num_variants: usize,
        num_annotations: usize,
        arguments: &ModelArguments,
    ) -> Self {
        let gaussians = (0..num_gaussians)
            .map(|_| Gaussian::new(num_variants, num_annotations))
            .collect();
        Self::with_gaussians(gaussians, num_annotations, arguments)
    }

    /// The model-report constructor: the Gaussians given, the hyperparameters left as they are.
    pub fn with_gaussians(
        gaussians: Vec<Gaussian>,
        num_annotations: usize,
        arguments: &ModelArguments,
    ) -> Self {
        let empirical_sigma = Matrix::identity(num_annotations)
            .times_scalar(200.0)
            .inverse()
            .unwrap_or_else(|_| Matrix::zeros(num_annotations, num_annotations));
        Self {
            gaussians,
            shrinkage: arguments.shrinkage,
            dirichlet: arguments.dirichlet,
            prior_counts: arguments.prior_counts,
            empirical_mu: vec![0.0; num_annotations],
            empirical_sigma,
            ready_for_evaluation: false,
            failed_to_converge: false,
        }
    }

    fn initialize_random_model(
        &mut self,
        data: &[&Datum],
        k_means_iterations: i32,
        random: &mut JavaRandom,
    ) -> Result<(), ModelError> {
        for gaussian in &mut self.gaussians {
            gaussian.initialize_random_mu(random);
        }
        self.k_means(data, k_means_iterations, random)?;
        let count = self.gaussians.len() as f64;
        for gaussian in &mut self.gaussians {
            gaussian.p_mixture_log10 = log10(1.0 / count);
            gaussian.sum_prob = 1.0 / count;
            gaussian.initialize_random_sigma(random);
            gaussian.hyper_a = self.prior_counts;
            gaussian.hyper_b = self.shrinkage;
            gaussian.hyper_lambda = self.dirichlet;
        }
        Ok(())
    }

    fn k_means(
        &mut self,
        data: &[&Datum],
        iterations: i32,
        random: &mut JavaRandom,
    ) -> Result<(), ModelError> {
        let mut assignment: Vec<Option<usize>> = vec![None; data.len()];
        let mut t = 0;
        while t < iterations {
            t += 1;
            for (index, datum) in data.iter().enumerate() {
                let mut min_distance = f64::MAX;
                let mut min_gaussian = None;
                for (g, gaussian) in self.gaussians.iter().enumerate() {
                    let distance = distance_squared(&datum.annotations, &gaussian.mu);
                    if distance < min_distance {
                        min_distance = distance;
                        min_gaussian = Some(g);
                    }
                }
                assignment[index] = min_gaussian;
            }
            for g in 0..self.gaussians.len() {
                let gaussian = &mut self.gaussians[g];
                gaussian.mu.iter_mut().for_each(|x| *x = 0.0);
                let mut assigned = 0i32;
                for (index, datum) in data.iter().enumerate() {
                    match assignment[index] {
                        None => {
                            return Err(ModelError::new(
                                NULL_POINTER,
                                "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.vqsr.\
                                 MultivariateGaussian.equals(Object)\" because \"datum.assignment\" \
                                 is null",
                            ))
                        }
                        Some(owner) if owner == g => {
                            assigned += 1;
                            for j in 0..gaussian.mu.len() {
                                gaussian.mu[j] += 1.0 * datum.annotations[j];
                            }
                        }
                        Some(_) => {}
                    }
                }
                if assigned != 0 {
                    for x in &mut gaussian.mu {
                        *x /= assigned as f64;
                    }
                } else {
                    gaussian.initialize_random_mu(random);
                }
            }
        }
        Ok(())
    }

    fn sum_hyper_lambda(&self) -> f64 {
        stream_sum(self.gaussians.iter().map(|g| g.hyper_lambda))
    }

    fn expectation_step(&mut self, data: &[&Datum]) -> Result<(), ModelError> {
        for g in 0..self.gaussians.len() {
            let sum_lambda = self.sum_hyper_lambda();
            self.gaussians[g].precompute_denominator_for_variational_bayes(sum_lambda)?;
        }
        for datum in data {
            let log10s: Vec<f64> = self
                .gaussians
                .iter()
                .map(|g| g.evaluate_datum_log10(&datum.annotations))
                .collect();
            let normalized = normalize_to_linear(&log10s);
            for (gaussian, p) in self.gaussians.iter_mut().zip(normalized) {
                let index = gaussian.p_var_index;
                gaussian.p_var[index] = p;
                gaussian.p_var_index += 1;
            }
        }
        Ok(())
    }

    fn maximization_step(&mut self, data: &[&Datum]) {
        let (mu, sigma) = (self.empirical_mu.clone(), self.empirical_sigma.clone());
        for gaussian in &mut self.gaussians {
            gaussian.maximize(
                data,
                &mu,
                &sigma,
                self.shrinkage,
                self.dirichlet,
                self.prior_counts,
            );
        }
    }

    /// `normalizePMixtureLog10`: the mixture weights renormalised, and how far they moved.
    pub fn normalize_p_mixture_log10(&mut self) -> f64 {
        let sum_pk = stream_sum(self.gaussians.iter().map(|g| g.sum_prob));
        let log10_sum_pk = log10(sum_pk);
        let mut weights: Vec<f64> = self
            .gaussians
            .iter()
            .map(|g| log10(g.sum_prob) - log10_sum_pk)
            .collect();
        normalize_log10_in_place(&mut weights);
        let mut sum_diff = 0.0;
        for (gaussian, weight) in self.gaussians.iter_mut().zip(weights) {
            sum_diff += (weight - gaussian.p_mixture_log10).abs();
            gaussian.p_mixture_log10 = weight;
        }
        sum_diff
    }

    fn evaluate_final_parameters(&mut self, data: &[&Datum]) {
        for gaussian in &mut self.gaussians {
            gaussian.evaluate_final_parameters(data);
        }
        self.normalize_p_mixture_log10();
    }

    fn precompute_for_evaluation(&mut self) -> Result<(), ModelError> {
        for gaussian in &mut self.gaussians {
            gaussian.precompute_denominator_for_evaluation()?;
        }
        self.ready_for_evaluation = true;
        Ok(())
    }

    /// `evaluateDatum`, which marginalises when any dimension is missing and so may draw.
    pub fn evaluate_datum(
        &self,
        datum: &mut Datum,
        random: &mut JavaRandom,
    ) -> Result<f64, ModelError> {
        if datum.is_null.iter().any(|null| *null) {
            return self.evaluate_datum_marginalized(datum, random);
        }
        let values: Vec<f64> = self
            .gaussians
            .iter()
            .map(|g| g.p_mixture_log10 + g.evaluate_datum_log10(&datum.annotations))
            .collect();
        nan_tolerant_log10_sum_log10(&values)
    }

    /// `evaluateDatumMarginalized`: twenty draws per missing dimension, WRITTEN INTO the datum, and
    /// the average of the linear probabilities.
    fn evaluate_datum_marginalized(
        &self,
        datum: &mut Datum,
        random: &mut JavaRandom,
    ) -> Result<f64, ModelError> {
        let mut draws = 0i32;
        let mut sum = 0.0;
        for i in 0..datum.annotations.len() {
            if datum.is_null[i] {
                for _ in 0..20 {
                    datum.annotations[i] = random.next_gaussian();
                    let values: Vec<f64> = self
                        .gaussians
                        .iter()
                        .map(|g| g.p_mixture_log10 + g.evaluate_datum_log10(&datum.annotations))
                        .collect();
                    sum += pow(10.0, nan_tolerant_log10_sum_log10(&values)?);
                    draws += 1;
                }
            }
        }
        Ok(log10(sum / draws as f64))
    }

    /// `evaluateDatumInOneDimension`, which reads the covariance's diagonal as a DEVIATION.
    fn evaluate_in_one_dimension(
        &self,
        datum: &Datum,
        i: usize,
    ) -> Result<Option<f64>, ModelError> {
        if datum.is_null[i] {
            return Ok(None);
        }
        let mut values = Vec::with_capacity(self.gaussians.len());
        for gaussian in &self.gaussians {
            let mut value = gaussian.p_mixture_log10;
            if gaussian.p_mixture_log10 != f64::NEG_INFINITY {
                value += normal_distribution_log10(
                    gaussian.mu[i],
                    gaussian.sigma.get(i, i),
                    datum.annotations[i],
                )?;
            }
            values.push(value);
        }
        nan_tolerant_log10_sum_log10(&values).map(Some)
    }
}

// ================================================================================================
// The engine.
// ================================================================================================

/// `VariantRecalibratorEngine.generateModel`: a mixture of `max_gaussians`, trained by VBEM.
pub fn generate_model(
    data: &[&Datum],
    max_gaussians: i32,
    arguments: &ModelArguments,
    random: &mut JavaRandom,
) -> Result<Mixture, ModelError> {
    if data.is_empty() {
        return Err(ModelError::new(NEGATIVE_FAILURE, "No data found."));
    }
    if max_gaussians <= 0 {
        return Err(ModelError::new(
            ILLEGAL_ARGUMENT,
            format!("maxGaussians must be a positive integer but found: {max_gaussians}"),
        ));
    }
    let mut model = Mixture::new(
        max_gaussians as usize,
        data.len(),
        data[0].annotations.len(),
        arguments,
    );
    model.initialize_random_model(data, arguments.k_means_iterations, random)?;
    model.normalize_p_mixture_log10();
    model.expectation_step(data)?;
    let mut iteration = 0;
    while iteration < arguments.max_iterations {
        iteration += 1;
        model.maximization_step(data);
        let change = model.normalize_p_mixture_log10();
        model.expectation_step(data)?;
        if iteration > 2 && change < MIN_PROB_CONVERGENCE {
            break;
        }
    }
    model.evaluate_final_parameters(data);
    Ok(model)
}

/// `VariantRecalibratorEngine.evaluateData` over the data at `indices` (every datum when `None`).
///
/// A model whose denominators cannot be computed, or a datum scored NaN, marks the model failed
/// and stops: the data after it keep whatever score they had.
pub fn evaluate_data(
    data: &mut [Datum],
    indices: Option<&[usize]>,
    model: &mut Mixture,
    contrastive: bool,
    random: &mut JavaRandom,
) -> Result<(), ModelError> {
    if !model.ready_for_evaluation && model.precompute_for_evaluation().is_err() {
        model.failed_to_converge = true;
        return Ok(());
    }
    let all: Vec<usize>;
    let indices = match indices {
        Some(indices) => indices,
        None => {
            all = (0..data.len()).collect();
            &all
        }
    };
    for &i in indices {
        let lod = model.evaluate_datum(&mut data[i], random)?;
        if lod.is_nan() {
            model.failed_to_converge = true;
            return Ok(());
        }
        let datum = &mut data[i];
        datum.lod = if contrastive {
            if datum.lod.is_infinite() {
                MIN_ACCEPTABLE_LOD_SCORE + random.next_double() * MIN_ACCEPTABLE_LOD_SCORE
            } else {
                datum.prior + datum.lod - lod
            }
        } else {
            lod
        };
    }
    Ok(())
}

/// `calculateWorstPerformingAnnotation`: the culprit of every datum.
pub fn worst_performing_annotation(
    data: &mut [Datum],
    indices: &[usize],
    good: &Mixture,
    bad: &Mixture,
) -> Result<(), ModelError> {
    for &index in indices {
        let datum = &mut data[index];
        let mut worst = -1i32;
        let mut min_prob = f64::MAX;
        let mut worst_value = -1.0;
        for i in 0..datum.annotations.len() {
            let good_prob = good.evaluate_in_one_dimension(datum, i)?;
            let bad_prob = bad.evaluate_in_one_dimension(datum, i)?;
            if let (Some(good_prob), Some(bad_prob)) = (good_prob, bad_prob) {
                let prob = good_prob - bad_prob;
                if prob < min_prob {
                    min_prob = prob;
                    worst = i as i32;
                    worst_value = datum.annotations[i];
                }
            }
        }
        datum.worst_annotation = worst;
        datum.worst_value = worst_value;
    }
    Ok(())
}

// ================================================================================================
// The model report.
// ================================================================================================

/// `writeModelReport`: the normalisation, the mixture weights and the two models' means and
/// covariances, every number at `%.16E`.
pub fn model_report(
    manager: &DataManager,
    good: &Mixture,
    bad: Option<&Mixture>,
    annotations: &[String],
) -> String {
    use gatk_engine::gatk_report::{Report, Sorting, Table, Value};
    const FORMAT: &str = "%.16E";
    let vector_table = |name: &str,
                        description: &str,
                        keys: &[String],
                        values: &[f64],
                        column: &str,
                        first: &str| {
        let mut table = Table::new(name, description, Sorting::DoNotSort);
        table.add_column(first, "");
        table.add_column(column, FORMAT);
        for (i, value) in values.iter().enumerate() {
            let key = i.to_string();
            table.set(&key, first, Value::Str(keys[i].clone()));
            table.set(&key, column, Value::Double(*value));
        }
        table
    };
    let mut report = Report::new();
    report.add_table(vector_table(
        "AnnotationMeans",
        "Mean for each annotation, used to normalize data",
        &manager.annotation_keys,
        &manager.mean_vector,
        "Mean",
        "Annotation",
    ));
    report.add_table(vector_table(
        "AnnotationStdevs",
        "Standard deviation for each annotation, used to normalize data",
        &manager.annotation_keys,
        &manager.variance_vector,
        "StandardDeviation",
        "Annotation",
    ));
    let pmix = |name: &str, model: &Mixture| {
        let names: Vec<String> = (0..model.gaussians.len()).map(|i| i.to_string()).collect();
        let values: Vec<f64> = model.gaussians.iter().map(|g| g.p_mixture_log10).collect();
        let mut table = Table::new(
            name,
            "Pmixture log 10 used to evaluate model",
            Sorting::DoNotSort,
        );
        table.add_column("Gaussian", "");
        table.add_column("pMixLog10", FORMAT);
        for (i, value) in values.iter().enumerate() {
            table.set(&names[i], "Gaussian", Value::Str(names[i].clone()));
            table.set(&names[i], "pMixLog10", Value::Double(*value));
        }
        table
    };
    report.add_table(pmix("GoodGaussianPMix", good));
    if let Some(bad) = bad {
        report.add_table(pmix("BadGaussianPMix", bad));
    }
    let means = |name: &str, description: &str, model: &Mixture| {
        let mut table = Table::new(name, description, Sorting::DoNotSort);
        table.add_column("Gaussian", "");
        for annotation in annotations {
            table.add_column(annotation, FORMAT);
        }
        for (i, gaussian) in model.gaussians.iter().enumerate() {
            let key = i.to_string();
            table.set(&key, "Gaussian", Value::Int(i as i64));
            for (j, annotation) in annotations.iter().enumerate() {
                table.set(&key, annotation, Value::Double(gaussian.mu[j]));
            }
        }
        table
    };
    let covariances = |name: &str, description: &str, model: &Mixture| {
        let mut table = Table::new(name, description, Sorting::DoNotSort);
        table.add_column("Gaussian", "");
        table.add_column("Annotation", "");
        for annotation in annotations {
            table.add_column(annotation, FORMAT);
        }
        let n = annotations.len();
        for (i, gaussian) in model.gaussians.iter().enumerate() {
            for j in 0..n {
                let key = (j + i * n).to_string();
                table.set(&key, "Gaussian", Value::Int(i as i64));
                table.set(&key, "Annotation", Value::Str(annotations[j].clone()));
                for (k, annotation) in annotations.iter().enumerate() {
                    table.set(&key, annotation, Value::Double(gaussian.sigma.get(j, k)));
                }
            }
        }
        table
    };
    report.add_table(means(
        "PositiveModelMeans",
        "Vector of annotation values to describe the (normalized) mean for each Gaussian in the \
         positive model",
        good,
    ));
    report.add_table(covariances(
        "PositiveModelCovariances",
        "Matrix to describe the (normalized) covariance for each Gaussian in the positive model; \
         covariance matrices are joined by row",
        good,
    ));
    if let Some(bad) = bad {
        report.add_table(means(
            "NegativeModelMeans",
            "Vector of annotation values to describe the (normalized) mean for each Gaussian in \
             the negative model",
            bad,
        ));
        report.add_table(covariances(
            "NegativeModelCovariances",
            "Matrix to describe the (normalized) covariance for each Gaussian in the negative \
             model; covariance matrices are joined by row",
            bad,
        ));
    }
    // `GATKReport` keeps its tables in a `TreeMap`, so they are written in name order.
    report.tables.sort_by(|a, b| a.name.cmp(&b.name));
    report.write()
}

/// One table of a model report as read back: its column names and its rows of cells.
#[derive(Debug, Clone, PartialEq)]
pub struct ReportTable {
    pub name: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl ReportTable {
    pub fn column(&self, name: &str) -> Option<usize> {
        self.columns.iter().position(|column| column == name)
    }
}

/// The tables of a `GATKReport` v1.1 whose cells hold no blank, which a model report's never do.
pub fn read_report_tables(text: &str) -> Vec<ReportTable> {
    let mut tables = Vec::new();
    let mut lines = text.lines().peekable();
    while let Some(line) = lines.next() {
        let Some(rest) = line.strip_prefix("#:GATKTable:") else {
            continue;
        };
        // The first `#:GATKTable` line carries the formats, the second the name.
        if rest
            .split(':')
            .next()
            .is_some_and(|field| field.parse::<usize>().is_ok())
        {
            continue;
        }
        let name = rest.split(':').next().unwrap_or_default().to_string();
        let columns: Vec<String> = lines
            .next()
            .unwrap_or_default()
            .split_whitespace()
            .map(str::to_string)
            .collect();
        let mut rows = Vec::new();
        while let Some(row) = lines.peek() {
            if row.trim().is_empty() || row.starts_with("#:") {
                break;
            }
            rows.push(row.split_whitespace().map(str::to_string).collect());
            lines.next();
        }
        tables.push(ReportTable {
            name,
            columns,
            rows,
        });
    }
    tables
}

/// `GMMFromTables`: a mixture rebuilt from a report's means, weights and covariances.
pub fn mixture_from_tables(
    means: &ReportTable,
    covariances: &ReportTable,
    weights: &ReportTable,
    num_annotations: usize,
    num_variants: usize,
    arguments: &ModelArguments,
) -> Mixture {
    let number = |text: &str| java_parse_double(text).unwrap_or(f64::NAN);
    let mut gaussians: Vec<Gaussian> = Vec::new();
    let mut annotation = 0;
    for (c, column) in means.columns.iter().enumerate() {
        if column == "Gaussian" {
            continue;
        }
        for (row, cells) in means.rows.iter().enumerate() {
            if gaussians.len() <= row {
                gaussians.push(Gaussian::new(num_variants, num_annotations));
            }
            if annotation < num_annotations {
                gaussians[row].mu[annotation] = number(&cells[c]);
            }
        }
        annotation += 1;
    }
    if let Some(c) = weights.column("pMixLog10") {
        for (row, cells) in weights.rows.iter().enumerate() {
            if let Some(gaussian) = gaussians.get_mut(row) {
                gaussian.p_mixture_log10 = number(&cells[c]);
            }
        }
    }
    let mut j = 0;
    for (c, column) in covariances.columns.iter().enumerate() {
        if column == "Gaussian" || column == "Annotation" {
            continue;
        }
        for (row, cells) in covariances.rows.iter().enumerate() {
            let g = row / num_annotations.max(1);
            let i = row % num_annotations.max(1);
            if let Some(gaussian) = gaussians.get_mut(g) {
                if j < num_annotations {
                    gaussian.sigma.set(i, j, number(&cells[c]));
                }
            }
        }
        j += 1;
    }
    Mixture::with_gaussians(gaussians, num_annotations, arguments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_diagonal_inverse_is_a_division() {
        let inverse = Matrix::identity(3).times_scalar(200.0).inverse().unwrap();
        assert_eq!(inverse.get(1, 1), 1.0 / 200.0);
        assert_eq!(inverse.get(0, 1), 0.0);
    }

    #[test]
    fn a_pivoted_determinant_keeps_its_sign() {
        let m = Matrix {
            a: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        };
        assert_eq!(m.det(), -1.0);
    }

    #[test]
    fn a_stream_sum_is_compensated() {
        let values = [1.0, 1e-16, 1e-16, 1e-16, 1e-16];
        assert!(stream_sum(values) > 1.0);
        assert_eq!(values.iter().sum::<f64>(), 1.0);
    }

    #[test]
    fn java_min_and_max_let_nan_win() {
        assert!(java_min(f64::NAN, 1.0).is_nan());
        assert!(java_max(1.0, f64::NAN).is_nan());
        assert_eq!(java_min(3.0, -4.0), -4.0);
    }
}
