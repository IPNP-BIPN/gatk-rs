//! Apache Commons Math 3's `BrentOptimizer`, as `OptimizationUtils.max` calls it.
//!
//! Golden-section search with parabolic interpolation over one variable. GATK fits the beta shape
//! of a somatic panel of normals through it, so its answers are part of that tool's output and the
//! two tolerances it stops on are part of them too.
//!
//! A function that is exactly a parabola is solved in one step by the interpolation, so the
//! tolerances never reach it. On anything else they decide where the search stops, and a loose
//! enough tolerance stops well short of the maximum.
//!
//! # The tolerances are the caller's
//!
//! `OptimizationUtils.max` has an overload that takes them, and that is the one the panel of
//! normals reaches: `OptimizationUtils.max(logLikelihood, 0.01, 100, 1, 0.01, 0.1, 100)` at
//! `CreateSomaticPanelOfNormals.java:251`, a relative tolerance of `0.01`, an absolute one of `0.1`
//! and a budget of 100 evaluations. The overload builds a fresh optimiser per call
//! (`OptimizationUtils.java:38`):
//!
//! ```java
//! final BrentOptimizer optimizer = new BrentOptimizer(relativeTolerance, absoluteTolerance);
//! ```
//!
//! The `brent-optimizer` golden (`crates/gatk-engine/tests/data/brent_optimizer.txt.gz`) holds
//! fifteen such answers and four refusals. Its `kinked-*` rows are the ones the tolerances move:
//! the same peak at 7 comes back as `7.00000000000277` at `1e-12` and as `7.28924288066983` at
//! `0.5`.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<PointValue, OptimizerError>` is "the answer or one of the three refusals below"; it
//!   is how this codebase models a method that throws. A `?` after a call that returns a `Result`
//!   rethrows: on an error the function returns it at once, on success the value is unwrapped.
//! - `Option<PointValue>` is "a point or nothing", Java's nullable reference with the null check
//!   enforced by the compiler; `Some(x)` is a value and `None` is the null.
//! - `F: Fn(f64) -> f64` is the function being optimised, Java's `DoubleUnaryOperator`; `|x| ...`
//!   is a lambda. A lambda that changes a captured variable (the evaluation counter below) is
//!   declared `let mut`.
//! - `f64` is Java's `double` and `usize` an unsigned integer used for counts and indices.

/// `BrentOptimizer.GOLDEN_SECTION`, which is `0.5 * (3 - sqrt(5))` evaluated in double precision.
///
/// The literal is the shortest decimal that round-trips to that product's exact bits. One ulp out
/// is enough to send the search somewhere else: written a digit longer, the symmetric case
/// converged on the other side of zero.
pub const GOLDEN_SECTION: f64 = 0.3819660112501051;

/// What the optimiser refuses, each with Apache Commons' own wording.
#[derive(Debug, Clone, PartialEq)]
pub enum OptimizerError {
    /// The evaluation budget, which is a refusal rather than a silent stop.
    TooManyEvaluations { max_evaluations: usize },
    /// A guess outside the interval.
    OutOfRange { value: f64, min: f64, max: f64 },
    /// An interval whose ends are the wrong way round.
    IntervalInverted { min: f64, max: f64 },
}

/// Apache Commons prints a whole number without a decimal point and anything else with one.
///
/// What the golden's `error` rows show: `-1 out of [0, 10] range` for a guess of `-1.0`, and
/// `10 is larger than, or equal to, the maximum (0)` for an inverted interval.
fn number(value: f64) -> String {
    if value == value.trunc() && value.abs() < 1e7 {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

impl OptimizerError {
    /// The exception's message, as the golden's `error` rows print it after the class name.
    pub fn message(&self) -> String {
        match self {
            OptimizerError::TooManyEvaluations { max_evaluations } => {
                format!("illegal state: maximal count ({max_evaluations}) exceeded: evaluations")
            }
            OptimizerError::OutOfRange { value, min, max } => format!(
                "{} out of [{}, {}] range",
                number(*value),
                number(*min),
                number(*max)
            ),
            OptimizerError::IntervalInverted { min, max } => format!(
                "{} is larger than, or equal to, the maximum ({})",
                number(*min),
                number(*max)
            ),
        }
    }
}

/// `UnivariatePointValuePair`: where the search stopped, and what the function was there.
///
/// `value` is the function's own value, not the negated one a maximisation searches on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointValue {
    pub point: f64,
    pub value: f64,
}

/// `Precision.equals(a, b)` at the default tolerance of one ulp, which for two finite doubles that
/// are not adjacent is plain equality.
///
/// The port is plain equality throughout, so it would part from the reference only on two
/// adjacent doubles among `x`, `v` and `w`, which the one-ulp tolerance calls equal.
fn precision_equals(a: f64, b: f64) -> bool {
    a == b
}

/// `BrentOptimizer.best`, which keeps the earlier of two equal points.
///
/// Either side may be absent (`previous` is, before the first step), and then the other is the
/// answer. `match` over the pair of `Option`s takes the first arm whose shape
/// fits, and the compiler checks that every combination is covered.
fn better(a: Option<PointValue>, b: Option<PointValue>, minimise: bool) -> Option<PointValue> {
    match (a, b) {
        (None, other) | (other, None) => other,
        (Some(left), Some(right)) => {
            let take_left = if minimise {
                left.value <= right.value
            } else {
                left.value >= right.value
            };
            Some(if take_left { left } else { right })
        }
    }
}

/// `SearchInterval`'s own validation, in its own order: the interval first, then the guess.
fn check_interval(min: f64, max: f64, start: f64) -> Result<(), OptimizerError> {
    if min >= max {
        return Err(OptimizerError::IntervalInverted { min, max });
    }
    if start < min || start > max {
        return Err(OptimizerError::OutOfRange {
            value: start,
            min,
            max,
        });
    }
    Ok(())
}

/// `OptimizationUtils.max(function, min, max, guess, relativeTolerance, absoluteTolerance,
/// maxEvaluations)`: `BrentOptimizer` with the given tolerances, `GoalType.MAXIMIZE`.
///
/// Maximising is minimising the negative, which the optimiser does internally, so the value that
/// comes back is the function's own.
pub fn maximize<F: Fn(f64) -> f64>(
    function: F,
    min: f64,
    max: f64,
    start: f64,
    relative_threshold: f64,
    absolute_threshold: f64,
    max_evaluations: usize,
) -> Result<PointValue, OptimizerError> {
    optimize(
        function,
        min,
        max,
        start,
        relative_threshold,
        absolute_threshold,
        max_evaluations,
        false,
    )
}

/// The same search minimising, which is what `BrentOptimizer` does natively.
pub fn minimize<F: Fn(f64) -> f64>(
    function: F,
    min: f64,
    max: f64,
    start: f64,
    relative_threshold: f64,
    absolute_threshold: f64,
    max_evaluations: usize,
) -> Result<PointValue, OptimizerError> {
    optimize(
        function,
        min,
        max,
        start,
        relative_threshold,
        absolute_threshold,
        max_evaluations,
        true,
    )
}

/// `BrentOptimizer.doOptimize`.
///
/// The two tolerances make one bound, `relative * |x| + absolute`, so the relative one is scaled by
/// where the search currently is. The stopping test compares the distance from the midpoint against
/// twice that bound less half the interval, which is Brent's own criterion.
///
/// The names are Brent's: `x` is the best point so far, `w` the second best, `v` the previous `w`;
/// `[a, b]` is the bracket, `d` the step just taken and `e` the step before it. The statements
/// keep the reference's order and temporaries, so the points visited are the reference's points.
#[allow(clippy::too_many_arguments)]
fn optimize<F: Fn(f64) -> f64>(
    function: F,
    min: f64,
    max: f64,
    start: f64,
    relative_threshold: f64,
    absolute_threshold: f64,
    max_evaluations: usize,
    minimise: bool,
) -> Result<PointValue, OptimizerError> {
    check_interval(min, max, start)?;

    let mut evaluations = 0usize;
    // `computeObjectiveValue` counts the call and refuses past the budget, so the budget belongs to
    // the search rather than to the function.
    let mut evaluate = |x: f64| -> Result<f64, OptimizerError> {
        evaluations += 1;
        if evaluations > max_evaluations {
            return Err(OptimizerError::TooManyEvaluations { max_evaluations });
        }
        let value = function(x);
        Ok(if minimise { value } else { -value })
    };

    // `check_interval` has already refused `min >= max`, so this always keeps `(min, max)`.
    // `let (a, b) = ...` binds both names from one pair.
    let (mut a, mut b) = if min < max { (min, max) } else { (max, min) };
    let mut x = start;
    let mut v = x;
    let mut w = x;
    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut fx = evaluate(x)?;
    let mut fv = fx;
    let mut fw = fx;

    // `evaluate` hands the search the negated value when maximising; `unsign` turns it back
    // before a point is recorded, so a returned `PointValue` carries the function's own value.
    let unsign = |value: f64| if minimise { value } else { -value };
    // The last two evaluated points and the best seen so far; the search answers with the best of
    // the three when it stops.
    let mut previous: Option<PointValue> = None;
    let mut current = Some(PointValue {
        point: x,
        value: unsign(fx),
    });
    let mut best = current;

    loop {
        let m = 0.5 * (a + b);
        let tol1 = relative_threshold * x.abs() + absolute_threshold;
        let tol2 = 2.0 * tol1;

        if (x - m).abs() <= tol2 - 0.5 * (b - a) {
            // Brent's own termination.
            return Ok(better(best, better(previous, current, minimise), minimise)
                .expect("at least one evaluated point"));
        }

        if e.abs() > tol1 {
            // Fit a parabola through the three best points so far.
            let mut r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
            if p > q * (a - x) && p < q * (b - x) && p.abs() < (0.5 * q * r).abs() {
                d = p / q;
                // The interpolated point is checked against the ends and the step replaced, not
                // the point: `u` is recomputed from `d` below either way.
                let interpolated = x + d;
                if interpolated - a < tol2 || b - interpolated < tol2 {
                    d = if x <= m { tol1 } else { -tol1 };
                }
            } else {
                e = if x < m { b - x } else { a - x };
                d = GOLDEN_SECTION * e;
            }
        } else {
            e = if x < m { b - x } else { a - x };
            d = GOLDEN_SECTION * e;
        }

        // The step is never smaller than the tolerance, so the search always moves.
        let u = if d.abs() < tol1 {
            if d >= 0.0 {
                x + tol1
            } else {
                x - tol1
            }
        } else {
            x + d
        };

        // `?`: a refusal from the budget ends the search with that error.
        let fu = evaluate(u)?;

        previous = current;
        current = Some(PointValue {
            point: u,
            value: unsign(fu),
        });
        best = better(best, better(previous, current, minimise), minimise);

        // Update the bracket and the three remembered points: `u` becomes the best if it is no
        // worse than `x`, otherwise it narrows the bracket and may replace `w` or `v`.
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || precision_equals(w, x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || precision_equals(v, x) || precision_equals(v, w) {
                v = u;
                fv = fu;
            }
        }
    }
}
