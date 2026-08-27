//! Apache Commons Math 3's `BrentOptimizer`, as `OptimizationUtils.max` calls it.
//!
//! Golden-section search with parabolic interpolation over one variable. GATK fits the beta shape
//! of a somatic panel of normals through it, so its answers are part of that tool's output and the
//! two tolerances it stops on are part of them too.
//!
//! A function that is exactly a parabola is solved in one step by the interpolation, so the
//! tolerances never reach it. On anything else they decide where the search stops, and a loose
//! enough tolerance stops well short of the maximum.

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
fn number(value: f64) -> String {
    if value == value.trunc() && value.abs() < 1e7 {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

impl OptimizerError {
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

/// Where the search stopped, and what the function was there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointValue {
    pub point: f64,
    pub value: f64,
}

/// `Precision.equals(a, b)` at the default tolerance of one ulp, which for two finite doubles that
/// are not adjacent is plain equality.
fn precision_equals(a: f64, b: f64) -> bool {
    a == b
}

/// `BrentOptimizer.best`, which keeps the earlier of two equal points.
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

/// `OptimizationUtils.max`: `BrentOptimizer` with the given tolerances, maximising.
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

    let (mut a, mut b) = if min < max { (min, max) } else { (max, min) };
    let mut x = start;
    let mut v = x;
    let mut w = x;
    let mut d: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut fx = evaluate(x)?;
    let mut fv = fx;
    let mut fw = fx;

    let unsign = |value: f64| if minimise { value } else { -value };
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

        let fu = evaluate(u)?;

        previous = current;
        current = Some(PointValue {
            point: u,
            value: unsign(fu),
        });
        best = better(best, better(previous, current, minimise), minimise);

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
