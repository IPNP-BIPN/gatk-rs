//! `org.hipparchus.stat.fitting.EmpiricalDistribution` with its default thousand bins, as
//! `VariantAnnotationsScorer.createScoreToCalibrationSensitivityConverter` uses it: `load` and
//! `cumulativeProbability`. Ported from Hipparchus 2.0 (Apache 2.0).
//!
//! A bin holding one value, or values all equal, is a point; any other is a Gaussian of the bin's
//! own mean and deviation, cut off at the bin's edges and scaled to the bin's mass.
use jmath::gamma::{erf, erfc};

const BINS: usize = 1000;

/// `Erf.X_CRIT`, the median of the normal over `sqrt(2)`.
const X_CRIT: f64 = 0.476_936_276_204_469_7;

/// `StreamingStatistics`: the count, and the mean and second moment `SecondMoment` keeps.
#[derive(Debug, Clone, Default)]
struct Statistics {
    n: u64,
    mean: f64,
    moment: f64,
}

impl Statistics {
    /// `SecondMoment.increment`.
    fn add(&mut self, value: f64) {
        if self.n < 1 {
            self.mean = 0.0;
            self.moment = 0.0;
        }
        self.n += 1;
        let n = self.n as f64;
        let deviation = value - self.mean;
        let scaled = deviation / n;
        self.mean += scaled;
        self.moment += (n - 1.0) * deviation * scaled;
    }

    /// `Variance.getResult`, bias-corrected.
    fn variance(&self) -> f64 {
        match self.n {
            0 => f64::NAN,
            1 => 0.0,
            n => self.moment / (n as f64 - 1.0),
        }
    }
}

/// `getKernel`: a point at the mean, or a Gaussian.
enum Kernel {
    Constant(f64),
    Normal { mean: f64, deviation: f64 },
}

impl Kernel {
    fn of(statistics: &Statistics) -> Kernel {
        if statistics.n < 2 || statistics.variance() == 0.0 {
            Kernel::Constant(statistics.mean)
        } else {
            Kernel::Normal {
                mean: statistics.mean,
                deviation: statistics.variance().sqrt(),
            }
        }
    }

    /// `NormalDistribution.cumulativeProbability`.
    fn normal_cumulative(mean: f64, deviation: f64, x: f64) -> f64 {
        let dev = x - mean;
        if dev.abs() > 40.0 * deviation {
            return if dev < 0.0 { 0.0 } else { 1.0 };
        }
        0.5 * erfc(-dev / (deviation * std::f64::consts::SQRT_2))
    }

    /// `Erf.erf(x1, x2)`, which takes the difference of whichever two tails lose the least.
    fn erf_between(x1: f64, x2: f64) -> f64 {
        if x1 > x2 {
            return -Self::erf_between(x2, x1);
        }
        if x1 < -X_CRIT {
            if x2 < 0.0 {
                erfc(-x2) - erfc(-x1)
            } else {
                erf(x2) - erf(x1)
            }
        } else if x2 > X_CRIT && x1 > 0.0 {
            erfc(x1) - erfc(x2)
        } else {
            erf(x2) - erf(x1)
        }
    }

    /// `NormalDistribution.probability(x0, x1)`.
    fn normal_probability(mean: f64, deviation: f64, x0: f64, x1: f64) -> f64 {
        let denominator = deviation * std::f64::consts::SQRT_2;
        let v0 = (x0 - mean) / denominator;
        let v1 = (x1 - mean) / denominator;
        0.5 * Self::erf_between(v0, v1)
    }
}

/// A loaded distribution.
#[derive(Debug, Clone)]
pub struct EmpiricalDistribution {
    min: f64,
    max: f64,
    delta: f64,
    bins: Vec<Statistics>,
    upper_bounds: Vec<f64>,
}

impl EmpiricalDistribution {
    /// `load(double[])`, whose values the caller has already checked are finite.
    pub fn load(values: &[f64]) -> EmpiricalDistribution {
        let min = values.iter().copied().fold(f64::NAN, f64::min);
        let max = values.iter().copied().fold(f64::NAN, f64::max);
        let mut out = EmpiricalDistribution {
            min,
            max,
            delta: (max - min) / BINS as f64,
            bins: vec![Statistics::default(); BINS],
            upper_bounds: vec![0.0; BINS],
        };
        for value in values {
            let bin = out.find_bin(*value);
            out.bins[bin].add(*value);
        }
        let n = values.len() as f64;
        out.upper_bounds[0] = out.bins[0].n as f64 / n;
        for i in 1..BINS - 1 {
            out.upper_bounds[i] = out.upper_bounds[i - 1] + out.bins[i].n as f64 / n;
        }
        out.upper_bounds[BINS - 1] = 1.0;
        out
    }

    /// `findBin`: `(int) ceil(...)`, which saturates, and a `NaN` of a zero width becomes 0.
    fn find_bin(&self, value: f64) -> usize {
        let raw = ((value - self.min) / self.delta).ceil();
        let index = (raw as i32).wrapping_sub(1);
        index.max(0).min(BINS as i32 - 1) as usize
    }

    /// `getUpperBounds()[i]`.
    fn bin_upper(&self, i: usize) -> f64 {
        if i == BINS - 1 {
            self.max
        } else {
            self.min + self.delta * (i as f64 + 1.0)
        }
    }

    /// `cumulativeProbability(x)`.
    pub fn cumulative_probability(&self, x: f64) -> f64 {
        if x < self.min {
            return 0.0;
        }
        if x >= self.max {
            return 1.0;
        }
        let bin = self.find_bin(x);
        let below = if bin == 0 {
            0.0
        } else {
            self.upper_bounds[bin - 1]
        };
        let mass = if bin == 0 {
            self.upper_bounds[0]
        } else {
            self.upper_bounds[bin] - self.upper_bounds[bin - 1]
        };
        match Kernel::of(&self.bins[bin]) {
            Kernel::Constant(mean) => {
                if x < mean {
                    below
                } else {
                    below + mass
                }
            }
            Kernel::Normal { mean, deviation } => {
                let lower = if bin == 0 {
                    self.min
                } else {
                    self.bin_upper(bin - 1)
                };
                let k_b = Kernel::normal_probability(mean, deviation, lower, self.bin_upper(bin));
                let within = (Kernel::normal_cumulative(mean, deviation, x)
                    - Kernel::normal_cumulative(mean, deviation, lower))
                    / k_b;
                below + mass * within
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn outside_the_values_the_answer_is_exact() {
        let distribution = EmpiricalDistribution::load(&[1.0, 2.0, 3.0]);
        assert_eq!(distribution.cumulative_probability(0.5), 0.0);
        assert_eq!(distribution.cumulative_probability(3.0), 1.0);
    }

    #[test]
    fn a_bin_holding_one_value_is_a_step() {
        let distribution = EmpiricalDistribution::load(&[1.0, 2.0, 3.0]);
        let below = distribution.cumulative_probability(1.5);
        assert_eq!(below, 1.0 / 3.0);
        assert_eq!(distribution.cumulative_probability(2.0), 2.0 / 3.0);
    }

    #[test]
    fn the_bins_inner_edges_are_the_left_open_ones() {
        let distribution = EmpiricalDistribution::load(&[0.0, 1000.0]);
        assert_eq!(distribution.find_bin(0.0), 0);
        assert_eq!(distribution.find_bin(1.0), 0);
        assert_eq!(distribution.find_bin(1.5), 1);
        assert_eq!(distribution.find_bin(1000.0), 999);
    }
}
