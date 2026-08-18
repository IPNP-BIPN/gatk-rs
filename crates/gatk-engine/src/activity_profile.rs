//! Ported from `org.broadinstitute.hellbender.utils.activityprofile.ActivityProfile` and
//! `BandPassActivityProfile` (GATK 4.6.2.0).
//!
//! What decides where an assembly region starts and stops, which is what `HaplotypeCaller` and
//! `Mutect2` assemble over. A boundary that moves by one base changes the haplotypes and therefore
//! the calls, so this is upstream of every variant those tools emit.
//!
//! Four behaviours decide the boundaries and none follows from "cut where the activity changes".
//!
//! # A probability is spread over a Gaussian, and the spread is added, not assigned
//!
//! `BandPassActivityProfile.processState` turns one state into `2 * filterSize + 1` states, each
//! carrying `prob * kernel[i]`, and `incorporateSingleState` **adds** each into whatever is already
//! at that position. So the probability at a site is the sum of the tails of its neighbours, and a
//! site that was never reported active can end up above the threshold.
//!
//! A state whose probability is exactly `0.0` skips the filter entirely and is added as-is, which
//! is the one case where the spread does not happen.
//!
//! # The filter size is chosen from the kernel, not from `sigma`
//!
//! `determineFilterSize` walks in from the edge of the full kernel while the values are at least
//! `MIN_PROB_TO_KEEP_IN_FILTER`, so the width depends on the kernel's *values*. Two sigmas that
//! look similar can give different widths, and the width feeds back into
//! `getMaxProbPropagationDistance`, which decides when a region is allowed to be popped at all.
//!
//! # The cut is a local minimum, searched backwards, with a strict inequality on one side
//!
//! `isMinimum` requires `p[i] <= p[i+1] && p[i] < p[i-1]`: less than or equal on the right, strictly
//! less on the left. On a plateau that asymmetry decides which end is chosen, and the search runs
//! from the far end towards `minRegionSize`, keeping the *last* strict improvement. Neither the
//! direction nor the asymmetry is arbitrary in its effect: they pick a different base.
//!
//! # A region cannot be popped until the profile is long enough
//!
//! `findEndOfRegion` refuses unless `stateList.size() >= maxRegionSize + maxProbPropagationDistance`,
//! so regions appear in bursts rather than as soon as their activity ends. `forceConversion`
//! bypasses it and also trims every state past the profile's span first.

use crate::interval::SimpleInterval;

/// `MIN_PROB_TO_KEEP_IN_FILTER`.
pub const MIN_PROB_TO_KEEP_IN_FILTER: f64 = 1e-5;

/// `BandPassActivityProfile.MAX_FILTER_SIZE`.
pub const MAX_FILTER_SIZE: i32 = 50;

/// `BandPassActivityProfile.DEFAULT_SIGMA`.
pub const DEFAULT_SIGMA: f64 = 17.0;

/// `MathUtils.ROOT_TWO_PI`, which is `Math.sqrt(2.0 * Math.PI)` computed at class initialisation.
///
/// Computed here too rather than written as a literal: a decimal literal is a transcription of a
/// value whose last bits are the point, and `sqrt` is one of the few operations IEEE 754 requires
/// to be correctly rounded, so computing it is exact where copying it is a guess.
fn root_two_pi() -> f64 {
    (2.0 * std::f64::consts::PI).sqrt()
}

/// `ActivityProfileState`, reduced to what the profile reads.
#[derive(Debug, Clone, PartialEq)]
pub struct ProfileState {
    pub contig: String,
    pub start: i32,
    pub is_active_prob: f64,
}

/// A popped region: `AssemblyRegion`'s span and its active flag, which is all this layer decides.
#[derive(Debug, Clone, PartialEq)]
pub struct PoppedRegion {
    pub span: SimpleInterval,
    pub is_active: bool,
    pub extension: i32,
}

/// `ActivityProfile`, with `BandPassActivityProfile`'s filter folded in as an option rather than a
/// subclass: the only thing the subclass changes is `processState` and the propagation distance.
pub struct ActivityProfile {
    max_prob_propagation_distance: i32,
    active_prob_threshold: f64,
    contig_length: i32,
    region_start: Option<i32>,
    region_stop: Option<i32>,
    contig: String,
    states: Vec<ProfileState>,
    /// `None` for the plain profile; `Some(kernel)` for the band pass one.
    kernel: Option<Vec<f64>>,
    filter_size: i32,
}

/// `MathUtils.normalDistribution(mean, sd, x)`.
///
/// Written as the reference writes it, including the order of operations: `exp(...)` divided by
/// `sd * ROOT_TWO_PI` rather than multiplied by a precomputed reciprocal. The division is a
/// separate rounding, so folding it in would change the last bits.
///
/// The `exp` is the **host** libm's, deliberately, and the choice is measured rather than assumed.
/// The reference is `Math.exp` (`MathUtils.normalDistribution`), not `StrictMath.exp`, so
/// `jmath::strict_math::exp` -- which is exact against `StrictMath` -- is faithful to a different
/// function here. Built with it, ten of the 266 kernel values the `activityprofile` suite pins move
/// by an ulp; built with the host `exp`, all 266 match the oracle bit for bit. Do not swap it
/// without redoing that measurement: `docs/numeric-functions-a-ported-call-site-reaches.md`.
pub fn normal_distribution(mean: f64, sd: f64, x: f64) -> f64 {
    assert!(sd >= 0.0, "sd: Standard deviation of normal must be >= 0");
    (-(x - mean) * (x - mean) / (2.0 * sd * sd)).exp() / (sd * root_two_pi())
}

/// `BandPassActivityProfile.makeKernel`.
pub fn make_kernel(filter_size: i32, sigma: f64) -> Vec<f64> {
    let band_size = 2 * filter_size + 1;
    let mut kernel: Vec<f64> = (0..band_size)
        .map(|i| normal_distribution(filter_size as f64, sigma, i as f64))
        .collect();
    // `normalizeSumToOne`: one sum, then one division per element, in index order.
    let sum: f64 = kernel.iter().sum();
    for value in &mut kernel {
        *value /= sum;
    }
    kernel
}

/// `BandPassActivityProfile.determineFilterSize`.
pub fn determine_filter_size(kernel: &[f64], min_prob_to_keep: f64) -> i32 {
    let middle = (kernel.len() - 1) / 2;
    let mut filter_end = middle;
    while filter_end > 0 {
        if kernel[filter_end - 1] < min_prob_to_keep {
            break;
        }
        filter_end -= 1;
    }
    (middle - filter_end) as i32
}

impl ActivityProfile {
    /// The plain profile, whose `processState` passes a state straight through.
    pub fn new(
        max_prob_propagation_distance: i32,
        active_prob_threshold: f64,
        contig: &str,
        contig_length: i32,
    ) -> Self {
        ActivityProfile {
            max_prob_propagation_distance,
            active_prob_threshold,
            contig_length,
            region_start: None,
            region_stop: None,
            contig: contig.to_string(),
            states: Vec::new(),
            kernel: None,
            filter_size: 0,
        }
    }

    /// `BandPassActivityProfile`, with `adaptiveFilterSize` as the reference's four-argument
    /// constructor sets it: true.
    pub fn band_pass(
        max_prob_propagation_distance: i32,
        active_prob_threshold: f64,
        max_filter_size: i32,
        sigma: f64,
        adaptive: bool,
        contig: &str,
        contig_length: i32,
    ) -> Self {
        assert!(sigma >= 0.0, "Sigma must be greater than or equal to 0");
        // The full kernel is built at the maximum width first, and its *values* choose the width
        // that is then rebuilt. Building once at the final width would give a different kernel,
        // because the normalisation divides by a different sum.
        let full = make_kernel(max_filter_size, sigma);
        let filter_size = if adaptive {
            determine_filter_size(&full, MIN_PROB_TO_KEEP_IN_FILTER)
        } else {
            max_filter_size
        };
        let mut profile = ActivityProfile::new(
            max_prob_propagation_distance,
            active_prob_threshold,
            contig,
            contig_length,
        );
        profile.kernel = Some(make_kernel(filter_size, sigma));
        profile.filter_size = filter_size;
        profile
    }

    /// `getMaxProbPropagationDistance`, which the band pass one extends by its filter size.
    pub fn max_prob_propagation_distance(&self) -> i32 {
        self.max_prob_propagation_distance + self.filter_size
    }

    pub fn kernel(&self) -> Option<&[f64]> {
        self.kernel.as_deref()
    }

    pub fn filter_size(&self) -> i32 {
        self.filter_size
    }

    pub fn size(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    /// `getEnd()`: the end of the region stop location, which is the last position **added** and
    /// not the last state, because the band pass filter writes past it.
    ///
    /// The traversal compares this against the next locus to decide whether to force a conversion,
    /// so taking the last state instead would stop it noticing a gap in the loci.
    pub fn end(&self) -> i32 {
        self.region_stop.unwrap_or(0)
    }

    /// `getSpan().size()`: the number of positions actually added, which is not the number of
    /// states, because the filter writes past the last added position.
    fn span_size(&self) -> usize {
        match (self.region_start, self.region_stop) {
            (Some(start), Some(stop)) => (stop - start + 1) as usize,
            _ => 0,
        }
    }

    /// `getLocForOffset`, which returns nothing off either end of the contig.
    fn loc_for_offset(&self, relative_start: i32, offset: i32) -> Option<i32> {
        let start = relative_start + offset;
        if start < 1 || start > self.contig_length {
            None
        } else {
            Some(start)
        }
    }

    /// `add`, which refuses a position that is not exactly one past the last one.
    pub fn add(&mut self, start: i32, is_active_prob: f64) {
        match self.region_stop {
            None => {
                self.region_start = Some(start);
                self.region_stop = Some(start);
            }
            Some(stop) => {
                assert_eq!(
                    stop,
                    start - 1,
                    "Bad add call to ActivityProfile: loc not immediately after last loc"
                );
                self.region_stop = Some(start);
            }
        }

        for state in self.process_state(start, is_active_prob) {
            self.incorporate(state);
        }
    }

    /// `BandPassActivityProfile.processState` over `ActivityProfile.processState`.
    ///
    /// The soft-clip branch of the parent is not here: it needs `ActivityProfileState.Type`, which
    /// only `HaplotypeCaller`'s activity calculation produces, and it is its own slice.
    fn process_state(&self, start: i32, is_active_prob: f64) -> Vec<(i32, f64)> {
        let Some(kernel) = &self.kernel else {
            return vec![(start, is_active_prob)];
        };
        // A probability of exactly zero is added unfiltered, so it lands on one position instead
        // of being spread over the band.
        if is_active_prob <= 0.0 {
            return vec![(start, is_active_prob)];
        }
        let mut states = Vec::new();
        for offset in -self.filter_size..=self.filter_size {
            if let Some(position) = self.loc_for_offset(start, offset) {
                states.push((
                    position,
                    is_active_prob * kernel[(offset + self.filter_size) as usize],
                ));
            }
        }
        states
    }

    /// `incorporateSingleState`: **adds** into an existing position, appends one past the end, and
    /// silently ignores anything before the profile's start.
    fn incorporate(&mut self, (start, prob): (i32, f64)) {
        let region_start = self.region_start.expect("add sets the start first");
        let position = start - region_start;
        assert!(
            position <= self.states.len() as i32,
            "Must add state contiguous to existing states"
        );
        if position < 0 {
            return;
        }
        let position = position as usize;
        if position < self.states.len() {
            self.states[position].is_active_prob += prob;
        } else {
            self.states.push(ProfileState {
                contig: self.contig.clone(),
                start,
                is_active_prob: prob,
            });
        }
    }

    /// The probabilities as the profile currently holds them, which is what a divergence in the
    /// kernel shows up in first.
    pub fn probabilities(&self) -> Vec<f64> {
        self.states.iter().map(|s| s.is_active_prob).collect()
    }

    /// `popReadyAssemblyRegions`.
    pub fn pop_ready_regions(
        &mut self,
        extension: i32,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Vec<PoppedRegion> {
        assert!(extension >= 0, "assemblyRegionExtension must be >= 0");
        assert!(min_region_size > 0, "minRegionSize must be >= 1");
        assert!(max_region_size > 0, "maxRegionSize must be >= 1");

        let mut regions = Vec::new();
        while let Some(region) = self.pop_next_ready_region(
            extension,
            min_region_size,
            max_region_size,
            force_conversion,
        ) {
            regions.push(region);
        }
        regions
    }

    /// `popNextReadyAssemblyRegion`.
    fn pop_next_ready_region(
        &mut self,
        extension: i32,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Option<PoppedRegion> {
        if self.states.is_empty() {
            return None;
        }

        // Flushing trims every state the filter wrote past the last position that was added, so a
        // forced pop cannot produce a region outside the interval being processed.
        if force_conversion {
            let span = self.span_size();
            if span < self.states.len() {
                self.states.truncate(span);
            }
        }

        let first_start = self.states[0].start;
        let is_active = self.states[0].is_active_prob > self.active_prob_threshold;
        let end_offset = self.find_end_of_region(
            is_active,
            min_region_size,
            max_region_size,
            force_conversion,
        )?;

        self.states.drain(0..=end_offset);
        if self.states.is_empty() {
            self.region_start = None;
            self.region_stop = None;
        } else {
            self.region_start = Some(self.states[0].start);
        }

        Some(PoppedRegion {
            span: SimpleInterval::new(&self.contig, first_start, first_start + end_offset as i32)
                .expect("a region built from two positions on the same contig"),
            is_active,
            extension,
        })
    }

    /// `findEndOfRegion`, returning the index of the region's last state.
    fn find_end_of_region(
        &self,
        is_active: bool,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Option<usize> {
        if !force_conversion
            && self.states.len() < max_region_size + self.max_prob_propagation_distance() as usize
        {
            // Not enough profile yet for the probabilities near the end to be final.
            return None;
        }

        let mut end = self.find_first_activity_boundary(is_active, max_region_size);
        if is_active && end == max_region_size {
            end = self.find_best_cut_site(end, min_region_size);
        }
        if end == 0 {
            // `endOfActiveRegion - 1` is -1, which the reference returns as "not found".
            return None;
        }
        Some(end - 1)
    }

    /// `findFirstActivityBoundary`.
    fn find_first_activity_boundary(&self, is_active: bool, max_region_size: usize) -> usize {
        let mut end = 0usize;
        while end < self.states.len() && end < max_region_size {
            if (self.states[end].is_active_prob > self.active_prob_threshold) != is_active {
                break;
            }
            end += 1;
        }
        end
    }

    /// `findBestCutSite`: the global minimum in `[minRegionSize - 1, end - 1]`, searched from the
    /// far end, keeping the last index that strictly improves on the best so far.
    fn find_best_cut_site(&self, end_of_active_region: usize, min_region_size: usize) -> usize {
        assert!(
            end_of_active_region >= min_region_size,
            "endOfActiveRegion must be >= minRegionSize"
        );
        let mut min_index = end_of_active_region - 1;
        let mut min_prob = f64::MAX;

        let mut index = min_index as i64;
        while index >= min_region_size as i64 - 1 {
            let current = self.states[index as usize].is_active_prob;
            if current < min_prob && self.is_minimum(index as usize) {
                min_prob = current;
                min_index = index as usize;
            }
            index -= 1;
        }
        min_index + 1
    }

    /// `isMinimum`, whose two comparisons are deliberately different: `<=` on the right and `<` on
    /// the left, so a plateau is a minimum only at its left-hand end.
    fn is_minimum(&self, index: usize) -> bool {
        if index == self.states.len() - 1 {
            return false;
        }
        if index < 1 {
            return false;
        }
        let here = self.states[index].is_active_prob;
        here <= self.states[index + 1].is_active_prob
            && here < self.states[index - 1].is_active_prob
    }
}
