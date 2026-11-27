//! `PersistenceOptimizer`, ported from
//! `org.broadinstitute.hellbender.tools.copynumber.utils.optimization` (GATK 4.7.0.0;
//! `PersistenceOptimizer.java` is unchanged since 4.6.2.0).
//!
//! A watershed over one-dimensional data: it returns every local minimum ordered by topological
//! persistence, which is how `ContaminationSegmenter` picks the changepoint candidates the kernel
//! segmenter then scores.
//!
//! # The ordering is `Double.compare`, and that is not `<`
//!
//! ```java
//! .sorted(Comparator.comparingDouble(i -> this.data[i]))
//! ```
//!
//! `comparingDouble` is `Double.compare`, which falls back to the **bit pattern** when neither
//! value is less than the other. So `-0.0` sorts below `0.0` and a `NaN` sorts above everything,
//! including positive infinity. The watershed starts from whichever point that ordering calls
//! lowest, so the answer changes with it: on `[0.0, -0.0, 0.0, -0.0]` the global minimum is index
//! 1 and one persistence comes out as **negative zero**, and data holding a `NaN` has that `NaN`
//! as its global maximum, which makes the global persistence `NaN`. The
//! `persistence_optimizer.txt.gz` golden pins both: its `persistence signed-zeroes` row is
//! `0.0,-0.0` and its `persistence with-nan` row is `NaN,NaN,1.5`.
//!
//! # A plateau is a minimum at its left end
//!
//! The sort is stable, so equal values keep their index order and the watershed reaches the
//! leftmost point of a plateau first. That point creates the component; the rest of the plateau
//! extends it. It is the documented behaviour, and it is a consequence of the sort rather than of
//! any test in the algorithm. The Java constructor's own Javadoc says so:
//!
//! ```java
//!  * Note that for a region where the data takes constant values that occurs either to the right
//!  * of a local maximum or at the beginning of the data, the leftmost point in the region
//!  * is considered a local minimum.
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `&[f64]` is a **borrowed, read-only view** of a `double[]`: the function reads it and the
//!   caller keeps it. `&mut [..]` is the writable form, which is how `colors` and `components`
//!   are updated in place, as the reference's `int[]` and `List` are.
//! - `Option<usize>` is "an index or nothing", checked by the compiler: `None` or `Some(c)`.
//! - `Result<Persistence, EmptyData>` is "the answer or the refusal"; `Err(..)` is the reference's
//!   `IllegalArgumentException`.
//! - `match (a, b) { (None, None) => .., (Some(left), None) => .., .. }` tests both neighbours at
//!   once and names the colour it found. The compiler checks that the four arms cover every case,
//!   which is the reference's four-way `if` chain.
//! - `|value: f64| -> i64 { .. }` is a closure, Java's lambda.
//! - `as i64` on a `u64` keeps the 64 bits and reads them as signed, which gives the `long` that
//!   `doubleToLongBits` returns.

use std::cmp::Ordering;

/// `Double.compare`, which orders `-0.0` below `0.0` and every `NaN` above `Infinity`.
///
/// ```java
/// if (d1 < d2) return -1;
/// if (d1 > d2) return 1;
/// long thisBits = Double.doubleToLongBits(d1);
/// ```
///
/// `doubleToLongBits` collapses every `NaN` to one pattern, so two `NaN`s compare equal, and the
/// comparison is on the **signed** long, which is what puts `-0.0` first.
///
/// The bits are reached only for equal values and for a `NaN`, as in the reference: there they
/// separate `-0.0` from `0.0`, and the canonical `NaN` pattern `0x7ff8000000000000` is above the
/// bits of every non-NaN double, positive infinity (`0x7ff0000000000000`) included. `to_bits` is
/// `doubleToRawLongBits`, which keeps a `NaN`'s payload, so the closure substitutes the canonical
/// pattern first.
pub fn java_compare(first: f64, second: f64) -> Ordering {
    if first < second {
        return Ordering::Less;
    }
    if first > second {
        return Ordering::Greater;
    }
    let bits = |value: f64| -> i64 {
        if value.is_nan() {
            0x7ff8_0000_0000_0000u64 as i64
        } else {
            value.to_bits() as i64
        }
    };
    bits(first).cmp(&bits(second))
}

/// `PersistenceOptimizer.Component`: a connected component, bounded by its two ends and named by
/// its lowest point. The reference creates one only at a local minimum, so `min_index` and
/// `min_value` are fixed at birth (they are `final` there) and only the two ends move.
struct Component {
    left_index: usize,
    right_index: usize,
    min_index: usize,
    min_value: f64,
}

/// `PersistenceOptimizer.ExtremaPair`: one minimum paired with the maximum that ends its
/// component. The reference also keeps `maxIndex`, commented "these could be used to find local
/// maxima, if desired", and never reads it; the port drops the field.
struct ExtremaPair {
    min_index: usize,
    persistence: f64,
}

impl ExtremaPair {
    /// `new ExtremaPair(data, datumIndex1, datumIndex2)`: the pair, whose ends are decided by value
    /// and, when the values are equal, by index.
    fn new(data: &[f64], first: usize, second: usize) -> ExtremaPair {
        // Three branches in the reference, of which the second and third have the same answer:
        // the value decides, and only when the two values are equal does the index. A `NaN` makes
        // both comparisons false, so it is the index that decides there too.
        let (min_index, max_index) = if data[first] > data[second] {
            (second, first)
        } else if data[first] < data[second] || first < second {
            (first, second)
        } else {
            (second, first)
        };
        ExtremaPair {
            min_index,
            // Not `abs`: the subtraction is the reference's, and on two zeroes of different signs
            // it gives -0.0.
            persistence: data[max_index] - data[min_index],
        }
    }
}

/// What `getMinimaIndices()` and `getPersistences()` return, as two parallel lists.
#[derive(Debug, Clone, PartialEq)]
pub struct Persistence {
    /// The indices, the global minimum first and the rest by decreasing persistence.
    pub minima_indices: Vec<usize>,
    /// The matching persistences, the whole range first.
    pub persistences: Vec<f64>,
}

/// The argument check, which is the only thing this refuses:
/// `Utils.validateArg(data.length > 0, "Data must contain at least one element.")`. A `null` array
/// (`Utils.nonNull`) cannot be expressed in the Rust signature, so it has no variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EmptyData;

impl EmptyData {
    /// The exception's message, as `validateArg` passes it through.
    pub fn message(&self) -> &'static str {
        "Data must contain at least one element."
    }

    /// The class `Utils.validateArg` throws.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }
}

/// `new PersistenceOptimizer(data)`, whose whole work happens in the constructor.
///
/// **How**: stably sort the indices by value, run the watershed, then prepend the global minimum
/// and the whole range, exactly as the reference's two `add(0, ..)` calls do:
///
/// ```java
/// minimaIndices.add(0, sortedIndices.get(0));     //add index of global minimum
/// ```
///
/// The prepended persistence is "last sorted value minus first sorted value", not a maximum over
/// the data, so the `Double.compare` ordering decides it too (a `NaN` sorts last and wins).
/// The reference's defensive copy of `data` is not needed: the borrow cannot be changed under it.
pub fn persistence_optimizer(data: &[f64]) -> Result<Persistence, EmptyData> {
    if data.is_empty() {
        return Err(EmptyData);
    }

    // The stable sort by `Double.compare`, which decides where the watershed starts. `sort_by` is
    // a stable sort, as Java's `Stream.sorted` is on an ordered stream, so equal values keep their
    // index order. `(0..n).collect()` is `IntStream.range(0, n).boxed().collect(toList())`.
    let mut sorted_indices: Vec<usize> = (0..data.len()).collect();
    sorted_indices.sort_by(|first, second| java_compare(data[*first], data[*second]));

    let pairs = find_extrema_pairs(data, &sorted_indices);

    // `pairs.iter().map(..).collect()` is the reference's `stream().map(..).collect(..)`.
    let mut minima_indices: Vec<usize> = pairs.iter().map(|pair| pair.min_index).collect();
    minima_indices.insert(0, sorted_indices[0]);
    let mut persistences: Vec<f64> = pairs.iter().map(|pair| pair.persistence).collect();
    persistences.insert(
        0,
        data[sorted_indices[sorted_indices.len() - 1]] - data[sorted_indices[0]],
    );

    Ok(Persistence {
        minima_indices,
        persistences,
    })
}

/// `findExtremaPairs`: the watershed, walked from the lowest point upward.
///
/// **How**: each point, in sorted order, either starts a component (no coloured neighbour: a
/// local minimum), joins the one neighbour that is coloured, or, with both neighbours coloured,
/// is a local maximum that closes the shallower component into a pair and merges the two. The
/// first and last points have one neighbour each, so the reference tests them first and skips the
/// two-sided test with `continue`; the port does the same.
fn find_extrema_pairs(data: &[f64], sorted_indices: &[usize]) -> Vec<ExtremaPair> {
    if data.len() == 1 {
        return Vec::new();
    }

    let mut components: Vec<Component> = Vec::with_capacity(data.len());
    // `NO_COLOR` is -1, so an `Option` says the same thing without a sentinel.
    let mut colors: Vec<Option<usize>> = vec![None; data.len()];
    let mut pairs: Vec<ExtremaPair> = Vec::with_capacity(data.len());

    for &index in sorted_indices {
        if index == 0 {
            match colors[index + 1] {
                None => create_component(data, &mut components, &mut colors, index),
                Some(color) => extend_component(&mut components, &mut colors, color, index),
            }
            continue;
        }
        if index == data.len() - 1 {
            match colors[index - 1] {
                None => create_component(data, &mut components, &mut colors, index),
                Some(color) => extend_component(&mut components, &mut colors, color, index),
            }
            continue;
        }

        match (colors[index - 1], colors[index + 1]) {
            (None, None) => create_component(data, &mut components, &mut colors, index),
            (Some(left), None) => extend_component(&mut components, &mut colors, left, index),
            (None, Some(right)) => extend_component(&mut components, &mut colors, right, index),
            (Some(left), Some(right)) => {
                // A local maximum: the pair keeps the minimum of whichever component is deeper,
                // compared with `<`, so a tie takes the LEFT one.
                if components[right].min_value < components[left].min_value {
                    pairs.push(ExtremaPair::new(data, components[left].min_index, index));
                } else {
                    pairs.push(ExtremaPair::new(data, components[right].min_index, index));
                }
                merge_components(&mut components, &mut colors, left, right, index);
            }
        }
    }

    // `Comparator.comparingDouble(persistence).reversed()`, on a stable sort, so equal
    // persistences keep the order the watershed found them in.
    pairs.sort_by(|first, second| java_compare(second.persistence, first.persistence));
    pairs
}

/// `createComponent`: a new colour, numbered by how many components exist, at a local minimum.
fn create_component(
    data: &[f64],
    components: &mut Vec<Component>,
    colors: &mut [Option<usize>],
    index: usize,
) {
    colors[index] = Some(components.len());
    components.push(Component {
        left_index: index,
        right_index: index,
        min_index: index,
        min_value: data[index],
    });
}

/// `extendComponent`: extend a component to a point beside it, which only moves the end it
/// actually touches. The reference's second test is `index - 1 == rightIndex`; the port writes
/// `right_index + 1 == index` because `index` is unsigned and `0 - 1` would not be `-1`. The two
/// agree on every index the reference can see.
fn extend_component(
    components: &mut [Component],
    colors: &mut [Option<usize>],
    component_index: usize,
    index: usize,
) {
    if index + 1 == components[component_index].left_index {
        components[component_index].left_index = index;
    } else if components[component_index].right_index + 1 == index {
        components[component_index].right_index = index;
    }
    colors[index] = Some(component_index);
}

/// Merge at a local maximum, keeping the component with the lower minimum and, on a tie, the lower
/// **colour**, which is the earlier component and not the earlier index.
fn merge_components(
    components: &mut [Component],
    colors: &mut [Option<usize>],
    left_color: usize,
    right_color: usize,
    index: usize,
) {
    let (keep, merge) = if components[left_color].min_value < components[right_color].min_value {
        (left_color, right_color)
    } else if components[left_color].min_value > components[right_color].min_value {
        (right_color, left_color)
    } else if left_color < right_color {
        (left_color, right_color)
    } else {
        (right_color, left_color)
    };

    let (merged_left, merged_right, merged_min) = (
        components[merge].left_index,
        components[merge].right_index,
        components[merge].min_index,
    );
    colors[merged_left] = Some(keep);
    colors[merged_right] = Some(keep);
    if components[keep].min_index > merged_min {
        components[keep].left_index = merged_left;
    } else {
        components[keep].right_index = merged_right;
    }
    // The maximum takes the colour of the point on its left, read AFTER the merge above, so which
    // component that is depends on what the merge just rewrote.
    colors[index] = colors[index - 1];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double_compare_is_not_the_ordinary_one() {
        assert_eq!(java_compare(-0.0, 0.0), Ordering::Less);
        assert_eq!(java_compare(f64::NAN, f64::INFINITY), Ordering::Greater);
        assert_eq!(java_compare(f64::NAN, f64::NAN), Ordering::Equal);
        assert_eq!(java_compare(1.0, 2.0), Ordering::Less);
    }

    #[test]
    fn a_persistence_can_be_negative_zero() {
        let answer = persistence_optimizer(&[0.0, -0.0, 0.0, -0.0]).expect("data");
        assert_eq!(answer.minima_indices, vec![1, 2]);
        assert_eq!(answer.persistences.len(), 2);
        assert!(answer.persistences[1].is_sign_negative() && answer.persistences[1] == 0.0);
    }

    #[test]
    fn a_nan_is_the_global_maximum() {
        let answer = persistence_optimizer(&[1.0, f64::NAN, 0.0, 2.0, 0.5]).expect("data");
        assert_eq!(answer.minima_indices, vec![2, 0, 4]);
        assert!(answer.persistences[0].is_nan());
        assert_eq!(answer.persistences[2], 1.5);
    }

    #[test]
    fn a_plateau_is_a_minimum_at_its_left_end() {
        let answer = persistence_optimizer(&[2.0, 2.0, 2.0, 5.0, 1.0]).expect("data");
        assert_eq!(answer.minima_indices, vec![4, 0]);
    }

    #[test]
    fn one_point_has_no_pairs_and_no_data_is_refused() {
        let answer = persistence_optimizer(&[1.0]).expect("data");
        assert_eq!(answer.minima_indices, vec![0]);
        assert_eq!(answer.persistences, vec![0.0]);
        assert_eq!(
            persistence_optimizer(&[]).unwrap_err().message(),
            "Data must contain at least one element."
        );
    }
}
