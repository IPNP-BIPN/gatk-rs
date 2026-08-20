//! `GatherTranches`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.GatherTranches` and `VQSLODTranche`
//! (GATK 4.6.2.0).
//!
//! The gather that ends a scattered `VariantRecalibrator` run: one VQSLOD tranche file per shard
//! in, one truth-sensitivity tranche file out.
//!
//! # A merged level's Ti/Tv is a ratio of sums
//!
//! ```java
//! double trancheKnownTransitions = (tranche.knownTiTv * tranche.numKnown) / (1 + tranche.knownTiTv);
//! sumKnownTransitions += trancheKnownTransitions;
//! sumKnownTransversions += (tranche.numKnown - trancheKnownTransitions);
//! ```
//!
//! Each shard's ratio is turned back into a count of transitions and a count of transversions, the
//! two are summed across shards, and the merged ratio is their quotient. A mean of the shards'
//! ratios would answer something else. A level whose shards have no known variants sums two zeroes
//! and the quotient is NaN, which the writer prints as `NaN`.
//!
//! # The sensitivity match is a walk that answers with the previous level
//!
//! The merged levels are visited in DESCENDING VQSLOD order, which is ascending sensitivity. For
//! each requested sensitivity the walk measures the distance to the current level and, the moment
//! that distance GROWS, writes the level BEFORE it. So the answer is the last level that was still
//! getting closer, and the first merged level is never an answer on its own: it is consumed as the
//! initial `currentTranche` and can only appear as somebody's `prevTranche`.
//!
//! The walk stops at the first target it cannot advance past, so asking for more sensitivities than
//! there are levels writes fewer rows than were asked for, without a word.

use gatk_engine::java_format::format_decimals;
use gatk_engine::tranches::{Mode, TrancheError, TruthSensitivityTranche, EXPECTED_COLUMN_COUNT};

/// `VQSLODTranche.CURRENT_VERSION`, which the INPUT files must carry.
pub const INPUT_VERSION: i32 = 6;

/// `TruthSensitivityTranche.CURRENT_VERSION`, which the OUTPUT carries. Not the same number.
pub const OUTPUT_VERSION: i32 = 5;

/// `GatherTranches`'s default truth-sensitivity levels.
pub const DEFAULT_TRUTH_SENSITIVITY_LEVELS: [f64; 2] = [90.0, 99.0];

/// One row of a VQSLOD tranche file.
#[derive(Debug, Clone, PartialEq)]
pub struct VqslodTranche {
    pub min_vqslod: f64,
    pub num_known: i64,
    pub known_titv: f64,
    pub num_novel: i64,
    pub novel_titv: f64,
    pub accessible_truth_sites: i32,
    pub calls_at_truth_sites: i32,
    pub model: Mode,
    pub name: String,
}

impl VqslodTranche {
    /// `getTruthSensitivity()`, zero rather than a division when nothing is accessible.
    pub fn truth_sensitivity(&self) -> f64 {
        if self.accessible_truth_sites > 0 {
            self.calls_at_truth_sites as f64 / (1.0 * self.accessible_truth_sites as f64)
        } else {
            0.0
        }
    }
}

/// `VQSLODTranche.readTranches`.
///
/// The version is checked on the comment line and nowhere else: a file whose columns are right and
/// whose version line says anything but six is refused before a single row is read.
pub fn read_vqslod_tranches(text: &str) -> Result<Vec<VqslodTranche>, TrancheError> {
    let mut header: Option<Vec<String>> = None;
    let mut tranches = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            if !line.contains("Version") {
                continue;
            }
            // `line.split("\\s+")[3]`, which is the number after `# Version number`.
            let words: Vec<&str> = line.split_whitespace().collect();
            let version = words
                .get(3)
                .and_then(|word| word.parse::<i32>().ok())
                .unwrap_or(-1);
            if version != INPUT_VERSION {
                return Err(TrancheError::VqslodVersion {
                    found: words.get(3).unwrap_or(&"").to_string(),
                    expected: INPUT_VERSION,
                });
            }
            continue;
        }
        let values: Vec<&str> = line.split(',').collect();
        let Some(names) = &header else {
            if values.len() != EXPECTED_COLUMN_COUNT {
                return Err(TrancheError::HeaderLength {
                    file: String::new(),
                    line: line.to_string(),
                });
            }
            header = Some(values.iter().map(|v| (*v).to_string()).collect());
            continue;
        };
        if names.len() != values.len() {
            return Err(TrancheError::RowLength {
                file: String::new(),
                header: names.len(),
                values: values.len(),
                line: line.to_string(),
            });
        }
        let get = |key: &str| -> Option<&str> {
            names
                .iter()
                .position(|name| name == key)
                .map(|index| values[index])
        };
        let number = |key: &str, default: f64| -> f64 {
            get(key).and_then(|v| v.parse().ok()).unwrap_or(default)
        };
        tranches.push(VqslodTranche {
            min_vqslod: number("minVQSLod", f64::NAN),
            num_known: number("numKnown", -1.0) as i64,
            known_titv: number("knownTiTv", -1.0),
            num_novel: number("numNovel", f64::NAN) as i64,
            novel_titv: number("novelTiTv", f64::NAN),
            accessible_truth_sites: number("accessibleTruthSites", -1.0) as i32,
            calls_at_truth_sites: number("callsAtTruthSites", -1.0) as i32,
            model: Mode::value_of(get("model"))?,
            name: get("filterName").unwrap_or_default().to_string(),
        });
    }
    Ok(tranches)
}

/// `mergeAndConvertTranches(list, mode)`: the shards of ONE VQSLOD level, summed.
pub fn merge_level(shards: &[VqslodTranche], mode: Mode) -> VqslodTranche {
    let index = shards[0].min_vqslod;
    let mut known = 0i64;
    let mut known_transitions = 0.0;
    let mut known_transversions = 0.0;
    let mut novel = 0i64;
    let mut novel_transitions = 0.0;
    let mut novel_transversions = 0.0;
    let mut accessible = 0i32;
    let mut called = 0i32;
    for shard in shards {
        known += shard.num_known;
        let transitions = (shard.known_titv * shard.num_known as f64) / (1.0 + shard.known_titv);
        known_transitions += transitions;
        known_transversions += shard.num_known as f64 - transitions;
        novel += shard.num_novel;
        let novel_ti = (shard.novel_titv * shard.num_novel as f64) / (1.0 + shard.novel_titv);
        novel_transitions += novel_ti;
        novel_transversions += shard.num_novel as f64 - novel_ti;
        accessible += shard.accessible_truth_sites;
        called += shard.calls_at_truth_sites;
    }
    VqslodTranche {
        min_vqslod: index,
        num_known: known,
        known_titv: known_transitions / known_transversions,
        num_novel: novel,
        novel_titv: novel_transitions / novel_transversions,
        accessible_truth_sites: accessible,
        calls_at_truth_sites: called,
        model: mode,
        // `"gathered" + indexVQSLOD`, which is Java's `Double.toString` of the level.
        name: format!(
            "gathered{}",
            gatk_engine::tsv_table::java_double_to_string(index)
        ),
    }
}

/// `mergeAndConvertTranches(map, tsLevels, mode)`: the whole gather.
///
/// `levels` is sorted here, as the reference sorts the caller's list in place.
pub fn merge_and_convert(
    shards: &[VqslodTranche],
    levels: &[f64],
    mode: Mode,
) -> Vec<TruthSensitivityTranche> {
    // A `TreeMap<Double, List<VQSLODTranche>>` walked by `descendingKeySet`.
    let mut keys: Vec<f64> = Vec::new();
    for shard in shards {
        if !keys.contains(&shard.min_vqslod) {
            keys.push(shard.min_vqslod);
        }
    }
    keys.sort_by(|left, right| right.partial_cmp(left).expect("no NaN level"));
    let merged: Vec<VqslodTranche> = keys
        .iter()
        .map(|key| {
            let level: Vec<VqslodTranche> = shards
                .iter()
                .filter(|shard| shard.min_vqslod == *key)
                .cloned()
                .collect();
            merge_level(&level, mode)
        })
        .collect();

    let mut targets = levels.to_vec();
    targets.sort_by(|left, right| left.partial_cmp(right).expect("no NaN level"));

    let mut gathered = Vec::new();
    if merged.is_empty() || targets.is_empty() {
        return gathered;
    }
    let mut target_index = 0;
    let mut target = targets[0];
    let mut sensitivity_delta = 100.0;
    let mut current = 0;
    while current + 1 < merged.len() {
        let previous_delta = sensitivity_delta;
        let previous = current;
        current += 1;
        sensitivity_delta = (target - merged[current].truth_sensitivity() * 100.0).abs();
        if sensitivity_delta > previous_delta {
            gathered.push(convert(&merged[previous], target, mode));
            target_index += 1;
            if target_index < targets.len() {
                target = targets[target_index];
                sensitivity_delta = (target - merged[current].truth_sensitivity() * 100.0).abs();
            } else {
                break;
            }
        }
        // The reference tests this INSIDE the loop, so the last level answers the target that was
        // still open when the levels ran out.
        if current + 1 >= merged.len() {
            gathered.push(convert(&merged[current], target, mode));
        }
    }
    gathered
}

fn convert(tranche: &VqslodTranche, target: f64, mode: Mode) -> TruthSensitivityTranche {
    TruthSensitivityTranche {
        target_truth_sensitivity: target,
        min_vqslod: tranche.min_vqslod,
        num_known: tranche.num_known,
        known_titv: tranche.known_titv,
        num_novel: tranche.num_novel,
        novel_titv: tranche.novel_titv,
        accessible_truth_sites: tranche.accessible_truth_sites,
        calls_at_truth_sites: tranche.calls_at_truth_sites,
        model: mode,
        name: String::new(),
    }
}

/// `TruthSensitivityTranche.printHeader()`, whose version is five however new the input was.
pub fn print_header() -> String {
    format!(
        "# Variant quality score tranches file\n# Version number {OUTPUT_VERSION}\n\
         targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,\
         accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
    )
}

/// `Tranche.tranchesString`: one row per tranche, each naming the band it closes.
///
/// The filter name is built from the PREVIOUS row's target sensitivity, or zero for the first, so
/// the names read as a partition of `0.00` to the last target.
pub fn tranches_string(tranches: &[TruthSensitivityTranche]) -> String {
    let mut text = String::new();
    let mut previous = 0.0;
    for tranche in tranches {
        text.push_str(&format!(
            "{},{},{},{},{},{},VQSRTranche{}{}to{},{},{},{},{}\n",
            format_decimals(tranche.target_truth_sensitivity, 2),
            tranche.num_known,
            tranche.num_novel,
            format_decimals(tranche.known_titv, 4),
            format_decimals(tranche.novel_titv, 4),
            format_decimals(tranche.min_vqslod, 4),
            tranche.model.name(),
            format_decimals(previous, 2),
            format_decimals(tranche.target_truth_sensitivity, 2),
            tranche.model.name(),
            tranche.accessible_truth_sites,
            tranche.calls_at_truth_sites,
            format_decimals(tranche.truth_sensitivity(), 4),
        ));
        previous = tranche.target_truth_sensitivity;
    }
    text
}

/// The whole run: the shards' text in, the gathered file out.
pub fn gather(shards: &[String], levels: &[f64], mode: Mode) -> Result<String, TrancheError> {
    let mut all = Vec::new();
    for text in shards {
        all.extend(read_vqslod_tranches(text)?);
    }
    let gathered = merge_and_convert(&all, levels, mode);
    Ok(format!("{}{}", print_header(), tranches_string(&gathered)))
}
