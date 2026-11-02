//! `TrainVariantAnnotationsModel`'s default backend without its interpreter: what
//! `isolation-forest.py` asks of scikit-learn 1.3.2 (BSD 3-Clause), which is `SimpleImputer` with
//! the median and `IsolationForest`.
//!
//! Every float the reference computes here is plain IEEE arithmetic except two: `numpy.log` of a
//! whole number, and `2 ** x`, which are libm's. See [`crate::libm_pow_two`] for the second, and
//! `docs/python-tools.md` for what was measured against the pinned image.
//!
//! * **The imputer** takes `numpy.ma.median`, a column's two middle values added and halved, or
//!   the one middle value added to itself and halved.
//! * **The forest** converts the imputed matrix to `float32`, and fits with `y` drawn from a
//!   `RandomState` of its own: the targets decide nothing but whether a node's impurity is zero,
//!   and they are written in the tree's `value` and `impurity` columns.
//! * **`BaseBagging._fit`** draws a seed per tree. Each seed gives the tree's `random_state`
//!   through a fresh `RandomState` and the tree's sample through another, whose drop-outs get
//!   weight zero rather than being removed, so a node's samples stay in row order.
//! * **`DepthFirstTreeBuilder` and `RandomSplitter`** are single-threaded per tree and draw from
//!   scikit-learn's own `our_rand_r`, seeded from the tree's `RandomState`. With `n_jobs` the
//!   trees run in other threads and land at the same index, so the thread count moves no byte.
//! * **`score_samples`** sums, per row and in tree order, the leaf's depth plus the average path
//!   length of its samples less one.
use crate::numpy_random::{sample_without_replacement, RandomState};

/// `sklearn.tree._utils.RAND_R_MAX`.
const RAND_R_MAX: u32 = 0x7fff_ffff;
/// `numpy.iinfo(numpy.int32).max`, `MAX_INT` in `_bagging.py`.
const MAX_INT: u64 = 0x7fff_ffff;
/// `numpy.euler_gamma`.
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
/// `numpy.finfo('double').eps`, `EPSILON` in `_tree.pyx`.
const EPSILON: f64 = f64::EPSILON;
/// `FEATURE_THRESHOLD` in `_splitter.pyx`.
const FEATURE_THRESHOLD: f64 = 1e-7;

/// `_tree.TREE_LEAF`.
const TREE_LEAF: i64 = -1;
/// `_tree.TREE_UNDEFINED`.
const TREE_UNDEFINED: i64 = -2;

/// The arguments of `IsolationForest(**hyperparameters)` this port reaches.
#[derive(Debug, Clone, PartialEq)]
pub struct Hyperparameters {
    pub random_state: u32,
    pub n_estimators: usize,
    /// Written in the pickle, and otherwise moves nothing: measured.
    pub n_jobs: Option<i64>,
}

/// One row of `Tree.nodes`.
#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    pub left_child: i64,
    pub right_child: i64,
    pub feature: i64,
    pub threshold: f64,
    pub impurity: f64,
    pub n_node_samples: i64,
    pub weighted_n_node_samples: f64,
    /// `Tree.value`, one output and one class per node.
    pub value: f64,
}

/// A fitted `ExtraTreeRegressor`.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    /// `get_params()['random_state']`.
    pub random_state: u32,
    pub max_depth: i64,
    pub nodes: Vec<Node>,
}

/// What scoring needs: a fitted forest, whether fitted here or read from a pickle.
#[derive(Debug, Clone, PartialEq)]
pub struct Scorer {
    pub annotation_names: Vec<String>,
    /// `SimpleImputer.statistics_`.
    pub statistics: Vec<f64>,
    /// `_max_samples`.
    pub max_samples: i64,
    /// Per tree: `feature` and `threshold` by node, children by node.
    pub trees: Vec<ScoringTree>,
    /// `_average_path_length_per_tree`.
    pub average_path_lengths: Vec<Vec<f64>>,
    /// `_decision_path_lengths`.
    pub decision_path_lengths: Vec<Vec<i64>>,
}

impl Hyperparameters {
    /// The JSON the script hands `IsolationForest(**kwargs)`, for the keyword arguments whose
    /// effect on the bytes has been measured. An unseeded forest has no bytes to match, and the
    /// rest change a fit this port does not carry, or the pickle's shape: the error names the one.
    pub fn from_json(text: &str) -> Result<Hyperparameters, String> {
        let parsed: serde_json::Value =
            serde_json::from_str(text).map_err(|error| format!("Hyperparameters JSON: {error}"))?;
        let object = parsed
            .as_object()
            .ok_or("Hyperparameters JSON that is not an object")?;
        let mut out = Hyperparameters {
            random_state: 0,
            n_estimators: 100,
            n_jobs: None,
        };
        let mut seeded = false;
        for (key, value) in object {
            let unsupported = || format!("IsolationForest({key}={value})");
            match key.as_str() {
                "random_state" => {
                    out.random_state = value
                        .as_u64()
                        .and_then(|seed| u32::try_from(seed).ok())
                        .ok_or_else(unsupported)?;
                    seeded = true;
                }
                "n_estimators" => {
                    out.n_estimators = value
                        .as_u64()
                        .filter(|count| (1..=100_000).contains(count))
                        .ok_or_else(unsupported)? as usize;
                }
                "n_jobs" => {
                    if !value.is_null() {
                        out.n_jobs =
                            Some(value.as_i64().filter(|jobs| *jobs != 0).ok_or_else(unsupported)?);
                    }
                }
                "verbose" if value.as_i64() == Some(0) => {}
                "max_samples" | "contamination" if value.as_str() == Some("auto") => {}
                "bootstrap" | "warm_start" if value.as_bool() == Some(false) => {}
                "max_features" if value.as_f64() == Some(1.0) && value.is_f64() => {}
                _ => return Err(unsupported()),
            }
        }
        if !seeded {
            return Err("IsolationForest without a random_state".to_string());
        }
        Ok(out)
    }
}

/// The columns of `Tree.nodes` that `apply` reads.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoringTree {
    pub left_child: Vec<i64>,
    pub right_child: Vec<i64>,
    pub feature: Vec<i64>,
    pub threshold: Vec<f64>,
}

/// The whole fit, which the pickle holds.
#[derive(Debug, Clone, PartialEq)]
pub struct Forest {
    pub hyperparameters: Hyperparameters,
    pub annotation_names: Vec<String>,
    pub statistics: Vec<f64>,
    pub n_samples: usize,
    pub max_samples: usize,
    pub max_depth: i64,
    pub seeds: Vec<i64>,
    pub trees: Vec<Tree>,
    pub average_path_lengths: Vec<Vec<f64>>,
    pub decision_path_lengths: Vec<Vec<i64>>,
}

impl Forest {
    pub fn scorer(&self) -> Scorer {
        Scorer {
            annotation_names: self.annotation_names.clone(),
            statistics: self.statistics.clone(),
            max_samples: self.max_samples as i64,
            trees: self
                .trees
                .iter()
                .map(|tree| ScoringTree {
                    left_child: tree.nodes.iter().map(|n| n.left_child).collect(),
                    right_child: tree.nodes.iter().map(|n| n.right_child).collect(),
                    feature: tree.nodes.iter().map(|n| n.feature).collect(),
                    threshold: tree.nodes.iter().map(|n| n.threshold).collect(),
                })
                .collect(),
            average_path_lengths: self.average_path_lengths.clone(),
            decision_path_lengths: self.decision_path_lengths.clone(),
        }
    }
}

/// `SimpleImputer(strategy='median').fit`: a column with no value at all is refused, which is
/// the shape assertion after `fit_transform` in the script.
pub fn median_statistics(rows: &[Vec<f64>]) -> Result<Vec<f64>, String> {
    let columns = rows.first().map_or(0, Vec::len);
    let mut out = Vec::with_capacity(columns);
    for column in 0..columns {
        let mut values: Vec<f64> = rows
            .iter()
            .map(|row| row[column])
            .filter(|value| !value.is_nan())
            .collect();
        if values.iter().any(|value| value.is_infinite()) {
            return Err("Input X contains infinity or a value too large for dtype('float64').".to_string());
        }
        if values.is_empty() {
            return Err("Shape of imputed annotations differs from shape of raw annotations; at \
                 least one feature is completely missing and hence dropped during imputation."
                .to_string());
        }
        values.sort_by(f64::total_cmp);
        let half = values.len() / 2;
        let high = values[half];
        let low = if values.len() % 2 == 1 {
            high
        } else {
            values[half - 1]
        };
        out.push((low + high) / 2.0);
    }
    Ok(out)
}

/// `SimpleImputer.transform`, then `_validate_data(dtype=float32)`: flat, by row.
fn imputed_float32(rows: &[Vec<f64>], statistics: &[f64]) -> Result<Vec<f32>, String> {
    let mut out = Vec::with_capacity(rows.len() * statistics.len());
    for row in rows {
        for (value, statistic) in row.iter().zip(statistics) {
            let imputed = if value.is_nan() { *statistic } else { *value };
            let narrow = imputed as f32;
            if !narrow.is_finite() {
                return Err(
                    "Input X contains infinity or a value too large for dtype('float32')."
                        .to_string(),
                );
            }
            out.push(narrow);
        }
    }
    Ok(out)
}

/// `_iforest._average_path_length`, one element at a time.
fn average_path_length(n: i64) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    if n == 2 {
        return 1.0;
    }
    let log = f64::from_bits(LOG_OF_WHOLE[n as usize - 2]);
    let n = n as f64;
    2.0 * (log + EULER_GAMMA) - 2.0 * (n - 1.0) / n
}

/// `numpy.log(k)` for the whole numbers `1..=256`, as bits: glibc's `log`, measured in the reference's
/// image. `max_samples='auto'` caps a tree at 256 samples, so nothing else is ever asked; the FMA and
/// the SSE2 variants of glibc agree on every whole number to 200000. fdlibm does not, by one unit in
/// the last place at 3, 48, 74 and 196 among these, and that is enough to change a scorer's bytes.
const LOG_OF_WHOLE: [u64; 256] = [
    0x0000000000000000, 0x3fe62e42fefa39ef, 0x3ff193ea7aad030b, 0x3ff62e42fefa39ef,
    0x3ff9c041f7ed8d33, 0x3ffcab0bfa2a2002, 0x3fff2272ae325a57, 0x4000a2b23f3bab73,
    0x400193ea7aad030b, 0x40026bb1bbb55516, 0x40032ee3b77f374c, 0x4003e116bcd39e7d,
    0x400485042b318c51, 0x40051cca16d7bba7, 0x4005aa16394d481f, 0x40062e42fefa39ef,
    0x4006aa6bc1fa7f7a, 0x40071f7b3a6b9186, 0x40078e360604b32c, 0x4007f7427b73e391,
    0x40085b2e946faeb1, 0x4008ba74773dc5c8, 0x4009157dfdd1b3f0, 0x40096ca77c922cf9,
    0x4009c041f7ed8d33, 0x400a1094eaf01acd, 0x400a5ddfb8038490, 0x400aa85ad6964a23,
    0x400af038cbcdfe1f, 0x400b35a6f90bd69b, 0x400b78ce48912b5a, 0x400bb9d3beb8c86b,
    0x400bf8d8f4d5b8d1, 0x400c35fc81b90df6, 0x400c715a530ff3c5, 0x400cab0bfa2a2002,
    0x400ce328ee3ac7bd, 0x400d19c6c5c341a8, 0x400d4ef968880dd6, 0x400d82d33b32720d,
    0x400db565458eac4e, 0x400de6bf542e3d2d, 0x400e16f01614ec11, 0x400e460536fc5443,
    0x400e740b76a3c9a4, 0x400ea10ebd90426c, 0x400ecd1a2f90f0b6, 0x400ef8383c50bb75,
    0x400f2272ae325a57, 0x400f4bd2b7ac1baf, 0x400f7460ff5100ff, 0x400f9c25aaaea949,
    0x400fc32868223cab, 0x400fe97077c2130c, 0x4010078259bafef3, 0x401019f5cb2a6c50,
    0x40102c15a1ad9a59, 0x40103de4c5c6464d, 0x40104f65f9c72979, 0x4010609bdc65328b,
    0x40107188eb126c0d, 0x4010822f8427dceb, 0x40109291e8e3181b, 0x4010a2b23f3bab73,
    0x4010b29293942975, 0x4010c234da4a23a6, 0x4010d19af128112d, 0x4010e0c6a0bbce39,
    0x4010efb99d941abb, 0x4010fe7589674121, 0x40110cfbf424d6c9, 0x40111b4e5cf4573f,
    0x4011296e33222fe6, 0x4011375cd6fcab1d, 0x4011451b9aa2075c, 0x401152abc2c0e812,
    0x4011600e874c323c, 0x40116d4514234e29, 0x40117a5089afb173, 0x40118731fd788044,
    0x401193ea7aad030b, 0x4011a07b02a69d65, 0x4011ace48d62e482, 0x4011b92809f665d4,
    0x4011c5465ef8a30a, 0x4011d1406ae9bd46, 0x4011dd1704923fd2, 0x4011e8cafb5d7160,
    0x4011f45d17ae8af6, 0x4011ffce1b312c10, 0x40120b1ec1255cbe, 0x4012164fbea76874,
    0x40122161c2f3d670, 0x40122c5577a7bf99, 0x4012372b80fdbce3, 0x401241e47e07a4f8,
    0x40124c8108e54bf8, 0x40125701b6f8746a, 0x4012616719161d2b, 0x40126bb1bbb55516,
    0x401275e2271bba31, 0x40127ff8df87c7bd, 0x401289f6655913ec, 0x401293db35369be2,
    0x40129da7c8333aa5, 0x4012a75c93f06593, 0x4012b0fa0abf475b, 0x4012ba809bc050c4,
    0x4012c3f0b3015506, 0x4012cd4ab99a4631, 0x4012d68f15c8a4a1, 0x4012dfbe2b09b38d,
    0x4012e8d85a33835c, 0x4012f1de018ce197, 0x4012facf7ce43d45, 0x401303ad25a58d8b,
    0x40130c7752ef47ae, 0x4013152e59a670b7, 0x40131dd28c89d653, 0x401326643c4479c9,
    0x40132ee3b77f374c, 0x401337514af1b34b, 0x40133fad417296e9, 0x401347f7e4072429,
    0x4013503179f229e7, 0x4013585a48c25f59, 0x4013607294602e42, 0x4013687a9f1af2b1,
    0x40137072a9b5b6cb, 0x4013785af37370b3, 0x40138033ba22c862, 0x401387fd3a296ae4,
    0x40138fb7ae8ef02c, 0x401397635107586b, 0x40139f0059fd2595, 0x4013a68f009b1577,
    0x4013ae0f7ad5809d, 0x4013b581fd7361f9, 0x4013bce6bc170b12, 0x4013c43de946885e,
    0x4013cb87b673b91e, 0x4013d2c454041e07, 0x4013d9f3f15861ce, 0x4013e116bcd39e7d,
    0x4013e82ce3e2625c, 0x4013ef3693017724, 0x4013f633f5c46dee, 0x4013fd2536dbf25a,
    0x4014040a801be72d, 0x40140ae3fa814e9a, 0x401411b1ce380053, 0x4014187422a02f50,
    0x40141f2b1e53c142, 0x401425d6e72b797a, 0x40142c77a243f8fa, 0x4014330d74029567,
    0x40143998801a085c, 0x40144018e98ef8b1, 0x4014468ed2bc5f18, 0x40144cfa5d57c782,
    0x4014535baa75708e, 0x401459b2da8c4a48, 0x401460000d79d57b, 0x401466436285e4a3,
    0x40146c7cf8663fb5, 0x401472aced422bc0, 0x401478d35eb5d75f, 0x40147ef069d5ad12,
    0x401485042b318c51, 0x40148b0ebed7ea48, 0x401491104058db1c, 0x40149708cac90484,
    0x40149cf878c47a87, 0x4014a2df64718710, 0x4014a8bda7835d30, 0x4014ae935b3cb89d,
    0x4014b46098726a3b, 0x4014ba25778dd234, 0x4014bfe2108f4858, 0x4014c5967b10734e,
    0x4014cb42ce468f2b, 0x4014d0e72104a3fc, 0x4014d68389bdacd0, 0x4014dc181e86afb2,
    0x4014e1a4f518c72b, 0x4014e72a22d31dae, 0x4014eca7bcbcdb63, 0x4014f21dd78706d7,
    0x4014f78c878e58de, 0x4014fcf3e0dd0421, 0x40150253f72c70ba, 0x401507acdde6ec36,
    0x40150cfea8294e62, 0x4015124968c49336, 0x4015178d323f6a38, 0x40151cca16d7bba7,
    0x40152200288423c1, 0x4015272f78f56469, 0x40152c581997cd86, 0x4015317a1b949c53,
    0x401536958fd351f0, 0x40153baa86fb016f, 0x401540b9117395a5, 0x401545c13f670efb,
    0x40154ac320c2b974, 0x40154fbec5385b2a, 0x401554b43c3f5b7d, 0x401559a39515e320,
    0x40155e8cdec1f53c, 0x40156370281281e3, 0x4015684d7fa071fb, 0x40156d24f3cfacd1,
    0x401571f692d0178c, 0x401576c26a9e8e99, 0x40157b888905d955, 0x40158048fb9f9802,
    0x40158503cfd52c43, 0x401589b912e09c44, 0x40158e68d1cd70a9, 0x4015931319798d6f,
    0x401597b7f69605e5, 0x40159c5775a7ebdf, 0x4015a0f1a3091a36, 0x4015a5868ae8facb,
    0x4015aa16394d481f, 0x4015aea0ba12ca9a, 0x4015b32618ee11ab, 0x4015b7a6616c28d5,
    0x4015bc219ef348c3, 0x4015c097dcc38483, 0x4015c50925f772fe, 0x4015c9758584d4c9,
    0x4015cddd063d365f, 0x4015d23fb2ce8eec, 0x4015d69d95c3dba8, 0x4015daf6b985b7f5,
    0x4015df4b285af236, 0x4015e39aec691d91, 0x4015e7e60fb5209a, 0x4015ec2c9c23c107,
    0x4015f06e9b7a2c7e, 0x4015f4ac175e7e8a, 0x4015f8e5195843cd, 0x4015fd19aad0fa89,
    0x40160149d5149079, 0x40160575a151de27, 0x4016099d189b1fbf, 0x40160dc043e66b67,
    0x401611df2c0e2544, 0x401615f9d9d17124, 0x40161a1055d4a1e2, 0x40161e22a8a1a697,
    0x40162230daa8759e, 0x4016263af43f7580, 0x40162a40fda3e3cc, 0x40162e42fefa39ef,
];

/// `our_rand_r`.
fn our_rand_r(seed: &mut u32) -> u32 {
    if *seed == 0 {
        *seed = 1;
    }
    *seed ^= *seed << 13;
    *seed ^= *seed >> 17;
    *seed ^= *seed << 5;
    // The `% (RAND_R_MAX + 1)` of a `uint32` whose modulus is 2^31.
    *seed % (RAND_R_MAX + 1)
}

/// `_splitter.pyx`'s and `_criterion.pyx`'s state over one tree.
struct Builder<'a> {
    x: &'a [f32],
    columns: usize,
    y: &'a [f64],
    samples: Vec<usize>,
    features: Vec<usize>,
    constant_features: Vec<usize>,
    rand_r_state: u32,
    weighted_n_samples: f64,
    // The criterion.
    sum_total: f64,
    sq_sum_total: f64,
    weighted_n_node: f64,
    sum_left: f64,
    sum_right: f64,
    weighted_left: f64,
    weighted_right: f64,
    start: usize,
    end: usize,
    position: usize,
}

struct Split {
    position: usize,
    feature: i64,
    threshold: f64,
    improvement: f64,
    impurity_left: f64,
    impurity_right: f64,
}

impl Builder<'_> {
    /// `RegressionCriterion.init`, with every weight one.
    fn node_reset(&mut self, start: usize, end: usize) {
        self.start = start;
        self.end = end;
        self.sq_sum_total = 0.0;
        self.sum_total = 0.0;
        self.weighted_n_node = 0.0;
        for p in start..end {
            let y = self.y[self.samples[p]];
            let w_y = 1.0 * y;
            self.sum_total += w_y;
            self.sq_sum_total += w_y * y;
            self.weighted_n_node += 1.0;
        }
        self.reset();
    }

    fn reset(&mut self) {
        self.sum_left = 0.0;
        self.sum_right = self.sum_total;
        self.weighted_left = 0.0;
        self.weighted_right = self.weighted_n_node;
        self.position = self.start;
    }

    fn reverse_reset(&mut self) {
        self.weighted_right = 0.0;
        self.weighted_left = self.weighted_n_node;
        self.position = self.end;
        self.sum_right = 0.0;
        self.sum_left = self.sum_total;
    }

    /// `RegressionCriterion.update`, from whichever end is nearer.
    fn update(&mut self, new_position: usize) {
        if new_position - self.position <= self.end - new_position {
            for p in self.position..new_position {
                let y = self.y[self.samples[p]];
                self.sum_left += 1.0 * y;
                self.weighted_left += 1.0;
            }
        } else {
            self.reverse_reset();
            for p in (new_position..self.end).rev() {
                let y = self.y[self.samples[p]];
                self.sum_left -= 1.0 * y;
                self.weighted_left -= 1.0;
            }
        }
        self.weighted_right = self.weighted_n_node - self.weighted_left;
        self.sum_right = self.sum_total - self.sum_left;
        self.position = new_position;
    }

    /// `MSE.node_impurity`.
    fn node_impurity(&self) -> f64 {
        let mean = self.sum_total / self.weighted_n_node;
        self.sq_sum_total / self.weighted_n_node - mean * mean
    }

    /// `MSE.children_impurity`.
    fn children_impurity(&self) -> (f64, f64) {
        let mut sq_sum_left = 0.0;
        for p in self.start..self.position {
            let y = self.y[self.samples[p]];
            sq_sum_left += 1.0 * y * y;
        }
        let sq_sum_right = self.sq_sum_total - sq_sum_left;
        let mean_left = self.sum_left / self.weighted_left;
        let mean_right = self.sum_right / self.weighted_right;
        (
            sq_sum_left / self.weighted_left - mean_left * mean_left,
            sq_sum_right / self.weighted_right - mean_right * mean_right,
        )
    }

    /// `Criterion.impurity_improvement`.
    fn impurity_improvement(&self, parent: f64, left: f64, right: f64) -> f64 {
        (self.weighted_n_node / self.weighted_n_samples)
            * (parent
                - self.weighted_right / self.weighted_n_node * right
                - self.weighted_left / self.weighted_n_node * left)
    }

    /// `node_split_random`, with `max_features` of one and a leaf of at least one sample.
    fn node_split(&mut self, parent: f64, n_constant: &mut usize) -> Split {
        let (start, end) = (self.start, self.end);
        let mut best = Split {
            position: end,
            feature: 0,
            threshold: 0.0,
            improvement: f64::NEG_INFINITY,
            impurity_left: f64::INFINITY,
            impurity_right: f64::INFINITY,
        };
        let n_known = *n_constant;
        let mut n_total = n_known;
        let (mut n_found, mut n_drawn, mut n_visited) = (0usize, 0usize, 0usize);
        let mut f_i = self.columns;
        let mut found = false;
        while f_i > n_total && (n_visited < 1 || n_visited <= n_found + n_drawn) {
            n_visited += 1;
            let range = (f_i - n_found - n_drawn) as u32;
            let mut f_j = n_drawn + (our_rand_r(&mut self.rand_r_state) % range) as usize;
            if f_j < n_known {
                self.features.swap(n_drawn, f_j);
                n_drawn += 1;
                continue;
            }
            f_j += n_found;
            let feature = self.features[f_j];

            // `find_min_max`.
            let mut low = self.x[self.samples[start] * self.columns + feature];
            let mut high = low;
            for p in start + 1..end {
                let value = self.x[self.samples[p] * self.columns + feature];
                if value < low {
                    low = value;
                } else if value > high {
                    high = value;
                }
            }
            if f64::from(high) <= f64::from(low) + FEATURE_THRESHOLD {
                self.features[f_j] = self.features[n_total];
                self.features[n_total] = feature;
                n_found += 1;
                n_total += 1;
                continue;
            }
            f_i -= 1;
            self.features.swap(f_i, f_j);

            // `rand_uniform`.
            let (low, high) = (f64::from(low), f64::from(high));
            let draw = f64::from(our_rand_r(&mut self.rand_r_state));
            let mut threshold = (high - low) * draw / f64::from(RAND_R_MAX) + low;
            if threshold == high {
                threshold = low;
            }

            // `partition_samples`.
            let mut p = start;
            let mut partition_end = end;
            while p < partition_end {
                if f64::from(self.x[self.samples[p] * self.columns + feature]) <= threshold {
                    p += 1;
                } else {
                    partition_end -= 1;
                    self.samples.swap(p, partition_end);
                }
            }
            if partition_end == start || partition_end == end {
                continue;
            }
            found = true;
            best.position = partition_end;
            best.feature = feature as i64;
            best.threshold = threshold;
        }
        if found && best.position < end {
            self.reset();
            self.update(best.position);
            let (left, right) = self.children_impurity();
            best.impurity_left = left;
            best.impurity_right = right;
            best.improvement = self.impurity_improvement(parent, left, right);
        }
        // The invariant on `features[:n_known]` that siblings and children read.
        self.features[..n_known].copy_from_slice(&self.constant_features[..n_known]);
        let newly = n_known + n_found;
        let copied: Vec<usize> = self.features[n_known..newly].to_vec();
        self.constant_features[n_known..newly].copy_from_slice(&copied);
        *n_constant = n_total;
        best
    }
}

/// `DepthFirstTreeBuilder.build` over the rows `kept` marks.
fn build_tree(
    x: &[f32],
    columns: usize,
    y: &[f64],
    kept: &[bool],
    random_state: u32,
    max_depth: i64,
) -> Tree {
    let mut random = RandomState::new(random_state);
    let rand_r_state = random.randint(0, u64::from(RAND_R_MAX)) as u32;
    let samples: Vec<usize> = (0..kept.len()).filter(|i| kept[*i]).collect();
    let n = samples.len();
    let mut builder = Builder {
        x,
        columns,
        y,
        weighted_n_samples: n as f64,
        samples,
        features: (0..columns).collect(),
        constant_features: vec![0; columns],
        rand_r_state,
        sum_total: 0.0,
        sq_sum_total: 0.0,
        weighted_n_node: 0.0,
        sum_left: 0.0,
        sum_right: 0.0,
        weighted_left: 0.0,
        weighted_right: 0.0,
        start: 0,
        end: 0,
        position: 0,
    };
    struct Pending {
        start: usize,
        end: usize,
        depth: i64,
        parent: Option<(usize, bool)>,
        impurity: f64,
        n_constant: usize,
    }
    let mut stack = vec![Pending {
        start: 0,
        end: n,
        depth: 0,
        parent: None,
        impurity: f64::INFINITY,
        n_constant: 0,
    }];
    let mut nodes: Vec<Node> = Vec::new();
    let mut deepest = 0;
    let mut first = true;
    while let Some(mut pending) = stack.pop() {
        let (start, end) = (pending.start, pending.end);
        builder.node_reset(start, end);
        let mut is_leaf = pending.depth >= max_depth || end - start < 2;
        if first {
            pending.impurity = builder.node_impurity();
            first = false;
        }
        is_leaf = is_leaf || pending.impurity <= EPSILON;
        let mut split = None;
        if !is_leaf {
            let found = builder.node_split(pending.impurity, &mut pending.n_constant);
            is_leaf = found.position >= end || found.improvement + EPSILON < 0.0;
            split = Some(found);
        }
        let id = nodes.len();
        let value = builder.sum_total / builder.weighted_n_node;
        if let Some((parent, is_left)) = pending.parent {
            if is_left {
                nodes[parent].left_child = id as i64;
            } else {
                nodes[parent].right_child = id as i64;
            }
        }
        let leaf_or = |number: f64| if is_leaf { TREE_UNDEFINED as f64 } else { number };
        nodes.push(Node {
            left_child: TREE_LEAF,
            right_child: TREE_LEAF,
            feature: if is_leaf {
                TREE_UNDEFINED
            } else {
                split.as_ref().map_or(TREE_UNDEFINED, |s| s.feature)
            },
            threshold: leaf_or(split.as_ref().map_or(0.0, |s| s.threshold)),
            impurity: pending.impurity,
            n_node_samples: (end - start) as i64,
            weighted_n_node_samples: builder.weighted_n_node,
            value,
        });
        deepest = deepest.max(pending.depth);
        if !is_leaf {
            let split = split.expect("a split on a node that is no leaf");
            stack.push(Pending {
                start: split.position,
                end,
                depth: pending.depth + 1,
                parent: Some((id, false)),
                impurity: split.impurity_right,
                n_constant: pending.n_constant,
            });
            stack.push(Pending {
                start,
                end: split.position,
                depth: pending.depth + 1,
                parent: Some((id, true)),
                impurity: split.impurity_left,
                n_constant: pending.n_constant,
            });
        }
    }
    Tree {
        random_state,
        max_depth: deepest,
        nodes,
    }
}

/// `Tree.compute_node_depths`: the root is at one.
fn node_depths(nodes: &[Node]) -> Vec<i64> {
    let mut depths = vec![0i64; nodes.len()];
    depths[0] = 1;
    for (id, node) in nodes.iter().enumerate() {
        if node.left_child != TREE_LEAF {
            depths[node.left_child as usize] = depths[id] + 1;
            depths[node.right_child as usize] = depths[id] + 1;
        }
    }
    depths
}

/// `train` in `isolation-forest.py`: the imputer and `IsolationForest.fit`.
pub fn fit(
    annotation_names: &[String],
    rows: &[Vec<f64>],
    hyperparameters: &Hyperparameters,
) -> Result<Forest, String> {
    let statistics = median_statistics(rows)?;
    let x = imputed_float32(rows, &statistics)?;
    let columns = statistics.len();
    let n_samples = rows.len();

    let seed = hyperparameters.random_state;
    let mut for_y = RandomState::new(seed);
    let y: Vec<f64> = (0..n_samples).map(|_| for_y.next_double()).collect();

    let max_samples = n_samples.min(256);
    let bits = max_samples.max(2);
    // `int(ceil(log2(max_samples)))`, which is exact on the integers this can be.
    let max_depth = i64::from(usize::BITS - (bits - 1).leading_zeros());

    let mut for_seeds = RandomState::new(seed);
    let seeds: Vec<i64> = (0..hyperparameters.n_estimators)
        .map(|_| for_seeds.randint(0, MAX_INT) as i64)
        .collect();

    let mut trees = Vec::with_capacity(seeds.len());
    for tree_seed in &seeds {
        let tree_seed = *tree_seed as u32;
        // `_make_estimator`'s `_set_random_states`.
        let tree_random_state = RandomState::new(tree_seed).randint(0, MAX_INT) as u32;
        // `_generate_bagging_indices`: the features first, all of them, which draws nothing.
        let mut for_indices = RandomState::new(tree_seed);
        sample_without_replacement(columns, columns, &mut for_indices);
        let indices = sample_without_replacement(n_samples, max_samples, &mut for_indices);
        let mut kept = vec![false; n_samples];
        for index in indices {
            kept[index] = true;
        }
        trees.push(build_tree(
            &x,
            columns,
            &y,
            &kept,
            tree_random_state,
            max_depth,
        ));
    }
    let average_path_lengths = trees
        .iter()
        .map(|tree| {
            tree.nodes
                .iter()
                .map(|node| average_path_length(node.n_node_samples))
                .collect()
        })
        .collect();
    let decision_path_lengths = trees.iter().map(|tree| node_depths(&tree.nodes)).collect();
    Ok(Forest {
        hyperparameters: hyperparameters.clone(),
        annotation_names: annotation_names.to_vec(),
        statistics,
        n_samples,
        max_samples,
        max_depth,
        seeds,
        trees,
        average_path_lengths,
        decision_path_lengths,
    })
}

/// `scorer_lambda`: the annotation names checked, the rows imputed, and `clf.score_samples`.
pub fn score_samples(
    scorer: &Scorer,
    annotation_names: &[String],
    rows: &[Vec<f64>],
) -> Result<Vec<f64>, String> {
    if annotation_names != scorer.annotation_names {
        return Err(format!(
            "Input annotation names ({annotation_names:?}) must be identical to those used to \
             train the scorer ({:?}).",
            scorer.annotation_names
        ));
    }
    let x = imputed_float32(rows, &scorer.statistics)?;
    let columns = scorer.statistics.len();
    let denominator = scorer.trees.len() as f64 * average_path_length(scorer.max_samples);
    let mut depths = vec![0.0f64; rows.len()];
    for (index, tree) in scorer.trees.iter().enumerate() {
        let lengths = &scorer.average_path_lengths[index];
        let decisions = &scorer.decision_path_lengths[index];
        for (row, depth) in depths.iter_mut().enumerate() {
            let mut node = 0usize;
            while tree.left_child[node] != TREE_LEAF {
                let value = x[row * columns + tree.feature[node] as usize];
                node = if f64::from(value) <= tree.threshold[node] {
                    tree.left_child[node] as usize
                } else {
                    tree.right_child[node] as usize
                };
            }
            *depth += decisions[node] as f64 + lengths[node] - 1.0;
        }
    }
    depths
        .iter()
        .map(|depth| {
            // `np.divide(..., out=np.ones_like(depths), where=denominator != 0)`.
            let ratio = if denominator != 0.0 {
                depth / denominator
            } else {
                1.0
            };
            crate::libm_pow_two::pow_two(-ratio).map(|score| -score).ok_or_else(|| {
                format!("A score of 2 ** {} leaves the path of pow that is ported.", -ratio)
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_median_of_an_even_column_is_halved_and_a_missing_value_is_skipped() {
        let rows = vec![vec![1.0], vec![f64::NAN], vec![4.0], vec![2.0], vec![9.0]];
        assert_eq!(median_statistics(&rows).unwrap(), vec![3.0]);
        let rows = vec![vec![1.0], vec![f64::NAN], vec![4.0]];
        assert_eq!(median_statistics(&rows).unwrap(), vec![2.5]);
        let rows = vec![vec![7.0], vec![f64::NAN], vec![8.0], vec![9.0]];
        assert_eq!(median_statistics(&rows).unwrap(), vec![8.0]);
    }

    #[test]
    fn a_column_with_nothing_in_it_is_refused() {
        assert!(median_statistics(&[vec![f64::NAN]]).is_err());
    }

    #[test]
    fn the_path_length_of_few_samples_is_a_constant() {
        assert_eq!(average_path_length(1), 0.0);
        assert_eq!(average_path_length(2), 1.0);
    }

    #[test]
    fn our_rand_r_turns_a_zero_seed_into_one() {
        let (mut zero, mut one) = (0u32, 1u32);
        assert_eq!(our_rand_r(&mut zero), our_rand_r(&mut one));
    }

    #[test]
    fn the_forest_is_a_function_of_its_seed() {
        let rows: Vec<Vec<f64>> = (0..50)
            .map(|i| vec![f64::from(i), f64::from(i * i % 17)])
            .collect();
        let names = vec!["A".to_string(), "B".to_string()];
        let hyper = Hyperparameters {
            random_state: 3,
            n_estimators: 5,
            n_jobs: None,
        };
        let forest = fit(&names, &rows, &hyper).unwrap();
        assert_eq!(forest, fit(&names, &rows, &hyper).unwrap());
        let scores = score_samples(&forest.scorer(), &names, &rows).unwrap();
        assert!(scores.iter().all(|s| *s < 0.0 && *s > -1.0));
        assert!(score_samples(&forest.scorer(), &["A".to_string()], &rows).is_err());
    }
}
