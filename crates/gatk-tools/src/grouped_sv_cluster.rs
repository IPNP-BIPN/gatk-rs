//! `GroupedSVCluster`: the stratification engine of [`crate::sv_stratify`] choosing which threshold
//! set of [`crate::sv_cluster`] a record is clustered under.
//!
//! The tool is the two of them wired together, and the wiring is where its own behaviour lives:
//! one stratum per record and one threshold set per stratum, a record matching two strata refused
//! outright, and a record matching none written straight out without being clustered at all.
//!
//! Reading the VCF, collapsing a cluster into a representative record and the order the records
//! come out in are not ported. The tool disables its own output index because that order is not
//! guaranteed, so what is compared here is which records ended up together and under which stratum.
//!
//! The reference is `org.broadinstitute.hellbender.tools.walkers.sv.GroupedSVCluster` and
//! `org.broadinstitute.hellbender.tools.sv.cluster.StratifiedClusteringTableParser` (GATK
//! 4.7.0.0).
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a value or an error"; [`GroupedError`] models the exceptions the Java
//!   throws, and `?` after a call returns its error to the caller at once;
//! - `Option<T>` is "a value or nothing", Java's nullable reference checked by the compiler;
//! - `&[T]` is a borrowed, read-only view of a list; the caller keeps ownership;
//! - `.iter().find(..)`, `.map(..)`, `.collect()` are Java's stream, `filter(..).findFirst()`,
//!   `map` and `collect`; a `|x| ...` is a lambda.

use crate::sv_cluster::{cluster, Algorithm, CallRecord, ClusteringParameters, Linkage};
use crate::sv_stratify::{Engine, StratifyError, Thresholds, DEFAULT_STRATUM};

/// The columns of the clustering configuration, in the order `ImmutableSet.of` keeps them, which is
/// the order a missing one is reported in.
pub const COLUMN_NAMES: &[&str] = &[
    "NAME",
    "RECIPROCAL_OVERLAP",
    "SIZE_SIMILARITY",
    "BREAKEND_WINDOW",
    "SAMPLE_OVERLAP",
];

/// One row of the clustering configuration: `StratifiedClusteringTableParser.StratumParameters`.
#[derive(Debug, Clone, PartialEq)]
pub struct StratumParameters {
    pub name: String,
    pub reciprocal_overlap: f64,
    pub size_similarity: f64,
    pub breakend_window: i32,
    pub sample_overlap: f64,
}

/// What this tool refuses, on top of what the stratification engine refuses on its own.
///
/// Each variant is one of the Java's `throw` or `Utils.validate` sites; [`GroupedError::message`]
/// is the text it carries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroupedError {
    // `"Missing column " + column`, for the first of `COLUMN_NAMES` that is absent.
    MissingColumn {
        column: String,
    },
    /// The same self-referential message the stratification parser prints: both halves read the
    /// count of columns that were actually found.
    ColumnCount {
        count: usize,
    },
    /// `onTraversalStart`: the stratification configuration defined nothing.
    NoStrata,
    /// `onTraversalStart`: `stratEngine.getStrata().size() == clusterEngineMap.size()` failed.
    GroupCountMismatch,
    /// `onTraversalStart`: a stratum the clustering configuration does not name.
    GroupNotFound {
        name: String,
    },
    /// A different message from `SVStratify`'s, which offers a flag instead.
    MultipleMatches {
        id: String,
        names: Vec<String>,
    },
    /// Anything the stratification engine itself refused first.
    Stratify(StratifyError),
}

impl GroupedError {
    /// The exception's message, character for character the Java's string concatenation.
    pub fn message(&self) -> String {
        match self {
            GroupedError::MissingColumn { column } => format!("Missing column {column}"),
            // `{count}` inside `format!` is the variable `count`, interpolated. The Java reads
            // `columns.columnCount()` for both halves, so both print the same number:
            // `"Expected " + columns.columnCount() + " columns but found "
            // + columns.columnCount()`.
            GroupedError::ColumnCount { count } => {
                format!("Expected {count} columns but found {count}")
            }
            GroupedError::NoStrata => "No strata defined with --stratify-config".to_string(),
            GroupedError::GroupCountMismatch => {
                "Stratification and clustering configurations have a different number of groups."
                    .to_string()
            }
            GroupedError::GroupNotFound { name } => {
                format!("Could not find group {name} in clustering configuration.")
            }
            GroupedError::MultipleMatches { id, names } => format!(
                "Record {id} matched multiple groups: {}. Groups must be mutually exclusive. \
                 Please modify the group configurations and/or tracks so that no variant can \
                 match more than one group.",
                names.join(", ")
            ),
            GroupedError::Stratify(error) => error.message(),
        }
    }
}

/// `StratifiedClusteringTableParser.tableParser`: every column must be present, and there must be
/// no others.
///
/// ```java
/// for (final String column : COLUMN_NAMES) {
///     if (!columns.contains(column)) {
///         throw exceptionFactory.apply("Missing column " + column);
///     }
/// }
/// if (columns.columnCount() != COLUMN_NAMES.size()) {
/// ```
///
/// The presence check runs first, so a table missing a column and carrying an extra one reports
/// the missing column, not the count.
pub fn check_columns(columns: &[String]) -> Result<(), GroupedError> {
    for name in COLUMN_NAMES {
        if !columns.iter().any(|column| column == name) {
            return Err(GroupedError::MissingColumn {
                column: (*name).to_string(),
            });
        }
    }
    if columns.len() != COLUMN_NAMES.len() {
        return Err(GroupedError::ColumnCount {
            count: columns.len(),
        });
    }
    Ok(())
}

/// One row becomes all three parameter sets. Unlike `SVCluster`, depth-only, PESR and mixed pairs
/// cannot be given different thresholds here, so which set a pair takes stops mattering.
///
/// `readStrataClusteringConfig` says so, `// Identical parameters for each linkage type`, and
/// passes the row's same four values, in the same order, to `createPesrParameters`,
/// `createMixedParameters` and `createDepthParameters`.
pub fn linkage_for(parameters: &StratumParameters, cluster_del_with_dup: bool) -> Linkage {
    let (overlap, similarity, window, sample) = (
        parameters.reciprocal_overlap,
        parameters.size_similarity,
        parameters.breakend_window,
        parameters.sample_overlap,
    );
    Linkage {
        depth: ClusteringParameters::depth(overlap, similarity, window, sample),
        mixed: ClusteringParameters::mixed(overlap, similarity, window, sample),
        pesr: ClusteringParameters::pesr(overlap, similarity, window, sample),
        cluster_del_with_dup,
    }
}

/// The engines a run holds, keyed by group name.
///
/// The Java side is a `HashMap`, so a clustering configuration naming the same group twice keeps
/// only the last row and is one engine short of its own row count.
#[derive(Debug, Clone, PartialEq)]
pub struct Engines {
    pub entries: Vec<(String, StratumParameters)>,
}

impl Engines {
    /// `readStrataClusteringConfig`'s `clusterEngineMap.put(parameters.name(), clusterEngine)`,
    /// row by row: a repeated name overwrites the earlier entry in place, as `put` does.
    ///
    /// The list keeps first-insertion order, which is not the `HashMap`'s; [`run`] recomputes the
    /// Java iteration order where it matters.
    pub fn new(configuration: &[StratumParameters]) -> Engines {
        let mut entries: Vec<(String, StratumParameters)> = Vec::new();
        for parameters in configuration {
            match entries
                .iter_mut()
                .find(|(name, _)| *name == parameters.name)
            {
                Some(entry) => entry.1 = parameters.clone(),
                None => entries.push((parameters.name.clone(), parameters.clone())),
            }
        }
        Engines { entries }
    }

    /// `clusterEngineMap.get(name)`, `None` where the Java returns null.
    pub fn get(&self, name: &str) -> Option<&StratumParameters> {
        self.entries
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, parameters)| parameters)
    }
}

/// `onTraversalStart`, in its own order: the strata must exist, the two configurations must hold
/// the same number of groups, and then every stratum must be named in the clustering one.
///
/// ```java
/// Utils.validate(stratEngine.getStrata().size() == clusterEngineMap.size(),
///         "Stratification and clustering configurations have a different number of groups.");
/// ```
///
/// `clusterEngineMap.size()` counts distinct names, which is why a duplicated row can fail the
/// count here (see [`Engines`]).
pub fn validate(engine: &Engine, engines: &Engines) -> Result<(), GroupedError> {
    if engine.strata.is_empty() {
        return Err(GroupedError::NoStrata);
    }
    if engine.strata.len() != engines.entries.len() {
        return Err(GroupedError::GroupCountMismatch);
    }
    for stratum in &engine.strata {
        if engines.get(&stratum.name).is_none() {
            return Err(GroupedError::GroupNotFound {
                name: stratum.name.clone(),
            });
        }
    }
    Ok(())
}

/// One cluster that was written out, as the stratum it belongs to and the ids it holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cluster {
    pub stratum: String,
    pub members: Vec<String>,
}

/// `applyRecord` over the whole input, then `onTraversalSuccess`.
///
/// The returned order is the strata in the Java `HashMap` iteration order of their names (the
/// order `clusterEngineMap.values()` flushes in), each stratum's clusters in the order
/// `sv_cluster::cluster` returns them, with the unmatched records last. That is NOT the order the
/// tool writes them in: it writes each unmatched record as it is read, emits clusters through
/// `addAndFlush` as they close and flushes the engines afterwards, and it disables its own output
/// index because of it.
///
/// ```java
/// if (stratifications.size() > 1) {
/// ```
///
/// A record matching two strata throws on that record, so records before it have already been
/// grouped; here the error is returned and nothing is.
pub fn run(
    records: &[CallRecord],
    engine: &Engine,
    engines: &Engines,
    thresholds: Thresholds,
    algorithm: Algorithm,
    cluster_del_with_dup: bool,
) -> Result<Vec<Cluster>, GroupedError> {
    validate(engine, engines)?;

    let mut grouped: Vec<(String, Vec<CallRecord>)> = Vec::new();
    let mut unmatched: Vec<Cluster> = Vec::new();
    for record in records {
        // The stratification engine reads a leaner record than the linkage does. `.clone()` copies
        // a string, since the new record must own its own; Java would share the reference.
        let stratified = crate::sv_stratify::CallRecord {
            id: record.id.clone(),
            sv_type: record.sv_type,
            contig_a: record.contig_a.clone(),
            position_a: record.position_a,
            contig_b: record.contig_b.clone(),
            position_b: record.position_b,
            length: record.length,
        };
        // `getMatches`, in the engine's order, which is the configuration's since 4.7.0.0.
        let names: Vec<String> = engine
            .matches(&stratified, thresholds)
            .map_err(GroupedError::Stratify)?
            .into_iter()
            .map(|stratum| stratum.name.clone())
            .collect();
        // `match` on the count: `0`, `1` and `_` (anything else). The compiler refuses a `match`
        // that leaves a case uncovered. The three arms are the Java's `isEmpty()`, "exactly one
        // match" and `size() > 1` branches; the unmatched record keeps
        // `SVStratify.DEFAULT_STRATUM`.
        match names.len() {
            0 => unmatched.push(Cluster {
                stratum: DEFAULT_STRATUM.to_string(),
                members: vec![record.id.clone()],
            }),
            1 => {
                let name = &names[0];
                match grouped.iter_mut().find(|(key, _)| key == name) {
                    Some(entry) => entry.1.push(record.clone()),
                    None => grouped.push((name.clone(), vec![record.clone()])),
                }
            }
            _ => {
                return Err(GroupedError::MultipleMatches {
                    id: record.id.clone(),
                    names,
                })
            }
        }
    }

    // `onTraversalSuccess` flushes `clusterEngineMap.values()`, a HashMap over the names, whatever
    // order the stratification engine keeps.
    let entries: Vec<(&crate::sv_stratify::Stratum, i32)> = engine
        .strata
        .iter()
        .map(|stratum| {
            (
                stratum,
                gatk_engine::java_hash::string_hash_code(&stratum.name),
            )
        })
        .collect();
    // `expect` aborts with that message where the call failed; it marks a failure the code
    // treats as impossible rather than one it reports.
    let flushed = gatk_engine::java_hash::hash_map_order(&entries).expect("a hashable stratum set");
    let mut out = Vec::new();
    for stratum in flushed {
        // `let Some(..) = .. else { continue }`: take the stratum's records, or skip a stratum
        // that no record matched (its engine flushes nothing).
        let Some((_, members)) = grouped.iter().find(|(key, _)| *key == stratum.name) else {
            continue;
        };
        let parameters = engines.get(&stratum.name).expect("a validated group");
        let linkage = linkage_for(parameters, cluster_del_with_dup);
        for members in cluster(members, &linkage, algorithm) {
            out.push(Cluster {
                stratum: stratum.name.clone(),
                members,
            });
        }
    }
    out.extend(unmatched);
    Ok(out)
}
