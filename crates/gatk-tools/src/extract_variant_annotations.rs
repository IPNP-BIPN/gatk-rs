//! `ExtractVariantAnnotations`: which variants a scalable-VQSR extraction keeps, and what it
//! writes for them.
//!
//! One prefix produces up to three files and none of them is a subset of another: a labelled
//! matrix, an unlabelled one, and a sites-only VCF. What is ported is the decision behind each
//! row: the filter and variant-type checks, what counts as the same variant as a resource record,
//! how an annotation becomes a number, and Algorithm R over the unlabelled reservoir.
//!
//! Reading a VCF and writing HDF5 are not ported. Neither is the random stream: the reservoir
//! takes the index it is given, because `java.util.Random` may not be transcribed here.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.ExtractVariantAnnotations`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.LabeledVariantAnnotationsWalker`,
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.VariantType` and
//! `org.broadinstitute.hellbender.tools.walkers.vqsr.scalable.data.LabeledVariantAnnotationsDatum`
//! in GATK 4.6.2.0.

use std::collections::{BTreeMap, BTreeSet};

// ================================================================================================
// The record and its type.
// ================================================================================================

/// One VCF record, reduced to what the extraction reads off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub reference: String,
    pub alternates: Vec<String>,
    /// The FILTER column, split; empty means the record passed.
    pub filters: Vec<String>,
    pub attributes: Vec<(String, String)>,
}

impl Record {
    pub fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// A record with no filters at all, which is what `isNotFiltered` asks about.
    pub fn is_not_filtered(&self) -> bool {
        self.filters.is_empty()
    }

    /// The end, which is the start plus the reference allele's length less one.
    pub fn end(&self) -> i32 {
        self.start + self.reference.len() as i32 - 1
    }
}

/// `VariantContext.Type`, as far as the extraction needs it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContextType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

/// An allele written in angle brackets, or a breakend.
pub fn is_symbolic(allele: &str) -> bool {
    allele.starts_with('<') || allele.contains('[') || allele.contains(']') || allele == "*"
}

/// `typeOfBiallelicVariant`, which is the whole of the type determination for one alternate.
fn biallelic_type(reference: &str, alternate: &str) -> ContextType {
    if is_symbolic(alternate) {
        return ContextType::Symbolic;
    }
    if reference.len() == alternate.len() {
        if alternate.len() == 1 {
            ContextType::Snp
        } else {
            ContextType::Mnp
        }
    } else {
        ContextType::Indel
    }
}

/// `determineType`: no alternates is no variation, and alternates that disagree are mixed.
pub fn context_type(record: &Record) -> ContextType {
    if record.alternates.is_empty() {
        return ContextType::NoVariation;
    }
    let mut kind: Option<ContextType> = None;
    for alternate in &record.alternates {
        let this = biallelic_type(&record.reference, alternate);
        match kind {
            None => kind = Some(this),
            Some(seen) if seen != this => return ContextType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(ContextType::NoVariation)
}

/// The two classes an extraction sorts variants into.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum VariantType {
    Snp,
    Indel,
}

impl VariantType {
    pub fn name(self) -> &'static str {
        match self {
            VariantType::Snp => "SNP",
            VariantType::Indel => "INDEL",
        }
    }
}

/// `getVariantType`, which puts an MNP with the SNPs and everything else with the indels.
///
/// A record with no alternate at all has no type and is refused rather than dropped.
pub fn variant_type(record: &Record) -> Option<VariantType> {
    match context_type(record) {
        ContextType::Snp | ContextType::Mnp => Some(VariantType::Snp),
        ContextType::Indel | ContextType::Mixed | ContextType::Symbolic => Some(VariantType::Indel),
        ContextType::NoVariation => None,
    }
}

/// `getAlleleSpecificVariantType`, which asks only whether the two alleles are the same length.
///
/// A spanning deletion would be a SNP by this rule, which is why it is filtered out upstream
/// rather than here.
pub fn allele_specific_variant_type(reference: &str, alternate: &str) -> VariantType {
    if reference.len() == alternate.len() {
        VariantType::Snp
    } else {
        VariantType::Indel
    }
}

/// `checkVariantType`, which maps the RESOURCE's context type onto one of the two classes and
/// asks the input for the same one.
///
/// An MNP resource record asks for a SNP, and a symbolic one asks for an indel.
pub fn check_variant_type(record: &Record, resource: &Record) -> bool {
    match context_type(resource) {
        ContextType::Snp | ContextType::Mnp => variant_type(record) == Some(VariantType::Snp),
        ContextType::Indel | ContextType::Mixed | ContextType::Symbolic => {
            variant_type(record) == Some(VariantType::Indel)
        }
        ContextType::NoVariation => false,
    }
}

// ================================================================================================
// The annotations.
// ================================================================================================

/// `GATKVCFConstants.ALLELE_SPECIFIC_PREFIX`, which is what makes an annotation name parsed per
/// alternate in allele-specific mode.
pub const ALLELE_SPECIFIC_PREFIX: &str = "AS_";

/// `decodeAnnotation`: the value one annotation takes on one record, as a double.
///
/// Every way of not having a number collapses to NaN, and they are not told apart: an annotation
/// the record does not carry, one that does not parse, and one that parses to an infinity are all
/// the same absence in the matrix.
pub fn decode_annotation(
    record: &Record,
    alternate: Option<&str>,
    name: &str,
    allele_specific: bool,
) -> f64 {
    let value = if allele_specific && name.starts_with(ALLELE_SPECIFIC_PREFIX) {
        let Some(text) = record.attribute(name) else {
            return f64::NAN;
        };
        let values: Vec<&str> = text.split(',').collect();
        if values.is_empty() || text.is_empty() {
            return f64::NAN;
        }
        // The index is into the ALTERNATES, which is the allele index less one.
        let alternate = alternate.expect("allele-specific mode is one alternate at a time");
        let Some(index) = record.alternates.iter().position(|a| a == alternate) else {
            // The allele is not one of this record's, which is a refusal rather than a NaN.
            return f64::NAN;
        };
        values.get(index).and_then(|v| parse_double(v))
    } else {
        record.attribute(name).and_then(parse_double)
    };
    match value {
        // An infinity is turned into a NaN, so the matrix carries no infinities at all.
        Some(value) if value.is_infinite() => f64::NAN,
        Some(value) => value,
        None => f64::NAN,
    }
}

/// `getAttributeAsDouble`, whose failures all become the default rather than a refusal.
fn parse_double(text: &str) -> Option<f64> {
    // `Double.parseDouble` accepts these three spellings, which Rust's parser also accepts.
    text.trim().parse::<f64>().ok()
}

// ================================================================================================
// Matching a resource record.
// ================================================================================================

/// `--resource-matching-strategy`, from the loosest to the strictest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchingStrategy {
    /// The DEFAULT: the start position and the variant class, and nothing about the alleles.
    StartPosition,
    /// One alternate in common, compared as written.
    StartPositionAndGivenRepresentation,
    /// One alternate in common after both are reduced to a minimal representation.
    StartPositionAndMinimalRepresentation,
}

/// The strategy a run uses when none is asked for.
pub const DEFAULT_MATCHING_STRATEGY: MatchingStrategy = MatchingStrategy::StartPosition;

/// `isAlleleInList`: whether `reference`/`alternate` names the same event as one of the
/// resource's alternates, once both pairs are trimmed.
///
/// The two records already share a start, so reconciling means removing the bases they pad with:
/// a common suffix first, then a common prefix.
pub fn is_allele_in_list(
    reference: &str,
    alternate: &str,
    resource_reference: &str,
    resource_alternates: &[String],
) -> bool {
    let (reference, alternate) = trim(reference, alternate);
    resource_alternates.iter().any(|resource_alternate| {
        let (resource_reference, resource_alternate) = trim(resource_reference, resource_alternate);
        reference == resource_reference && alternate == resource_alternate
    })
}

/// The minimal representation of one reference and one alternate: the shared suffix goes first,
/// and then the shared prefix, and at least one base is always left on each side.
fn trim<'a>(reference: &'a str, alternate: &'a str) -> (String, String) {
    let mut reference: Vec<u8> = reference.as_bytes().to_vec();
    let mut alternate: Vec<u8> = alternate.as_bytes().to_vec();
    while reference.len() > 1 && alternate.len() > 1 && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }
    let mut front = 0;
    while reference.len() - front > 1
        && alternate.len() - front > 1
        && reference[front] == alternate[front]
    {
        front += 1;
    }
    (
        String::from_utf8(reference[front..].to_vec()).expect("ascii"),
        String::from_utf8(alternate[front..].to_vec()).expect("ascii"),
    )
}

/// `isMatchingVariant`: whether one resource record labels one input record.
///
/// The three preconditions run before the strategy does, and a resource record that is filtered,
/// carries no alternate, or is of the other variant class never labels anything whatever the
/// strategy says.
pub fn is_matching_variant(
    record: &Record,
    resource: &Record,
    alternate: Option<&str>,
    trust_all_polymorphic: bool,
    resource_is_polymorphic: bool,
    resource_has_genotypes: bool,
    strategy: MatchingStrategy,
) -> bool {
    if !resource.is_not_filtered()
        || resource.alternates.is_empty()
        || !check_variant_type(record, resource)
    {
        return false;
    }
    if !(trust_all_polymorphic || !resource_has_genotypes || resource_is_polymorphic) {
        return false;
    }
    match strategy {
        MatchingStrategy::StartPosition => true,
        MatchingStrategy::StartPositionAndGivenRepresentation => record
            .alternates
            .iter()
            .any(|a| resource.alternates.contains(a)),
        MatchingStrategy::StartPositionAndMinimalRepresentation => match alternate {
            Some(alternate) => is_allele_in_list(
                &record.reference,
                alternate,
                &resource.reference,
                &resource.alternates,
            ),
            None => record.alternates.iter().any(|a| {
                is_allele_in_list(
                    &record.reference,
                    a,
                    &resource.reference,
                    &resource.alternates,
                )
            }),
        },
    }
}

// ================================================================================================
// The extraction.
// ================================================================================================

/// One resource file and the labels its tag carries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Resource {
    pub name: String,
    /// The tag's key-and-value pairs, of which only those whose value is exactly `true` label.
    pub tags: Vec<(String, String)>,
    pub records: Vec<Record>,
    pub has_genotypes: bool,
    pub is_polymorphic: bool,
}

impl Resource {
    /// The labels this resource applies, sorted, which is what a `TreeSet` gives.
    ///
    /// The value has to be the STRING `true`: a tag written `training=false` labels nothing, and
    /// so a run whose only resource is tagged that way extracts nothing at all.
    pub fn labels(&self) -> BTreeSet<String> {
        self.tags
            .iter()
            .filter(|(_, value)| value == "true")
            .map(|(key, _)| key.clone())
            .collect()
    }
}

/// The label the matrix carries for the variant class, which is why it may not come from a
/// resource: a run's own `snp` column would collide with it.
pub const SNP_LABEL: &str = "snp";
pub const TRAINING_LABEL: &str = "training";
pub const CALIBRATION_LABEL: &str = "calibration";

/// The refusal a resource tagged with the reserved label produces.
pub fn check_resource_labels(labels: &BTreeSet<String>) -> Result<(), String> {
    if labels.contains(SNP_LABEL) {
        return Err(format!(
            "Bad input: The resource label \"{SNP_LABEL}\" is reserved for labeling variant types."
        ));
    }
    Ok(())
}

/// What one input record contributes: one entry per row it becomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extracted {
    /// In allele-specific mode a single alternate; otherwise every alternate that passed.
    pub alternates: Vec<String>,
    pub variant_type: VariantType,
    /// Sorted, and empty for an unlabelled row.
    pub labels: BTreeSet<String>,
}

/// The arguments the extraction reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    /// `--mode`, which may name both.
    pub modes: BTreeSet<VariantType>,
    /// `--ignore-filter`, each naming one filter.
    pub ignored_filters: BTreeSet<String>,
    /// `--ignore-all-filters`.
    pub ignore_all_filters: bool,
    /// The opposite of `--do-not-trust-all-polymorphic`.
    pub trust_all_polymorphic: bool,
    pub strategy: MatchingStrategy,
    /// Whether any requested annotation is declared `Number=A`, which switches the WHOLE run.
    pub allele_specific: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            modes: [VariantType::Snp, VariantType::Indel].into_iter().collect(),
            ignored_filters: BTreeSet::new(),
            ignore_all_filters: false,
            trust_all_polymorphic: true,
            strategy: DEFAULT_MATCHING_STRATEGY,
            allele_specific: false,
        }
    }
}

/// Whether a record's filters let it through.
///
/// `--ignore-filter` has to name EVERY filter the record carries, not just one of them.
pub fn passes_filters(record: &Record, arguments: &Arguments) -> bool {
    arguments.ignore_all_filters
        || record.is_not_filtered()
        || record
            .filters
            .iter()
            .all(|filter| arguments.ignored_filters.contains(filter))
}

/// `extractVariantMetadata`: the rows one input record becomes, which may be none.
///
/// `extract_unlabeled` is on exactly when a reservoir was asked for, and it is what decides
/// whether a row with no label survives this function at all.
pub fn extract_variant_metadata(
    record: &Record,
    resources: &[Resource],
    arguments: &Arguments,
    extract_unlabeled: bool,
) -> Vec<Extracted> {
    if !passes_filters(record, arguments) {
        return Vec::new();
    }
    if !arguments.allele_specific {
        let Some(kind) = variant_type(record) else {
            return Vec::new();
        };
        if !arguments.modes.contains(&kind) {
            return Vec::new();
        }
        let labels = matching_labels(record, None, resources, arguments);
        if extract_unlabeled || !labels.is_empty() {
            return vec![Extracted {
                alternates: record.alternates.clone(),
                variant_type: kind,
                labels,
            }];
        }
        return Vec::new();
    }
    record
        .alternates
        .iter()
        // A spanning deletion would be called a SNP by the allele-specific rule, so it is dropped
        // here instead.
        .filter(|alternate| *alternate != "*")
        .filter(|alternate| {
            arguments
                .modes
                .contains(&allele_specific_variant_type(&record.reference, alternate))
        })
        .map(|alternate| Extracted {
            alternates: vec![alternate.clone()],
            variant_type: allele_specific_variant_type(&record.reference, alternate),
            labels: matching_labels(record, Some(alternate), resources, arguments),
        })
        .filter(|extracted| extract_unlabeled || !extracted.labels.is_empty())
        .collect()
}

/// `findMatchingResourceLabels`: every label of every resource with a matching record at the
/// input's start.
fn matching_labels(
    record: &Record,
    alternate: Option<&str>,
    resources: &[Resource],
    arguments: &Arguments,
) -> BTreeSet<String> {
    let mut labels = BTreeSet::new();
    for resource in resources {
        // The resource is queried at the START position alone, so a record that merely overlaps
        // is not consulted.
        let matches = resource
            .records
            .iter()
            .filter(|candidate| {
                candidate.contig == record.contig && candidate.start == record.start
            })
            .any(|candidate| {
                is_matching_variant(
                    record,
                    candidate,
                    alternate,
                    arguments.trust_all_polymorphic,
                    resource.is_polymorphic,
                    resource.has_genotypes,
                    // Allele-specific mode forces the strictest strategy, whatever was asked for.
                    if arguments.allele_specific {
                        MatchingStrategy::StartPositionAndMinimalRepresentation
                    } else {
                        arguments.strategy
                    },
                )
            });
        if matches {
            labels.extend(resource.labels());
        }
    }
    labels
}

// ================================================================================================
// The two matrices.
// ================================================================================================

/// One row of a matrix: where it came from, what it is, what labelled it, and its numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference: String,
    pub alternates: Vec<String>,
    pub variant_type: VariantType,
    pub labels: BTreeSet<String>,
    pub annotations: Vec<f64>,
}

/// The annotation columns, SORTED BY NAME whatever order they were asked for in, and deduplicated.
pub fn sorted_annotation_names(requested: &[String]) -> Vec<String> {
    let sorted: BTreeSet<String> = requested.iter().cloned().collect();
    sorted.into_iter().collect()
}

/// The label columns, sorted, with the reserved one always first because it sorts before the
/// others a resource can supply.
pub fn sorted_labels(resource_labels: &BTreeSet<String>) -> Vec<String> {
    let mut labels: BTreeSet<String> = resource_labels.clone();
    labels.insert(SNP_LABEL.to_string());
    labels.into_iter().collect()
}

/// One row built from one extracted entry.
pub fn row(record: &Record, extracted: &Extracted, names: &[String], allele_specific: bool) -> Row {
    let alternate = if allele_specific {
        Some(extracted.alternates[0].as_str())
    } else {
        None
    };
    Row {
        contig: record.contig.clone(),
        start: record.start,
        end: record.end(),
        reference: record.reference.clone(),
        alternates: extracted.alternates.clone(),
        variant_type: extracted.variant_type,
        labels: extracted.labels.clone(),
        annotations: names
            .iter()
            .map(|name| decode_annotation(record, alternate, name, allele_specific))
            .collect(),
    }
}

/// Algorithm R, with the random index handed in.
///
/// The reservoir fills in order and is then overwritten in place, which is why it is NOT in
/// genomic order: a record that arrives late lands wherever the index says.
#[derive(Debug, Clone, PartialEq)]
pub struct Reservoir {
    pub maximum: usize,
    pub rows: Vec<Row>,
    /// How many unlabelled records have been offered, which is what the index is drawn against.
    pub seen: usize,
}

impl Reservoir {
    pub fn new(maximum: usize) -> Reservoir {
        Reservoir {
            maximum,
            rows: Vec::new(),
            seen: 0,
        }
    }

    /// One unlabelled record, and the index the random stream would have produced for it.
    ///
    /// `index` is `rng.nextInt(seen)` and is consulted ONLY once the reservoir is full, so a
    /// caller that has not filled it yet may pass anything.
    pub fn offer(&mut self, rows: Vec<Row>, index: usize) {
        if self.seen < self.maximum {
            self.rows.extend(rows);
        } else if index < self.maximum {
            // `set` replaces the whole record at that slot, not one row of it.
            self.rows[index] = rows.into_iter().next().expect("at least one row");
        }
        self.seen += 1;
    }
}

/// A whole extraction: the labelled rows in the order they were seen, and the reservoir beside
/// them.
#[derive(Debug, Clone, PartialEq)]
pub struct Extraction {
    pub names: Vec<String>,
    pub labeled: Vec<Row>,
    pub unlabeled: Option<Reservoir>,
}

impl Extraction {
    /// Whether a matrix file is written at all.
    ///
    /// An extraction that kept nothing writes NO labelled matrix, which is a warning rather than
    /// a refusal, and still writes its sites-only VCF.
    pub fn writes_labeled_matrix(&self) -> bool {
        !self.labeled.is_empty()
    }

    /// The unlabelled matrix is written only when a reservoir was asked for at all.
    pub fn writes_unlabeled_matrix(&self) -> bool {
        self.unlabeled.as_ref().is_some_and(|r| !r.rows.is_empty())
    }
}

/// One pass over the input: every record, split into the labelled rows and the reservoir.
///
/// `random` is `rng.nextInt(bound)`, called once per unlabelled record after the reservoir fills.
pub fn extract(
    records: &[Record],
    resources: &[Resource],
    arguments: &Arguments,
    requested: &[String],
    maximum_unlabeled: usize,
    random: &mut dyn FnMut(usize) -> usize,
) -> Extraction {
    let names = sorted_annotation_names(requested);
    let mut labeled = Vec::new();
    let mut reservoir = (maximum_unlabeled > 0).then(|| Reservoir::new(maximum_unlabeled));
    for record in records {
        let metadata = extract_variant_metadata(record, resources, arguments, reservoir.is_some());
        if metadata.is_empty() {
            continue;
        }
        for extracted in metadata.iter().filter(|e| !e.labels.is_empty()) {
            labeled.push(row(record, extracted, &names, arguments.allele_specific));
        }
        let Some(reservoir) = reservoir.as_mut() else {
            continue;
        };
        let unlabeled: Vec<Row> = metadata
            .iter()
            .filter(|e| e.labels.is_empty())
            .map(|extracted| row(record, extracted, &names, arguments.allele_specific))
            .collect();
        if unlabeled.is_empty() {
            continue;
        }
        // The counter advances per RECORD that had an unlabelled row, not per row.
        let index = if reservoir.seen < reservoir.maximum {
            0
        } else {
            random(reservoir.seen)
        };
        reservoir.offer(unlabeled, index);
    }
    Extraction {
        names,
        labeled,
        unlabeled: reservoir,
    }
}

/// The reference allele of every row, in order, which is one of the two allele datasets.
pub fn reference_alleles(rows: &[Row]) -> Vec<String> {
    rows.iter().map(|row| row.reference.clone()).collect()
}

/// The alternate alleles of every row, FLATTENED.
///
/// A multiallelic record contributes both of its alternates, so this array is longer than the
/// reference one beside it and the two cannot be read side by side.
pub fn alternate_alleles(rows: &[Row]) -> Vec<String> {
    rows.iter()
        .flat_map(|row| row.alternates.iter().cloned())
        .collect()
}

/// The `snp` column, which is written for every run whatever the resources say.
pub fn snp_column(rows: &[Row]) -> Vec<bool> {
    rows.iter()
        .map(|row| row.variant_type == VariantType::Snp)
        .collect()
}

/// One label's column.
pub fn label_column(rows: &[Row], label: &str) -> Vec<bool> {
    rows.iter().map(|row| row.labels.contains(label)).collect()
}

/// The INFO column of the sites-only VCF, which carries the labels as bare flags and nothing else.
pub fn sites_only_info(labels: &BTreeSet<String>) -> String {
    if labels.is_empty() {
        ".".to_string()
    } else {
        labels.iter().cloned().collect::<Vec<_>>().join(";")
    }
}

/// The attributes of a record, as the dump's INFO column writes them.
pub fn attributes(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
    pairs
        .iter()
        .map(|(key, value)| ((*key).to_string(), (*value).to_string()))
        .collect()
}

/// A convenience for the suffixes the two matrices are written under.
pub const ANNOTATIONS_HDF5_SUFFIX: &str = ".annot.hdf5";
pub const UNLABELED_TAG: &str = ".unlabeled";

/// The paths the matrix's datasets live at, which is what makes `--omit-alleles-in-hdf5`
/// observable: those two go and the rest stay.
pub const ALLELES_REF_PATH: &str = "/alleles/ref";
pub const ALLELES_ALT_PATH: &str = "/alleles/alt";
pub const INTERVALS_PATH: &str = "/intervals";

/// Which datasets a matrix carries.
pub fn datasets(omit_alleles: bool) -> BTreeMap<&'static str, bool> {
    [
        (INTERVALS_PATH, true),
        (ALLELES_REF_PATH, !omit_alleles),
        (ALLELES_ALT_PATH, !omit_alleles),
    ]
    .into_iter()
    .collect()
}
