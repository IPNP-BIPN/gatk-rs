//! `VariantRecalibrator`: a callset's annotations modelled by two Gaussian mixtures, every variant
//! scored by their log-odds, and the tranches a truth sensitivity target cuts out of the scores.
//!
//! The model is [`gatk_tools::variant_recalibrator_model`] and the tranches
//! [`gatk_tools::variant_recalibrator`]. What is here is the `MultiVariantWalker` around them, in
//! the order the reference touches things, because that order is also the order of the draws
//! from the one random stream the whole run shares:
//!
//! * **`onTraversalStart`**: the Rscript check when `--rscript-file` is given, the tranches file
//!   opened, the resources read into training sets and refused when none is a training or a truth
//!   set, the serialised model read and its annotations matched against the command line, the
//!   recal writer's header, and then four hundred `nextDouble` draws for a replicate list nothing
//!   reads again;
//! * **the traversal**: the inputs merged in `MergingIterator`'s order, variants QUEUED by start
//!   and turned into data only when the start moves on, so an exception while decoding one is
//!   wrapped by the walker naming the record whose arrival flushed the queue. Each datum's
//!   annotations are decoded (and jittered) as it is built, and its labels come from the resource
//!   records that START where the variant starts;
//! * **`onTraversalSuccess`**, up to `--max-attempts` times over the SAME data, since the reference
//!   retries without resetting anything: the normalisation, the positive model, the negative model
//!   over the worst scores, the contrastive scores, the model report, the culprits, the tranches
//!   (which sort the data by score in place) and the recal file (which sorts them back by
//!   coordinate, stably), and last the plotting script.
use super::*;
use gatk_tools::variant_recalibrator as tranches;
use gatk_tools::variant_recalibrator_model as model;
use htsjdk_vcf::variant::VariantContext;

/// A resource, its tags, and its records.
struct Resource {
    set: model::TrainingSet,
    records: Vec<VariantContext>,
}

impl Resource {
    /// `featureContext.getValues(source, start)`: the records overlapping the variant that start
    /// where it starts, in file order.
    fn at<'a>(&'a self, contig: &str, start: i64, stop: i64) -> Vec<&'a VariantContext> {
        records_at(&self.records, contig, start, stop)
    }
}

fn records_at<'a>(
    records: &'a [VariantContext],
    contig: &str,
    start: i64,
    stop: i64,
) -> Vec<&'a VariantContext> {
    records
        .iter()
        .filter(|record| {
            record.contig == contig
                && record.start == start
                && record.start <= stop
                && record.stop >= start
        })
        .collect()
}

fn thrown(error: model::ModelError) -> Thrown {
    Thrown {
        failure: if error.is_user() {
            Failure::User
        } else {
            Failure::Other
        },
        exception: error.class,
        message: Some(error.message),
    }
}

/// `VariantContext.getType`, in the classes `checkVariationClass` asks about.
fn is_snp_or_mnp(record: &VariantContext) -> bool {
    use gatk_tools::remove_nearby_indels::{variant_type, VariantType};
    matches!(variant_type(record), VariantType::Snp | VariantType::Mnp)
}

fn is_indel_like(record: &VariantContext) -> bool {
    use gatk_tools::remove_nearby_indels::{variant_type, VariantType};
    matches!(
        variant_type(record),
        VariantType::Indel | VariantType::Mixed | VariantType::Symbolic
    )
}

/// `checkVariationClass(evalVC, mode)`.
fn of_mode(record: &VariantContext, mode: tranches::Mode) -> bool {
    match mode {
        tranches::Mode::Snp => is_snp_or_mnp(record),
        tranches::Mode::Indel => is_indel_like(record),
        tranches::Mode::Both => true,
    }
}

/// `checkVariationClass(evalVC, allele, mode)`, where a spanning deletion counts as a SNP.
fn allele_of_mode(
    record: &VariantContext,
    allele: &htsjdk_vcf::allele::Allele,
    mode: tranches::Mode,
) -> bool {
    match mode {
        tranches::Mode::Snp => record.reference().len() == allele.len(),
        tranches::Mode::Indel => record.reference().len() != allele.len() || allele.is_symbolic(),
        tranches::Mode::Both => true,
    }
}

/// `checkVariationClass(evalVC, trainVC)`: the training record's type decides the class asked.
fn same_class(eval: &VariantContext, train: &VariantContext) -> bool {
    use gatk_tools::remove_nearby_indels::{variant_type, VariantType};
    match variant_type(train) {
        VariantType::Snp | VariantType::Mnp => of_mode(eval, tranches::Mode::Snp),
        VariantType::Indel | VariantType::Mixed | VariantType::Symbolic => {
            of_mode(eval, tranches::Mode::Indel)
        }
        VariantType::NoVariation => false,
    }
}

/// `isPolymorphicInSamples`: some called allele is not the reference.
fn is_polymorphic_in_samples(record: &VariantContext) -> bool {
    let (mut called, mut reference) = (0, 0);
    for genotype in record.genotypes.iter() {
        for allele in &genotype.alleles {
            if allele.is_no_call() {
                continue;
            }
            called += 1;
            if allele.is_reference() {
                reference += 1;
            }
        }
    }
    record.alleles.len() > 1 && called != reference
}

/// `isValidVariant`.
fn is_valid_training_variant(eval: &VariantContext, train: &VariantContext, trust: bool) -> bool {
    !train.is_filtered()
        && train.alleles.len() > 1
        && same_class(eval, train)
        && (trust || train.genotypes.iter().next().is_none() || is_polymorphic_in_samples(train))
}

/// `Allele.getBaseString`, which a reference allele carries without its star.
fn bases(allele: &htsjdk_vcf::allele::Allele) -> String {
    allele.base_string()
}

/// `GATKVariantContextUtils.isTransition` on a biallelic SNP.
fn is_transition(record: &VariantContext) -> bool {
    let reference = bases(record.reference()).to_ascii_uppercase();
    let alternate = record
        .alternate_alleles()
        .first()
        .map(|allele| bases(allele).to_ascii_uppercase())
        .unwrap_or_default();
    matches!(
        (reference.as_bytes().first(), alternate.as_bytes().first()),
        (Some(b'A'), Some(b'G'))
            | (Some(b'G'), Some(b'A'))
            | (Some(b'C'), Some(b'T'))
            | (Some(b'T'), Some(b'C'))
    )
}

/// What an INFO key holds on a record, as the decoder reads it.
fn attribute<'a>(record: &'a VariantContext, key: &str) -> model::AttributeValue<'a> {
    use htsjdk_vcf::variant::Value;
    let Some((_, value)) = record.attributes.iter().find(|(name, _)| name == key) else {
        return model::AttributeValue::Absent;
    };
    match value {
        Value::Str(text) => {
            if text.contains(',') {
                model::AttributeValue::Many(text.split(',').collect())
            } else {
                model::AttributeValue::One(text)
            }
        }
        Value::List(values) => model::AttributeValue::Many(
            values
                .iter()
                .map(|value| match value {
                    Value::Str(text) => text.as_str(),
                    _ => ".",
                })
                .collect(),
        ),
        Value::Missing => model::AttributeValue::One("."),
        Value::Bool(true) => model::AttributeValue::One("true"),
        Value::Bool(false) => model::AttributeValue::Absent,
        // A decoded record keeps its INFO values as text; the other forms never reach here.
        _ => model::AttributeValue::One("."),
    }
}

/// Everything `addVariantDatum` needs that is not the record.
struct DatumContext<'a> {
    keys: &'a [String],
    arguments: &'a model::ModelArguments,
    mode: tranches::Mode,
    resources: &'a [Resource],
    trust_all_polymorphic: bool,
    ignore_all_filters: bool,
    ignored_filters: &'a [String],
}

/// `addVariantDatum`: nothing for a record the filters or the mode exclude, one datum for a site,
/// and one per alternate in allele-specific mode.
///
/// `site` is the input variant the feature context belongs to, which is the record itself unless
/// this is an aggregate one.
fn add_variant_datum(
    data: &mut Vec<model::Datum>,
    record: &VariantContext,
    site: &VariantContext,
    is_input: bool,
    context: &DatumContext,
    random: &mut gatk_engine::java_random::JavaRandom,
) -> Result<(), Thrown> {
    let filters = record.filters.clone().unwrap_or_default();
    let passes = context.ignore_all_filters
        || filters.is_empty()
        || filters
            .iter()
            .all(|filter| context.ignored_filters.contains(filter));
    if !passes {
        return Ok(());
    }
    let use_as = context.arguments.use_as_annotations;
    if !use_as {
        if of_mode(record, context.mode) {
            add_datum(data, record, site, is_input, None, context, random)?;
        }
        return Ok(());
    }
    for (index, allele) in record.alternate_alleles().iter().enumerate() {
        if allele.display_string() == "*" || !allele_of_mode(record, allele, context.mode) {
            continue;
        }
        add_datum(data, record, site, is_input, Some(index), context, random)?;
    }
    Ok(())
}

/// `addDatum`: the annotations decoded, the labels read off the resources, the prior as log odds.
fn add_datum(
    data: &mut Vec<model::Datum>,
    record: &VariantContext,
    site: &VariantContext,
    is_input: bool,
    allele: Option<usize>,
    context: &DatumContext,
    random: &mut gatk_engine::java_random::JavaRandom,
) -> Result<(), Thrown> {
    let mut annotations = Vec::with_capacity(context.keys.len());
    let mut is_null = Vec::with_capacity(context.keys.len());
    for key in context.keys {
        let value = model::decode_annotation(
            key,
            attribute(record, key),
            allele,
            context.arguments,
            random,
        )
        .map_err(thrown)?;
        is_null.push(value.is_nan());
        annotations.push(value);
    }
    let mut datum = model::Datum::new(annotations, is_null);
    let reference = bases(record.reference());
    let alternate = allele.map(|index| record.alternate_alleles()[index].display_string());
    if allele.is_some() {
        datum.reference_allele = Some(reference.clone());
        datum.alternate_allele = alternate.clone();
    }
    if is_input {
        datum.loc = Some((record.contig.clone(), record.start, record.stop));
    }
    datum.is_snp = matches!(
        gatk_tools::remove_nearby_indels::variant_type(record),
        gatk_tools::remove_nearby_indels::VariantType::Snp
    ) && record.alleles.len() == 2;
    datum.is_transition = datum.is_snp && is_transition(record);
    datum.is_aggregate = !is_input;

    // `parseTrainingSets`, over the resource records starting where the FEATURE CONTEXT's
    // variant starts.
    datum.prior = 2.0;
    for resource in context.resources {
        for train in resource.at(&site.contig, site.start, site.stop) {
            if let Some(alternate) = &alternate {
                let others: Vec<String> = train
                    .alternate_alleles()
                    .iter()
                    .map(|allele| allele.display_string())
                    .collect();
                match is_allele_in_list(&reference, alternate, &bases(train.reference()), &others) {
                    Ok(true) => {}
                    Ok(false) => continue,
                    Err(()) => {
                        return Err(Thrown::non_user(
                            "java.lang.IllegalStateException",
                            format!(
                                "Reference allele mismatch at position {}:{} : ",
                                train.contig, train.start
                            ),
                        ))
                    }
                }
            }
            if is_valid_training_variant(record, train, context.trust_all_polymorphic) {
                datum.is_known = datum.is_known || resource.set.is_known;
                datum.at_truth_site = datum.at_truth_site || resource.set.is_truth;
                datum.at_training_site = datum.at_training_site || resource.set.is_training;
                datum.prior = model::java_max(datum.prior, resource.set.prior);
            }
            datum.at_anti_training_site =
                datum.at_anti_training_site || resource.set.is_anti_training;
        }
    }
    datum.prior = model::prior_log_odds(datum.prior);
    data.push(datum);
    Ok(())
}

/// The input's text read as a VCF, the way the port reads every driving file.
fn read_feature_vcf(path: &str) -> Result<htsjdk_vcf::reader::VcfFile, Thrown> {
    let (_, text) = open_feature_input(path)?;
    htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })
}

/// The run's arguments, read once.
struct Settings {
    mode: tranches::Mode,
    arguments: model::ModelArguments,
    max_attempts: i32,
    targets: Vec<f64>,
    scatter: bool,
    vqslod_tranches: Vec<f64>,
}

fn double_argument(parser: &Parser, name: &str, default: f64) -> f64 {
    scalar(parser, name)
        .and_then(|value| model::java_parse_double(&value))
        .unwrap_or(default)
}

fn int_argument(parser: &Parser, name: &str, default: i32) -> i32 {
    scalar(parser, name)
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

/// `VQSLOD_TRANCHES`' default, built by the reference's own accumulating loops.
fn default_vqslod_tranches() -> Vec<f64> {
    let mut out = Vec::new();
    let mut i = 10.0;
    while i > 5.0 {
        out.push(i);
        i -= 0.1;
    }
    let mut i = 5.0;
    while i > -5.0 {
        out.push(i);
        i -= 0.01;
    }
    let mut i = -5.0;
    while i > -10.0 {
        out.push(i);
        i -= 0.1;
    }
    out
}

fn settings(parser: &Parser) -> Settings {
    let mode = match scalar(parser, "mode").as_deref() {
        Some("INDEL") => tranches::Mode::Indel,
        Some("BOTH") => tranches::Mode::Both,
        _ => tranches::Mode::Snp,
    };
    let list = |name: &str, default: Vec<f64>| {
        let values = arguments(parser, name);
        let from_definition: Vec<f64> = parser
            .definitions()
            .iter()
            .find(|definition| definition.long_name() == name)
            .and_then(|definition| match &definition.value {
                gatk_barclay::Value::List(values) => Some(
                    values
                        .iter()
                        .filter_map(|value| match value {
                            gatk_barclay::Value::Double(number) => Some(*number),
                            other => model::java_parse_double(&other.to_java_string()),
                        })
                        .collect(),
                ),
                _ => None,
            })
            .unwrap_or_default();
        if !from_definition.is_empty() {
            from_definition
        } else if !values.is_empty() {
            values
                .iter()
                .filter_map(|value| model::java_parse_double(value))
                .collect()
        } else {
            default
        }
    };
    Settings {
        mode,
        arguments: model::ModelArguments {
            use_as_annotations: flag(parser, "use-allele-specific-annotations"),
            max_gaussians: int_argument(parser, "max-gaussians", 8),
            max_negative_gaussians: int_argument(parser, "max-negative-gaussians", 2),
            max_iterations: int_argument(parser, "max-iterations", 150),
            k_means_iterations: int_argument(parser, "k-means-iterations", 100),
            std_threshold: double_argument(parser, "standard-deviation-threshold", 10.0),
            shrinkage: double_argument(parser, "shrinkage", 1.0),
            dirichlet: double_argument(parser, "dirichlet", 0.001),
            prior_counts: double_argument(parser, "prior-counts", 20.0),
            max_training: int_argument(parser, "maximum-training-variants", 2_500_000),
            min_bad_variants: int_argument(parser, "minimum-bad-variants", 1000),
            bad_lod_cutoff: double_argument(parser, "bad-lod-score-cutoff", -5.0),
            mq_cap: int_argument(parser, "mq-cap-for-logit-jitter-transform", 0),
            mq_jitter: double_argument(parser, "mq-jitter", 0.05),
        },
        max_attempts: int_argument(parser, "max-attempts", 1),
        targets: list("truth-sensitivity-tranche", vec![100.0, 99.9, 99.0, 90.0]),
        scatter: flag(parser, "output-tranches-for-scatter"),
        vqslod_tranches: list("vqslod-tranche", default_vqslod_tranches()),
    }
}

/// The serialised model's tables, and the order the command line's annotations take from it.
struct LoadedModel {
    tables: Vec<model::ReportTable>,
    order: Vec<usize>,
}

impl LoadedModel {
    fn table(&self, name: &str) -> Result<&model::ReportTable, Thrown> {
        self.tables
            .iter()
            .find(|table| table.name == name)
            .ok_or_else(|| {
                Thrown::non_user(
                    "java.lang.IllegalArgumentException",
                    format!("Table {name} not found in the report"),
                )
            })
    }
}

pub fn variant_recalibrator(parser: &Parser) -> Outcome {
    let _ = resolve_read_filters(parser, "VariantRecalibrator")?;
    let settings = settings(parser);
    let paths = arguments(parser, "variant");
    if paths.is_empty() {
        return Err(Thrown::command_line(
            "Argument variant was missing: Argument 'variant' is required",
        ));
    }
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let tranches_path = argument(parser, "tranches-file").ok_or_else(|| {
        Thrown::command_line(
            "Argument tranches-file was missing: Argument 'tranches-file' is required",
        )
    })?;

    // `initializeDrivingVariants`: a repeated input is refused before any is opened.
    let names: Vec<String> = paths.iter().map(|path| absolute_input_name(path)).collect();
    for (index, name) in names.iter().enumerate() {
        if names[..index].contains(name) {
            return Err(Thrown::user(format!(
                "Bad input: Feature inputs must be unique: {name}"
            )));
        }
    }
    let mut inputs = Vec::new();
    let mut texts = Vec::new();
    for path in &paths {
        let (codec, text) = open_feature_input(path)?;
        let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
            failure: Failure::User,
            exception: "htsjdk.tribble.TribbleException",
            message: Some(failure.error.message()),
        })?;
        texts.push((codec, text));
        inputs.push(file);
    }
    let dictionaries: Vec<SamHeader> = texts.iter().map(|(_, text)| vcf_dictionary(text)).collect();
    let names_of = |header: &SamHeader| -> Vec<(String, i32)> {
        header
            .sequences
            .iter()
            .map(|sequence| (sequence.name.clone(), sequence.length))
            .collect()
    };
    if dictionaries
        .iter()
        .any(|header| names_of(header) != names_of(&dictionaries[0]))
    {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "--variant inputs whose sequence dictionaries differ are merged into one by GATK, \
             which this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let VariantWalkerStart { intervals, .. } = variant_walker_validation(
        parser,
        paths[0].clone(),
        texts[0].1.clone(),
        texts[0].0,
        dictionaries[0].clone(),
    )?;
    let contigs: Vec<(String, i32)> = names_of(&dictionaries[0]);

    // The resources and the aggregate inputs are feature inputs too, opened at startup.
    let resource_values: Vec<(String, Vec<(String, String)>)> = parser
        .definitions()
        .iter()
        .find(|definition| definition.long_name() == "resource")
        .map(|definition| {
            let one = |value: &gatk_barclay::Value| match value {
                gatk_barclay::Value::Tagged {
                    value, attributes, ..
                } => Some((value.clone(), attributes.clone())),
                gatk_barclay::Value::Str(text) => Some((text.clone(), Vec::new())),
                _ => None,
            };
            match &definition.value {
                gatk_barclay::Value::List(values) => values.iter().filter_map(one).collect(),
                other => one(other).into_iter().collect(),
            }
        })
        .unwrap_or_default();
    let mut resource_files = Vec::new();
    for (path, attributes) in &resource_values {
        resource_files.push((path.clone(), attributes.clone(), read_feature_vcf(path)?));
    }
    let mut aggregates = Vec::new();
    for path in arguments(parser, "aggregate") {
        aggregates.push((path.clone(), read_feature_vcf(&path)?));
    }

    // `onTraversalStart`.
    let keys: Vec<String> = arguments(parser, "use-annotation");
    let mut manager = model::DataManager::new(&keys);
    let rscript = argument(parser, "rscript-file");
    let dont_run_rscript = flag(parser, "dont-run-rscript");
    if rscript.is_some() && !dont_run_rscript {
        // The pinned container has no Rscript on its PATH, and neither does the port.
        return Err(Thrown::user(
            "Rscript not found in environment path. Fix executor or run with --dont-run-rscript \
             argument to generate rscript file without running.",
        ));
    }
    let mut ignored: Vec<String> = arguments(parser, "ignore-filter");
    ignored.sort();
    ignored.dedup();
    std::fs::write(&tranches_path, b"").map_err(|_| {
        Thrown::user(format!(
            "Couldn't write file {} because exception {}",
            java_absolute_path(&tranches_path),
            tranches_path
        ))
    })?;
    let mut resources = Vec::new();
    for (_, attributes, file) in resource_files {
        let set = model::TrainingSet::from_attributes(&attributes).map_err(thrown)?;
        resources.push(Resource {
            set,
            records: file.records,
        });
    }
    if !resources.iter().any(|resource| resource.set.is_training) {
        return Err(Thrown::command_line(
            "No training set found! Please provide sets of known polymorphic loci marked with the \
             training=true feature input tag. For example, -resource \
             hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf",
        ));
    }
    if !resources.iter().any(|resource| resource.set.is_truth) {
        return Err(Thrown::command_line(
            "No truth set found! Please provide sets of known polymorphic loci marked with the \
             truth=true feature input tag. For example, -resource \
             hapmap,VCF,known=false,training=true,truth=true,prior=12.0 hapmapFile.vcf",
        ));
    }
    let loaded = match argument(parser, "input-model") {
        None => None,
        Some(path) => {
            let text = std::fs::read_to_string(&path).map_err(|_| {
                Thrown::user(format!(
                    "Couldn't read file. Error was: File: ({path}) with exception: {path}"
                ))
            })?;
            let tables = model::read_report_tables(&text);
            let means = tables
                .iter()
                .find(|table| table.name == "AnnotationMeans")
                .cloned()
                .ok_or_else(|| {
                    Thrown::non_user(
                        "java.lang.IllegalArgumentException",
                        "Table AnnotationMeans not found in the report",
                    )
                })?;
            let mut order = Vec::new();
            for row in &means.rows {
                for (j, key) in manager.annotation_keys.iter().enumerate() {
                    if row.first() == Some(key) {
                        order.push(j);
                    }
                }
            }
            if order.len() != means.rows.len() || order.len() != manager.annotation_keys.len() {
                return Err(Thrown::command_line(
                    "Annotations specified on the command line do not match annotations in the \
                     model report.",
                ));
            }
            let vector = |name: &str| -> Vec<(String, f64)> {
                tables
                    .iter()
                    .find(|table| table.name == name)
                    .map(|table| {
                        table
                            .rows
                            .iter()
                            .map(|row| {
                                (
                                    row[0].clone(),
                                    model::java_parse_double(&row[1]).unwrap_or(f64::NAN),
                                )
                            })
                            .collect()
                    })
                    .unwrap_or_default()
            };
            manager.set_normalization(&vector("AnnotationMeans"), &vector("AnnotationStdevs"));
            Some(LoadedModel { tables, order })
        }
    };

    let mut random = gatk_random();
    for _ in 0..400 {
        random.next_double();
    }

    // The traversal: the inputs merged, each cut to `-L`, each record queued by start.
    let mut reached: Vec<Vec<&VariantContext>> = Vec::new();
    for (input, file) in paths.iter().zip(&inputs) {
        reached.push(variants_in_traversal(
            &file.records,
            intervals.as_deref(),
            input,
        )?);
    }
    let contig_index = |name: &str| {
        contigs
            .iter()
            .position(|(contig, _)| contig == name)
            .map_or(-1, |index| index as i64)
    };
    let keys_of: Vec<Vec<(i64, i64)>> = reached
        .iter()
        .map(|records| {
            records
                .iter()
                .map(|record| (contig_index(&record.contig), record.start))
                .collect()
        })
        .collect();
    let failed = |failure: gatk_tools::combine_gvcfs::Failure| match failure {
        gatk_tools::combine_gvcfs::Failure::User(message) => Thrown::user(message),
        gatk_tools::combine_gvcfs::Failure::Runtime { class, message } => {
            Thrown::non_user(java_class_name(&class), message)
        }
    };
    let order = gatk_tools::combine_gvcfs::merging_order(&keys_of).map_err(failed)?;
    let context = DatumContext {
        keys: &manager.annotation_keys.clone(),
        arguments: &settings.arguments,
        mode: settings.mode,
        resources: &resources,
        trust_all_polymorphic: flag(parser, "trust-all-polymorphic"),
        ignore_all_filters: flag(parser, "ignore-all-filters"),
        ignored_filters: &ignored,
    };
    let mut data: Vec<model::Datum> = Vec::new();
    let mut queue: Vec<&VariantContext> = Vec::new();
    let consume = |queue: &mut Vec<&VariantContext>,
                   data: &mut Vec<model::Datum>,
                   random: &mut gatk_engine::java_random::JavaRandom|
     -> Result<(), Thrown> {
        for record in queue.iter() {
            add_variant_datum(data, record, record, true, &context, random)?;
        }
        if let Some(first) = queue.first() {
            // Every aggregate input's records at the first queued variant's start, input by input.
            for (_, file) in &aggregates {
                for record in records_at(&file.records, &first.contig, first.start, first.stop) {
                    add_variant_datum(data, record, first, false, &context, random)?;
                }
            }
        }
        queue.clear();
        Ok(())
    };
    for (input, index) in order {
        let record = reached[input][index];
        let changed = queue
            .first()
            .is_some_and(|first| first.start != record.start || first.contig != record.contig);
        if changed {
            consume(&mut queue, &mut data, &mut random).map_err(|_| {
                Thrown::non_user(
                    "org.broadinstitute.hellbender.exceptions.GATKException",
                    format!(
                        "Exception thrown at {}:{} {}",
                        record.contig,
                        record.start,
                        java_variant_context_string(record, &paths[input])
                    ),
                )
            })?;
        }
        queue.push(record);
    }
    consume(&mut queue, &mut data, &mut random)?;
    manager.active = (0..data.len()).collect();
    manager.data = data;

    // `onTraversalSuccess`, attempt after attempt over the same data.
    let output_model = argument(parser, "output-model");
    let mut last_error = None;
    for attempt in 1..=settings.max_attempts.max(1) {
        match train_and_write(
            parser,
            &settings,
            &mut manager,
            loaded.as_ref(),
            output_model.as_deref(),
            rscript.as_deref(),
            &contigs,
            &mut random,
        ) {
            Ok(outputs) => {
                std::fs::write(&tranches_path, outputs.tranches.as_bytes()).map_err(|error| {
                    Thrown::non_user(PORT_FAILURE, format!("{tranches_path}: {error}"))
                })?;
                write_variant_output(parser, &output, &outputs.recal)?;
                if let Some((path, text)) = outputs.rscript {
                    std::fs::write(&path, text.as_bytes()).map_err(|error| {
                        Thrown::non_user(PORT_FAILURE, format!("{path}: {error}"))
                    })?;
                }
                // `onTraversalSuccess` returns `true`, which `Main` prints.
                return Ok(Some("true".to_string()));
            }
            Err(error) => {
                if attempt >= settings.max_attempts.max(1) {
                    return Err(error);
                }
                last_error = Some(error);
            }
        }
    }
    Err(last_error.unwrap_or_else(|| Thrown::non_user(PORT_FAILURE, "no attempt was made")))
}

/// What one successful attempt leaves behind.
struct Outputs {
    tranches: String,
    recal: String,
    rscript: Option<(String, String)>,
}

/// One pass of `onTraversalSuccess`'s loop body.
#[allow(clippy::too_many_arguments)]
fn train_and_write(
    parser: &Parser,
    settings: &Settings,
    manager: &mut model::DataManager,
    loaded: Option<&LoadedModel>,
    output_model: Option<&str>,
    rscript: Option<&str>,
    contigs: &[(String, i32)],
    random: &mut gatk_engine::java_random::JavaRandom,
) -> Result<Outputs, Thrown> {
    let arguments = &settings.arguments;
    manager
        .normalize(
            loaded.is_none(),
            loaded.map(|loaded| loaded.order.as_slice()),
            arguments.std_threshold,
            random,
        )
        .map_err(thrown)?;
    let mut positive = manager.training_data(arguments, random);
    let num_annotations = manager.annotation_keys.len();

    let (mut good, mut bad, mut negative) = match loaded {
        Some(loaded) => {
            let mut good = model::mixture_from_tables(
                loaded.table("PositiveModelMeans")?,
                loaded.table("PositiveModelCovariances")?,
                loaded.table("GoodGaussianPMix")?,
                num_annotations,
                positive.len(),
                arguments,
            );
            let active = manager.active.clone();
            model::evaluate_data(&mut manager.data, Some(&active), &mut good, false, random)
                .map_err(thrown)?;
            let negative = manager.select_worst(arguments.bad_lod_cutoff);
            let bad = model::mixture_from_tables(
                loaded.table("NegativeModelMeans")?,
                loaded.table("NegativeModelCovariances")?,
                loaded.table("BadGaussianPMix")?,
                num_annotations,
                negative.len(),
                arguments,
            );
            (good, bad, negative)
        }
        None => {
            let training: Vec<&model::Datum> = positive.iter().map(|i| &manager.data[*i]).collect();
            let mut good =
                model::generate_model(&training, arguments.max_gaussians, arguments, random)
                    .map_err(thrown)?;
            let active = manager.active.clone();
            model::evaluate_data(&mut manager.data, Some(&active), &mut good, false, random)
                .map_err(thrown)?;
            if good.failed_to_converge {
                if let Some(path) = output_model {
                    let keys = manager.annotation_keys.clone();
                    write_model(path, &model::model_report(manager, &good, None, &keys))?;
                }
                return Err(Thrown {
                    failure: Failure::User,
                    exception: model::POSITIVE_FAILURE,
                    message: Some(
                        "Positive training model failed to converge.  One or more annotations \
                         (usually MQ) may have insufficient variance.  Please consider lowering \
                         the maximum number of Gaussians allowed for use in the model (via \
                         --max-gaussians 4, for example)."
                            .to_string(),
                    ),
                });
            }
            let negative = manager.select_worst(arguments.bad_lod_cutoff);
            let worst: Vec<&model::Datum> = negative.iter().map(|i| &manager.data[*i]).collect();
            let bad = model::generate_model(
                &worst,
                arguments
                    .max_negative_gaussians
                    .min(arguments.max_gaussians),
                arguments,
                random,
            )
            .map_err(thrown)?;
            if bad.failed_to_converge {
                return Err(Thrown {
                    failure: Failure::User,
                    exception: model::NEGATIVE_FAILURE,
                    message: Some(
                        "NaN LOD value assigned. Clustering with this few variants and these \
                         annotations is unsafe. Please consider raising the number of variants \
                         used to train the negative model (via --minimum-bad-variants 5000, for \
                         example)."
                            .to_string(),
                    ),
                });
            }
            (good, bad, negative)
        }
    };

    // `dropAggregateData`: out of the list, still in the training lists that hold them.
    let data = &manager.data;
    let mut active = manager.active.clone();
    active.retain(|i| !data[*i].is_aggregate);
    manager.active = active.clone();

    model::evaluate_data(&mut manager.data, Some(&active), &mut bad, true, random)
        .map_err(thrown)?;
    if let Some(path) = output_model {
        let keys = manager.annotation_keys.clone();
        write_model(
            path,
            &model::model_report(manager, &good, Some(&bad), &keys),
        )?;
    }
    model::worst_performing_annotation(&mut manager.data, &active, &good, &bad).map_err(thrown)?;

    // `findTranches` sorts the list by score IN PLACE, stably.
    let calls_at_truth = active
        .iter()
        .filter(|i| manager.data[**i].at_truth_site)
        .count();
    let data = &manager.data;
    active.sort_by(|a, b| model::double_compare(data[*a].lod, data[*b].lod));
    let tranche_data: Vec<tranches::Datum> = active
        .iter()
        .map(|i| {
            let datum = &manager.data[*i];
            tranches::Datum {
                lod: datum.lod,
                is_known: datum.is_known,
                at_truth_site: datum.at_truth_site,
                is_snp: datum.is_snp,
                is_transition: datum.is_transition,
            }
        })
        .collect();
    let tranches_text = if !settings.scatter {
        let found = tranches::find_tranches(
            &tranche_data,
            &settings.targets,
            calls_at_truth,
            settings.mode,
        )
        .map_err(Thrown::user)?;
        tranches::truth_sensitivity_file(&found)
    } else {
        let found =
            tranches::find_vqslod_tranches(&tranche_data, &settings.vqslod_tranches, settings.mode);
        tranches::vqslod_file(&found)
    };

    // `writeOutRecalibrationTable`: the list sorted back by coordinate, stably.
    let contig_index = |name: &str| {
        contigs
            .iter()
            .position(|(contig, _)| contig == name)
            .map_or(i64::MAX, |index| index as i64)
    };
    let data = &manager.data;
    active.sort_by_key(|i| {
        let (contig, start, end) = data[*i].loc.clone().unwrap_or_default();
        (contig_index(&contig), start, end)
    });
    manager.active = active;
    let recal = recal_file(parser, manager, contigs)?;

    let rscript = match rscript {
        None => None,
        Some(path) => {
            let mut evaluation = manager.evaluation_data();
            let plotted = model::random_data_for_plotting(
                1000,
                &mut positive,
                &mut negative,
                &mut evaluation,
                random,
            );
            let text = visualization_script(path, manager, &plotted, &mut good, &mut bad, random)
                .map_err(thrown)?;
            Some((path.to_string(), text))
        }
    };
    Ok(Outputs {
        tranches: tranches_text,
        recal,
        rscript,
    })
}

fn write_model(path: &str, report: &str) -> Result<(), Thrown> {
    std::fs::write(path, report.as_bytes()).map_err(|error| {
        Thrown::user(format!(
            "Couldn't write file {path} because Exception writing to report output.  with \
             exception {error}"
        ))
    })
}

/// The recal file: one record per datum, at its locus, with its score, its culprit and its labels.
fn recal_file(
    parser: &Parser,
    manager: &model::DataManager,
    contigs: &[(String, i32)],
) -> Result<String, Thrown> {
    use htsjdk_vcf::allele::Allele;
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
    use htsjdk_vcf::variant::Value;
    let compound = |id: &str, number: Cardinality, line_type: LineType, description: &str| {
        HeaderLine::Compound {
            key: "INFO".to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: description.to_string(),
            extra: Vec::new(),
        }
    };
    let mut lines = vec![
        compound(
            "END",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Stop position of the interval",
        ),
        compound(
            "VQSLOD",
            Cardinality::Fixed(1),
            LineType::Float,
            "Log odds of being a true variant versus being false under the trained gaussian mixture model",
        ),
        compound(
            "culprit",
            Cardinality::Fixed(1),
            LineType::String,
            "The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out",
        ),
        compound(
            "POSITIVE_TRAIN_SITE",
            Cardinality::Fixed(0),
            LineType::Flag,
            "This variant was used to build the positive training set of good variants",
        ),
        compound(
            "NEGATIVE_TRAIN_SITE",
            Cardinality::Fixed(0),
            LineType::Flag,
            "This variant was used to build the negative training set of bad variants",
        ),
        HeaderLine::Filter {
            id: "PASS".to_string(),
            description: "Site contains at least one allele that passes filters".to_string(),
        },
    ];
    lines.extend(default_tool_vcf_header_lines(parser, "VariantRecalibrator"));
    for (index, (name, length)) in contigs.iter().enumerate() {
        lines.push(HeaderLine::Contig {
            index: index as i32,
            fields: vec![
                ("ID".to_string(), name.clone()),
                ("length".to_string(), length.to_string()),
            ],
        });
    }
    let header = update_header_contig_lines(
        parser,
        VcfHeader {
            lines,
            samples: Vec::new(),
        },
    )?;
    let placeholder = || {
        vec![
            Allele::from_str("N", true).expect("a reference allele"),
            Allele::from_str("<VQSR>", false).expect("a symbolic allele"),
        ]
    };
    let mut records = Vec::with_capacity(manager.data.len());
    for datum in manager.active.iter().map(|i| &manager.data[*i]) {
        let Some((contig, start, end)) = &datum.loc else {
            continue;
        };
        let alleles = match (&datum.reference_allele, &datum.alternate_allele) {
            (Some(reference), Some(alternate)) if settings_as(manager) => vec![
                Allele::from_str(reference, true).expect("a reference allele"),
                Allele::from_str(alternate, false).expect("an alternate allele"),
            ],
            _ => placeholder(),
        };
        let mut record = VariantContext::new(contig, *start, alleles);
        record.stop = *end;
        record
            .attributes
            .push(("END".to_string(), Value::Int(*end)));
        record.attributes.push((
            "VQSLOD".to_string(),
            Value::Str(gatk_engine::java_format::format_decimals(datum.lod, 4)),
        ));
        let culprit = if datum.worst_annotation != -1 {
            manager.annotation_keys[datum.worst_annotation as usize].clone()
        } else {
            "NULL".to_string()
        };
        record
            .attributes
            .push(("culprit".to_string(), Value::Str(culprit)));
        if datum.at_training_site {
            record
                .attributes
                .push(("POSITIVE_TRAIN_SITE".to_string(), Value::Bool(true)));
        }
        if datum.at_anti_training_site {
            record
                .attributes
                .push(("NEGATIVE_TRAIN_SITE".to_string(), Value::Bool(true)));
        }
        records.push(record);
    }
    htsjdk_vcf::vcf_file::write_vcf(&header, &records)
        .map_err(|error| Thrown::user(format!("{error:?}")))
}

/// Whether the data are allele-specific, which every datum's alternate answers alike.
fn settings_as(manager: &model::DataManager) -> bool {
    manager
        .data
        .first()
        .is_some_and(|datum| datum.alternate_allele.is_some())
}

/// `createVisualizationScript`: for each pair of annotations, the model's contrastive score over a
/// 61 by 61 grid, then the plotted data, as R vectors.
fn visualization_script(
    path: &str,
    manager: &model::DataManager,
    plotted: &[usize],
    good: &mut model::Mixture,
    bad: &mut model::Mixture,
    random: &mut gatk_engine::java_random::JavaRandom,
) -> Result<String, model::ModelError> {
    use gatk_engine::java_format::format_decimals as f4;
    let fmt = |value: f64| f4(value, 4);
    let keys = &manager.annotation_keys;
    let mut s = String::new();
    s.push_str("library(ggplot2)\nlibrary(tools)\nlibrary(grid)\n");
    for line in [
        "vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)",
        "arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {",
        "dots <- list(...)",
        "n <- length(dots)",
        "if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}",
        "if(is.null(nrow)) { nrow = ceiling(n/ncol)}",
        "if(is.null(ncol)) { ncol = ceiling(n/nrow)}",
        "grid.newpage()",
        "pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )",
        "ii.p <- 1",
        "for(ii.row in seq(1, nrow)){",
        "ii.table.row <- ii.row ",
        "if(as.table) {ii.table.row <- nrow - ii.table.row + 1}",
        "for(ii.col in seq(1, ncol)){",
        "ii.table <- ii.p",
        "if(ii.p > n) break",
        "print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))",
        "ii.p <- ii.p + 1",
        "}",
        "}",
        "}",
    ] {
        s.push_str(line);
        s.push('\n');
    }
    s.push_str(&format!("outputPDF <- \"{path}.pdf\"\n"));
    s.push_str("pdf(outputPDF)\n");
    let n = keys.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let (mut min1, mut max1, mut min2, mut max2) = (100.0, -100.0, 100.0, -100.0);
            for index in plotted {
                let datum = &manager.data[*index];
                min1 = model::java_min(min1, datum.annotations[i]);
                max1 = model::java_max(max1, datum.annotations[i]);
                min2 = model::java_min(min2, datum.annotations[j]);
                max2 = model::java_max(max2, datum.annotations[j]);
            }
            let mut fake: Vec<model::Datum> = Vec::new();
            let mut a1 = min1;
            while a1 <= max1 {
                let mut a2 = min2;
                while a2 <= max2 {
                    let mut datum = model::Datum::new(vec![0.0; n], vec![true; n]);
                    datum.prior = 0.0;
                    datum.annotations[i] = a1;
                    datum.annotations[j] = a2;
                    datum.is_null[i] = false;
                    datum.is_null[j] = false;
                    fake.push(datum);
                    a2 += (max2 - min2) / 60.0;
                }
                a1 += (max1 - min1) / 60.0;
            }
            model::evaluate_data(&mut fake, None, good, false, random)?;
            model::evaluate_data(&mut fake, None, bad, true, random)?;
            s.push_str("surface <- c(");
            for datum in &fake {
                s.push_str(&format!(
                    "{}, {}, {}, ",
                    fmt(manager.denormalize(datum.annotations[i], i)),
                    fmt(manager.denormalize(datum.annotations[j], j)),
                    fmt(model::java_min(4.0, model::java_max(-4.0, datum.lod)))
                ));
            }
            s.push_str("NA,NA,NA)\n");
            s.push_str("s <- matrix(surface,ncol=3,byrow=T)\n");
            s.push_str("data <- c(");
            for index in plotted {
                let datum = &manager.data[*index];
                let training = if datum.at_anti_training_site {
                    -1
                } else if datum.at_training_site {
                    1
                } else {
                    0
                };
                s.push_str(&format!(
                    "{}, {}, {}, {}, {},",
                    fmt(manager.denormalize(datum.annotations[i], i)),
                    fmt(manager.denormalize(datum.annotations[j], j)),
                    fmt(if datum.lod < 0.0 { -1.0 } else { 1.0 }),
                    training,
                    if datum.is_known { 1 } else { -1 }
                ));
            }
            s.push_str("NA,NA,NA,NA,1)\n");
            s.push_str("d <- matrix(data,ncol=5,byrow=T)\n");
            let (k1, k2) = (&keys[i], &keys[j]);
            let surface = format!("sf.{k1}.{k2}");
            let frame = format!("df.{k1}.{k2}");
            let white = "+theme(panel.background = element_rect(fill = \"white\"), panel.grid.minor = element_line(colour = \"white\"), panel.grid.major = element_line(colour = \"white\"))";
            s.push_str(&format!(
                "{surface} <- data.frame(x=s[,1], y=s[,2], lod=s[,3])\n"
            ));
            s.push_str(&format!(
                "{frame} <- data.frame(x=d[,1], y=d[,2], retained=d[,3], training=d[,4], novelty=d[,5])\n"
            ));
            s.push_str(&format!("dummyData <- {frame}[1,]\n"));
            s.push_str("dummyData$x <- NaN\n");
            s.push_str("dummyData$y <- NaN\n");
            s.push_str(&format!(
                "p <- ggplot(data={surface}, aes(x=x, y=y)) {white}\n"
            ));
            s.push_str(&format!("p1 = p +ggtitle(\"model PDF\") + labs(x=\"{k1}\", y=\"{k2}\") + geom_tile(aes(fill = lod)) + scale_fill_gradient(high=\"green\", low=\"red\", space=\"Lab\")\n"));
            s.push_str(&format!(
                "p <- qplot(x,y,data={frame}, color=retained, alpha=I(1/7),legend=FALSE) {white}\n"
            ));
            s.push_str("q <- geom_point(aes(x=x,y=y,color=retained),data=dummyData, alpha=1.0, na.rm=TRUE)\n");
            s.push_str(&format!("p2 = p + q + labs(x=\"{k1}\", y=\"{k2}\") + scale_colour_gradient(name=\"outcome\", high=\"black\", low=\"red\",breaks=c(-1,1),guide=\"legend\",labels=c(\"filtered\",\"retained\"))\n"));
            s.push_str(&format!("p <- qplot(x,y,data={frame}[{frame}$training != 0,], color=training, alpha=I(1/7)) {white}\n"));
            s.push_str("q <- geom_point(aes(x=x,y=y,color=training),data=dummyData, alpha=1.0, na.rm=TRUE)\n");
            s.push_str(&format!("p3 = p + q + labs(x=\"{k1}\", y=\"{k2}\") + scale_colour_gradient(high=\"green\", low=\"purple\",breaks=c(-1,1),guide=\"legend\", labels=c(\"neg\", \"pos\"))\n"));
            s.push_str(&format!(
                "p <- qplot(x,y,data={frame}, color=novelty, alpha=I(1/7)) {white}\n"
            ));
            s.push_str("q <- geom_point(aes(x=x,y=y,color=novelty),data=dummyData, alpha=1.0, na.rm=TRUE)\n");
            s.push_str(&format!("p4 = p + q + labs(x=\"{k1}\", y=\"{k2}\") + scale_colour_gradient(name=\"novelty\", high=\"blue\", low=\"red\",breaks=c(-1,1),guide=\"legend\", labels=c(\"novel\",\"known\"))\n"));
            s.push_str("arrange(p1, p2, p3, p4, ncol=2)\n");
        }
    }
    s.push_str("dev.off()\n");
    s.push_str("if (exists(\"compactPDF\")) {\n");
    s.push_str("compactPDF(outputPDF)\n");
    s.push_str("}\n");
    Ok(s)
}
