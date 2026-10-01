//! `ExtractVariantAnnotations`: the annotations of the variants a set of resources labels, written
//! as an HDF5 matrix, a reservoir of unlabelled ones beside it, and a sites-only VCF of the
//! labelled sites.
//!
//! The decisions are [`gatk_tools::extract_variant_annotations`]'s and the HDF5 file is
//! [`gatk_tools::hdf5_writer`]'s. What is here is the walker, in the reference's order:
//!
//! * **`onTraversalStart`**: the output names checked, the resources' labels gathered (the
//!   reserved `snp` refused), allele-specific mode switched on when any requested annotation is
//!   declared `Number=A` (a requested annotation the header does not declare is the reference's
//!   `NullPointerException`), the matching strategy forced to the minimal representation in that
//!   mode, and the VCF's header written;
//! * **the traversal**: each record's alternates checked for filters and type, labelled by the
//!   resource records that START where it starts, and either added to the labelled data and the
//!   VCF or offered to the reservoir by Algorithm R, whose draws come from
//!   `new Random(--reservoir-sampling-random-seed)`;
//! * **`afterNthPass`**: the labelled matrix (skipped with a warning when nothing was labelled),
//!   then the reservoir's, which refuses an empty reservoir AFTER the first file is written.
use super::*;
use gatk_tools::extract_variant_annotations as extract;
use gatk_tools::hdf5_writer::Hdf5File;
use std::collections::BTreeSet;

/// One extracted record: its alternates, type and labels, as `Triple`s per datum.
type Metadata = Vec<extract::Extracted>;

/// A resource: its labels and its records, each with the genotype facts the polymorphism check
/// reads.
struct Resource {
    labels: BTreeSet<String>,
    records: Vec<(extract::Record, bool, bool)>,
}

/// The port's reduced record, from a decoded one.
fn reduced(record: &htsjdk_vcf::variant::VariantContext) -> extract::Record {
    use htsjdk_vcf::variant::Value;
    let attributes = record
        .attributes
        .iter()
        .filter_map(|(key, value)| {
            let text = match value {
                Value::Str(text) => text.clone(),
                Value::List(values) => values
                    .iter()
                    .map(|value| match value {
                        Value::Str(text) => text.clone(),
                        other => other.format().unwrap_or_default(),
                    })
                    .collect::<Vec<_>>()
                    .join(","),
                Value::Bool(false) => return None,
                other => other.format().unwrap_or_default(),
            };
            Some((key.clone(), text))
        })
        .collect();
    extract::Record {
        contig: record.contig.clone(),
        start: record.start as i32,
        reference: record.reference().base_string(),
        alternates: record
            .alternate_alleles()
            .iter()
            .map(|allele| allele.display_string())
            .collect(),
        filters: record.filters.clone().unwrap_or_default(),
        attributes,
    }
}

/// `hasGenotypes` and `isPolymorphicInSamples` of a resource record.
fn genotype_facts(record: &htsjdk_vcf::variant::VariantContext) -> (bool, bool) {
    let has = record.genotypes.iter().next().is_some();
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
    (has, record.alleles.len() > 1 && called != reference)
}

/// The labels the resources give one record, or one of its alternates.
fn matching_labels(
    record: &extract::Record,
    alternate: Option<&str>,
    resources: &[Resource],
    arguments: &extract::Arguments,
) -> BTreeSet<String> {
    let mut labels = BTreeSet::new();
    for resource in resources {
        for (candidate, has_genotypes, polymorphic) in &resource.records {
            if candidate.contig != record.contig || candidate.start != record.start {
                continue;
            }
            if extract::is_matching_variant(
                record,
                candidate,
                alternate,
                arguments.trust_all_polymorphic,
                *polymorphic,
                *has_genotypes,
                arguments.strategy,
            ) {
                labels.extend(resource.labels.iter().cloned());
            }
        }
    }
    labels
}

/// `extractVariantMetadata`.
fn metadata(
    record: &extract::Record,
    resources: &[Resource],
    arguments: &extract::Arguments,
    extract_unlabeled: bool,
) -> Result<Metadata, Thrown> {
    if !extract::passes_filters(record, arguments) {
        return Ok(Vec::new());
    }
    if !arguments.allele_specific {
        let Some(kind) = extract::variant_type(record) else {
            return Err(Thrown::non_user(
                "java.lang.IllegalStateException",
                "Encountered unknown variant type: NO_VARIATION",
            ));
        };
        if arguments.modes.contains(&kind) {
            let labels = matching_labels(record, None, resources, arguments);
            if extract_unlabeled || !labels.is_empty() {
                return Ok(vec![extract::Extracted {
                    alternates: record.alternates.clone(),
                    variant_type: kind,
                    labels,
                }]);
            }
        }
        return Ok(Vec::new());
    }
    Ok(record
        .alternates
        .iter()
        .filter(|alternate| *alternate != "*")
        .filter(|alternate| {
            arguments
                .modes
                .contains(&extract::allele_specific_variant_type(
                    &record.reference,
                    alternate,
                ))
        })
        .map(|alternate| extract::Extracted {
            alternates: vec![alternate.clone()],
            variant_type: extract::allele_specific_variant_type(&record.reference, alternate),
            labels: matching_labels(record, Some(alternate), resources, arguments),
        })
        .filter(|extracted| extract_unlabeled || !extracted.labels.is_empty())
        .collect())
}

/// `LabeledVariantAnnotationsData.writeHDF5`: the rows of every record, flattened in order.
fn matrix(
    records: &[Vec<extract::Row>],
    names: &[String],
    labels: &[String],
    allele_specific: bool,
    omit_alleles: bool,
) -> Vec<u8> {
    let rows: Vec<&extract::Row> = records.iter().flatten().collect();
    let mut file = Hdf5File::new();
    // `writeIntervals`: contigs indexed in the order they first appear.
    let mut contigs: Vec<String> = Vec::new();
    let mut index = Vec::with_capacity(rows.len());
    for row in &rows {
        let at = match contigs.iter().position(|contig| *contig == row.contig) {
            Some(at) => at,
            None => {
                contigs.push(row.contig.clone());
                contigs.len() - 1
            }
        };
        index.push(at as f64);
    }
    file.double_matrix(
        "/intervals/transposed_index_start_end",
        &[
            index,
            rows.iter().map(|row| row.start as f64).collect(),
            rows.iter().map(|row| row.end as f64).collect(),
        ],
    );
    file.string_array("/intervals/indexed_contig_names", &contigs);
    if !omit_alleles {
        file.string_array(
            "/alleles/ref",
            &rows
                .iter()
                .map(|row| row.reference.clone())
                .collect::<Vec<_>>(),
        );
        let alternates: Vec<String> = rows
            .iter()
            .map(|row| {
                if allele_specific {
                    row.alternates[0].clone()
                } else {
                    row.alternates.join(",")
                }
            })
            .collect();
        file.string_array("/alleles/alt", &alternates);
    }
    file.string_array("/annotations/names", names);
    // `writeChunkedDoubleMatrix`, which needs one chunk for any matrix this corpus makes.
    file.double_array("/annotations/num_rows", &[rows.len() as f64]);
    file.double_array("/annotations/num_columns", &[names.len() as f64]);
    file.double_array("/annotations/num_chunks", &[1.0]);
    file.double_matrix(
        "/annotations/chunk_0",
        &rows
            .iter()
            .map(|row| row.annotations.clone())
            .collect::<Vec<_>>(),
    );
    file.double_array(
        "/labels/snp",
        &rows
            .iter()
            .map(|row| {
                if row.variant_type == extract::VariantType::Snp {
                    1.0
                } else {
                    0.0
                }
            })
            .collect::<Vec<_>>(),
    );
    for label in labels {
        file.double_array(
            &format!("/labels/{label}"),
            &rows
                .iter()
                .map(|row| if row.labels.contains(label) { 1.0 } else { 0.0 })
                .collect::<Vec<_>>(),
        );
    }
    file.to_bytes()
}

pub fn extract_variant_annotations(parser: &Parser) -> Outcome {
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup(parser, "ExtractVariantAnnotations")?;
    let prefix = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let file = htsjdk_vcf::reader::read_vcf(&text).map_err(|failure| Thrown {
        failure: Failure::User,
        exception: "htsjdk.tribble.TribbleException",
        message: Some(failure.error.message()),
    })?;

    // The resources are feature inputs, opened at startup.
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
    let mut resources = Vec::new();
    for (path, attributes) in &resource_values {
        let (_, resource_text) = open_feature_input(path)?;
        let resource_file =
            htsjdk_vcf::reader::read_vcf(&resource_text).map_err(|failure| Thrown {
                failure: Failure::User,
                exception: "htsjdk.tribble.TribbleException",
                message: Some(failure.error.message()),
            })?;
        resources.push((path.clone(), attributes.clone(), resource_file));
    }

    // `onTraversalStart`.
    let mut ignored: BTreeSet<String> = arguments(parser, "ignore-filter").into_iter().collect();
    ignored.remove("null");
    let modes: BTreeSet<extract::VariantType> = {
        let listed = arguments(parser, "mode");
        let listed = if listed.is_empty() {
            match scalar(parser, "mode") {
                Some(value) => value
                    .trim_matches(|c| c == '[' || c == ']')
                    .split(", ")
                    .map(str::to_string)
                    .collect(),
                None => vec!["SNP".to_string(), "INDEL".to_string()],
            }
        } else {
            listed
        };
        listed
            .iter()
            .filter_map(|mode| match mode.as_str() {
                "SNP" => Some(extract::VariantType::Snp),
                "INDEL" => Some(extract::VariantType::Indel),
                _ => None,
            })
            .collect()
    };
    let gzip = !flag(parser, "do-not-gzip-vcf-output");
    let annotations_path = format!("{prefix}{}", extract::ANNOTATIONS_HDF5_SUFFIX);
    let vcf_path = format!("{prefix}{}", if gzip { ".vcf.gz" } else { ".vcf" });
    for path in [&annotations_path, &vcf_path] {
        let parent = std::path::Path::new(path)
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or(std::path::Path::new("."));
        if !parent.is_dir() {
            return Err(Thrown::user(format!(
                "Couldn't write file {} because The output file could not be written.",
                java_absolute_path(path)
            )));
        }
    }
    let mut labels: BTreeSet<String> = BTreeSet::new();
    let mut loaded = Vec::new();
    for (_, attributes, resource_file) in resources {
        let own: BTreeSet<String> = attributes
            .iter()
            .filter(|(_, value)| value == "true")
            .map(|(key, _)| key.clone())
            .collect();
        labels.extend(own.iter().cloned());
        loaded.push(Resource {
            labels: own,
            records: resource_file
                .records
                .iter()
                .map(|record| {
                    let (has, polymorphic) = genotype_facts(record);
                    (reduced(record), has, polymorphic)
                })
                .collect(),
        });
    }
    extract::check_resource_labels(&labels).map_err(|message| Thrown {
        failure: Failure::User,
        exception: "org.broadinstitute.hellbender.exceptions.UserException$BadInput",
        message: Some(message),
    })?;

    // `isAlleleSpecificAnnotationRequested`: the first requested name the header does not declare
    // is a `NullPointerException`, unless a `Number=A` one came first.
    let requested: Vec<String> = arguments(parser, "annotation");
    let mut distinct: Vec<String> = Vec::new();
    for name in &requested {
        if !distinct.contains(name) {
            distinct.push(name.clone());
        }
    }
    let mut allele_specific = false;
    for name in &distinct {
        let line = file.header.lines.iter().find_map(|line| match line {
            htsjdk_vcf::header::HeaderLine::Compound {
                key, id, number, ..
            } if key == "INFO" && id == name => Some(*number),
            _ => None,
        });
        match line {
            None => {
                return Err(Thrown::non_user(
                    "java.lang.NullPointerException",
                    "Cannot invoke \"htsjdk.variant.vcf.VCFInfoHeaderLine.getCountType()\" because \
                     the return value of \"htsjdk.variant.vcf.VCFHeader.getInfoHeaderLine(String)\" \
                     is null",
                ))
            }
            Some(htsjdk_vcf::header::Cardinality::A) => {
                allele_specific = true;
                break;
            }
            Some(_) => {}
        }
    }
    let strategy = match scalar(parser, "resource-matching-strategy").as_deref() {
        Some("START_POSITION_AND_GIVEN_REPRESENTATION") => {
            extract::MatchingStrategy::StartPositionAndGivenRepresentation
        }
        Some("START_POSITION_AND_MINIMAL_REPRESENTATION") => {
            extract::MatchingStrategy::StartPositionAndMinimalRepresentation
        }
        _ => extract::MatchingStrategy::StartPosition,
    };
    let arguments = extract::Arguments {
        modes,
        ignored_filters: ignored,
        ignore_all_filters: flag(parser, "ignore-all-filters"),
        trust_all_polymorphic: !flag(parser, "do-not-trust-all-polymorphic"),
        strategy: if allele_specific {
            extract::MatchingStrategy::StartPositionAndMinimalRepresentation
        } else {
            strategy
        },
        allele_specific,
    };
    let names = extract::sorted_annotation_names(&requested);
    let sorted_labels: Vec<String> = labels.iter().cloned().collect();

    // The VCF's header, written before the first record.
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
    let mut lines: Vec<HeaderLine> = sorted_labels
        .iter()
        .map(|label| HeaderLine::Compound {
            key: "INFO".to_string(),
            id: label.clone(),
            number: Cardinality::Fixed(0),
            line_type: LineType::Flag,
            description: format!("This site was labeled as {label} according to resources"),
            extra: Vec::new(),
        })
        .collect();
    lines.push(HeaderLine::Filter {
        id: "PASS".to_string(),
        description: "Site contains at least one allele that passes filters".to_string(),
    });
    let dictionary = vcf_dictionary(&text);
    for (index, sequence) in dictionary.sequences.iter().enumerate() {
        lines.push(HeaderLine::Contig {
            index: index as i32,
            fields: vec![
                ("ID".to_string(), sequence.name.clone()),
                ("length".to_string(), sequence.length.to_string()),
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
    let mut header = header;
    header.lines.extend(default_tool_vcf_header_lines(
        parser,
        "ExtractVariantAnnotations",
    ));

    let maximum = number_or(parser, "maximum-number-of-unlabeled-variants", 0).max(0) as usize;
    let seed = number_or(parser, "reservoir-sampling-random-seed", 0) as i64;
    let mut random = gatk_engine::java_random::JavaRandom::new(seed);

    // The traversal.
    let traversed = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let mut labeled: Vec<Vec<extract::Row>> = Vec::new();
    let mut reservoir: Vec<Vec<extract::Row>> = Vec::new();
    let mut unlabeled_index = 0usize;
    let mut written: Vec<htsjdk_vcf::variant::VariantContext> = Vec::new();
    for record in traversed {
        let reduced_record = reduced(record);
        let found = metadata(&reduced_record, &loaded, &arguments, maximum > 0)?;
        if found.is_empty() {
            continue;
        }
        let labeled_found: Vec<&extract::Extracted> =
            found.iter().filter(|e| !e.labels.is_empty()).collect();
        if !labeled_found.is_empty() {
            labeled.push(
                labeled_found
                    .iter()
                    .map(|e| extract::row(&reduced_record, e, &names, allele_specific))
                    .collect(),
            );
            // `writeExtractedVariantToVCF`: the reference and the labelled alternates, flagged
            // with the union of their labels.
            let mut alleles = vec![record.reference().clone()];
            for extracted in &labeled_found {
                for alternate in &extracted.alternates {
                    if let Some(allele) = record
                        .alternate_alleles()
                        .iter()
                        .find(|allele| allele.display_string() == *alternate)
                    {
                        alleles.push(allele.clone());
                    }
                }
            }
            let mut out =
                htsjdk_vcf::variant::VariantContext::new(&record.contig, record.start, alleles);
            out.stop = record.stop;
            let union: BTreeSet<&String> = labeled_found.iter().flat_map(|e| &e.labels).collect();
            for label in union {
                out.attributes
                    .push((label.clone(), htsjdk_vcf::variant::Value::Bool(true)));
            }
            written.push(out);
        }
        if maximum > 0 {
            let unlabeled: Vec<extract::Row> = found
                .iter()
                .filter(|e| e.labels.is_empty())
                .map(|e| extract::row(&reduced_record, e, &names, allele_specific))
                .collect();
            if !unlabeled.is_empty() {
                if unlabeled_index < maximum {
                    reservoir.push(unlabeled);
                } else {
                    let j = random.next_int_bound(unlabeled_index as i32) as usize;
                    if j < maximum {
                        reservoir[j] = unlabeled;
                    }
                }
                unlabeled_index += 1;
            }
        }
    }

    // `afterNthPass`: the labelled matrix, the reservoir's, and the VCF closed.
    let omit = flag(parser, "omit-alleles-in-hdf5");
    if !labeled.is_empty() {
        write_file(
            &annotations_path,
            &matrix(&labeled, &names, &sorted_labels, allele_specific, omit),
        )?;
    }
    let finish = |written: &[htsjdk_vcf::variant::VariantContext]| -> Result<(), Thrown> {
        let text = htsjdk_vcf::vcf_file::write_vcf(&header, written)
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
        write_variant_output(parser, &vcf_path, &text)
    };
    if maximum > 0 {
        if reservoir.is_empty() {
            finish(&written)?;
            return Err(Thrown::non_user(
                "org.broadinstitute.hellbender.exceptions.GATKException",
                "No unlabeled variants were present in the input VCF. Consider setting the \
                 maximum-number-of-unlabeled-variants argument to 0.",
            ));
        }
        write_file(
            &format!(
                "{prefix}{}{}",
                extract::UNLABELED_TAG,
                extract::ANNOTATIONS_HDF5_SUFFIX
            ),
            &matrix(&reservoir, &names, &sorted_labels, allele_specific, omit),
        )?;
    }
    finish(&written)?;
    Ok(None)
}
