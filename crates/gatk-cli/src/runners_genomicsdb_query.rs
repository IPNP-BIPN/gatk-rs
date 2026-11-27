//! A `gendb://` workspace as the driving variants of `GenotypeGVCFs`, `GnarlyGenotyper` and
//! `SelectVariants`: `FeatureDataSource.getGenomicsDBFeatureReader` with the VCF codec.
//!
//! The text the native query streams is [`gatk_tools::genomicsdb_query`]'s; what is here is the
//! startup around it, in `GATKTool.onStartup`'s order (the reference, then the reader, then the
//! dictionaries and `-L`), and the queries `FeatureDataSource` makes: one per traversal interval
//! over each array it overlaps, or every array in `ChrArrayFolderComparator` order without `-L`.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is a value or nothing; `?` after one returns nothing from the function when
//!   there is no value, and after a `Result` it returns the error, Java's `throws`;
//! - `let Some(x) = y else { return ...; };` unwraps the value or takes the early exit;
//! - `|what: String| ...` is a closure, Java's lambda; one declared `let mut` may change what
//!   it captured (here, the reference reader it reads bases through);
//! - `x.max(0) as u32` clamps at zero and then converts; `as` never fails, so the clamp is what
//!   keeps a negative count from wrapping to a huge unsigned one.
use super::*;
use gatk_tools::genomicsdb_query::{ExportOptions, Workspace};

/// `IOUtils.getGenomicsDBPath`: what follows `gendb://`, or `scheme://...` for `gendb.scheme://`,
/// and `None` for an input that is not a GenomicsDB URI.
///
/// ```java
/// Pattern.compile("^" + GENOMIC_DB_URI_SCHEME + "(\\.?)(.*)(://)(.*)");
/// ```
///
/// ```java
/// if (!matcher.group(1).isEmpty()) { // path has a period after gendb, so it is a URI
///     if (!matcher.group(2).isEmpty()) { //path has a scheme, so it is valid URI for GenomicsDB
///         genomicsdbPath = matcher.group(2) + matcher.group(3) + matcher.group(4);
///     }
/// } else if (matcher.group(2).isEmpty()) {
///     genomicsdbPath = matcher.group(4);
/// }
/// ```
///
/// So `gendb.://x` (a dot and no scheme) and `gendbx://x` (no dot and a scheme) are both not
/// GenomicsDB paths; the tests below pin each case.
pub(super) fn genomicsdb_path(input: &str) -> Option<String> {
    // ^gendb(\.?)(.*)(://)(.*), whose greedy second group ends at the LAST `://`.
    let rest = input.strip_prefix("gendb")?;
    let (dotted, rest) = match rest.strip_prefix('.') {
        Some(rest) => (true, rest),
        None => (false, rest),
    };
    let at = rest.rfind("://")?;
    let (scheme, path) = (&rest[..at], &rest[at + 3..]);
    match (dotted, scheme.is_empty()) {
        (true, false) => Some(format!("{scheme}://{path}")),
        (false, true) => Some(path.to_string()),
        _ => None,
    }
}

/// `getGenomicsDBOptions` as each of the three tools overrides it, and the check the
/// `GenomicsDBOptions` constructor makes on the result.
///
/// ```java
/// genomicsdbArgs.maxDiploidAltAllelesThatCanBeGenotyped - 1 < genotypeCalcArgs.maxAlternateAlleles
/// ```
///
/// `SelectVariants` calls `new GenomicsDBOptions(referenceArguments.getReferencePath(),
/// genomicsdbArgs)`, which passes a fresh `GenotypeCalculationArgumentCollection`, so its check
/// and its genotype count use the defaults, `DEFAULT_MAX_ALTERNATE_ALLELES = 6` and
/// `DEFAULT_MAX_GENOTYPE_COUNT = 1024`, whatever the command line says.
fn export_options(parser: &Parser, tool: &str) -> Result<ExportOptions, Thrown> {
    let max_alternate_alleles = number_or(parser, "max-alternate-alleles", 6);
    let max_genotype_count = number_or(parser, "max-genotype-count", 1024);
    let (produce_gt, max_diploid) = match tool {
        // `genomicsdbArgs.maxDiploidAltAllelesThatCanBeGenotyped = maxAlternateAlleles + 1`.
        "GnarlyGenotyper" => (true, max_alternate_alleles + 1),
        // `genomicsdbArgs.callGenotypes = true;` before the options are built.
        "GenotypeGVCFs" => (
            true,
            number_or(parser, "genomicsdb-max-alternate-alleles", 50),
        ),
        // `new GenomicsDBOptions(reference, genomicsdbArgs)`: the default genotype arguments.
        _ => (
            flag(parser, "call-genotypes"),
            number_or(parser, "genomicsdb-max-alternate-alleles", 50),
        ),
    };
    let (checked_alternates, max_genotype_count) = if tool == "SelectVariants" {
        (6, 1024)
    } else {
        (max_alternate_alleles, max_genotype_count)
    };
    if max_diploid - 1 < checked_alternates {
        return Err(bad_input(
            "GenomicsDB max alternate alleles (genomicsdb-max-alternate-alleles) must be at least \
             one greater than genotype calculation max alternate alleles (max-alternate-alleles), \
             accounting for the non-ref allele"
                .to_string(),
        ));
    }
    Ok(ExportOptions {
        produce_gt,
        max_diploid_alt_alleles: max_diploid.max(0) as u32,
        max_genotype_count: max_genotype_count.max(0) as u32,
    })
}

/// A variant walker's startup over a GenomicsDB workspace: the header and dictionary the reader
/// reports, the usual validation, and the text of every query the traversal makes.
pub(super) fn genomicsdb_startup(
    parser: &Parser,
    tool: &str,
    input: String,
) -> Result<VariantWalkerStart, Thrown> {
    // The port's own refusal, for a workspace this reader does not handle; never GATK's text.
    let limitation = |what: String| {
        Thrown::non_user(
            PORT_LIMITATION,
            format!("{what} This message is the port's own and not GATK's."),
        )
    };
    // The caller only routes a `gendb` input here, so the URI parses; `expect` panics otherwise.
    let path = genomicsdb_path(&input).expect("a GenomicsDB URI");
    if path.contains("://") {
        return Err(limitation(format!(
            "GenomicsDB workspace {input} is not on a local file system, which this port's \
             GenomicsDB reader does not read."
        )));
    }
    if flag(parser, "genomicsdb-use-bcf-codec") {
        return Err(limitation(
            "--genomicsdb-use-bcf-codec streams the workspace through BCF2Codec, whose records \
             this port's GenomicsDB reader has not measured."
                .to_string(),
        ));
    }
    // `initializeReference` runs before `initializeFeatures` opens the reader.
    let _ = reference_dictionary(parser)?;
    // `getGenomicsDBOptions()` is an argument of `addToFeatureSources`, so its `BadInput` comes
    // before the reader is opened and before the reference check inside it.
    let options = export_options(parser, tool)?;
    // `throw new UserException.MissingReference("You must provide a reference if you want to load
    // from GenomicsDB");` (`FeatureDataSource.getFeatureReader`).
    let Some(reference) = argument(parser, "reference") else {
        return Err(Thrown {
            failure: Failure::User,
            exception: "org.broadinstitute.hellbender.exceptions.UserException$MissingReference",
            message: Some(
                "You must provide a reference if you want to load from GenomicsDB".into(),
            ),
        });
    };
    let directory = java_absolute_path(&path);
    if !std::path::Path::new(&directory).exists() {
        // The `FeatureInput` is named, its path made an absolute `gendb://` URI:
        // `throw new UserException("GenomicsDB workspace " + path + " does not exist");`
        return Err(Thrown::user(format!(
            "GenomicsDB workspace drivingVariantFile:gendb://{directory} does not exist"
        )));
    }
    let workspace = Workspace::open(std::path::Path::new(&directory)).map_err(limitation)?;
    let header = workspace.header().map_err(limitation)?;
    let dictionary = vcf_dictionary(&header);
    let mut start = variant_walker_validation(
        parser,
        input,
        header,
        gatk_tools::feature_codec::Codec::Vcf,
        dictionary,
    )?;

    // `transformTraversalIntervals`: `--merge-input-intervals` spans each contig's intervals,
    // `IntervalUtils.getSpanningIntervals`. Only `GenotypeGVCFs` and `GnarlyGenotyper` declare the
    // flag. The spans are kept in first-seen contig order, which for traversal intervals already
    // sorted by the dictionary is the order the reference sorts them back into.
    let traversal: Option<Vec<gatk_engine::interval::SimpleInterval>> =
        start.intervals.as_ref().map(|list| {
            if !flag(parser, "merge-input-intervals") {
                return list.clone();
            }
            let mut spanning: Vec<gatk_engine::interval::SimpleInterval> = Vec::new();
            for interval in list {
                match spanning
                    .iter_mut()
                    .find(|span| span.contig == interval.contig)
                {
                    Some(span) => {
                        span.start = span.start.min(interval.start);
                        span.end = span.end.max(interval.end);
                    }
                    None => spanning.push(interval.clone()),
                }
            }
            spanning
        });
    let mut source =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let mut bases = ReferenceBases(&mut source);
    let mut records = String::new();
    // One query's records appended to the text. The closure borrows `workspace`, `options` and
    // `bases` for as long as it lives; `records` is passed in rather than captured so the loop
    // below can keep appending to it.
    let mut query = |array, begin, end, records: &mut String| {
        workspace
            .stream(array, begin, end, &options, &mut bases, records)
            .map_err(limitation)
    };
    match &traversal {
        // `GenomicsDBFeatureReader.iterator`: every array over its whole span.
        None => {
            for array in workspace.arrays() {
                query(array, array.start, array.end, &mut records)?;
            }
        }
        // `GenomicsDBFeatureReader.query`: each array of the contig the interval overlaps,
        // clipped to it.
        Some(list) => {
            for interval in list {
                // `as i64` widens the interval's coordinates to the arrays' type; it cannot lose
                // a value.
                let (begin, end) = (interval.start as i64, interval.end as i64);
                for array in workspace.arrays().iter().filter(|array| {
                    array.contig == interval.contig && begin <= array.end && end >= array.start
                }) {
                    query(
                        array,
                        begin.max(array.start),
                        end.min(array.end),
                        &mut records,
                    )?;
                }
            }
        }
    }
    start.text.push_str(&records);
    Ok(start)
}

#[cfg(test)]
mod tests {
    use super::genomicsdb_path;

    #[test]
    fn genomicsdb_uris_resolve_as_ioutils_does() {
        assert_eq!(genomicsdb_path("gendb://ws").as_deref(), Some("ws"));
        assert_eq!(
            genomicsdb_path("gendb:///abs/ws").as_deref(),
            Some("/abs/ws")
        );
        assert_eq!(
            genomicsdb_path("gendb.gs://b/ws").as_deref(),
            Some("gs://b/ws")
        );
        assert_eq!(genomicsdb_path("gendb.://ws"), None);
        assert_eq!(genomicsdb_path("gendbx://ws"), None);
        assert_eq!(genomicsdb_path("ws.g.vcf"), None);
    }
}
