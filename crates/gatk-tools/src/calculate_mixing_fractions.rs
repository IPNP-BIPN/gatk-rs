//! `CalculateMixingFractions`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.CalculateMixingFractions`
//! (GATK 4.6.2.0).
//!
//! One bucket per sample, filled at singleton het SNPs from the pileup of a pooled bam, and
//! normalised at the end by the sum of every bucket's alt fraction.
//!
//! # The row order is a hash bucket order
//!
//! ```java
//! private final Map<String, AltAndTotalReadCounts> sampleCounts = new HashMap<>();
//! ...
//! sampleCounts.entrySet().stream().map(e -> new MixingFraction(e.getKey(), ...))
//! ```
//!
//! The rows come out in the order a `HashMap` iterates its keys, which is neither the header's
//! sample order nor alphabetical: the golden's three samples are written `zebra`, `mike`, `alpha`
//! for a header that lists `zebra`, `alpha`, `mike`. [`gatk_engine::java_hash::hash_map_order`] is
//! what reproduces it, and it is the same layout the allele maps needed.
//!
//! # A sample with no counted site poisons every row
//!
//! An empty bucket's alt fraction is `0/0`, which is NaN, and the normalizer is the **sum** of
//! every sample's fraction. One sample the traversal never counted therefore makes every mixing
//! fraction in the table NaN, including the samples that were counted perfectly well. A run
//! restricted by `-L` to one site does it, and so does a run with no `-I` at all.
//!
//! A sum of zero does the same by a different route: every fraction is `0/0` at the end, so a file
//! whose sites carry no alt read at all is NaN rather than zero.
//!
//! # Singleton is either of two tests, and the first one wins
//!
//! ```java
//! vc.isBiallelic() && vc.isSNP() && ((vc.hasAttribute(AC) && AC[0] == 1)
//!     || vc.getGenotypes().stream().filter(Genotype::isHet).count() == 1)
//! ```
//!
//! `AC=2` with one het is still counted, because the second test is reached when the first fails;
//! `AC=1` with two hets is counted as well, because the first test never looks at the genotypes.
//! And a record that passes with no het at all is then dropped, since the sample it would be
//! attributed to comes from `findFirst()` over the het genotypes.
//!
//! # The pileup walks each read to the variant's start
//!
//! `AlignmentStateMachine` is stepped until its genome position reaches the site, and the read
//! contributes only if it lands exactly there. A read whose alignment passes the position inside a
//! deletion **does** land there and contributes a deletion, which counts towards the total and
//! never towards the alt. A read that ends before the site contributes nothing.
//!
//! Vendor-failed reads are skipped by the tool itself rather than by a read filter.

use std::collections::HashMap;

use gatk_engine::alignment_state::AlignmentStateMachine;
use gatk_engine::java_hash::{hash_map_order, string_hash_code, HashOrderError};
use gatk_engine::tsv_table::{java_double_to_string, write_table};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK CalculateMixingFractions";

/// `MixingFractionTableColumn`, in the order the enum declares.
pub const COLUMNS: [&str; 2] = ["SAMPLE", "MIXING_FRACTION"];

/// The `0x200` flag the tool tests itself.
const VENDOR_QUALITY_CHECK_FAILED: u16 = 0x200;

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CalculateMixingFractionsError {
    /// The output could not be created, thrown after the whole traversal.
    CouldNotCreateOutputFile { path: String },
}

impl CalculateMixingFractionsError {
    /// A plain `UserException`, not one of its subclasses.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    /// The message, whose format string ends in a full stop after the path.
    pub fn message(&self) -> String {
        match self {
            CalculateMixingFractionsError::CouldNotCreateOutputFile { path } => {
                format!("Encountered an IO exception while trying to create output file {path}.")
            }
        }
    }
}

/// One sample's bucket: `AltAndTotalReadCounts`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AltAndTotalReadCounts {
    pub alt: i64,
    pub total: i64,
}

impl AltAndTotalReadCounts {
    /// `getAltFraction()`, which is `0/0` and therefore NaN for a bucket nothing was added to.
    pub fn alt_fraction(&self) -> f64 {
        self.alt as f64 / self.total as f64
    }
}

/// `isBiallelicSingletonHetSnp`: two tests, of which the second is only reached when the first
/// fails.
pub fn is_biallelic_singleton_het_snp(variant: &VariantContext) -> bool {
    if !is_biallelic(variant) || !is_snp(variant) {
        return false;
    }
    if let Some(allele_count) = first_allele_count(variant) {
        if allele_count == 1.0 {
            return true;
        }
    }
    variant.genotypes.iter().filter(|g| is_het(g)).count() == 1
}

/// `isBiallelic()`: the reference and exactly one alternate.
fn is_biallelic(variant: &VariantContext) -> bool {
    variant.alleles.len() == 2
}

/// `isSNP()`, which for a biallelic record is both alleles being one base and neither symbolic.
fn is_snp(variant: &VariantContext) -> bool {
    crate::remove_nearby_indels::variant_type(variant)
        == crate::remove_nearby_indels::VariantType::Snp
}

/// `getArrayAttribute(vc, AC)[0]`, read as a double the way `VariantContextGetters` does.
///
/// A decoded record carries `AC` as whatever the codec made of it, so every shape a number can
/// arrive in is read the same way here: only the **first** element is looked at, which is why an
/// `AC=1,0` on a multi-allelic record would pass a test the biallelic check has already failed.
fn first_allele_count(variant: &VariantContext) -> Option<f64> {
    let value = variant
        .attributes
        .iter()
        .find(|(key, _)| key == "AC")
        .map(|(_, value)| value)?;
    first_number(value)
}

fn first_number(value: &Value) -> Option<f64> {
    match value {
        Value::Int(number) => Some(*number as f64),
        Value::Double(number) => Some(*number),
        Value::Str(text) => text.split(',').next()?.trim().parse::<f64>().ok(),
        Value::List(values) => values.first().and_then(first_number),
        Value::Bool(_) | Value::Missing => None,
    }
}

/// `Genotype.isHet()`, which is `getType() == HET` and therefore false for a call carrying a
/// no-call allele: `./1` is MIXED rather than het, while `1/2` is het.
pub fn is_het(genotype: &Genotype) -> bool {
    let mut saw_no_call = false;
    let mut observed: Option<&htsjdk_vcf::allele::Allele> = None;
    let mut saw_multiple = false;
    for allele in &genotype.alleles {
        if allele.is_no_call() {
            saw_no_call = true;
        } else if observed.is_none() {
            observed = Some(allele);
        } else if Some(allele) != observed {
            saw_multiple = true;
        }
    }
    if saw_no_call || observed.is_none() {
        return false;
    }
    saw_multiple
}

/// The sample a counted record is attributed to: `findFirst()` over the het genotypes, so the
/// header's order decides when there is more than one.
pub fn variant_sample(variant: &VariantContext) -> Option<String> {
    variant
        .genotypes
        .iter()
        .find(|genotype| is_het(genotype))
        .map(|genotype| genotype.sample_name.clone())
}

/// The alt and total counts of one site, from the reads the walker handed to `apply`.
///
/// Each read is walked to the variant's start; a read that does not land exactly there is dropped,
/// and a read that lands there inside a deletion counts towards the total alone.
pub fn site_counts(reads: &[BamRecord], start: i32, alt_base: u8) -> AltAndTotalReadCounts {
    let mut counts = AltAndTotalReadCounts::default();
    for read in reads {
        if read.flags & VENDOR_QUALITY_CHECK_FAILED != 0 {
            continue;
        }
        let mut machine = AlignmentStateMachine::new(read);
        // The reference's loop ends on a null step, and a malformed read would have been refused
        // by the walker's filters long before this point.
        while let Ok(Some(_)) = machine.step_forward_on_genome() {
            if machine.genome_position() >= start {
                break;
            }
        }
        if machine.genome_position() != start {
            continue;
        }
        let offset = machine.read_offset();
        let base = if offset < 0 || offset as usize >= read.read_bases.len() {
            // A deletion, whose pileup element reports the deletion base rather than a read base.
            b'D'
        } else {
            read.read_bases[offset as usize]
        };
        counts.total += 1;
        if base == alt_base {
            counts.alt += 1;
        }
    }
    counts
}

/// `onTraversalSuccess`: the rows, in the order the `HashMap` iterates them.
///
/// `samples` is the header's sample list, which is the insertion order of the map.
pub fn mixing_fractions(
    samples: &[String],
    counts: &HashMap<String, AltAndTotalReadCounts>,
) -> Result<Vec<(String, f64)>, HashOrderError> {
    let normalizer: f64 = samples
        .iter()
        .map(|sample| {
            counts
                .get(sample)
                .copied()
                .unwrap_or_default()
                .alt_fraction()
        })
        .sum();
    let entries: Vec<(String, i32)> = samples
        .iter()
        .map(|sample| (sample.clone(), string_hash_code(sample)))
        .collect();
    Ok(hash_map_order(&entries)?
        .into_iter()
        .map(|sample| {
            let fraction = counts
                .get(&sample)
                .copied()
                .unwrap_or_default()
                .alt_fraction();
            (sample, fraction / normalizer)
        })
        .collect())
}

/// The whole output file.
pub fn table(rows: &[(String, f64)]) -> String {
    let rendered: Vec<Vec<String>> = rows
        .iter()
        .map(|(sample, fraction)| vec![sample.clone(), java_double_to_string(*fraction)])
        .collect();
    write_table(&COLUMNS, &rendered, &[])
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn genotype(sample: &str, call: &[usize], alleles: &[Allele]) -> Genotype {
        Genotype::new(
            sample,
            call.iter().map(|index| alleles[*index].clone()).collect(),
        )
    }

    fn variant(
        alternates: &[&str],
        allele_count: Option<&str>,
        calls: &[(&str, &[usize])],
    ) -> VariantContext {
        let mut alleles = vec![Allele::create(b"A", true).expect("a reference")];
        for alternate in alternates {
            alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
        }
        let mut context = VariantContext::new("chr1", 20, alleles);
        if let Some(count) = allele_count {
            context
                .attributes
                .push(("AC".to_string(), Value::Str(count.to_string())));
        }
        let alleles = context.alleles.clone();
        context.genotypes = calls
            .iter()
            .map(|(sample, call)| genotype(sample, call, &alleles))
            .collect();
        context
    }

    #[test]
    fn the_first_test_wins_and_the_second_is_only_reached_when_it_fails() {
        // AC=1 with one het: the first test.
        assert!(is_biallelic_singleton_het_snp(&variant(
            &["C"],
            Some("1"),
            &[("zebra", &[0, 1]), ("alpha", &[0, 0])]
        )));
        // AC=2 with one het: the first fails, the second passes.
        assert!(is_biallelic_singleton_het_snp(&variant(
            &["C"],
            Some("2"),
            &[("zebra", &[0, 0]), ("alpha", &[0, 1])]
        )));
        // AC=1 with two hets: counted, because the first test never looks at the genotypes.
        assert!(is_biallelic_singleton_het_snp(&variant(
            &["C"],
            Some("1"),
            &[("zebra", &[0, 1]), ("alpha", &[0, 1])]
        )));
        // AC=2 with two hets: neither test passes.
        assert!(!is_biallelic_singleton_het_snp(&variant(
            &["C"],
            Some("2"),
            &[("zebra", &[0, 1]), ("alpha", &[0, 1])]
        )));
        // Multi-allelic and indel are refused before either test.
        assert!(!is_biallelic_singleton_het_snp(&variant(
            &["C", "G"],
            Some("1"),
            &[("zebra", &[0, 1])]
        )));
        assert!(!is_biallelic_singleton_het_snp(&variant(
            &["ACC"],
            Some("1"),
            &[("zebra", &[0, 1])]
        )));
    }

    #[test]
    fn a_singleton_with_no_het_has_nobody_to_attribute_it_to() {
        let hom_var = variant(&["C"], Some("1"), &[("zebra", &[1, 1]), ("alpha", &[0, 0])]);
        assert!(is_biallelic_singleton_het_snp(&hom_var));
        assert_eq!(variant_sample(&hom_var), None);
    }

    #[test]
    fn the_sample_is_the_first_het_in_genotype_order() {
        let two_hets = variant(&["C"], Some("1"), &[("zebra", &[0, 1]), ("alpha", &[0, 1])]);
        assert_eq!(variant_sample(&two_hets).as_deref(), Some("zebra"));
    }

    fn counts(pairs: &[(&str, i64, i64)]) -> HashMap<String, AltAndTotalReadCounts> {
        pairs
            .iter()
            .map(|(sample, alt, total)| {
                (
                    sample.to_string(),
                    AltAndTotalReadCounts {
                        alt: *alt,
                        total: *total,
                    },
                )
            })
            .collect()
    }

    #[test]
    fn the_rows_come_out_in_hash_order_and_not_the_headers() {
        let samples: Vec<String> = ["zebra", "alpha", "mike"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let rows = mixing_fractions(
            &samples,
            &counts(&[("zebra", 2, 4), ("alpha", 1, 2), ("mike", 1, 1)]),
        )
        .expect("no treeified bucket");
        let order: Vec<&str> = rows.iter().map(|(sample, _)| sample.as_str()).collect();
        assert_eq!(order, vec!["zebra", "mike", "alpha"]);
        assert_eq!(rows[0].1, 0.25);
        assert_eq!(rows[1].1, 0.5);
        assert_eq!(rows[2].1, 0.25);
    }

    #[test]
    fn one_uncounted_sample_makes_every_row_nan() {
        let samples: Vec<String> = ["zebra", "alpha", "mike"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let rows =
            mixing_fractions(&samples, &counts(&[("zebra", 2, 4)])).expect("no treeified bucket");
        assert!(rows.iter().all(|(_, fraction)| fraction.is_nan()));
        assert_eq!(
            table(&rows),
            "SAMPLE\tMIXING_FRACTION\nzebra\tNaN\nmike\tNaN\nalpha\tNaN\n"
        );
    }

    #[test]
    fn a_sum_of_zero_is_nan_as_well() {
        let samples: Vec<String> = ["zebra"].iter().map(|s| s.to_string()).collect();
        let rows =
            mixing_fractions(&samples, &counts(&[("zebra", 0, 4)])).expect("no treeified bucket");
        assert!(rows[0].1.is_nan());
    }
}
