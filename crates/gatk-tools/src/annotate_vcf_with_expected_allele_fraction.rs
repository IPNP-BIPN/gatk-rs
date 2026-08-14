//! `AnnotateVcfWithExpectedAlleleFraction`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.AnnotateVcfWithExpectedAlleleFraction`
//! (GATK 4.6.2.0).
//!
//! One Float INFO field per record: the dot product of each sample's genotype weight with its
//! mixing fraction. Two of the three lines that produce it are not what they look like.
//!
//! # The weights and the fractions are in different orders
//!
//! ```java
//! mixingFractionsInSampleOrder = inputHeader.getSampleNamesInOrder().stream()
//!         .mapToDouble(mixingfractionsMap::get).toArray();
//! ...
//! final double[] weights = vc.getGenotypes().stream().mapToDouble(g -> weight(g)).toArray();
//! final double expected = MathUtils.sum(MathArrays.ebeMultiply(weights, mixingFractionsInSampleOrder));
//! ```
//!
//! `getSampleNamesInOrder()` is **sorted**, while `vc.getGenotypes()` is the VCF's **column**
//! order, and the two arrays are multiplied element by element. The golden settles it: a file whose
//! columns are `zebra, alpha, mike` with fractions `zebra=0.3, alpha=0.2, mike=0.1` annotates a
//! zebra-only het as `AF_EXP=0.100`, which is zebra's weight against **alpha's** fraction. The tool
//! computes the number it computes only when the columns happen to be in sorted order.
//!
//! A port that paired the two by sample name would produce 0.150 here: the arithmetic would be
//! right and the file would be wrong.
//!
//! # The default tool header lines never reach the file
//!
//! ```java
//! final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
//! headerLines.addAll(getDefaultToolVCFHeaderLines());
//! vcfWriter.writeHeader(vcfHeader);
//! ```
//!
//! The set is added to **after** the header was built from it, and the constructor copies. So this
//! tool's output carries no `##source=` and no `##GATKCommandLine`, where its sibling
//! [`crate::annotate_vcf_with_bam_depth`], whose two statements are the other way round, carries
//! both. The two goldens differ by exactly those lines.
//!
//! # The weight is 1.0, 0.5 or nothing
//!
//! `isHomVar` then `isHet`, with everything else falling through to zero: a no-call, a half-call
//! and a hom ref all weigh zero, and a `1/2` call weighs 0.5 like any other het.
//!
//! # Nothing checks the fractions
//!
//! They are not required to sum to one, so `AF_EXP` is not bounded by one either. A sample of the
//! VCF that the table does not name is a `NullPointerException` out of unboxing a null, and a
//! sample named twice is an `IllegalStateException` out of `Collectors.toMap`.

use htsjdk_vcf::variant::{format_vcf_double, Genotype, VariantContext};

use crate::calculate_mixing_fractions::is_het;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK AnnotateVcfWithExpectedAlleleFraction";

/// `EXPECTED_ALLELE_FRACTION_NAME`.
pub const AF_EXP: &str = "AF_EXP";

/// The INFO line the tool adds.
pub const AF_EXP_HEADER_LINE: &str =
    "##INFO=<ID=AF_EXP,Number=1,Type=Float,Description=\"expected allele fraction in pooled bam\">";

/// What this tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExpectedAlleleFractionError {
    /// A sample of the VCF that the table does not name: unboxing a null.
    SampleMissingFromTable,
    /// A sample the table names twice, out of `Collectors.toMap`.
    DuplicateSample {
        sample: String,
        first: String,
        second: String,
    },
}

impl ExpectedAlleleFractionError {
    pub fn class(&self) -> &'static str {
        match self {
            ExpectedAlleleFractionError::SampleMissingFromTable => "java.lang.NullPointerException",
            ExpectedAlleleFractionError::DuplicateSample { .. } => {
                "java.lang.IllegalStateException"
            }
        }
    }

    /// The message, which for the missing sample is the word `null`: nothing words that refusal.
    pub fn message(&self) -> String {
        match self {
            ExpectedAlleleFractionError::SampleMissingFromTable => "null".to_string(),
            ExpectedAlleleFractionError::DuplicateSample {
                sample,
                first,
                second,
            } => format!("Duplicate key {sample} (attempted merging values {first} and {second})"),
        }
    }
}

/// `weight(genotype)`: hom var, het, or nothing at all.
pub fn weight(genotype: &Genotype) -> f64 {
    if is_hom_var(genotype) {
        1.0
    } else if is_het(genotype) {
        0.5
    } else {
        0.0
    }
}

/// `isHomVar()`, which is `getType() == HOM_VAR`: every allele called, all the same, and not the
/// reference. A call carrying a no-call is MIXED and weighs nothing.
fn is_hom_var(genotype: &Genotype) -> bool {
    if genotype.alleles.is_empty() {
        return false;
    }
    if genotype.alleles.iter().any(|allele| allele.is_no_call()) {
        return false;
    }
    let first = &genotype.alleles[0];
    !first.is_reference() && genotype.alleles.iter().all(|allele| allele == first)
}

/// The mixing fractions in `getSampleNamesInOrder()` order, which is sorted.
///
/// The table itself is a `HashMap`, so the order it was written in is lost; what survives is the
/// sorted sample list of the **VCF header**, which is why a table whose rows are in another order
/// changes nothing and a VCF whose columns are in another order changes everything.
pub fn fractions_in_sample_order(
    table: &[(String, f64)],
    header_samples: &[String],
) -> Result<Vec<f64>, ExpectedAlleleFractionError> {
    for (index, (sample, value)) in table.iter().enumerate() {
        if let Some((_, earlier)) = table[..index].iter().find(|(name, _)| name == sample) {
            return Err(ExpectedAlleleFractionError::DuplicateSample {
                sample: sample.clone(),
                first: trim_zeroes(*earlier),
                second: trim_zeroes(*value),
            });
        }
    }
    let mut sorted: Vec<String> = header_samples.to_vec();
    sorted.sort();
    sorted
        .iter()
        .map(|sample| {
            table
                .iter()
                .find(|(name, _)| name == sample)
                .map(|(_, value)| *value)
                .ok_or(ExpectedAlleleFractionError::SampleMissingFromTable)
        })
        .collect()
}

/// `Double.toString` as the duplicate-key message prints it.
fn trim_zeroes(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// `MathUtils.sum(MathArrays.ebeMultiply(weights, fractions))`: element by element, left to right.
///
/// The weights are in the record's column order and the fractions in sorted sample order, and this
/// function pairs them by **position**, as the reference does. The mismatch is the behaviour.
pub fn expected_allele_fraction(weights: &[f64], fractions: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (weight, fraction) in weights.iter().zip(fractions) {
        sum += weight * fraction;
    }
    sum
}

/// The whole of `apply` for one record: the value, and the string the writer puts in the INFO
/// column.
pub fn annotation(variant: &VariantContext, fractions: &[f64]) -> String {
    let weights: Vec<f64> = variant.genotypes.iter().map(weight).collect();
    format_vcf_double(expected_allele_fraction(&weights, fractions))
}

/// The metadata lines of the output header: the input's, plus `AF_EXP`, and **not** the default
/// tool lines, which are added to the set after the header has already been built from it.
pub fn header_lines(input_lines: &[String]) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for line in input_lines {
        if !lines.contains(line) {
            lines.push(line.clone());
        }
    }
    let added = AF_EXP_HEADER_LINE.to_string();
    if !lines.contains(&added) {
        lines.push(added);
    }
    lines
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn genotype(sample: &str, call: &str, alleles: &[Allele]) -> Genotype {
        Genotype::new(
            sample,
            call.split(['/', '|'])
                .map(|index| match index.parse::<usize>() {
                    Ok(at) => alleles[at].clone(),
                    Err(_) => Allele::no_call(),
                })
                .collect(),
        )
    }

    /// The golden's own record shapes, with the columns in the file's order.
    fn record(calls: &[&str], alternates: &[&str]) -> VariantContext {
        let mut alleles = vec![Allele::create(b"A", true).expect("a reference")];
        for alternate in alternates {
            alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
        }
        let mut variant = VariantContext::new("chr1", 20, alleles.clone());
        variant.genotypes = ["zebra", "alpha", "mike"]
            .iter()
            .zip(calls)
            .map(|(sample, call)| genotype(sample, call, &alleles))
            .collect();
        variant
    }

    /// The golden's `fractions` table and header, in the orders they are written in.
    fn fractions() -> Vec<f64> {
        let table = vec![
            ("alpha".to_string(), 0.2),
            ("mike".to_string(), 0.1),
            ("zebra".to_string(), 0.3),
        ];
        let header = vec!["zebra".to_string(), "alpha".to_string(), "mike".to_string()];
        fractions_in_sample_order(&table, &header).expect("every sample is named")
    }

    #[test]
    fn the_fractions_come_out_sorted_whatever_the_columns_are() {
        assert_eq!(fractions(), vec![0.2, 0.1, 0.3]);
    }

    #[test]
    fn a_zebra_only_het_is_paired_with_alphas_fraction() {
        // The golden's first record: 0.100, not the 0.150 that pairing by name would give.
        assert_eq!(
            annotation(&record(&["0/1", "0/0", "0/0"], &["C"]), &fractions()),
            "0.100"
        );
    }

    #[test]
    fn the_weights_are_one_a_half_and_nothing() {
        let fractions = fractions();
        // Hom var and het together.
        assert_eq!(
            annotation(&record(&["1/1", "0/1", "0/0"], &["C"]), &fractions),
            "0.250"
        );
        // Every sample hom var, which is the sum of the fractions.
        assert_eq!(
            annotation(&record(&["1/1", "1/1", "1/1"], &["C"]), &fractions),
            "0.600"
        );
        // A no-call, a half-call and a hom ref, none of which weigh anything.
        assert_eq!(
            annotation(&record(&["./.", "./1", "0/0"], &["C"]), &fractions),
            "0.00"
        );
        // A 1/2 call is a het.
        assert_eq!(
            annotation(&record(&["1/2", "0/0", "0/0"], &["C", "G"]), &fractions),
            "0.100"
        );
    }

    #[test]
    fn a_sample_missing_from_the_table_is_a_null_unboxing() {
        let table = vec![("zebra".to_string(), 0.5), ("alpha".to_string(), 0.5)];
        let header = vec!["zebra".to_string(), "alpha".to_string(), "mike".to_string()];
        let refused = fractions_in_sample_order(&table, &header).expect_err("mike is missing");
        assert_eq!(refused.class(), "java.lang.NullPointerException");
        assert_eq!(refused.message(), "null");
    }

    #[test]
    fn a_sample_named_twice_is_a_duplicate_key() {
        let table = vec![
            ("zebra".to_string(), 0.5),
            ("zebra".to_string(), 0.4),
            ("alpha".to_string(), 0.05),
            ("mike".to_string(), 0.05),
        ];
        let header = vec!["zebra".to_string(), "alpha".to_string(), "mike".to_string()];
        let refused = fractions_in_sample_order(&table, &header).expect_err("zebra twice");
        assert_eq!(refused.class(), "java.lang.IllegalStateException");
        assert_eq!(
            refused.message(),
            "Duplicate key zebra (attempted merging values 0.5 and 0.4)"
        );
    }

    #[test]
    fn the_header_carries_no_source_line() {
        let input = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=200>".to_string(),
        ];
        let lines = header_lines(&input);
        assert!(lines.contains(&AF_EXP_HEADER_LINE.to_string()));
        assert!(
            !lines.iter().any(|line| line.starts_with("##source=")),
            "the default tool lines are added after the header is built"
        );
    }
}
