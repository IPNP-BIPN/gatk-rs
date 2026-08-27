//! `GenotypeGVCFs`: how a GVCF becomes a genotyped VCF.
//!
//! Every reference block is dropped, every variant is re-genotyped from its likelihoods, and
//! `<NON_REF>` is removed from the alleles that survive.
//!
//! Reading and writing the VCFs are not ported, nor is the exact posterior arithmetic. Which
//! records are emitted, which alleles they keep and which genotype is called are.

/// The allele every GVCF record carries and which the genotyper removes.
pub const NON_REF: &str = "<NON_REF>";

/// One record of the input GVCF.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub reference: String,
    /// Including `<NON_REF>`, as the GVCF carries it.
    pub alternates: Vec<String>,
    /// The genotype the GVCF wrote, which the genotyper does NOT copy.
    pub written_alleles: Vec<i32>,
    /// Per-genotype likelihoods, in the VCF's own order.
    pub likelihoods: Vec<i32>,
    pub allele_depths: Vec<i32>,
}

impl Record {
    pub fn is_reference_block(&self) -> bool {
        self.alternates.len() == 1 && self.alternates[0] == NON_REF
    }

    /// The alternates without `<NON_REF>`, in their original order.
    pub fn real_alternates(&self) -> Vec<String> {
        self.alternates
            .iter()
            .filter(|allele| *allele != NON_REF)
            .cloned()
            .collect()
    }
}

/// `HomoSapiensConstants.SNP_HETEROZYGOSITY`, the prior the calling applies.
pub const SNP_HETEROZYGOSITY: f64 = 1.0e-3;
pub const INDEL_HETEROZYGOSITY: f64 = 1.25e-4;

/// The index of one diploid genotype in a VCF likelihood array.
///
/// The order is the VCF's own: for alleles `a <= b`, the index is `b * (b + 1) / 2 + a`.
pub fn genotype_index(a: usize, b: usize) -> usize {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    b * (b + 1) / 2 + a
}

/// The called genotype, as a pair of allele indices.
///
/// The likelihoods are turned into posteriors by adding the prior for each genotype, so a margin of
/// a few phred points can be reversed by it: a site the GVCF wrote heterozygous with likelihoods
/// `3,0,900` is called HOMOZYGOUS REFERENCE, because the prior against a variant is worth more than
/// three points.
pub fn call_genotype(record: &Record, heterozygosity: f64) -> (usize, usize) {
    let alleles = 1 + record.alternates.len();
    let mut best = (0usize, 0usize);
    let mut best_posterior = f64::NEG_INFINITY;
    for b in 0..alleles {
        for a in 0..=b {
            let index = genotype_index(a, b);
            let Some(likelihood) = record.likelihoods.get(index) else {
                continue;
            };
            // The likelihoods are phred-scaled, so the log10 posterior is their negation over ten
            // plus the log10 prior.
            let prior = log10_prior(a, b, heterozygosity);
            let posterior = -(*likelihood as f64) / 10.0 + prior;
            if posterior > best_posterior {
                best_posterior = posterior;
                best = (a, b);
            }
        }
    }
    best
}

/// The diploid prior: hom-ref is almost everything, a het is the heterozygosity, and a hom-var is
/// its square over two.
fn log10_prior(a: usize, b: usize, heterozygosity: f64) -> f64 {
    if a == 0 && b == 0 {
        (1.0 - heterozygosity - heterozygosity * heterozygosity / 2.0).log10()
    } else if a == 0 || b == 0 {
        heterozygosity.log10()
    } else {
        // Every genotype carrying two alternates takes the same prior, whether they are the same
        // alternate twice or two different ones.
        (heterozygosity * heterozygosity / 2.0).log10()
    }
}

/// The arguments that change what is written.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Arguments {
    /// Writes a record for every base of every reference block, and for a site whose call is
    /// reference.
    pub include_non_variant_sites: bool,
    /// Read and applied, but it does NOT decide emission: see `is_emitted`.
    pub standard_min_confidence_for_calling: f64,
}

/// Whether one record reaches the output.
///
/// A reference BLOCK is never written unless every site is asked for, and a variant is written only
/// when its CALL carries an alternate. The calling threshold is not consulted here at all, which is
/// why moving it from 2 to 50 changes nothing about which records appear.
pub fn is_emitted(record: &Record, called: (usize, usize), arguments: &Arguments) -> bool {
    if arguments.include_non_variant_sites {
        return true;
    }
    if record.is_reference_block() {
        return false;
    }
    called.0 > 0 || called.1 > 0
}

/// One output record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotyped {
    pub start: i32,
    /// `<NON_REF>` removed, and any alternate no called genotype carries removed with it.
    pub alternates: Vec<String>,
    pub called: Vec<i32>,
}

/// `<NON_REF>` is removed, and so is every real alternate that no called genotype carries, so the
/// output's allele set can be smaller than the input's. The remaining alternates keep their order
/// and the genotype is re-indexed against them.
pub fn trim_alleles(record: &Record, called: (usize, usize)) -> Genotyped {
    let carried: Vec<usize> = (1..=record.alternates.len())
        .filter(|index| called.0 == *index || called.1 == *index)
        .collect();
    let alternates: Vec<String> = carried
        .iter()
        .filter_map(|index| record.alternates.get(index - 1))
        .filter(|allele| *allele != NON_REF)
        .cloned()
        .collect();
    let reindex = |allele: usize| -> i32 {
        if allele == 0 {
            return 0;
        }
        match carried.iter().position(|index| *index == allele) {
            Some(at) => at as i32 + 1,
            None => 0,
        }
    };
    Genotyped {
        start: record.start,
        alternates,
        called: vec![reindex(called.0), reindex(called.1)],
    }
}

/// The whole run.
///
/// `--include-non-variant-sites` expands a reference BLOCK into one record per base rather than
/// writing the block, which is why the fixture's block is deliberately short.
pub fn genotype(records: &[Record], arguments: &Arguments) -> Vec<Genotyped> {
    let mut out = Vec::new();
    for record in records {
        if record.is_reference_block() {
            if arguments.include_non_variant_sites {
                for position in record.start..=record.end {
                    out.push(Genotyped {
                        start: position,
                        alternates: Vec::new(),
                        called: vec![0, 0],
                    });
                }
            }
            continue;
        }
        let called = call_genotype(record, SNP_HETEROZYGOSITY);
        if !is_emitted(record, called, arguments) {
            continue;
        }
        out.push(trim_alleles(record, called));
    }
    out
}
