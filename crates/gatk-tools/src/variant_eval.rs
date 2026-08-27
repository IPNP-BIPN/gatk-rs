//! `VariantEval`: how a call set is counted against a comparison one.
//!
//! The output is a GATKReport of one table per evaluation module, each stratified by whatever
//! stratifiers were asked for. What a row counts depends on the stratification as much as on the
//! data.
//!
//! Reading the VCFs and writing the report are not ported. Which strata a run has, which stratum a
//! record falls in, and what each module counts are.

/// One record of the evaluation set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub contig: String,
    pub position: i32,
    pub reference: String,
    pub alternates: Vec<String>,
    /// The one sample's alleles, as indices into reference-then-alternates.
    pub alleles: Vec<i32>,
}

/// `VariantContext.Type`, as the counters distinguish them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Snp,
    Mnp,
    Insertion,
    Deletion,
    Mixed,
    Symbolic,
    NoVariation,
}

impl Record {
    /// The tool's own reading of a record's type. A MULTIALLELIC site of single bases is still a
    /// SNP: it is counted ONCE, not once per alternate.
    pub fn variant_type(&self) -> VariantType {
        if self.alternates.is_empty() {
            return VariantType::NoVariation;
        }
        if self.alternates.iter().any(|allele| allele.starts_with('<')) {
            return VariantType::Symbolic;
        }
        let reference = self.reference.len();
        let same_length = self
            .alternates
            .iter()
            .all(|allele| allele.len() == reference);
        if same_length {
            return if reference == 1 {
                VariantType::Snp
            } else {
                VariantType::Mnp
            };
        }
        let all_longer = self
            .alternates
            .iter()
            .all(|allele| allele.len() > reference);
        let all_shorter = self
            .alternates
            .iter()
            .all(|allele| allele.len() < reference);
        if all_longer {
            VariantType::Insertion
        } else if all_shorter {
            VariantType::Deletion
        } else {
            VariantType::Mixed
        }
    }

    /// An indel's length, SIGNED: an insertion is positive and a deletion negative.
    pub fn indel_length(&self) -> Option<i32> {
        match self.variant_type() {
            VariantType::Insertion | VariantType::Deletion => {
                Some(self.alternates[0].len() as i32 - self.reference.len() as i32)
            }
            _ => None,
        }
    }

    pub fn is_het(&self) -> bool {
        self.alleles.len() == 2 && self.alleles[0] != self.alleles[1]
    }

    pub fn is_hom_var(&self) -> bool {
        self.alleles.len() == 2 && self.alleles[0] == self.alleles[1] && self.alleles[0] > 0
    }

    /// The base substitution, for a biallelic SNP alone.
    pub fn substitution(&self) -> Option<(u8, u8)> {
        if self.variant_type() != VariantType::Snp || self.alternates.len() != 1 {
            return None;
        }
        Some((
            self.reference.as_bytes()[0].to_ascii_uppercase(),
            self.alternates[0].as_bytes()[0].to_ascii_uppercase(),
        ))
    }

    /// A transition is a purine-to-purine or pyrimidine-to-pyrimidine change.
    pub fn is_transition(&self) -> bool {
        matches!(
            self.substitution(),
            Some((b'A', b'G')) | Some((b'G', b'A')) | Some((b'C', b'T')) | Some((b'T', b'C'))
        )
    }

    pub fn is_transversion(&self) -> bool {
        self.substitution().is_some() && !self.is_transition()
    }
}

/// The novelty stratum a record falls in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Novelty {
    Known,
    Novel,
}

/// `Novelty`'s decision, which reads the dbSNP track and NOT `--comp`.
///
/// The same file given as a comparison leaves every site novel; given as dbSNP it splits them.
/// That is the whole of the difference between the two arguments for this stratifier.
pub fn novelty(record: &Record, dbsnp_positions: &[i32]) -> Novelty {
    if dbsnp_positions.contains(&record.position) {
        Novelty::Known
    } else {
        Novelty::Novel
    }
}

/// The rows one run's table has, in the order the report writes them.
///
/// The standard stratifiers apply unless turned off, and they contribute the `all` row and the two
/// novelty ones. A stratifier MULTIPLIES the rows rather than adding a column.
pub fn strata(standard: bool, extra: &[Vec<String>]) -> Vec<Vec<String>> {
    let novelty: Vec<Vec<String>> = if standard {
        vec![
            vec!["all".to_string()],
            vec!["known".to_string()],
            vec!["novel".to_string()],
        ]
    } else {
        vec![vec!["all".to_string()]]
    };
    let mut out = novelty;
    for stratifier in extra {
        let mut multiplied = Vec::new();
        for row in &out {
            for value in stratifier {
                let mut next = row.clone();
                next.push(value.clone());
                multiplied.push(next);
            }
        }
        out = multiplied;
    }
    out
}

/// What `CountVariants` counts for one set of records.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Counts {
    pub called_loci: i32,
    pub snps: i32,
    pub mnps: i32,
    pub insertions: i32,
    pub deletions: i32,
    pub complex: i32,
    pub symbolic: i32,
    pub mixed: i32,
    pub hets: i32,
    pub hom_var: i32,
}

/// `CountVariants.update1`.
pub fn count_variants(records: &[Record]) -> Counts {
    let mut counts = Counts::default();
    for record in records {
        counts.called_loci += 1;
        match record.variant_type() {
            VariantType::Snp => counts.snps += 1,
            VariantType::Mnp => counts.mnps += 1,
            VariantType::Insertion => counts.insertions += 1,
            VariantType::Deletion => counts.deletions += 1,
            VariantType::Mixed => counts.mixed += 1,
            VariantType::Symbolic => counts.symbolic += 1,
            VariantType::NoVariation => {}
        }
        if record.is_het() {
            counts.hets += 1;
        } else if record.is_hom_var() {
            counts.hom_var += 1;
        }
    }
    counts
}

/// What `TiTvVariantEvaluator` counts: the same records under a different question.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TiTv {
    pub transitions: i32,
    pub transversions: i32,
}

impl TiTv {
    /// The ratio the report prints, which is zero rather than infinite when nothing transverted.
    pub fn ratio(&self) -> f64 {
        if self.transversions == 0 {
            0.0
        } else {
            self.transitions as f64 / self.transversions as f64
        }
    }
}

pub fn ti_tv(records: &[Record]) -> TiTv {
    let mut out = TiTv::default();
    for record in records {
        if record.is_transition() {
            out.transitions += 1;
        } else if record.is_transversion() {
            out.transversions += 1;
        }
    }
    out
}

/// What the tool refuses about its module and stratifier names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvalError {
    /// One message for BOTH namespaces: an unknown evaluator and an unknown stratifier are refused
    /// with the same wording.
    ModuleNotFound { name: String },
}

impl EvalError {
    pub fn message(&self) -> String {
        match self {
            EvalError::ModuleNotFound { name } => format!(
                "Module {name} could not be found; please check that you have specified the class \
                 name correctly"
            ),
        }
    }
}

/// The name check, which runs before a record is read.
pub fn check_module(name: &str, known: &[String]) -> Result<(), EvalError> {
    if known.iter().any(|candidate| candidate == name) {
        Ok(())
    } else {
        Err(EvalError::ModuleNotFound {
            name: name.to_string(),
        })
    }
}
