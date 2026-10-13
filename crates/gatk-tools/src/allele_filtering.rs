//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AlleleFiltering`,
//! `AlleleFilteringHC`, `AlleleFilteringMutect` and `OccurrenceMatrix` (GATK 4.6.2.0):
//! `--flow-filter-alleles`, the haplotypes dropped before genotyping for carrying an allele too
//! weak or too strand-biased to keep.
//!
//! The alleles (events) of the haplotypes are grouped: two alleles no haplotype carries together
//! are linked when they start at most three bases apart, or at most twenty apart and inserting
//! either into the reference haplotype gives two sequences that differ in one hmer's length only.
//! Within each connected group, until nothing more goes: every allele still carried by an active
//! haplotype is scored by marginalising the read likelihoods to it against all the others (its
//! inverse): a quality (the het or hom-alt PL over the hom-inverse one, the smallest over samples;
//! for Mutect2, ten times the somatic log odds) and a strand odds ratio against all reads. The
//! weakest by quality under the threshold, then the strand-biased ones, are the candidates, and the
//! first one's haplotypes are disabled. A lone allele is removed only on the stringent quality of
//! one, unless lone alleles may go.
//!
//! # Orders the reference takes from hash tables
//!
//! The active haplotypes are a `HashSet`, so the alleles are scored in the order of the haplotypes'
//! `Arrays.hashCode` over their bases, and that order breaks ties between equal qualities. The
//! haplotypes kept are put back through a `HashSet` and a `HashMap`, so the matrix that comes out
//! has its haplotypes in that hash order, whether or not any was removed. The two alleles of a
//! marginalised matrix are a `HashMap`'s keys too, which decides which of them Mutect2's prior
//! takes as the first.

use std::collections::BTreeSet;

use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::event_map::Event;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::java_hash::{JavaHashCode, JavaHashMap};

/// The thresholds and the lone-allele switch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlleleFilterArguments {
    /// `--flow-filter-alleles-qual-threshold`.
    pub qual_threshold: f32,
    /// `--flow-filter-alleles-sor-threshold`.
    pub sor_threshold: f32,
    /// `--flow-filter-lone-alleles`.
    pub filter_lone_alleles: bool,
}

impl Default for AlleleFilterArguments {
    fn default() -> Self {
        AlleleFilterArguments {
            qual_threshold: 30.0,
            sor_threshold: 3.0,
            filter_lone_alleles: false,
        }
    }
}

/// How an allele's quality against its inverse is scored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterModel {
    /// `AlleleFilteringHC`: the samples' genotype PLs at this ploidy.
    Germline { ploidy: usize },
    /// `AlleleFilteringMutect`: the somatic log odds, with these prior pseudocounts.
    Somatic {
        ref_pseudocount: f64,
        alt_pseudocount: f64,
    },
}

/// What filtering leaves: the matrix over the haplotypes kept (in their hash order), the starts of
/// the alleles that were weighed where a candidate was found, and how many haplotypes went.
pub struct Filtered {
    pub likelihoods: AlleleLikelihoods<BamRecord, Haplotype>,
    pub suspicious_locations: BTreeSet<i32>,
    pub filtered_haplotype_count: i32,
}

/// A haplotype as a hash key: its index in the matrix, hashed as `Haplotype.hashCode`.
#[derive(Debug, Clone, PartialEq)]
struct HaplotypeKey {
    index: usize,
    hash: i32,
}

impl JavaHashCode for HaplotypeKey {
    fn java_hash_code(&self) -> i32 {
        self.hash
    }
}

/// `filterAlleles(readLikelihoods, activeWindowStart, suspiciousLocations)`.
///
/// The haplotypes' event maps are built first, over `reference` at `ref_loc`, as
/// `buildEventMapsForHaplotypes` does just before.
pub fn filter_alleles(
    likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
    reference: &[u8],
    ref_loc: &gatk_engine::interval::SimpleInterval,
    max_mnp_distance: i32,
    arguments: &AlleleFilterArguments,
    model: FilterModel,
) -> Result<Filtered, String> {
    let haplotypes: Vec<Haplotype> = (0..likelihoods.number_of_alleles())
        .map(|a| likelihoods.get_allele(a).expect("an allele").clone())
        .collect();
    let mut events: Vec<Vec<Event>> = Vec::with_capacity(haplotypes.len());
    for h in &haplotypes {
        let map = gatk_engine::event_map::EventMap::from_haplotype(
            h,
            reference,
            ref_loc,
            max_mnp_distance,
        )
        .map_err(|e| e.message())?;
        events.push(map.events().cloned().collect());
    }
    let keys: Vec<HaplotypeKey> = haplotypes
        .iter()
        .enumerate()
        .map(|(index, h)| HaplotypeKey {
            index,
            hash: h.java_hash_code(),
        })
        .collect();

    // `haplotypeAlleleMap`, a `HashMap` filled in the matrix's order, is the occurrence matrix's
    // row order; its columns are the alleles in order of first sighting along those rows.
    let mut haplotype_allele_map: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::new();
    for key in &keys {
        haplotype_allele_map.insert(key.clone(), ());
    }
    haplotype_allele_map
        .check()
        .map_err(|e| format!("port limitation: {e:?}"))?;
    let rows: Vec<usize> = haplotype_allele_map.keys().map(|k| k.index).collect();
    let mut columns: Vec<Event> = Vec::new();
    for &row in &rows {
        for event in &events[row] {
            if !columns.contains(event) {
                columns.push(event.clone());
            }
        }
    }
    let occurs = |row: usize, column: usize| events[row].contains(&columns[column]);
    let mut non_co_occurring: Vec<(usize, usize)> = Vec::new();
    for i in 0..columns.len() {
        for j in i + 1..columns.len() {
            if !rows.iter().any(|&r| occurs(r, i) && occurs(r, j)) {
                non_co_occurring.push((i, j));
            }
        }
    }
    let distance = |(i, j): &(usize, usize)| (columns[*i].start() - columns[*j].start()).abs();
    let close: Vec<(usize, usize)> = non_co_occurring
        .iter()
        .copied()
        .filter(|pair| distance(pair) <= 3)
        .collect();
    let reference = haplotypes.iter().find(|h| h.is_reference());
    let mut linked: Vec<(usize, usize)> = Vec::new();
    for pair in non_co_occurring.iter().filter(|pair| distance(pair) <= 20) {
        let reference = reference.ok_or("java.lang.NullPointerException")?;
        let insert = |event: &Event| {
            insert_allele(
                reference,
                event.ref_allele(),
                event.alt_allele(),
                event.start(),
            )
            .ok_or("java.lang.NullPointerException".to_string())
        };
        let left = insert(&columns[pair.0])?;
        let right = insert(&columns[pair.1])?;
        if equal_up_to_hmer_change(&left, &right) {
            linked.push(*pair);
        }
    }
    linked.extend(close);
    let sets = connected_sets(columns.len(), &linked);

    let mut disabled: Vec<bool> = vec![false; haplotypes.len()];
    let mut suspicious = BTreeSet::new();
    for set in &sets {
        let in_set = |event: &Event| set.iter().any(|&c| columns[c] == *event);
        let mut active: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::copy_of(keys.iter());
        loop {
            let mut active_alleles: Vec<Event> = Vec::new();
            for key in active.keys() {
                for event in events[key.index].iter().filter(|e| in_set(e)) {
                    if !active_alleles.contains(event) {
                        active_alleles.push(event.clone());
                    }
                }
            }
            let is_active = |index: usize| active.keys().any(|k| k.index == index);
            let carriers = |event: &Event| -> Vec<usize> {
                (0..haplotypes.len())
                    .filter(|&h| is_active(h) && events[h].contains(event))
                    .collect()
            };
            let mut rpls = Vec::with_capacity(active_alleles.len());
            let mut sors = Vec::with_capacity(active_alleles.len());
            for event in &active_alleles {
                let enabled: Vec<usize> = (0..haplotypes.len()).filter(|&h| is_active(h)).collect();
                let marginal = Marginal::new(likelihoods, &events, &enabled, event)?;
                rpls.push(marginal.quality(model)?);
                sors.push(marginal.sor()?);
            }
            let candidates = identify_bad_alleles(
                &rpls,
                &sors,
                f64::from(arguments.qual_threshold),
                f64::from(arguments.sor_threshold),
            );
            let stringent = identify_bad_alleles(&rpls, &sors, 1.0, f64::from(i32::MAX));
            if !candidates.is_empty() && !set.is_empty() {
                for event in &active_alleles {
                    suspicious.insert(event.start());
                }
            }
            let remove = (!candidates.is_empty() && active_alleles.len() > 1)
                || (active_alleles.len() == 1 && !stringent.is_empty())
                || (!candidates.is_empty() && arguments.filter_lone_alleles);
            if !remove {
                break;
            }
            if !stringent.is_empty() && candidates.is_empty() {
                return Err("org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException: The thresholds for stringent allele filtering should always be higher than for the relaxed one".to_string());
            }
            let candidate = &active_alleles[candidates[0]];
            for h in carriers(candidate) {
                disabled[h] = true;
                active.remove(&keys[h]);
            }
        }
    }

    // `eventualAlleles`, a `HashSet`, then `removeAllelesToSubset`'s `HashMap`: the haplotypes kept
    // come out in that map's order.
    let mut eventual: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::new();
    for key in keys.iter().filter(|k| !disabled[k.index]) {
        eventual.insert(key.clone(), ());
    }
    let mut subset: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::new();
    for key in eventual.keys() {
        subset.insert(key.clone(), ());
    }
    subset
        .check()
        .map_err(|e| format!("port limitation: {e:?}"))?;
    let new_to_old: Vec<(Haplotype, Vec<Haplotype>)> = subset
        .keys()
        .map(|k| {
            (
                haplotypes[k.index].clone(),
                vec![haplotypes[k.index].clone()],
            )
        })
        .collect();
    let kept = likelihoods
        .marginalize(&new_to_old)
        .map_err(|e| format!("{e:?}"))?;
    let filtered_haplotype_count = (haplotypes.len() - new_to_old.len()) as i32;
    Ok(Filtered {
        likelihoods: kept,
        suspicious_locations: suspicious,
        filtered_haplotype_count,
    })
}

/// One allele against its inverse: the matrix marginalised to the two, its alleles in the order of
/// the reference's `HashMap` over them.
struct Marginal {
    likelihoods: AlleleLikelihoods<BamRecord, Allele>,
    /// Which of the two is the allele, and which its inverse.
    allele_index: usize,
    inverse_index: usize,
}

impl Marginal {
    fn new(
        likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
        events: &[Vec<Event>],
        enabled: &[usize],
        event: &Event,
    ) -> Result<Marginal, String> {
        let alt = event.alt_allele().clone();
        let alt_bases = alt.display_string().into_bytes();
        // `InverseAllele.of(alt, true)`: a reference allele, unequal to the alternate, hashed over
        // its own and the alternate's.
        let inverse = Allele::create(&alt_bases, true).map_err(|e| format!("{e:?}"))?;
        let alt_hash = allele_hash(&alt_bases, false);
        let inverse_hash = allele_hash(&alt_bases, true)
            .wrapping_mul(31)
            .wrapping_add(alt_hash);
        let carrying: Vec<Haplotype> = enabled
            .iter()
            .filter(|&&h| events[h].contains(event))
            .map(|&h| likelihoods.get_allele(h).expect("an allele").clone())
            .collect();
        let others: Vec<Haplotype> = enabled
            .iter()
            .filter(|&&h| !events[h].contains(event))
            .map(|&h| likelihoods.get_allele(h).expect("an allele").clone())
            .collect();
        let order =
            gatk_engine::java_hash::hash_map_order(&[(0usize, alt_hash), (1usize, inverse_hash)])
                .map_err(|e| format!("port limitation: {e:?}"))?;
        let mut new_to_old: Vec<(Allele, Vec<Haplotype>)> = Vec::new();
        for which in &order {
            if *which == 0 {
                new_to_old.push((alt.clone(), carrying.clone()));
            } else {
                new_to_old.push((inverse.clone(), others.clone()));
            }
        }
        let allele_index = order.iter().position(|w| *w == 0).expect("present");
        let inverse_index = 1 - allele_index;
        Ok(Marginal {
            likelihoods: likelihoods
                .marginalize(&new_to_old)
                .map_err(|e| format!("{e:?}"))?,
            allele_index,
            inverse_index,
        })
    }

    /// `getAlleleLikelihoodVsInverse`.
    fn quality(&self, model: FilterModel) -> Result<i32, String> {
        match model {
            FilterModel::Germline { ploidy } => {
                let mut smallest: Option<i32> = None;
                for sample in 0..self.likelihoods.number_of_samples() {
                    let count = self.likelihoods.sample_evidence_count(sample);
                    let row = |a: usize| -> Vec<f64> {
                        (0..count)
                            .map(|r| self.likelihoods.value(sample, a, r))
                            .collect()
                    };
                    // The genotyping allele list is (inverse, allele).
                    let rows = vec![row(self.inverse_index), row(self.allele_index)];
                    let gls = gatk_engine::genotype_likelihood_calculator::
                        compute_log10_genotype_likelihoods(ploidy, &rows);
                    let pls = gatk_engine::genotype_likelihood_calculator::gls_to_pls(&gls);
                    let het = pls[1].wrapping_sub(pls[0]);
                    let hom = *pls.get(2).ok_or_else(|| {
                        "java.lang.ArrayIndexOutOfBoundsException: Index 2 out of bounds for length 2"
                            .to_string()
                    })?;
                    let value = het.min(hom.wrapping_sub(pls[0]));
                    smallest = Some(smallest.map_or(value, |s| s.min(value)));
                }
                smallest.ok_or_else(|| "java.util.NoSuchElementException".to_string())
            }
            FilterModel::Somatic {
                ref_pseudocount,
                alt_pseudocount,
            } => {
                // `combinedLikelihoodMatrix`: every sample's evidence in one, in the two alleles'
                // order.
                let mut matrix: Vec<Vec<f64>> = vec![Vec::new(), Vec::new()];
                for sample in 0..self.likelihoods.number_of_samples() {
                    for r in 0..self.likelihoods.sample_evidence_count(sample) {
                        for (a, row) in matrix.iter_mut().enumerate() {
                            row.push(self.likelihoods.value(sample, a, r));
                        }
                    }
                }
                let pseudocounts = |n: usize| -> Vec<f64> {
                    (0..n)
                        .map(|i| {
                            if i == 0 {
                                ref_pseudocount
                            } else {
                                alt_pseudocount
                            }
                        })
                        .collect()
                };
                let empty = matrix[0].is_empty();
                let with_all = if empty {
                    0.0
                } else {
                    crate::mutect2::log_evidence(&matrix, &pseudocounts(2))
                        .map_err(|e| format!("{e:?}"))?
                };
                let without: Vec<Vec<f64>> = vec![matrix[self.inverse_index].clone()];
                let with_out = if empty {
                    0.0
                } else {
                    crate::mutect2::log_evidence(&without, &pseudocounts(1))
                        .map_err(|e| format!("{e:?}"))?
                };
                Ok((10.0 * (-with_all + with_out)) as i32)
            }
        }
    }

    /// `getAlleleSOR`: `getContingencyTableWrtAll` against the inverse, minimum count one.
    fn sor(&self) -> Result<f64, String> {
        let mut table = [[0i32; 2]; 2];
        for sample in 0..self.likelihoods.number_of_samples() {
            let mut sample_table = [0i32; 4];
            let evidence = self.likelihoods.sample_evidence(sample).unwrap_or(&[]);
            for best in self
                .likelihoods
                .best_alleles_breaking_ties_for_sample(sample, None)
            {
                if !best.is_informative() {
                    continue;
                }
                let Some(allele) = &best.allele else { continue };
                let read = &evidence[best.evidence_index];
                let index = self.likelihoods.index_of_allele(allele);
                let offset = if index == Some(self.inverse_index) {
                    0
                } else if index == Some(self.allele_index) {
                    2
                } else {
                    continue;
                };
                let forward = read.flags & 0x10 == 0;
                sample_table[offset + usize::from(!forward)] += 1;
            }
            if sample_table.iter().sum::<i32>() > 1 {
                table[0][0] += sample_table[0];
                table[0][1] += sample_table[1];
                table[1][0] += sample_table[2];
                table[1][1] += sample_table[3];
            }
        }
        for column in 0..2 {
            table[0][column] += table[1][column];
        }
        Ok(gatk_annotation::strand_bias::calculate_sor(table))
    }
}

/// `SimpleAllele.hashCode`: `Arrays.hashCode(bases) * 31 + Boolean.hashCode(isRef)`.
fn allele_hash(bases: &[u8], is_ref: bool) -> i32 {
    let mut result: i32 = 1;
    for &base in bases {
        result = result.wrapping_mul(31).wrapping_add(i32::from(base as i8));
    }
    result
        .wrapping_mul(31)
        .wrapping_add(if is_ref { 1231 } else { 1237 })
}

/// `identifyBadAlleles`: indices into the alleles, the weakest by quality first while above the
/// negated threshold, then those over the strand odds ratio threshold in the same order.
fn identify_bad_alleles(
    rpls: &[i32],
    sors: &[f64],
    qual_threshold: f64,
    sor_threshold: f64,
) -> Vec<usize> {
    let mut sorted: Vec<usize> = (0..rpls.len()).collect();
    sorted.sort_by_key(|&i| rpls[i]);
    let threshold = -qual_threshold;
    let mut result = Vec::new();
    for &i in sorted.iter().rev() {
        if f64::from(rpls[i]) <= threshold {
            break;
        }
        result.push(i);
    }
    for &i in sorted.iter().rev() {
        if sors[i] > sor_threshold && !result.contains(&i) {
            result.push(i);
        }
    }
    result
}

/// `ConnectivityInspector.connectedSets()` over the columns: each component, found from the
/// columns in order.
fn connected_sets(columns: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); columns];
    for &(a, b) in edges {
        if a != b {
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
    }
    let mut seen = vec![false; columns];
    let mut sets = Vec::new();
    for start in 0..columns {
        if seen[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = std::collections::VecDeque::from([start]);
        seen[start] = true;
        while let Some(v) = queue.pop_front() {
            component.push(v);
            for &w in &adjacency[v] {
                if !seen[w] {
                    seen[w] = true;
                    queue.push_back(w);
                }
            }
        }
        sets.push(component);
    }
    sets
}

/// `Haplotype.insertAllele(ref, alt, location)`: the haplotype's bases with the allele put in at
/// the reference location, or `None` where the reference returns null.
fn insert_allele(
    haplotype: &Haplotype,
    reference: &Allele,
    alt: &Allele,
    location: i32,
) -> Option<Vec<u8>> {
    let bases = haplotype.bases();
    let start = haplotype.genome_location()?.start;
    if reference.len() > alt.len() && location == start - 1 {
        let deletion = reference.len() - alt.len();
        if deletion > bases.len() {
            return None;
        }
        return Some(bases[deletion..].to_vec());
    }
    let cigar = haplotype.cigar()?;
    let (index, op) =
        gatk_engine::read_utils::read_index_for_reference_coordinate(start, cigar, location);
    if index == gatk_engine::read_utils::READ_INDEX_NOT_FOUND
        || !op.is_some_and(|op| op.consumes_read_bases())
    {
        return None;
    }
    let index = index as usize;
    if index + reference.len() > bases.len() {
        return None;
    }
    let mut out = bases[..index].to_vec();
    out.extend_from_slice(alt.display_string().as_bytes());
    out.extend_from_slice(&bases[index + reference.len()..]);
    Some(out)
}

/// `BaseUtils.equalUpToHmerChange`: the same runs of the same bases, at most one of a different
/// length.
fn equal_up_to_hmer_change(first: &[u8], second: &[u8]) -> bool {
    let runs = |bases: &[u8]| -> Vec<(u8, usize)> {
        let mut out: Vec<(u8, usize)> = Vec::new();
        for &b in bases {
            match out.last_mut() {
                Some((base, length)) if *base == b => *length += 1,
                _ => out.push((b, 1)),
            }
        }
        out
    };
    let (r1, r2) = (runs(first), runs(second));
    let mut found = false;
    let shared = r1.len().min(r2.len());
    for i in 0..shared {
        if r1[i].0 != r2[i].0 {
            return false;
        }
        if r1[i].1 == r2[i].1 {
            continue;
        }
        if found {
            return false;
        }
        found = true;
    }
    r1.len() == r2.len()
}
