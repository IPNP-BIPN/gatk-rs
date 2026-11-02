//! `PSScorer` (GATK 4.6.2.0): the third step of PathSeq, with no alignment of its own.
//!
//! # Every sum is taken in a hash table's order
//!
//! A score is a sum of doubles, so the order is in the bytes. A hit's taxa come out of a `HashSet`
//! into the ancestors' descendant scores, each input is its own partition with its own map, and the
//! normalization walks a `HashMap` built from `collectAsMap`'s Scala map. That last order is a Java
//! table's, except within a bucket, where it is the Scala table's: slots by the same spreading,
//! chains ascending by hash, so nothing in it depends on the order the keys arrived.

use std::collections::HashMap;

use gatk_engine::java_hash::JavaHashMap;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};

use crate::pathseq_filter::matches_less_deletions_of;
use crate::pathseq_taxonomy::PsTree;

/// `SCORE_GENOME_LENGTH_UNITS`.
const SCORE_GENOME_LENGTH_UNITS: f64 = 1e6;
/// `PSTaxonomyConstants.ROOT_ID`.
const ROOT_ID: i32 = 1;

/// What the scorer refuses with, and which of the tool's refusals it is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScoreError {
    BadInput(String),
    Other {
        class: &'static str,
        message: String,
    },
}

/// `PSTaxonomyDatabase`, read.
pub struct TaxonomyDatabase {
    pub tree: PsTree,
    pub accession_to_tax_id: JavaHashMap<String, i32>,
}

/// Kryo's reads, the mirror of `gatk_engine::kryo::Output`.
struct Input<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Input<'_> {
    fn byte(&mut self) -> Result<u8, String> {
        let byte = *self.bytes.get(self.at).ok_or("Buffer underflow.")?;
        self.at += 1;
        Ok(byte)
    }

    fn int(&mut self) -> Result<i32, String> {
        let mut value = 0u32;
        for _ in 0..4 {
            value = (value << 8) | u32::from(self.byte()?);
        }
        Ok(value as i32)
    }

    fn long(&mut self) -> Result<i64, String> {
        let mut value = 0u64;
        for _ in 0..8 {
            value = (value << 8) | u64::from(self.byte()?);
        }
        Ok(value as i64)
    }

    /// `Input.readString`: the ASCII form ends at a set high bit, the other is `count + 1` in six
    /// bits and then seven, followed by that many UTF-8 characters.
    fn string(&mut self) -> Result<Option<String>, String> {
        let first = self.byte()?;
        if first & 0x80 == 0 {
            let mut text = vec![first];
            loop {
                let byte = self.byte()?;
                text.push(byte & 0x7f);
                if byte & 0x80 != 0 {
                    return Ok(Some(String::from_utf8_lossy(&text).into_owned()));
                }
            }
        }
        let mut count = u64::from(first & 0x3f);
        if first & 0x40 != 0 {
            let mut shift = 6;
            loop {
                let byte = self.byte()?;
                count |= u64::from(byte & 0x7f) << shift;
                shift += 7;
                if byte & 0x80 == 0 {
                    break;
                }
            }
        }
        match count {
            0 => Ok(None),
            1 => Ok(Some(String::new())),
            count => {
                let start = self.at;
                for _ in 0..count - 1 {
                    let lead = self.byte()?;
                    let extra = match lead >> 4 {
                        12 | 13 => 1,
                        14 => 2,
                        _ => 0,
                    };
                    for _ in 0..extra {
                        self.byte()?;
                    }
                }
                Ok(Some(
                    String::from_utf8_lossy(&self.bytes[start..self.at]).into_owned(),
                ))
            }
        }
    }

    fn integer_string(&mut self) -> Result<i32, String> {
        let text = self.string()?.ok_or("null")?;
        text.parse()
            .map_err(|_| format!("For input string: \"{text}\""))
    }
}

/// `readTaxonomyDatabase`, whose nodes only need their name, rank, parent and length here.
pub fn read_taxonomy_database(bytes: &[u8]) -> Result<TaxonomyDatabase, String> {
    let mut input = Input { bytes, at: 0 };
    let root = input.integer_string()?;
    let mut tree = PsTree::new(root);
    let size = input.int()?;
    let mut nodes = Vec::with_capacity(size.max(0) as usize);
    for _ in 0..size {
        let id = input.int()?;
        let name = input.string()?;
        let rank = input.string()?;
        let parent = input.int()?;
        let length = input.long()?;
        for _ in 0..input.int()? {
            input.integer_string()?;
        }
        nodes.push((id, name, rank, parent, length));
    }
    for (id, name, rank, parent, length) in nodes {
        if id == root {
            tree.set_root(name, rank, parent, length);
        } else {
            tree.add_node(
                id,
                name.as_deref().unwrap_or_default(),
                parent,
                length,
                rank.as_deref().unwrap_or_default(),
            );
        }
    }
    let mut accession_to_tax_id = JavaHashMap::new();
    for _ in 0..input.int()? {
        let key = input.string()?.ok_or("null")?;
        accession_to_tax_id.insert(key, input.integer_string()?);
    }
    Ok(TaxonomyDatabase {
        tree,
        accession_to_tax_id,
    })
}

/// `PSPathogenTaxonScore`.
#[derive(Debug, Clone, Default)]
struct TaxonScore {
    self_score: f64,
    descendent_score: f64,
    score_normalized: f64,
    total_reads: i32,
    unambiguous_reads: i32,
    reference_length: i64,
    kingdom: i32,
}

/// A template's hits: `PSPathogenAlignmentHit`.
struct AlignmentHit {
    tax_ids: Vec<i32>,
    mates: i32,
}

/// `SAMSequenceRecord.truncateSequenceName`: up to the first whitespace.
fn truncated(name: &str) -> &str {
    let end = name
        .find([' ', '\t', '\n', '\r', '\u{c}'])
        .unwrap_or(name.len());
    &name[..end]
}

/// `TextCigarCodec.decode`.
pub fn text_cigar(text: &str) -> Result<Vec<CigarElement>, ScoreError> {
    let invalid = || ScoreError::Other {
        class: "java.lang.IllegalArgumentException",
        message: format!("Malformed CIGAR string: {text}"),
    };
    if text == "*" {
        return Ok(Vec::new());
    }
    let mut elements = Vec::new();
    let mut length = 0u32;
    let mut digits = false;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length.wrapping_mul(10).wrapping_add(u32::from(byte - b'0'));
            digits = true;
            continue;
        }
        if !digits {
            return Err(invalid());
        }
        let op = match byte {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            other => {
                return Err(ScoreError::Other {
                    class: "java.lang.IllegalArgumentException",
                    message: format!("Unrecognized CigarOperator: {other}"),
                })
            }
        };
        elements.push(CigarElement { length, op });
        length = 0;
        digits = false;
    }
    if digits {
        return Err(invalid());
    }
    Ok(elements)
}

/// `getValidAlternateHits`.
fn alternate_hits(
    record: &BamRecord,
    tag: &[u8; 2],
    indices: (usize, usize, usize),
    min_identity_bases: f64,
    hits: &mut Vec<(i32, String)>,
) -> Result<(), ScoreError> {
    let (contig, cigar, mismatches) = indices;
    let name = String::from_utf8_lossy(tag).into_owned();
    let Some(value) = record.tags.get(Tag::new(tag)) else {
        return Ok(());
    };
    let TagValue::Str(value) = value else {
        return Ok(());
    };
    let expected = contig.max(cigar).max(mismatches) + 1;
    // `String.split`, which drops trailing empty strings.
    let split = |text: &str, on: char| {
        let mut tokens: Vec<String> = text.split(on).map(str::to_string).collect();
        while tokens.len() > 1 && tokens.last().is_some_and(String::is_empty) {
            tokens.pop();
        }
        if tokens.len() == 1 && tokens[0].is_empty() && text.is_empty() {
            return tokens;
        }
        tokens
    };
    for token in split(value, ';') {
        let values = split(&token, ',');
        if values.len() < expected {
            return Err(ScoreError::BadInput(format!(
                "Error parsing {name} tag: expected at least {expected} values per alignment but \
                 found {}",
                values.len()
            )));
        }
        let text = &values[mismatches];
        let count: i32 = text.parse().map_err(|_| ScoreError::Other {
            class: "java.lang.NumberFormatException",
            message: format!("For input string: \"{text}\""),
        })?;
        let elements = text_cigar(&values[cigar])?;
        let matches = matches_less_deletions_of(&elements, count);
        if f64::from(matches) >= min_identity_bases {
            hits.push((matches, truncated(&values[contig]).to_string()));
        }
    }
    Ok(())
}

/// `getValidHits`, as the `Collectors.toSet` it returns: a `HashSet`'s order.
fn valid_hits(
    record: &BamRecord,
    contigs: &[String],
    database: &TaxonomyDatabase,
    min_identity: f64,
    identity_margin: f64,
) -> Result<JavaHashMap<i32, ()>, ScoreError> {
    let mut out = JavaHashMap::new();
    if gatk_engine::read::is_unmapped(record) {
        return Ok(out);
    }
    let Some(TagValue::Int(nm)) = record.tags.get(Tag::new(b"NM")) else {
        return Err(ScoreError::BadInput(
            "SAM flag indicates a read is mapped, but the NM tag is absent".to_string(),
        ));
    };
    let min_identity_bases = min_identity * record.read_bases.len() as f64;
    let matches = matches_less_deletions_of(&record.cigar.elements, *nm as i32);
    let mut hits = Vec::new();
    if f64::from(matches) >= min_identity_bases {
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| contigs.get(index))
            .map_or("", String::as_str);
        hits.push((matches, truncated(contig).to_string()));
    }
    alternate_hits(record, b"XA", (0, 2, 3), min_identity_bases, &mut hits)?;
    alternate_hits(record, b"SA", (0, 3, 5), min_identity_bases, &mut hits)?;
    let Some(best) = hits.iter().map(|hit| hit.0).max() else {
        return Ok(out);
    };
    let margin = (1.0 - identity_margin) * f64::from(best);
    for (matches, accession) in &hits {
        if f64::from(*matches) >= margin {
            if let Some(tax_id) = database.accession_to_tax_id.get(accession) {
                out.insert(*tax_id, ());
            }
        }
    }
    Ok(out)
}

/// `mapGroupedReadsToTax` over one template: sets the tag, and answers its hits.
fn map_template(
    template: &mut [BamRecord],
    contigs: &[String],
    database: &TaxonomyDatabase,
    min_identity: f64,
    identity_margin: f64,
) -> Result<AlignmentHit, ScoreError> {
    let mates = template.len() as i32;
    let mut sets = Vec::with_capacity(template.len());
    for record in template.iter() {
        sets.push(valid_hits(
            record,
            contigs,
            database,
            min_identity,
            identity_margin,
        )?);
    }
    let tax_ids: Vec<i32> = if mates > 1 {
        let mut counts: JavaHashMap<i32, i64> = JavaHashMap::new();
        for set in &sets {
            for id in set.keys() {
                // `groupingBy`'s `computeIfAbsent`.
                match counts.get_mut(id) {
                    Some(count) => *count += 1,
                    None => counts.insert_at_head(*id, 1),
                }
            }
        }
        counts
            .iter()
            .filter(|(_, count)| **count == i64::from(mates))
            .map(|(id, _)| *id)
            .collect()
    } else {
        sets.iter().flat_map(|set| set.keys().copied()).collect()
    };
    if !tax_ids.is_empty() {
        let joined: Vec<String> = tax_ids.iter().map(i32::to_string).collect();
        for record in template.iter_mut() {
            record
                .tags
                .insert(Tag::new(b"YP"), TagValue::Str(joined.join(",")));
        }
    }
    Ok(AlignmentHit { tax_ids, mates })
}

fn path_of(tree: &PsTree, id: i32) -> Vec<i32> {
    tree.path_of(id)
}

fn entry<'a>(
    scores: &'a mut HashMap<i32, TaxonScore>,
    tree: &PsTree,
    id: i32,
) -> &'a mut TaxonScore {
    scores.entry(id).or_insert_with(|| TaxonScore {
        reference_length: tree.length_of(id),
        ..TaxonScore::default()
    })
}

/// `computeTaxScores` over one partition.
fn compute_tax_scores(
    hits: &[AlignmentHit],
    tree: &PsTree,
    divide_by_genome_length: bool,
) -> HashMap<i32, TaxonScore> {
    let mut scores = HashMap::new();
    for hit in hits {
        let mut set: JavaHashMap<i32, ()> = JavaHashMap::copy_of(hit.tax_ids.iter());
        set.retain(|id, _| tree.has_node(*id) && tree.length_of(*id) != 0);
        let taxa: Vec<i32> = set.keys().copied().collect();
        if taxa.is_empty() {
            continue;
        }
        // `getLCA`: the lowest node of the first path that every other path has.
        let paths: Vec<Vec<i32>> = taxa.iter().map(|id| path_of(tree, *id)).collect();
        let lowest = paths[0]
            .iter()
            .copied()
            .find(|node| paths.iter().all(|path| path.contains(node)))
            .expect("the root is on every path");
        for id in path_of(tree, lowest) {
            entry(&mut scores, tree, id).unambiguous_reads += hit.mates;
        }
        let mut on_paths = std::collections::BTreeSet::new();
        for (id, path) in taxa.iter().zip(&paths) {
            let mut score = f64::from(hit.mates) / taxa.len() as f64;
            if divide_by_genome_length {
                score *= SCORE_GENOME_LENGTH_UNITS / tree.length_of(*id) as f64;
            }
            for node in path {
                on_paths.insert(*node);
                let info = entry(&mut scores, tree, *node);
                if node == id {
                    info.self_score += score;
                } else {
                    info.descendent_score += score;
                }
            }
        }
        for id in on_paths {
            entry(&mut scores, tree, id).total_reads += hit.mates;
        }
    }
    scores
}

/// The key order of `new HashMap<>(collectAsMap())`: Scala's `mutable.HashMap` after its `sizeHint`,
/// copied into a Java table sized for the map.
fn collect_as_map_order(keys: Vec<i32>) -> Vec<i32> {
    let spread = |key: i32| key ^ ((key as u32) >> 16) as i32;
    let size = keys.len();
    let table_for = |capacity: i64| -> i64 {
        let floor = (capacity - 1).max(4);
        ((1i64 << (63 - floor.leading_zeros())) * 2).min(1 << 30)
    };
    let mut table = 16i64.max(table_for(((size + 1) as f64 / 0.75) as i64));
    for count in 0..size as i64 {
        if count + 1 >= (table as f64 * 0.75) as i64 {
            table *= 2;
        }
    }
    let mut scala = keys;
    scala.sort_by_key(|key| (i64::from(spread(*key)) & (table - 1), spread(*key)));
    let mut java: JavaHashMap<i32, ()> =
        JavaHashMap::with_capacity(JavaHashMap::<i32, ()>::copy_capacity(size));
    for key in scala {
        java.insert(key, ());
    }
    java.keys().copied().collect()
}

/// `Double.toString`.
fn java_double(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// `computeNormalizedScores` and `writeScoresFile`.
fn scores_file(
    mut scores: HashMap<i32, TaxonScore>,
    tree: &PsTree,
    not_normalized_by_kingdom: bool,
) -> String {
    let order = collect_as_map_order(scores.keys().copied().collect());
    let mut sums: HashMap<i32, f64> = HashMap::new();
    for id in &order {
        let self_score = scores[id].self_score;
        let kingdom = if not_normalized_by_kingdom {
            ROOT_ID
        } else {
            path_of(tree, *id)
                .into_iter()
                .find(|node| {
                    let rank = tree.rank_of(*node).unwrap_or_default();
                    rank == "kingdom" || rank == "superkingdom" || *node == ROOT_ID
                })
                .unwrap_or(0)
        };
        *sums.entry(kingdom).or_insert(0.0) += self_score;
        scores.get_mut(id).expect("a known taxon").kingdom = kingdom;
    }
    for id in &order {
        let score = &scores[id];
        let sum = sums[&score.kingdom];
        let normalized = if sum == 0.0 {
            0.0
        } else {
            100.0 * score.self_score / sum
        };
        for node in path_of(tree, *id) {
            scores.get_mut(&node).expect("an ancestor").score_normalized += normalized;
        }
    }
    let mut text = String::from(
        "tax_id\ttaxonomy\ttype\tname\tkingdom\tscore\tscore_normalized\treads\tunambiguous\treference_length\n",
    );
    for id in &order {
        let score = &scores[id];
        let mut path: Vec<&str> = path_of(tree, *id)
            .into_iter()
            .map(|node| tree.name_of(node).unwrap_or("null"))
            .collect();
        path.reverse();
        let line = format!(
            "{id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            path.join("|"),
            tree.rank_of(*id).unwrap_or("null"),
            tree.name_of(*id).unwrap_or("null"),
            tree.name_of(score.kingdom).unwrap_or("null"),
            java_double(score.self_score + score.descendent_score),
            java_double(score.score_normalized),
            score.total_reads,
            score.unambiguous_reads,
            score.reference_length,
        );
        text.push_str(&line.replace(' ', "_"));
        text.push('\n');
    }
    text
}

/// `PSScoreMetrics`.
#[derive(Debug, Clone, Default)]
pub struct ScoreMetrics {
    pub mapped: i64,
    pub unmapped: i64,
}

impl MetricBean for ScoreMetrics {
    fn class_name(&self) -> &str {
        "org.broadinstitute.hellbender.tools.spark.pathseq.loggers.PSScoreMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        &["MAPPED_READS", "UNMAPPED_READS"]
    }

    fn values(&self) -> Vec<Value> {
        vec![Value::Long(self.mapped), Value::Long(self.unmapped)]
    }
}

/// `PSScoreFileLogger.close`'s file.
pub fn metrics_file(metrics: &ScoreMetrics, headers: &[String]) -> String {
    let mut file = MetricsFile::default();
    for line in headers {
        file.add_header(line);
    }
    file.add_metric(metrics);
    file.write()
}

/// `writeMissingReferenceAccessions`: a `HashSet`, printed with nothing between its names.
pub fn missing_accessions(contigs: &[String], database: &TaxonomyDatabase) -> String {
    let mut set: JavaHashMap<String, ()> = JavaHashMap::new();
    for name in contigs {
        if !database.accession_to_tax_id.contains_key(name) {
            set.insert(name.clone(), ());
        }
    }
    set.keys().map(String::as_str).collect()
}

/// `primaryReads`.
pub fn primary_reads(records: Vec<BamRecord>) -> Vec<BamRecord> {
    records
        .into_iter()
        .filter(|record| record.flags & (0x100 | 0x800) == 0)
        .collect()
}

/// What `scoreReads` and the logger hand back.
pub struct Scored {
    pub reads: Vec<BamRecord>,
    /// How many of `reads`, from the start, came from the paired input.
    pub paired: usize,
    pub scores: String,
    pub metrics: ScoreMetrics,
}

/// `scoreReads`: each input is its own partition, and `reduceByKey` adds the two maps.
#[allow(clippy::too_many_arguments)]
pub fn score_reads(
    paired: Option<Vec<BamRecord>>,
    unpaired: Option<Vec<BamRecord>>,
    contigs: &[String],
    database: &TaxonomyDatabase,
    min_identity: f64,
    identity_margin: f64,
    divide_by_genome_length: bool,
    not_normalized_by_kingdom: bool,
) -> Result<Scored, ScoreError> {
    let mut reads = Vec::new();
    let mut partitions = Vec::new();
    if let Some(mut records) = paired {
        let mut hits = Vec::new();
        let mut chunks = records.chunks_mut(2);
        for template in &mut chunks {
            if template.len() < 2 || template[0].read_name != template[1].read_name {
                return Err(ScoreError::BadInput(format!(
                    "Found an unpaired read but expected all reads to be paired: {}",
                    template[0].read_name
                )));
            }
            hits.push(map_template(
                template,
                contigs,
                database,
                min_identity,
                identity_margin,
            )?);
        }
        reads.extend(records);
        partitions.push(compute_tax_scores(
            &hits,
            &database.tree,
            divide_by_genome_length,
        ));
    }
    let paired_reads = reads.len();
    if let Some(mut records) = unpaired {
        let mut hits = Vec::new();
        for template in records.chunks_mut(1) {
            hits.push(map_template(
                template,
                contigs,
                database,
                min_identity,
                identity_margin,
            )?);
        }
        reads.extend(records);
        partitions.push(compute_tax_scores(
            &hits,
            &database.tree,
            divide_by_genome_length,
        ));
    }
    let mut totals: HashMap<i32, TaxonScore> = HashMap::new();
    for partition in partitions {
        for (id, score) in partition {
            match totals.get_mut(&id) {
                None => {
                    totals.insert(id, score);
                }
                Some(total) => {
                    total.self_score += score.self_score;
                    total.descendent_score += score.descendent_score;
                    total.total_reads += score.total_reads;
                    total.unambiguous_reads += score.unambiguous_reads;
                }
            }
        }
    }
    let scores = scores_file(totals, &database.tree, not_normalized_by_kingdom);
    let mapped = reads
        .iter()
        .filter(|read| read.tags.get(Tag::new(b"YP")).is_some())
        .count() as i64;
    let metrics = ScoreMetrics {
        mapped,
        unmapped: reads.len() as i64 - mapped,
    };
    Ok(Scored {
        reads,
        paired: paired_reads,
        scores,
        metrics,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_split_keeps_inner_empties_and_a_cigar_refuses_a_bare_operator() {
        assert!(text_cigar("4M1I").is_ok());
        assert!(text_cigar("M").is_err());
        assert!(text_cigar("4").is_err());
    }

    #[test]
    fn order_is_the_java_table_when_no_bucket_is_shared() {
        assert_eq!(collect_as_map_order(vec![9, 3, 17]), vec![17, 3, 9]);
    }
}
