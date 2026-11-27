//! `PSScorer` (GATK 4.7.0.0; `PSScorer.java` is unchanged since 4.6.2.0): the third step of
//! PathSeq, with no alignment of its own.
//!
//! It reads the Kryo taxonomy database `PathSeqBuildReferenceTaxonomy` wrote, turns each read
//! template's alignments (the primary one plus the `XA` and `SA` alternates) into a set of taxon
//! ids, tags the reads with them as `YP`, and sums per-taxon scores into the scores table.
//!
//! # Every sum is taken in a hash table's order
//!
//! A score is a sum of doubles, so the order is in the bytes. A hit's taxa come out of a `HashSet`
//! into the ancestors' descendant scores, each input is its own partition with its own map, and the
//! normalization walks a `HashMap` built from `collectAsMap`'s Scala map. That last order is a Java
//! table's, except within a bucket, where it is the Scala table's: slots by the same spreading,
//! chains ascending by hash, so nothing in it depends on the order the keys arrived.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, E>` is "a `T` or the error `E`", which is how a Java `throw` is modelled; `Ok(..)`
//!   is the return and `Err(..)` the throw. A `?` after a call passes an `Err` straight up to the
//!   caller, like an exception propagating through a method that declares it.
//! - `Option<T>` is "a `T` or nothing", Java's nullable reference made explicit: `Some(x)` or
//!   `None`. `let Some(x) = y else { .. }` unwraps it or runs the block, which must leave.
//! - `JavaHashMap` reproduces `java.util.HashMap`'s iteration order (bucket by bucket, and the
//!   capacity it was built with); a `JavaHashMap<K, ()>` is a `HashSet<K>`. Rust's own `HashMap`
//!   is used only where no order is ever read out of it.
//! - `&[T]` is a read-only borrowed view of an array or list, and `&mut [T]` a writable one: the
//!   function may change the elements but the caller keeps the list.
//! - `x as i32`, `as u32` and friends are casts that truncate or wrap silently, like Java's
//!   narrowing casts; `wrapping_mul` and `wrapping_add` are Java's overflowing `int` arithmetic.
//! - `.iter().map(..).collect()` is Java's stream, map and collect; `|x| ..` is a lambda.

use std::collections::HashMap;

use gatk_engine::java_hash::JavaHashMap;
use htsjdk_bam::cigar::{CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};

use crate::pathseq_filter::matches_less_deletions_of;
use crate::pathseq_taxonomy::PsTree;

/// `PSScorer.SCORE_GENOME_LENGTH_UNITS`: scores per million reference bases.
const SCORE_GENOME_LENGTH_UNITS: f64 = 1e6;
/// `PSTaxonomyConstants.ROOT_ID`.
const ROOT_ID: i32 = 1;

/// What the scorer refuses with, and which of the tool's refusals it is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ScoreError {
    /// `UserException.BadInput`, carrying the message without its "Bad input: " prefix.
    BadInput(String),
    /// Any other exception, named by its Java class, such as the `IllegalArgumentException` of a
    /// malformed CIGAR or the `NumberFormatException` of a mismatch count.
    Other {
        class: &'static str,
        message: String,
    },
}

/// `PSTaxonomyDatabase`, read: the tree and the contig-name-to-taxon map.
///
/// The map is only ever looked up here (`valid_hits`, `missing_accessions`), so the order it was
/// filled in reaches no output of this tool.
pub struct TaxonomyDatabase {
    pub tree: PsTree,
    pub accession_to_tax_id: JavaHashMap<String, i32>,
}

/// Kryo's reads, the mirror of `gatk_engine::kryo::Output`.
///
/// A cursor over the whole file: `bytes` is the file and `at` the next byte to read. The `'a` says
/// the cursor borrows the bytes and cannot outlive them. Every read that runs off the end answers
/// `Err("Buffer underflow.")`, the message Kryo's own underflow carries.
struct Input<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Input<'_> {
    /// One byte. `ok_or` turns the `None` of a read past the end into the error, and `?` returns
    /// it.
    fn byte(&mut self) -> Result<u8, String> {
        let byte = *self.bytes.get(self.at).ok_or("Buffer underflow.")?;
        self.at += 1;
        Ok(byte)
    }

    /// `Input.readInt`: four bytes, most significant first (the order this loop shifts them in).
    /// Assembled unsigned and cast, so a high bit set becomes a negative `int` as in Java.
    fn int(&mut self) -> Result<i32, String> {
        let mut value = 0u32;
        for _ in 0..4 {
            value = (value << 8) | u32::from(self.byte()?);
        }
        Ok(value as i32)
    }

    /// `Input.readLong`: the same, over eight bytes.
    fn long(&mut self) -> Result<i64, String> {
        let mut value = 0u64;
        for _ in 0..8 {
            value = (value << 8) | u64::from(self.byte()?);
        }
        Ok(value as i64)
    }

    /// `Input.readString`: the ASCII form ends at a set high bit, the other is `count + 1` in six
    /// bits and then seven, followed by that many UTF-8 characters.
    ///
    /// A count of zero is a written `null` (hence the `Option`), and one is the empty string.
    fn string(&mut self) -> Result<Option<String>, String> {
        let first = self.byte()?;
        // A first byte with its high bit clear starts the ASCII form: bytes up to and including
        // the one whose high bit marks the end, that bit cleared.
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
                // The count is of characters, not bytes, so each lead byte says how many
                // continuation bytes follow it (`0b110x` one, `0b1110` two); the bytes are then
                // decoded in one piece.
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

    /// A number written as a string, which is how the tree writes its root and every child:
    /// `root = Integer.valueOf(input.readString());` (`PSTree.java`). A `null` or a non-number
    /// is the refusal `Integer.valueOf` would throw.
    fn integer_string(&mut self) -> Result<i32, String> {
        let text = self.string()?.ok_or("null")?;
        text.parse()
            .map_err(|_| format!("For input string: \"{text}\""))
    }
}

/// `readTaxonomyDatabase`, whose nodes only need their name, rank, parent and length here.
///
/// The field order is `PSTree`'s and `PSTreeNode`'s Kryo constructors:
///
/// ```java
/// name = input.readString();
/// rank = input.readString();
/// parent = input.readInt();
/// length = input.readLong();
/// final int numChildren = input.readInt();
/// ```
///
/// The reference reads the nodes straight into `new HashMap<>(treeSize)`; this one rebuilds the
/// tree with `add_node` and skips the children, which `add_node` re-derives from the parents.
/// Nothing in this tool iterates the tree or a node's children, only looks nodes up, so neither
/// order is needed.
pub fn read_taxonomy_database(bytes: &[u8]) -> Result<TaxonomyDatabase, String> {
    let mut input = Input { bytes, at: 0 };
    let root = input.integer_string()?;
    let mut tree = PsTree::new(root);
    let size = input.int()?;
    // All nodes are read first, because the root has to be set rather than added.
    let mut nodes = Vec::with_capacity(size.max(0) as usize);
    for _ in 0..size {
        let id = input.int()?;
        let name = input.string()?;
        let rank = input.string()?;
        let parent = input.int()?;
        let length = input.long()?;
        // The children: parsed, so a malformed one refuses as it would upstream, and dropped.
        for _ in 0..input.int()? {
            input.integer_string()?;
        }
        nodes.push((id, name, rank, parent, length));
    }
    for (id, name, rank, parent, length) in nodes {
        if id == root {
            tree.set_root(name, rank, parent, length);
        } else {
            // `as_deref().unwrap_or_default()` is the stored name, or "" for a written `null`.
            tree.add_node(
                id,
                name.as_deref().unwrap_or_default(),
                parent,
                length,
                rank.as_deref().unwrap_or_default(),
            );
        }
    }
    // `PSTaxonomyDatabase`'s map: a count, then each key and its value, both strings.
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

/// `PSPathogenTaxonScore`: one taxon's running totals. `Default` is every field zero, the Java
/// field initialisers.
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

/// A template's hits: `PSPathogenAlignmentHit`, the taxon ids in the order `mapGroupedReadsToTax`
/// listed them and the number of reads (1 or 2) they stand for.
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

/// `TextCigarCodec.decode` (htsjdk 5.0.0), for the CIGARs inside `XA` and `SA`.
///
/// `"*"` is the empty CIGAR, an operator with no digits before it is "Malformed CIGAR string",
/// and an unknown operator is `CigarOperator.characterToEnum`'s refusal, which concatenates the
/// `byte`, so it prints the number (`"Unrecognized CigarOperator: " + b`). Rust's `{other}` on a
/// `u8` prints the number too.
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
        // The reference's `length * 10 + digit` is `int` arithmetic, which wraps.
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
    // Trailing digits with no operator. The reference reads past the array here instead (its
    // inner loop indexes `cigarBytes[i]` with no bound check), so its exception differs.
    if digits {
        return Err(invalid());
    }
    Ok(elements)
}

/// `getValidAlternateHits`: the alternate alignments in an `XA` or `SA` tag that reach
/// `min_identity_bases`, appended to `hits` as (matches, truncated contig name).
///
/// `indices` is the reference's `contigIndex, cigarIndex, numMismatchesIndex`: `(0, 2, 3)` for
/// `XA`, `(0, 3, 5)` for `SA`. Each alternate is
///
/// ```java
/// final String[] tagTokens = tagValue.split(";");
/// for (final String tok : tagTokens) {
///     final String[] subtokens = tok.split(",");
/// ```
fn alternate_hits(
    record: &BamRecord,
    tag: &[u8; 2],
    indices: (usize, usize, usize),
    min_identity_bases: f64,
    hits: &mut Vec<(i32, String)>,
) -> Result<(), ScoreError> {
    let (contig, cigar, mismatches) = indices;
    let name = String::from_utf8_lossy(tag).into_owned();
    // No tag, or a tag that is not a string, contributes nothing.
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
        // `Integer.valueOf(subtokens[numMismatchesIndex])`, then the CIGAR, in that order, so a
        // bad count refuses before a bad CIGAR does.
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
///
/// The primary alignment, then the `XA` alternates, then the `SA` ones, each kept if its matches
/// (less deletions) reach `min_identity` of the read length; then only those within
/// `identity_margin` of the best survive, mapped to taxa:
///
/// ```java
/// .map(hit -> taxonomyDatabase.accessionToTaxId.getOrDefault(hit.getAccession(), null))
///         .filter(Objects::nonNull).collect(Collectors.toSet());
/// ```
///
/// That set's order is what the `YP` tag of an unpaired read prints, which is why the set is a
/// `JavaHashMap` filled in the hits' order.
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
        // `getAssignedContig()`: the reference index looked up in the dictionary. `try_from`
        // turns the `-1` of no contig into `None`, and the chain answers "" for it.
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| contigs.get(index))
            .map_or("", String::as_str);
        hits.push((matches, truncated(contig).to_string()));
    }
    alternate_hits(record, b"XA", (0, 2, 3), min_identity_bases, &mut hits)?;
    alternate_hits(record, b"SA", (0, 3, 5), min_identity_bases, &mut hits)?;
    // `max()` is `None` on no hits, the reference's `if (hits.isEmpty()) return emptySet()`.
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
///
/// A pair keeps only the taxa both mates hit, in the order of the `groupingBy` map the reference
/// counts them in:
///
/// ```java
/// taxIds.collect(Collectors.groupingBy(e -> e, Collectors.counting()));
/// ```
///
/// An unpaired read keeps its set as it came. Both reads get the same `YP`, and a template with no
/// hit gets none.
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

/// `PSTree.getPathOf`: the node, its parent, and so on up to the root.
fn path_of(tree: &PsTree, id: i32) -> Vec<i32> {
    tree.path_of(id)
}

/// `getOrAddScoreInfo`: the taxon's score, created with the tree's length on first use.
///
/// The `'a` ties the returned `&mut` to the map it points into, so the caller can change the
/// score in place, as the Java does through the object reference.
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
///
/// Per template: drop taxa absent from the tree or of length zero; credit the lowest common
/// ancestor's path with unambiguous reads; split the template's reads evenly over its taxa (over
/// genome length too, with `--divide-by-genome-length`), adding each share to the taxon's own
/// score and to every ancestor's descendant score; and credit every node on any of those paths
/// once with total reads.
///
/// The shares are summed in the `HashSet`'s order, which is why the hit's taxa go through a
/// `JavaHashMap` first:
///
/// ```java
/// final Set<Integer> hitTaxIds = new HashSet<>(hit.taxIDs);
/// ```
///
/// The result is keyed in a plain `HashMap` because its order is discarded: `reduceByKey` and then
/// `collect_as_map_order` decide the order the table is walked in.
fn compute_tax_scores(
    hits: &[AlignmentHit],
    tree: &PsTree,
    divide_by_genome_length: bool,
) -> HashMap<i32, TaxonScore> {
    let mut scores = HashMap::new();
    for hit in hits {
        // `retain` is `removeAll(hitInvalidTaxIds)`: removing keeps the survivors' order.
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
        // `hitPathNodes`, a `HashSet` upstream. It only adds the same integer to each member,
        // so its order cannot show, and a sorted set does as well.
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

/// The key order of `new HashMap<>(collectAsMap())`: Scala's `mutable.HashMap` after its
/// `sizeHint`, copied into a Java table sized for the map.
///
/// How: work out the Scala table's size, sort the keys by (slot, spread hash) to get the order
/// that table iterates in, then insert them in that order into a `JavaHashMap` built with the
/// capacity `new HashMap<>(map)` asks for, and read its keys back. Inserting in the Scala order
/// matters only where two keys share a Java bucket, since a bucket keeps its insertion order.
fn collect_as_map_order(keys: Vec<i32>) -> Vec<i32> {
    // `key ^ (key >>> 16)`: the unsigned shift is done on the `u32` view of the key.
    let spread = |key: i32| key ^ ((key as u32) >> 16) as i32;
    let size = keys.len();
    // The smallest power of two above `max(capacity - 1, 4)`, capped at 2^30. This and the
    // growth loop below are the port's reading of Scala's table sizing; there is no Scala source
    // in the reference trees to quote.
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
///
/// Three passes over the same key order: `assignKingdoms` sums each taxon's own score into its
/// kingdom's total (the first `kingdom`, `superkingdom` or root on its path, or always the root
/// with `--not-normalized-by-kingdom`); each taxon's share of that total, times 100, is added to
/// every node on its path; then one line per taxon, spaces turned to underscores:
///
/// ```java
/// printStream.println(line.replace(" ", "_"));
/// ```
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
        // Root first: the path is read leaf to root and reversed, as `Collections.reverse`.
        // A `null` name prints as "null", which is Java's string concatenation of a null.
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

/// `PSScoreMetrics`: reads with a `YP` tag and reads without.
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
///
/// ```java
/// unknownSequences.stream().forEach(file::print);
/// ```
///
/// `print`, not `println`, so the names run together with no separator and no final newline.
pub fn missing_accessions(contigs: &[String], database: &TaxonomyDatabase) -> String {
    let mut set: JavaHashMap<String, ()> = JavaHashMap::new();
    for name in contigs {
        if !database.accession_to_tax_id.contains_key(name) {
            set.insert(name.clone(), ());
        }
    }
    set.keys().map(String::as_str).collect()
}

/// `PSUtils.primaryReads`: drops secondary (0x100) and supplementary (0x800) records.
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
///
/// The paired input is taken two records at a time (`groupPairedReadsPartition`), and refuses a
/// record whose neighbour is missing or has another name; the unpaired input one at a time. The
/// reads come back paired first, as `union` orders them, with their `YP` tags set.
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
    // `if let Some(x) = y` runs the block only when the input was given.
    if let Some(mut records) = paired {
        let mut hits = Vec::new();
        // `chunks_mut(2)` hands out writable two-record windows of the list, the last one
        // shorter when the count is odd.
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
    // `reduceByKey(PSPathogenTaxonScore::add)`. With at most two partitions each sum has two
    // terms, so the order they are added in cannot change a double.
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
    // `PSScoreFileLogger.logReadCounts`: mapped is "has a `YP` tag".
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
