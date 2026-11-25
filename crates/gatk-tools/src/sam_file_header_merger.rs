//! The header `htsjdk.samtools.SamFileHeaderMerger` builds, as `GATKSparkTool.initializeReads` asks
//! for it with `mergeDictionaries` true (htsjdk 5.0.0, GATK 4.7.0.0).
//!
//! Every input goes through it, even alone. The merged header starts from `new SAMFileHeader()`, so
//! it has no `@HD` of the inputs' but the sort order, which two inputs leave `unsorted`, and its
//! group order is `none`. Read groups and programs come out sorted by id. Two records of one id
//! that differ are a collision, and the later one takes a new id from a running counter, while the
//! reads keep the id they came with: the Spark source translates nothing.
//!
//! # The order of the merged header's lines
//!
//! The writer fixes the line order (`@HD`, `@SQ`, `@RG`, `@PG`, `@CO`); this file decides what is
//! inside each:
//!
//! - `@HD` is `VN:1.6` from the constructor, then `GO:none`, then `SO`, because the constructor
//!   calls `setGroupOrder` before `setSortOrder` and the attributes keep insertion order;
//! - `@SQ` is the first input's dictionary when every input's is the same as it, otherwise the
//!   merged one (below);
//! - `@RG` and `@PG` are sorted by their (possibly renamed) id with `String.compareTo`, so input
//!   order is lost and `A.1` sorts after `A`; the attributes of each keep the order of the first
//!   record that carried that id and those values;
//! - `@CO` is every input's comments, input after input, never deduplicated.
//!
//! # How a collision is renamed
//!
//! Two records with one id are the same record when `Map.equals` says their attributes are equal,
//! which ignores attribute order; they collapse into one, kept under the id. A record with the
//! same id and different attributes is a collision. The first set of attributes seen for an id
//! keeps the id (unless an earlier level of programs already took it), and each later one
//! becomes `recordId + "." + positiveFourDigitBase36Str(recordCounter++)`.
//!
//! The counter is **shared**: it counts every record merged so far, colliding or not, and records
//! are merged id by id, in the order each id first appears. So the suffix is not "the second A".
//! Two files that each hold one `@RG ID:A` with different samples give `A` and `A.1`; had the
//! first file also held an `@RG ID:B` before its `A`, they would give `A`, `A.2` and `B`. It is
//! written in upper-case base 36 with as many digits as it needs, despite the function's name.
//! Read groups and programs each start the counter at zero. Programs are merged level by level
//! along their `PP` chains, and a `PP` that pointed at a renamed program is rewritten to the new
//! id.
//!
//! # When two dictionaries are "the same"
//!
//! `SequenceUtil.assertSequenceDictionariesEqual` against the first input: the same number of
//! contigs, and at each index `SAMSequenceRecord.isSameSequence`. That test is looser than
//! names: a length of zero matches any length, two records that both carry `M5` are compared on
//! the checksum alone (their names are not looked at), and without `M5` a name matches its
//! counterpart's `AN` aliases. Failing it is not a refusal: the dictionaries are merged instead.
//! The merge matches contigs by name alone, so two inputs that disagree on a contig's length are
//! merged without complaint and the first input's record is the one kept. It refuses a pair of
//! contigs in opposite orders in two inputs, and in the Java also a merged list that
//! `new SAMSequenceDictionary(resultingDict)` will not hold, which this port does not check (see
//! `merge_sequences`).
//!
//! # Reading this file without reading Rust
//!
//! - `&[SamHeader]` is a borrowed, read-only view of a list; `&mut` is the writable kind, and
//!   only one may exist at a time, which is why the shared state of a merge (`taken`,
//!   `translations`, `counter`) is passed down as `&mut` rather than held in fields.
//! - `Result<T, String>` is "a `T`, or the message of the `SAMException` the reference throws";
//!   `?` after a call rethrows its error to this function's caller.
//! - `Option<T>` is a value or nothing, Java's nullable reference checked by the compiler.
//!   `if let Some(x) = y` runs its block only when there is a value.
//! - `.iter().map(..).collect()`, `.filter`, `.any`, `.all`, `.position` and `.partition` are
//!   Java stream operations; `.partition` splits one list in two by a predicate, keeping order.
//! - `|x| ...` is a closure, Java's lambda. A function's last expression is its return value.

use htsjdk_bam::header::{Attributes, ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use std::collections::{HashMap, HashSet};

/// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH`.
const UNKNOWN_SEQUENCE_LENGTH: i32 = 0;
/// `SAMProgramRecord.PREVIOUS_PROGRAM_GROUP_ID_TAG`.
const PREVIOUS_PROGRAM: &str = "PP";

/// What the merge made.
pub struct Merged {
    pub header: SamHeader,
    /// Whether the dictionaries differed, so that a read's contig indices are the merged ones only
    /// through its contig's name.
    pub merged_dictionaries: bool,
}

/// One header record and the input it came from, which `HeaderRecordAndFileHeader` holds.
///
/// The Java keys its translation tables by the header object in an `IdentityHashMap`; `file` is
/// the input's position instead, which names the same thing without holding a reference.
#[derive(Clone)]
struct Entry {
    id: String,
    attributes: Attributes,
    file: usize,
}

/// `AbstractSAMHeaderRecord.toString`, on a `SAMFileHeader`.
///
/// ```java
/// return getClass().getSimpleName() + this.mAttributes.toString();
/// ```
///
/// `SAMFileHeader` does not override it, so it is the class name and the `@HD` attributes as a
/// `LinkedHashMap` prints them, `SAMFileHeader{VN=1.6, SO=coordinate}`: no read groups, programs
/// or contigs. It is what the duplicate-id messages name an input by.
pub fn header_to_string(header: &SamHeader) -> String {
    let attributes: Vec<String> = header
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    format!("SAMFileHeader{{{}}}", attributes.join(", "))
}

/// `new SamFileHeaderMerger(identifySortOrder(headers), headers, true).getMergedHeader()`, or the
/// message of the `SAMException` it throws.
///
/// **How**: the constructor's steps in its order: dictionary, programs, read groups, group
/// order, sort order, comments. The order matters only for which error wins when an input has
/// two faults: a dictionary-order refusal comes before a duplicate `@PG` id, which comes before a
/// duplicate `@RG` id.
pub fn merge(headers: &[SamHeader]) -> Result<Merged, String> {
    // `getSequenceDictionary(headers)` and, on its `SequenceListsDifferException`,
    // `mergeSequenceDictionaries(headers)`. `match` on the `Option` takes the first arm when the
    // dictionaries are all the same and the second when they are not.
    let (sequences, merged_dictionaries) = match same_dictionary(headers) {
        Some(sequences) => (sequences, false),
        None => (merged_dictionary(headers)?, true),
    };
    // `new SAMFileHeader()` with the dictionary set: `..SamHeader::default()` fills every field
    // not named from the default header, which carries the constructor's `VN:1.6`.
    let mut header = SamHeader {
        sequences,
        ..SamHeader::default()
    };
    header.programs = merged_programs(headers)?
        .into_iter()
        .map(|(id, attributes)| ProgramRecord { id, attributes })
        .collect();
    header.read_groups = merged_read_groups(headers)?
        .into_iter()
        .map(|(id, attributes)| ReadGroup { id, attributes })
        .collect();
    header.set_group_order("none");
    // `GATKSparkTool.identifySortOrder` (GATKSparkTool.java:596), not the `ReadsPathDataSource`
    // method of the same name, which keeps a sort order all the inputs agree on. Here any second
    // input gives `unsorted`: `if (headers.size() > 1){`. `[only]` is a slice pattern that
    // matches a list of exactly one element and names it.
    //
    // The Java reads the single input's order through `getSortOrder()`, which answers `unsorted`
    // for no `SO` (as `unwrap_or` does). A non-conforming `SO` is already `unknown` by then,
    // because `SAMFileHeader.setAttribute` maps it while the header is parsed; this copies the
    // value it is given through unchanged, so it depends on the reader having done the same.
    let order = match headers {
        [only] => only.attributes.get("SO").unwrap_or("unsorted"),
        _ => "unsorted",
    };
    header.set_sort_order(order);
    // `addComment` prefixes `@CO\t` only when it is missing, and a comment read from an input
    // already has it, so copying the strings as they are is the same.
    for input in headers {
        header.comments.extend(input.comments.iter().cloned());
    }
    Ok(Merged {
        header,
        merged_dictionaries,
    })
}

/// `getSequenceDictionary`: the first dictionary, if each of the others is the same as it.
///
/// `None` stands for the `SequenceListsDifferException` the Java throws and `merge` catches; the
/// exception's own message is never seen, so it is not built. Every input is compared with the
/// FIRST, not with its neighbour, as `SequenceUtil.assertSequenceDictionariesEqual(sequences,
/// currentSequences)` is called with `sequences` fixed at the first header's.
fn same_dictionary(headers: &[SamHeader]) -> Option<Vec<SequenceRecord>> {
    // `?` on an `Option` returns `None` at once when there is no first header.
    let first = headers.first()?;
    let same = |other: &SamHeader| {
        first.sequences.len() == other.sequences.len()
            && first
                .sequences
                .iter()
                .zip(&other.sequences)
                .all(|(a, b)| is_same_sequence(a, b))
    };
    // `bool.then(|| x)` is `Some(x)` when true and `None` when false.
    headers[1..]
        .iter()
        .all(same)
        .then(|| first.sequences.clone())
}

/// `SAMSequenceRecord.isSameSequence`, for two records at one index.
///
/// The Java also compares `mSequenceIndex`, which is implied here: the caller pairs the records
/// by position. The checks run in the Java's order, and the `M5` branch returns without looking
/// at names:
///
/// ```java
/// if (this.getAttribute(SAMSequenceRecord.MD5_TAG) != null
///         && that.getAttribute(SAMSequenceRecord.MD5_TAG) != null) {
/// ```
///
/// so `chr1` and `1` with one checksum are the same sequence, and so are two equal names whose
/// checksums differ only in leading zeros or letter case.
fn is_same_sequence(a: &SequenceRecord, b: &SequenceRecord) -> bool {
    if a.length != UNKNOWN_SEQUENCE_LENGTH
        && b.length != UNKNOWN_SEQUENCE_LENGTH
        && a.length != b.length
    {
        return false;
    }
    if let (Some(a), Some(b)) = (a.attributes.get("M5"), b.attributes.get("M5")) {
        // Compared as `BigInteger`s of base sixteen: `new BigInteger(..., 16)` then `equals`.
        // Dropping leading zeros and lower-casing gives the same answer for unsigned hex text.
        // It differs where `BigInteger` does: that constructor accepts a leading `+` or `-` and
        // throws `NumberFormatException` on a non-hex checksum, which this compares as text.
        let number = |text: &str| text.trim_start_matches('0').to_ascii_lowercase();
        return number(a) == number(b);
    }
    // `getAlternativeSequenceNames()`: the `AN` attribute split on `,`. `is_some_and` is false
    // when there is no `AN`, which is the Java's empty set.
    let alternatives = |record: &SequenceRecord, name: &str| {
        record
            .attributes
            .get("AN")
            .is_some_and(|names| names.split(',').any(|alternative| alternative == name))
    };
    a.name == b.name || alternatives(a, &b.name) || alternatives(b, &a.name)
}

/// `mergeSequenceDictionaries`, from an empty dictionary.
///
/// The Java then builds `createSequenceMapping`, old index to merged index per input. Nothing
/// here needs it: a read names its contig, and the caller resolves the name in the merged
/// dictionary (`Merged::merged_dictionaries`).
fn merged_dictionary(headers: &[SamHeader]) -> Result<Vec<SequenceRecord>, String> {
    let mut merged = Vec::new();
    for header in headers {
        merged = merge_sequences(merged, &header.sequences)?;
    }
    Ok(merged)
}

/// `mergeSequences`: a sequence the dictionary lacks goes in front of the next one it has.
///
/// **How**: walk the new input's contigs. One the merged list lacks is held back; one it has
/// flushes everything held into the merged list just before it, so a new contig lands next to
/// the known contig that followed it in its own input. Contigs still held at the end are
/// appended. A known contig found at a position before the previous known one is refused:
///
/// ```java
/// } else if (prevloc > loc) {
/// ```
///
/// Matching is by name only (`getIndexOfSequenceName`), never by length or `M5`, so the record
/// already in the merged list wins and the new input's version of it is dropped.
///
/// **Not ported**: the Java returns `new SAMSequenceDictionary(resultingDict)`, whose
/// `addSequence` throws `IllegalArgumentException` on a name already in its map, and the map
/// holds every record's `AN` aliases too. Inputs `[chr1 AN:1, chr2]` then `[1]` merge to
/// `chr1, chr2, 1` here; the reference refuses the `1` as already present. No golden covers it.
fn merge_sequences(
    mut merged: Vec<SequenceRecord>,
    from: &[SequenceRecord],
) -> Result<Vec<SequenceRecord>, String> {
    let mut holder: Vec<SequenceRecord> = Vec::new();
    // `prevloc` and `previouslyMerged` in one: `None` is the Java's `-1` and `null`.
    let mut previous: Option<(usize, &str)> = None;
    for sequence in from {
        match merged.iter().position(|known| known.name == sequence.name) {
            None => holder.push(sequence.clone()),
            Some(at) => {
                if let Some((_, name)) = previous.filter(|(before, _)| *before > at) {
                    return Err(format!(
                        "Cannot merge sequence dictionaries because sequence {} and {name} are in \
                         different orders in two input sequence dictionaries.",
                        sequence.name
                    ));
                }
                // `resultingDict.addAll(loc, holder)`: `splice(at..at, ..)` inserts at `at`
                // without removing anything, and `drain(..)` empties the holder as it hands its
                // records over, which is the Java's `holder.clear()`.
                let held = holder.len();
                merged.splice(at..at, holder.drain(..));
                previous = Some((at + held, &sequence.name));
            }
        }
    }
    merged.append(&mut holder);
    Ok(merged)
}

/// Each input's records, refusing a file that holds an id twice.
///
/// The first loop of `mergeReadGroups` and of `mergeProgramGroups`, shared. A repeated id
/// within one input is a `SAMException` and not a collision; across inputs it is a collision.
/// `records` is the accessor that picks read groups or programs out of a header, `impl Fn` being
/// "any closure with this signature". The lifetime `'a` says the borrowed ids and attributes it
/// returns live as long as `headers` does.
fn entries<'a>(
    headers: &'a [SamHeader],
    kind: &str,
    records: impl Fn(&'a SamHeader) -> Vec<(&'a str, &'a Attributes)>,
) -> Result<Vec<Entry>, String> {
    let mut all = Vec::new();
    for (file, header) in headers.iter().enumerate() {
        // `idsThatAreAlreadyTaken`, cleared per input in the Java. `insert` returns false when
        // the id was already there, which is `!add(...)`.
        let mut seen = HashSet::new();
        for (id, attributes) in records(header) {
            if !seen.insert(id) {
                return Err(format!(
                    "Input file: {} contains more than one {kind} with the same id ({id})",
                    header_to_string(header)
                ));
            }
            all.push(Entry {
                id: id.to_string(),
                attributes: attributes.clone(),
                file,
            });
        }
    }
    Ok(all)
}

/// `mergeReadGroups`: one `mergeHeaderRecords` over every input's read groups, with a fresh
/// taken-set and the counter at zero, then sorted by id.
fn merged_read_groups(headers: &[SamHeader]) -> Result<Vec<(String, Attributes)>, String> {
    let groups = entries(headers, "RG", |header| {
        header
            .read_groups
            .iter()
            .map(|group| (group.id.as_str(), &group.attributes))
            .collect()
    })?;
    let mut translations = vec![HashMap::new(); headers.len()];
    // The read-group translations are built and then dropped: `Merged` carries no mapping, so
    // nothing here plays the part of `getReadGroupId`.
    let mut result = merge_records(&groups, &mut HashSet::new(), &mut translations, &mut 0);
    // `Collections.sort(result, RECORD_ID_COMPARATOR)`. Both sorts are stable. `String::cmp`
    // compares UTF-8 bytes and `compareTo` UTF-16 units; the two orders agree except between a
    // character above U+FFFF and one in U+E000..U+FFFF.
    result.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(result)
}

/// `mergeProgramGroups`: the roots first, then whatever follows a program just merged in its own
/// file, its `PP` following the id that program took.
///
/// **How**: level by level. Each round merges `current` with `mergeHeaderRecords`, applies the
/// growing per-input translation table to `current` (ids) and to `left` (ids and `PP`s), then
/// moves into the next `current` every program in `left` whose `PP` names, in the same input, a
/// program just merged. Unlike read groups, the taken-set and the counter run on across rounds,
/// so a program in a later round can be renamed against an id an earlier round took. Programs
/// whose `PP` names nothing reachable are left over and refused together, with the reference's
/// message (including its space before the newline).
fn merged_programs(headers: &[SamHeader]) -> Result<Vec<(String, Attributes)>, String> {
    let programs = entries(headers, "PG", |header| {
        header
            .programs
            .iter()
            .map(|program| (program.id.as_str(), &program.attributes))
            .collect()
    })?;
    // The roots: programs with no `PP`, in input order, as the Java's iterator-remove loop
    // leaves them. `into_iter()` hands the list's records over rather than borrowing them.
    let (mut current, mut left): (Vec<Entry>, Vec<Entry>) = programs
        .into_iter()
        .partition(|entry| entry.attributes.get(PREVIOUS_PROGRAM).is_none());
    let mut taken = HashSet::new();
    let mut translations = vec![HashMap::new(); headers.len()];
    let mut counter = 0;
    let mut overall = Vec::new();
    while !current.is_empty() {
        overall.extend(merge_records(
            &current,
            &mut taken,
            &mut translations,
            &mut counter,
        ));
        current = translate_ids(current, &translations, false);
        left = translate_ids(left, &translations, true);
        // The Java requires the same file header object and
        // `ppIdOfRecordLeftToProcess.equals(idJustProcessed)`: a `PP` only ever links within its
        // own input, compared after both sides were translated.
        let (next, rest): (Vec<Entry>, Vec<Entry>) = left.into_iter().partition(|entry| {
            let previous = entry.attributes.get(PREVIOUS_PROGRAM);
            current
                .iter()
                .any(|done| done.file == entry.file && previous == Some(done.id.as_str()))
        });
        current = next;
        left = rest;
    }
    if !left.is_empty() {
        let mut message = format!(
            "{} program groups weren't processed. Do their PP ids point to existing PGs? \n",
            left.len()
        );
        for entry in &left {
            // `getProgramName()` and `getPreviousProgramGroupId()` concatenated into a string,
            // where an absent attribute prints as Java's `null`.
            let attribute = |key| entry.attributes.get(key).unwrap_or("null");
            message.push_str(&format!(
                "@PG ID:{} PN:{} PP:{}\n",
                entry.id,
                attribute("PN"),
                attribute(PREVIOUS_PROGRAM)
            ));
        }
        return Err(message);
    }
    overall.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(overall)
}

/// `translateIds`: each record's id, and when asked its `PP`, looked up in its own input's table.
///
/// An id the table does not hold is kept, as the Java keeps the original record when
/// `translations.get(id)` is null. A translated `PP` is set in place, so it keeps its position
/// among the attributes (`setAttribute` on the `LinkedHashMap` copy), which is what the written
/// `@PG` line shows.
fn translate_ids(
    entries: Vec<Entry>,
    translations: &[HashMap<String, String>],
    translate_previous: bool,
) -> Vec<Entry> {
    entries
        .into_iter()
        .map(|mut entry| {
            let table = &translations[entry.file];
            if let Some(id) = table.get(&entry.id) {
                entry.id = id.clone();
            }
            if translate_previous {
                let previous = entry.attributes.get(PREVIOUS_PROGRAM).map(str::to_string);
                if let Some(previous) = previous.and_then(|previous| table.get(&previous)) {
                    entry.attributes.set(PREVIOUS_PROGRAM, previous);
                }
            }
            entry
        })
        .collect()
}

/// `mergeHeaderRecords`: records of one id and one set of attributes collapse, and an id already
/// taken is remapped.
///
/// **How**: the Java bins records in `idToRecord`, a `LinkedHashMap` of id to a `LinkedHashMap`
/// of record to the inputs holding it: by id in first-seen order, then by record equality in
/// first-seen order. `by_id` is the same two-level binning as lists searched linearly, which
/// keeps both first-seen orders without a hash. Then every bin, in that order, takes the id if
/// it is free and otherwise the first free `id.<base36 counter>`, and every input that held it
/// records old id to new id.
///
/// **Why the counter moves on a free id too**: the Java writes `++recordCounter` in the free
/// branch and `recordCounter++` in the colliding one, so the suffix is the number of records
/// merged before it, not a per-id ordinal (module note).
fn merge_records(
    entries: &[Entry],
    taken: &mut HashSet<String>,
    translations: &mut [HashMap<String, String>],
    counter: &mut i32,
) -> Vec<(String, Attributes)> {
    // `Map.equals`, which does not see the order of the attributes. `SAMReadGroupRecord.equals`
    // and `SAMProgramRecord.equals` are `attributesEqual(that)` plus the id, and the id is the
    // outer bin here. Sorting the pairs gives one key per set of attributes regardless of order.
    let sorted = |attributes: &Attributes| {
        let mut pairs: Vec<(String, String)> = attributes
            .iter()
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect();
        pairs.sort();
        pairs
    };
    // One inner bin: its comparison key, the attributes of the FIRST record that made it (the
    // `LinkedHashMap` key the Java keeps, so the output has that record's attribute order), and
    // the inputs that held it.
    type Same<'a> = (Vec<(String, String)>, &'a Attributes, Vec<usize>);
    let mut by_id: Vec<(&str, Vec<Same>)> = Vec::new();
    for entry in entries {
        let at = match by_id.iter().position(|(id, _)| *id == entry.id) {
            Some(at) => at,
            None => {
                by_id.push((&entry.id, Vec::new()));
                by_id.len() - 1
            }
        };
        let key = sorted(&entry.attributes);
        let same = &mut by_id[at].1;
        match same.iter_mut().find(|(known, _, _)| *known == key) {
            Some((_, _, files)) => files.push(entry.file),
            None => same.push((key, &entry.attributes, vec![entry.file])),
        }
    }
    let mut result = Vec::new();
    for (id, same) in by_id {
        for (_, attributes, files) in same {
            // `taken.insert` is the Java's `!contains` followed by `add`. `if` is an expression
            // here, and `break candidate` makes the `loop` evaluate to `candidate`.
            let new_id = if taken.insert(id.to_string()) {
                *counter += 1;
                id.to_string()
            } else {
                // A `while` loop on `idsThatAreAlreadyTaken.contains(...)` over
                // `newId = recordId + "." + positiveFourDigitBase36Str(recordCounter++)`: format
                // with the old value, then increment, and skip any suffix already taken.
                loop {
                    let candidate = format!("{id}.{}", positive_four_digit_base36(*counter));
                    *counter += 1;
                    if taken.insert(candidate.clone()) {
                        break candidate;
                    }
                }
            };
            for file in files {
                translations[file].insert(id.to_string(), new_id.clone());
            }
            result.push((new_id, attributes.clone()));
        }
    }
    result
}

/// `positiveFourDigitBase36Str`, of no fixed width despite its name.
///
/// Digits `0-9` then `A-Z`, most significant first, no padding: 35 is `Z` and 36 is `10`.
fn positive_four_digit_base36(mut left_over: i32) -> String {
    if left_over == 0 {
        return "0".to_string();
    }
    let mut digits = Vec::new();
    while left_over > 0 {
        digits.push(b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[(left_over % 36) as usize] as char);
        left_over /= 36;
    }
    digits.iter().rev().collect()
}
