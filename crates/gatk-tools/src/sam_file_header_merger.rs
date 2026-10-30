//! The header `htsjdk.samtools.SamFileHeaderMerger` builds, as `GATKSparkTool.initializeReads` asks
//! for it with `mergeDictionaries` true.
//!
//! Every input goes through it, even alone. The merged header starts from `new SAMFileHeader()`, so
//! it has no `@HD` of the inputs' but the sort order, which two inputs leave `unsorted`, and its
//! group order is `none`. Read groups and programs come out sorted by id. Two records of one id
//! that differ are a collision, and the later one takes a new id from a running counter, while the
//! reads keep the id they came with: the Spark source translates nothing.

use htsjdk_bam::header::{Attributes, ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use std::collections::{HashMap, HashSet};

/// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH`.
const UNKNOWN_SEQUENCE_LENGTH: i32 = 0;
const PREVIOUS_PROGRAM: &str = "PP";

/// What the merge made.
pub struct Merged {
    pub header: SamHeader,
    /// Whether the dictionaries differed, so that a read's contig indices are the merged ones only
    /// through its contig's name.
    pub merged_dictionaries: bool,
}

/// One header record and the input it came from, which `HeaderRecordAndFileHeader` holds.
#[derive(Clone)]
struct Entry {
    id: String,
    attributes: Attributes,
    file: usize,
}

/// `AbstractSAMHeaderRecord.toString`, on a `SAMFileHeader`.
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
pub fn merge(headers: &[SamHeader]) -> Result<Merged, String> {
    let (sequences, merged_dictionaries) = match same_dictionary(headers) {
        Some(sequences) => (sequences, false),
        None => (merged_dictionary(headers)?, true),
    };
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
    // `identifySortOrder`.
    let order = match headers {
        [only] => only.attributes.get("SO").unwrap_or("unsorted"),
        _ => "unsorted",
    };
    header.set_sort_order(order);
    for input in headers {
        header.comments.extend(input.comments.iter().cloned());
    }
    Ok(Merged {
        header,
        merged_dictionaries,
    })
}

/// `getSequenceDictionary`: the first dictionary, if each of the others is the same as it.
fn same_dictionary(headers: &[SamHeader]) -> Option<Vec<SequenceRecord>> {
    let first = headers.first()?;
    let same = |other: &SamHeader| {
        first.sequences.len() == other.sequences.len()
            && first
                .sequences
                .iter()
                .zip(&other.sequences)
                .all(|(a, b)| is_same_sequence(a, b))
    };
    headers[1..]
        .iter()
        .all(same)
        .then(|| first.sequences.clone())
}

/// `SAMSequenceRecord.isSameSequence`, for two records at one index.
fn is_same_sequence(a: &SequenceRecord, b: &SequenceRecord) -> bool {
    if a.length != UNKNOWN_SEQUENCE_LENGTH
        && b.length != UNKNOWN_SEQUENCE_LENGTH
        && a.length != b.length
    {
        return false;
    }
    if let (Some(a), Some(b)) = (a.attributes.get("M5"), b.attributes.get("M5")) {
        // Compared as `BigInteger`s of base sixteen.
        let number = |text: &str| text.trim_start_matches('0').to_ascii_lowercase();
        return number(a) == number(b);
    }
    let alternatives = |record: &SequenceRecord, name: &str| {
        record
            .attributes
            .get("AN")
            .is_some_and(|names| names.split(',').any(|alternative| alternative == name))
    };
    a.name == b.name || alternatives(a, &b.name) || alternatives(b, &a.name)
}

/// `mergeSequenceDictionaries`, from an empty dictionary.
fn merged_dictionary(headers: &[SamHeader]) -> Result<Vec<SequenceRecord>, String> {
    let mut merged = Vec::new();
    for header in headers {
        merged = merge_sequences(merged, &header.sequences)?;
    }
    Ok(merged)
}

/// `mergeSequences`: a sequence the dictionary lacks goes in front of the next one it has.
fn merge_sequences(
    mut merged: Vec<SequenceRecord>,
    from: &[SequenceRecord],
) -> Result<Vec<SequenceRecord>, String> {
    let mut holder: Vec<SequenceRecord> = Vec::new();
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
fn entries<'a>(
    headers: &'a [SamHeader],
    kind: &str,
    records: impl Fn(&'a SamHeader) -> Vec<(&'a str, &'a Attributes)>,
) -> Result<Vec<Entry>, String> {
    let mut all = Vec::new();
    for (file, header) in headers.iter().enumerate() {
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

fn merged_read_groups(headers: &[SamHeader]) -> Result<Vec<(String, Attributes)>, String> {
    let groups = entries(headers, "RG", |header| {
        header
            .read_groups
            .iter()
            .map(|group| (group.id.as_str(), &group.attributes))
            .collect()
    })?;
    let mut translations = vec![HashMap::new(); headers.len()];
    let mut result = merge_records(&groups, &mut HashSet::new(), &mut translations, &mut 0);
    result.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(result)
}

/// `mergeProgramGroups`: the roots first, then whatever follows a program just merged in its own
/// file, its `PP` following the id that program took.
fn merged_programs(headers: &[SamHeader]) -> Result<Vec<(String, Attributes)>, String> {
    let programs = entries(headers, "PG", |header| {
        header
            .programs
            .iter()
            .map(|program| (program.id.as_str(), &program.attributes))
            .collect()
    })?;
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

/// `translateIds`.
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
fn merge_records(
    entries: &[Entry],
    taken: &mut HashSet<String>,
    translations: &mut [HashMap<String, String>],
    counter: &mut i32,
) -> Vec<(String, Attributes)> {
    // `Map.equals`, which does not see the order of the attributes.
    let sorted = |attributes: &Attributes| {
        let mut pairs: Vec<(String, String)> = attributes
            .iter()
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect();
        pairs.sort();
        pairs
    };
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
            let new_id = if taken.insert(id.to_string()) {
                *counter += 1;
                id.to_string()
            } else {
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
