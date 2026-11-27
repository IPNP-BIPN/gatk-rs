//! A side feature file read alongside a sorted traversal instead of held whole.
//!
//! `featureContext.getValues(side)` at a locus is every record of the side file overlapping it,
//! in the file's own order. A traversal walks its loci in dictionary order, so a side file sorted
//! the same way can be read forward with it: only the records reaching the current locus are kept.
//! The cursor only answers what holding the whole file answers when the file IS in that order, so
//! it stops ([`Stop::Disorder`]) the moment it cannot tell, and the caller answers from the whole
//! file instead.

use std::collections::HashMap;

/// Why the cursor cannot answer.
#[derive(Debug)]
pub enum Stop<E> {
    /// The side file's reader refused.
    Input(E),
    /// The loci or the side file are not in the order the cursor relies on.
    Disorder,
}

/// One record of the side file: contig, start, stop and what the caller keeps of it.
pub type Feature<T> = (String, i64, i64, T);

/// A record whose contig is resolved to its rank in the traversal's order.
type Ranked<T> = (usize, i64, i64, T);

/// The side file read forward, a locus at a time.
pub struct ForwardOverlaps<T, E, I: Iterator<Item = Result<Feature<T>, E>>> {
    ranks: HashMap<String, usize>,
    input: I,
    /// The next record not yet taken: its rank, start, stop and payload.
    ahead: Option<Ranked<T>>,
    /// The current contig's records starting at or before the locus, in file order.
    window: Vec<(i64, i64, T)>,
    contig: Option<usize>,
    last_locus: Option<(usize, i64)>,
    last_record: Option<(usize, i64)>,
}

impl<T, E, I: Iterator<Item = Result<Feature<T>, E>>> ForwardOverlaps<T, E, I> {
    /// `contigs` in the traversal's order, which is the order the loci come in.
    pub fn new<'a>(contigs: impl IntoIterator<Item = &'a str>, input: I) -> Self {
        ForwardOverlaps {
            ranks: contigs
                .into_iter()
                .enumerate()
                .map(|(rank, name)| (name.to_string(), rank))
                .collect(),
            input,
            ahead: None,
            window: Vec::new(),
            contig: None,
            last_locus: None,
            last_record: None,
        }
    }

    /// The next record on a contig the traversal knows; a record on any other contig overlaps no
    /// locus and is passed over.
    fn pull(&mut self) -> Result<Option<Ranked<T>>, Stop<E>> {
        if let Some(ahead) = self.ahead.take() {
            return Ok(Some(ahead));
        }
        for item in self.input.by_ref() {
            let (contig, start, stop, payload) = item.map_err(Stop::Input)?;
            let Some(&rank) = self.ranks.get(&contig) else {
                continue;
            };
            if self.last_record.is_some_and(|last| (rank, start) < last) {
                return Err(Stop::Disorder);
            }
            self.last_record = Some((rank, start));
            return Ok(Some((rank, start, stop, payload)));
        }
        Ok(None)
    }

    /// The first record in file order overlapping `contig:position`. Loci must come in the order
    /// of `contigs` and of position.
    pub fn first(&mut self, contig: &str, position: i64) -> Result<Option<&T>, Stop<E>> {
        let rank = *self.ranks.get(contig).ok_or(Stop::Disorder)?;
        if self.last_locus.is_some_and(|last| (rank, position) < last) {
            return Err(Stop::Disorder);
        }
        self.last_locus = Some((rank, position));
        if self.contig != Some(rank) {
            self.window.clear();
            self.contig = Some(rank);
        }
        while let Some((at, start, stop, payload)) = self.pull()? {
            if at < rank {
                continue;
            }
            if at > rank || start > position {
                self.ahead = Some((at, start, stop, payload));
                break;
            }
            self.window.push((start, stop, payload));
        }
        // A record ending before this locus ends before every later one.
        self.window.retain(|(_, stop, _)| *stop >= position);
        Ok(self.window.first().map(|(_, _, payload)| payload))
    }

    /// The rest of the side file read through, so a refusal past the last locus still stops it.
    pub fn finish(mut self) -> Result<(), E> {
        for item in self.input.by_ref() {
            item?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The answer holding the whole file gives: the lowest file index among the records
    /// overlapping the locus.
    fn whole(records: &[Feature<usize>], contig: &str, position: i64) -> Option<usize> {
        records
            .iter()
            .filter(|(c, start, stop, _)| c == contig && *start <= position && *stop >= position)
            .map(|(_, _, _, index)| *index)
            .min()
    }

    #[test]
    fn answers_what_the_whole_file_answers_when_sorted() {
        let mut seed: u64 = 7;
        let mut next = |bound: u64| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (seed >> 33) % bound
        };
        let contigs = ["c1", "c2", "c3"];
        for _ in 0..500 {
            // A sorted file, with records on a contig the traversal does not know mixed in.
            let mut records: Vec<Feature<usize>> = Vec::new();
            for (rank, contig) in ["c1", "x", "c2", "c3"].iter().enumerate() {
                let mut start = 1;
                for _ in 0..next(12) {
                    start += next(4) as i64;
                    let stop = start + next(5) as i64 - 1;
                    records.push((contig.to_string(), start, stop, rank * 1000 + records.len()));
                }
            }
            let indexed: Vec<Feature<usize>> = records
                .iter()
                .enumerate()
                .map(|(i, (c, s, e, _))| (c.clone(), *s, *e, i))
                .collect();
            let mut cursor =
                ForwardOverlaps::new(contigs, indexed.clone().into_iter().map(Ok::<_, ()>));
            for contig in contigs {
                let mut position = 0;
                for _ in 0..next(10) {
                    position += next(6) as i64;
                    let got = cursor.first(contig, position).unwrap().copied();
                    assert_eq!(got, whole(&indexed, contig, position));
                }
            }
            assert!(cursor.finish().is_ok());
        }
    }

    #[test]
    fn stops_when_it_cannot_tell() {
        let unsorted = vec![("c1".to_string(), 10, 10, 0), ("c1".to_string(), 5, 5, 1)];
        let mut cursor = ForwardOverlaps::new(["c1"], unsorted.into_iter().map(Ok::<_, ()>));
        assert!(matches!(cursor.first("c1", 20), Err(Stop::Disorder)));

        let later = vec![("c1".to_string(), 1, 1, 0)];
        let mut cursor = ForwardOverlaps::new(["c1"], later.into_iter().map(Ok::<_, ()>));
        cursor.first("c1", 5).unwrap();
        assert!(matches!(cursor.first("c1", 4), Err(Stop::Disorder)));
        assert!(matches!(cursor.first("c9", 4), Err(Stop::Disorder)));

        let refused = vec![Ok(("c1".to_string(), 1, 1, 0)), Err("bad")];
        let cursor = ForwardOverlaps::new(["c1"], refused.into_iter());
        assert_eq!(cursor.finish(), Err("bad"));
    }
}
