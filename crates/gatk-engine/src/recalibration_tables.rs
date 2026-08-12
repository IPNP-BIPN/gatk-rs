//! `NestedIntegerArray` and `RecalibrationTables`, ported from
//! `org.broadinstitute.hellbender.utils.collections` and
//! `org.broadinstitute.hellbender.utils.recalibration` (GATK 4.6.2.0).
//!
//! The container a recalibration table is: four sparse arrays of [`RecalDatum`] indexed by covariate
//! key. Two are "special", because the report writes them differently, and every shape comes from
//! the covariates' own `maximumKeyValue` rather than from a constant here.
//!
//! | table | shape |
//! |---|---|
//! | read group | `numReadGroups` x 3 |
//! | quality score | `numReadGroups` x 94 x 3 |
//! | one per additional covariate | `numReadGroups` x 94 x (`maximumKeyValue` + 1) x 3 |
//!
//! # Which lookup answers nothing and which one fails
//!
//! The reference has five ways to read a value and they do not agree.
//!
//! ```java
//! public T get(final int... keys) {
//!     final int numNestedDimensions = numDimensions - 1;
//!     for( int i = 0; i < numNestedDimensions; i++ ) {
//!         if ( keys[i] >= dimensions[i] ) return null;
//!         ...
//!     }
//!     return (T)myData[keys[numNestedDimensions]];   // <- no check
//! }
//! ```
//!
//! The loop stops one short of the last dimension, so the **last key is never bounds-checked** and
//! an out-of-range one is an `ArrayIndexOutOfBoundsException` where every other out-of-range key is
//! a null. The specialised `get2Keys`, `get3Keys` and `get4Keys` do check it, so
//! `get(1, 2, 5)` throws where `get3Keys(1, 2, 5)` answers null. They are documented as a
//! performance specialisation of the same function, and they are not the same function.
//!
//! `get1Key` checks nothing at all: its comment says the bounds check is done in the caller, and no
//! caller does it. And a **negative key** is refused by nothing anywhere, because every test is
//! `>=` against the dimension, so it reaches the index and throws in all five.
//!
//! This port keeps every one of those distinctions, because a table read the wrong way is a
//! recalibration that silently uses the wrong datum.
//!
//! # Combining shares objects, and `safeCombine` is not safe
//!
//! `combineTables` walks the right table's leaves and, where the left table has nothing, **stores
//! the right table's object itself**. The two tables then share that datum. `safeCombine` allocates
//! a new set of tables and combines both arguments into it, so the first combine moves the left
//! table's datums into the new one and the second combine mutates them: after
//! `safeCombine(one, two)` the new table's datum is the same object as `one`'s and holds the sum of
//! both. The port reproduces this with `Rc<RefCell<RecalDatum>>`, because copying would leave the
//! caller holding something the reference does not.

use std::cell::RefCell;
use std::rc::Rc;

use crate::covariates::{CovariateKind, StandardCovariateList};
use crate::recal_datum::{EventType, RecalDatum, RecalDatumError};

/// The reference's `NUM_DIMENSIONS_TO_PREALLOCATE`.
///
/// It decides nothing about the values a table can hold; it decides only which branches exist before
/// anything is written. It is kept because `getAllLeaves` walks the tree and a differently
/// preallocated one would walk a different shape, even though every leaf is null either way.
const NUM_DIMENSIONS_TO_PREALLOCATE: usize = 2;

/// A datum in a table, shared rather than copied. See the module note on `safeCombine`.
pub type SharedDatum = Rc<RefCell<RecalDatum>>;

/// What the array refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NestedArrayError {
    /// `Utils.validateArg` in the constructor.
    NoDimensions,
    /// `put` with the wrong number of keys.
    WrongKeyCount { expected: usize, provided: usize },
    /// `put` with a nested key past its dimension. The maximum in the message is the dimension
    /// minus one, which is what the reference prints.
    KeyTooLarge {
        key: i32,
        dimension: usize,
        maximum: i32,
    },
    /// Where the reference throws `ArrayIndexOutOfBoundsException`: an unchecked last key, or any
    /// negative key.
    IndexOutOfBounds { index: i32, length: usize },
    /// `combineTables` on two tables of different shapes.
    DifferentShapes { left: String, right: String },
}

impl NestedArrayError {
    pub fn message(&self) -> String {
        match self {
            NestedArrayError::NoDimensions => {
                "There must be at least one dimension to an NestedIntegerArray".to_string()
            }
            NestedArrayError::WrongKeyCount { expected, provided } => format!(
                "Exactly {expected} keys should be passed to this NestedIntegerArray but {provided} were provided"
            ),
            NestedArrayError::KeyTooLarge {
                key,
                dimension,
                maximum,
            } => format!("Key {key} is too large for dimension {dimension} (max is {maximum})"),
            NestedArrayError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            NestedArrayError::DifferentShapes { left, right } => {
                format!("Table1 {left} not equal to {right}")
            }
        }
    }
}

/// One node of the tree: either more branches or the leaf row of values.
#[derive(Debug, Clone, PartialEq)]
enum Node {
    /// A branch that has not been created. Only reachable past the preallocated dimensions.
    Absent,
    Branch(Vec<Node>),
    Leaf(Option<SharedDatum>),
}

/// `NestedIntegerArray<RecalDatum>`.
#[derive(Debug, Clone, PartialEq)]
pub struct NestedIntegerArray {
    dimensions: Vec<usize>,
    data: Vec<Node>,
}

impl NestedIntegerArray {
    /// The constructor, which preallocates the first two dimensions and no more.
    pub fn new(dimensions: &[usize]) -> Result<NestedIntegerArray, NestedArrayError> {
        if dimensions.is_empty() {
            return Err(NestedArrayError::NoDimensions);
        }
        let to_preallocate = dimensions.len().min(NUM_DIMENSIONS_TO_PREALLOCATE);
        let data = build(dimensions, 0, to_preallocate);
        Ok(NestedIntegerArray {
            dimensions: dimensions.to_vec(),
            data,
        })
    }

    pub fn dimensions(&self) -> &[usize] {
        &self.dimensions
    }

    /// `get(int...)`: the varargs form, whose **last key is not bounds-checked**.
    ///
    /// Every nested key past its dimension answers `None`; the last one, and any negative key,
    /// is [`NestedArrayError::IndexOutOfBounds`].
    pub fn get(&self, keys: &[i32]) -> Result<Option<SharedDatum>, NestedArrayError> {
        let nested = self.dimensions.len() - 1;
        let mut node = &self.data;
        for (i, key) in keys.iter().take(nested).enumerate() {
            if *key >= self.dimensions[i] as i32 {
                return Ok(None);
            }
            match index(node, *key)? {
                Node::Absent => return Ok(None),
                Node::Branch(children) => node = children,
                // Unreachable while the key count matches the dimensions.
                Node::Leaf(_) => return Ok(None),
            }
        }
        match index(node, keys[nested])? {
            Node::Leaf(value) => Ok(value.clone()),
            _ => Ok(None),
        }
    }

    /// `get1Key(key0)`, which checks **nothing**: not the dimension, not the sign.
    pub fn get1_key(&self, key0: i32) -> Result<Option<SharedDatum>, NestedArrayError> {
        match index(&self.data, key0)? {
            Node::Leaf(value) => Ok(value.clone()),
            _ => Ok(None),
        }
    }

    /// `get2Keys`, which checks **both** keys, unlike [`NestedIntegerArray::get`].
    pub fn get2_keys(&self, key0: i32, key1: i32) -> Result<Option<SharedDatum>, NestedArrayError> {
        self.get_checked(&[key0, key1])
    }

    /// `get3Keys`, which checks all three.
    pub fn get3_keys(
        &self,
        key0: i32,
        key1: i32,
        key2: i32,
    ) -> Result<Option<SharedDatum>, NestedArrayError> {
        self.get_checked(&[key0, key1, key2])
    }

    /// `get4Keys`, which checks all four.
    pub fn get4_keys(
        &self,
        key0: i32,
        key1: i32,
        key2: i32,
        key3: i32,
    ) -> Result<Option<SharedDatum>, NestedArrayError> {
        self.get_checked(&[key0, key1, key2, key3])
    }

    /// The specialised getters' shared body: every key tested against its dimension, then the walk.
    ///
    /// A negative key still reaches the index, because the test is `>=` and not a range check.
    fn get_checked(&self, keys: &[i32]) -> Result<Option<SharedDatum>, NestedArrayError> {
        for (i, key) in keys.iter().enumerate() {
            if *key >= self.dimensions[i] as i32 {
                return Ok(None);
            }
        }
        let mut node = &self.data;
        for key in &keys[..keys.len() - 1] {
            match index(node, *key)? {
                Node::Absent => return Ok(None),
                Node::Branch(children) => node = children,
                Node::Leaf(_) => return Ok(None),
            }
        }
        match index(node, keys[keys.len() - 1])? {
            Node::Leaf(value) => Ok(value.clone()),
            _ => Ok(None),
        }
    }

    /// `put(value, keys...)`.
    ///
    /// The key count is checked, every nested key is checked against its dimension with a message
    /// naming both, and **the last key is not checked**, exactly as in `get`.
    pub fn put(&mut self, value: SharedDatum, keys: &[i32]) -> Result<(), NestedArrayError> {
        if keys.len() != self.dimensions.len() {
            return Err(NestedArrayError::WrongKeyCount {
                expected: self.dimensions.len(),
                provided: keys.len(),
            });
        }
        let nested = self.dimensions.len() - 1;
        for (i, key) in keys.iter().take(nested).enumerate() {
            if *key >= self.dimensions[i] as i32 {
                return Err(NestedArrayError::KeyTooLarge {
                    key: *key,
                    dimension: i,
                    maximum: self.dimensions[i] as i32 - 1,
                });
            }
        }

        let dimensions = self.dimensions.clone();
        let mut node = &mut self.data;
        for (i, key) in keys.iter().take(nested).enumerate() {
            let length = node.len();
            let slot = usize::try_from(*key)
                .ok()
                .filter(|slot| *slot < length)
                .ok_or(NestedArrayError::IndexOutOfBounds {
                    index: *key,
                    length,
                })?;
            // Past the preallocated dimensions the branch is made on demand.
            if let Node::Absent = node[slot] {
                node[slot] = Node::Branch(build(
                    &dimensions,
                    i + 1,
                    dimensions.len().min(NUM_DIMENSIONS_TO_PREALLOCATE),
                ));
            }
            node = match &mut node[slot] {
                Node::Branch(children) => children,
                _ => return Ok(()),
            };
        }
        let length = node.len();
        let slot = usize::try_from(keys[nested])
            .ok()
            .filter(|slot| *slot < length)
            .ok_or(NestedArrayError::IndexOutOfBounds {
                index: keys[nested],
                length,
            })?;
        node[slot] = Node::Leaf(Some(value));
        Ok(())
    }

    /// `getAllValues()`: every datum, in the tree walk's order rather than in insertion order.
    pub fn all_values(&self) -> Vec<SharedDatum> {
        self.all_leaves()
            .into_iter()
            .map(|(_, value)| value)
            .collect()
    }

    /// `getAllLeaves()`: every datum with the keys that reach it.
    pub fn all_leaves(&self) -> Vec<(Vec<i32>, SharedDatum)> {
        let mut out = Vec::new();
        fill_leaves(&self.data, &mut Vec::new(), &mut out);
        out
    }
}

/// The value at one index of a node's children, or the reference's index error.
fn index(node: &[Node], key: i32) -> Result<&Node, NestedArrayError> {
    usize::try_from(key)
        .ok()
        .and_then(|slot| node.get(slot))
        .ok_or(NestedArrayError::IndexOutOfBounds {
            index: key,
            length: node.len(),
        })
}

/// `preallocateArray`: branches down to `to_preallocate`, leaves below that.
fn build(dimensions: &[usize], dimension: usize, to_preallocate: usize) -> Vec<Node> {
    let width = dimensions[dimension];
    if dimension == dimensions.len() - 1 {
        return vec![Node::Leaf(None); width];
    }
    if dimension + 1 >= to_preallocate {
        // Not preallocated: the branch is made by `put` when something is written under it.
        return vec![Node::Absent; width];
    }
    (0..width)
        .map(|_| Node::Branch(build(dimensions, dimension + 1, to_preallocate)))
        .collect()
}

fn fill_leaves(node: &[Node], path: &mut Vec<i32>, out: &mut Vec<(Vec<i32>, SharedDatum)>) {
    for (key, child) in node.iter().enumerate() {
        match child {
            Node::Absent | Node::Leaf(None) => continue,
            Node::Leaf(Some(value)) => {
                path.push(key as i32);
                out.push((path.clone(), Rc::clone(value)));
                path.pop();
            }
            Node::Branch(children) => {
                path.push(key as i32);
                fill_leaves(children, path, out);
                path.pop();
            }
        }
    }
}

/// `RecalUtils.combineTables(table1, table2)`.
///
/// Where the left table already holds a datum the right one's is **combined into it**; where it does
/// not, the right one's **object** is stored, so the two tables share it afterwards. See the module
/// note.
pub fn combine_tables(
    table1: &mut NestedIntegerArray,
    table2: &NestedIntegerArray,
) -> Result<(), CombineError> {
    if table1.dimensions() != table2.dimensions() {
        let join = |dimensions: &[usize]| {
            dimensions
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        return Err(CombineError::Nested(NestedArrayError::DifferentShapes {
            left: join(table1.dimensions()),
            right: join(table2.dimensions()),
        }));
    }
    for (keys, value) in table2.all_leaves() {
        match table1.get(&keys).map_err(CombineError::Nested)? {
            Some(mine) => {
                let other = value.borrow().clone();
                mine.borrow_mut()
                    .combine(&other)
                    .map_err(CombineError::Datum)?;
            }
            None => table1.put(value, &keys).map_err(CombineError::Nested)?,
        }
    }
    Ok(())
}

/// Either the array's refusal or the datum's.
#[derive(Debug, Clone, PartialEq)]
pub enum CombineError {
    Nested(NestedArrayError),
    Datum(RecalDatumError),
}

impl CombineError {
    pub fn message(&self) -> String {
        match self {
            CombineError::Nested(error) => error.message(),
            CombineError::Datum(error) => error.message(),
        }
    }
}

/// `RecalibrationTables`: the four tables, in the order the report writes them.
#[derive(Debug, Clone, PartialEq)]
pub struct RecalibrationTables {
    /// Every table, special first, which is the order `numTables` and `getTable` use.
    pub all_tables: Vec<NestedIntegerArray>,
    /// Which covariate each table belongs to, parallel to `all_tables`.
    pub kinds: Vec<CovariateKind>,
    pub num_read_groups: usize,
    pub qual_dimension: usize,
}

/// `EventType.values().length`.
pub const EVENT_DIMENSION: usize = 3;

impl RecalibrationTables {
    /// The one-argument constructor: the read group count comes from the covariate.
    pub fn new(
        covariates: &StandardCovariateList,
    ) -> Result<RecalibrationTables, NestedArrayError> {
        let num_read_groups = (covariates.read_group.maximum_key_value() + 1).max(0) as usize;
        RecalibrationTables::with_read_groups(covariates, num_read_groups)
    }

    /// The two-argument constructor, which is what `safeCombine` uses to make a new set.
    pub fn with_read_groups(
        covariates: &StandardCovariateList,
        num_read_groups: usize,
    ) -> Result<RecalibrationTables, NestedArrayError> {
        let qual_dimension = (covariates.quality_score.maximum_key_value() + 1) as usize;

        let mut all_tables = vec![
            NestedIntegerArray::new(&[num_read_groups, EVENT_DIMENSION])?,
            NestedIntegerArray::new(&[num_read_groups, qual_dimension, EVENT_DIMENSION])?,
        ];
        let mut kinds = vec![CovariateKind::ReadGroup, CovariateKind::QualityScore];
        for kind in covariates.additional_covariates() {
            let width = (covariates.maximum_key_value(kind) + 1) as usize;
            all_tables.push(NestedIntegerArray::new(&[
                num_read_groups,
                qual_dimension,
                width,
                EVENT_DIMENSION,
            ])?);
            kinds.push(kind);
        }

        Ok(RecalibrationTables {
            all_tables,
            kinds,
            num_read_groups,
            qual_dimension,
        })
    }

    pub fn num_tables(&self) -> usize {
        self.all_tables.len()
    }

    pub fn read_group_table(&self) -> &NestedIntegerArray {
        &self.all_tables[0]
    }

    pub fn read_group_table_mut(&mut self) -> &mut NestedIntegerArray {
        &mut self.all_tables[0]
    }

    pub fn quality_score_table(&self) -> &NestedIntegerArray {
        &self.all_tables[1]
    }

    pub fn quality_score_table_mut(&mut self) -> &mut NestedIntegerArray {
        &mut self.all_tables[1]
    }

    /// The tables of the additional covariates, in list order.
    pub fn additional_tables(&self) -> &[NestedIntegerArray] {
        &self.all_tables[2..]
    }

    /// `getTableForCovariate`.
    pub fn table_for_covariate(&self, kind: CovariateKind) -> Option<&NestedIntegerArray> {
        self.kinds
            .iter()
            .position(|seen| *seen == kind)
            .map(|index| &self.all_tables[index])
    }

    /// `makeQualityScoreTable`: the same shape, and a **different** table.
    pub fn make_quality_score_table(&self) -> Result<NestedIntegerArray, NestedArrayError> {
        NestedIntegerArray::new(&[self.num_read_groups, self.qual_dimension, EVENT_DIMENSION])
    }

    /// `isEmpty()`: no table holds a datum.
    pub fn is_empty(&self) -> bool {
        self.all_tables
            .iter()
            .all(|table| table.all_values().is_empty())
    }

    /// `combine(toMerge)`: table by table, in place.
    ///
    /// The reference checks the table **count** and leaves the shapes to `combineTables`.
    pub fn combine(&mut self, to_merge: &RecalibrationTables) -> Result<(), CombineError> {
        if self.num_tables() != to_merge.num_tables() {
            return Err(CombineError::Nested(NestedArrayError::DifferentShapes {
                left: self.num_tables().to_string(),
                right: to_merge.num_tables().to_string(),
            }));
        }
        for (mine, theirs) in self.all_tables.iter_mut().zip(&to_merge.all_tables) {
            combine_tables(mine, theirs)?;
        }
        Ok(())
    }

    /// `safeCombine(left, right)`, which is **not** safe: see the module note. It mutates the datums
    /// of both arguments, because the new table holds the same objects.
    pub fn safe_combine(
        covariates: &StandardCovariateList,
        left: &RecalibrationTables,
        right: &RecalibrationTables,
    ) -> Result<RecalibrationTables, CombineError> {
        let mut combined = RecalibrationTables::with_read_groups(covariates, left.num_read_groups)
            .map_err(CombineError::Nested)?;
        combined.combine(left)?;
        combined.combine(right)?;
        Ok(combined)
    }
}

/// The event type a table's last dimension is indexed by, for readers of this module.
pub fn event_index(event: EventType) -> usize {
    event.ordinal()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn datum(observations: i64) -> SharedDatum {
        Rc::new(RefCell::new(
            RecalDatum::new(observations, 1.0, 30).unwrap(),
        ))
    }

    #[test]
    fn the_varargs_getter_does_not_check_its_last_key_and_the_specialised_one_does() {
        let mut table = NestedIntegerArray::new(&[3, 4, 5]).unwrap();
        table.put(datum(1000), &[1, 2, 3]).unwrap();

        assert!(table.get(&[1, 2, 3]).unwrap().is_some());
        assert!(table.get3_keys(1, 2, 3).unwrap().is_some());
        // The divergence, in one pair of lines.
        assert_eq!(
            table.get(&[1, 2, 5]).unwrap_err(),
            NestedArrayError::IndexOutOfBounds {
                index: 5,
                length: 5
            }
        );
        assert_eq!(table.get3_keys(1, 2, 5).unwrap(), None);
        // A nested key is null both ways.
        assert_eq!(table.get(&[3, 0, 0]).unwrap(), None);
        assert_eq!(table.get3_keys(3, 0, 0).unwrap(), None);
    }

    #[test]
    fn a_negative_key_is_checked_by_nothing() {
        let table = NestedIntegerArray::new(&[3, 4, 5]).unwrap();
        assert!(table.get(&[-1, 0, 0]).is_err());
        assert!(table.get3_keys(-1, 0, 0).is_err());
        let flat = NestedIntegerArray::new(&[2]).unwrap();
        assert!(flat.get1_key(-1).is_err());
        // And get1Key does not check the dimension either.
        assert!(flat.get1_key(2).is_err());
    }

    #[test]
    fn put_names_the_dimension_and_its_maximum() {
        let mut table = NestedIntegerArray::new(&[3, 4, 5]).unwrap();
        assert_eq!(
            table.put(datum(1), &[3, 0, 0]).unwrap_err().message(),
            "Key 3 is too large for dimension 0 (max is 2)"
        );
        assert_eq!(
            table.put(datum(1), &[1, 2]).unwrap_err().message(),
            "Exactly 3 keys should be passed to this NestedIntegerArray but 2 were provided"
        );
        // And the last key is an index error rather than a message, like `get`'s.
        assert_eq!(
            table.put(datum(1), &[0, 0, 5]).unwrap_err(),
            NestedArrayError::IndexOutOfBounds {
                index: 5,
                length: 5
            }
        );
    }

    #[test]
    fn combining_shares_the_datum_it_did_not_have() {
        let mut left = NestedIntegerArray::new(&[2, 3]).unwrap();
        let mut right = NestedIntegerArray::new(&[2, 3]).unwrap();
        let only_right = datum(500);
        right.put(Rc::clone(&only_right), &[1, 1]).unwrap();
        combine_tables(&mut left, &right).unwrap();
        // The same object, not a copy: a change to one is a change to both.
        assert!(Rc::ptr_eq(
            &left.get(&[1, 1]).unwrap().unwrap(),
            &only_right
        ));
    }

    #[test]
    fn the_leaves_come_back_in_tree_order() {
        let mut table = NestedIntegerArray::new(&[2, 3, 4, 5]).unwrap();
        table.put(datum(10), &[1, 2, 3, 4]).unwrap();
        table.put(datum(20), &[0, 0, 0, 0]).unwrap();
        let keys: Vec<Vec<i32>> = table
            .all_leaves()
            .into_iter()
            .map(|(keys, _)| keys)
            .collect();
        // Inserted last, walked first.
        assert_eq!(keys, vec![vec![0, 0, 0, 0], vec![1, 2, 3, 4]]);
    }

    #[test]
    fn an_empty_table_has_no_leaves_however_much_of_it_is_preallocated() {
        let table = NestedIntegerArray::new(&[2, 3, 4, 5]).unwrap();
        assert!(table.all_values().is_empty());
        assert!(table.all_leaves().is_empty());
        assert_eq!(
            NestedIntegerArray::new(&[]).unwrap_err().message(),
            "There must be at least one dimension to an NestedIntegerArray"
        );
    }
}
