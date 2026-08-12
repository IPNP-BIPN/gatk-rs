//! Conformance for `NestedIntegerArray` and `RecalibrationTables` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/RecalibrationTablesDump.java`. Most of this suite is
//! about which lookup answers nothing and which one fails, because a table read the wrong way is a
//! recalibration that silently uses the wrong datum.
//!
//! # What this suite is for
//!
//!  * **the varargs `get` does not bounds-check its last key** and the specialised `getNKeys` do, so
//!    the two disagree on exactly that case;
//!  * **`get1Key` checks nothing at all**;
//!  * **a negative key is checked by nothing anywhere**, because every test is `>=`;
//!  * **`put` bounds-checks only the nested dimensions**, with two exact messages;
//!  * **the shapes come from the covariates**, not from any constant;
//!  * **combining shares objects**, and `safeCombine` mutates both its arguments' datums.

use std::cell::RefCell;
use std::rc::Rc;

use gatk_corpus as corpus;
use gatk_engine::covariates::{RecalibrationArguments, StandardCovariateList};
use gatk_engine::recal_datum::RecalDatum;
use gatk_engine::recalibration_tables::{
    combine_tables, NestedIntegerArray, RecalibrationTables, SharedDatum,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/recalibration_tables.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn field(text: &str, kind: &str, name: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("the golden has no {kind} {name}"))[1]
        .to_string()
}

/// The covariates the dump computed its shapes from, built from the read groups it printed.
fn covariates(text: &str) -> StandardCovariateList {
    let groups: Vec<String> = rows(text, "readgroups")[0][0]
        .split(',')
        .map(|id| id.to_string())
        .collect();
    StandardCovariateList::new(&RecalibrationArguments::default(), &groups).unwrap()
}

fn datum(observations: i64, mismatches: f64, quality: i8) -> SharedDatum {
    Rc::new(RefCell::new(
        RecalDatum::new(observations, mismatches, quality).unwrap(),
    ))
}

/// `RecalDatum.toString()`, which is how the golden writes a datum.
fn text_of(value: &Option<SharedDatum>) -> String {
    match value {
        None => "null".to_string(),
        Some(datum) => datum.borrow_mut().to_text(),
    }
}

/// One lookup's outcome as the golden writes it: the datum, `null`, or `E:<exception>:<message>`.
fn outcome(
    result: Result<Option<SharedDatum>, gatk_engine::recalibration_tables::NestedArrayError>,
) -> String {
    match result {
        Ok(value) => text_of(&value),
        Err(error) => format!("E:ArrayIndexOutOfBoundsException:{}", error.message()),
    }
}

#[test]
fn the_shapes_come_from_the_covariates() {
    let text = golden();
    let covariates = covariates(&text);
    let tables = RecalibrationTables::new(&covariates).unwrap();

    assert_eq!(
        field(&text, "const", "numTables"),
        tables.num_tables().to_string()
    );
    assert_eq!(
        field(&text, "const", "numReadGroups"),
        tables.num_read_groups.to_string()
    );
    assert_eq!(
        field(&text, "const", "qualDimension"),
        tables.qual_dimension.to_string()
    );
    assert_eq!(
        field(&text, "const", "isEmpty"),
        tables.is_empty().to_string()
    );

    let shape = |dimensions: &[usize]| {
        dimensions
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    assert_eq!(
        field(&text, "shape", "readGroup"),
        shape(tables.read_group_table().dimensions())
    );
    assert_eq!(
        field(&text, "shape", "qualityScore"),
        shape(tables.quality_score_table().dimensions())
    );
    for (index, kind) in covariates.additional_covariates().into_iter().enumerate() {
        assert_eq!(
            field(&text, "shape", kind.parsed_name()),
            shape(tables.additional_tables()[index].dimensions()),
            "{}",
            kind.parsed_name()
        );
    }
    // A freshly made quality score table has the same shape and is a different table.
    assert_eq!(
        shape(tables.make_quality_score_table().unwrap().dimensions()),
        shape(tables.quality_score_table().dimensions())
    );
    assert_eq!(
        field(&text, "const", "madeQualityScoreIsQualityScore"),
        "false"
    );
}

/// The five lookups, which do not agree, over every key the golden asks for.
#[test]
fn every_lookup_is_the_reference() {
    let text = golden();
    let expected = |label: &str| field(&text, "get", label);

    let mut table = NestedIntegerArray::new(&[3, 4, 5]).unwrap();
    table.put(datum(1000, 10.0, 30), &[1, 2, 3]).unwrap();

    assert_eq!(outcome(table.get(&[1, 2, 3])), expected("in-range-varargs"));
    assert_eq!(
        outcome(table.get3_keys(1, 2, 3)),
        expected("in-range-3keys")
    );
    assert_eq!(outcome(table.get(&[1, 2, 4])), expected("unset-varargs"));
    assert_eq!(outcome(table.get3_keys(1, 2, 4)), expected("unset-3keys"));
    assert_eq!(
        outcome(table.get(&[0, 0, 0])),
        expected("unset-branch-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(0, 0, 0)),
        expected("unset-branch-3keys")
    );

    // The divergence.
    assert_eq!(
        outcome(table.get(&[1, 2, 5])),
        expected("last-key-too-big-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(1, 2, 5)),
        expected("last-key-too-big-3keys")
    );

    assert_eq!(
        outcome(table.get(&[3, 0, 0])),
        expected("first-key-too-big-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(3, 0, 0)),
        expected("first-key-too-big-3keys")
    );
    assert_eq!(
        outcome(table.get(&[0, 4, 0])),
        expected("second-key-too-big-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(0, 4, 0)),
        expected("second-key-too-big-3keys")
    );

    assert_eq!(
        outcome(table.get(&[-1, 0, 0])),
        expected("negative-first-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(-1, 0, 0)),
        expected("negative-first-3keys")
    );
    assert_eq!(
        outcome(table.get(&[1, 2, -1])),
        expected("negative-last-varargs")
    );
    assert_eq!(
        outcome(table.get3_keys(1, 2, -1)),
        expected("negative-last-3keys")
    );

    let mut flat = NestedIntegerArray::new(&[2]).unwrap();
    flat.put(datum(1000, 10.0, 30), &[0]).unwrap();
    assert_eq!(outcome(flat.get1_key(0)), expected("flat-1key-in-range"));
    assert_eq!(outcome(flat.get1_key(1)), expected("flat-1key-unset"));
    assert_eq!(outcome(flat.get1_key(2)), expected("flat-1key-too-big"));
    assert_eq!(outcome(flat.get1_key(-1)), expected("flat-1key-negative"));
    assert_eq!(outcome(flat.get(&[2])), expected("flat-varargs-too-big"));

    let mut two = NestedIntegerArray::new(&[2, 3]).unwrap();
    two.put(datum(1000, 10.0, 30), &[1, 2]).unwrap();
    assert_eq!(outcome(two.get2_keys(1, 2)), expected("two-2keys"));
    assert_eq!(
        outcome(two.get2_keys(1, 3)),
        expected("two-2keys-last-too-big")
    );
    assert_eq!(
        outcome(two.get(&[1, 3])),
        expected("two-varargs-last-too-big")
    );

    let mut four = NestedIntegerArray::new(&[2, 3, 4, 5]).unwrap();
    four.put(datum(1000, 10.0, 30), &[1, 2, 3, 4]).unwrap();
    assert_eq!(outcome(four.get4_keys(1, 2, 3, 4)), expected("four-4keys"));
    assert_eq!(
        outcome(four.get4_keys(1, 2, 3, 5)),
        expected("four-4keys-last-too-big")
    );
    assert_eq!(
        outcome(four.get(&[1, 2, 3, 5])),
        expected("four-varargs-last-too-big")
    );
    assert_eq!(
        outcome(four.get(&[0, 0, 0, 0])),
        expected("four-unset-branch")
    );
    assert_eq!(
        outcome(four.get4_keys(0, 0, 0, 0)),
        expected("four-4keys-unset-branch")
    );
}

/// `put`'s checks and its two exact messages.
#[test]
fn every_insertion_is_the_reference() {
    let text = golden();
    let expected = |label: &str| field(&text, "put", label);
    let refused = |result: Result<(), gatk_engine::recalibration_tables::NestedArrayError>,
                   exception: &str| match result {
        Ok(()) => "null".to_string(),
        Err(error) => format!("E:{exception}:{}", error.message()),
    };

    let mut table = NestedIntegerArray::new(&[3, 4, 5]).unwrap();
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[1, 2]),
            "IllegalArgumentException"
        ),
        expected("wrong-key-count-too-few")
    );
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[1, 2, 3, 4]),
            "IllegalArgumentException"
        ),
        expected("wrong-key-count-too-many")
    );
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[3, 0, 0]),
            "IllegalArgumentException"
        ),
        expected("first-key-too-big")
    );
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[0, 4, 0]),
            "IllegalArgumentException"
        ),
        expected("second-key-too-big")
    );
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[0, 0, 5]),
            "ArrayIndexOutOfBoundsException"
        ),
        expected("last-key-too-big")
    );
    assert_eq!(
        refused(
            table.put(datum(1, 0.0, 30), &[-1, 0, 0]),
            "ArrayIndexOutOfBoundsException"
        ),
        expected("negative-first-key")
    );

    // Overwriting is silent.
    let mut small = NestedIntegerArray::new(&[2, 2]).unwrap();
    small.put(datum(1000, 10.0, 30), &[0, 0]).unwrap();
    small.put(datum(7, 1.0, 20), &[0, 0]).unwrap();
    assert_eq!(text_of(&small.get(&[0, 0]).unwrap()), expected("overwrite"));

    assert_eq!(
        format!(
            "E:IllegalArgumentException:{}",
            NestedIntegerArray::new(&[]).unwrap_err().message()
        ),
        expected("no-dimensions")
    );
    // A zero-length first dimension makes an array nothing fits in.
    let empty = NestedIntegerArray::new(&[0, 3]).unwrap();
    assert_eq!(
        empty.all_values().len().to_string(),
        expected("zero-length-dimension")
    );
}

/// The tree walk, whose order is the tree's and not the insertion order.
#[test]
fn the_traversal_is_the_reference() {
    let text = golden();
    let value = |label: &str| field(&text, "values", label);

    let mut table = NestedIntegerArray::new(&[2, 3, 4, 5]).unwrap();
    assert_eq!(table.all_values().len().to_string(), value("empty-four"));
    assert_eq!(
        table.all_leaves().len().to_string(),
        value("empty-four-leaves")
    );

    table.put(datum(10, 1.0, 30), &[0, 0, 0, 0]).unwrap();
    table.put(datum(20, 2.0, 30), &[1, 2, 3, 4]).unwrap();
    table.put(datum(30, 3.0, 30), &[0, 1, 0, 1]).unwrap();
    assert_eq!(table.all_values().len().to_string(), value("three-values"));

    let expected_leaves: Vec<Vec<&str>> = rows(&text, "leaf")
        .into_iter()
        .filter(|row| row[0] == "three-values")
        .collect();
    let ours = table.all_leaves();
    assert_eq!(ours.len(), expected_leaves.len());
    for ((keys, datum), row) in ours.iter().zip(&expected_leaves) {
        let keys: Vec<String> = keys.iter().map(|key| key.to_string()).collect();
        assert_eq!(keys.join(","), row[1]);
        assert_eq!(datum.borrow_mut().to_text(), row[2]);
    }

    // The order, which is the tree's.
    let order: Vec<String> = table
        .all_values()
        .iter()
        .map(|datum| datum.borrow_mut().to_text())
        .collect();
    assert_eq!(order.join(";"), value("three-values-order"));

    let mut flat = NestedIntegerArray::new(&[3]).unwrap();
    flat.put(datum(5, 0.0, 30), &[2]).unwrap();
    let leaves = flat.all_leaves();
    let expected: Vec<Vec<&str>> = rows(&text, "leaf")
        .into_iter()
        .filter(|row| row[0] == "flat")
        .collect();
    assert_eq!(leaves.len(), expected.len());
    assert_eq!(leaves[0].0, vec![2]);
    assert_eq!(leaves[0].1.borrow_mut().to_text(), expected[0][2]);
}

/// Combining, and the object sharing that makes `safeCombine` unsafe.
#[test]
fn combining_shares_objects_like_the_reference() {
    let text = golden();

    let mut left = NestedIntegerArray::new(&[2, 3]).unwrap();
    let mut right = NestedIntegerArray::new(&[2, 3]).unwrap();
    left.put(datum(1000, 10.0, 30), &[0, 0]).unwrap();
    right.put(datum(2000, 20.0, 20), &[0, 0]).unwrap();
    let only_right = datum(500, 5.0, 25);
    right.put(Rc::clone(&only_right), &[1, 1]).unwrap();
    left.put(datum(300, 3.0, 35), &[1, 2]).unwrap();

    combine_tables(&mut left, &right).unwrap();

    let expected: Vec<Vec<&str>> = rows(&text, "combine")
        .into_iter()
        .filter(|row| row[0] == "merged")
        .collect();
    let leaves = left.all_leaves();
    assert_eq!(leaves.len(), expected.len());
    for ((keys, datum), row) in leaves.iter().zip(&expected) {
        let keys: Vec<String> = keys.iter().map(|key| key.to_string()).collect();
        assert_eq!(keys.join(","), row[1]);
        assert_eq!(datum.borrow_mut().to_text(), row[2]);
    }

    // The datum the left table did not have is the right table's own object.
    assert_eq!(
        Rc::ptr_eq(&left.get(&[1, 1]).unwrap().unwrap(), &only_right).to_string(),
        field(&text, "shared", "only-right")
    );
    // And the one it did have is not.
    assert_eq!(
        Rc::ptr_eq(
            &left.get(&[0, 0]).unwrap().unwrap(),
            &right.get(&[0, 0]).unwrap().unwrap()
        )
        .to_string(),
        field(&text, "shared", "combined")
    );

    // Mismatched shapes, with both in the message.
    let mut narrow = NestedIntegerArray::new(&[2, 3]).unwrap();
    let wide = NestedIntegerArray::new(&[2, 4]).unwrap();
    assert_eq!(
        combine_tables(&mut narrow, &wide).unwrap_err().message(),
        rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == "combine-different-shapes")
            .unwrap()[2]
    );
}

/// `safeCombine`, which mutates the datums of both arguments because it holds their objects.
#[test]
fn safe_combine_is_not_safe() {
    let text = golden();
    let covariates = covariates(&text);

    let mut one = RecalibrationTables::new(&covariates).unwrap();
    let mut two = RecalibrationTables::new(&covariates).unwrap();
    let ones_datum = datum(100, 1.0, 30);
    one.read_group_table_mut()
        .put(Rc::clone(&ones_datum), &[0, 0])
        .unwrap();
    two.read_group_table_mut()
        .put(datum(200, 2.0, 30), &[0, 0])
        .unwrap();
    two.read_group_table_mut()
        .put(datum(300, 3.0, 30), &[1, 1])
        .unwrap();

    let combined = RecalibrationTables::safe_combine(&covariates, &one, &two).unwrap();

    let expected: Vec<Vec<&str>> = rows(&text, "combine")
        .into_iter()
        .filter(|row| row[0] == "safeCombine")
        .collect();
    let leaves = combined.read_group_table().all_leaves();
    assert_eq!(leaves.len(), expected.len());
    for ((keys, datum), row) in leaves.iter().zip(&expected) {
        let keys: Vec<String> = keys.iter().map(|key| key.to_string()).collect();
        assert_eq!(keys.join(","), row[1]);
        assert_eq!(datum.borrow_mut().to_text(), row[2]);
    }
    assert_eq!(
        combined.is_empty().to_string(),
        field(&text, "const", "safeCombineIsEmpty")
    );

    // The new table holds the left table's object, so the left table's datum changed under it.
    assert_eq!(
        Rc::ptr_eq(
            &combined.read_group_table().get(&[0, 0]).unwrap().unwrap(),
            &one.read_group_table().get(&[0, 0]).unwrap().unwrap()
        )
        .to_string(),
        field(&text, "shared", "safeCombine-left")
    );
    assert!(Rc::ptr_eq(
        &combined.read_group_table().get(&[0, 0]).unwrap().unwrap(),
        &ones_datum
    ));
    assert_eq!(ones_datum.borrow().num_observations(), 300);
}
