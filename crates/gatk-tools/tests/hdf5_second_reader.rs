//! hdf5-pure, a reader that does not share our reading of the HDF5 specification, opens what
//! `hdf5_writer` writes (IPNP-BIPN/gatk-rs#1294).
//!
//! Until this test the only independent reader of the port's HDF5 was h5py in the covering-array
//! harness, which needs Docker. This one runs in `cargo test`: every file is walked from the root
//! by hdf5-pure, which verifies each object header's `lookup3` checksum, and the tree it finds
//! (groups, datasets, datatype, shape, every value bit for bit) must be the tree that was written,
//! and must agree with what the port's own `hdf5_reader` makes of the same bytes.

use std::collections::BTreeMap;

use gatk_engine::interval::SimpleInterval;
use gatk_tools::hdf5_reader::Hdf5Reader;
use gatk_tools::hdf5_writer::Hdf5File;

/// A dataset as hdf5-pure reads it.
#[derive(Debug, PartialEq)]
enum Leaf {
    /// Shape and the doubles' bits, so a NaN or a -0.0 compares exactly.
    Doubles(Vec<u64>, Vec<u64>),
    Strings(Vec<String>),
}

/// Walk the file from its root with hdf5-pure: every group and dataset by absolute path.
fn walk(bytes: &[u8]) -> (Vec<String>, BTreeMap<String, Leaf>) {
    let file = hdf5_pure::File::from_bytes(bytes.to_vec()).expect("hdf5-pure opens the file");
    let mut groups = Vec::new();
    let mut datasets = BTreeMap::new();
    let mut pending = vec![(String::new(), file.root())];
    while let Some((path, group)) = pending.pop() {
        groups.push(if path.is_empty() {
            "/".to_string()
        } else {
            path.clone()
        });
        for name in group.datasets().expect("hdf5-pure lists the datasets") {
            let full = format!("{path}/{name}");
            let dataset = group.dataset(&name).expect("hdf5-pure opens the dataset");
            let shape = dataset.shape().expect("a shape");
            let leaf = match dataset.dtype().expect("a datatype") {
                hdf5_pure::DType::F64 => Leaf::Doubles(
                    shape,
                    dataset
                        .read_f64()
                        .expect("doubles")
                        .iter()
                        .map(|v| v.to_bits())
                        .collect(),
                ),
                hdf5_pure::DType::String => {
                    Leaf::Strings(dataset.read_string().expect("fixed-length strings"))
                }
                other => panic!("{full}: a datatype the writer never writes: {other}"),
            };
            datasets.insert(full, leaf);
        }
        for name in group.groups().expect("hdf5-pure lists the groups") {
            let child = group.group(&name).expect("hdf5-pure opens the group");
            pending.push((format!("{path}/{name}"), child));
        }
    }
    groups.sort();
    (groups, datasets)
}

/// hdf5-pure's tree against the port's own reader, dataset by dataset.
fn agrees_with_ours(bytes: &[u8], datasets: &BTreeMap<String, Leaf>) {
    let ours = Hdf5Reader::open(bytes.to_vec()).expect("the port's reader opens the file");
    for (path, leaf) in datasets {
        let ours = match leaf {
            Leaf::Doubles(..) => ours.doubles(path).map(|(shape, values)| {
                Leaf::Doubles(shape, values.iter().map(|v| v.to_bits()).collect())
            }),
            Leaf::Strings(_) => ours.strings(path).map(Leaf::Strings),
        };
        assert_eq!(
            ours.as_ref(),
            Some(leaf),
            "{path}: the two readers disagree"
        );
    }
}

#[test]
fn every_shape_the_writer_has() {
    let doubles = [
        0.0,
        -0.0,
        1.5,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::MIN_POSITIVE,
        f64::from_bits(1), // the smallest subnormal
        f64::MAX,
        -1e-300,
    ];
    let mut file = Hdf5File::new();
    file.double_array("/one", &[42.0]);
    file.double_array("/values/doubles", &doubles);
    file.double_array("/values/empty", &[]);
    file.double_matrix(
        "/values/matrix",
        &[vec![1.0, 2.0, 3.0], vec![4.0, f64::NAN, -0.0]],
    );
    file.double_matrix("/values/row", &[doubles.to_vec()]);
    let strings: Vec<String> = ["", "a", "QD", "a longer one, with a comma", "chr1\tchr2\n"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    file.string_array("/text/strings", &strings);
    file.string_array("/text/single", &["only".to_string()]);
    // Deep and wide: a group per level, and more children than one header's first chunk.
    file.double_array("/a/b/c/d/e/leaf", &[1.0]);
    for i in 0..40 {
        file.double_array(&format!("/wide/child_{i:02}"), &[i as f64]);
    }
    let bytes = file.to_bytes();
    let (groups, datasets) = walk(&bytes);

    let mut expected_groups = vec!["/", "/a", "/a/b", "/a/b/c", "/a/b/c/d", "/a/b/c/d/e"];
    expected_groups.extend(["/text", "/values", "/wide"]);
    expected_groups.sort();
    assert_eq!(groups, expected_groups);
    assert_eq!(datasets.len(), 7 + 1 + 40);

    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    assert_eq!(datasets["/one"], Leaf::Doubles(vec![1], bits(&[42.0])));
    assert_eq!(
        datasets["/values/doubles"],
        Leaf::Doubles(vec![doubles.len() as u64], bits(&doubles))
    );
    assert_eq!(
        datasets["/values/empty"],
        Leaf::Doubles(vec![0], Vec::new())
    );
    assert_eq!(
        datasets["/values/matrix"],
        Leaf::Doubles(vec![2, 3], bits(&[1.0, 2.0, 3.0, 4.0, f64::NAN, -0.0]))
    );
    assert_eq!(
        datasets["/values/row"],
        Leaf::Doubles(vec![1, doubles.len() as u64], bits(&doubles))
    );
    assert_eq!(datasets["/text/strings"], Leaf::Strings(strings));
    assert_eq!(
        datasets["/text/single"],
        Leaf::Strings(vec!["only".to_string()])
    );
    assert_eq!(
        datasets["/wide/child_39"],
        Leaf::Doubles(vec![1], bits(&[39.0]))
    );
    agrees_with_ours(&bytes, &datasets);
}

/// What CollectReadCounts writes, through the function the tool calls.
#[test]
fn collect_read_counts_file() {
    let sequences = vec![("chr1".to_string(), 1000), ("chr2".to_string(), 500)];
    let intervals = vec![
        SimpleInterval::new("chr1", 1, 100).unwrap(),
        SimpleInterval::new("chr1", 101, 200).unwrap(),
        SimpleInterval::new("chr2", 1, 500).unwrap(),
    ];
    let file =
        gatk_tools::collect_read_counts::write_hdf5(&sequences, "sample", &intervals, &[3, 0, 17]);
    let bytes = file.to_bytes();
    let (groups, datasets) = walk(&bytes);
    assert!(groups.contains(&"/intervals".to_string()), "{groups:?}");
    assert_eq!(
        datasets["/sample_metadata/sample_name"],
        Leaf::Strings(vec!["sample".to_string()])
    );
    assert_eq!(
        datasets["/counts/values"],
        Leaf::Doubles(
            vec![1, 3],
            [3.0f64, 0.0, 17.0].iter().map(|v| v.to_bits()).collect()
        )
    );
    agrees_with_ours(&bytes, &datasets);
}

/// The control: rename the one link in the root's object header without fixing its checksum. A
/// reader that skipped the checksum would find `/walues` and read it; hdf5-pure must refuse it,
/// so a file that opens above was checked, and not merely tolerated.
#[test]
fn a_header_whose_checksum_is_wrong_is_refused() {
    let mut file = Hdf5File::new();
    file.double_array("/values", &[1.0, 2.0]);
    let mut bytes = file.to_bytes();
    let name = bytes
        .windows(6)
        .position(|w| w == b"values")
        .expect("the link's name");
    bytes[name] = b'w';
    let read = hdf5_pure::File::from_bytes(bytes)
        .and_then(|f| f.dataset("/walues").and_then(|d| d.read_f64()));
    let error = format!("{read:?}");
    assert!(error.contains("hecksum"), "{error}");
}
