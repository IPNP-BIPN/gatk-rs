//! Against a measurement directory of `tools/python-tools/measure_vets.sh`: the port's fit, scores
//! and scorer beside the files the reference left there.
//!
//! The directory is a run of `broadinstitute/gatk:4.6.2.0` and is not committed; goldens come from
//! CI. Without `GATK_RS_VETS_RUN` this test has nothing to compare and says so by passing.
use gatk_tools::h5py_scores;
use gatk_tools::hdf5_reader::Hdf5Reader;
use gatk_tools::isolation_forest::{fit, score_samples, Hyperparameters};
use gatk_tools::scorer_pickle;

/// The reference's bytes with the random digits of its temporary script's name zeroed.
fn masked(mut bytes: Vec<u8>) -> Vec<u8> {
    let needle = b"/tmp/isolation-forest.";
    let at = bytes
        .windows(needle.len())
        .position(|window| window == needle)
        .expect("a script name")
        + needle.len();
    let mut i = at;
    while bytes[i].is_ascii_digit() {
        bytes[i] = b'0';
        i += 1;
    }
    bytes
}

fn matrix(reader: &Hdf5Reader) -> (Vec<String>, Vec<Vec<f64>>) {
    let names = reader.strings("/annotations/names").expect("names");
    let (_, values) = reader.doubles("/annotations/chunk_0").expect("a chunk");
    let rows = values.chunks(names.len()).map(<[f64]>::to_vec).collect();
    (names, rows)
}

#[test]
fn the_port_matches_a_measured_run() {
    let Ok(run) = std::env::var("GATK_RS_VETS_RUN") else {
        return;
    };
    let reader =
        Hdf5Reader::open(std::fs::read(format!("{run}/extract.annot.hdf5")).unwrap()).unwrap();
    let (names, rows) = matrix(&reader);
    let label = |name: &str| reader.doubles(&format!("/labels/{name}")).unwrap().1;
    let (training, calibration, snp) = (label("training"), label("calibration"), label("snp"));
    let hyper = Hyperparameters {
        random_state: 0,
        n_estimators: 100,
        n_jobs: None,
    };
    for (tag, wanted) in [("snp", 1.0), ("indel", 0.0)] {
        let pick = |labels: &[f64]| -> Vec<Vec<f64>> {
            rows.iter()
                .enumerate()
                .filter(|(i, _)| labels[*i] == 1.0 && snp[*i] == wanted)
                .map(|(_, row)| row.clone())
                .collect()
        };
        let train = pick(&training);
        let forest = fit(&names, &train, &hyper).unwrap();

        let reference_pickle = std::fs::read(format!("{run}/train.{tag}.scorer.pkl")).unwrap();
        let reference_scorer =
            scorer_pickle::read(&reference_pickle).expect("the reference's scorer");
        assert_eq!(
            reference_scorer.statistics, forest.statistics,
            "{tag} statistics"
        );
        assert_eq!(reference_scorer.max_samples, forest.max_samples as i64);
        for (index, tree) in forest.scorer().trees.iter().enumerate() {
            assert_eq!(*tree, reference_scorer.trees[index], "{tag} tree {index}");
        }
        assert_eq!(
            reference_scorer.decision_path_lengths,
            forest.decision_path_lengths
        );
        assert_eq!(
            reference_scorer.average_path_lengths, forest.average_path_lengths,
            "{tag} paths"
        );

        let scores = score_samples(&forest.scorer(), &names, &train).unwrap();
        let file = h5py_scores::write(&scores);
        let reference = std::fs::read(format!("{run}/train.{tag}.trainingScores.hdf5")).unwrap();
        assert!(file == reference, "{tag}: training scores differ");
        let calibrated = pick(&calibration);
        let file =
            h5py_scores::write(&score_samples(&reference_scorer, &names, &calibrated).unwrap());
        let reference = std::fs::read(format!("{run}/train.{tag}.calibrationScores.hdf5")).unwrap();
        assert!(file == reference, "{tag}: calibration scores differ");

        let ours = scorer_pickle::write(&forest, scorer_pickle::SCRIPT_NAME);
        let theirs = masked(reference_pickle);
        let first = ours.iter().zip(&theirs).position(|(a, b)| a != b);
        assert!(
            first.is_none() && ours.len() == theirs.len(),
            "{tag}: pickle differs at {first:?}, lengths {} and {}",
            ours.len(),
            theirs.len()
        );
    }
}

/// `2**ys` in the reference's numpy against the port's `pow`, on the points a probe of the image
/// wrote: `ys.bin` and `nppow.bin` are little-endian doubles.
#[test]
fn the_port_pow_matches_a_measured_numpy() {
    let Ok(directory) = std::env::var("GATK_RS_LIBM_DIR") else {
        return;
    };
    let doubles = |name: &str| -> Vec<f64> {
        std::fs::read(format!("{directory}/{name}"))
            .unwrap()
            .chunks_exact(8)
            .map(|b| f64::from_le_bytes(b.try_into().unwrap()))
            .collect()
    };
    let (ys, wanted) = (doubles("ys.bin"), doubles("nppow.bin"));
    assert!(ys.len() >= 300_000 && ys.len() == wanted.len());
    let differ = ys
        .iter()
        .zip(&wanted)
        .filter(|(y, want)| {
            gatk_tools::libm_pow_two::pow_two(**y).map(f64::to_bits) != Some(want.to_bits())
        })
        .count();
    assert_eq!(differ, 0, "of {}", ys.len());
}

/// The doubles of the port's `score.scores.hdf5` against the reference's: Java's HDF5 library holds
/// the time and lays the file out its own way, so its values are what can be compared.
#[test]
fn the_port_scores_match_a_measured_run() {
    let (Ok(run), Ok(port)) = (
        std::env::var("GATK_RS_VETS_RUN"),
        std::env::var("GATK_RS_VETS_PORT"),
    ) else {
        return;
    };
    let scores = |directory: &str| -> Vec<u64> {
        let bytes = std::fs::read(format!("{directory}/score.scores.hdf5")).unwrap();
        let reader = Hdf5Reader::open(bytes).unwrap();
        let (_, values) = reader.doubles("/data/scores").unwrap();
        values.iter().map(|x| x.to_bits()).collect()
    };
    assert!(scores(&run) == scores(&port), "the scores differ");
}
