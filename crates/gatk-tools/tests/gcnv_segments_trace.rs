//! The segmenter's intermediates against the frozen reference's, bit for bit.
//!
//! `tools/python-tools/pp_trace.py` runs `segment_gcnv_calls.py` in `broadinstitute/gatk:4.7.0.0`
//! under `gcnv_frozen_env.sh` and writes every call of every compiled pytensor function, inputs and
//! outputs, and every segment quality, as raw float64 bits. `GCNV_PP_TRACE` names the directory
//! holding `s<index>/trace.tsv`, `GCNV_PP_RUN` the run whose shards it read. Without them the test
//! has nothing to compare and passes.

use gatk_tools::gcnv_segments::*;

struct Record {
    key: String,
    shape: Vec<usize>,
    values: Vec<String>,
}

impl Record {
    fn floats(&self) -> Vec<f64> {
        self.values
            .iter()
            .map(|v| f64::from_bits(u64::from_str_radix(v, 16).unwrap()))
            .collect()
    }
    fn bits(&self) -> Vec<String> {
        self.values.clone()
    }
}

fn bits(values: &[f64]) -> Vec<String> {
    values
        .iter()
        .map(|v| format!("{:016x}", v.to_bits()))
        .collect()
}

fn read_trace(path: &str) -> Vec<Record> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| {
            let mut fields = line.splitn(3, '\t');
            let key = fields.next().unwrap().to_string();
            let shape = fields
                .next()
                .unwrap()
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();
            let values = fields
                .next()
                .unwrap_or("")
                .split(' ')
                .filter(|s| !s.is_empty())
                .map(str::to_string)
                .collect();
            Record { key, shape, values }
        })
        .collect()
}

fn rows(record: &Record) -> Vec<Vec<f64>> {
    let width = *record.shape.last().unwrap();
    record.floats().chunks(width).map(<[f64]>::to_vec).collect()
}

#[test]
fn every_compiled_call_and_every_quality_is_the_references() {
    let (Ok(trace_dir), Ok(run)) = (std::env::var("GCNV_PP_TRACE"), std::env::var("GCNV_PP_RUN"))
    else {
        return;
    };
    let mut checked = [0usize; 5];
    for sample in 0..6 {
        let path = format!("{trace_dir}/s{sample}/trace.tsv");
        if !std::path::Path::new(&path).exists() {
            continue;
        }
        let records = read_trace(&path);
        let get = |i: usize, key: &str| {
            assert_eq!(records[i].key, key, "sample {sample} record {i}");
            &records[i]
        };
        let mut qualities = Vec::new();
        let mut i = 0;
        while i < records.len() {
            let key = records[i].key.clone();
            if key == "f0.in0" {
                let pi = get(i, "f0.in0");
                let pi_kc: Vec<Vec<f64>> = rows(pi);
                let stay = get(i + 1, "f0.in1").floats();
                let lq = rows(get(i + 2, "f0.in2"));
                let (prior, trans) = hmm_specs(&pi_kc, &stay, &lq);
                assert_eq!(
                    bits(&prior),
                    get(i + 4, "f0.out0").bits(),
                    "sample {sample} prior"
                );
                assert_eq!(
                    bits(&trans),
                    get(i + 5, "f0.out1").bits(),
                    "sample {sample} trans"
                );
                checked[0] += 1;
                i += 6;
            } else if key == "f1.in0" {
                let hmm = Hmm {
                    states: records[i].shape[0],
                    log_prior: get(i, "f1.in0").floats(),
                    log_trans: get(i + 1, "f1.in1").floats(),
                    log_emission: get(i + 2, "f1.in2").floats(),
                };
                let fb = forward_backward(&hmm);
                assert_eq!(
                    bits(&fb.log_posterior),
                    get(i + 3, "f1.out0").bits(),
                    "posterior"
                );
                assert_eq!(
                    bits(&[fb.log_data_likelihood]),
                    get(i + 4, "f1.out1").bits(),
                    "likelihood"
                );
                assert_eq!(bits(&fb.alpha), get(i + 5, "f1.out2").bits(), "alpha");
                assert_eq!(bits(&fb.beta), get(i + 6, "f1.out3").bits(), "beta");
                checked[1] += 1;
                i += 7;
                // The Viterbi call takes the same three inputs.
                let expected: Vec<String> = get(i + 3, "f2.out0").values.clone();
                let path: Vec<String> = viterbi(&hmm).iter().map(usize::to_string).collect();
                assert_eq!(path, expected, "sample {sample} viterbi");
                checked[2] += 1;
                i += 4;
            } else if key == "f3.in0" {
                let logp = constrained_path_logp(
                    &get(i, "f3.in0").floats(),
                    &get(i + 1, "f3.in1").floats(),
                    &get(i + 2, "f3.in2").floats(),
                    &get(i + 3, "f3.in3").floats(),
                    get(i + 4, "f3.in4").floats()[0],
                );
                assert_eq!(bits(&[logp]), get(i + 5, "f3.out0").bits(), "constrained");
                checked[3] += 1;
                i += 6;
            } else {
                if key.starts_with("q.") {
                    qualities.push(records[i].values[0].clone());
                }
                i += 1;
            }
        }
        let inputs = read_sample_inputs(
            &[format!("{run}/gcnv/g-model")],
            &[format!("{run}/gcnv/g-calls")],
            &format!("{run}/ploidy/p-calls"),
            sample,
        )
        .unwrap();
        let mut ours = Vec::new();
        for segment in segments(&inputs).unwrap() {
            let q = segment.qualities;
            if segment.num_points > 1 {
                ours.extend(bits(&q));
            } else {
                ours.extend(bits(&q[..1]));
            }
        }
        assert_eq!(ours, qualities, "sample {sample} qualities");
        checked[4] += qualities.len();
    }
    eprintln!(
        "hmm specs {}, forward-backward {}, viterbi {}, constrained {}, qualities {}",
        checked[0], checked[1], checked[2], checked[3], checked[4]
    );
}
