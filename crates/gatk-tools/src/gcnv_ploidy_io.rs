//! What `cohort_determine_ploidy_and_depth.py` writes, byte for byte: Python's `repr` of a float,
//! pandas' `to_csv` of a float column, gcnvkernel's `write_ndarray_to_tsv` and its JSON, and the
//! per-sample tables of `io_ploidy.SamplePloidyWriter`.

/// Python's `repr(float)` (and numpy's `str` of a float64, which pandas' `to_csv` writes): the
/// shortest digits that read back, in exponent form below 1e-4 or from 1e16 on.
pub fn python_repr(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0"
        } else {
            "0.0"
        }
        .to_string();
    }
    // Rust's `{:e}` is the same shortest round-trip digits, as `d.ddde<exp>`.
    let scientific = format!("{value:e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("an integer exponent");
    let negative = mantissa.starts_with('-');
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let sign = if negative { "-" } else { "" };
    if !(-4..16).contains(&exponent) {
        let mut out = format!("{sign}{}", &digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        out.push_str(&format!(
            "e{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        ));
        return out;
    }
    let point = exponent + 1;
    if point <= 0 {
        format!("{sign}0.{}{digits}", "0".repeat((-point) as usize))
    } else if point as usize >= digits.len() {
        format!(
            "{sign}{digits}{}.0",
            "0".repeat(point as usize - digits.len())
        )
    } else {
        let (whole, fraction) = digits.split_at(point as usize);
        format!("{sign}{whole}.{fraction}")
    }
}

/// `io_commons.write_ndarray_to_tsv` of a float64 vector: the shape and dtype comments, the
/// extra comment lines, `VALUE_0`, then one value per line.
pub fn ndarray_tsv(values: &[f64], extra_comments: &[String]) -> String {
    let shape = if values.len() == 1 {
        "(1,)".to_string()
    } else {
        format!("({},)", values.len())
    };
    let mut out = format!("@CO\tshape:{shape}\n@CO\tdtype:float64\n");
    for line in extra_comments {
        out.push('@');
        out.push_str(line);
        out.push('\n');
    }
    out.push_str("VALUE_0\n");
    for value in values {
        out.push_str(&python_repr(*value));
        out.push('\n');
    }
    out
}

/// `json.dump(dict, indent=1)` of the sorted float settings, and the newline gcnvkernel adds.
pub fn json_floats(entries: &[(&str, f64)]) -> String {
    let mut sorted = entries.to_vec();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    let body: Vec<String> = sorted
        .iter()
        .map(|(key, value)| format!(" \"{key}\": {}", python_repr(*value)))
        .collect();
    format!("{{\n{}\n}}\n", body.join(",\n"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floats_are_written_as_python_writes_them() {
        assert_eq!(python_repr(0.15295824923253837), "0.15295824923253837");
        assert_eq!(python_repr(-8.486688839284168), "-8.486688839284168");
        assert_eq!(python_repr(2.0), "2.0");
        assert_eq!(python_repr(0.0001), "0.0001");
        assert_eq!(python_repr(0.00001), "1e-05");
        assert_eq!(python_repr(1e16), "1e+16");
        assert_eq!(python_repr(1.5e300), "1.5e+300");
        assert_eq!(python_repr(123456789012345.6), "123456789012345.6");
        assert_eq!(python_repr(-2.5e-7), "-2.5e-07");
    }

    #[test]
    fn the_ploidy_config_is_the_references_json() {
        assert_eq!(
            json_floats(&[
                ("psi_s_scale", 1e-4),
                ("mean_bias_sd", 1.0),
                ("mapping_error_rate", 0.3),
                ("psi_j_scale", 1e-3)
            ]),
            "{\n \"mapping_error_rate\": 0.3,\n \"mean_bias_sd\": 1.0,\n \"psi_j_scale\": 0.001,\n \
             \"psi_s_scale\": 0.0001\n}\n"
        );
    }
}
