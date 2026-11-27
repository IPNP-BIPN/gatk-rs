//! What `cohort_determine_ploidy_and_depth.py` writes, byte for byte: Python's `repr` of a float,
//! pandas' `to_csv` of a float column, gcnvkernel's `write_ndarray_to_tsv` and its JSON, and the
//! numbers in the per-sample tables of `io_ploidy.SamplePloidyWriter`. (Those tables are put
//! together in `gcnv_ploidy.rs`, which formats their floats with [`python_repr`]; the
//! gcnvkernel sources are under `src/main/python/org/broadinstitute/hellbender/gcnvkernel` in
//! GATK 4.7.0.0.) CPython, numpy and pandas have no local source, so what is said about them
//! here is held to the test rows at the end of this file.
//!
//! # Reading this file without reading Rust
//!
//! - `&[f64]` is a borrowed, read-only view of a `double[]`; the caller keeps the list.
//! - `format!(..)` is Java's `String.format`, with `{}` for each argument and `{name}` for a
//!   variable in scope. `{:e}` is exponent form and `{:02}` zero-pads to two digits.
//! - `.expect("..")` unwraps a value that cannot be missing, and stops the program with that
//!   text if it is.

/// Python's `repr(float)` (and numpy's `str` of a float64, which pandas' `to_csv` writes): the
/// shortest digits that read back, in exponent form below 1e-4 or from 1e16 on.
///
/// **How**: take Rust's shortest round-trip digits in exponent form, then lay them out as
/// CPython's `repr` does: exponent form when the decimal exponent is below -4 or at least 16
/// (so `0.0001` stays `0.0001` and `0.00001` is `1e-05`), plain decimals otherwise, with `.0`
/// added to a whole number.
///
/// **Why it is written out**: Rust's own `{}` for an `f64` never switches to exponent form and
/// writes `2` for two, and its `{:e}` writes `1e-5`, not Python's two-digit `1e-05`.
pub fn python_repr(value: f64) -> String {
    // Python spells the special values `nan`, `inf`, `-inf` and keeps the sign of zero.
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
    // `split_once` returns the text before and after the first `e`, as a pair.
    let (mantissa, exponent) = scientific.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("an integer exponent");
    let negative = mantissa.starts_with('-');
    // The significant digits alone: the sign and the decimal point are dropped.
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let sign = if negative { "-" } else { "" };
    // `(-4..16)` is the half-open range -4 to 15, so this is "below -4 or 16 and up": Python's
    // exponent form, `d.ddd` then `e`, a sign, and at least two exponent digits.
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
    // Plain decimals. `point` is how many digits stand before the decimal point: none (pad with
    // zeros after `0.`), all of them or more (pad with zeros and add `.0`), or some.
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
///
/// The header it writes (`io_commons.py`, GATK 4.7.0.0):
///
/// ```text
/// if write_shape_info:
///     header += compose_sam_comment(io_consts.shape_key_value, repr(shape)) + '\n'
///     header += compose_sam_comment(io_consts.type_key_value, str(dtype)) + '\n'
/// if extra_comment_lines is not None:
///     header += '\n'.join(comment + comment_line for comment_line in extra_comment_lines) + '\n'
/// ```
///
/// `compose_sam_comment` is `@`, `CO`, a tab, the key, `:` and the value, so the first lines are
/// `@CO\tshape:(5,)` and `@CO\tdtype:float64`. The column name of a vector is
/// `output_column_prefix + "0"`, `VALUE_0`, and the values are pandas' `to_csv` with no index and
/// no header.
///
/// An empty `extra_comments` stands for Python's `None` (the global model files pass none). An
/// empty Python LIST would still write the `'\n'` above, an empty line, which no caller asks for.
pub fn ndarray_tsv(values: &[f64], extra_comments: &[String]) -> String {
    // `repr` of a one-element tuple keeps its trailing comma, `(1,)`, so both branches print the
    // same form; a per-sample file is `np.atleast_1d` of one value and always `(1,)`.
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
///
/// `io_commons.write_dict_to_json_file`:
///
/// ```text
/// filtered_dict = {k: v for k, v in sorted(dict_to_write.items()) if k not in ignored_keys}
/// with open(output_file, 'w') as fp:
///     json.dump(filtered_dict, fp, indent=1)
///     fp.write('\n')  # json.dump does not add newline to end of file
/// ```
///
/// So the keys come out in sorted order, one per line indented by one space, `"key": value`, the
/// lines joined by `,` and a newline. The test below holds that shape.
pub fn json_floats(entries: &[(&str, f64)]) -> String {
    // `to_vec` copies the borrowed list so it can be sorted; `a.0` is a pair's first half, the
    // key. Rust's `str` order is by code point, as Python's `sorted` of string keys is.
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
