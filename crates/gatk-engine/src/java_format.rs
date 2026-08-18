//! `java.util.Formatter`'s `%f`, for the call sites that print a double into an output file.
//!
//! This is not a general formatter. It is the one conversion two ported tools need and cannot get
//! from Rust's own: **Java rounds HALF_UP on the decimal expansion of the double, and Rust rounds
//! half-to-even on the same expansion**. They disagree whenever the first dropped digit is a 5 that
//! terminates the expansion, which is a value a ratio can easily take: `0.0625` prints `0.063` in
//! Java and `0.062` in Rust.
//!
//! # It rounds the shortest digits, not the exact expansion
//!
//! `Formatter` does not see the double's exact decimal expansion. It takes the digits
//! `Double.toString` would produce -- the shortest that round-trip -- and pads or rounds THOSE to the
//! requested scale. So a value whose shortest form ends in a five rounds **up** even where its exact
//! expansion would round down: `2.675` is `2.67499999999999982...` exactly, and Java prints `2.68`
//! at two places where C's printf prints `2.67`. And a large value is padded with zeros rather than
//! expanded: `1e300` is a one and three hundred zeros, not the
//! `1000000000000000052504760255204420248704...` the double actually is.
//!
//! This port worked from the exact expansion until the `string-format` golden was measured, which is
//! what #262 was filed about. It goes through [`crate::tsv_table::java_double_to_string`] now, which
//! is the same `Double.toString` the reference's formatter starts from.
//!
//! The three non-finite spellings are Java's too: `NaN`, `Infinity`, `-Infinity`, with no sign on
//! `NaN`. `ClipReads` reaches all three through `(100.0 * n) / 0` when its `--read` argument names
//! no read at all.
//!
//! It lives here rather than in either caller because two implementations of one rounding rule is
//! how a rounding rule drifts.

/// `String.format("%.Nf", value)`, half-up on the decimal expansion.
pub fn format_decimals(value: f64, places: usize) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    // The digits `Double.toString` produces, laid out as a plain decimal with one more fraction
    // digit than the scale asks for, then a half-up rounding of THOSE.
    let text = plain_decimal(value.abs(), places + 1);
    let (whole, fraction) = text.split_once('.').expect("a decimal point");
    let mut digits: Vec<u8> = whole
        .bytes()
        .chain(fraction.bytes().take(places))
        .map(|b| b - b'0')
        .collect();
    let round_up = fraction.as_bytes()[places] >= b'5';
    if round_up {
        let mut index = digits.len();
        loop {
            if index == 0 {
                digits.insert(0, 1);
                break;
            }
            index -= 1;
            if digits[index] == 9 {
                digits[index] = 0;
            } else {
                digits[index] += 1;
                break;
            }
        }
    }
    let split = digits.len() - places;
    let whole: String = digits[..split].iter().map(|d| (d + b'0') as char).collect();
    let fraction: String = digits[split..].iter().map(|d| (d + b'0') as char).collect();
    let sign = if value.is_sign_negative() { "-" } else { "" };
    // `%.0f` prints no point at all, which is `Formatter`'s doing and not the digits'.
    if places == 0 {
        return format!("{sign}{whole}");
    }
    format!("{sign}{whole}.{fraction}")
}

/// `Double.toString(value)` laid out as a plain decimal, with at least `places` fraction digits.
///
/// `Double.toString` answers `1.0E300` for a large value and `4.9E-324` for a tiny one; the exponent
/// is applied by moving the point and padding with zeros, which is what `FormattedFloatingDecimal`
/// does and is why a large value prints as a one and three hundred zeros.
fn plain_decimal(value: f64, places: usize) -> String {
    let shortest = crate::tsv_table::java_double_to_string(value);
    let (mantissa, exponent) = match shortest.split_once('E') {
        Some((mantissa, exponent)) => (mantissa, exponent.parse::<i32>().expect("an exponent")),
        None => (shortest.as_str(), 0),
    };
    let point = mantissa
        .find('.')
        .expect("Double.toString always has a point");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    // Where the point sits once the exponent has moved it.
    let position = point as i32 + exponent;

    let (mut whole, mut fraction) = if position <= 0 {
        ("0".to_string(), "0".repeat((-position) as usize) + &digits)
    } else if position as usize >= digits.len() {
        (
            digits.clone() + &"0".repeat(position as usize - digits.len()),
            String::new(),
        )
    } else {
        (
            digits[..position as usize].to_string(),
            digits[position as usize..].to_string(),
        )
    };
    if whole.is_empty() {
        whole = "0".to_string();
    }
    // `Double.toString`'s trailing `.0` is a digit like any other, and the padding is the
    // formatter's.
    while fraction.len() < places {
        fraction.push('0');
    }
    format!("{whole}.{fraction}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_formatter_rounds_half_up_where_rust_would_round_to_even() {
        assert_eq!(format_decimals(0.0625, 3), "0.063");
        assert_eq!(format!("{:.3}", 0.0625_f64), "0.062", "which is not Java's");
    }

    #[test]
    fn the_non_finite_spellings_are_javas() {
        assert_eq!(format_decimals(f64::NAN, 2), "NaN");
        assert_eq!(format_decimals(f64::INFINITY, 2), "Infinity");
        assert_eq!(format_decimals(f64::NEG_INFINITY, 2), "-Infinity");
    }

    #[test]
    fn the_percentages_clip_reads_prints() {
        assert_eq!(format_decimals(100.0 * 2.0 / 7.0, 2), "28.57");
        assert_eq!(format_decimals(100.0 * 10.0 / 65.0, 2), "15.38");
        assert_eq!(format_decimals(100.0 * 7.0 / 7.0, 2), "100.00");
        let nothing_examined = 0.0_f64;
        assert_eq!(
            format_decimals(100.0 * nothing_examined / nothing_examined, 2),
            "NaN"
        );
    }

    #[test]
    fn a_carry_runs_the_whole_way_up() {
        assert_eq!(format_decimals(9.999, 2), "10.00");
        assert_eq!(format_decimals(0.999, 2), "1.00");
    }
}
