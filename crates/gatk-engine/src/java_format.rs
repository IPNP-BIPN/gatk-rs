//! `java.util.Formatter`'s `%f`, for the call sites that print a double into an output file.
//!
//! This is not a general formatter. It is the one conversion two ported tools need and cannot get
//! from Rust's own: **Java rounds HALF_UP on the decimal expansion of the double, and Rust rounds
//! half-to-even on the same expansion**. They disagree whenever the first dropped digit is a 5 that
//! terminates the expansion, which is a value a ratio can easily take: `0.0625` prints `0.063` in
//! Java and `0.062` in Rust.
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
    // The exact decimal expansion of the double, then a half-up rounding of it.
    let text = format!("{:.*}", 30, value.abs());
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
    format!("{sign}{whole}.{fraction}")
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
