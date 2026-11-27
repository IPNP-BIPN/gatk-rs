//! `java.util.Formatter`'s `%f`, for the call sites that print a double into an output file.
//!
//! This is not a general formatter. It is the one conversion the ported tools that print a double
//! with `%.Nf` or `%.NE` need and cannot get from Rust's own: **Java rounds HALF_UP on the digits
//! `Double.toString` produces, and Rust's `{:.N}` rounds half-to-even on the double's exact
//! expansion**. They disagree on a tie, which is a value a ratio can easily take: `0.0625` prints
//! `0.063` in Java and `0.062` in Rust (`string_format.txt.gz`, `format one-sixteenth 3=0.063`).
//! They also disagree where the shortest digits end in a 5 the exact expansion does not reach,
//! which the next section is about.
//!
//! # It rounds the shortest digits, not the exact expansion
//!
//! `Formatter` does not see the double's exact decimal expansion. It takes the digits
//! `Double.toString` would produce -- the shortest that round-trip -- and pads or rounds THOSE to
//! the requested scale. So a value whose shortest form ends in a five rounds **up** even where its
//! exact expansion would round down: `2.675` is `2.67499999999999982...` exactly, and Java prints
//! `2.68` at two places where C's printf prints `2.67`. And a large value is padded with zeros
//! rather than expanded: `1e300` is a one and three hundred zeros, not the
//! `1000000000000000052504760255204420248704...` the double actually is.
//!
//! This port worked from the exact expansion until the `string-format` golden was measured, which
//! is what #262 was filed about. It goes through [`crate::tsv_table::java_double_to_string`] now,
//! which is the same `Double.toString` the reference's formatter starts from.
//!
//! The three non-finite spellings are Java's too: `NaN`, `Infinity`, `-Infinity`, with no sign on
//! `NaN`. `ClipReads` reaches all three through `(100.0 * n) / 0` when its `--read` argument names
//! no read at all.
//!
//! It lives here rather than in any one caller because two implementations of one rounding rule is
//! how a rounding rule drifts.
//!
//! # Reading this file without reading Rust
//!
//! - The digits are handled as bytes. `b'0'` is the byte value of the character `0` (Java's
//!   `(byte) '0'`), so `b - b'0'` turns a digit character into its number and `(d + b'0') as char`
//!   turns it back; `Vec<u8>` is a growable `byte[]`, Java's `ArrayList<Byte>` without the boxing.
//! - `text.split_once('.')` is "the part before the first `.` and the part after", or nothing if
//!   there is no `.`; `.expect("...")` takes the value and panics with that message if there is
//!   none, which here marks a case `Double.toString` cannot produce.
//! - `a.bytes().chain(b).map(..).collect()` is a Java stream: the bytes of `a` followed by those of
//!   `b`, each mapped, gathered into a new list. `.take(n)` keeps only the first `n`.
//! - `s[..n]` and `s[n..]` are the first `n` elements and the rest (Java's `substring(0, n)` and
//!   `substring(n)`); `x as i32` is a cast that, between integer types, truncates or wraps like
//!   Java's `(int)` rather than failing.
//! - `format!("{sign}{whole}")` is `String.format` with the variables named in the braces.

/// `String.format("%.Nf", value)` (`java.util.Formatter`'s `%f`), half-up on the digits
/// `Double.toString` would produce, not on the exact decimal expansion.
///
/// **How**: take `Double.toString` of the magnitude, lay it out as a plain decimal with one
/// fraction digit more than `places` (padding with zeros), keep `places` digits and round up when
/// the next one is 5 or more, carrying leftwards. The sign is put back last.
///
/// **Why not the exact expansion**: `2.675` is `2.67499999999999982...` exactly, so a rounding of
/// the expansion gives `2.67`; its `Double.toString` is `2.675`, and Java prints `2.68`
/// (`string_format.txt.gz`, `tostring two-point-six-seven-five 2.675` and
/// `format two-point-six-seven-five 2=2.68`). The tests below use `0.0625`, an exact tie, which
/// shows the half-up but not the choice of digits; the golden is what shows that.
///
/// The sign comes from `is_sign_negative`, which is true for `-0.0`, so negative zero prints
/// `-0.00` (`format negative-zero 2=-0.00`). `NaN` prints with no sign.
pub fn format_decimals(value: f64, places: usize) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    // The digits `Double.toString` produces, laid out as a plain decimal with one more fraction
    // digit than the scale asks for, then a half-up rounding of THOSE. The magnitude is rounded,
    // so a negative tie rounds away from zero (`format negative 2=-2.68`), which is HALF_UP.
    let text = plain_decimal(value.abs(), places + 1);
    let (whole, fraction) = text.split_once('.').expect("a decimal point");
    let mut digits: Vec<u8> = whole
        .bytes()
        .chain(fraction.bytes().take(places))
        .map(|b| b - b'0')
        .collect();
    let round_up = fraction.as_bytes()[places] >= b'5';
    // Only the first dropped digit is looked at: the digits after it cannot change a HALF_UP
    // decision on a decimal string. The carry walks left over nines, and past the first digit a
    // new leading 1 is inserted, so `9.999` at two places is `10.00`.
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

/// `String.format("%.NE", value)`: scientific notation over the SAME shortest digits `%f` starts
/// from, rounded half-up to `places` digits after the point and padded with zeros past them.
///
/// It is not the exact expansion rounded: `FormattedFloatingDecimal` never sees more digits than
/// `Double.toString` has, so a value whose shortest form holds sixteen significant digits prints a
/// seventeenth that is always `0` at `%.16E`. `VariantRecalibrator`'s model report is written this
/// way, `7.0276012385529370E-01` where the exact value continues `...53703...`.
///
/// The exponent is signed and at least two digits wide, and a non-finite value is Java's own
/// upper-cased spelling, since `%E` upper-cases the whole conversion.
pub fn format_scientific_upper(value: f64, places: usize) -> String {
    if value.is_nan() {
        return "NAN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "INFINITY" } else { "-INFINITY" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    if value == 0.0 {
        let zeros = "0".repeat(places);
        return if places == 0 {
            format!("{sign}0E+00")
        } else {
            format!("{sign}0.{zeros}E+00")
        };
    }
    let shortest = crate::tsv_table::java_double_to_string(value.abs());
    // `split_once` answers an `Option`: `Some((before, after))` when `Double.toString` used its
    // `E` form (`1.0E-7`), `None` for the plain form (`0.5`), and `match` must handle
    // both, which the compiler checks.
    let (mantissa, exponent) = match shortest.split_once('E') {
        Some((mantissa, exponent)) => (mantissa, exponent.parse::<i32>().expect("an exponent")),
        None => (shortest.as_str(), 0),
    };
    let point = mantissa
        .find('.')
        .expect("Double.toString always has a point") as i32;
    let all: Vec<u8> = mantissa
        .bytes()
        .filter(|b| *b != b'.')
        .map(|b| b - b'0')
        .collect();
    // The value is 0.<all> times ten to the (point + exponent); leading zeros move the exponent.
    let leading = all.iter().take_while(|digit| **digit == 0).count();
    let mut digits: Vec<u8> = all[leading..].to_vec();
    while digits.last() == Some(&0) && digits.len() > 1 {
        digits.pop();
    }
    let mut scientific = point + exponent - leading as i32 - 1;
    let wanted = places + 1;
    if digits.len() > wanted {
        let round_up = digits[wanted] >= 5;
        digits.truncate(wanted);
        if round_up {
            let mut index = wanted;
            loop {
                if index == 0 {
                    digits.insert(0, 1);
                    digits.pop();
                    scientific += 1;
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
    }
    while digits.len() < wanted {
        digits.push(0);
    }
    let rendered: String = digits.iter().map(|d| (d + b'0') as char).collect();
    let body = if places == 0 {
        rendered
    } else {
        format!("{}.{}", &rendered[..1], &rendered[1..])
    };
    let exponent_sign = if scientific < 0 { '-' } else { '+' };
    format!("{sign}{body}E{exponent_sign}{:02}", scientific.abs())
}

/// `Double.toString(value)` laid out as a plain decimal, with at least `places` fraction digits.
///
/// `Double.toString` answers `1.0E300` for a large value and `4.9E-324` for a tiny one; the
/// exponent is applied by moving the point and padding with zeros, which is what
/// `FormattedFloatingDecimal` does and is why a large value prints as a one and three hundred
/// zeros.
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

    // Three layouts: the point lands before every digit (`1.0E-7` gives `0` and `00000010`),
    // after all of them (`1.0E21` gives the digits padded with zeros and no fraction), or among
    // them. `if` is an expression here, so each branch's last value is what the pair is set to.
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
