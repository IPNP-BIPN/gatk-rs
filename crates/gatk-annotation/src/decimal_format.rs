//! `java.text.DecimalFormat`, for the two patterns `AllelePseudoDepth` emits through.
//!
//! Written from the class's documented behaviour and from measurement of the pinned oracle, not
//! from its source: `java.text` is GPL2-with-classpath-exception, so transcribing it into this
//! crate is the thing htsjdk-rs decision 0014 refuses. What is reproduced here is behaviour that
//! the Javadoc specifies (`HALF_EVEN`, pattern semantics) plus three facts that it does not, each
//! established by running the reference and comparing 5,699,818 formatted values against an
//! independent rule.
//!
//! # It rounds the shortest decimal form, not the value
//!
//! This is the fact everything else follows from, and it is not in the documentation. Formatting
//! `0.1` to forty fraction digits gives `0.1`, not `0.1000000000000000055511151231257827…`, which
//! is what the double actually is. So the digits being rounded are the ones `Double.toString`
//! produces — the shortest decimal that round-trips — and the true binary value is consulted only
//! to break a tie.
//!
//! Two visible consequences, both measured:
//!
//!  * **past 2^53 the tail is zeros.** `-9.81505222706934e219` formats as its sixteen significant
//!    digits followed by two hundred zeros, where the exact value continues `…087125181525995…`;
//!  * **a value can have fewer fraction digits than the pattern allows.** `5.985315667820839e13`
//!    at four fraction digits gives `59853156678208.39`, because the shortest form has run out:
//!    the exact value's `…3906` is not available to be printed.
//!
//! # The tie is broken by the true value, except in one place
//!
//! When the first dropped digit is `5` and nothing follows it in the shortest form, the shortest
//! form *looks* like an exact tie. Usually it is not, and which way the true value sits decides:
//! `0.155` rounds **down** (the double is `0.15499999999999999888…`) while `0.165` rounds **up**
//! (`0.16500000000000000777…`), and no rule that looks only at the digit string gets both right.
//! When the shortest form **is** the value, as for `0.125`, it is a real tie and `HALF_EVEN` sends
//! it to the even neighbour.
//!
//! The exception is when the rounding position falls before the first significant digit, which
//! happens when the whole value is at or below half of the last place the pattern can show. There
//! the two patterns disagree with each other:
//!
//! ```text
//! pattern   value   rounding position   result   the double is
//! #.##      0.005   index 0             0.01     above the tie
//! #.####    5e-5    index 0             0        above the tie
//! ```
//!
//! Same shape, one decade apart, opposite answers. Swept across every decade from `5e-1` to
//! `5e-12`: patterns with **two or fewer** fraction digits follow the value, patterns with three or
//! more round to even — and with no preceding digit, even means down. That boundary is where
//! `DecimalFormat`'s internal fast path stops applying, so the pattern decides which rule runs.
//! It is reproduced here because it is reachable: a pseudo-fraction of exactly `5e-5` is a legal
//! output, and the reference would print `0` for it.
//!
//! # Where this stops being the reference, measured
//!
//! **Below 2^53 there is no divergence at all**, on 903,121 values formatted through both patterns
//! by the pinned oracle and by this module. That covers everything `AllelePseudoDepth` can produce
//! by a wide margin: a pseudo-fraction lies in `[0, 1]` and a pseudo-depth is a sum of posterior
//! counts over the reads at one site.
//!
//! Above 2^53 the reference stops printing the shortest form, and what it prints instead is not one
//! thing. On a sweep of 493 such values it gave the shortest form 472 times, the double's exact
//! value 9 times, and neither 12 times — `2^60` comes out as its exact value rounded to eighteen
//! significant digits. Those are branches inside Java 17's pre-Schubfach `FloatingDecimal`, not a
//! rule, and reproducing them would mean transcribing GPL2 source or fitting an implementation to
//! measurements. Both are refused; see `docs/licence-compatibility-risk.md`.
//!
//! There used to be a second limit here, at sixteen significant digits, and it was this port's
//! fault rather than Java's. See [`closest_of_that_length`].
//!
//! # What is not modelled
//!
//! Locale is pinned to `en-US` by the conformance harness (`-Duser.language=en
//! -Duser.country=US`), so the decimal separator is `.` and the digits are ASCII. Grouping is off
//! for both patterns. `minimumIntegerDigits` is 1, so a pure fraction keeps its leading `0`, and
//! `minimumFractionDigits` is 0, so trailing zeros are dropped.

use std::cmp::Ordering;

/// One `DecimalFormat`, identified by the only part of its pattern that varies here.
///
/// `#.##` and `#.####` differ from each other in exactly one way that the pattern language
/// expresses — how many fraction digits they will show — and in one way it does not, which is the
/// tie rule at the underflow boundary described above.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DecimalFormat {
    max_fraction_digits: usize,
}

/// `AllelePseudoDepth.DEPTH_FORMAT`, the pattern `#.##`.
pub const DEPTH_FORMAT: DecimalFormat = DecimalFormat::new(2);

/// `AllelePseudoDepth.FRACTION_FORMAT`, the pattern `#.####`.
pub const FRACTION_FORMAT: DecimalFormat = DecimalFormat::new(4);

impl DecimalFormat {
    pub const fn new(max_fraction_digits: usize) -> Self {
        Self {
            max_fraction_digits,
        }
    }

    /// `DecimalFormat.format(double)` under `en-US`.
    pub fn format(&self, value: f64) -> String {
        if value.is_nan() {
            // The symbol is `NaN`, and a NaN's sign is not printed.
            return "NaN".to_string();
        }
        let negative = value.is_sign_negative();
        if value.is_infinite() {
            // U+221E, and the sign comes from the negative subpattern.
            return if negative { "-∞" } else { "∞" }.to_string();
        }

        let body = if value == 0.0 {
            "0".to_string()
        } else {
            let (digits, exp10) = shortest_decimal(value.abs());
            self.round(&digits, exp10, value.abs())
        };
        // A negative value that rounds away to nothing still prints its sign: `-0`.
        if negative {
            format!("-{body}")
        } else {
            body
        }
    }

    /// Rounds `0.<digits> * 10^exp10` and renders it.
    fn round(&self, digits: &str, exp10: i32, magnitude: f64) -> String {
        // How many of `digits` survive. A digit at index `i` sits at fraction position
        // `i + 1 - exp10`, so keeping `max_fraction_digits` fraction digits keeps indices below
        // `max_fraction_digits + exp10`.
        let keep = self.max_fraction_digits as i32 + exp10;

        if keep < 0 {
            // Everything is smaller than half the last place: nothing can survive, and there is
            // not even a tie to consider.
            return "0".to_string();
        }

        let mut kept: Vec<u8>;
        let mut exp10 = exp10;

        if keep == 0 {
            // The rounding position falls before the first digit. The result is either zero or one
            // unit in the last place, and this is the boundary where the two patterns disagree.
            let first = digits.as_bytes()[0];
            let round_up = match first.cmp(&b'5') {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => {
                    if digits.len() > 1 {
                        // Something nonzero follows the 5, so the value is past the tie whatever
                        // the true binary value is.
                        true
                    } else if self.max_fraction_digits <= 2 {
                        // The fast path: the true value decides, and an exact tie goes to even,
                        // which with no preceding digit means down.
                        compare_to_exact(digits, exp10, magnitude) == Ordering::Greater
                    } else {
                        // The general path: a lone 5 is treated as a tie, and the digit before it
                        // is absent, which counts as even. So down.
                        false
                    }
                }
            };
            if !round_up {
                return "0".to_string();
            }
            kept = vec![1];
            exp10 += 1;
        } else {
            let keep = keep as usize;
            if keep >= digits.len() {
                // The pattern allows more fraction digits than the shortest form has. Nothing is
                // dropped, and nothing is invented to fill the gap either.
                kept = digits.bytes().map(|b| b - b'0').collect();
            } else {
                kept = digits.bytes().take(keep).map(|b| b - b'0').collect();
                let rest = &digits[keep..];
                let first_dropped = rest.as_bytes()[0];
                let more_follows = rest.bytes().skip(1).any(|b| b != b'0');
                let round_up = match first_dropped.cmp(&b'5') {
                    Ordering::Greater => true,
                    Ordering::Less => false,
                    Ordering::Equal if more_follows => true,
                    // A lone 5 at the cut. Only the true binary value can say whether this is a
                    // tie at all, and only if it is does `HALF_EVEN` look at the digit before it.
                    Ordering::Equal => match compare_to_exact(digits, exp10, magnitude) {
                        // The true value is past the printed 5, so it is not a tie: round up.
                        Ordering::Greater => true,
                        // It is short of it: round down.
                        Ordering::Less => false,
                        // A real tie, and only now does the preceding digit matter.
                        Ordering::Equal => kept.last().is_some_and(|d| d % 2 == 1),
                    },
                };
                if round_up {
                    carry(&mut kept, &mut exp10);
                }
            }
        }

        // `minimumFractionDigits` is 0, so trailing fraction zeros go. The integer part keeps its
        // places, which is what stops 10 becoming 1.
        while kept.len() > 1 && kept.last() == Some(&0) && kept.len() as i32 > exp10 {
            kept.pop();
        }
        render(&kept, exp10)
    }
}

/// Adds one unit to the last kept digit, propagating.
fn carry(kept: &mut Vec<u8>, exp10: &mut i32) {
    for index in (0..kept.len()).rev() {
        if kept[index] == 9 {
            kept[index] = 0;
        } else {
            kept[index] += 1;
            return;
        }
    }
    // Every digit was a 9: the number gained a place.
    kept.insert(0, 1);
    *exp10 += 1;
}

/// Writes `0.<kept> * 10^exp10` without grouping and without scientific notation.
fn render(kept: &[u8], exp10: i32) -> String {
    let mut out = String::new();
    if exp10 <= 0 {
        // `minimumIntegerDigits` is 1, so a pure fraction gets a leading zero.
        out.push_str("0.");
        for _ in 0..(-exp10) {
            out.push('0');
        }
        for digit in kept {
            out.push((b'0' + digit) as char);
        }
    } else {
        let integer_len = exp10 as usize;
        for index in 0..integer_len {
            // Past the end of the shortest form the integer part is padded with zeros. This is
            // where a value beyond 2^53 loses its true tail.
            out.push((b'0' + kept.get(index).copied().unwrap_or(0)) as char);
        }
        if kept.len() > integer_len {
            out.push('.');
            for digit in &kept[integer_len..] {
                out.push((b'0' + digit) as char);
            }
        }
    }
    out
}

/// The shortest decimal that round-trips, as `(digits, exp10)` with the value `0.<digits> * 10^exp10`.
///
/// Rust's `{:e}` produces the shortest round-tripping form, which is what Java 17's
/// `Double.toString` is *intended* to produce. The two are not guaranteed identical: the
/// pre-JDK-19 algorithm occasionally emits one digit more than necessary. An extra digit lands at
/// the seventeenth significant place, so it can only change an answer for a value whose rounding
/// position is out there too — far beyond the four fraction digits used here for anything of
/// ordinary magnitude. The conformance suite is what would catch it.
fn shortest_decimal(value: f64) -> (String, i32) {
    let text = format!("{value:e}");
    let (mantissa, exponent) = text.split_once('e').expect("scientific form");
    let exponent: i32 = exponent.parse().expect("exponent");
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let digits = digits.trim_end_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    // `{:e}` writes one digit before the point, so `0.<digits>` needs one more.
    closest_of_that_length(digits, exponent + 1, value)
}

/// Among the decimals of the length Rust chose, the one nearest the double, ties to even.
///
/// "Shortest, then nearest, then even" is the specification every modern shortest-representation
/// algorithm implements — Ryu, Schubfach, and therefore Java 19 and later. Rust's formatter gets
/// the length right and does not always get the last digit right: `6.985838094673373e14` is
/// exactly `698583809467337.25`, so the two sixteen-digit forms `…337.2` and `…337.3` are
/// **equidistant** and both round-trip. Java prints the even one; Rust prints the other.
///
/// Fixing it here is implementing to the specification, not to Java 17: the reference agrees with
/// the rule, and where it stops agreeing — above 2^53, where its digit generation predates
/// Schubfach — no rule of this kind reaches it anyway.
///
/// The cost is two string parses per call. The exact expansion, which is the expensive part, is
/// computed only when a neighbour also round-trips, which needs the double to be a short decimal
/// and is rare.
fn closest_of_that_length(digits: &str, exp10: i32, value: f64) -> (String, i32) {
    let length = digits.len();
    let Ok(number) = digits.parse::<u128>() else {
        // Beyond what a u128 holds there is no neighbour to consider: the length is already past
        // anything a double can distinguish.
        return (digits.to_string(), exp10);
    };
    // The scale that turns the digit string back into the value it names.
    let scale = exp10 - length as i32;
    let magnitude = value.abs();

    // Reaching the exact expansion needs the double to sit exactly on a midpoint, which needs it
    // to *be* a short decimal. Testing that first is what keeps this affordable: without the gate
    // the expansion ran on nearly every value, because at seventeen digits some twenty neighbouring
    // decimals round-trip. Reading a midpoint back as a double does not help either — it sits
    // inside the same rounding interval and parses to the same double.
    if !may_be_a_short_decimal(magnitude) {
        return (digits.to_string(), exp10);
    }

    for lower in [number.wrapping_sub(1), number] {
        let upper = lower + 1;
        let (lower_text, upper_text) = (lower.to_string(), upper.to_string());
        // Only a neighbour of the *same* length competes. A shorter one would have been chosen
        // already, since Rust's answer is the shortest that round-trips.
        if lower_text.len() != length || upper_text.len() != length {
            continue;
        }
        let both_round_trip = [&lower_text, &upper_text].iter().all(|candidate| {
            format!("{candidate}e{scale}")
                .parse::<f64>()
                .is_ok_and(|parsed| parsed == magnitude)
        });
        if !both_round_trip {
            continue;
        }
        // The midpoint of two consecutive n-digit decimals is that pair's lower half with a 5
        // appended, so the exact comparison this crate already does answers which side the double
        // falls on.
        let chosen = match compare_to_exact(&format!("{lower_text}5"), exp10, magnitude) {
            Ordering::Less => lower_text,
            Ordering::Greater => upper_text,
            // Equidistant, and this is the case Rust gets wrong.
            Ordering::Equal => {
                if lower % 2 == 0 {
                    lower_text
                } else {
                    upper_text
                }
            }
        };
        return (chosen, exp10);
    }
    (digits.to_string(), exp10)
}

/// Where the double's exact value sits relative to its own shortest decimal form.
///
/// `Less` means the shortest form was rounded **up** to reach it, so a apparent tie is really below
/// the halfway point; `Greater` means it was rounded down; `Equal` means the decimal is the value
/// and the tie is real.
///
/// Every double is a finite decimal, so this is an exact comparison rather than an estimate. The
/// expansion is computed by the schoolbook route: `m * 2^e` with `e` negative is
/// `m * 5^-e / 10^-e`, and multiplying a decimal by five is one pass over its digits.
fn compare_to_exact(digits: &str, exp10: i32, magnitude: f64) -> Ordering {
    let (exact_digits, exact_exp10) = exact_decimal(magnitude);
    // Both are written as `0.<digits> * 10^exp`, with no leading zero, so the exponent orders them
    // whenever it differs. It does differ for a value like 1e23, whose shortest form is `1` at
    // exponent 24 while the double itself is 0.99999999999999991611392 at exponent 23.
    if exp10 != exact_exp10 {
        return exact_exp10.cmp(&exp10);
    }
    let shortest = digits.as_bytes();
    for index in 0..shortest.len().max(exact_digits.len()) {
        let ours = exact_digits.get(index).copied().unwrap_or(0);
        let theirs = shortest.get(index).map_or(0, |b| b - b'0');
        match ours.cmp(&theirs) {
            Ordering::Equal => {}
            other => return other,
        }
    }
    Ordering::Equal
}

/// Whether the double could be a decimal short enough for an equidistant pair to exist at all.
///
/// A necessary condition, and a cheap one. Write the value as `odd * 2^power` with `odd` odd. When
/// `power` is negative the exact decimal is `odd * 5^-power` over `10^-power`, and `5^-power` is
/// odd too, so nothing cancels: the expansion has at least `digits(odd) + floor(0.699 * -power)`
/// significant digits. When `power` is positive the value is `odd * 2^power`, at least
/// `digits(odd) + floor(0.301 * power)` digits by the same argument. A tie between two forms of at
/// most eighteen digits cannot happen once either exceeds nineteen.
///
/// False positives are harmless — the exact comparison then runs and finds no tie. False negatives
/// would be a bug, which is why the bound is the pessimistic one.
fn may_be_a_short_decimal(value: f64) -> bool {
    let bits = value.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let (mantissa, exponent) = if biased == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1u64 << 52), biased - 1075)
    };
    if mantissa == 0 {
        return true;
    }
    let shift = mantissa.trailing_zeros() as i32;
    let power = exponent + shift;
    let odd_digits = (mantissa >> shift).to_string().len() as u32;
    if power >= 0 {
        return odd_digits * 1000 + power as u32 * 301 <= 19_000;
    }
    odd_digits * 1000 + (-power) as u32 * 699 <= 19_000
}

/// The double's exact decimal expansion, as `(digits, exp10)` with the value `0.<digits> * 10^exp10`.
fn exact_decimal(value: f64) -> (Vec<u8>, i32) {
    let bits = value.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    // Subnormals have no implicit leading one and a fixed exponent.
    let (mantissa, exponent) = if biased == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1u64 << 52), biased - 1075)
    };

    // Least significant digit first, which is the direction a carry travels.
    let mut digits: Vec<u8> = Vec::with_capacity(800);
    let mut rest = mantissa;
    if rest == 0 {
        digits.push(0);
    }
    while rest > 0 {
        digits.push((rest % 10) as u8);
        rest /= 10;
    }

    // `m * 2^e` for e >= 0 is repeated doubling; for e < 0 it is `m * 5^-e` with the decimal point
    // moved -e places left, because `2^e = 5^-e / 10^-e`.
    let (factor, shift) = if exponent >= 0 {
        (2u8, 0)
    } else {
        (5u8, -exponent)
    };
    let repetitions = if exponent >= 0 { exponent } else { -exponent };
    for _ in 0..repetitions {
        let mut carry = 0u8;
        for digit in digits.iter_mut() {
            let product = *digit * factor + carry;
            *digit = product % 10;
            carry = product / 10;
        }
        while carry > 0 {
            digits.push(carry % 10);
            carry /= 10;
        }
    }

    // Back to most significant first, and to the `0.<digits>` convention.
    digits.reverse();
    let exp10 = digits.len() as i32 - shift;
    (digits, exp10)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The values measured against the oracle, kept as unit tests so a refactor that breaks the
    /// rule fails here rather than three commits later in a conformance run.
    #[test]
    fn ties_follow_the_true_value_and_not_the_digit_string() {
        // The double is below the printed 5, so down, even though the digit before it is odd.
        assert_eq!(DEPTH_FORMAT.format(0.155), "0.15");
        // The double is above it, so up, even though the digit before it is even.
        assert_eq!(DEPTH_FORMAT.format(0.165), "0.17");
        assert_eq!(DEPTH_FORMAT.format(0.145), "0.14");
        assert_eq!(DEPTH_FORMAT.format(0.175), "0.17");
        assert_eq!(DEPTH_FORMAT.format(0.185), "0.18");
        assert_eq!(DEPTH_FORMAT.format(0.195), "0.2");
        assert_eq!(DEPTH_FORMAT.format(1.005), "1");
        assert_eq!(DEPTH_FORMAT.format(2.675), "2.67");
        assert_eq!(DEPTH_FORMAT.format(8.835), "8.84");
    }

    #[test]
    fn an_exact_tie_goes_to_the_even_neighbour() {
        // These are exactly representable, so the tie is real and HALF_EVEN applies.
        assert_eq!(DEPTH_FORMAT.format(0.125), "0.12");
        assert_eq!(DEPTH_FORMAT.format(0.375), "0.38");
        assert_eq!(DEPTH_FORMAT.format(0.625), "0.62");
        assert_eq!(DEPTH_FORMAT.format(0.875), "0.88");
        assert_eq!(DEPTH_FORMAT.format(1.125), "1.12");
        assert_eq!(FRACTION_FORMAT.format(0.015625), "0.0156");
    }

    /// The boundary the two patterns disagree on, which is the one place the pattern itself
    /// changes the rule rather than only the number of digits.
    #[test]
    fn the_underflow_tie_depends_on_the_pattern() {
        // Two fraction digits: the true value decides, and 0.005 is above the tie.
        assert_eq!(DEPTH_FORMAT.format(0.005), "0.01");
        // Four: a lone 5 is a tie with no preceding digit, and even means down — although this
        // double is also above the tie.
        assert_eq!(FRACTION_FORMAT.format(5e-5), "0");
        assert_eq!(FRACTION_FORMAT.format(1.5e-4), "0.0001");
        assert_eq!(FRACTION_FORMAT.format(2.5e-4), "0.0003");
        assert_eq!(FRACTION_FORMAT.format(3.5e-4), "0.0003");
    }

    #[test]
    fn the_shortest_form_is_what_gets_rounded() {
        // Beyond 2^53 the tail is zeros, not the double's true digits.
        assert!(DEPTH_FORMAT
            .format(-9.815_052_227_069_34e219)
            .starts_with("-9815052227069340000"));
        // And the shortest form can run out before the pattern does.
        assert_eq!(
            FRACTION_FORMAT.format(5.985_315_667_820_839e13),
            "59853156678208.39"
        );
    }

    #[test]
    fn zero_and_the_non_finite_symbols() {
        assert_eq!(DEPTH_FORMAT.format(0.0), "0");
        assert_eq!(DEPTH_FORMAT.format(-0.0), "-0");
        // A value too small to show keeps its sign.
        assert_eq!(DEPTH_FORMAT.format(-1e-9), "-0");
        assert_eq!(DEPTH_FORMAT.format(f64::NAN), "NaN");
        assert_eq!(DEPTH_FORMAT.format(f64::INFINITY), "∞");
        assert_eq!(DEPTH_FORMAT.format(f64::NEG_INFINITY), "-∞");
    }

    #[test]
    fn trailing_zeros_go_and_the_leading_zero_stays() {
        assert_eq!(DEPTH_FORMAT.format(1.0), "1");
        assert_eq!(DEPTH_FORMAT.format(0.5), "0.5");
        assert_eq!(DEPTH_FORMAT.format(10.0), "10");
        assert_eq!(DEPTH_FORMAT.format(1000.0), "1000");
        assert_eq!(DEPTH_FORMAT.format(0.999_9), "1");
        assert_eq!(FRACTION_FORMAT.format(0.999_9), "0.9999");
        assert_eq!(DEPTH_FORMAT.format(12_345.678_9), "12345.68");
        assert_eq!(FRACTION_FORMAT.format(1.0 / 3.0), "0.3333");
        assert_eq!(FRACTION_FORMAT.format(2.0 / 3.0), "0.6667");
    }

    /// An equidistant pair of shortest forms, which Rust resolves one way and the specification the
    /// other. This used to be a recorded divergence and is now a passing case.
    #[test]
    fn an_equidistant_shortest_form_goes_to_the_even_digit() {
        // 698583809467337.25 exactly, so `…337.2` and `…337.3` are equidistant and both
        // round-trip. Rust's formatter gives the odd one; the reference and the specification give
        // the even one.
        assert_eq!(
            DEPTH_FORMAT.format(6.985_838_094_673_373e14),
            "698583809467337.2"
        );
        assert_eq!(
            DEPTH_FORMAT.format(5.936_134_122_025_243e14),
            "593613412202524.2"
        );
    }

    /// The one limit that is left, and it is Java's rather than this port's.
    #[test]
    fn above_two_to_the_fifty_three_the_reference_leaves_the_shortest_form() {
        // The reference prints `-200971145216768832`, the double's exact value. This gives the
        // shortest form padded with zeros, and no rule of the kind above reaches the difference.
        assert_eq!(
            FRACTION_FORMAT.format(-2.009_711_452_167_688_3e17),
            "-200971145216768830"
        );
    }

    #[test]
    fn the_exact_expansion_is_exact() {
        // 0.1 is famously not 0.1, and the comparison has to see that.
        assert_eq!(compare_to_exact("1", 0, 0.1), Ordering::Greater);
        // 0.5 is exactly itself.
        assert_eq!(compare_to_exact("5", 0, 0.5), Ordering::Equal);
        // 1e23's shortest form is a whole decade away from the double it names.
        assert_eq!(compare_to_exact("1", 24, 1e23), Ordering::Less);
    }
}
