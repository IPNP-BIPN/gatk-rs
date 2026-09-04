//! Barclay's argument value model.
//!
//! Ported from `org.broadinstitute.barclay.argparser` (Barclay 5.0.0, the version GATK 4.6.2.0's
//! `build.gradle` pins), and from `joptsimple` (jopt-simple 5.0.3) for the grammar underneath it.
//!
//! This is the layer a covering-array vector is interpreted by. The unified CLI dispatcher is out
//! of scope, as ROADMAP G1.8 says; what is here is which vectors are accepted, what value each
//! field ends up holding, and which exception the rejected ones produce, because a port that
//! disagreed about any of those would be answering a different question from the one the array
//! asked.
//!
//! # Provenance
//!
//! Barclay is BSD 3-Clause (`LICENSE.txt` at tag `5.0.0` of `broadinstitute/barclay`, "Copyright
//! (c) 2009-2016, GATK Authors"), and jopt-simple 5.0.3 is MIT (the licence header of every file
//! in its sources jar). Both are permissive and compatible with this crate's Apache 2.0. They are
//! two libraries and the split matters: **the grammar is jopt-simple's, not Barclay's**, so the
//! rules about what a token is are cited to `joptsimple` below and the rules about what a value
//! means are cited to `barclay`.
//!
//! # `optional()` is not what decides whether an argument is optional
//!
//! ```java
//! this.defaultValueAsString = convertDefaultValueToString();
//! this.isOptional = argumentAnnotation.optional() || !this.defaultValueAsString.equals(NULL_ARGUMENT_STRING);
//! ```
//!
//! A field declared `optional = false` that was **initialised** is optional anyway, because its
//! default renders as something other than the string `"null"`. And `convertDefaultValueToString`
//! maps an *empty collection* back to `"null"`, so an initialised-but-empty `List` is required
//! while an initialised-and-non-empty one is not. Two fields that look identical in the
//! annotation differ by what the constructor left in them.
//!
//! # `"null"` is a value, and it means three different things
//!
//! On a collection it clears the collection, and warns rather than fails if it is not the first
//! value; on a non-optional argument it throws; and on a scalar whose **field** is primitive it
//! throws a different exception. Note which test the last one uses:
//! `getUnderlyingField().getType().isPrimitive()`, the raw field, not
//! `getUnderlyingFieldClass()`, which has already boxed `int` into `Integer`.
//!
//! # The bounds check calls a null out of range, and checks the recommended range against the hard one
//!
//! ```java
//! private boolean isValueOutOfRange(final Double value) {
//!     return value == null || getMinValue() != Double.NEGATIVE_INFINITY && value < getMinValue()
//!             || getMaxValue() != Double.POSITIVE_INFINITY && value > getMaxValue();
//! }
//! ```
//!
//! Two consequences, both measured in the golden. `--bounded-int null` is an
//! `OutOfRangeArgumentValue` rather than an accepted null. And `checkArgumentRange` uses this same
//! method for the **recommended** range, which compares against `minValue`/`maxValue`: for an
//! argument with a recommended range and no hard range those are infinities, so the warning can
//! only fire on a null; for an argument with both, the hard check has already thrown. The
//! recommended-range warning is close to unreachable, and this port reproduces that rather than
//! repairing it.
//!
//! The message tells the two apart in a way nothing documents: it formats the bounds as integers
//! only when the *value* is an `Integer`, so a rejected `0` reports `allowed range [1, 10]` and a
//! rejected `null` on the same argument reports `allowed range [1.0, 10.0]`.

use std::collections::BTreeSet;
use std::fmt::Write as _;

/// The string Barclay treats as the absence of a value, and as the rendering of an absent default.
///
/// `NULL_ARGUMENT_STRING`. It is a *value* on the command line, not a token the grammar knows.
pub const NULL_ARGUMENT_STRING: &str = "null";

/// A value a field can hold, in the shapes this model distinguishes.
///
/// The Java type is kept rather than flattened to a string for the same reason
/// `gatk_annotation::AnnotationValue` keeps it: the range check asks whether the value is a
/// `Number`, and the out-of-range message asks whether it is an `Integer`. Both questions have no
/// answer once everything is a string.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// A Java `null`. What a field holds before anything is assigned to it.
    Null,
    Str(String),
    Int(i32),
    /// A `Long`, which is a different field from an `Int` because its rendering is its own: a
    /// value past an int's range renders as itself rather than wrapping.
    Int64(i64),
    Double(f64),
    Bool(bool),
    /// An enum constant, by name. `Enum.valueOf` is case-**sensitive**, which the golden shows.
    Enum(String),
    /// A `Collection`, which is what an `@Argument` on a `List` field holds.
    List(Vec<Value>),
    /// A value whose type implements `TaggedArgument`: the value itself, plus the logical name and
    /// attributes `populateArgumentTags` wrote onto it.
    ///
    /// The tag is written **before** the value is stored, and it is written even when there is no
    /// tag: `populateArgumentTags(null)` sets the tag to null and the attributes to an empty map,
    /// so an untagged taggable argument is distinguishable from one nobody touched.
    Tagged {
        value: String,
        tag: Option<String>,
        attributes: Vec<(String, String)>,
    },
}

impl Value {
    /// `String.valueOf(value)`, which is how a field is rendered and how `AbstractCollection`
    /// renders its elements.
    pub fn to_java_string(&self) -> String {
        match self {
            Value::Null => "null".to_string(),
            Value::Str(text) => text.clone(),
            Value::Int(number) => number.to_string(),
            Value::Int64(number) => number.to_string(),
            Value::Double(number) => java_double_to_string(*number),
            Value::Bool(flag) => flag.to_string(),
            Value::Enum(name) => name.clone(),
            // The dump's `TaggedPath.toString` is the value it was constructed from; the tag
            // travels beside it rather than in it.
            Value::Tagged { value, .. } => value.clone(),
            // `AbstractCollection.toString`: square brackets, `", "` between elements.
            Value::List(values) => {
                let parts: Vec<String> = values.iter().map(Value::to_java_string).collect();
                format!("[{}]", parts.join(", "))
            }
        }
    }

    /// `((Number) value).doubleValue()`, or `None` where the reference has a null or a non-number.
    fn as_number(&self) -> Option<f64> {
        match self {
            Value::Int(number) => Some(f64::from(*number)),
            // A long widens to a double the way `Number.doubleValue()` widens it, which loses
            // precision past 2^53 exactly as the reference does.
            Value::Int64(number) => Some(*number as f64),
            Value::Double(number) => Some(*number),
            _ => None,
        }
    }
}

/// `Double.toString`, for the values this model reaches.
///
/// **Not the general algorithm.** Java's is specified as "the smallest number of digits that
/// uniquely distinguishes the argument value from adjacent values of type double", with a decimal
/// point always present and `E` notation outside `[1e-3, 1e7)`. Rust's `{}` for `f64` has the same
/// shortest-round-trip property, so the two agree on the finite non-scientific values a bound or a
/// bounded field can carry here, and this function adds Java's trailing `.0`. Anything outside
/// that range is refused rather than guessed: the general `Double.toString` is `FloatingDecimal`,
/// which is JDK source, and htsjdk-rs decision 0013 refused to transcribe it.
pub fn java_double_to_string(value: f64) -> String {
    assert!(
        value.is_finite() && (value == 0.0 || (1e-3..1e7).contains(&value.abs())),
        "java_double_to_string is ported only for finite values Java prints without an exponent; \
         {value} needs FloatingDecimal, which is not portable"
    );
    let rendered = format!("{value}");
    if rendered.contains('.') {
        rendered
    } else {
        // Java always prints "at least one digit ... after the decimal point".
        format!("{rendered}.0")
    }
}

/// What kind of value a field's *boxed* class accepts.
///
/// This is `getUnderlyingFieldClass()`, which boxes: a primitive `int` field reports `Integer`, so
/// it is a `Number` for the range check even though `getUnderlyingField().getType().isPrimitive()`
/// is what the null check asks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValueClass {
    Integer,
    Double,
    Text,
    Boolean,
    /// A type built from a `String` that implements `TaggedArgument`, which is the only thing
    /// `getValuePopulatedWithTags` asks of it.
    Tagged,
    /// An enum, with its simple name and its constants in declaration order: both appear in the
    /// message a bad value produces.
    Enum {
        simple_name: &'static str,
        constants: &'static [&'static str],
    },
    /// A `Float`, whose grammar is `Float.valueOf`'s and not `str::parse`'s, and whose refusal
    /// names `Float` where a double's names `Double`.
    Float,
    /// A `Long`, which is neither an `Integer` with a wider range nor a `Float` with a narrower
    /// grammar: `Long.valueOf` refuses a leading space where `Float.valueOf` trims one, and
    /// refuses a hexadecimal literal where the float takes `0x1p3`.
    Long,
    /// A class built from a `String` by a constructor that accepts every string.
    ///
    /// `File`, `GATKPath` and `FeatureInput` are all of this shape: a bad path is not a bad value,
    /// so no message exists for one. What separates them is the name their refusals would carry
    /// and whether they implement `TaggedArgument`, which is what decides whether a tag on the
    /// argument is a tag or an error.
    Constructed {
        simple_name: &'static str,
        taggable: bool,
    },
}

impl ValueClass {
    /// `Number.class.isAssignableFrom(getUnderlyingFieldClass())`.
    fn is_number(&self) -> bool {
        matches!(
            self,
            ValueClass::Integer | ValueClass::Double | ValueClass::Float | ValueClass::Long
        )
    }

    /// `getUnderlyingFieldClass().getSimpleName()`, which appears in two messages: the failure of
    /// the `String` constructor and the failure of `Enum.valueOf`.
    pub fn simple_name(&self) -> &'static str {
        match self {
            ValueClass::Integer => "Integer",
            ValueClass::Double => "Double",
            ValueClass::Text => "String",
            ValueClass::Boolean => "Boolean",
            ValueClass::Tagged => "TaggedPath",
            ValueClass::Enum { simple_name, .. } => simple_name,
            ValueClass::Float => "Float",
            ValueClass::Long => "Long",
            ValueClass::Constructed { simple_name, .. } => simple_name,
        }
    }

    /// `constructFromString`: `Enum.valueOf` for an enum, otherwise the `String` constructor.
    ///
    /// The failure of the `String` constructor arrives as an `InvocationTargetException` and is
    /// reported as "Failure constructing '<class>' from the string '<value>'."; the failure of
    /// `Enum.valueOf` is an `IllegalArgumentException` and is reported with the allowed values.
    fn construct_from_string(&self, text: &str, argument_name: &str) -> Result<Value, Error> {
        match self {
            ValueClass::Integer => text.parse::<i32>().map(Value::Int).map_err(|_| {
                Error::bad_argument_value_with_message(
                    argument_name,
                    text,
                    &format!("Failure constructing 'Integer' from the string '{text}'."),
                )
            }),
            ValueClass::Double => text.parse::<f64>().map(Value::Double).map_err(|_| {
                Error::bad_argument_value_with_message(
                    argument_name,
                    text,
                    &format!("Failure constructing 'Double' from the string '{text}'."),
                )
            }),
            ValueClass::Long => java_long(text).map(Value::Int64).ok_or_else(|| {
                Error::bad_argument_value_with_message(
                    argument_name,
                    text,
                    &format!("Failure constructing 'Long' from the string '{text}'."),
                )
            }),
            ValueClass::Float => java_float(text)
                .map(|value| Value::Double(f64::from(value)))
                .ok_or_else(|| {
                    Error::bad_argument_value_with_message(
                        argument_name,
                        text,
                        &format!("Failure constructing 'Float' from the string '{text}'."),
                    )
                }),
            ValueClass::Constructed { taggable, .. } => {
                if *taggable {
                    Ok(Value::Tagged {
                        value: text.to_string(),
                        tag: None,
                        attributes: Vec::new(),
                    })
                } else {
                    Ok(Value::Str(text.to_string()))
                }
            }
            ValueClass::Text => Ok(Value::Str(text.to_string())),
            ValueClass::Tagged => Ok(Value::Tagged {
                value: text.to_string(),
                tag: None,
                attributes: Vec::new(),
            }),
            // A `Boolean` field is a flag, and its values have already been through
            // `StrictBooleanConverter` in the grammar, so only "true" and "false" reach here.
            ValueClass::Boolean => Ok(Value::Bool(text == "true")),
            ValueClass::Enum {
                simple_name,
                constants,
            } => {
                if constants.contains(&text) {
                    Ok(Value::Enum(text.to_string()))
                } else {
                    // `getEnumOptions` renders "Possible values: {A, B}", and the caller wraps that
                    // in "Allowed values are %s", so the message says "Allowed values are Possible
                    // values: {...}" and ends in a space. Both are the reference's.
                    Err(Error::bad_argument_value_with_message(
                        argument_name,
                        text,
                        &format!(
                            "'{text}' is not a valid value for {simple_name}. Allowed values are Possible values: {{{}}} ",
                            constants.join(", ")
                        ),
                    ))
                }
            }
        }
    }
}

/// The `@Argument` annotation, as far as the value model reads it.
///
/// The bounds are `double` in the annotation whatever the field type is, and their defaults are the
/// infinities rather than absent, which is what `hasBoundedRange` tests for.
#[derive(Debug, Clone)]
pub struct Annotation {
    pub full_name: &'static str,
    pub short_name: &'static str,
    pub doc: &'static str,
    pub optional: bool,
    pub mutex: &'static [&'static str],
    pub min_value: f64,
    pub max_value: f64,
    pub min_recommended_value: f64,
    pub max_recommended_value: f64,
    /// `suppressFileExpansion()`, which the constructor refuses on anything but a collection.
    pub suppress_file_expansion: bool,
}

impl Default for Annotation {
    /// The annotation's own defaults, from `Argument.java`.
    fn default() -> Self {
        Self {
            full_name: "",
            short_name: "",
            doc: "Undocumented option",
            optional: false,
            mutex: &[],
            min_value: f64::NEG_INFINITY,
            max_value: f64::INFINITY,
            min_recommended_value: f64::NEG_INFINITY,
            max_recommended_value: f64::INFINITY,
            suppress_file_expansion: false,
        }
    }
}

/// A rejected command line: the Java exception class, and its message.
///
/// The class is carried because the reference's callers distinguish them (`MissingArgument` and
/// `BadArgumentValue` are both `CommandLineException`, and `OutOfRangeArgumentValue` is a
/// `BadArgumentValue`), and because a dump can report it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub class: &'static str,
    pub message: String,
}

const EXCEPTION: &str = "org.broadinstitute.barclay.argparser.CommandLineException";

impl Error {
    /// `new CommandLineException(msg)`.
    pub fn command_line(message: impl Into<String>) -> Self {
        Self {
            class: EXCEPTION,
            message: message.into(),
        }
    }

    /// `new CommandLineException.MissingArgument(arg, message)`.
    fn missing_argument(argument: &str, message: &str) -> Self {
        Self {
            class: "org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument",
            message: format!("Argument {argument} was missing: {message}"),
        }
    }

    /// `new CommandLineException.BadArgumentValue(message)`, the one-argument form, which prefixes
    /// rather than naming an argument.
    fn bad_argument_value(message: &str) -> Self {
        Self {
            class: "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue",
            message: format!("Illegal argument value: {message}"),
        }
    }

    /// `new CommandLineException.BadArgumentValue(arg, value, message)`.
    fn bad_argument_value_with_message(argument: &str, value: &str, message: &str) -> Self {
        Self {
            class: "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue",
            message: format!("Argument {argument} has a bad value: {value}. {message}"),
        }
    }

    /// `new CommandLineException.OutOfRangeArgumentValue(name, min, max, value)`, which is a
    /// `BadArgumentValue` with a range message built from the two bounds.
    fn out_of_range(name: &str, min: f64, max: f64, value: &Value) -> Self {
        // `getValueString`: a null renders as the four characters, not as an empty string.
        let value_string = value.to_java_string();
        // `asInt` is `value instanceof Integer`, so a **null** value formats the bounds as
        // doubles even on an Integer argument. That is the whole of the difference between
        // "allowed range [1, 10]." and "allowed range [1.0, 10.0].".
        let as_int = matches!(value, Value::Int(_));
        let render = |bound: f64| -> String {
            if as_int {
                // `Integer.toString((int) Math.rint(v))`.
                format!("{}", rint(bound) as i32)
            } else {
                java_double_to_string(bound)
            }
        };
        let has_min = min != f64::NEG_INFINITY;
        let has_max = max != f64::INFINITY;
        let range = if has_min && has_max {
            format!("allowed range [{}, {}].", render(min), render(max))
        } else if has_min {
            format!("minimum allowed value {}", render(min))
        } else if has_max {
            format!("maximum allowed value {}", render(max))
        } else {
            // `ShouldNeverReachHereException`: an unbounded range cannot produce this exception,
            // because `hasBoundedRange()` gates it.
            unreachable!("Unbounded range should never result in this exception")
        };
        Self {
            class:
                "org.broadinstitute.barclay.argparser.CommandLineException$OutOfRangeArgumentValue",
            message: format!("Argument {name} has a bad value: {value_string}. {range}"),
        }
    }
}

/// `Math.rint`: round half to **even**, unlike `round`.
fn rint(value: f64) -> f64 {
    let rounded = value.round();
    if (value - value.trunc()).abs() == 0.5 && rounded % 2.0 != 0.0 {
        rounded - value.signum()
    } else {
        rounded
    }
}

/// `Long.valueOf`, which is NOT `Float.valueOf` with a narrower grammar.
///
/// It refuses a leading space where the float trims one, refuses a hexadecimal literal where the
/// float takes `0x1p3`, and takes a leading plus like the float does. A value past a long's range
/// is a refusal and not a saturation, which is the other half of what separates it from the float.
pub fn java_long(text: &str) -> Option<i64> {
    if text.is_empty() {
        return None;
    }
    let (sign, digits) = match text.as_bytes()[0] {
        b'-' => (-1i128, &text[1..]),
        b'+' => (1i128, &text[1..]),
        _ => (1i128, text),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    let mut value: i128 = 0;
    for byte in digits.bytes() {
        value = value * 10 + i128::from(byte - b'0');
        if value > i128::from(u64::MAX) {
            return None;
        }
    }
    let signed = sign * value;
    if signed > i128::from(i64::MAX) || signed < i128::from(i64::MIN) {
        return None;
    }
    Some(signed as i64)
}

/// `Float.valueOf`, whose grammar is `FloatingDecimal.readJavaFormatString`'s.
///
/// It is not `str::parse`, and the two disagree in both directions. Java takes a trailing type
/// suffix (`1.5f`, `1.5d`), a hexadecimal literal with a binary exponent (`0x1p3` is eight),
/// leading and trailing whitespace and a leading plus; it spells its infinity and its not-a-number
/// with capitals, refusing `inf` and `nan` where Rust accepts both. All ten of those spellings are
/// in the `tool-argument-value-classes` golden.
///
/// A value out of a float's range is an infinity rather than a refusal, which is `Float.valueOf`'s
/// own rounding and not this function's.
pub fn java_float(text: &str) -> Option<f32> {
    // `readJavaFormatString` trims the string before looking at it, which is why a value with a
    // leading space parses and one with an underscore does not.
    let trimmed = text.trim_matches(|c: char| c.is_ascii_whitespace());
    if trimmed.is_empty() {
        return None;
    }
    let (sign, body) = match trimmed.as_bytes()[0] {
        b'-' => (-1.0f32, &trimmed[1..]),
        b'+' => (1.0f32, &trimmed[1..]),
        _ => (1.0f32, trimmed),
    };
    if body == "NaN" {
        return Some(f32::NAN);
    }
    if body == "Infinity" {
        return Some(sign * f32::INFINITY);
    }
    // The trailing type suffix is part of the grammar and is dropped before the digits are read.
    let body = match body.as_bytes().last() {
        Some(b'f' | b'F' | b'd' | b'D') if !body.eq_ignore_ascii_case("infinity") => {
            &body[..body.len() - 1]
        }
        _ => body,
    };
    if body.is_empty() {
        return None;
    }
    let lower = body.to_ascii_lowercase();
    if let Some(rest) = lower.strip_prefix("0x") {
        return java_hexadecimal_float(rest).map(|value| sign * value);
    }
    // What is left is a decimal literal, which Rust reads the same way EXCEPT that it also
    // accepts `inf` and `nan`: both are refused above by their spelling, and neither can reach
    // here, since a body of `inf` has no digit and `str::parse` is only asked for one that does.
    if !body.bytes().any(|byte| byte.is_ascii_digit()) {
        return None;
    }
    body.parse::<f32>().ok().map(|value| sign * value)
}

/// The hexadecimal half of the grammar: `0x<hex>[.<hex>]p<decimal exponent>`.
///
/// The binary exponent is REQUIRED, which is what separates `0x1p3` from `0x1`.
fn java_hexadecimal_float(text: &str) -> Option<f32> {
    let (digits, exponent) = text.split_once('p')?;
    let exponent: i32 = exponent.parse().ok()?;
    let (whole, fraction) = match digits.split_once('.') {
        Some((whole, fraction)) => (whole, fraction),
        None => (digits, ""),
    };
    if whole.is_empty() && fraction.is_empty() {
        return None;
    }
    let mut value = 0.0f64;
    for digit in whole.bytes() {
        value = value * 16.0 + f64::from((digit as char).to_digit(16)?);
    }
    let mut scale = 1.0f64 / 16.0;
    for digit in fraction.bytes() {
        value += f64::from((digit as char).to_digit(16)?) * scale;
        scale /= 16.0;
    }
    Some((value * 2.0f64.powi(exponent)) as f32)
}

/// `NamedArgumentDefinition`: one `@Argument` field, its rules, and its current value.
#[derive(Debug, Clone)]
pub struct Definition {
    pub annotation: Annotation,
    /// The field's declared name, used as the long name when `fullName` is empty.
    pub field_name: &'static str,
    /// `getUnderlyingFieldClass()`, which is boxed. For a collection it is the element class.
    pub class: ValueClass,
    /// Whether the field is a `Collection`.
    pub is_collection: bool,
    /// `getUnderlyingField().getType().isPrimitive()`, the **unboxed** question, which only the
    /// null check asks.
    pub field_is_primitive: bool,
    /// The value the constructor left in the field.
    pub value: Value,
    /// `getDescriptorForControllingPlugin()`, `None` for an argument the tool declares itself.
    ///
    /// A controlled argument exists in the parser whether or not anybody selected its plugin: the
    /// descriptor registers every implementation it discovers. What decides whether it is *usable*
    /// is [`Parser::validate_plugin_argument_values`].
    pub controlled_by: Option<PluginControl>,
    /// `defaultValueAsString`, computed once at construction from that initial value.
    default_value_as_string: String,
    /// `hasBeenSet`.
    has_been_set: bool,
}

impl Definition {
    /// The `NamedArgumentDefinition` constructor, as far as it concerns values.
    ///
    /// `defaultValueAsString` is computed **here**, from the initial value, and never again: a
    /// later assignment does not change whether the argument is optional.
    pub fn new(
        annotation: Annotation,
        field_name: &'static str,
        class: ValueClass,
        is_collection: bool,
        field_is_primitive: bool,
        initial: Value,
    ) -> Self {
        let default_value_as_string = convert_default_value_to_string(&initial, is_collection);
        Self {
            annotation,
            controlled_by: None,
            field_name,
            class,
            is_collection,
            field_is_primitive,
            value: initial,
            default_value_as_string,
            has_been_set: false,
        }
    }

    /// `getLongName()`: the full name, or the field's own name when the annotation has none.
    pub fn long_name(&self) -> &str {
        if self.annotation.full_name.is_empty() {
            self.field_name
        } else {
            self.annotation.full_name
        }
    }

    /// `getArgumentAliases()`: the short name first, if it exists, then the long name.
    pub fn argument_aliases(&self) -> Vec<&str> {
        let mut aliases = Vec::with_capacity(2);
        if !self.annotation.short_name.is_empty() {
            aliases.push(self.annotation.short_name);
        }
        aliases.push(self.long_name());
        aliases
    }

    /// `getArgumentAliasDisplayString()`, joined with `/`. This is what an error message names the
    /// argument by, and it is **not** the long name: it is `S/optional-string`.
    pub fn alias_display_string(&self) -> String {
        self.argument_aliases().join("/")
    }

    /// `getDefaultValueAsString()`.
    pub fn default_value_as_string(&self) -> &str {
        &self.default_value_as_string
    }

    /// `isOptional()`: the annotation **or** an initialised field. See the module note.
    pub fn is_optional(&self) -> bool {
        self.annotation.optional || self.default_value_as_string != NULL_ARGUMENT_STRING
    }

    /// `isFlag()`: a boolean-valued argument, which may appear with no value at all.
    pub fn is_flag(&self) -> bool {
        self.class == ValueClass::Boolean
    }

    /// `getHasBeenSet()`.
    pub fn has_been_set(&self) -> bool {
        self.has_been_set
    }

    /// `hasBoundedRange()`.
    pub fn has_bounded_range(&self) -> bool {
        self.annotation.min_value != f64::NEG_INFINITY || self.annotation.max_value != f64::INFINITY
    }

    /// `hasRecommendedRange()`.
    pub fn has_recommended_range(&self) -> bool {
        self.annotation.max_recommended_value != f64::INFINITY
            || self.annotation.min_recommended_value != f64::NEG_INFINITY
    }

    /// `isValueOutOfRange(value)`, including its first clause: a null is out of range.
    pub fn is_value_out_of_range(&self, value: Option<f64>) -> bool {
        match value {
            None => true,
            Some(number) => {
                (self.annotation.min_value != f64::NEG_INFINITY
                    && number < self.annotation.min_value)
                    || (self.annotation.max_value != f64::INFINITY
                        && number > self.annotation.max_value)
            }
        }
    }

    /// `checkArgumentRange(value)`.
    ///
    /// The recommended-range branch only warns, and this port keeps it as a returned flag rather
    /// than a log line: what it would print is `logger.warn`, which no golden can see, and what
    /// matters is that it does not throw. It is also nearly unreachable, for the reason in the
    /// module note.
    pub fn check_argument_range(&self, value: &Value) -> Result<bool, Error> {
        // "Only validate numeric types because we have already ensured at constructor time that
        // only numeric types have bounds."
        if !self.class.is_number() {
            return Ok(false);
        }
        let as_double = value.as_number();
        if self.has_bounded_range() && self.is_value_out_of_range(as_double) {
            return Err(Error::out_of_range(
                self.long_name(),
                self.annotation.min_value,
                self.annotation.max_value,
                value,
            ));
        }
        // The same test, against the same hard bounds. Not a typo here: a transcription.
        Ok(self.has_recommended_range() && self.is_value_out_of_range(as_double))
    }

    /// `setArgumentValues(parser, stream, preprocessedValues)`.
    ///
    /// `append` is the parser's `APPEND_TO_COLLECTIONS` option, whose default is "replace". It is
    /// the only thing that stops a collection from being cleared before the first value.
    pub fn set_argument_values(
        &mut self,
        values: &[String],
        append: bool,
        surrogates: &TagSurrogates,
        files: &dyn FileSource,
    ) -> Result<(), Error> {
        if self.is_collection {
            self.set_collection_values(values, append, surrogates, files)?;
        } else {
            self.set_scalar_value(values, surrogates)?;
        }
        self.has_been_set = true;
        Ok(())
    }

    /// `getNormalizedTagValuePair`: a value that is a surrogate key stands for a (tag, value) pair,
    /// and one that is not stands for itself with no tag.
    fn normalized_tag_value_pair(
        &self,
        text: &str,
        surrogates: &TagSurrogates,
    ) -> (Option<String>, String) {
        match surrogates.get(text) {
            Some((tag, value)) => (Some(tag.clone()), value.clone()),
            None => (None, text.to_string()),
        }
    }

    /// `getValuePopulatedWithTags`: build the value, then write the tag onto it — in that order,
    /// and onto **every** value an expansion file produced.
    ///
    /// A tag on a field whose type does not implement `TaggedArgument` is refused here rather than
    /// during preprocessing, and the message names the argument as `shortName/fullName`, so an
    /// argument with no short name reports a leading slash.
    fn value_populated_with_tags(&self, tag: Option<&str>, text: &str) -> Result<Value, Error> {
        let value = self.class.construct_from_string(text, self.long_name())?;
        match (&value, tag) {
            (Value::Tagged { value: raw, .. }, _) => {
                let (name, attributes) = match tag {
                    Some(tag_string) => {
                        let parsed = parse_tag(self.long_name(), tag_string)?;
                        (Some(parsed.0), parsed.1)
                    }
                    // `populateArgumentTags(null)` is still called: the tag is set to null and the
                    // attributes to an empty map.
                    None => (None, Vec::new()),
                };
                Ok(Value::Tagged {
                    value: raw.clone(),
                    tag: name,
                    attributes,
                })
            }
            (_, Some(tag_string)) => Err(Error::command_line(format!(
                "The argument: \"{}/{}\" does not accept tags: \"{tag_string}\"",
                self.annotation.short_name, self.annotation.full_name
            ))),
            (_, None) => Ok(value),
        }
    }

    fn set_collection_values(
        &mut self,
        values: &[String],
        append: bool,
        surrogates: &TagSurrogates,
        files: &dyn FileSource,
    ) -> Result<(), Error> {
        let mut collected: Vec<Value> = match (&self.value, append) {
            // "if this is a collection then we only want to clear it once at the beginning, before
            // we process any of the values, unless we're in APPEND_TO_COLLECTIONS mode". So the
            // field's declared contents are discarded the moment the user names the argument.
            (Value::List(existing), true) => existing.clone(),
            _ => Vec::new(),
        };

        for value in values {
            if value == NULL_ARGUMENT_STRING {
                // A "null" that is not the first value warns and clobbers; the warning is a
                // `logger.warn` and nothing observable, so only the clobbering is modelled.
                if !self.is_optional() {
                    return Err(Error::command_line(format!(
                        "Non \"null\" value must be provided for '{}'",
                        self.alias_display_string()
                    )));
                }
                collected.clear();
            } else {
                let (tag, raw) = self.normalized_tag_value_pair(value, surrogates);
                // `expandFromExpansionFile` is reached only from here, so expansion is a
                // **collection-only** mechanism: the identical value on a scalar is a value.
                let expanded = if self.annotation.suppress_file_expansion {
                    vec![raw]
                } else if EXPANSION_FILE_EXTENSIONS
                    .iter()
                    .any(|extension| raw.ends_with(extension))
                {
                    load_collection_list_file(&raw, files)?
                } else {
                    vec![raw]
                };
                for expanded_value in expanded {
                    let actual = self.value_populated_with_tags(tag.as_deref(), &expanded_value)?;
                    self.check_argument_range(&actual)?;
                    collected.push(actual);
                }
            }
        }
        self.value = Value::List(collected);
        Ok(())
    }

    fn set_scalar_value(
        &mut self,
        values: &[String],
        surrogates: &TagSurrogates,
    ) -> Result<(), Error> {
        // Two occurrences of a scalar are an error, not "the last one wins"; and so is one
        // occurrence carrying two values.
        if self.has_been_set || values.len() > 1 {
            return Err(Error::bad_argument_value(&format!(
                "Argument '{}' cannot be specified more than once.",
                self.alias_display_string()
            )));
        }

        if self.is_flag() && values.is_empty() {
            // A flag with no value is `true`. This is the only place a value appears from nowhere.
            self.value = Value::Bool(true);
            return Ok(());
        }

        let text = &values[0];
        let value = if text == NULL_ARGUMENT_STRING {
            // The **raw field**, not the boxed class: an `int` refuses, an `Integer` accepts.
            if self.field_is_primitive {
                return Err(Error::bad_argument_value(&format!(
                    "Argument '{}' is not a nullable argument type.",
                    self.alias_display_string()
                )));
            }
            Value::Null
        } else {
            let (tag, raw) = self.normalized_tag_value_pair(text, surrogates);
            self.value_populated_with_tags(tag.as_deref(), &raw)?
        };
        // Reached with the null too, which is what makes `--bounded-int null` out of range.
        self.check_argument_range(&value)?;
        self.value = value;
        Ok(())
    }

    /// `validateValues(parser)`, minus the plugin machinery: the mutex check, then the required
    /// check that a set mutex partner satisfies.
    pub fn validate_values(&self, provided_mutex_arguments: &[String]) -> Result<(), Error> {
        if self.has_been_set && !provided_mutex_arguments.is_empty() {
            return Err(Error::command_line(format!(
                "Argument '{}' cannot be used in conjunction with argument(s) {}",
                self.long_name(),
                provided_mutex_arguments.join(" ")
            )));
        }
        if !self.is_optional() {
            let missing = if self.is_collection {
                matches!(&self.value, Value::List(values) if values.is_empty())
            } else {
                !self.has_been_set
            };
            if missing && provided_mutex_arguments.is_empty() {
                return Err(Error::missing_argument(
                    self.long_name(),
                    &self.arg_required_error_message(),
                ));
            }
        }
        Ok(())
    }

    /// `getArgRequiredErrorMessage()`, whose second form renders the mutex targets as a
    /// `LinkedHashSet`'s `toString`.
    fn arg_required_error_message(&self) -> String {
        if self.annotation.mutex.is_empty() {
            format!("Argument '{}' is required", self.long_name())
        } else {
            format!(
                "Argument '{}' is required unless one of [{}] are provided",
                self.long_name(),
                self.annotation.mutex.join(", ")
            )
        }
    }
}

/// `convertDefaultValueToString`, which maps an **empty collection** to the same string as an
/// uninitialised field, and everything else to its `toString`.
fn convert_default_value_to_string(initial: &Value, is_collection: bool) -> String {
    match initial {
        Value::Null => NULL_ARGUMENT_STRING.to_string(),
        Value::List(values) if is_collection && values.is_empty() => {
            NULL_ARGUMENT_STRING.to_string()
        }
        other => other.to_java_string(),
    }
}

/// Where an expansion file's lines come from.
///
/// The reference opens a `FileReader` on the value. The port goes through a source so a
/// conformance suite can hand it the same bytes the dump wrote rather than rebuilding a file to
/// match a description of one; [`Filesystem`] is the behaviour the reference has.
pub trait FileSource {
    /// The file's contents, or [`IoError`] where the reference gets an `IOException`.
    fn read(&self, path: &str) -> Result<String, IoError>;
}

/// What the reference catches as an `IOException`. It carries nothing, because the message the
/// parser builds from it names the path rather than the cause.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IoError;

/// The default source: the filesystem, as `new FileReader(collectionListFile)` reads it.
pub struct Filesystem;

impl FileSource for Filesystem {
    fn read(&self, path: &str) -> Result<String, IoError> {
        std::fs::read_to_string(path).map_err(|_| IoError)
    }
}

/// `ARGUMENT_FILE_COMMENT`.
const ARGUMENT_FILE_COMMENT: &str = "#";

/// `SpecialArgumentsCollection.ARGUMENTS_FILE_FULLNAME`.
///
/// The parser does **not** declare it. A tool that wants argument files has to hold a
/// `SpecialArgumentsCollection`, which GATK's `CommandLineProgram` does; a tool that does not gets
/// "not a recognized option" for `--arguments_file`, which is the correct answer for it.
pub const ARGUMENTS_FILE_FULLNAME: &str = "arguments_file";

/// `EXPANSION_FILE_EXTENSIONS`.
const EXPANSION_FILE_EXTENSIONS: [&str; 2] = [".list", ".args"];

/// `loadCollectionListFile`: trim every line, drop the empty ones, drop the comments.
///
/// The `@`-prefixed warning is a `messageStream.println` and changes nothing about the result, so
/// it is not modelled: the file is expanded either way.
fn load_collection_list_file(path: &str, files: &dyn FileSource) -> Result<Vec<String>, Error> {
    match files.read(path) {
        Ok(text) => Ok(text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .filter(|line| !line.starts_with(ARGUMENT_FILE_COMMENT))
            .map(str::to_string)
            .collect()),
        // `new CommandLineException("I/O error loading list file:" + file, e)`. Note the missing
        // space after the colon, which is the reference's.
        Err(IoError) => Err(Error::command_line(format!(
            "I/O error loading list file:{path}"
        ))),
    }
}

/// `loadArgumentsFile`: every non-comment non-blank line, split on runs of whitespace.
///
/// `StringUtils.split(line)` collapses runs, so alignment inside a file is free — and an argument
/// value containing a space cannot be written in one at all.
fn load_arguments_file(path: &str, files: &dyn FileSource) -> Result<Vec<String>, Error> {
    match files.read(path) {
        Ok(text) => Ok(text
            .lines()
            // The comment test is on the **raw** line and the blank test is on the trimmed one, so
            // an indented `#` is not a comment while an indented empty line is blank.
            .filter(|line| !line.starts_with(ARGUMENT_FILE_COMMENT) && !line.trim().is_empty())
            .flat_map(|line| line.split_whitespace().map(str::to_string))
            .collect()),
        // Note the missing space after the colon, and that it is a different message from the
        // expansion-file one.
        Err(IoError) => Err(Error::command_line(format!(
            "I/O error loading arguments file:{path}"
        ))),
    }
}

/// `TaggedArgumentParser`: the rewrite that happens **before** jopt-simple sees the command line.
///
/// `--argument:logical_name,key=value raw_value` becomes `--argument` plus a *surrogate key*, and
/// the map remembers what the key stands for. Two things follow from the key's shape,
/// `option_string + ':' + raw_value`:
///
///  * it is the uniqueness test. `tagSurrogates.put` returning non-null is "duplicated on the
///    command line", so the same option with the same tag and the same value twice is an error
///    rather than two values, while the same tag with two different values is fine;
///  * the value jopt-simple parses is synthetic. Nothing downstream can look at it and see a path.
#[derive(Debug, Default)]
pub struct TagSurrogates {
    /// Insertion-ordered, because the only operation is "was this key already here".
    entries: Vec<(String, (String, String))>,
}

impl TagSurrogates {
    /// `getTaggedOptionForSurrogate`: the (tag, raw value) pair, or `None` for an ordinary value.
    pub fn get(&self, key: &str) -> Option<&(String, String)> {
        self.entries
            .iter()
            .find(|(candidate, _)| candidate == key)
            .map(|(_, pair)| pair)
    }

    fn put(&mut self, key: String, tag: String, value: String) -> Result<(), Error> {
        if self.get(&key).is_some() {
            // `rawOptionString` here is the option string **without** the prefix, which is why the
            // message reads `"tagged-collection:one a.bam"` and not `"--tagged-collection:one ..."`.
            let (option, raw) = key.rsplit_once(':').expect("the key carries its separator");
            return Err(Error::bad_argument_value(&format!(
                "The argument value: \"{option} {raw}\" was duplicated on the command line"
            )));
        }
        self.entries.push((key, (tag, value)));
        Ok(())
    }
}

/// `TaggedArgumentParser.USAGE`.
const TAG_USAGE: &str =
    "Tagged arguments must be of the form argument_name or argument_name:logical_name(,key=value)*";

/// `ParsedArgument.of(longArgName, rawTagValue)`: the logical name, then `key=value` pairs.
///
/// The three refusals are three *different* exception shapes, and two of them do not name the
/// argument at all: `BadArgumentValue("")` renders as `Argument  has a bad value:` with two
/// spaces, which the golden shows. Transcribed rather than tidied.
pub fn parse_tag(
    long_argument_name: &str,
    raw_tag_value: &str,
) -> Result<(String, Vec<(String, String)>), Error> {
    // `split(",", -1)`: the limit keeps trailing empty tokens, which is what makes a trailing
    // comma an "empty tag or attribute" rather than a silently ignored one.
    let tokens: Vec<&str> = raw_tag_value.split(',').collect();

    // "first token is required to be a name"
    if tokens[0].contains('=') {
        return Err(Error::bad_argument_value(&format!(
            "Missing tag name for argument: {raw_tag_value}"
        )));
    }
    if tokens.iter().any(|token| token.is_empty()) {
        return Err(Error::bad_argument_value_with_message(
            long_argument_name,
            raw_tag_value,
            &format!("Empty tag or attribute encountered. {TAG_USAGE}"),
        ));
    }

    let mut attributes: Vec<(String, String)> = Vec::new();
    for token in &tokens[1..] {
        let pair: Vec<&str> = token.split('=').collect();
        if pair.len() != 2 || pair[0].is_empty() || pair[1].is_empty() {
            // `BadArgumentValue("", rawTagValue, USAGE)`: the empty argument name is the
            // reference's, and it is what puts two spaces in the message.
            return Err(Error::bad_argument_value_with_message(
                "",
                raw_tag_value,
                TAG_USAGE,
            ));
        }
        if attributes.iter().any(|(key, _)| key == pair[0]) {
            return Err(Error::bad_argument_value_with_message(
                "",
                raw_tag_value,
                &format!("Duplicate key {}\n{TAG_USAGE}", pair[0]),
            ));
        }
        attributes.push((pair[0].to_string(), pair[1].to_string()));
    }
    Ok((tokens[0].to_string(), attributes))
}

/// `CommandLineParserUtilities.getAllFields` plus `createArgumentDefinitions`: how a nested
/// declaration becomes one flat list of arguments.
///
/// This is how `-L`, `-XL` and the read-filter arguments reach a tool. They are not declared on
/// the tool: they live on `@ArgumentCollection` objects it holds, and the parser flattens those
/// into a single namespace where nothing records which object an argument came from.
///
/// Two orderings are decided here and neither is stated anywhere in Barclay.
///
/// **Subclass before superclass.** `getAllFields` adds `clazz.getDeclaredFields()` and then climbs
/// to `getSuperclass()`, so a subclass's own fields are registered *first*. That order is the order
/// values are propagated in and the order `validateArgumentValues` reports a missing required
/// argument in, so which of two missing arguments a user is told about depends on which class
/// declared it.
///
/// **Depth-first, at the position of the field.** A collection declared between two `@Argument`
/// fields inserts all of its arguments between them.
pub enum FieldDecl {
    /// A field carrying `@Argument`.
    Argument(Box<Definition>),
    /// A field carrying `@ArgumentCollection`, holding the object's own declaration.
    Collection(ClassDecl),
    /// A field carrying `@ArgumentCollection` that the constructor left null. Not a value the
    /// parser can do anything with: it is a refusal, raised while the definitions are built.
    UninitialisedCollection { field_name: &'static str },
    /// A field carrying both `@Argument` and `@ArgumentCollection`, which is refused before either
    /// is looked at.
    BothAnnotations { field_name: &'static str },
}

/// One class's declaration: its own fields, in order, and the class it extends.
pub struct ClassDecl {
    /// The class's own name, which appears in the message an uninitialised collection produces.
    pub name: &'static str,
    /// `clazz.getDeclaredFields()`, in declaration order.
    pub declared: Vec<FieldDecl>,
    /// `clazz.getSuperclass()`, walked **after** the declared fields.
    pub superclass: Option<Box<ClassDecl>>,
}

impl ClassDecl {
    pub fn new(name: &'static str, declared: Vec<FieldDecl>) -> Self {
        Self {
            name,
            declared,
            superclass: None,
        }
    }

    /// The class this one extends. Its fields are registered after this class's own.
    pub fn extending(mut self, superclass: ClassDecl) -> Self {
        self.superclass = Some(Box::new(superclass));
        self
    }
}

impl Error {
    /// `new CommandLineException.CommandLineParserInternalException(msg)`. Raised while the
    /// definitions are built, so a tool with one of these is unusable before any command line
    /// exists.
    fn parser_internal(message: String) -> Self {
        Self {
            class: "org.broadinstitute.barclay.argparser.CommandLineException$CommandLineParserInternalException",
            message,
        }
    }
}

/// `createArgumentDefinitions`, flattened into the list the parser holds.
///
/// The duplicate check is `inArgumentMap`, which tests **every** alias of the new definition
/// against every alias already registered, so a short name colliding with somebody else's long
/// name is the same failure as two identical long names.
pub fn create_argument_definitions(class: &ClassDecl) -> Result<Vec<Definition>, Error> {
    let mut definitions: Vec<Definition> = Vec::new();
    collect(class, &mut definitions)?;
    Ok(definitions)
}

fn collect(class: &ClassDecl, into: &mut Vec<Definition>) -> Result<(), Error> {
    for field in &class.declared {
        match field {
            FieldDecl::BothAnnotations { field_name } => {
                return Err(Error::parser_internal(format!(
                    "Field {field_name}: Only one of @Argument, @ArgumentCollection or \
                     @PositionalArguments can be used"
                )))
            }
            FieldDecl::UninitialisedCollection { field_name } => {
                return Err(Error::parser_internal(format!(
                    "The ArgumentCollection field '{field_name}' in '{}' must have an initial value",
                    class.name
                )))
            }
            FieldDecl::Argument(definition) => {
                let aliases: Vec<String> = definition
                    .argument_aliases()
                    .iter()
                    .map(|alias| alias.to_string())
                    .collect();
                let clash = into.iter().any(|existing| {
                    existing
                        .argument_aliases()
                        .iter()
                        .any(|alias| aliases.iter().any(|new| new == alias))
                });
                if clash {
                    // The message names the alias display string, not the field, so two different
                    // fields in two different classes report the same text.
                    return Err(Error::parser_internal(format!(
                        "{} has already been used.",
                        definition.alias_display_string()
                    )));
                }
                into.push((**definition).clone());
            }
            // Depth-first, here rather than after the loop.
            FieldDecl::Collection(nested) => collect(nested, into)?,
        }
    }
    // The superclass last, which is what puts a base class's required argument after the
    // subclass's in every message that reports one.
    if let Some(superclass) = &class.superclass {
        collect(superclass, into)?;
    }
    Ok(())
}

/// What makes an argument a plugin's rather than a tool's.
///
/// The reference carries a whole `CommandLinePluginDescriptor` here. Two of its members are what
/// `validatePluginArgumentValues` actually reads, and they are the two kept:
///
///  * the **predecessor class**, whose `getSimpleName()` is what the error names — not the argument
///    that would have selected it, which is the surprising half;
///  * `isDependentArgumentAllowed`, which every descriptor in GATK answers by looking at the
///    argument that names plugins. It is modelled here as that argument's long name, so the
///    predicate is "the selector's values contain the predecessor's simple name".
///
/// The discovery half — `ClassFinder` scanning a package for implementations and instantiating
/// them — is not modelled. It decides *which* definitions exist, not what any of them means, and
/// it belongs with the tool that has plugins rather than with the argument layer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PluginControl {
    /// The plugin class's simple name.
    pub predecessor: &'static str,
    /// The long name of the argument that names plugins for this descriptor.
    pub selector: &'static str,
}

/// `validateAndResolvePlugins`, as GATK's read-filter descriptor implements it.
///
/// The interface method is the descriptor's own, and Barclay only promises to call it; what it
/// does is the descriptor's business. `GATKReadFilterPluginDescriptor` uses it to refuse a name
/// that matches no filter, which is why an unknown `--read-filter` is rejected there rather than
/// by the argument layer, and why the message is GATK's wording rather than Barclay's.
pub struct PluginResolution {
    /// The long name of the argument that names plugins.
    pub selector: &'static str,
    /// The names that resolve to something.
    pub known: Vec<String>,
    /// The message's prefix, up to the offending name. GATK's is
    /// `"Unrecognized read filter name: "`.
    pub unrecognized_prefix: &'static str,
}

/// `CommandLineArgumentParser` over a set of definitions.
pub struct Parser {
    definitions: Vec<Definition>,
    append_to_collections: bool,
    /// The descriptor's own `validateAndResolvePlugins`, if the caller has one.
    plugin_resolution: Option<PluginResolution>,
    /// The plugins the tool handed the descriptor as defaults, which count as selected.
    default_plugins: Vec<String>,
    /// `pluginDescriptor.validateAndResolvePlugins()`, which is the descriptor's OWN validation
    /// and not the parser's.
    ///
    /// It lives here rather than in the caller because of WHEN it runs: `validateArgumentValues`
    /// calls it before it walks the definitions, so a command line that both names a filter twice
    /// and breaks a mutex reports the filter. A port that ran it after parsing reported the mutex
    /// (#1070).
    #[allow(clippy::type_complexity)]
    plugin_validation: Option<Box<dyn Fn(&Parser) -> Result<(), Error>>>,
    /// `argumentsFilesLoadedAlready`.
    ///
    /// Parser state rather than per-call state, because the recursion is the same parser calling
    /// itself. Every file *named* in a pass goes in here, including ones that were skipped for
    /// already being present, which is what stops a file that includes itself and a pair of files
    /// that include each other.
    arguments_files_loaded_already: Vec<String>,
}

impl Parser {
    /// `new CommandLineArgumentParser(callerArguments)`, with the definitions supplied rather than
    /// discovered by reflection.
    pub fn new(definitions: Vec<Definition>) -> Self {
        Self {
            definitions,
            append_to_collections: false,
            plugin_resolution: None,
            default_plugins: Vec::new(),
            plugin_validation: None,
            arguments_files_loaded_already: Vec::new(),
        }
    }

    /// The descriptor's own `validateAndResolvePlugins`, which runs after the plugin trim and
    /// before the required check.
    pub fn with_plugin_resolution(mut self, resolution: PluginResolution) -> Self {
        self.plugin_resolution = Some(resolution);
        self
    }

    /// The descriptor's own validation, which runs where the reference runs it.
    pub fn with_plugin_validation(
        mut self,
        validate: impl Fn(&Parser) -> Result<(), Error> + 'static,
    ) -> Self {
        self.plugin_validation = Some(Box::new(validate));
        self
    }

    /// The plugin instances the tool handed the descriptor, whose arguments are allowed with no
    /// `--read-filter` on the command line.
    ///
    /// `GATKReadFilterPluginDescriptor` is constructed with the tool's `getDefaultReadFilters()`,
    /// and `isDependentArgumentAllowed` answers for a default the same way it answers for a filter
    /// the command line named. That is why a plain walker command line accepts a default filter's
    /// argument and refuses everybody else's, and it is measured either way in the
    /// `plugin-argument-ownership` golden.
    pub fn with_default_plugins(mut self, names: Vec<String>) -> Self {
        self.default_plugins = names;
        self
    }

    /// `CommandLineParserOptions.APPEND_TO_COLLECTIONS`, whose absence is "replace".
    pub fn with_append_to_collections(mut self) -> Self {
        self.append_to_collections = true;
        self
    }

    /// `new CommandLineArgumentParser(callerArguments)` over a nested declaration, which it
    /// flattens with [`create_argument_definitions`].
    pub fn from_class(class: &ClassDecl) -> Result<Self, Error> {
        Ok(Self::new(create_argument_definitions(class)?))
    }

    pub fn definitions(&self) -> &[Definition] {
        &self.definitions
    }

    /// The definition an alias names, if any. `namedArgumentsDefinitionsByAlias`.
    fn index_of_alias(&self, alias: &str) -> Option<usize> {
        self.definitions
            .iter()
            .position(|definition| definition.argument_aliases().contains(&alias))
    }

    /// `parseArguments(messageStream, args)`.
    ///
    /// Three phases, in this order, because the order decides which error a doubly-wrong command
    /// line reports: the grammar (jopt-simple), then the values (`propagateParsedValues`), then
    /// `validateArgumentValues`.
    pub fn parse_arguments(&mut self, argv: &[&str]) -> Result<(), Error> {
        self.parse_arguments_with(argv, &Filesystem)
    }

    /// The same, with the expansion files read from a source the caller names.
    ///
    /// The reference has no such seam: `loadCollectionListFile` opens a `FileReader`. It is here so
    /// a conformance suite can hand the port the bytes the dump wrote instead of rebuilding files
    /// to match a description of them, which is a second fixture.
    pub fn parse_arguments_with(
        &mut self,
        argv: &[&str],
        files: &dyn FileSource,
    ) -> Result<(), Error> {
        // `parser.parse(tagParser.preprocessTaggedOptions(args))`: the rewrite happens **first**,
        // so jopt-simple never sees a tag and every error about one comes from before or after it.
        let (preprocessed, surrogates) = self.preprocess_tagged_options(argv)?;
        let borrowed: Vec<&str> = preprocessed.iter().map(String::as_str).collect();
        let (parsed, positionals) = self.tokenize(&borrowed)?;

        // `--arguments_file` is checked here, between the grammar and the values, and it is the
        // only argument that changes the command line rather than a field.
        if let Some(expanded) = self.expand_from_argument_file(&parsed, files)? {
            // The file's arguments come **first**: `newArgs.addAll(Arrays.asList(args))` appends
            // the original command line to the expansion rather than the other way round, wherever
            // `--arguments_file` itself sat. That is why a collection ends up with the file's
            // values before the user's and a scalar given in both is a duplicate.
            let mut next: Vec<String> = expanded;
            next.extend(argv.iter().map(|arg| (*arg).to_string()));
            let borrowed: Vec<&str> = next.iter().map(String::as_str).collect();
            // The tag surrogates of this pass are dropped by construction: `preprocess_tagged_
            // options` builds a fresh map per call, which is what `resetTagSurrogates` does by
            // hand in the reference.
            return self.parse_arguments_with(&borrowed, files);
        }

        // `for (final OptionSpec<?> optSpec : parsedArguments.asMap().keySet())`: the map is over
        // the specs as they were registered, which is field order, and `has` filters it to the
        // ones actually given. So values are propagated in **declaration** order, not in the order
        // the user wrote them.
        for index in 0..self.definitions.len() {
            let Some(values) = parsed.iter().find(|(i, _)| *i == index).map(|(_, v)| v) else {
                continue;
            };
            let append = self.append_to_collections;
            self.definitions[index].set_argument_values(values, append, &surrogates, files)?;
        }

        if !positionals.is_empty() {
            // `stringValues.stream().collect(Collectors.joining("{", ",", "}"))`. The three
            // arguments of `joining` are (delimiter, prefix, suffix), so this is a delimiter of
            // `{`, a prefix of `,` and a suffix of `}`: one value renders as `,maybe}` and two as
            // `,a{b}`. Transcribed, not repaired.
            let mut joined = String::from(",");
            for (position, value) in positionals.iter().enumerate() {
                if position > 0 {
                    joined.push('{');
                }
                let _ = write!(joined, "{value}");
            }
            joined.push('}');
            return Err(Error::bad_argument_value(&format!(
                "Positional arguments were provided '{joined}' but no positional argument is \
                 defined for this tool."
            )));
        }

        self.validate_argument_values()
    }

    /// `expandFromArgumentFile`, plus the guard in `parseArguments` that decides whether to recurse.
    ///
    /// `Ok(None)` means "do not recurse", which covers both "the argument was not given" and "every
    /// file it named has already been read". The second is how the recursion terminates: the
    /// argument is **not** removed from the command line, so the next pass sees it again and
    /// expands it to nothing.
    fn expand_from_argument_file(
        &mut self,
        parsed: &[(usize, Vec<String>)],
        files: &dyn FileSource,
    ) -> Result<Option<Vec<String>>, Error> {
        let Some(index) = self.index_of_alias(ARGUMENTS_FILE_FULLNAME) else {
            return Ok(None);
        };
        let Some((_, argfiles)) = parsed.iter().find(|(i, _)| *i == index) else {
            return Ok(None);
        };

        // `.distinct().filter(not already loaded)`, in that order.
        let mut seen: Vec<&String> = Vec::new();
        let mut expanded: Vec<String> = Vec::new();
        for file in argfiles {
            if seen.contains(&file) {
                continue;
            }
            seen.push(file);
            if self.arguments_files_loaded_already.contains(file) {
                continue;
            }
            expanded.extend(load_arguments_file(file, files)?);
        }
        // Every file *named*, not every file read: a file skipped for being already loaded is
        // still recorded, which costs nothing and is what the reference does.
        for file in argfiles {
            self.arguments_files_loaded_already.push(file.clone());
        }

        if expanded.is_empty() {
            Ok(None)
        } else {
            Ok(Some(expanded))
        }
    }

    /// `TaggedArgumentParser.preprocessTaggedOptions`.
    ///
    /// Every token that looks like an option is inspected for a `:`. Without one, the token passes
    /// through and the `=` check runs on it. With one, the option name is split off, the following
    /// token is consumed as the value, and the pair is replaced by the bare option name and a
    /// surrogate key.
    ///
    /// The refusals here happen **before** anything knows which arguments exist, which is why
    /// `--:tumour` is "Zero length argument name" rather than "not a recognized option".
    #[allow(clippy::type_complexity)]
    fn preprocess_tagged_options(
        &self,
        argv: &[&str],
    ) -> Result<(Vec<String>, TagSurrogates), Error> {
        let mut out: Vec<String> = Vec::with_capacity(argv.len());
        let mut surrogates = TagSurrogates::default();
        let mut cursor = 0;

        while cursor < argv.len() {
            let token = argv[cursor];
            cursor += 1;
            let (prefix, rest) = if let Some(rest) = token.strip_prefix("--") {
                ("--", rest)
            } else if token.starts_with('-') && token != "-" {
                ("-", &token[1..])
            } else {
                out.push(token.to_string());
                continue;
            };

            let Some(separator) = rest.find(':') else {
                // `detectAndRejectHybridSyntax`, which is where the `=` refusal lives.
                reject_hybrid_syntax(rest)?;
                out.push(format!("{prefix}{rest}"));
                continue;
            };

            let option_name = &rest[..separator];
            reject_hybrid_syntax(option_name)?;
            let Some(value) = argv.get(cursor) else {
                return Err(Error::command_line(format!(
                    "No argument value found for tagged argument: {rest}"
                )));
            };
            if option_name.is_empty() {
                return Err(Error::command_line(format!(
                    "Zero length argument name found in tagged argument: {rest}"
                )));
            }
            let tag = &rest[separator + 1..];
            if tag.is_empty() {
                return Err(Error::command_line(format!(
                    "Zero length tag name found in tagged argument: {rest}"
                )));
            }
            if looks_like_an_option(value) {
                // The value slot holds another option, so there is nothing to consume. The message
                // is the same one the end-of-list case produces.
                return Err(Error::command_line(format!(
                    "No argument value found for tagged argument: {rest}"
                )));
            }
            cursor += 1;

            // `makeSurrogateKey`: the option string **without its prefix**, a colon, and the raw
            // value. It is never parsed; it is only compared.
            let key = format!("{rest}:{value}");
            surrogates.put(key.clone(), tag.to_string(), (*value).to_string())?;
            out.push(format!("{prefix}{option_name}"));
            out.push(key);
        }

        Ok((out, surrogates))
    }

    /// The grammar: jopt-simple 5.0.3 as `BarclayOptionParser` configures it, plus the one check
    /// Barclay does before the parser sees the tokens.
    ///
    /// Only the surface these definitions reach is ported. Short-argument **clustering** is off,
    /// which is the whole reason `BarclayOptionParser` exists; abbreviations are off, because
    /// Barclay constructs it with `allowAbbreviations = false`.
    #[allow(clippy::type_complexity)]
    fn tokenize(&self, argv: &[&str]) -> Result<(Vec<(usize, Vec<String>)>, Vec<String>), Error> {
        let mut given: Vec<(usize, Vec<String>)> = Vec::new();
        let mut positionals: Vec<String> = Vec::new();
        let mut cursor = 0;

        while cursor < argv.len() {
            let token = argv[cursor];
            let name = if let Some(rest) = token.strip_prefix("--") {
                rest
            } else if token.len() > 1 && token.starts_with('-') && !looks_like_number(token) {
                &token[1..]
            } else {
                positionals.push(token.to_string());
                cursor += 1;
                continue;
            };

            // The `=` refusal has already run in `preprocess_tagged_options`, which is where the
            // reference puts it; repeating it here would be a second implementation of one rule.

            let Some(index) = self.index_of_alias(name) else {
                // jopt-simple's `UnrecognizedOptionException`, whose message Barclay re-wraps in a
                // plain `CommandLineException`.
                return Err(Error::command_line(format!(
                    "{name} is not a recognized option"
                )));
            };
            let definition = &self.definitions[index];
            cursor += 1;

            let value = if definition.is_flag() {
                // `withOptionalArg().withValuesConvertedBy(new StrictBooleanConverter())`, and
                // `OptionalArgumentOptionSpec.detectOptionArgument`: the next token is taken only
                // if it does not look like an option **and** the converter accepts it. So
                // `--flag maybe` leaves "maybe" on the command line as a positional argument
                // rather than reporting a bad boolean.
                match argv.get(cursor) {
                    Some(next) if !looks_like_an_option(next) => match strict_boolean(next) {
                        Some(converted) => {
                            cursor += 1;
                            Some(converted)
                        }
                        None => None,
                    },
                    _ => None,
                }
            } else {
                // `withRequiredArg()`: the next token is the value whatever it looks like, and its
                // absence is jopt-simple's `OptionMissingRequiredArgumentException`.
                match argv.get(cursor) {
                    Some(next) => {
                        cursor += 1;
                        Some((*next).to_string())
                    }
                    None => {
                        return Err(Error::command_line(format!(
                            "Option {} requires an argument",
                            definition.alias_display_string()
                        )))
                    }
                }
            };

            match given.iter_mut().find(|(i, _)| *i == index) {
                Some((_, values)) => values.extend(value),
                None => given.push((index, value.into_iter().collect())),
            }
        }

        Ok((given, positionals))
    }

    /// `isDependentArgumentAllowed`: did the command line name this plugin?
    /// Whether a plugin the parser knows about is SELECTED, by a command line or by the tool's
    /// own defaults. Public because the expanded command line asks the same question.
    pub fn plugin_is_selected(&self, control: &PluginControl) -> bool {
        self.is_dependent_argument_allowed(control)
    }

    fn is_dependent_argument_allowed(&self, control: &PluginControl) -> bool {
        // A default counts as named: the descriptor was constructed with the tool's own filters,
        // and it does not distinguish those from the ones the command line asked for.
        if self
            .default_plugins
            .iter()
            .any(|name| name == control.predecessor)
        {
            return true;
        }
        match self.value_of(control.selector) {
            Some(Value::List(values)) => values
                .iter()
                .any(|value| value.to_java_string() == control.predecessor),
            _ => false,
        }
    }

    /// `validatePluginArgumentValues()`, which runs **before** the required check and rewrites the
    /// list it checks.
    ///
    /// Three cases, and the first is the one worth having: a controlled argument whose plugin
    /// nobody selected and which nobody gave is **dropped from the list**. Not made optional —
    /// removed, so a *required* argument belonging to an unselected plugin does not fire. Every
    /// GATK read filter has arguments in that state on every run.
    fn validate_plugin_argument_values(&self) -> Result<Vec<&Definition>, Error> {
        let mut actual: Vec<&Definition> = Vec::new();
        for definition in &self.definitions {
            let Some(control) = &definition.controlled_by else {
                actual.push(definition);
                continue;
            };
            let allowed = self.is_dependent_argument_allowed(control);
            if definition.has_been_set() {
                if !allowed {
                    // The message names the predecessor **class**, not the argument that selects
                    // it, and builds the argument's name as `getShortName() + "/" + getLongName()`
                    // with no guard, so one with no short name reports a leading slash.
                    return Err(Error::command_line(format!(
                        "Argument \"{}/{}\" is only valid when the argument \"{}\" is specified",
                        definition.annotation.short_name,
                        definition.long_name(),
                        control.predecessor
                    )));
                }
                actual.push(definition);
            } else if allowed {
                // Kept so the required check can still see it: a selected plugin's required
                // argument does fire.
                actual.push(definition);
            }
        }
        // "finally, give each plugin a chance to trim down any unseen instances from its own
        // list": after the trim, before the required check.
        if let Some(resolution) = &self.plugin_resolution {
            if let Some(Value::List(values)) = self.value_of(resolution.selector) {
                for value in values {
                    let name = value.to_java_string();
                    if !resolution.known.contains(&name) {
                        return Err(Error::command_line(format!(
                            "{}{name}",
                            resolution.unrecognized_prefix
                        )));
                    }
                }
            }
        }
        Ok(actual)
    }

    /// `validateArgumentValues()`: every definition, in declaration order, so the first missing
    /// required argument reported is the first one declared.
    fn validate_argument_values(&self) -> Result<(), Error> {
        let definitions = self.validate_plugin_argument_values()?;
        // The descriptor's own validation, BEFORE the per-definition walk: this is where the
        // reference decides that a read filter is both enabled and disabled, and a command line
        // that also breaks a mutex reports the filter rather than the mutex.
        if let Some(validate) = &self.plugin_validation {
            validate(self)?;
        }
        for definition in definitions {
            // The partners that were actually given, by long name, in the annotation's order.
            let provided: Vec<String> = definition
                .annotation
                .mutex
                .iter()
                .filter_map(|target| self.index_of_alias(target))
                .filter(|index| self.definitions[*index].has_been_set())
                .map(|index| self.definitions[index].long_name().to_string())
                .collect();
            definition.validate_values(&provided)?;
        }
        Ok(())
    }

    /// The value a field ended up holding, by long name.
    pub fn value_of(&self, long_name: &str) -> Option<&Value> {
        self.definitions
            .iter()
            .find(|definition| definition.long_name() == long_name)
            .map(|definition| &definition.value)
    }
}

/// `detectAndRejectHybridSyntax`: an option **name** may not contain `=`.
///
/// The comment above it in the reference explains why this is a refusal rather than a convenience:
/// jopt-simple would accept `-O=value`, but a value containing a `:` would then be read as tagging
/// syntax and fail somewhere else, so the spelling is refused everywhere instead of working
/// sometimes.
fn reject_hybrid_syntax(option_name: &str) -> Result<(), Error> {
    if option_name.contains('=') {
        return Err(Error::command_line(format!(
            "Can't parse option name containing an embedded '=' ({option_name})"
        )));
    }
    Ok(())
}

/// `OptionParser.looksLikeAnOption`: a token that starts with a dash and is longer than one
/// character.
fn looks_like_an_option(token: &str) -> bool {
    token.len() > 1 && token.starts_with('-')
}

/// A lone `-` is not an option, and neither is a negative number: jopt-simple's short-option
/// recognition is what decides, and Barclay disables clustering so a single dash introduces one
/// name.
fn looks_like_number(token: &str) -> bool {
    token[1..].parse::<f64>().is_ok()
}

/// `StrictBooleanConverter.convert`, which is case-insensitive and accepts the single letters.
///
/// Its rejection is a `ValueConversionException`, which in the optional-argument path is not an
/// error at all: it only means the token was not this option's value.
fn strict_boolean(value: &str) -> Option<String> {
    let lower = value.to_ascii_lowercase();
    match lower.as_str() {
        "true" | "t" => Some("true".to_string()),
        "false" | "f" => Some("false".to_string()),
        _ => None,
    }
}

/// The long names of a set of definitions, sorted, for a caller that wants to report on them.
pub fn long_names(definitions: &[Definition]) -> BTreeSet<String> {
    definitions
        .iter()
        .map(|definition| definition.long_name().to_string())
        .collect()
}
