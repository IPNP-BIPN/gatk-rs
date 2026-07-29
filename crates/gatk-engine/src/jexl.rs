//! Ported from `org.apache.commons.jexl2` 2.1.1, the version GATK 4.6.2.0 pins: the `Parser`
//! number-literal rules, `JexlArithmetic`, and the `Interpreter` nodes an expression over read
//! attributes can reach.
//!
//! This exists for one filter, `JexlExpressionReadTagValueFilter`, but the semantics are not the
//! filter's: they are a general-purpose expression language's, and the surprising parts all come
//! from how it coerces types it was never given.
//!
//! # What the filter hands it
//!
//! Every attribute arrives as a **String**, because the filter's context calls
//! `read.getAttributeAsString(name)`. So `NM > 3` is a String against an Integer, and which of
//! `JexlArithmetic.compare`'s six branches runs is decided by the *literal's* type rather than by
//! the attribute's:
//!
//!  * `NM > 3`: `3` is an `Integer`, `isNumberable` is true, so both sides go through `toLong` and
//!    the comparison is integral. `toLong("30")` is `Long.parseLong`, which **throws** on a
//!    non-numeric attribute rather than answering false;
//!  * `NM > 3.5`: `3.5` is a `Float`, not a Double, because `setReal` defaults to `Float`. So
//!    `isFloatingPoint` is true and both sides go through `toDouble`, where an empty string is
//!    `NaN` and a non-numeric string throws;
//!  * `NM > '3'`: neither side is a number, so it is a **lexicographic** String comparison, and
//!    `"30" > "3"` is true for a different reason than `30 > 3` is.
//!
//! # The upstream bug this reproduces
//!
//! `JexlArithmetic.toLong(Object)` reads:
//!
//! ```java
//! } else if (val instanceof Double) {
//!     if (!Double.isNaN(((Double) val).doubleValue())) {
//!         return 0;
//!     } else {
//!         return ((Double) val).longValue();
//!     }
//! }
//! ```
//!
//! The test is inverted: a Double that is *not* NaN coerces to `0`, and a NaN one coerces to
//! `(long) NaN`, which is also `0`. Every Double is therefore `0L` in an integral context. It is
//! reproduced here rather than corrected, and pinned by the conformance suite.
//!
//! # What this refuses rather than guesses
//!
//! Arbitrary-precision literals (`1h`, `1.5b`) and results beyond `i128` are refused with
//! [`JexlError::Unsupported`]. A wrong answer would be indistinguishable from a right one in a
//! filter that returns a boolean; a refusal is not.

use std::collections::HashMap;

/// A JEXL value. The variants are the Java classes, because the class is what the arithmetic
/// dispatches on.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    /// `java.lang.Integer`.
    Int(i32),
    /// `java.lang.Long`.
    Long(i64),
    /// `java.lang.Float`. The default type of a real literal, which is why it is separate.
    Float(f32),
    /// `java.lang.Double`.
    Double(f64),
    Str(String),
}

/// What an evaluation can refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JexlError {
    /// The expression did not parse.
    Parse(String),
    /// `ArithmeticException`, which `setLenient(false)` makes the reference throw rather than
    /// swallow. Includes `NumberFormatException` from the coercions, which reaches the caller the
    /// same way.
    Arithmetic(String),
    /// `NumberFormatException`, which is a distinct variant because `add` catches it and
    /// concatenates while every other operator lets it out as an `ArithmeticException`.
    NumberFormat(String),
    /// A construct this port refuses rather than approximating.
    Unsupported(String),
}

impl Value {
    /// `JexlArithmetic.isFloatingPoint`: the *object* is a Float or a Double. A numeric-looking
    /// String is not, which is the whole reason `NM > 3` and `NM > 3.5` take different branches.
    fn is_floating_point(&self) -> bool {
        matches!(self, Value::Float(_) | Value::Double(_))
    }

    /// `JexlArithmetic.isNumberable`: Integer, Long, Byte, Short or Character. Not Float, not
    /// Double, not String.
    fn is_numberable(&self) -> bool {
        matches!(self, Value::Int(_) | Value::Long(_))
    }

    fn is_string(&self) -> bool {
        matches!(self, Value::Str(_))
    }
}

/// `JexlArithmetic.toBoolean` under `setLenient(false)`: a null operand throws.
fn to_boolean(value: &Value) -> Result<bool, JexlError> {
    match value {
        Value::Null => Err(JexlError::Arithmetic("null operand".into())),
        Value::Bool(b) => Ok(*b),
        Value::Int(i) => Ok(*i != 0),
        Value::Long(l) => Ok(*l != 0),
        Value::Float(f) => Ok(*f != 0.0),
        Value::Double(d) => Ok(*d != 0.0),
        Value::Str(s) => Ok(s == "true"),
    }
}

/// `JexlArithmetic.toDouble`.
fn to_double(value: &Value) -> Result<f64, JexlError> {
    match value {
        Value::Null => Err(JexlError::Arithmetic("null operand".into())),
        Value::Double(d) => Ok(*d),
        // `Double.parseDouble(String.valueOf(val))` rather than `doubleValue()`, so that a Float
        // goes through its own `toString` and widens the way the printed value reads rather than
        // the way the bits do: 6.4f becomes 6.4, not 6.400000095367432.
        Value::Float(f) => format_float(*f)
            .parse::<f64>()
            .map_err(|_| JexlError::Arithmetic(format!("Double coercion: {f}"))),
        Value::Int(i) => Ok(*i as f64),
        Value::Long(l) => Ok(*l as f64),
        Value::Bool(b) => Ok(if *b { 1.0 } else { 0.0 }),
        Value::Str(s) => {
            let trimmed = s.trim();
            if trimmed.is_empty() {
                Ok(f64::NAN)
            } else {
                parse_java_double(trimmed)
                    .ok_or_else(|| JexlError::NumberFormat(format!("For input string: \"{s}\"")))
            }
        }
    }
}

/// `JexlArithmetic.toLong`, inverted Double test included. See the module doc.
fn to_long(value: &Value) -> Result<i64, JexlError> {
    match value {
        Value::Null => Err(JexlError::Arithmetic("null operand".into())),
        // The upstream inversion: not-NaN yields 0, and NaN yields (long) NaN, which is 0 too.
        Value::Double(_) => Ok(0),
        Value::Float(f) => Ok(*f as i64),
        Value::Int(i) => Ok(*i as i64),
        Value::Long(l) => Ok(*l),
        Value::Bool(b) => Ok(if *b { 1 } else { 0 }),
        Value::Str(s) => {
            if s.is_empty() {
                Ok(0)
            } else {
                // `Long.parseLong`, which does not trim and does not accept a decimal point.
                parse_java_long(s)
                    .ok_or_else(|| JexlError::NumberFormat(format!("For input string: \"{s}\"")))
            }
        }
    }
}

/// `JexlArithmetic.toString`: a NaN Double prints as the empty string, everything else through
/// `String.valueOf`.
fn to_jexl_string(value: &Value) -> Result<String, JexlError> {
    match value {
        Value::Null => Err(JexlError::Arithmetic("null operand".into())),
        Value::Double(d) if d.is_nan() => Ok(String::new()),
        Value::Double(d) => Ok(format_double(*d)),
        Value::Float(f) => Ok(format_float(*f)),
        Value::Int(i) => Ok(i.to_string()),
        Value::Long(l) => Ok(l.to_string()),
        Value::Bool(b) => Ok(b.to_string()),
        Value::Str(s) => Ok(s.clone()),
    }
}

/// `Long.parseLong`: optional sign, then digits, no whitespace, no separators.
fn parse_java_long(text: &str) -> Option<i64> {
    let bytes = text.as_bytes();
    let (start, negative) = match bytes.first() {
        Some(b'-') => (1, true),
        Some(b'+') => (1, false),
        _ => (0, false),
    };
    if start >= bytes.len() {
        return None;
    }
    if !bytes[start..].iter().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let magnitude: i64 = text[start..].parse().ok()?;
    Some(if negative { -magnitude } else { magnitude })
}

/// `Double.parseDouble`, restricted to what a read attribute can hold: an optional sign, digits, a
/// decimal point and an exponent. Java also accepts hexadecimal floats and a trailing `d`/`f`,
/// which are refused here rather than mis-parsed.
fn parse_java_double(text: &str) -> Option<f64> {
    let ok = text
        .chars()
        .all(|c| c.is_ascii_digit() || matches!(c, '+' | '-' | '.' | 'e' | 'E'));
    if !ok {
        return None;
    }
    text.parse::<f64>().ok()
}

/// `Float.toString`. Java prints the shortest decimal that round-trips, which is what Rust's
/// `{}` does for `f32` as well.
fn format_float(value: f32) -> String {
    if value == value.trunc() && value.is_finite() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// `Double.toString`, same rule.
fn format_double(value: f64) -> String {
    if value == value.trunc() && value.is_finite() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// `JexlArithmetic.compare(left, right, operator)`: the branch order is the behaviour.
fn compare(left: &Value, right: &Value, operator: &str) -> Result<i32, JexlError> {
    if matches!(left, Value::Null) || matches!(right, Value::Null) {
        return Err(JexlError::Arithmetic(format!(
            "Object comparison:({left:?} {operator} {right:?})"
        )));
    }
    if left.is_floating_point() || right.is_floating_point() {
        let lhs = to_double(left)?;
        let rhs = to_double(right)?;
        // NaN sorts below everything, and two NaNs are equal. Java's own comparators say the same
        // and Rust's `partial_cmp` says neither, so the branch is written out.
        return Ok(if lhs.is_nan() {
            if rhs.is_nan() {
                0
            } else {
                -1
            }
        } else if rhs.is_nan() {
            1
        } else if lhs < rhs {
            -1
        } else if lhs > rhs {
            1
        } else {
            0
        });
    }
    if left.is_numberable() || right.is_numberable() {
        let lhs = to_long(left)?;
        let rhs = to_long(right)?;
        return Ok(match lhs.cmp(&rhs) {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Greater => 1,
            std::cmp::Ordering::Equal => 0,
        });
    }
    if left.is_string() || right.is_string() {
        // `String.compareTo`, which compares UTF-16 code units. For the ASCII a tag can hold this
        // is byte order.
        let lhs = to_jexl_string(left)?;
        let rhs = to_jexl_string(right)?;
        return Ok(match compare_java_strings(&lhs, &rhs) {
            d if d < 0 => -1,
            d if d > 0 => 1,
            _ => 0,
        });
    }
    if operator == "==" {
        return Ok(if left == right { 0 } else { -1 });
    }
    Err(JexlError::Arithmetic(format!(
        "Object comparison:({left:?} {operator} {right:?})"
    )))
}

/// `String.compareTo`: the difference of the first differing char, else of the lengths. The sign
/// is all the caller uses, but the reference returns the difference and so does this.
fn compare_java_strings(left: &str, right: &str) -> i32 {
    for (a, b) in left.chars().zip(right.chars()) {
        if a != b {
            return a as i32 - b as i32;
        }
    }
    left.chars().count() as i32 - right.chars().count() as i32
}

/// `JexlArithmetic.equals`: identity, then null, then Boolean, then `compare(..., "==")`.
fn jexl_equals(left: &Value, right: &Value) -> Result<bool, JexlError> {
    if left == right {
        return Ok(true);
    }
    if matches!(left, Value::Null) || matches!(right, Value::Null) {
        return Ok(false);
    }
    if matches!(left, Value::Bool(_)) || matches!(right, Value::Bool(_)) {
        return Ok(to_boolean(left)? == to_boolean(right)?);
    }
    Ok(compare(left, right, "==")? == 0)
}

/// `lessThan`, `greaterThan` and the two inclusive ones, whose null handling differs from
/// `equals`: a null operand is *false* rather than an error, and `left == right` short-circuits to
/// true for the inclusive pair and false for the strict pair.
fn relational(left: &Value, right: &Value, operator: &str) -> Result<bool, JexlError> {
    let inclusive = operator == "<=" || operator == ">=";
    if left == right {
        return Ok(inclusive);
    }
    if matches!(left, Value::Null) || matches!(right, Value::Null) {
        return Ok(false);
    }
    let ordering = compare(left, right, operator)?;
    Ok(match operator {
        "<" => ordering < 0,
        ">" => ordering > 0,
        "<=" => ordering <= 0,
        ">=" => ordering >= 0,
        _ => unreachable!("relational called with {operator}"),
    })
}

// ---------------------------------------------------------------------------------------------
// The expression tree, and the recursive-descent parser that builds it.
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Node {
    Literal(Value),
    Identifier(String),
    Not(Box<Node>),
    Negate(Box<Node>),
    And(Box<Node>, Box<Node>),
    Or(Box<Node>, Box<Node>),
    Binary(&'static str, Box<Node>, Box<Node>),
    /// `empty(x)`, a JEXL built-in rather than a method call.
    Empty(Box<Node>),
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(Value),
    Str(String),
    Ident(String),
    Symbol(&'static str),
}

fn tokenize(text: &str) -> Result<Vec<Token>, JexlError> {
    let bytes: Vec<char> = text.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        if c.is_ascii_digit() {
            let start = i;
            while i < bytes.len()
                && (bytes[i].is_ascii_digit()
                    || bytes[i] == '.'
                    || matches!(
                        bytes[i],
                        'x' | 'X'
                            | 'e'
                            | 'E'
                            | 'l'
                            | 'L'
                            | 'h'
                            | 'H'
                            | 'b'
                            | 'B'
                            | 'd'
                            | 'D'
                            | 'f'
                            | 'F'
                    )
                    || (bytes[i].is_ascii_hexdigit() && text[start..].starts_with("0x")))
            {
                i += 1;
            }
            let literal: String = bytes[start..i].iter().collect();
            tokens.push(Token::Number(parse_number_literal(&literal)?));
            continue;
        }
        if c == '\'' || c == '"' {
            let quote = c;
            i += 1;
            let mut value = String::new();
            while i < bytes.len() && bytes[i] != quote {
                // JEXL's own escape handling, which only ever unescapes the quote in use.
                if bytes[i] == '\\' && i + 1 < bytes.len() && bytes[i + 1] == quote {
                    value.push(quote);
                    i += 2;
                    continue;
                }
                value.push(bytes[i]);
                i += 1;
            }
            if i >= bytes.len() {
                return Err(JexlError::Parse(format!("unterminated string in {text}")));
            }
            i += 1;
            tokens.push(Token::Str(value));
            continue;
        }
        if c.is_alphabetic() || c == '_' {
            let start = i;
            while i < bytes.len() && (bytes[i].is_alphanumeric() || bytes[i] == '_') {
                i += 1;
            }
            tokens.push(Token::Ident(bytes[start..i].iter().collect()));
            continue;
        }
        let two: String = bytes[i..(i + 2).min(bytes.len())].iter().collect();
        let symbol: Option<&'static str> = match two.as_str() {
            "==" => Some("=="),
            "!=" => Some("!="),
            "<=" => Some("<="),
            ">=" => Some(">="),
            "&&" => Some("&&"),
            "||" => Some("||"),
            _ => None,
        };
        if let Some(symbol) = symbol {
            tokens.push(Token::Symbol(symbol));
            i += 2;
            continue;
        }
        let one: Option<&'static str> = match c {
            '<' => Some("<"),
            '>' => Some(">"),
            '!' => Some("!"),
            '+' => Some("+"),
            '-' => Some("-"),
            '*' => Some("*"),
            '/' => Some("/"),
            '%' => Some("%"),
            '(' => Some("("),
            ')' => Some(")"),
            _ => None,
        };
        match one {
            Some(symbol) => {
                tokens.push(Token::Symbol(symbol));
                i += 1;
            }
            None => return Err(JexlError::Parse(format!("unexpected character {c:?}"))),
        }
    }
    Ok(tokens)
}

/// `ASTNumberLiteral.setNatural` and `setReal`.
///
/// The default classes are the surprise: a whole number is an **Integer** (widening to Long then
/// BigInteger only if it does not fit), and a real is a **Float**, not a Double. Both decide which
/// branch of `compare` a comparison against a String attribute takes.
fn parse_number_literal(text: &str) -> Result<Value, JexlError> {
    let last = text.chars().last().unwrap_or('0');
    let real = text.contains('.')
        || ((text.contains('e') || text.contains('E')) && !text.starts_with("0x"));
    if real {
        let body = &text[..text.len() - 1];
        return match last {
            'b' | 'B' => Err(JexlError::Unsupported(format!(
                "BigDecimal literal {text}: this port refuses what it cannot reproduce exactly"
            ))),
            'd' | 'D' => body
                .parse::<f64>()
                .map(Value::Double)
                .map_err(|_| JexlError::Parse(format!("bad real literal {text}"))),
            'f' | 'F' => body
                .parse::<f32>()
                .map(Value::Float)
                .map_err(|_| JexlError::Parse(format!("bad real literal {text}"))),
            _ => text
                .parse::<f32>()
                .map(Value::Float)
                .map_err(|_| JexlError::Parse(format!("bad real literal {text}"))),
        };
    }
    let (body, base) =
        if let Some(hex) = text.strip_prefix("0x").or_else(|| text.strip_prefix("0X")) {
            (hex, 16)
        } else if text.len() > 1 && text.starts_with('0') {
            (&text[1..], 8)
        } else {
            (text, 10)
        };
    match last {
        'l' | 'L' => i64::from_str_radix(&body[..body.len() - 1], base)
            .map(Value::Long)
            .map_err(|_| JexlError::Parse(format!("bad long literal {text}"))),
        'h' | 'H' => Err(JexlError::Unsupported(format!(
            "BigInteger literal {text}: this port refuses what it cannot reproduce exactly"
        ))),
        _ => {
            if let Ok(value) = i32::from_str_radix(body, base) {
                Ok(Value::Int(value))
            } else if let Ok(value) = i64::from_str_radix(body, base) {
                Ok(Value::Long(value))
            } else {
                Err(JexlError::Unsupported(format!(
                    "integer literal {text} does not fit a Long"
                )))
            }
        }
    }
}

struct Parser {
    tokens: Vec<Token>,
    position: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.position)
    }

    fn eat_symbol(&mut self, symbol: &str) -> bool {
        if matches!(self.peek(), Some(Token::Symbol(s)) if *s == symbol) {
            self.position += 1;
            return true;
        }
        false
    }

    fn parse_or(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_and()?;
        while self.eat_symbol("||") {
            left = Node::Or(Box::new(left), Box::new(self.parse_and()?));
        }
        Ok(left)
    }

    fn parse_and(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_equality()?;
        while self.eat_symbol("&&") {
            left = Node::And(Box::new(left), Box::new(self.parse_equality()?));
        }
        Ok(left)
    }

    fn parse_equality(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_relational()?;
        loop {
            let operator = if self.eat_symbol("==") {
                "=="
            } else if self.eat_symbol("!=") {
                "!="
            } else {
                return Ok(left);
            };
            left = Node::Binary(operator, Box::new(left), Box::new(self.parse_relational()?));
        }
    }

    fn parse_relational(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_additive()?;
        loop {
            let operator = if self.eat_symbol("<=") {
                "<="
            } else if self.eat_symbol(">=") {
                ">="
            } else if self.eat_symbol("<") {
                "<"
            } else if self.eat_symbol(">") {
                ">"
            } else {
                return Ok(left);
            };
            left = Node::Binary(operator, Box::new(left), Box::new(self.parse_additive()?));
        }
    }

    fn parse_additive(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_multiplicative()?;
        loop {
            let operator = if self.eat_symbol("+") {
                "+"
            } else if self.eat_symbol("-") {
                "-"
            } else {
                return Ok(left);
            };
            left = Node::Binary(
                operator,
                Box::new(left),
                Box::new(self.parse_multiplicative()?),
            );
        }
    }

    fn parse_multiplicative(&mut self) -> Result<Node, JexlError> {
        let mut left = self.parse_unary()?;
        loop {
            let operator = if self.eat_symbol("*") {
                "*"
            } else if self.eat_symbol("/") {
                "/"
            } else if self.eat_symbol("%") {
                "%"
            } else {
                return Ok(left);
            };
            left = Node::Binary(operator, Box::new(left), Box::new(self.parse_unary()?));
        }
    }

    fn parse_unary(&mut self) -> Result<Node, JexlError> {
        if self.eat_symbol("!") {
            return Ok(Node::Not(Box::new(self.parse_unary()?)));
        }
        if self.eat_symbol("-") {
            return Ok(Node::Negate(Box::new(self.parse_unary()?)));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Node, JexlError> {
        match self.peek().cloned() {
            Some(Token::Number(value)) => {
                self.position += 1;
                Ok(Node::Literal(value))
            }
            Some(Token::Str(value)) => {
                self.position += 1;
                Ok(Node::Literal(Value::Str(value)))
            }
            Some(Token::Ident(name)) => {
                self.position += 1;
                match name.as_str() {
                    "true" => Ok(Node::Literal(Value::Bool(true))),
                    "false" => Ok(Node::Literal(Value::Bool(false))),
                    "null" => Ok(Node::Literal(Value::Null)),
                    "empty" if self.eat_symbol("(") => {
                        let inner = self.parse_or()?;
                        if !self.eat_symbol(")") {
                            return Err(JexlError::Parse("expected )".into()));
                        }
                        Ok(Node::Empty(Box::new(inner)))
                    }
                    _ => Ok(Node::Identifier(name)),
                }
            }
            Some(Token::Symbol("(")) => {
                self.position += 1;
                let inner = self.parse_or()?;
                if !self.eat_symbol(")") {
                    return Err(JexlError::Parse("expected )".into()));
                }
                Ok(inner)
            }
            other => Err(JexlError::Parse(format!("unexpected token {other:?}"))),
        }
    }
}

/// One compiled expression.
#[derive(Debug, Clone)]
pub struct Expression {
    root: Node,
}

/// `JexlEngine.createExpression`.
pub fn create_expression(text: &str) -> Result<Expression, JexlError> {
    let tokens = tokenize(text)?;
    let mut parser = Parser {
        tokens,
        position: 0,
    };
    let root = parser.parse_or()?;
    if parser.position != parser.tokens.len() {
        return Err(JexlError::Parse(format!("trailing tokens in {text}")));
    }
    Ok(Expression { root })
}

/// The attributes an expression can read. `get` returns `None` where the reference returns null,
/// which is what makes an absent tag an error under `setLenient(false)` rather than a false.
pub type Context = HashMap<String, String>;

impl Expression {
    /// `Expression.evaluate`.
    pub fn evaluate(&self, context: &Context) -> Result<Value, JexlError> {
        evaluate(&self.root, context)
    }
}

fn evaluate(node: &Node, context: &Context) -> Result<Value, JexlError> {
    match node {
        Node::Literal(value) => Ok(value.clone()),
        // `GATKReadJexlContext.get` returns `getAttributeAsString`, which is null for an absent
        // tag. Null is a value here and only becomes an error when an operator touches it.
        Node::Identifier(name) => Ok(match context.get(name) {
            Some(value) => Value::Str(value.clone()),
            None => Value::Null,
        }),
        Node::Not(inner) => Ok(Value::Bool(!to_boolean(&evaluate(inner, context)?)?)),
        Node::Negate(inner) => {
            let value = evaluate(inner, context)?;
            match value {
                Value::Int(i) => Ok(Value::Int(-i)),
                Value::Long(l) => Ok(Value::Long(-l)),
                Value::Float(f) => Ok(Value::Float(-f)),
                Value::Double(d) => Ok(Value::Double(-d)),
                other => Ok(Value::Long(-to_long(&other)?)),
            }
        }
        // Short-circuiting, as the Interpreter's ASTAndNode and ASTOrNode do.
        Node::And(left, right) => {
            if !to_boolean(&evaluate(left, context)?)? {
                return Ok(Value::Bool(false));
            }
            Ok(Value::Bool(to_boolean(&evaluate(right, context)?)?))
        }
        Node::Or(left, right) => {
            if to_boolean(&evaluate(left, context)?)? {
                return Ok(Value::Bool(true));
            }
            Ok(Value::Bool(to_boolean(&evaluate(right, context)?)?))
        }
        Node::Empty(inner) => {
            let value = evaluate(inner, context)?;
            Ok(Value::Bool(match value {
                Value::Null => true,
                Value::Str(s) => s.is_empty(),
                _ => false,
            }))
        }
        Node::Binary(operator, left, right) => {
            let left = evaluate(left, context)?;
            let right = evaluate(right, context)?;
            match *operator {
                "==" => Ok(Value::Bool(jexl_equals(&left, &right)?)),
                "!=" => Ok(Value::Bool(!jexl_equals(&left, &right)?)),
                "<" | ">" | "<=" | ">=" => Ok(Value::Bool(relational(&left, &right, operator)?)),
                "+" | "-" | "*" | "/" | "%" => arithmetic(operator, &left, &right),
                other => Err(JexlError::Unsupported(format!("operator {other}"))),
            }
        }
    }
}

/// `JexlArithmetic.add`, `subtract`, `multiply`, `divide` and `mod`.
///
/// Two things about their shape are the behaviour rather than the implementation:
///
///  * **the branch is `isFloatingPointNumber`, not `isFloatingPoint`.** That predicate is textual
///    for a String: any attribute containing `.`, `e` or `E` sends the operation down the double
///    path, so `NM + 1` is integral for `NM=30` and floating for `NM=3.0`;
///  * **only `add` concatenates, and only from a `catch`.** The reference tries the numeric path
///    and falls back to `toString(left).concat(toString(right))` when a coercion throws
///    `NumberFormatException`. Every other operator lets that exception out. So `NM + 1` on a
///    non-numeric tag is a string, and `NM - 1` on the same tag is an error.
///
/// The integral path goes through `BigInteger`, not `long`, so it does not wrap. This port carries
/// it in `i128` and refuses beyond that rather than wrapping silently.
fn arithmetic(operator: &str, left: &Value, right: &Value) -> Result<Value, JexlError> {
    if matches!(left, Value::Null) && matches!(right, Value::Null) {
        return Err(JexlError::Arithmetic("null operands".into()));
    }
    let numeric = numeric_arithmetic(operator, left, right);
    match (numeric, operator) {
        (Ok(value), _) => Ok(value),
        // The `catch (NumberFormatException)` that only `add` has.
        (Err(JexlError::NumberFormat(_)), "+") => Ok(Value::Str(format!(
            "{}{}",
            to_jexl_string(left)?,
            to_jexl_string(right)?
        ))),
        // Every other operator lets the NumberFormatException out unchanged, which is why
        // `NM + 1` on a non-numeric tag is a string and `NM - 1` on the same tag is an error.
        (Err(other), _) => Err(other),
    }
}

fn numeric_arithmetic(operator: &str, left: &Value, right: &Value) -> Result<Value, JexlError> {
    if is_floating_point_number(left) || is_floating_point_number(right) {
        let lhs = to_double(left)?;
        let rhs = to_double(right)?;
        let result = match operator {
            "+" => lhs + rhs,
            "-" => lhs - rhs,
            "*" => lhs * rhs,
            "/" => {
                if rhs == 0.0 {
                    return Err(JexlError::Arithmetic("/".into()));
                }
                lhs / rhs
            }
            "%" => {
                if rhs == 0.0 {
                    return Err(JexlError::Arithmetic("%".into()));
                }
                lhs % rhs
            }
            other => return Err(JexlError::Unsupported(format!("operator {other}"))),
        };
        return Ok(Value::Double(result));
    }
    let lhs = to_big_integer(left)?;
    let rhs = to_big_integer(right)?;
    let result = match operator {
        "+" => lhs.checked_add(rhs),
        "-" => lhs.checked_sub(rhs),
        "*" => lhs.checked_mul(rhs),
        "/" => {
            if rhs == 0 {
                return Err(JexlError::Arithmetic("/".into()));
            }
            lhs.checked_div(rhs)
        }
        "%" => {
            if rhs == 0 {
                return Err(JexlError::Arithmetic("%".into()));
            }
            lhs.checked_rem(rhs)
        }
        other => return Err(JexlError::Unsupported(format!("operator {other}"))),
    };
    let result = result.ok_or_else(|| {
        JexlError::Unsupported("result outside i128: BigInteger is unbounded upstream".into())
    })?;
    // `narrowBigInteger`: an Integer if it fits, else a Long, else the BigInteger itself.
    Ok(match i32::try_from(result) {
        Ok(narrowed) => Value::Int(narrowed),
        Err(_) => match i64::try_from(result) {
            Ok(narrowed) => Value::Long(narrowed),
            Err(_) => {
                return Err(JexlError::Unsupported(
                    "result needs a BigInteger, which this port refuses".into(),
                ))
            }
        },
    })
}

/// `JexlArithmetic.toBigInteger`.
///
/// Unlike `toLong`, this one has no inverted test: a non-NaN Double becomes `new
/// BigInteger(val.toString())`, which *throws* on anything with a decimal point, and a NaN Double
/// becomes zero. It also trims before testing for empty and then parses **untrimmed**, so a string
/// of spaces is zero while `" 1"` throws.
fn to_big_integer(value: &Value) -> Result<i128, JexlError> {
    match value {
        Value::Null => Err(JexlError::Arithmetic("null operand".into())),
        Value::Double(d) if d.is_nan() => Ok(0),
        Value::Double(d) => parse_big_integer(&format_double(*d)),
        Value::Float(f) => parse_big_integer(&format_float(*f)),
        Value::Int(i) => Ok(*i as i128),
        Value::Long(l) => Ok(*l as i128),
        Value::Bool(_) => Err(JexlError::Arithmetic("BigInteger coercion: Boolean".into())),
        Value::Str(s) => {
            if s.trim().is_empty() {
                Ok(0)
            } else {
                parse_big_integer(s)
            }
        }
    }
}

/// `new BigInteger(String)`: sign then digits, nothing else, and a `NumberFormatException`
/// otherwise. That exception is what `add` catches to concatenate.
fn parse_big_integer(text: &str) -> Result<i128, JexlError> {
    let bytes = text.as_bytes();
    let (start, negative) = match bytes.first() {
        Some(b'-') => (1, true),
        Some(b'+') => (1, false),
        _ => (0, false),
    };
    if start >= bytes.len() || !bytes[start..].iter().all(|b| b.is_ascii_digit()) {
        return Err(JexlError::NumberFormat(format!(
            "For input string: \"{text}\""
        )));
    }
    let magnitude: i128 = text[start..]
        .parse()
        .map_err(|_| JexlError::Unsupported(format!("{text} does not fit an i128")))?;
    Ok(if negative { -magnitude } else { magnitude })
}

/// `JexlArithmetic.isFloatingPointNumber`: a Float, a Double, or a String containing `.`, `e` or
/// `E`. Note that the String test is textual, so `"1e"` counts.
fn is_floating_point_number(value: &Value) -> bool {
    match value {
        Value::Float(_) | Value::Double(_) => true,
        Value::Str(s) => s.contains('.') || s.contains('e') || s.contains('E'),
        _ => false,
    }
}
