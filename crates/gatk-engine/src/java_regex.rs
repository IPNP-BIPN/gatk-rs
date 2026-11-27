//! A `java.util.regex` subset, enough for the patterns GATK compiles from command lines.
//!
//! **Pinned by measurement, not by source.** `java.util.regex` is JDK code (GPL2) and is never
//! quoted here. The behaviour below is pinned by the golden `tests/data/java_regex.txt.gz`,
//! written by `tools/readfilter-conformance/JavaRegexDump.java`, which records for each pattern
//! whether it compiles, its refusal's description, index and message, and what `find()` answers
//! on a fixed set of inputs. `tests/java_regex.rs` replays every row, and its `BOUNDARY` lists are
//! the patterns the reference compiles and this refuses.
//!
//! GATK hands user strings to `Pattern.compile` in several places and then calls **`find()`**, not
//! `matches()`: `Utils.filterCollectionByExpressions`, which is what `-se` and `-xl-se` go through,
//! and `ClipReads`' sequence clipper. An unanchored search is what those arguments mean, and it is
//! why `-se s1` selects a sample named `xs10`. Reproducing the selection means reproducing the
//! search, so the matcher is here rather than assumed away. The golden's `find  s1  xs10  true`
//! row is that selection.
//!
//! # What this supports, and what it does not
//!
//! Literals, `.`, `^`, `$`, character classes with ranges and negation, the `\d \D \w \W \s \S`
//! shorthands, escaped metacharacters, groups, alternation, and the `* + ? {n} {n,} {n,m}`
//! quantifiers in their greedy, reluctant and possessive forms. That is the whole of what a sample
//! expression or a base sequence can usefully be.
//!
//! It does **not** support back references, look-around, named groups, POSIX or Unicode property
//! classes, class intersection, flags, or `\b` and the other boundary matchers. Compiling one of
//! those is an error here rather than a silent mismatch, which is the only safe way to be a subset:
//! a pattern this cannot represent must not quietly match differently from the reference.
//!
//! # Where the reference's own behaviour shows through
//!
//!  * matching is **backtracking**, left to right, first-alternative-first, exactly as Java's is,
//!    so a greedy quantifier gives back one repetition at a time and alternation prefers the
//!    branch written first. Both are observable through `find()`'s leftmost result;
//!  * `.` does not match a line terminator, `^` matches only at the start of the input, and `$`
//!    matches at the end **or before a final line terminator**, which is what Java does without
//!    `MULTILINE`;
//!  * a repetition that consumes nothing stops the loop rather than spinning, which is Java's own
//!    guard and what keeps `(a*)*` from hanging;
//!  * and the compile failure carries the message Java's does, `<description> near index <i>`
//!    followed by the pattern and a caret under that index. `Pattern.compile("[")` throws
//!    `Unclosed character class near index 0`, and SelectVariants lets it out unwrapped, so the
//!    text reaches the user as the regex engine wrote it. The descriptions produced are the ones
//!    listed on [`PatternSyntaxError::description`]; the golden's `compile ... error` rows pin
//!    each of them but the two this port adds for constructs outside the subset (and the `'+'`
//!    and `'?'` forms of the dangling-meta one).
//!
//! # Reading this file without reading Rust
//!
//! - `enum Node` is a tree of the compiled pattern: each value is one of the listed variants, and
//!   a variant can carry fields (`Repeat` carries the repeated node and its bounds). `Box<Node>`
//!   is a node stored behind a pointer, which a recursive type needs; Java's object references
//!   are always that.
//! - `match` picks the arm for the variant it is given, and the compiler refuses a `match` that
//!   misses one, so every node kind is handled in `match_node`.
//! - `&dyn Fn(usize) -> bool` is a function value, Java's `IntPredicate`. The matcher passes "the
//!   rest of the match" down as one (a continuation), which is how backtracking is written here
//!   without an explicit stack: a node calls it once per place it could end, and stops at the
//!   first that answers `true`.
//! - `Option<char>` is "a character or nothing"; `let Some(c) = x else { ... };` takes the
//!   character or leaves early. `Result<T, PatternSyntaxError>` is "a value or the exception", and
//!   `?` after a call passes the exception up, Java's `throws` at one call site.
//! - `char` is a Unicode code point, not Java's UTF-16 `char`. The input is split into code points
//!   with `.chars()`, so a character outside the Basic Multilingual Plane is one position here
//!   where a Java string index counts two; see [`PatternSyntaxError::index`].

/// A compiled pattern: `Pattern.compile`'s result, without flags.
#[derive(Debug, Clone)]
pub struct Pattern {
    /// The parsed tree.
    node: Node,
    /// The text compiled, `Pattern.pattern()`.
    source: String,
}

/// What `Pattern.compile` throws, with the three parts Java's `getMessage` joins with newlines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatternSyntaxError {
    /// The description, `getDescription()`. The ones the reference was measured to produce, each
    /// with a golden row: `Unclosed character class`, `Unclosed group`, `Unmatched closing ')'`,
    /// `Dangling meta character '*'` (the `'+'` and `'?'` forms share its code and have no row),
    /// `Illegal repetition`, `Illegal repetition range`, `Unclosed counted closure`,
    /// `Unexpected internal error` (a trailing backslash), `Illegal character range` and
    /// `Unknown inline modifier`.
    ///
    /// Two are this port's own, for constructs it refuses: `Unsupported escape sequence` (a
    /// back reference, an octal or control escape, `\b` and every other letter escape it does
    /// not know) and `Unsupported class union` (a `[` inside a class). Where the reference
    /// compiles the pattern, `tests/java_regex.rs` lists it in `BOUNDARY`; where it refuses it
    /// too, its description has not been measured and may differ.
    pub description: String,
    /// The index the caret points at, `getIndex()`. Negative where Java has nothing to point at:
    /// a bare `)` reports -1.
    ///
    /// Counted in code points (the parser reads `.chars()`), where a Java string index counts
    /// UTF-16 units. No golden row has a character outside the Basic Multilingual Plane, so the
    /// two are not compared where they could differ.
    pub index: i32,
    /// `getPattern()`: the text that was compiled.
    pub pattern: String,
}

impl PatternSyntaxError {
    /// `getMessage()`, whose three parts are each conditional.
    ///
    /// Java appends the "near index" clause only for a non-negative index, and the caret line only
    /// when that index is inside the pattern. A bare `)` reports -1 and gets neither; a trailing
    /// backslash reports the length and gets the clause without the caret.
    pub fn message(&self) -> String {
        let mut out = self.description.clone();
        if self.index >= 0 {
            out.push_str(&format!(" near index {}", self.index));
        }
        out.push('\n');
        out.push_str(&self.pattern);
        let length = self.pattern.chars().count() as i32;
        if self.index >= 0 && self.index < length {
            out.push('\n');
            out.push_str(&" ".repeat(self.index as usize));
            out.push('^');
        }
        out
    }

    /// The exception's class name, for callers that report a refusal as the reference names it.
    pub fn java_class(&self) -> &'static str {
        "java.util.regex.PatternSyntaxException"
    }
}

/// One node of the compiled pattern.
#[derive(Debug, Clone)]
enum Node {
    /// The branches of a `|`, tried in the order they were written.
    Alternation(Vec<Node>),
    /// Nodes matched one after another. Empty matches the empty string, which is what `a|` and
    /// `{2}` compile to in part (both are golden rows that compile).
    Sequence(Vec<Node>),
    /// A quantified node. `max` is `usize::MAX` for an unbounded `*`, `+` or `{n,}`.
    Repeat {
        node: Box<Node>,
        min: usize,
        max: usize,
        kind: RepeatKind,
    },
    /// One character, compared exactly (no case folding: no flags are supported).
    Literal(char),
    /// `.`, which is every character but a line terminator.
    Any,
    /// `[...]`, or a shorthand such as `\d` on its own, which is compiled as a one-member class.
    Class {
        negated: bool,
        members: Vec<ClassMember>,
    },
    /// `^`, the start of the input only (no `MULTILINE`).
    Start,
    /// `$`, the end of the input or just before a final line terminator.
    End,
}

/// The three quantifier forms: `*`, `*?` and `*+` (and likewise for every quantifier).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RepeatKind {
    /// As many as possible, giving back one at a time.
    Greedy,
    /// As few as possible, taking one more at a time.
    Reluctant,
    /// `*+`, which takes as much as it can and never gives any back.
    Possessive,
}

/// One item inside a class.
#[derive(Debug, Clone, Copy)]
enum ClassMember {
    Single(char),
    Range(char, char),
    /// One of the `\d \w \s` shorthands, with its negation folded in.
    Shorthand(Shorthand, bool),
}

/// The three predefined classes. The upper-case forms are the same with negation set.
#[derive(Debug, Clone, Copy)]
enum Shorthand {
    Digit,
    Word,
    Space,
}

/// Java's line terminators, which `.` refuses and `$` allows itself to stand before: `\n`, `\r`,
/// U+0085, U+2028 and U+2029 (and `\r\n` as a pair, which `match_node` handles for `$`). The
/// golden's `a.c` on `a\rc` row is `false`, which is why `\r` is here and not only `\n`.
fn is_line_terminator(character: char) -> bool {
    matches!(character, '\n' | '\r' | '\u{85}' | '\u{2028}' | '\u{2029}')
}

impl Shorthand {
    /// Whether `character` is in the class, before any negation.
    fn contains(&self, character: char) -> bool {
        match self {
            // Java's \d is ASCII digits only without UNICODE_CHARACTER_CLASS, and so are \w and \s.
            // The golden pins two: `\d` on U+0663 and `\w` on U+00E9 are both `false`.
            Shorthand::Digit => character.is_ascii_digit(),
            Shorthand::Word => character.is_ascii_alphanumeric() || character == '_',
            Shorthand::Space => matches!(character, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r'),
        }
    }
}

/// A recursive-descent parser over the pattern's code points. `'a` says it borrows them.
struct Parser<'a> {
    pattern: &'a [char],
    /// The next character to read; also what most refusal indices are measured from.
    at: usize,
    /// The pattern as given, for the refusal's third line.
    source: String,
}

impl<'a> Parser<'a> {
    /// A refusal with the given description and index.
    fn error(&self, description: &str, index: i32) -> PatternSyntaxError {
        PatternSyntaxError {
            description: description.to_string(),
            index,
            pattern: self.source.clone(),
        }
    }

    /// The next character, or `None` at the end, without consuming it.
    fn peek(&self) -> Option<char> {
        self.pattern.get(self.at).copied()
    }

    /// A whole alternation, which is what a group and the pattern itself both are.
    fn parse_alternation(&mut self) -> Result<Node, PatternSyntaxError> {
        let mut branches = vec![self.parse_sequence()?];
        while self.peek() == Some('|') {
            self.at += 1;
            branches.push(self.parse_sequence()?);
        }
        Ok(if branches.len() == 1 {
            branches.pop().expect("one branch")
        } else {
            Node::Alternation(branches)
        })
    }

    /// Atoms and their quantifiers up to a `|`, a `)` or the end.
    fn parse_sequence(&mut self) -> Result<Node, PatternSyntaxError> {
        let mut nodes: Vec<Node> = Vec::new();
        while let Some(character) = self.peek() {
            if character == '|' || character == ')' {
                break;
            }
            // Measured (golden `{2}` and `|{2}` rows, `find` true on every input): `{2}` with
            // nothing before it compiles, and matches everywhere. Java reads
            // the brace as a quantifier over the EMPTY atom rather than as a literal or an error,
            // which is why it is the one metacharacter that may open an expression. `*` in the
            // same place is a dangling meta character.
            let atom = if character == '{' {
                Node::Sequence(Vec::new())
            } else {
                self.parse_atom()?
            };
            let node = self.parse_quantifier(atom)?;
            nodes.push(node);
        }
        Ok(Node::Sequence(nodes))
    }

    /// One atom: an anchor, `.`, a group, a class, an escape or a literal. `start` is kept because
    /// several refusals are measured from the atom's first character.
    fn parse_atom(&mut self) -> Result<Node, PatternSyntaxError> {
        let start = self.at;
        let character = self.peek().expect("a character");
        self.at += 1;
        match character {
            '^' => Ok(Node::Start),
            '$' => Ok(Node::End),
            '.' => Ok(Node::Any),
            '(' => {
                if self.peek() == Some('?') {
                    // Non-capturing groups, look-around and named groups all start `(?`; only the
                    // first is representable here and telling them apart is not worth guessing.
                    // Java reads `(?` as an inline modifier and names it that way; the
                    // constructs it introduces are all outside this subset either way.
                    // Golden: `(?` is refused at index 2.
                    return Err(self.error("Unknown inline modifier", self.at as i32 + 1));
                }
                // A group only groups: it is not captured, because nothing here reads a capture.
                let inner = self.parse_alternation()?;
                if self.peek() != Some(')') {
                    // Measured: the index is the character AFTER the `(` (golden `(` is 1,
                    // `abc(` is 4).
                    return Err(self.error("Unclosed group", start as i32 + 1));
                }
                self.at += 1;
                Ok(inner)
            }
            // The character before the `)`, so a leading one is -1 (golden `)`, `a)`, `abc)`).
            ')' => Err(self.error("Unmatched closing ')'", start as i32 - 1)),
            '[' => self.parse_class(start),
            '\\' => self.parse_escape(start),
            // A quantifier with nothing to quantify (golden `*a` at 0, `a**` at 2).
            '*' | '+' | '?' => Err(self.error(
                &format!("Dangling meta character '{character}'"),
                start as i32,
            )),
            other => Ok(Node::Literal(other)),
        }
    }

    /// A `\` and what follows it, either a shorthand class or a literal.
    fn parse_escape(&mut self, start: usize) -> Result<Node, PatternSyntaxError> {
        let Some(character) = self.peek() else {
            // Measured: a trailing backslash reports the length, and gets no caret line.
            return Err(self.error("Unexpected internal error", self.at as i32));
        };
        self.at += 1;
        match character {
            'd' => Ok(shorthand_class(Shorthand::Digit, false)),
            'D' => Ok(shorthand_class(Shorthand::Digit, true)),
            'w' => Ok(shorthand_class(Shorthand::Word, false)),
            'W' => Ok(shorthand_class(Shorthand::Word, true)),
            's' => Ok(shorthand_class(Shorthand::Space, false)),
            'S' => Ok(shorthand_class(Shorthand::Space, true)),
            'n' => Ok(Node::Literal('\n')),
            'r' => Ok(Node::Literal('\r')),
            't' => Ok(Node::Literal('\t')),
            'f' => Ok(Node::Literal('\u{c}')),
            'a' => Ok(Node::Literal('\u{7}')),
            'e' => Ok(Node::Literal('\u{1b}')),
            // Back references and octal escapes, outside the subset.
            '0'..='9' => Err(self.error("Unsupported escape sequence", start as i32)),
            // Every other letter escape (`\b`, `\p`, `\x`, `\Q`, ...), which in Java either means
            // something this does not support or is itself an error. `other if ..` is a match
            // arm with a condition.
            other if other.is_ascii_alphabetic() => {
                Err(self.error("Unsupported escape sequence", start as i32))
            }
            // An escaped non-letter is that character: `\.`, `\[`, `\\`.
            other => Ok(Node::Literal(other)),
        }
    }

    /// `[...]` after the `[`, with `start` its index for the unclosed-class refusal.
    fn parse_class(&mut self, start: usize) -> Result<Node, PatternSyntaxError> {
        let mut negated = false;
        if self.peek() == Some('^') {
            negated = true;
            self.at += 1;
        }
        let mut members: Vec<ClassMember> = Vec::new();
        // A `]` first in the class is a literal, which is Java's rule as well as POSIX's.
        let mut first = true;
        loop {
            // Measured: an unclosed class points at its `[` (golden `[` 0, `a[` 1, `abc[` 3).
            let Some(character) = self.peek() else {
                return Err(self.error("Unclosed character class", start as i32));
            };
            if character == ']' && !first {
                self.at += 1;
                return Ok(Node::Class { negated, members });
            }
            first = false;
            self.at += 1;
            let low = if character == '\\' {
                match self.parse_escape(self.at - 1)? {
                    Node::Literal(literal) => Some(literal),
                    Node::Class {
                        negated: inner_negated,
                        members: inner,
                    } => {
                        // A shorthand inside a class contributes its whole membership.
                        for member in inner {
                            members.push(match member {
                                ClassMember::Shorthand(kind, flip) => {
                                    ClassMember::Shorthand(kind, flip != inner_negated)
                                }
                                other => other,
                            });
                        }
                        None
                    }
                    _ => None,
                }
            } else if character == '[' {
                // Java reads a nested `[` as the start of a union, which this does not represent.
                return Err(self.error("Unsupported class union", self.at as i32 - 1));
            } else {
                Some(character)
            };
            let Some(low) = low else {
                continue;
            };
            // A `-` is a range only with something other than `]` after it, so `[a-]` and `[-a]`
            // hold a literal `-` (both golden rows compile).
            if self.peek() == Some('-')
                && self
                    .pattern
                    .get(self.at + 1)
                    .copied()
                    .is_some_and(|c| c != ']')
            {
                self.at += 1;
                let high = self.peek().expect("a range end");
                self.at += 1;
                if high == '\\' {
                    // Measured: `[a-\d]` is "Illegal character range" pointing at the letter, one
                    // past the backslash. An escape that stands for a single character would be a
                    // valid range end in the reference; this subset has none of those, so the
                    // whole case lands here.
                    return Err(self.error("Illegal character range", self.at as i32));
                }
                // Golden `[z-a]` points at the `a`, index 3.
                if high < low {
                    return Err(self.error("Illegal character range", self.at as i32 - 1));
                }
                members.push(ClassMember::Range(low, high));
            } else {
                members.push(ClassMember::Single(low));
            }
        }
    }

    /// The quantifier after `atom`, if any, and the form (`?` reluctant, `+` possessive) after it.
    fn parse_quantifier(&mut self, atom: Node) -> Result<Node, PatternSyntaxError> {
        let Some(character) = self.peek() else {
            return Ok(atom);
        };
        let (min, max) = match character {
            '*' => (0, usize::MAX),
            '+' => (1, usize::MAX),
            '?' => (0, 1),
            '{' => return self.parse_braced_quantifier(atom),
            _ => return Ok(atom),
        };
        self.at += 1;
        Ok(Node::Repeat {
            node: Box::new(atom),
            min,
            max,
            kind: self.parse_repeat_kind(),
        })
    }

    /// `{n}`, `{n,}` or `{n,m}`, with `self.at` on the `{`.
    fn parse_braced_quantifier(&mut self, atom: Node) -> Result<Node, PatternSyntaxError> {
        let open = self.at;
        self.at += 1;
        let mut low = String::new();
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            low.push(self.peek().expect("a digit"));
            self.at += 1;
        }
        if low.is_empty() {
            // Measured: `{a}`, `{}` and `a{,2}` are all "Illegal repetition" at the character
            // after the brace, while `{2` unterminated is an "Unclosed counted closure".
            return Err(self.error("Illegal repetition", open as i32 + 1));
        }
        // A count too large for a `usize` is refused rather than wrapped; no golden row has one,
        // so the description and index here are this port's choice.
        let min: usize = low
            .parse()
            .map_err(|_| self.error("Illegal repetition", open as i32))?;
        let max = match self.peek() {
            Some('}') => min,
            Some(',') => {
                self.at += 1;
                let mut high = String::new();
                while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    high.push(self.peek().expect("a digit"));
                    self.at += 1;
                }
                if high.is_empty() {
                    usize::MAX
                } else {
                    high.parse()
                        .map_err(|_| self.error("Illegal repetition", open as i32))?
                }
            }
            _ => return Err(self.error("Unclosed counted closure", self.at as i32)),
        };
        if self.peek() != Some('}') {
            return Err(self.error("Unclosed counted closure", self.at as i32));
        }
        // Measured: a range whose low is above its high reports the CLOSING BRACE (golden
        // `a{2,1}` at 5, `ab{3,2}` at 6).
        if min > max {
            return Err(self.error("Illegal repetition range", self.at as i32));
        }
        self.at += 1;
        Ok(Node::Repeat {
            node: Box::new(atom),
            min,
            max,
            kind: self.parse_repeat_kind(),
        })
    }

    /// The suffix after a quantifier: `?` reluctant, `+` possessive, otherwise greedy.
    fn parse_repeat_kind(&mut self) -> RepeatKind {
        match self.peek() {
            Some('?') => {
                self.at += 1;
                RepeatKind::Reluctant
            }
            Some('+') => {
                self.at += 1;
                RepeatKind::Possessive
            }
            _ => RepeatKind::Greedy,
        }
    }
}

/// A shorthand as a node. The negation goes on the member, not the class, so the same member can
/// be copied into an enclosing `[...]` with its sense intact (`[^\W\d]` is a golden row).
fn shorthand_class(kind: Shorthand, negated: bool) -> Node {
    Node::Class {
        negated: false,
        members: vec![ClassMember::Shorthand(kind, negated)],
    }
}

impl Pattern {
    /// `Pattern.compile`, refusing what this subset cannot represent rather than mismatching it.
    pub fn compile(pattern: &str) -> Result<Pattern, PatternSyntaxError> {
        let characters: Vec<char> = pattern.chars().collect();
        let mut parser = Parser {
            pattern: &characters,
            at: 0,
            source: pattern.to_string(),
        };
        // `pattern.chars().collect()` splits the text into code points once, so the parser can
        // index it like an array.
        let node = parser.parse_alternation()?;
        if parser.at < characters.len() {
            // Only a `)` can stop the parse without consuming the input, and Java points at the
            // character before it.
            return Err(parser.error("Unmatched closing ')'", parser.at as i32 - 1));
        }
        Ok(Pattern {
            node,
            source: pattern.to_string(),
        })
    }

    /// `Pattern.pattern()`: the text compiled.
    pub fn source(&self) -> &str {
        &self.source
    }

    /// `Matcher.find()`: whether the pattern occurs anywhere, which is what GATK asks.
    ///
    /// Every start from `0` to the length inclusive is tried, leftmost first, and the first that
    /// matches ends the search; the end position is the one that lets an empty pattern match an
    /// empty input. Only the yes or no is returned, which is all `filterCollectionByExpressions`
    /// reads.
    pub fn find(&self, text: &str) -> bool {
        let characters: Vec<char> = text.chars().collect();
        (0..=characters.len()).any(|start| match_node(&self.node, &characters, start, &|_| true))
    }

    /// `Matcher.matches()`: whether it occupies the whole input.
    pub fn matches(&self, text: &str) -> bool {
        let characters: Vec<char> = text.chars().collect();
        match_node(&self.node, &characters, 0, &|at| at == characters.len())
    }
}

/// One node against `text` from `at`, handing every possible end to `continuation`.
///
/// This is the backtracking Java does: the continuation is what comes after, so a greedy repeat
/// that consumed too much can give a repetition back and ask again.
///
/// Each arm consumes what its node matches and calls `continuation` with the position after it;
/// `&&` and `.any(..)` stop at the first success, which is the "first alternative first" order.
/// `text.get(at)` is `None` past the end, so no arm reads out of bounds.
fn match_node(node: &Node, text: &[char], at: usize, continuation: &dyn Fn(usize) -> bool) -> bool {
    match node {
        Node::Literal(literal) => text.get(at) == Some(literal) && continuation(at + 1),
        Node::Any => text
            .get(at)
            .is_some_and(|c| !is_line_terminator(*c) && continuation(at + 1)),
        Node::Class { negated, members } => text.get(at).is_some_and(|c| {
            let inside = members.iter().any(|member| match member {
                ClassMember::Single(single) => single == c,
                ClassMember::Range(low, high) => low <= c && c <= high,
                ClassMember::Shorthand(kind, flip) => kind.contains(*c) != *flip,
            });
            inside != *negated && continuation(at + 1)
        }),
        Node::Start => at == 0 && continuation(at),
        Node::End => {
            // Java's `$` without MULTILINE: the end, or just before a line terminator that ends it.
            let end = at == text.len()
                || (at + 1 == text.len() && is_line_terminator(text[at]))
                || (at + 2 == text.len() && text[at] == '\r' && text[at + 1] == '\n');
            end && continuation(at)
        }
        Node::Sequence(nodes) => match_sequence(nodes, text, at, continuation),
        Node::Alternation(branches) => branches
            .iter()
            .any(|branch| match_node(branch, text, at, continuation)),
        Node::Repeat {
            node,
            min,
            max,
            kind,
        } => match_repeat(node, *min, *max, *kind, text, at, continuation),
    }
}

/// A sequence: the first node, with "the rest of the sequence, then `continuation`" as its
/// continuation. `split_first` is the head and the tail of a slice, or `None` when it is empty.
fn match_sequence(
    nodes: &[Node],
    text: &[char],
    at: usize,
    continuation: &dyn Fn(usize) -> bool,
) -> bool {
    match nodes.split_first() {
        None => continuation(at),
        Some((head, rest)) => match_node(head, text, at, &|next| {
            match_sequence(rest, text, next, continuation)
        }),
    }
}

/// A repeat with `min` still owed and up to `max` allowed, from `at`.
///
/// The mandatory repetitions are matched first, one at a time. After that a greedy repeat tries
/// one more before trying to stop, a reluctant one tries to stop first, and both offer every
/// alternative on failure. A possessive repeat is different: it takes the repeated node's first
/// match at each step and never revisits them, which is why `^.*+1$` fails on `s1` while `^.*1$`
/// matches (both golden rows).
fn match_repeat(
    node: &Node,
    min: usize,
    max: usize,
    kind: RepeatKind,
    text: &[char],
    at: usize,
    continuation: &dyn Fn(usize) -> bool,
) -> bool {
    if kind == RepeatKind::Possessive {
        // Take as many as possible and never give one back, which is what `*+` means.
        let mut end = at;
        let mut taken = 0;
        while taken < max {
            // `Cell` lets the closure record where the node ended although the closure only has
            // read access to its surroundings; answering `true` stops at the first match.
            let next = std::cell::Cell::new(None);
            match_node(node, text, end, &|reached| {
                next.set(Some(reached));
                true
            });
            match next.get() {
                // A repetition that matched nothing ends the loop, as below.
                Some(reached) if reached != end => {
                    end = reached;
                    taken += 1;
                }
                _ => break,
            }
        }
        return taken >= min && continuation(end);
    }

    if min > 0 {
        return match_node(node, text, at, &|next| {
            match_repeat(
                node,
                min - 1,
                max.saturating_sub(1),
                kind,
                text,
                next,
                continuation,
            )
        });
    }
    if max == 0 {
        return continuation(at);
    }
    // A repetition that consumed nothing would loop for ever, so it stops the loop instead: the
    // guard is the reference's own. Golden `^(a*)*$` rows: true on `aaa` and on the empty input.
    let another = |at: usize| {
        match_node(node, text, at, &|next| {
            next != at && match_repeat(node, 0, max - 1, kind, text, next, continuation)
        })
    };
    match kind {
        RepeatKind::Reluctant => continuation(at) || another(at),
        _ => another(at) || continuation(at),
    }
}

/// `Utils.filterCollectionByExpressions`, which is what `-se` and `-xl-se` are.
///
/// The result keeps the **source** order, not the expression order, since it accumulates into a
/// `LinkedHashSet` while walking the values. Every expression is compiled up front, so an
/// uncompilable one refuses before any value has been looked at, and matching is `find()` unless
/// `exact_match` is asked for, which turns the whole thing into equality.
///
/// The reference's `LinkedHashSet` also drops a value that occurs twice in `values`; this `Vec`
/// keeps both. The callers pass a VCF's sample names.
///
/// `.collect::<Result<Vec<_>, _>>()?` compiles every expression and stops at the first refusal,
/// as `compilePatterns`' loop does when `Pattern.compile` throws.
pub fn filter_collection_by_expressions(
    values: &[String],
    expressions: &[String],
    exact_match: bool,
) -> Result<Vec<String>, PatternSyntaxError> {
    let patterns = if exact_match {
        Vec::new()
    } else {
        expressions
            .iter()
            .map(|expression| Pattern::compile(expression))
            .collect::<Result<Vec<_>, _>>()?
    };
    let mut filtered: Vec<String> = Vec::new();
    for value in values {
        // Equality first, whether or not an exact match was asked for: an expression that is also a
        // value selects that value before any pattern is consulted.
        let selected = expressions.contains(value)
            || (!exact_match && patterns.iter().any(|pattern| pattern.find(value)));
        if selected {
            filtered.push(value.clone());
        }
    }
    Ok(filtered)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn find(pattern: &str, text: &str) -> bool {
        Pattern::compile(pattern).expect("compiled").find(text)
    }

    #[test]
    fn find_is_a_search_and_matches_is_not() {
        assert!(find("s1", "xs10"));
        assert!(!Pattern::compile("s1").expect("compiled").matches("xs10"));
        assert!(Pattern::compile("s1").expect("compiled").matches("s1"));
    }

    #[test]
    fn the_anchors_are_the_whole_input_without_multiline() {
        assert!(find("^s1$", "s1"));
        assert!(!find("^s1$", "xs10"));
        assert!(find("^NA", "NA12878"));
        // `$` stands before a final line terminator, which is Java's rule.
        assert!(find("^s1$", "s1\n"));
        assert!(find("^s1$", "s1\r\n"));
        assert!(!find("^s1$", "s1\nx"));
    }

    #[test]
    fn a_dot_is_every_character_but_a_line_terminator() {
        assert!(find(".", "a"));
        assert!(!find(".", "\n"));
        assert!(find("a.c", "abc"));
    }

    #[test]
    fn classes_carry_ranges_negation_and_shorthands() {
        assert!(find("[a-c]", "b"));
        assert!(!find("[a-c]", "d"));
        assert!(find("[^a-c]", "d"));
        assert!(find("[]]", "]"));
        assert!(find("\\d", "s1"));
        assert!(!find("\\D", "1"));
        assert!(find("[\\d_]", "_"));
    }

    #[test]
    fn quantifiers_backtrack_the_way_java_does() {
        assert!(find("^a*b$", "aaab"));
        assert!(find("^a+b$", "ab"));
        assert!(!find("^a+b$", "b"));
        assert!(find("^ab?$", "a"));
        assert!(find("^a{2,3}$", "aaa"));
        assert!(!find("^a{2,3}$", "aaaa"));
        assert!(find("^a{2}$", "aa"));
        // Greedy first, then giving one back: without backtracking this would fail.
        assert!(find("^.*1$", "NA12891"));
        // Possessive takes everything and keeps it, so the same pattern fails.
        assert!(!find("^.*+1$", "NA12891"));
        assert!(find("^.*?1$", "NA12891"));
    }

    #[test]
    fn alternation_and_groups_are_tried_in_order() {
        assert!(find("^(NA|s)1$", "s1"));
        assert!(find("^(NA|s)12878$", "NA12878"));
        assert!(!find("^(NA|s)1$", "x1"));
        assert!(find("(ab)+c", "ababc"));
    }

    #[test]
    fn a_repetition_that_consumes_nothing_stops_rather_than_spinning() {
        assert!(find("^(a*)*$", "aaa"));
        assert!(find("^(a*)*$", ""));
    }

    #[test]
    fn an_unclosed_class_is_the_reference_s_message() {
        let error = Pattern::compile("[").expect_err("refused");
        assert_eq!(error.description, "Unclosed character class");
        assert_eq!(error.index, 0);
        assert_eq!(
            error.message(),
            "Unclosed character class near index 0\n[\n^"
        );
    }

    #[test]
    fn what_this_subset_cannot_represent_is_refused_rather_than_mismatched() {
        assert!(Pattern::compile("(?:a)").is_err());
        assert!(Pattern::compile("\\b").is_err());
        assert!(Pattern::compile("(a").is_err());
        assert!(Pattern::compile("a)").is_err());
        assert!(Pattern::compile("*a").is_err());
    }

    #[test]
    fn the_filter_keeps_the_source_order_and_compiles_before_it_looks() {
        let values: Vec<String> = ["tumor", "s1", "NA12891", "xs10", "s0", "NA12878"]
            .iter()
            .map(|value| value.to_string())
            .collect();
        let matched = filter_collection_by_expressions(&values, &["^NA".to_string()], false)
            .expect("matched");
        assert_eq!(matched, vec!["NA12891".to_string(), "NA12878".to_string()]);

        // An exact match is equality, and equality is checked before any pattern is.
        let exact =
            filter_collection_by_expressions(&values, &["^NA".to_string()], true).expect("matched");
        assert!(exact.is_empty());

        assert!(filter_collection_by_expressions(&values, &["[".to_string()], false).is_err());
    }
}
