//! A `java.util.regex` subset, enough for the patterns GATK compiles from command lines.
//!
//! GATK hands user strings to `Pattern.compile` in several places and then calls **`find()`**, not
//! `matches()`: `Utils.filterCollectionByExpressions`, which is what `-se` and `-xl-se` go through,
//! and `ClipReads`' sequence clipper. An unanchored search is what those arguments mean, and it is
//! why `-se s1` selects a sample named `xs10`. Reproducing the selection means reproducing the
//! search, so the matcher is here rather than assumed away.
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
//!    text reaches the user as the regex engine wrote it. Only the descriptions listed in
//!    [`PatternSyntaxError`] are produced; a golden pins the unclosed-class one.

/// A compiled pattern.
#[derive(Debug, Clone)]
pub struct Pattern {
    node: Node,
    source: String,
}

/// What `Pattern.compile` throws, with the three parts Java's `getMessage` joins with newlines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatternSyntaxError {
    /// The description, e.g. `Unclosed character class`.
    pub description: String,
    /// The index the caret points at, in characters.
    pub index: usize,
    pub pattern: String,
}

impl PatternSyntaxError {
    /// `getMessage()`: the description, the pattern, and a caret line.
    pub fn message(&self) -> String {
        format!(
            "{} near index {}\n{}\n{}^",
            self.description,
            self.index,
            self.pattern,
            " ".repeat(self.index)
        )
    }

    pub fn java_class(&self) -> &'static str {
        "java.util.regex.PatternSyntaxException"
    }
}

#[derive(Debug, Clone)]
enum Node {
    /// The branches of a `|`, tried in the order they were written.
    Alternation(Vec<Node>),
    Sequence(Vec<Node>),
    Repeat {
        node: Box<Node>,
        min: usize,
        max: usize,
        kind: RepeatKind,
    },
    Literal(char),
    /// `.`, which is every character but a line terminator.
    Any,
    Class {
        negated: bool,
        members: Vec<ClassMember>,
    },
    Start,
    End,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RepeatKind {
    Greedy,
    Reluctant,
    /// `*+`, which takes as much as it can and never gives any back.
    Possessive,
}

#[derive(Debug, Clone, Copy)]
enum ClassMember {
    Single(char),
    Range(char, char),
    /// One of the `\d \w \s` shorthands, with its negation folded in.
    Shorthand(Shorthand, bool),
}

#[derive(Debug, Clone, Copy)]
enum Shorthand {
    Digit,
    Word,
    Space,
}

/// Java's line terminators, which `.` refuses and `$` allows itself to stand before.
fn is_line_terminator(character: char) -> bool {
    matches!(character, '\n' | '\r' | '\u{85}' | '\u{2028}' | '\u{2029}')
}

impl Shorthand {
    fn contains(&self, character: char) -> bool {
        match self {
            // Java's \d is ASCII digits only without UNICODE_CHARACTER_CLASS, and so are \w and \s.
            Shorthand::Digit => character.is_ascii_digit(),
            Shorthand::Word => character.is_ascii_alphanumeric() || character == '_',
            Shorthand::Space => matches!(character, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r'),
        }
    }
}

struct Parser<'a> {
    pattern: &'a [char],
    at: usize,
    source: String,
}

impl<'a> Parser<'a> {
    fn error(&self, description: &str, index: usize) -> PatternSyntaxError {
        PatternSyntaxError {
            description: description.to_string(),
            index,
            pattern: self.source.clone(),
        }
    }

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

    fn parse_sequence(&mut self) -> Result<Node, PatternSyntaxError> {
        let mut nodes: Vec<Node> = Vec::new();
        while let Some(character) = self.peek() {
            if character == '|' || character == ')' {
                break;
            }
            let atom = self.parse_atom()?;
            let node = self.parse_quantifier(atom)?;
            nodes.push(node);
        }
        Ok(Node::Sequence(nodes))
    }

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
                    return Err(self.error("Unsupported group", self.at));
                }
                let inner = self.parse_alternation()?;
                if self.peek() != Some(')') {
                    return Err(self.error("Unclosed group", self.pattern.len()));
                }
                self.at += 1;
                Ok(inner)
            }
            ')' => Err(self.error("Unmatched closing ')'", start)),
            '[' => self.parse_class(start),
            '\\' => self.parse_escape(start),
            '*' | '+' | '?' => {
                Err(self.error(&format!("Dangling meta character '{character}'"), start))
            }
            other => Ok(Node::Literal(other)),
        }
    }

    /// A `\` and what follows it, either a shorthand class or a literal.
    fn parse_escape(&mut self, start: usize) -> Result<Node, PatternSyntaxError> {
        let Some(character) = self.peek() else {
            return Err(self.error("Unexpected internal error", start));
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
            '0'..='9' => Err(self.error("Unsupported escape sequence", start)),
            other if other.is_ascii_alphabetic() => {
                Err(self.error("Unsupported escape sequence", start))
            }
            other => Ok(Node::Literal(other)),
        }
    }

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
            let Some(character) = self.peek() else {
                return Err(self.error("Unclosed character class", start));
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
                return Err(self.error("Unsupported class union", self.at - 1));
            } else {
                Some(character)
            };
            let Some(low) = low else {
                continue;
            };
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
                    return Err(self.error("Unsupported escaped range", self.at - 1));
                }
                if high < low {
                    return Err(self.error("Illegal character range", self.at - 1));
                }
                members.push(ClassMember::Range(low, high));
            } else {
                members.push(ClassMember::Single(low));
            }
        }
    }

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

    fn parse_braced_quantifier(&mut self, atom: Node) -> Result<Node, PatternSyntaxError> {
        let open = self.at;
        self.at += 1;
        let mut low = String::new();
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            low.push(self.peek().expect("a digit"));
            self.at += 1;
        }
        if low.is_empty() {
            return Err(self.error("Unclosed counted closure", open + 1));
        }
        let min: usize = low
            .parse()
            .map_err(|_| self.error("Illegal repetition", open))?;
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
                        .map_err(|_| self.error("Illegal repetition", open))?
                }
            }
            _ => return Err(self.error("Unclosed counted closure", self.at)),
        };
        if self.peek() != Some('}') {
            return Err(self.error("Unclosed counted closure", self.at));
        }
        self.at += 1;
        Ok(Node::Repeat {
            node: Box::new(atom),
            min,
            max,
            kind: self.parse_repeat_kind(),
        })
    }

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
        let node = parser.parse_alternation()?;
        if parser.at < characters.len() {
            // Only a `)` can stop the parse without consuming the input.
            return Err(parser.error("Unmatched closing ')'", parser.at));
        }
        Ok(Pattern {
            node,
            source: pattern.to_string(),
        })
    }

    pub fn source(&self) -> &str {
        &self.source
    }

    /// `Matcher.find()`: whether the pattern occurs anywhere, which is what GATK asks.
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
            let next = std::cell::Cell::new(None);
            match_node(node, text, end, &|reached| {
                next.set(Some(reached));
                true
            });
            match next.get() {
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
    // guard is the reference's own.
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
