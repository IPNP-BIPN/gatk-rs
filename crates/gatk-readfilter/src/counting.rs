//! Ported from `org.broadinstitute.hellbender.engine.filters.CountingReadFilter` (GATK 4.7.0.0;
//! `CountingReadFilter.java` is unchanged since 4.6.2.0).
//!
//! This is the first piece of *output text* in this crate. Every GATK tool that reads reads ends
//! its run by printing this block, so the summary is part of the bytes a run produces, and the
//! byte-identity claim covers it.
//!
//! It is also what makes the order of a conjunction observable. `fromList` wraps every filter of
//! the tool's list in its own `CountingReadFilter` and joins them with the counting `and`, which
//! short-circuits, so only the first failing filter of the list increments. That holds between
//! the filters of the list and not inside one: `WellformedReadFilter` is eight filters joined by
//! the plain `ReadFilter.and`, wrapped as a single leaf, so its count cannot say which of the
//! eight rejected a read, and [`parse`] builds it as one leaf.
//!
//! # Where the bytes hide
//!
//! Three formatting details that no reasonable person would guess and that a diff catches at once:
//!
//! - a leaf's line ends `"<name> \n"`, with a space before the newline;
//! - a composite's line ends `"<name>\n"`, without it, **except** when its count is zero, in which
//!   case it takes the leaf's branch and the space comes back;
//! - the flattened summary ends without a trailing newline at all.
//!
//! The three as the Java writes them, in that order (the leaf's `getSummaryLineForLevel`, the
//! non-zero branch of `CountingBinopReadFilter`'s, and the end of the flattened summary):
//!
//! ```java
//! " read(s) filtered by: " + getName() + " \n";
//! " read(s) filtered by: " + getName() + "\n"
//! " total reads filtered out of "+ this.totalCount+ " reads processed");
//! ```
//!
//! # Reading this file without reading Rust
//!
//! - `enum Node` is a closed set of shapes, each carrying its own fields; `match` on it must
//!   name every shape, which the compiler checks. Java spreads the same four shapes over four
//!   subclasses;
//! - `Box<T>` is a heap-allocated `T` the node owns; a tree of nodes needs it, since a node
//!   cannot contain itself directly. `Box<dyn Fn(&BamRecord) -> bool>` is an owned closure,
//!   Java's `Predicate<GATKRead>`;
//! - `&mut self` lets `test` update the counts; only one such borrow can exist at a time;
//! - `Option<T>` is a value or nothing; `?` after one returns nothing from the function.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

/// A node of a `CountingReadFilter` tree.
enum Node {
    /// `CountingReadFilter` wrapping a plain filter.
    ///
    /// The predicate is boxed rather than a function pointer because the header-dependent filters
    /// have to carry the header with them, and the engine wraps those exactly like any other.
    Leaf {
        name: String,
        filter: Box<dyn Fn(&BamRecord) -> bool>,
    },
    /// `CountingAndReadFilter`, which also counts the reads it saw.
    And(Box<Counting>, Box<Counting>),
    /// `CountingOrReadFilter`.
    Or(Box<Counting>, Box<Counting>),
    /// `CountingNegateReadFilter`.
    Not(Box<Counting>),
}

/// A `CountingReadFilter`: a filter that remembers how many reads it rejected.
pub struct Counting {
    node: Node,
    filtered_count: u64,
    /// Only `CountingAndReadFilter` keeps this, and only its flattened summary reads it.
    total_count: u64,
}

impl Counting {
    /// `new CountingReadFilter(readFilter)`: a plain filter, counted under its own name.
    ///
    /// `impl Fn(&BamRecord) -> bool + 'static` accepts any closure or function of that shape
    /// that owns everything it captures, so it can be stored for the whole run.
    pub fn leaf(name: &str, filter: impl Fn(&BamRecord) -> bool + 'static) -> Counting {
        Counting {
            node: Node::Leaf {
                name: name.to_string(),
                filter: Box::new(filter),
            },
            filtered_count: 0,
            total_count: 0,
        }
    }

    /// A composite with both counters at zero, as each counting subclass starts.
    fn binop(node: Node) -> Counting {
        Counting {
            node,
            filtered_count: 0,
            total_count: 0,
        }
    }

    /// `CountingReadFilter.and`: a `CountingAndReadFilter` over this and `other`. `self` is
    /// taken by value, so the two children move into the new node.
    pub fn and(self, other: Counting) -> Counting {
        Counting::binop(Node::And(Box::new(self), Box::new(other)))
    }

    /// `CountingReadFilter.or`: a `CountingOrReadFilter`.
    pub fn or(self, other: Counting) -> Counting {
        Counting::binop(Node::Or(Box::new(self), Box::new(other)))
    }

    /// `CountingReadFilter.negate`: a `CountingNegateReadFilter`.
    pub fn negate(self) -> Counting {
        Counting::binop(Node::Not(Box::new(self)))
    }

    /// `CountingReadFilter.fromList`: a left-nested chain of ANDs, the shape the engine builds
    /// from `--read-filter` arguments. An empty list is `ALLOW_ALL_READS`, not an error.
    ///
    /// ```java
    /// CountingReadFilter compositeFilter = new CountingReadFilter(readFilters.get(0));
    /// for (int i = 1; i < readFilters.size(); i++) {
    ///     compositeFilter = compositeFilter.and(new CountingReadFilter(readFilters.get(i)));
    /// }
    /// ```
    pub fn from_list(filters: Vec<(&str, crate::ReadFilter)>) -> Counting {
        let mut iter = filters.into_iter();
        let Some((name, filter)) = iter.next() else {
            return Counting::leaf("AllowAllReadsReadFilter", crate::allow_all_reads);
        };
        let mut composite = Counting::leaf(name, filter);
        for (name, filter) in iter {
            composite = composite.and(Counting::leaf(name, filter));
        }
        composite
    }

    /// `getFilteredCount`.
    pub fn filtered_count(&self) -> u64 {
        self.filtered_count
    }

    /// `getName`: a leaf's own name, `(a AND b)`, `(a OR b)` or `NOT a`, built recursively.
    pub fn name(&self) -> String {
        match &self.node {
            Node::Leaf { name, .. } => name.clone(),
            Node::And(lhs, rhs) => format!("({} AND {})", lhs.name(), rhs.name()),
            Node::Or(lhs, rhs) => format!("({} OR {})", lhs.name(), rhs.name()),
            Node::Not(delegate) => format!("NOT {}", delegate.name()),
        }
    }

    /// `test`: the node's verdict, counting a rejection at this node.
    ///
    /// Every node counts its own rejections, so a read that fails a leaf inside an AND increments
    /// both the leaf and each AND above it that it fails.
    pub fn test(&mut self, read: &BamRecord) -> bool {
        let accept = match &mut self.node {
            Node::Leaf { filter, .. } => filter(read),
            // `lhs.test(read) && rhs.test(read)`: short-circuiting, so a read rejected by the left
            // never reaches the right and never increments its counter.
            Node::And(lhs, rhs) => lhs.test(read) && rhs.test(read),
            Node::Or(lhs, rhs) => lhs.test(read) || rhs.test(read),
            Node::Not(delegate) => !delegate.test(read),
        };
        if !accept {
            self.filtered_count += 1;
        }
        // Only the AND node counts what it saw, and it counts it whether or not the read passed.
        if matches!(self.node, Node::And(_, _)) {
            self.total_count += 1;
        }
        accept
    }

    /// `CountingReadFilter.getSummaryLine`.
    pub fn summary_line(&self) -> String {
        self.summary_line_for_level(0)
    }

    /// `getSummaryLineForLevel`, dispatched by shape where the Java dispatches by subclass.
    fn summary_line_for_level(&self, indent_level: usize) -> String {
        let indent = "  ".repeat(indent_level);
        match &self.node {
            // A leaf, and a negation, take the base implementation: name then a space then the
            // newline.
            Node::Leaf { .. } | Node::Not(_) => {
                format!(
                    "{indent}{} read(s) filtered by: {} \n",
                    self.filtered_count,
                    self.name()
                )
            }
            Node::And(lhs, rhs) => {
                // At the top level an all-AND tree prints flattened, with no line for the
                // composite itself. Anything else in the tree falls back to the nested form.
                if indent_level == 0 {
                    if let Some(simplified) = self.simplified_summary() {
                        return simplified;
                    }
                }
                self.binop_summary(indent_level, lhs, rhs)
            }
            Node::Or(lhs, rhs) => self.binop_summary(indent_level, lhs, rhs),
        }
    }

    /// `CountingBinopReadFilter.getSummaryLineForLevel`, shared by AND and OR.
    fn binop_summary(&self, indent_level: usize, lhs: &Counting, rhs: &Counting) -> String {
        let indent = "  ".repeat(indent_level);
        if self.filtered_count == 0 {
            // The zero branch is written as a literal and keeps the space before the newline that
            // the non-zero branch drops.
            return format!("{indent}0 read(s) filtered by: {} \n", self.name());
        }
        let mut out = format!(
            "{indent}{} read(s) filtered by: {}\n",
            self.filtered_count,
            self.name()
        );
        // A child that rejected nothing is left out entirely.
        if lhs.filtered_count > 0 {
            out.push_str(&lhs.summary_line_for_level(indent_level + 1));
        }
        if rhs.filtered_count > 0 {
            out.push_str(&rhs.summary_line_for_level(indent_level + 1));
        }
        out
    }

    /// `getSummaryLineForLevelAllAndsSimplified`, or `None` where the Java returns "".
    ///
    /// The walk is an explicit stack pushing right then left, which is what puts the leaves in
    /// source order. Any OR or NOT anywhere in the tree abandons the whole flattening.
    fn simplified_summary(&self) -> Option<String> {
        let Node::And(lhs, rhs) = &self.node else {
            return None;
        };
        let mut out = String::new();
        // `Stack<CountingReadFilter>`: a `Vec` used as a stack, `push` and `pop` at its end.
        let mut unread: Vec<&Counting> = vec![rhs, lhs];
        while let Some(current) = unread.pop() {
            match &current.node {
                Node::And(l, r) => {
                    unread.push(r);
                    unread.push(l);
                }
                Node::Or(_, _) | Node::Not(_) => return None,
                Node::Leaf { .. } => out.push_str(&current.summary_line_for_level(0)),
            }
        }
        // No trailing newline: this is the last line a tool prints.
        out.push_str(&format!(
            "{} total reads filtered out of {} reads processed",
            self.filtered_count, self.total_count
        ));
        Some(out)
    }
}

/// Rebuild a tree from the name the reference gives it, so both sides compose the same thing.
///
/// The grammar is exactly what `getName` produces: `(a AND b)`, `(a OR b)`, `NOT a`, and a bare
/// class name for a leaf. Parsing the reference's own output rather than keeping a second list of
/// compositions is what stops the two from drifting.
pub fn parse(spec: &str, header: &SamHeader) -> Option<Counting> {
    let spec = spec.trim();
    if let Some(rest) = spec.strip_prefix("NOT ") {
        return Some(parse(rest, header)?.negate());
    }
    // `strip_prefix` gives the rest of the string, or nothing when the prefix is not there.
    if let Some(inner) = spec.strip_prefix('(').and_then(|s| s.strip_suffix(')')) {
        // Split on the operator at depth zero, which is the one this node owns.
        let bytes = inner.as_bytes();
        let mut depth = 0;
        for i in 0..bytes.len() {
            match bytes[i] {
                b'(' => depth += 1,
                b')' => depth -= 1,
                _ if depth == 0 => {
                    for (token, is_and) in [(" AND ", true), (" OR ", false)] {
                        if inner[i..].starts_with(token) {
                            let lhs = parse(&inner[..i], header)?;
                            let rhs = parse(&inner[i + token.len()..], header)?;
                            return Some(if is_and { lhs.and(rhs) } else { lhs.or(rhs) });
                        }
                    }
                }
                _ => {}
            }
        }
        return None;
    }
    // The header-dependent leaves carry a copy of the header, because the engine gives them one
    // through setHeader before wrapping them and they answer differently without it.
    let owned = header.clone();
    Some(match spec {
        "WellformedReadFilter" => Counting::leaf(spec, move |read| {
            crate::with_header::wellformed(read, &owned)
        }),
        "AlignmentAgreesWithHeaderReadFilter" => Counting::leaf(spec, move |read| {
            crate::with_header::alignment_agrees_with_header(read, &owned)
        }),
        // An unknown name gives nothing, and so does the whole parse.
        _ => Counting::leaf(spec, crate::by_name(spec)?),
    })
}
