# Commenting standard

This programme reimplements GATK, Picard and htsjdk in Rust and claims **byte-identical** output.
A reader who wants to check that claim has to compare two things: what the Java does and what the
Rust does. Most of them can read the Java. Not all of them can read the Rust, and none of them
should have to read it twice to find out *why* a line is written the way it is.

So the code is commented for someone who does not know Rust.

## The three questions

Every item that is not self-evident answers three questions, in this order:

- **What** it computes, in the domain's own words. "The index of the largest value", not "the loop
  variable after the loop".
- **How** it computes it, when the method is not obvious from the name. One or two sentences.
- **Why** it is written *this way and not the obvious way*. This is the one that carries the value
  of the whole exercise, because in a byte-identity port the obvious way is usually wrong.

The **why** is the point. A comment that says `// add one to i` is worse than no comment. A comment
that says `// strictly greater, so the first of two equal maxima wins, and that tie is reachable`
is why the file exists.

## Levels

**Module** (`//!` at the top of the file). Which reference class this is, which version, and the
handful of things about it that would surprise a careful reader. This level already exists across
the codebase and is not changed by this standard.

**Item** (`///` above a function, type or constant). The three questions. Always name the reference
symbol it ports, so a reader can put the two side by side.

**Inline** (`//` inside a body). Two kinds, both wanted:

- *Intent*: what this step is for, and what the reference does at the same point.
- *Mechanics*: what a piece of Rust syntax means, **where the idiom is not guessable from Java**.
  `.iter().enumerate().skip(n)` deserves a sentence. `let x = 1;` does not.

## Mechanics worth explaining, once per file at most

These are the constructs a Java reader trips on. Explain them the first time a file uses one, and
not again:

| Rust | What it means to a Java reader |
|---|---|
| `&[T]` | a borrowed, read-only view of an array; the caller keeps ownership |
| `&mut T` | a borrowed, writable reference; only one may exist at a time |
| `Option<T>` | "a value or nothing", checked by the compiler; Java's `null` without the crash |
| `Result<T, E>` | "a value or an error"; how this codebase models a reference exception |
| `let Some(x) = y else { ... }` | unwrap or take the early exit; there is no Java equivalent |
| `.iter().map(..).collect()` | Java's stream, map, collect |
| `?` | propagate the error to the caller; Java's `throws` at a single call site |
| final expression, no `return` | the last expression of a body is its value |
| `f64`, `i32`, `i64`, `usize` | `double`, `int`, `long`, and an unsigned index type Java lacks |

## What not to write

- Do not restate the type. The compiler already checked it and the comment will rot.
- Do not narrate control flow that reads plainly. `// loop over the samples` above
  `for sample in samples` is noise.
- Do not explain the same idiom twice in one file.
- Do not write a comment you cannot defend from the reference. If the reason is "I think this is
  faster", say that; if the reason is a line of Java, quote the line.

## Quoting the reference

Where a comment turns on what the Java does, quote the Java, in a fenced block for an item comment
and inline for a one-liner. The quote is the evidence. A paraphrase invites a reader to trust the
port's reading of the reference, which is exactly what the port is not entitled to ask.

## Density, and why it is measured rather than mandated

`tools/audit/comment_density.py` reports, per file, the ratio of comment lines to code lines. It
does **not** fail a file for being below the bar: a generated table of five thousand constants
should not be commented line by line, and a mechanical accessor needs nothing.

What it does fail is **regression**. A file listed in `tools/audit/commented.txt` has been brought
to this standard deliberately, and the guard fails if its density drops. That way the work is
ratchet-shaped: nothing already done can quietly come undone, and the remaining files are a
measured list rather than a feeling.

Run it with no arguments for the full table:

```
python3 tools/audit/comment_density.py
```

## Order of work

The list in `tools/audit/commented.txt` grows one tranche at a time. The order is by how much a
reader needs the file to understand anything else:

1. the engine's numeric primitives (`math_utils`, `histogram`, `java_hash`, `fragment`)
2. the likelihood matrix and its allele and sample axes
3. the annotations, which are the widest surface and the least surprising individually
4. the readers and writers, where the surprises are about file formats rather than arithmetic
5. everything else

`htsjdk-rs` and `picard-rs` carry their own copies of this file and their own lists.
