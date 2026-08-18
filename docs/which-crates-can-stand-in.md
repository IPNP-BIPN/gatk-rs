# Which crates can stand in for the reference, measured rather than argued

Prompted by a read of [awesome-rust](https://github.com/rust-unofficial/awesome-rust). Rather than
reason about which crates could replace ported code, each candidate was run against the pinned oracle
container over a corpus and compared answer by answer. The probes live outside the tree; the numbers
are below, and they are the reason each verdict is what it is.

The rule they serve is `docs/when-a-dependency-is-cheaper-than-a-port.md`: depend where the bytes are
unambiguous, port where the semantics are the reference's. This note is that rule's evidence.

## Regular expressions, which is what `-se`, `-xl-se` and ClipReads compile

47 patterns by 27 inputs, boolean `find()`, against `java.util.regex`.

| engine | patterns accepted differently | find() answers differing | compile messages differing |
|---|---|---|---|
| `regex` 1.x | 4 | 17 | 7/7 |
| `regex-lite` | 8 | 10 | 7/7 |
| `fancy-regex` | 1 | 16 | 7/7 |
| `regress` | 2 | 74 | 7/7 |
| `onig`, default syntax | 2 | 17 | 6/7 |
| `onig`, `Syntax::java()` | 1 | 12 | 7/7 |
| `pcre2` | 1 | 21 | 7/7 |
| **the port, `gatk_engine::java_regex`** | **13, all of them refusals** | **0** | **1/7, the index only** |

The causes, all of them silent:

- `$` before a final line terminator: Java's `^s1$` matches `"s1\n"`, `regex` and `fancy-regex` do not;
- `.` and `\r`: Java excludes four line terminators, the crates exclude only `\n`;
- `\d`, `\w`, `\p{Alpha}`, `\p{Lower}`: Java is ASCII without `UNICODE_CHARACTER_CLASS`, the crates are Unicode. Oniguruma's `ONIG_OPTION_*_IS_ASCII` bits are not reachable through `rust-onig` alongside `Syntax::java()`;
- **possessive quantifiers**: `regex` reads `.*+` as `(.*)+` and flips a match rather than refusing to compile.

None of them produces Java's `PatternSyntaxException` text, which SelectVariants lets reach the user unwrapped. Fixing the three configurable causes needs a parser, because the tokens must not be rewritten inside a character class, and that parser is most of what the port already is.

**Verdict: keep the port.** A crate would replace about 120 lines of matcher and add a dependency on the byte path.

## Compression, which is the deflate of X.4 (#75)

Five real inputs (a golden, 136 KB of a covariates dump, 1 MB of zeros, 500 KB pseudo-random, 200 KB of DNA text) by levels 0 to 9, raw and zlib-wrapped, against `java.util.zip.Deflater`, compared as the md5 of the bytes.

`java.util.zip.Deflater` **is** zlib, so the question is not which crate compresses well but which backend is zlib.

| backend | byte-identical |
|---|---|
| `flate2`, miniz_oxide (the default) | 1/100 |
| `libdeflater` | 2/100 |
| `flate2`, `zlib-rs` backend | 12/100 |
| `flate2`, C zlib backend | **92/100** |
| C zlib through Python, as a control | 92/100 |

Why each one misses:

- **miniz_oxide** is not zlib. It reimplements miniz, with its own lazy-matching and block-splitting heuristics: same format, other bytes;
- **libdeflater** is a different algorithm and a better one. Level 1 on the covariates dump: 13266 bytes against zlib's 14873. Compressing better is compressing differently;
- **`zlib-rs`** descends from **zlib-ng**, not from stock zlib. zlib-ng replaced level 1 with `deflate_quick` and retuned the middle levels: 18574 bytes where zlib writes 14873. A Rust rewrite of zlib is not a rewrite of *this* zlib;
- **C zlib** matches, and all eight remaining misses are at **level 0**, on the large inputs only. Java emits stored blocks of 65535 bytes (1 MB of zeros: 16 headers, 1000080 bytes); zlib 1.2.12 emits about thirty (1000150). `deflate_stored` was rewritten between zlib versions and the JDK bundles its own, so byte-identity depends on the zlib **version** as well, and level 0 is where the versions part.

htsjdk writes BGZF at level 5 by default, and level 5 matches exactly.

And the port beats all of them: **`gkl-deflate` (htsjdk-rs) is 45/45** against `java.util.zip.Deflater` over the same corpus at levels 1 to 9, which is the whole of what BGZF asks for. It also carries the Intel GKL flavour, checked against hashes the real library produced in the pinned container. #75 is open only for igzip at levels 1 and 2.

**Verdict, revised twice.** A first reading said no crate could carry the write path, measured against pure-Rust backends only. A second said the choice was between linking C zlib and porting `deflate`; the port already exists. What the crates add is nothing, and what they cost is a dependency on the byte path. A first run over four tiny inputs had miniz_oxide at 36/36, which is the trap worth naming twice: small inputs mostly agree.

## BGZF writers, since the format is a second layer of decisions

`bgzf` 0.5 from Fulcrum Genomics, against htsjdk's `BlockCompressedOutputStream` over the same five inputs at levels 1, 5, 6 and 9, compared as the md5 of the whole file:

| writer | byte-identical |
|---|---|
| `bgzf` 0.5 (Fulcrum Genomics) | 0/20 |
| the port, `htsjdk-bgzf` | **20/20** |

It compresses with libdeflater, and it compresses *better*: 1 MB of zeros at level 1 comes to 1666 bytes against htsjdk's 5067, the DNA text at level 5 to 55306 against 60588. Compressing better is writing other bytes, so the BAM is not the same BAM.

**Where it is legitimate: reading.** What an inflate returns is defined by the stream and not by the implementation, so a third-party BGZF *reader* is plumbing under the rule in `docs/when-a-dependency-is-cheaper-than-a-port.md`. A writer never is.

## Collections, which several ports emulate by hand

- `indexmap` reproduces `LinkedHashSet` and `LinkedHashMap` iteration order exactly, **including re-insertion keeping the original position**. Verdict: a legitimate dependency, and it would replace order-preserving maps written by hand in a dozen modules;
- `BTreeSet<String>` does **not** reproduce `TreeSet<String>`: Java orders by UTF-16 code units and Rust by UTF-8 bytes, and they disagree for every supplementary character (2 of 289 pairs in the probe). `gatk_engine::java_hash::compare_strings` is the fix, and one port was found using the wrong order: #260.

## Special functions, against Apache Commons Math

| function | `statrs` exact | `puruspe` exact |
|---|---|---|
| `Gamma.logGamma` | 1/9 | 0/9 |
| `Gamma.digamma` | 0/9 | not offered |
| `Gamma.regularizedGammaP` | 4/9 | 1/9 |
| `NormalDistribution.cumulativeProbability` | 1/7 | not offered |
| `Beta.regularizedBeta` | 2/5 | 0/5 |

Worst observed distance 490000 ULP. **Verdict: keep the ports.** Rust's own `f64::exp` and `f64::ln` did match `Math.exp`/`Math.log` and `FastMath.exp`/`FastMath.log` at 0 ULP on the five points tried, which is worth a wider probe before anything rests on it.

## Float formatting, which every report and every annotation goes through

`Double.toString` against Rust's `{}` (18 of 27 differ) and against `ryu` (13 of 27 differ): Java switches to exponent notation at different thresholds and always writes a `.0`. `String.format("%.3f")` against Rust's `{:.3}` differs on 4 of 27, including `Infinity` versus `inf`.

The port, `gatk_engine::java_format::format_decimals`, differed from Java on 2 of 27 when this was
written: `1e300` and `Double.MAX_VALUE`, where Java writes the shortest representation's digits padded
with zeros and the port wrote the exact binary value. Filed as #262 and **since fixed**: the
`string-format` golden showed the difference was wider than those two magnitudes, reaching `2.675` and
`1.005` in the ordinary range, and the conversion now rounds the digits `Double.toString` produces.
`Double.toString` itself still differs on fifteen values of a 1059-value corpus, which the
`double-to-string` golden pins and #399 tracks.

## Not probed, with reasons

- `tch` and the ML stack: settled in `docs/the-ml-surface-cannot-be-bit-identical.md`, byte-identity unreachable;
- `wgpu`, `cudarc`, `rust-gpu`: nothing to compare until a kernel exists, feeds #76;
- `rayon`: a scheduling question rather than an output question, belongs with #85's determinism gates;
- `noodles-vcf`: reading a VCF is a decision and not plumbing, per the rule above.

## Tooling worth adopting, where there is no byte path to risk

`cargo-fuzz` with `arbitrary` and `cargo-llvm-cov` for #84, `proptest` for port invariants, `criterion` or `divan` for milestone S, `cargo-nextest` for the 29-job matrix, and `cargo-deny` to mechanise the exact-pin invariant the repository holds by hand.



