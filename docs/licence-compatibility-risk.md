# R11: the reference implementation's licence may be incompatible with the port

The plan's risk register lists ten risks. This is the eleventh, and it is **critical**. It was
not anticipated, it has already fired once, and the failure it produces is not a wrong number —
it is a licence violation in a public repository.

## What happened

Two decisions in htsjdk-rs, in the order they were written:

- **0013**: `FloatingDecimal` cannot be ported. It is a `java.base` class, `java.base` is GPL2,
  and the OpenJDK Assembly Exception grants permission to **link**, not to translate and
  relicense. The 112 remaining `FormatUtil` divergences and one SAM-text float divergence are
  therefore reclassified from "not yet ported" to **licence-blocked**.
- **0014**: an audit run immediately afterwards found that `crates/jmath/src/exp.rs` was
  *already* an operation-by-operation transcription of
  `src/hotspot/cpu/x86/macroAssembler_x86_exp.cpp`, whose header reads **"GNU General Public
  License version 2 only"** with no Classpath Exception. It had been merged and published under
  MIT. It is removed; `Math.exp` is now unported.

The gap between writing 0013 and finding the violation was minutes. The gap between committing
the violation and finding it was eleven commits.

## Why it is critical rather than medium

The three ported projects are permissively licensed: htsjdk MIT, Picard MIT, GATK Apache 2.0.
That is what made the licence question feel settled. **It is not settled, because all three run
on a GPL2 JVM, and the port has to reproduce the JVM's behaviour wherever that behaviour reaches
an output byte.**

Places where it already does, or predictably will:

| behaviour | JDK source | status |
|---|---|---|
| `Math.exp`, `Math.pow` intrinsics | HotSpot, **GPL2 only** | exp withdrawn, pow never written |
| `Double.toString` / `DecimalFormat` rounding | `java.base`, GPL2+CPE | licence-blocked, decision 0013 |
| `String.hashCode` iteration order | `java.base`, GPL2+CPE | reachable wherever a `HashSet` orders output |
| `Arrays.sort` tie-breaking | `java.base`, GPL2+CPE | reachable in every sort-order claim |
| `TreeMap` iteration | `java.base`, GPL2+CPE | already relied on by `Histogram` |

The last two are worth dwelling on. htsjdk's `Histogram` iterates a `TreeMap`, and every
statistic depends on that order. The port reproduces the *order* — sorted by key — which is a
documented property of a sorted map and not a transcription of anything. That is fine. It would
**not** be fine to transcribe `TreeMap`'s balancing to reproduce some order-dependent detail,
and the distinction between "implement to a documented property" and "translate the source" is
now the load-bearing one.

## The rule

**A symbol is portable into this program only if it comes from the pinned htsjdk, Picard or
GATK clones.** Anything reached through `java.lang`, `java.util`, `java.text`, `java.math` or
HotSpot is GPL2 and is not portable, whatever the Classpath Exception says about linking.

Where such behaviour is observable in an output byte, there are four options and only four:

1. **Establish an independent property and implement to it.** This is what htsjdk-rs decision
   0006 did for `Math.log`: it asked whether the intrinsic was correctly rounded, found that it
   was, and implemented correct rounding. The result is exact and owes nothing to the GPL2
   source. Where such a property exists this is the right answer, and it is worth looking for it
   before concluding a function is blocked.
2. **Obtain permission.** For a small self-contained function this is cheap to ask.
3. **Change the oracle.** From JDK 19, `Double.toString` is Schubfach and *is*
   shortest-round-trip, which Rust already produces — so cause A of decision 0011 disappears
   against a JDK 19+ oracle. This makes the oracle's JDK version a **licensing** decision as
   well as a fidelity one, which the plan does not currently treat it as.
4. **Quarantine.** Report the affected values as bio-identical, with the exact list committed.

## The label is more dangerous than the block (2026-08-04)

The table above lists `DecimalFormat` rounding as licence-blocked. G1.9.3 needed it, went looking
for option 1, and found that **most of what the label covered was never blocked**.

Two things had been run together:

- *the shortest decimal representation of a double*, which is a genuinely hard problem and is what
  `FloatingDecimal` exists for. Rust's `{:e}` already supplies it, so nothing needs porting;
- *the exact decimal expansion of a double*, which htsjdk-rs decision 0013 said also required
  `FloatingDecimal`. It does not. Every finite double **is** a finite decimal: `m * 2^e` is `m`
  doubled `e` times when `e` is non-negative, and `m * 5^-e` with the point moved when it is not,
  because `2^e = 5^-e / 10^-e`. Multiplying a decimal digit string by five is one pass over its
  digits. Thirty lines of schoolbook arithmetic on the bits of an IEEE 754 double, and a
  translation of nothing.

The second is what `DigitList.shouldRoundUp` needs to break a tie, and supplying it closed **44 of
the 112** quarantined `FormatUtil` divergences with none introduced (htsjdk-rs #72, decision 0026).
The 68 that remain are one cause, and it is the first bullet: Java 17's pre-Schubfach digit
generation, which option 3 above already covers.

### How wide that one cause is (2026-08-19)

The 68 are the divergences the *suites* reach. A corpus that goes looking finds the same failure mode
across the whole exponent range: `Double.toString` measured over 1059 deterministic values -- 27 named
ones, 32 powers of two, and a thousand bit patterns from splitmix64 seeded at one -- differs on
**fifteen**, 1.4 per cent, and in every one **Java emits more digits than the shortest**. Three shapes:

| | the reference | this port |
|---|---|---|
| the smallest subnormal | `4.9E-324` | `5.0E-324` |
| `1e23` | `9.999999999999999E22` | `1.0E23` |
| eleven values needing sixteen digits | `7.2911220195563975E-304` | `7.291122019556398E-304` |

Two things follow. It is **not confined above 2^53**: the smallest subnormal is in the list, so a
claim that the region is unreachable has to be made per call site rather than per magnitude. And every
one of the fifteen **parses back to the same double**, so this is a rendering difference and not a
value difference -- which is why nothing caught it until a golden printed the renderings themselves.

The fifteen are committed as the `double-to-string` suite and named in its conformance test with both
renderings and an asserted count, which is option 4 carried out rather than proposed: a sixteenth
value diverging, one of the fifteen changing, and one of them being fixed all fail a test. Option 1
stays closed for this one -- htsjdk-rs decision 0013 blocks `FloatingDecimal` itself, and a clean-room
implementation of the *specification* would agree with this port rather than with the oracle, which
is what Android's `RealToString` demonstrated. Option 3 is the only route that removes them, and it is
the same decision the 68 wait on.

**The lesson is about the label, not the licence.** Once work is marked licence-blocked it stops
being examined, and the mark covers whatever was nearby when it was applied. Option 1 was written
down and correct; what failed was that nobody re-ran it after the classification. Anything
currently carrying that status is worth one more look, and the question to ask is narrow: *which
specific fact does the output depend on, and is that fact reachable without the source?*

The same re-examination reopened `AllelePseudoDepth` itself (G1.9): it was refused on `Math.exp`,
and both values it emits leave through a formatter at two and four decimals, about twelve orders of
magnitude coarser than the 1 ulp that was being worried about.

**What this does not reach.** `Math.exp` bit-for-bit stays out of scope, and for a different reason
that this finding does not weaken: there is no specification to implement against, so recovering
its bits by black-box measurement would be reverse engineering toward a functional copy rather than
implementing to a property. See htsjdk-rs #71. Option 1 needs a property to exist; where none does,
the option is not available.

## Mitigation now in place

`htsjdk-rs/tools/audit/provenance.py` checks every `Ported from` claim in the tree against an
explicit allow/deny list of source licences, and CI runs it on every push. It is verified in
both directions: clean on the current tree, and it catches the exact file that caused this.

The same guard belongs in picard-rs and gatk-rs before either ports a symbol whose source is not
obviously in the pinned clones.

## What this does to the estimate

Nothing directly, and something indirectly. No line count changes. But `Math.pow`'s 2,220 lines
were sized as a large porting job in decision 0007. They are not a porting job at all — they are
unportable, and the honest entry for `pow` is now "quarantined, pending one of the four options
above" rather than "large".

The same reclassification may apply to other work currently sized as difficult. Anything whose
difficulty is "reproduce the JVM exactly" should be re-examined under this rule before it is
scheduled.
