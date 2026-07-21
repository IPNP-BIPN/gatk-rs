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
