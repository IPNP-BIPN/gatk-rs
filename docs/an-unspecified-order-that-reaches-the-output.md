# An unspecified iteration order that reaches the output

`AlignmentContextIteratorBuilder` collects the header's sample names with `Collectors.toSet()`.
`LocusIteratorByState` creates one per-sample state manager per element of that set, in iteration
order, and concatenates their elements in the same order to build every pileup.

So the element order of a multi-sample pileup is a `java.util.HashSet`'s iteration order over the
sample names. Any tool that prints a pileup, or that reduces over its elements in a way that is not
commutative, inherits that order in its output bytes.

## Why this is awkward

The order is **deterministic but unspecified**. `HashMap`'s documentation states that it makes no
guarantee about iteration order. So there is no contract to implement: the only description of the
behaviour is the implementation, and the implementation is `java.util`, which is GPL2.

That is the case htsjdk-rs decision 0013 already refused once, for `FloatingDecimal`, and which
`licence-compatibility-risk.md` records as this programme's critical risk. The OpenJDK Assembly
Exception grants permission to **link**, not to translate. Transcribing `HashMap` into an Apache 2.0
repository is the failure mode that document exists to prevent, and the provenance guard caught the
first version of `crates/gatk-engine/src/java_hash.rs`, which claimed exactly that provenance.

## What was done instead

The order is treated as an **observable of the pinned oracle**, in the same way a command line or a
fixture BAM is.

- `String.hashCode` is *specified*: its Javadoc gives the value as
  `s[0]*31^(n-1) + ... + s[n-1]` in `int` arithmetic. Computing that is implementing a published
  contract, and the guard classifies it as such.
- The bucket layout is *not* specified, so the conformance golden records, for each probed set of
  names, both the order the reference produced and each name's `String.hashCode`. The golden is the
  definition of the behaviour; the Rust code is a hypothesis that reproduces it. Where the two
  disagree, the measurement is right.

## What this costs

A standing obligation rather than a one-off. Any shape of sample-name set the suite does not probe
is unverified, and a divergence there would surface as a reordered pileup rather than as an error.
The probe therefore covers a set large enough to cross the growth point and a name whose hash is
negative, and the function refuses outright once a bucket grows past the length where nothing has
been measured.

Two consequences follow for the rest of the programme:

1. **This pattern will recur.** Any GATK behaviour that depends on the iteration order of a
   `HashSet` or `HashMap` over strings lands here. Each such case gets a probe, and none of them may
   be answered by reading OpenJDK.
2. **`Random` is the same problem, and worse.** Downsampling draws from `Utils.getRandomGenerator`,
   which is a `java.util.Random`. Its algorithm *is* specified in its Javadoc, unusually, so that
   one can be implemented from the contract. The G1.5 downsampling box stays open until that is
   done deliberately rather than incidentally.
