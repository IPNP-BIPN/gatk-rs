# When a dependency is cheaper than a port, and when it is not

`gatk-engine` depends on `noodles-fasta` for indexed FASTA access and on `noodles-bam` for parsing
the `.bai`. These are the first third-party implementations of a file format in this programme,
and the rule they establish is worth stating once, because applying it wrongly would quietly
weaken every claim downstream.

## The rule

**Depend where the bytes are unambiguous. Port where the semantics are the reference's.**

A `.fai` index is five numbers a line, and seeking to `offset + (position / line_bases) *
line_width + (position % line_bases)` has one right answer. Two implementations of that either
agree or one of them is broken, and the conformance suite says which. Porting htsjdk's copy of it a
second time buys no measurable property.

What a GATK tool actually sees is a different question, and it is not the file's bytes.
`CachingIndexedFastaSequenceFile` defaults to `preserveCase = false` and `preserveIUPAC = false`,
so every query comes back upper-cased with every IUPAC ambiguity code replaced by `N`. Measured in
`crates/gatk-engine/tests/data/reference_query.txt.gz`:

| in the file | what the reference returns |
|---|---|
| `acgtNNNNacgt` | `ACGTNNNNACGT` |
| `ACGTRYKMSWBD` | `ACGTNNNNNNNN` |

Soft-masking is erased and ambiguity codes are flattened. A port that returned what any FASTA
reader gives, `noodles` included, would differ from the reference at every soft-masked or ambiguous
position in a genome, and those are not rare: roughly half of the human reference is soft-masked.

So the dependency provides the plumbing and the transformation on top of it is ported and measured.
The suite compares the *whole* answer, which is what makes the split safe: if `noodles` and htsjdk
ever disagree about a line boundary, an offset or an edge, the golden fails and the difference gets
ported rather than inherited.

## The second application: the `.bai`

`ReadsDataSource` splits the same way, and the line falls in a place worth naming. The `.bai`
bytes are parsed by `noodles-bam`: a bin's chunk list is what the format says it is. Everything
that decides *which records come back* is ported into `crates/gatk-engine/src/reads.rs`, because
each of those is htsjdk's or GATK's and not the format's:

| ported | why it is not plumbing |
|---|---|
| `GenomicIndexUtil.regionToBins` | 1-based in, decremented before shifting; htsjdk's own comment calls this "suspicious" and keeps it |
| `LinearIndex.getMinimumOffset` | out-of-range windows mean *no* constraint, not an empty result |
| `Chunk.optimizeChunkList` | drops chunks below the minimum offset, coalesces chunks whose pointers exactly touch |
| `QueryInterval.optimizeIntervals` | merges **abutting** intervals, so two adjacent `-L` arguments return a spanning read once rather than twice |
| `BAMQueryMultipleIntervalsIteratorFilter` | stateful and single-pass: the interval index only advances, and the traversal *stops* once every interval is behind the record |
| `AbstractBAMFileIndex.getStartOfLastLinearBin` | the **last** reference's last entry, not the largest entry |

The filter is where the argument for porting is strongest. It special-cases an unmapped read that
carries its mate's coordinate to `end = start`, because `getAlignmentEnd()` is `0` for anything
with the unmapped flag; without it, every mate-placed unmapped read would sort before every
interval and be invisible to every query. A generic "does this record overlap this interval"
would drop them silently, and a caller counting reads over a region would be wrong by however many
half-mapped pairs the region holds.

Records are decompressed by `htsjdk-bgzf` and decoded by `htsjdk-bam`, this programme's own ports,
rather than by `noodles-bam`'s reader: what a record *is* is htsjdk's decision, and that is the
"reading is itself a decision" case below.

## Where this does not apply

- **The write path.** Byte-identity lives there: which deflate level, which tag integer width,
  which order. htsjdk-rs ports those, and a dependency would replace a measured property with a
  hope. `noodles` is not on the write path of any of the three repositories.
- **Formats whose reading is itself a decision.** htsjdk's SAM reader refuses `RNAME is not
  specified but flags indicate mapped`; its BAM tag codec picks integer widths from the *value*
  rather than the declared type. A reader that accepts more, or normalises differently, silently
  changes which records exist. Those stay ported.
- **CRAM.** It is tempting for exactly the reason it is dangerous: the sub-project is large, and a
  dependency would swap a byte-identity claim for a bio-identity one without anything failing. If
  it is ever taken, it is taken explicitly and the status drops to bio-identical in
  `docs/STATUS.md`.

## Pinning

`noodles-fasta = "=0.61.0"` and `noodles-bam = "=0.92.0"`, exact versions rather than caret ranges.
A byte-identity claim cannot float its dependencies: a patch release that changed an edge case
would change the port's answers with nothing in the diff to show for it. The same rule as
`rev = "..."` for the two sibling repositories.

`noodles-bam` drags `noodles-sam` and `rayon` into the build for what is, here, a `.bai` parser.
That is a build cost and not a correctness one, and it is the honest place to note that the day
this crate needs anything else from `noodles-bam` the answer is still no: its record reader is on
the wrong side of the line above.

## Licence

`noodles` is MIT, `gatk-rs` is Apache 2.0, and MIT is compatible with inclusion in an Apache-2.0
work. `tools/audit/provenance.py` checks ported *symbols*, which a dependency is not; this document
is the record for the dependency itself.
