# When a dependency is cheaper than a port, and when it is not

`gatk-engine` depends on `noodles-fasta` for indexed FASTA access. This is the first third-party
implementation of a file format in this programme, and the rule it establishes is worth stating
once, because applying it wrongly would quietly weaken every claim downstream.

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

`noodles-fasta = "=0.61.0"`, an exact version rather than a caret range. A byte-identity claim
cannot float its dependencies: a patch release that changed an edge case would change the port's
answers with nothing in the diff to show for it. The same rule as `rev = "..."` for the two sibling
repositories.

## Licence

`noodles` is MIT, `gatk-rs` is Apache 2.0, and MIT is compatible with inclusion in an Apache-2.0
work. `tools/audit/provenance.py` checks ported *symbols*, which a dependency is not; this document
is the record for the dependency itself.
