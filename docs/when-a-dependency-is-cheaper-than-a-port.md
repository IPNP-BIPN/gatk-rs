# When a dependency is cheaper than a port, and when it is not

`gatk-engine` had two third-party readers of a file format: `noodles-fasta` for indexed FASTA
access and `noodles-bam` for parsing the `.bai`. They were the first in this programme, and the
rule they established is worth stating once, because applying it wrongly would quietly weaken
every claim downstream.

Both are gone. The engine now reaches a FASTA and a `.bai` through `htsjdk-bam`, this programme's
own port, and no third-party implementation of a file format is left in any of the three
repositories. The rule below still stands; the two applications of it did not survive being
measured, and that is the more useful half of this document.

## The rule

**Depend where the bytes are unambiguous. Port where the semantics are the reference's.**

A `.fai` index is five numbers a line, and seeking to `offset + (position / line_bases) *
line_width + (position % line_bases)` has one right answer. Two implementations of that either
agree or one of them is broken, and the conformance suite says which.

What a GATK tool actually sees is a different question, and it is not the file's bytes.
`CachingIndexedFastaSequenceFile` defaults to `preserveCase = false` and `preserveIUPAC = false`,
so every query comes back upper-cased with every IUPAC ambiguity code replaced by `N`. Measured in
`crates/gatk-engine/tests/data/reference_query.txt.gz`:

| in the file | what the reference returns |
|---|---|
| `acgtNNNNacgt` | `ACGTNNNNACGT` |
| `ACGTRYKMSWBD` | `ACGTNNNNNNNN` |

Soft-masking is erased and ambiguity codes are flattened. A port that returned what any FASTA
reader gives would differ from the reference at every soft-masked or ambiguous position in a
genome, and those are not rare: roughly half of the human reference is soft-masked. That half of
the split was right and has not changed: the transformation is ported here and measured against
the reference's own answers.

## Why the FASTA half was taken back

The premise that failed is "the bytes are unambiguous". A query is not answered by reading bases
until the right ones arrive: `getSubsequenceAt` **seeks**. The first base's byte offset is computed
from the `.fai`'s bases-per-line and bytes-per-line columns, and every line boundary the query
crosses is a jump over a terminator whose length is the difference between those two columns.
Nothing ever looks for a newline. A CRLF file is therefore read correctly only because the index
says the terminator is two bytes, and a `.fai` that disagrees with its own file is read wrongly and
silently.

That is not a hypothesis. The first version of htsjdk-rs's suite wrote the `.fai` by hand, got the
CRLF offset one byte wrong, and the oracle answered `chr1:6-7` with a terminator byte among the
bases, reported as an answer and not as an error. So the file's bytes do not determine the answer;
the index's arithmetic does, and that arithmetic is htsjdk's. A reader that scans for newlines
agrees on every well-formed file and disagrees on the rest, which is the shape of divergence this
programme exists to refuse.

Two of the bounds are htsjdk's own and one of them looks wrong: the malformed-query test is
`start > stop + 1`, so an **empty** query is legal and answers with no bases, while `start > stop`
by two is refused. Past the end of a contig is refused by the *index*, so it is the `.fai`'s size
column that decides and not the file's length.

The measurement lives in htsjdk-rs's `indexed-fasta` suite: sixteen answered queries and three
refusals over three files, with the `.fai` written by `FastaSequenceIndexCreator` and carried in
the golden, so the creator's arithmetic is pinned together with the reader's.

## The second application: the `.bai`

`ReadsDataSource` splits the same way, and the line falls in a place worth naming. Everything that
decides *which records come back* is ported into `crates/gatk-engine/src/reads.rs`, because each of
those is htsjdk's or GATK's and not the format's:

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

The parse underneath them was the part a dependency answered for, and a bin's chunk list really is
what the format says it is: this half was not taken back on correctness. It was taken back on
price. `noodles-bam` dragged `noodles-sam` and `rayon` into the build for what was, here, a `.bai`
parser; `htsjdk_bam::index::read_bai` was already written and already measured, by the
`textual-index` suite, which parses six `.bai` files and reprints each one through
`TextualBAMIndexWriter`'s format. Dropping both dependencies removes 43 packages from
`Cargo.lock`, and the reader that answers is now the one whose write side is compared byte for
byte.

Records are decompressed by `htsjdk-bgzf` and decoded by `htsjdk-bam`, this programme's own ports,
because what a record *is* is htsjdk's decision, and that is the "reading is itself a decision"
case below.

## Where this does not apply

- **The write path.** Byte-identity lives there: which deflate level, which tag integer width,
  which order. htsjdk-rs ports those, and a dependency would replace a measured property with a
  hope. No third-party format implementation is on the read or the write path of any of the three
  repositories.
- **Formats whose reading is itself a decision.** htsjdk's SAM reader refuses `RNAME is not
  specified but flags indicate mapped`; its BAM tag codec picks integer widths from the *value*
  rather than the declared type. A reader that accepts more, or normalises differently, silently
  changes which records exist. Those stay ported.
- **CRAM.** It is tempting for exactly the reason it is dangerous: the sub-project is large, and a
  dependency would swap a byte-identity claim for a bio-identity one without anything failing. If
  it is ever taken, it is taken explicitly and the status drops to bio-identical in
  `docs/STATUS.md`.

## Pinning

The two that existed were pinned as `noodles-fasta = "=0.66.0"` and `noodles-bam = "=0.94.0"`,
exact versions rather than caret ranges, and the next one would be too. A byte-identity claim
cannot float its dependencies: a patch release that changed an edge case would change the port's
answers with nothing in the diff to show for it. That is the same rule as `rev = "..."` for the two
sibling repositories, which is how `htsjdk-bam` itself arrives here.

## Licence

`noodles` is MIT and `gatk-rs` is Apache 2.0, which was compatible while it lasted.
`tools/audit/provenance.py` checks ported *symbols*, which a dependency is not, so a document like
this one is the only record a dependency leaves. That is a further reason to prefer the port where
the choice is close: a port is audited by a gate, and a dependency is audited by whoever remembers
to read the file.
