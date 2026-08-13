# A truncated alignment that looks finished

`GATKVariantContextUtils.leftAlignAndTrim` slides an indel left through the reference. It does not
read the whole contig: it takes a slice, ten bases before the variant to start with, and widens it
while the shift keeps landing on the slice's edge.

```java
for (int leadingBases = Math.min(maxLeadingBases, 10); leadingBases <= maxLeadingBases; leadingBases = Math.min(2*leadingBases, maxLeadingBases)) {
    ...
    } else if (shifts.getLeft() == variantOffsetInRef && leadingBases < maxLeadingBases) {
        continue;
    }
```

The second half of that condition is the bound: once the window is at `--max-leading-bases`, the
loop stops widening and returns whatever the last pass produced. So a record can come out **shifted
as far as the argument allowed and not as far as the reference would allow**, and it is returned
the same way a fully aligned record is.

The `left-align-and-trim` suite holds one deletion at four windows:

| `--max-leading-bases` | where it lands | left aligned |
|---|---|---|
| 2 | `chr1:15-16` | no |
| 10 | `chr1:49-50` | no |
| 20 | `chr1:39-40` | no |
| 1000 | `chr1:30-31` | yes |

All four are `VariantContext`s with nothing to tell them apart. Two runs of the same tool over the
same data with different `--max-leading-bases` produce different, equally silent files, and a
downstream comparison of two such files reports a difference in the data rather than in the
arguments.

## Why the port does not fix it

The bound is deliberate: the whole point of the widening loop is to avoid reading megabases of
reference for a variant that will move three bases. Removing it would change output bytes, and
byte-identity with GATK 4.6.2.0 is the thing this programme exists to prove. The truncation is
reproduced exactly, and `left_align_and_trim` returns the reference's record.

What the port does not reproduce is the **silence**, because the silence is not in the bytes.
`left_align_and_trim_reporting` returns the same record with an `Alignment` beside it:

- `Complete`, the shift stopped before the edge, or one more base of window would not move it;
- `WindowExhausted`, the window is what stopped it, and a larger one would move it further;
- `NotAttempted`, a non-indel or a window of zero or less.

Nothing about the record changes. The conformance suite asserts that on all nineteen calls of the
golden: the record from the reporting function is the record from the plain one, and the report
classifies the three truncated calls, the four never attempted, and the rest.

## The measurement inside the fix

Telling `Complete` from `WindowExhausted` is not `shifts.getLeft() == variantOffsetInRef`. A window
of exactly the distance the record wants to walk lands **on** the edge and is complete: the
deletion above with `--max-leading-bases 7` reaches `chr1:10-11`, the same place the window of 1000
reaches, with the shift equal to the offset. The first version of this reporting got that wrong and
its own unit test caught it.

The honest test is to offer one more base and see whether the shift grows. Left alignment slides
base by base, so a record that can still move moves by at least one when given one more base, and a
record already at the start of the contig cannot be given one. That is one extra pass over one
slightly wider slice, and only in the ambiguous case.

## What this costs

A caller that wants the guarantee has to ask for it, and then decide what to do with the answer.
Nothing in the tool chain asks yet: `LeftAlignAndTrimVariants` is not ported, and when it is, it
will pass `min(maxLeadingBases, distanceToLastVariant - 1)`, which truncates for its own reasons
and produces `WindowExhausted` routinely and correctly. The report is there so that a future caller
can distinguish "as left as it goes" from "as left as we paid for", which the reference cannot.
