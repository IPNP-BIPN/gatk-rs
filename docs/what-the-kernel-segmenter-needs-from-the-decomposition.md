# What the kernel segmenter needs from the decomposition

`CalculateContamination` reaches a singular value decomposition, and #230 asked what has to be true
of it: byte identity through the decomposition, or only through what the segmenter does with it.

The `kernel-segmentation` golden measured the decomposition itself. This is the follow-up question,
and it is answerable from the reference's source plus a small amount of algebra, without another
oracle run.

## The decomposition leaves the segmenter as five integers

`KernelSegmenter.findChangepoints` returns `List<Integer>`. `ContaminationSegmenter` turns those
indices into `SimpleInterval`s by looking up sites at `changepoints.get(n) + 1`, and
`ContaminationModel` consumes the intervals. Nothing else in the chain touches the decomposition:
`SingularValueDecomposition` appears once in the whole tool, at `KernelSegmenter.java:224`.

So the tool's output depends on the decomposition **only through the changepoint indices**. No
singular value and no `U` entry is carried into a printed number.

## Every cost sees the reduced matrix only through row inner products

The reduced observation matrix is

```
Z = K_reduced * (U * diag(1 / (sqrt(s) + 1e-10)))
```

and two functions consume it. `calculateKernelApproximationDiagonal` takes `||Z_i||^2`.
`calculateSegmentCost` accumulates `D` from those norms and `V` from `Z_tau . W`, where `W` is a
running sum of rows of `Z`. `calculateWindowCosts` slides the same recurrences.

Both are functions of `Z Z^T` alone. That has a consequence worth stating precisely, because it
retires the objection this port raised when the golden landed:

**The arbitrary basis does not matter.** The `U` columns spanning the null space are not determined
by the matrix -- any orthonormal basis of it is as valid -- and the earlier reading was that this is
fatal, since those are exactly the columns multiplied by `1e10`. It is not fatal. A different basis
is `U -> U Q` with `Q` orthogonal inside the degenerate block, the scaling is *constant* inside that
block, so `Q` commutes with it and `Z -> Z Q`. Then `Z Q (Z Q)^T = Z Z^T`. The freedom cancels
exactly, before any rounding.

Measured on the golden's own decomposition of the two-step series, with a random orthonormal basis
of its two-dimensional null space and a stand-in reduced matrix:

| perturbation | relative change in `Z Z^T` |
|---|---|
| rotate the null-space basis | **5.6e-16** -- rounding, nothing else |
| render one exactly-zero singular value as `1.65e-17` | **96 per cent** |

## What is not invariant is the rendering of the near-zero singular values

The second row is the real hazard, and it is not a basis question. Singular values are
mathematically unique, so any conforming implementation agrees on them -- but agreement at the
rounding level is not enough here, because `1 / (sqrt(s) + 1e-10)` at `s ~ 0` is not a continuous
function of `s` in any useful sense: `0` gives `1e10`, `1.65e-17` gives `2.4e8`, a factor of 40, and
the Gram matrix moves by a factor of 1600 in that column's contribution.

The golden shows the reference doing **both within one decomposition**: `1.65e-17` for the fourth
singular value of the two-step series and exact `0.0` for the fifth and sixth. Where its
bidiagonalisation stops producing digits and starts producing zeros is an artefact of the algorithm,
not a property of the matrix, and a different implementation will draw that line somewhere else.

## So the decision #230 asked for

**Byte identity through `U` is not required and not meaningful** -- proven above, not estimated.

**Byte identity through the singular values is required if `Z` is required**, and reproducing them
means reproducing commons-math's rounding at the rank deficiency, which is what a transcription
buys and a crate does not.

**Whether the changepoints survive anyway is an empirical question, not an argument.** They are
integers chosen by local minima and a backward merge, so a 96 per cent change in the Gram matrix may
still leave them where they were, and may not. That measurement needs a second implementation run on
the same matrices and has to be done on real pileups rather than on the four synthetic series, since
the flat and step series are exactly the rank-deficient cases.

Until it is done, the honest position for `CalculateContamination` is the one the golden already
pins: the decomposition is measured, the changepoints are measured, and the tool waits under G2.3
rather than being ported with a claim nobody has checked.
