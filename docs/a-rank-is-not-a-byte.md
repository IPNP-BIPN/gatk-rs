# A rank is not a byte

`CreateReadCountPanelOfNormals` writes a panel whose last component is a basis of eigensamples,
computed by a singular value decomposition. The dump reported how many there were. The number moved:

    committed golden (produced on CI)   eigensamples=8
    a later CI run, same image          eigensamples=6
    this laptop, twice, same container   eigensamples=7

The three numbers are three answers to the same question about the same fixture, and none of them
is wrong. The count is the decomposition's RANK, which is the number of singular values a solver
judged to be non-zero, and that judgement is a comparison against a tolerance made on values that
are themselves the output of an iterative, distributed solver. Nothing in the tool's input decides
it; the arithmetic path does.

This is the third hazard of this kind the programme has hit, and they are worth naming together:

  * **log4j's clock**, which is time and was masked;
  * **BWA's index image**, whose differing bytes are in-process POINTERS
    (`docs/pointers-that-reach-the-output.md`);
  * and a **rank**, which is a decision about a number rather than the number itself.

The first two are things the reference writes down that are not about the data. This one is
different and more insidious: it IS about the data, it is stable enough to look reproducible, and it
is what a reader would naturally record. The suite carried it as a golden line for weeks and only a
re-run on a differently loaded runner showed it moving.

## What the suite reports instead

Two facts that do not move:

  * `eigensamples-at-most`, the number of samples the panel was built from, which is the cap the
    rank is taken under and is a count of the input;
  * `eigensamples-positive`, whether the basis is empty at all, which is the DECISION the tool
    makes when every singular value is zero and is what the refusal path downstream depends on.

The singular values themselves were never in the golden, and their count is the same rank, so what
is reported of them is whether the panel hands them over: a panel with no eigensamples refuses,
which is a decision rather than a number.

## The rule this leaves

A number a golden holds has to be a count of something the input decides, or a value the reference
computes by a path that has no tolerance in it. A rank, a cluster count, an iteration count and a
convergence flag are none of those. Where one of them is the interesting fact, measure the BOUND it
is taken under and the branch it selects, and say in the dump why the number itself is absent.
