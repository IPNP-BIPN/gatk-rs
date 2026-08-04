# Commenting gatk-rs to the standard: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring all 70 Rust source files in `gatk-rs` to the commenting standard in `docs/COMMENTING.md`, so a reader who knows Java and biology but not Rust can check the port against the reference without learning Rust first.

**Architecture:** No behaviour changes. Every task edits comments only, then proves it changed nothing by running the existing test suite unchanged, and records the file's new comment density in `tools/audit/commented.txt` so the CI ratchet protects it. One task per coherent module group; each task ends in a commit and a green ratchet, so tasks can be reviewed and merged independently and in any order.

**Tech Stack:** Rust 1.97.1 (pinned in `rust-toolchain.toml`), `cargo test` / `cargo clippy` / `cargo fmt`, Python 3 for `tools/audit/comment_density.py`, GitHub Actions for CI.

## Global Constraints

- Rust toolchain is pinned to **1.97.1**. Do not change `rust-toolchain.toml` or any dependency version.
- **No behaviour changes in any task of this plan.** Comments, doc comments and the density list only. If a task tempts you to fix a bug, stop and open a separate issue; a comment change that also changes code cannot be reviewed as either.
- **Never use em-dashes** (—) in any output: code comments, commit messages, PR bodies, documentation. Use commas, colons, parentheses or separate sentences.
- Commit messages end with the trailer `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`. PR bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- The pinned reference clone at `gatk/` is **read-only**. Read it, never edit it. Do not run `git` commands from inside it: `cd /Users/benjamin/GATK` first, every time.
- Every comment must be defensible from the reference. If the reason is a line of Java, quote the line. If the reason is a measurement, name the golden row. Do not write a comment whose justification is a guess.
- One PR per task group, merged only on green CI.
- `docs/COMMENTING.md` is the standard this plan implements. Read it before Task 1 and keep it open.

---

## File Structure

No files are created except the plan's own tooling change in Task 1. The work is spread over the existing tree:

| Path | Responsibility | Task |
|---|---|---|
| `tools/audit/comment_density.py` | measure density, ratchet against regression | 1 |
| `tools/audit/commented.txt` | the list of files at the standard, with floors | every task |
| `crates/gatk-engine/src/histogram.rs` | GATK's run-length histogram; **partially done**, at 0.609 | 2 |
| `crates/gatk-engine/src/{java_hash,java_random,well19937c,permutation,fisher_exact,mann_whitney}.rs` | numeric and Java-semantics primitives | 3 |
| `crates/gatk-engine/src/{allele_likelihoods,allele_list}.rs` | the likelihood matrix and its two axes | 4 |
| `crates/gatk-engine/src/{interval,interval_args,locus_shards,features,feature_intervals,variant_source,variant_getters}.rs` | intervals and feature sources | 5 |
| `crates/gatk-engine/src/{read,read_utils,read_group,reads,read_states,read_pileup,pileup,alignment_state,downsampling}.rs` | reads and pileups | 6 |
| `crates/gatk-engine/src/{cigar_builder,cigar_utils,clipping}.rs` | cigar arithmetic and clipping | 7 |
| `crates/gatk-engine/src/{context,context_iterator,locus_iterator,reference,assembly_region,assembly_region_iterator,assembly_region_walker,activity_profile}.rs` | contexts, regions, traversal | 8 |
| `crates/gatk-engine/src/{jexl,lib}.rs` | the JEXL expression engine and the crate root | 9 |
| `crates/gatk-annotation/src/{lib,info_annotation,chromosome_counts,coverage,raw_gt_count,sample_list,original_alignment}.rs` | the annotation interface and the counting family | 10 |
| `crates/gatk-annotation/src/{rank_sum,per_allele,strand_bias,read_grouping}.rs` | rank sums, medians, strand bias | 11 |
| `crates/gatk-annotation/src/allele_specific_{rank_sum,strand_bias,site_statistics}.rs` | the allele-specific family | 12 |
| `crates/gatk-annotation/src/{flow,pedigree,fragment_counts,site_statistics,tandem_repeat,heterozygosity,mapping_quality,depth_per_allele}.rs` | the rest of the annotations | 13 |
| `crates/gatk-readfilter/src/{lib,counting}.rs`, `crates/gatk-tools/src/*.rs`, `crates/gatk-corpus/src/lib.rs` | read filters, walkers, the corpus reader | 14 |
| `ROADMAP.md` | the running count | 15 |

Test files under `crates/*/tests/` are **out of scope for this plan**. They are read by whoever is checking a golden, and their density is low by design: a fixture table explains itself. A separate plan may take them later.

---

### Task 1: Close the ratchet's back door

`comment_density.py --record` currently rewrites every listed floor from the current tree. Someone who deletes the comments in a finished file and runs `--record` lowers its floor silently, which defeats the whole guard. Close that before adding 67 more files to the list.

**Files:**
- Modify: `tools/audit/comment_density.py`
- Test: `tools/audit/test_comment_density.py` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `comment_density.py --record` refuses to lower an existing floor unless `--allow-lower` is passed. Every later task calls `--record` and relies on that refusal. Exit code 1 on refusal, 0 on success.

- [ ] **Step 1: Write the failing test**

Create `tools/audit/test_comment_density.py`:

```python
"""The ratchet must not have a back door. `--record` may raise a floor and may add a file, but
lowering one is a mistake often enough that it has to be asked for explicitly."""

import pathlib
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
SCRIPT = HERE / "comment_density.py"


def run(args, cwd):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def make_tree(root, comment_lines, code_lines):
    """A throwaway repository holding one Rust file with a known comment ratio."""
    crate = root / "crates" / "demo" / "src"
    crate.mkdir(parents=True)
    body = "\n".join(["// a comment"] * comment_lines + ["let x = 1;"] * code_lines)
    (crate / "demo.rs").write_text(body + "\n", encoding="utf-8")
    return crate / "demo.rs"


def test_record_refuses_to_lower_a_floor(tmp_path):
    make_tree(tmp_path, comment_lines=10, code_lines=10)
    listed = tmp_path / "tools" / "audit" / "commented.txt"
    listed.parent.mkdir(parents=True)
    listed.write_text("crates/demo/src/demo.rs 1.000\n", encoding="utf-8")

    # Gut the comments: the ratio falls from 1.0 to 0.0.
    (tmp_path / "crates" / "demo" / "src" / "demo.rs").write_text(
        "let x = 1;\n" * 10, encoding="utf-8"
    )

    result = run(["--record", "--root", str(tmp_path)], cwd=tmp_path)
    assert result.returncode == 1, result.stdout + result.stderr
    assert "would lower" in result.stdout
    # The file on disk is unchanged, so the floor still protects the file.
    assert "1.000" in listed.read_text(encoding="utf-8")


def test_record_allows_lowering_when_asked(tmp_path):
    make_tree(tmp_path, comment_lines=10, code_lines=10)
    listed = tmp_path / "tools" / "audit" / "commented.txt"
    listed.parent.mkdir(parents=True)
    listed.write_text("crates/demo/src/demo.rs 1.000\n", encoding="utf-8")
    (tmp_path / "crates" / "demo" / "src" / "demo.rs").write_text(
        "let x = 1;\n" * 10, encoding="utf-8"
    )

    result = run(["--record", "--allow-lower", "--root", str(tmp_path)], cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "0.000" in listed.read_text(encoding="utf-8")


def test_record_raises_a_floor_without_asking(tmp_path):
    make_tree(tmp_path, comment_lines=10, code_lines=10)
    listed = tmp_path / "tools" / "audit" / "commented.txt"
    listed.parent.mkdir(parents=True)
    listed.write_text("crates/demo/src/demo.rs 0.500\n", encoding="utf-8")

    result = run(["--record", "--root", str(tmp_path)], cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "1.000" in listed.read_text(encoding="utf-8")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /Users/benjamin/GATK && python3 -m pytest tools/audit/test_comment_density.py -v`

Expected: FAIL. The first two tests fail on the unrecognised `--root` argument; `comment_density.py` currently hardcodes its root and has no `--allow-lower`.

If `pytest` is not installed, run `python3 -m pip install --user pytest` first. Do not add pytest to any Cargo manifest.

- [ ] **Step 3: Add `--root` and `--allow-lower`**

In `tools/audit/comment_density.py`, replace the module-level constants and the `--record` branch.

Replace:

```python
ROOT = pathlib.Path(__file__).resolve().parents[2]
CRATES = ROOT / "crates"
LISTED = pathlib.Path(__file__).resolve().parent / "commented.txt"
```

with:

```python
# The repository root. Overridable so the tests can point the script at a throwaway tree; every
# real invocation leaves it alone.
DEFAULT_ROOT = pathlib.Path(__file__).resolve().parents[2]
```

Change `sources()` to take the root:

```python
def sources(root: pathlib.Path) -> list[pathlib.Path]:
    """Every Rust source under `crates/`, excluding build output, sorted for a stable report."""
    return sorted(
        p for p in (root / "crates").rglob("*.rs") if "target" not in p.parts
    )
```

Change `load_floors()` to take the list path:

```python
def load_floors(listed: pathlib.Path) -> dict[str, float]:
    """The recorded floor for each file that has been brought to the standard.

    Format is one `<path> <ratio>` pair per line, `#` starting a comment. Paths are relative to the
    repository root so the file reads as a checklist.
    """
    floors: dict[str, float] = {}
    if not listed.exists():
        return floors
    for raw in listed.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        name, _, value = line.partition(" ")
        floors[name.strip()] = float(value.strip())
    return floors
```

In `main()`, add the two arguments and resolve the paths:

```python
    parser.add_argument("--check", action="store_true", help="fail on regression")
    parser.add_argument("--record", action="store_true", help="rewrite the floors")
    parser.add_argument(
        "--allow-lower",
        action="store_true",
        help="with --record, permit a floor to drop; without it, a drop is refused",
    )
    parser.add_argument(
        "--root",
        default=str(DEFAULT_ROOT),
        help="repository root; only the tests pass this",
    )
    args = parser.parse_args()

    root = pathlib.Path(args.root).resolve()
    listed = root / "tools" / "audit" / "commented.txt"
```

Then replace every use of `ROOT`, `CRATES` and `LISTED` in `main()` with `root`, `root / "crates"` and `listed`, and pass them into `sources` and `load_floors`.

Finally, guard the `--record` branch. Replace the loop that builds `lines` with:

```python
    if args.record:
        # A floor may rise freely: a file gained explanations. A floor may only fall on request,
        # because the usual reason for a fall is that someone deleted comments and did not notice.
        lowered = [
            (name, floors[name], measured[name][2])
            for name in sorted(floors)
            if name in measured and measured[name][2] < floors[name] - 0.02
        ]
        if lowered and not args.allow_lower:
            for name, floor, current in lowered:
                print(f"{name}: recording would lower the floor {floor:.3f} to {current:.3f}")
            print()
            print("if that is deliberate, pass --allow-lower and say why in the commit message")
            return 1

        lines = [
            "# Files brought to the standard in docs/COMMENTING.md, with the comment-to-code",
            "# ratio each had when it was done. tools/audit/comment_density.py --check fails if",
            "# one of them drops below its recorded floor, so the work is a ratchet: nothing",
            "# already explained can quietly come undone.",
            "#",
            "# Regenerate with: python3 tools/audit/comment_density.py --record",
            "",
        ]
        for name in sorted(floors):
            if name not in measured:
                print(f"listed file no longer exists: {name}", file=sys.stderr)
                return 1
            lines.append(f"{name} {measured[name][2]:.3f}")
        listed.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"recorded {len(floors)} floors")
        return 0
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /Users/benjamin/GATK && python3 -m pytest tools/audit/test_comment_density.py -v`

Expected: PASS, 3 passed.

- [ ] **Step 5: Verify the real tree still checks clean**

Run:

```bash
cd /Users/benjamin/GATK
python3 tools/audit/comment_density.py --check
python3 tools/audit/comment_density.py | head -5
```

Expected: `3 files at or above their recorded comment density`, then the summary. If `--check` fails, something else regressed; stop and investigate before continuing.

- [ ] **Step 6: Commit**

```bash
cd /Users/benjamin/GATK
git checkout -b docs/commenting-tranche-1
git add tools/audit/comment_density.py tools/audit/test_comment_density.py
git commit -m "$(cat <<'EOF'
Close the ratchet's back door

--record rewrote every listed floor from the current tree, so deleting the
comments in a finished file and re-recording lowered its floor silently.
That defeats the guard, which exists precisely to catch that.

A floor may now rise freely, because a file that gained explanations needs
no permission. A floor may only fall when --allow-lower is passed, so the
fall is a decision someone made and can be asked about in review.

--root is added for the tests alone, so they can run against a throwaway
tree rather than the repository.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Finish `histogram.rs`

It is on the list at 0.609 because it was patched section by section rather than worked through. Its two neighbours are at 1.93 and 1.15. Bring it to parity so the list means one thing.

**Files:**
- Modify: `crates/gatk-engine/src/histogram.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new. Behaviour is unchanged and no signature moves.

- [ ] **Step 1: Read the reference**

Read `gatk/src/main/java/org/broadinstitute/hellbender/utils/Histogram.java` and `CompressedDataList.java` in full. They are 160 and 90 lines. You are looking for the parts the current comments do not yet explain: `get`, `add(Histogram)`, `with_bin_size`'s precision derivation, `Display for Histogram`, and `format_fixed`.

- [ ] **Step 2: Comment the uncommented items**

Every `pub fn` and every non-trivial block in `histogram.rs` gets the three questions. The pattern to follow is the one already in the file for `binned_value` and `median`. Specific things that must be explained because they are not guessable:

- `with_bin_size`: the precision is `Math.round(-Math.log10(binSize))`, so a bin size that is not a power of a tenth still gets a whole number of decimals, and a bin size above one gives a negative exponent that `String.format` would reject. No caller in the annotations does either.
- `add_histogram`: the bin sizes must match exactly, compared as floats. Two histograms built with `0.1` always match; one built with `1.0/10.0` would not be guaranteed to.
- `Display for Histogram`: an empty histogram renders as the four characters `NaN`, because the reference returns `Double.toString(Double.NaN)`. Say where that reaches a record: `AS_RAW_BaseQRankSum=|NaN` in the allele-specific rank-sum golden.
- `format_fixed`: it is a second copy of the half-up decimal rule that `gatk-annotation` also has. Explain why the duplication is deliberate: this crate sits below that one, and a histogram's rendering must not depend on an annotation crate.
- `get`: returns absent rather than zero for a bin nothing landed in, which is a different fact from a bin with a count of zero.

- [ ] **Step 3: Verify nothing changed**

Run:

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine --lib
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

Expected: tests pass with the same count as before (34 at the time of writing), clippy silent, no rustdoc warnings. A rustdoc warning here almost always means a broken `[link]` in a doc comment you just wrote.

- [ ] **Step 4: Record the new floor and check**

```bash
cd /Users/benjamin/GATK
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
grep histogram tools/audit/commented.txt
```

Expected: `--record` succeeds (the floor rose, which needs no permission), `--check` passes, and the printed ratio is at least 1.0. If it is below 1.0, the file is not yet at parity: go back to Step 2.

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src/histogram.rs tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Finish histogram.rs

It went on the list at 0.609 because it was patched section by section
rather than worked through, while its two neighbours are at 1.93 and 1.15.
A list whose entries mean different things is not a list.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: The engine's remaining numeric primitives

**Files:**
- Modify: `crates/gatk-engine/src/java_hash.rs`
- Modify: `crates/gatk-engine/src/java_random.rs`
- Modify: `crates/gatk-engine/src/well19937c.rs`
- Modify: `crates/gatk-engine/src/permutation.rs`
- Modify: `crates/gatk-engine/src/fisher_exact.rs`
- Modify: `crates/gatk-engine/src/mann_whitney.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read each file's reference**

Each file's module header (`//!` at the top) names the Java class it ports. Open that class in `gatk/src/main/java/...` before commenting the file. Do not comment a file whose reference you have not read: a comment that guesses at the why is worse than no comment, and this standard exists to stop exactly that.

Known hooks, which must appear in the comments because they are the reason each file is not the obvious implementation:

- `java_hash.rs`: `hash_map_order` reproduces Java's `HashMap` iteration, which is observable wherever a port takes a key set as an array. `string_hash_code` is exact by specification, which is what makes the reproduction possible. The treeify threshold is refused rather than guessed.
- `java_random.rs` and `well19937c.rs`: a seeded generator's stream is part of the output. `well19937c.rs` already carries a refusal note; make sure a reader sees what is refused and why before they see any code.
- `permutation.rs`: the identity permutation is a distinct case in the reference and skips work; say what that means for a caller.
- `fisher_exact.rs`: `REL_ERR = 1.0 - 10e-7`, which is `1 - 1e-6` and not `1 - 1e-7`, and R's own algorithm is the model rather than a textbook Fisher test. `pow10` names decision 0007.
- `mann_whitney.rs`: the ranks are `f32`, not `f64`. That single fact changes reported Z scores and is the most important comment in the crate. Also: the continuity correction is dropped when there are no ties, and the alternate series is the **first** argument so swapping the two flips the sign of every rank-sum annotation.

- [ ] **Step 2: Comment the six files**

Follow `crates/gatk-engine/src/math_utils.rs` as the worked example. Its shape is: module header (already present, leave it), then per item the three questions, then inline comments where the Rust idiom is not guessable from Java and where the reference does something surprising.

Explain each Rust idiom at most once per file. The table in `docs/COMMENTING.md` lists the ones worth explaining.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine --lib
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

Expected: same test count as before the task, clippy silent, no doc warnings.

- [ ] **Step 4: Run the conformance suites that cover these files**

```bash
cd /Users/benjamin/GATK
cargo test -p gatk-engine --test mann_whitney --test fisher_exact 2>/dev/null || cargo test -p gatk-engine
```

Expected: pass. These files back oracle-backed suites; if a suite fails, you changed behaviour, which this plan forbids. Revert and try again.

- [ ] **Step 5: Add the six files to the list, record, check**

```bash
cd /Users/benjamin/GATK
for f in java_hash java_random well19937c permutation fisher_exact mann_whitney; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

Expected: `recorded 10 floors`, then `10 files at or above their recorded comment density`.

- [ ] **Step 6: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the engine's numeric primitives

Six files whose implementations are not the obvious ones, and whose
comments now say why: the HashMap iteration order that is observable
through a key set, the seeded generators whose stream is part of the
output, R's Fisher algorithm rather than a textbook one, and the
Mann-Whitney ranks that are 32-bit floats and change every reported Z
score because of it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: The likelihood matrix and its axes

The widest-reaching type in the crate. Every annotation reads it, so this is the file a new reader hits first and understands least.

**Files:**
- Modify: `crates/gatk-engine/src/allele_likelihoods.rs`
- Modify: `crates/gatk-engine/src/allele_list.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the reference**

`gatk/src/main/java/org/broadinstitute/hellbender/utils/genotyper/AlleleLikelihoods.java`, `IndexedAlleleList.java`, `IndexedSampleList.java`. `AlleleLikelihoods` is long; read `bestAllelesBreakingTies`, `searchBestAllele`, `marginalize` and `groupEvidence` closely and skim the rest.

- [ ] **Step 2: Comment both files**

The facts that must be explained, because each has already caught the port at least once:

- `search_best_allele` breaks a tie by **keeping the first index**, so the allele order decides which allele a read supporting two equally is attributed to. That is how `HashMap` order reaches `AD`.
- `is_informative` is an **absolute** difference against a threshold, so summing two reads' likelihoods into a fragment makes the threshold easier to clear.
- `marginalize` takes the **maximum** over the old alleles, not the sum, and preserves `is_natural_log`.
- `group_by_fragment` sums the log likelihoods of a group, and the fragment it builds may hold fewer reads than the sum covered. Point at the golden row that shows it.
- `value(sample, allele, evidence)` is indexed `[sample][allele][evidence]`, which is not the order the constructor's argument names suggest at a glance.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo test -p gatk-annotation
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

Expected: all pass. `gatk-annotation`'s suites are the real check here: they exercise the matrix through every annotation.

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in allele_likelihoods allele_list; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the likelihood matrix and its axes

The widest-reaching type in the crate, and the one a reader hits first and
understands least. The comments now carry the four facts that have each
caught this port at least once: the tie broken by keeping the first index,
which is how HashMap order reaches AD; the informativeness threshold being
an absolute difference, which is why summing a pair makes it easier to
clear; marginalize taking a maximum and not a sum; and the index order of
the value table.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Intervals and feature sources

**Files:**
- Modify: `crates/gatk-engine/src/interval.rs`
- Modify: `crates/gatk-engine/src/interval_args.rs`
- Modify: `crates/gatk-engine/src/locus_shards.rs`
- Modify: `crates/gatk-engine/src/features.rs`
- Modify: `crates/gatk-engine/src/feature_intervals.rs`
- Modify: `crates/gatk-engine/src/variant_source.rs`
- Modify: `crates/gatk-engine/src/variant_getters.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references named in each module header**

`SimpleInterval`, `IntervalArgumentCollection`, `IntervalUtils`, `ShardBoundary`, `FeatureDataSource`, `VariantContextGetters`. Read `IntervalUtils.parseIntervalArguments` and its error paths carefully: `interval_args.rs` is mostly a catalogue of the reference's exceptions and each variant needs its Java counterpart named.

- [ ] **Step 2: Comment the seven files**

Facts that must appear:

- an interval parse is ambiguous when a contig name contains a colon, and the reference resolves it by trying the whole string as a contig first;
- commas are stripped from positions, and a trailing `+` runs to the end of the contig;
- adjacency merges only under the `ALL` rule, not under the default;
- `feature_intervals.rs` dispatches BED, IntervalList and VCF, and only the VCF codec decides by reading the file's bytes rather than its name;
- each `interval_args.rs` error variant names the exact reference exception class it stands for, since that class is what a conformance dump records.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in interval interval_args locus_shards features feature_intervals variant_source variant_getters; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the intervals and the feature sources

An interval is what -L means, and a tool given the wrong one reads the
wrong data and every number it produces is wrong in a way no downstream
comparison can attribute. The comments now say where the parse is
ambiguous, how the reference resolves it, and which exception class each
refusal stands for, since that class is what a conformance dump records.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Reads and pileups

**Files:**
- Modify: `crates/gatk-engine/src/read.rs`
- Modify: `crates/gatk-engine/src/read_utils.rs`
- Modify: `crates/gatk-engine/src/read_group.rs`
- Modify: `crates/gatk-engine/src/reads.rs`
- Modify: `crates/gatk-engine/src/read_states.rs`
- Modify: `crates/gatk-engine/src/read_pileup.rs`
- Modify: `crates/gatk-engine/src/pileup.rs`
- Modify: `crates/gatk-engine/src/alignment_state.rs`
- Modify: `crates/gatk-engine/src/downsampling.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

`GATKRead`, `ReadUtils`, `ReadsDataSource`, `AlignmentStateMachine`, `ReadPileup`, `PileupElement`, `LocusIteratorByState`, `ReservoirDownsampler` / `PositionalDownsampler`.

- [ ] **Step 2: Comment the nine files**

Facts that must appear:

- the downsampler's reservoir is seeded, so which reads survive is reproducible and is part of the output; name the generator;
- `reads.rs` merges abutting intervals so a read spanning a boundary is returned once, and a zero end swallows what follows on the same contig;
- the alignment state machine's position advances differently for each cigar operator, and the deletion case is the one that surprises;
- `read_utils.rs` has soft-clip-aware starts and ends that are not the alignment start and end, and which one a caller wants is the usual bug.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in read read_utils read_group reads read_states read_pileup pileup alignment_state downsampling; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the reads and the pileups

Nine files where the surprises are about which coordinate is meant. The
comments now distinguish the soft-clip-aware start from the alignment
start at every point one is chosen, say why the downsampler's reservoir is
seeded and therefore part of the output, and record that abutting
intervals merge so a read spanning a boundary is returned once.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Cigar arithmetic and clipping

**Files:**
- Modify: `crates/gatk-engine/src/cigar_builder.rs`
- Modify: `crates/gatk-engine/src/cigar_utils.rs`
- Modify: `crates/gatk-engine/src/clipping.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

`CigarBuilder`, `CigarUtils`, `ReadClipper`, `ClippingOp`. `cigar_builder.rs` is at 0.145, the lowest in the crate, and is the file most in need of this task.

- [ ] **Step 2: Comment the three files**

Facts that must appear, all already witnessed by the crate's own unit tests:

- consecutive identical operators merge, so a builder's output is not the sequence of pushes;
- a deletion after an insertion moves **before** it, which is a normalisation the reference performs and a naive builder would not;
- deletions at either end are removed and their lengths counted, and the count is what a caller uses to fix the alignment start;
- a completely soft-clipped cigar is refused, and so is a soft clip in the middle.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in cigar_builder cigar_utils clipping; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the cigar arithmetic and the clipping

cigar_builder.rs was the least explained file in the crate at 0.145, and
it is the one whose output least resembles its input: consecutive
operators merge, a deletion after an insertion moves before it, and
deletions at either end are removed and counted so the caller can fix the
alignment start. All three are normalisations the reference performs and a
naive builder would not.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Contexts, regions and traversal

**Files:**
- Modify: `crates/gatk-engine/src/context.rs`
- Modify: `crates/gatk-engine/src/context_iterator.rs`
- Modify: `crates/gatk-engine/src/locus_iterator.rs`
- Modify: `crates/gatk-engine/src/reference.rs`
- Modify: `crates/gatk-engine/src/assembly_region.rs`
- Modify: `crates/gatk-engine/src/assembly_region_iterator.rs`
- Modify: `crates/gatk-engine/src/assembly_region_walker.rs`
- Modify: `crates/gatk-engine/src/activity_profile.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

`ReferenceContext`, `FeatureContext`, `AlignmentContext`, `ReferenceDataSource`, `AssemblyRegion`, `AssemblyRegionIterator`, `AssemblyRegionWalker`, `ActivityProfile` and `BandPassActivityProfile`.

- [ ] **Step 2: Comment the eight files**

Facts that must appear:

- a reference context carries a **window** that may be wider than the interval, and an annotation reading it must say which it wants; `REF_BASES` pads with `N` on the right only and can therefore be off-centre;
- soft masking and ambiguity codes do not survive a reference query, which the crate's own test asserts;
- the activity profile's band pass filter has a fixed kernel whose width is derived from a sigma, and the derivation rounds;
- an assembly region's extended span is not its span, and which one a caller wants is the usual bug here too.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in context context_iterator locus_iterator reference assembly_region assembly_region_iterator assembly_region_walker activity_profile; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the contexts, the regions and the traversal

Eight files where the recurring question is which span is meant: an
interval or the window around it, a region or its extended form. The
comments now say which at every point one is chosen, and record that soft
masking and ambiguity codes do not survive a reference query.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: JEXL and the crate root

**Files:**
- Modify: `crates/gatk-engine/src/jexl.rs`
- Modify: `crates/gatk-engine/src/lib.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the reference**

`VariantContextUtils.match`, `JexlEngine` usage in GATK, and the Apache Commons JEXL 2 grammar for the subset GATK exercises. `jexl.rs` is 720 lines at 0.200, the second-least explained file in the crate.

- [ ] **Step 2: Comment both files**

`jexl.rs` is an expression parser and evaluator, so the comments carry more weight than usual: a reader cannot check a parser against a grammar without knowing which production each function is. Name the production for each parsing function, and say what the reference does with an undefined attribute, since that is the case every filter expression hits.

`lib.rs` is the crate root: its module list is the crate's table of contents, and each `pub mod` deserves a one-line description of what lives there.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-engine
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-engine 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in jexl lib; do
  echo "crates/gatk-engine/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit and open the PR for Tasks 1 to 9**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-engine/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the JEXL engine and the crate root

A reader cannot check a parser against a grammar without knowing which
production each function is, so each parsing function now names its
production, and the undefined-attribute case is spelled out because every
filter expression hits it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
git push -u origin docs/commenting-tranche-1
gh workflow run "CI" --ref docs/commenting-tranche-1
```

Wait for CI to go green, then open the PR. Body must state: no behaviour changed, the test counts before and after are identical, and the list grew from 3 files to 40.

---

### Task 10: The annotation interface and the counting family

**Files:**
- Modify: `crates/gatk-annotation/src/lib.rs`
- Modify: `crates/gatk-annotation/src/info_annotation.rs`
- Modify: `crates/gatk-annotation/src/chromosome_counts.rs`
- Modify: `crates/gatk-annotation/src/coverage.rs`
- Modify: `crates/gatk-annotation/src/raw_gt_count.rs`
- Modify: `crates/gatk-annotation/src/sample_list.rs`
- Modify: `crates/gatk-annotation/src/original_alignment.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Start a new branch**

```bash
cd /Users/benjamin/GATK
git checkout main && git pull
git checkout -b docs/commenting-tranche-2
```

- [ ] **Step 2: Read the references**

`InfoFieldAnnotation`, `ChromosomeCounts`, `Coverage`, `RawGtCount`, `SampleList`, `OriginalAlignment` in `gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/annotator/`.

- [ ] **Step 3: Comment the seven files**

The three facts the crate root already states must be carried down into `info_annotation.rs` item by item, because they are what makes `AnnotationValue` an enum rather than a string:

- the Java type of the value put in the map is observable: `Coverage` puts a `String` from `String.format("%d", depth)`, `CountNs` puts a boxed `Long`, `ChromosomeCounts` puts an `Integer` or an `ArrayList` depending on the alternate count, and the encoder renders each differently in the edge cases;
- an annotation with nothing to say returns an **empty map**, not a zero, so the key is absent from the record;
- `getKeyNames()` is the declaration order and nothing else, because the encoder sorts.

`to_java_string` returns `None` for a `Double` and must say why: `Double.toString` is its own algorithm, not a format string, and producing a plausible rendering would be inventing a golden.

- [ ] **Step 4: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-annotation
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-annotation 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 5: Record and check**

```bash
cd /Users/benjamin/GATK
for f in lib info_annotation chromosome_counts coverage raw_gt_count sample_list original_alignment; do
  echo "crates/gatk-annotation/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 6: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-annotation/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the annotation interface and the counting family

AnnotationValue is an enum rather than a string because the Java type of
the value put in the map is observable, and the comments now say so at
each variant. Two more facts carried down from the crate root into the
items that implement them: an annotation with nothing to say returns an
empty map and not a zero, and getKeyNames is the declaration order and
nothing else, because the encoder sorts.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: Rank sums, medians and strand bias

**Files:**
- Modify: `crates/gatk-annotation/src/rank_sum.rs`
- Modify: `crates/gatk-annotation/src/per_allele.rs`
- Modify: `crates/gatk-annotation/src/strand_bias.rs`
- Modify: `crates/gatk-annotation/src/read_grouping.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

`RankSumTest` and its four subclasses, `PerAlleleAnnotation` and its four, `StrandBiasTest`, `FisherStrand`, `StrandOddsRatio`, `StrandBiasBySample`, `UniqueAltReadCount`, `BaseQualityHistogram`, `ReferenceBases`.

- [ ] **Step 2: Comment the four files**

Facts that must appear:

- `format_decimals` rounds **half up on the decimal expansion**, which is Java's rule and not Rust's half-to-even. `0.0625` prints `0.063` in Java and `0.062` in Rust, and a Z score can be exactly that;
- `INVALID_ELEMENT_FROM_READ` is negative infinity, a **value** and not an absence, and reads carrying it are dropped at a different point from reads with no value at all;
- `FisherStrand` uses a minimum count of 2 and `StrandOddsRatio` uses 0, so the two disagree on a site whose only sample has one or two reads;
- the genotype `SB` field wins over the matrix: if any genotype carries it, the likelihoods are never consulted;
- `AS_UNIQ_ALT_READ_COUNT` counts distinct `(start, fragmentLength)` pairs, so a hundred PCR duplicates count once and two genuinely distinct fragments sharing both count once too.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-annotation
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-annotation 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in rank_sum per_allele strand_bias read_grouping; do
  echo "crates/gatk-annotation/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-annotation/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the rank sums, the medians and the strand bias

Four files whose numbers differ from the obvious implementation in ways a
type signature cannot show: a formatter that rounds half up on the decimal
expansion where Rust rounds half to even, an invalid element that is a
value and not an absence, and two strand-bias annotations with different
minimum counts that therefore disagree on a shallow site.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: The allele-specific family

**Files:**
- Modify: `crates/gatk-annotation/src/allele_specific_rank_sum.rs`
- Modify: `crates/gatk-annotation/src/allele_specific_strand_bias.rs`
- Modify: `crates/gatk-annotation/src/allele_specific_site_statistics.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

The `allelespecific/` package: `AS_RankSumTest`, `AS_StrandBiasTest`, `StrandBiasUtils`, `AS_QualByDepth`, `AS_RMSMappingQuality`, `AS_InbreedingCoeff`, `AlleleSpecificAnnotationData`.

- [ ] **Step 2: Comment the three files**

These three module headers already carry the findings. This task's work is to push them **down into the items**, so a reader who jumps to a function does not have to scroll up. Every one of these must appear at the item that implements it:

- the direct `annotate()` path of an allele-specific annotation is not allele-specific: it pools every alternate into one series;
- the rank-sum raw string starts with its delimiter because the reference's slot is skipped but kept, while the strand-bias raw string puts the delimiter between entries and includes the reference. Two families, two conventions, one delimiter;
- an entry that is present but empty renders as nothing, and is not the same as an absent entry, which renders as `0,0`;
- `AS_SOR` computes a value for the reference against itself and never prints it, while `AS_FS` filters the reference out first;
- `AS_QUAL` is read with `getAttributeAsList`, which does **not** split a comma-separated string;
- `AS_RMSMappingQuality` collects into a container nothing pre-populates, so an allele with no read is skipped entirely, while the parsed form is pre-populated with nulls and writes the missing value for the same state.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-annotation
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-annotation 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in allele_specific_rank_sum allele_specific_strand_bias allele_specific_site_statistics; do
  echo "crates/gatk-annotation/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-annotation/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Push the allele-specific findings down into the items

The three module headers already carried them. A reader who jumps to a
function had to scroll up to learn that the direct path is not
allele-specific, that the two families use the same delimiter with
opposite conventions, and that a present-but-empty entry is not an absent
one. Each now sits at the item that implements it.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 13: The rest of the annotations

**Files:**
- Modify: `crates/gatk-annotation/src/flow.rs`
- Modify: `crates/gatk-annotation/src/pedigree.rs`
- Modify: `crates/gatk-annotation/src/fragment_counts.rs`
- Modify: `crates/gatk-annotation/src/site_statistics.rs`
- Modify: `crates/gatk-annotation/src/tandem_repeat.rs`
- Modify: `crates/gatk-annotation/src/heterozygosity.rs`
- Modify: `crates/gatk-annotation/src/mapping_quality.rs`
- Modify: `crates/gatk-annotation/src/depth_per_allele.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Read the references**

The flow package, `PossibleDeNovo`, `TransmittedSingleton`, `MendelianViolation`, `OrientationBiasReadCounts`, `FragmentDepthPerAlleleBySample`, `QualByDepth`, `GenotypeSummaries`, `TandemRepeat`, `GATKVariantContextUtils`'s repeat helpers, `ExcessHet`, `InbreedingCoeff`, `HeterozygosityCalculator`, `GenotypeUtils`, `RMSMappingQuality`, `MappingQualityZero`, `DepthPerAlleleBySample`.

- [ ] **Step 2: Comment the eight files**

As in Task 12, the headers carry the findings and this task pushes them into the items. The ones that must land at an item:

- `QualByDepth` above 35 is randomised and this port refuses the branch; the refusal must be explained where the error variant is defined, not only in the header;
- `TransmittedSingleton` reads the **child's** depth three times under three names; the comment goes on the three lines;
- an absent `DP` is minus one, not zero, and the default threshold is zero, so no-depth and zero-depth take different branches;
- `findRepeatedSubstring` cannot see a partial trailing repeat because `Arrays.copyOfRange` pads with zero bytes, and returns a unit of one zero byte on an empty input;
- the heterozygosity calculator adds the hom-ref mass to the reference's count **inside** the loop over alternates, so a triallelic site counts it twice;
- `MQ` drops mapping quality 255 from both numerator and denominator, so an all-unavailable matrix writes the four characters `NaN`;
- `AD`'s marginalisation order is a `HashMap`'s, and the tie in `searchBestAllele` is how that order reaches the output.

- [ ] **Step 3: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test -p gatk-annotation
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps -p gatk-annotation 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 4: Record and check**

```bash
cd /Users/benjamin/GATK
for f in flow pedigree fragment_counts site_statistics tandem_repeat heterozygosity mapping_quality depth_per_allele; do
  echo "crates/gatk-annotation/src/$f.rs 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
```

- [ ] **Step 5: Commit and open the PR for Tasks 10 to 13**

```bash
cd /Users/benjamin/GATK
git add crates/gatk-annotation/src tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the remaining annotations

Eight files, and the same move as the last task: the findings were in the
headers and are now at the items that implement them. The randomised QD
branch is explained where its refusal is defined, the three depth tests
that all read the child are commented on the three lines, and the
HashMap order that reaches AD is named at the marginalisation rather than
three modules away.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
git push -u origin docs/commenting-tranche-2
gh workflow run "CI" --ref docs/commenting-tranche-2
```

Wait for green, open the PR, merge.

---

### Task 14: Read filters, walkers and the corpus reader

**Files:**
- Modify: `crates/gatk-readfilter/src/lib.rs`
- Modify: `crates/gatk-readfilter/src/counting.rs`
- Modify: `crates/gatk-tools/src/lib.rs`
- Modify: `crates/gatk-tools/src/locus_walker.rs`
- Modify: `crates/gatk-tools/src/read_walker.rs`
- Modify: `crates/gatk-tools/src/interval_walker.rs`
- Modify: `crates/gatk-tools/src/print_reads.rs`
- Modify: `crates/gatk-corpus/src/lib.rs`
- Modify: `tools/audit/commented.txt`

**Interfaces:**
- Consumes: the `--record` refusal from Task 1.
- Produces: nothing new.

- [ ] **Step 1: Start a new branch**

```bash
cd /Users/benjamin/GATK
git checkout main && git pull
git checkout -b docs/commenting-tranche-3
```

- [ ] **Step 2: Read the references**

`engine/filters/` (32 files, 55 filters), `CountingReadFilter`, `LocusWalker`, `ReadWalker`, `IntervalWalker`, `PrintReads`.

- [ ] **Step 3: Comment the eight files**

`gatk-readfilter/src/lib.rs` is 1091 lines holding 55 filters. Each filter is short, so the item comment is one or two sentences: which reads it keeps, and the one thing about it that is not obvious from its name. Several filters disagree with their names; those are the ones to spend words on.

`counting.rs` must explain that the counts are per filter **and** per composite, and that the composite's count is not the sum of its parts because a read stops at the first rejection.

The walkers must say in which order the callbacks fire and what is guaranteed about the context each receives, since a walker's contract is the only thing a tool author can rely on.

- [ ] **Step 4: Verify nothing changed**

```bash
cd /Users/benjamin/GATK
cargo fmt --all
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo doc --no-deps --workspace 2>&1 | grep -i "warning" || echo "no doc warnings"
```

- [ ] **Step 5: Record, check, and confirm the whole tree is listed**

```bash
cd /Users/benjamin/GATK
for f in crates/gatk-readfilter/src/lib.rs crates/gatk-readfilter/src/counting.rs \
         crates/gatk-tools/src/lib.rs crates/gatk-tools/src/locus_walker.rs \
         crates/gatk-tools/src/read_walker.rs crates/gatk-tools/src/interval_walker.rs \
         crates/gatk-tools/src/print_reads.rs crates/gatk-corpus/src/lib.rs; do
  echo "$f 0.0" >> tools/audit/commented.txt
done
python3 tools/audit/comment_density.py --record
python3 tools/audit/comment_density.py --check
python3 tools/audit/comment_density.py | head -4
```

Expected: the summary line reads `70 of 122 files brought to the standard`. The 52 not listed are the test files, which this plan puts out of scope. If the number is not 70, a source file was missed. Find it with:

```bash
comm -23 <(ls crates/*/src/*.rs | sort) <(grep -v '^#' tools/audit/commented.txt | awk 'NF {print $1}' | sort)
```

- [ ] **Step 6: Commit**

```bash
cd /Users/benjamin/GATK
git add crates tools/audit/commented.txt
git commit -m "$(cat <<'EOF'
Explain the read filters, the walkers and the corpus reader

Fifty-five filters in one file, each getting the one sentence that says
which reads it keeps and the one thing that is not obvious from its name.
Several of them disagree with their names, and those got the words. The
counting wrapper now records that a composite's count is not the sum of
its parts, because a read stops at the first rejection.

Every source file in the repository is now on the list, 70 of them. The 52
files not listed are tests, which are out of scope by design: a fixture
table explains itself.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 15: Close the count on the roadmap

**Files:**
- Modify: `ROADMAP.md`
- Modify: `docs/COMMENTING.md`

**Interfaces:**
- Consumes: the finished list from Task 14.
- Produces: nothing.

- [ ] **Step 1: Update the roadmap count**

In `ROADMAP.md`, find the section "Explaining the code as it is written" and replace the second bullet:

```markdown
- [ ] the tranches, in the order a reader needs them: the engine's numeric primitives, then the
      likelihood matrix, then the annotations, then the readers and writers. **3 of 122 files** in
      this repository, 332 across the three
```

with:

```markdown
- [x] every source file in this repository: **70 of 70**, the remaining 52 being tests, which are
      out of scope because a fixture table explains itself
- [ ] `htsjdk-rs` (69 source files) and `picard-rs` (53), each with its own copy of the standard
      and its own list
```

- [ ] **Step 2: Record the outcome in the standard**

In `docs/COMMENTING.md`, replace the "Order of work" section's closing line with the measured outcome, so the document says what happened rather than what was intended:

```markdown
`gatk-rs` finished this order in three tranches. The overall density went from 0.266 to the figure
`tools/audit/comment_density.py` now reports; the per-file floors are in
`tools/audit/commented.txt`. `htsjdk-rs` and `picard-rs` carry their own copies of this file and
their own lists.
```

Run `python3 tools/audit/comment_density.py | head -3` and paste the real overall figure into that paragraph rather than leaving it as a reference to the tool.

- [ ] **Step 3: Verify**

```bash
cd /Users/benjamin/GATK
python3 tools/audit/comment_density.py --check
python3 tools/conformance/generate_ci.py
git diff --stat .github/workflows/ci.yml
```

Expected: check passes; the generated CI is unchanged, since this task touches no manifest. If `ci.yml` shows a diff, something else changed it and that diff belongs in another commit.

- [ ] **Step 4: Commit, push, PR**

```bash
cd /Users/benjamin/GATK
git add ROADMAP.md docs/COMMENTING.md
git commit -m "$(cat <<'EOF'
Close the count

67 of 67 source files. The roadmap now carries the outcome rather than the
intention, and the standard records the density the tree actually reached.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
)"
git push -u origin docs/commenting-tranche-3
gh workflow run "CI" --ref docs/commenting-tranche-3
```

Wait for green, open the PR, merge.

---

## Self-review

**Spec coverage.** The request was "add ultra-detailed comments to all the code, what/how/why, understandable without knowing Rust, useful to an LLM". Task 1 protects the measurement; Tasks 2 to 14 cover all 70 source files in `gatk-rs` by name (38 engine files in Tasks 2 and 5 to 9 plus the two already done, 22 annotation files in Tasks 10 to 13, 8 in Task 14), with no file appearing twice and none omitted (verified against `ls crates/*/src/*.rs` in Task 14 Step 5); Task 15 records the outcome. The "without knowing Rust" requirement is met by the mechanics-comment rule in `docs/COMMENTING.md` and by the idiom table, and every task names the specific Rust constructs its files use. Test files are explicitly out of scope with a stated reason, which is a narrowing and is called out here rather than hidden.

**Gap, stated rather than papered over:** `htsjdk-rs` (69 source files) and `picard-rs` (53) are **not** covered by this plan. They are separate repositories with separate CI and separate merge queues, and a single plan spanning three repositories could not produce a reviewable PR. Each needs its own plan, and Task 15 puts them on the roadmap so they are visible rather than forgotten.

**Placeholder scan.** No step says "add appropriate comments" without naming what must be explained: every commenting task lists the specific facts that must appear, each traceable to a line of Java or a golden row. The one instruction that is deliberately open is "read the reference named in the module header", which is an action, not a placeholder: the file names its own reference and the plan cannot enumerate 55 read filters' quirks without reading them.

**Type consistency.** `comment_density.py`'s three flags are used identically everywhere: `--record` in every task, `--check` in every task and in CI, `--allow-lower` only in Task 1's test. `--root` is introduced in Task 1 and used only by the test file created in the same task. The list file is `tools/audit/commented.txt` in every task. The commit trailer is identical in all fifteen commit messages.

---

## Execution handoff

Plan complete and saved to `docs/superpowers/plans/2026-08-01-commenting-gatk-rs.md`. Two execution options:

**1. Subagent-Driven (recommended)** - a fresh subagent per task, review between tasks, fast iteration. Suits this plan well: each task is independent, the ratchet catches regressions mechanically, and a fresh reader per task is exactly the audience the comments are for.

**2. Inline Execution** - execute tasks in this session, batch execution with checkpoints for review.

Which approach?
