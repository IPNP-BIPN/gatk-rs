# Three readings of a tool's arguments, and which one a command line obeys

Asking "what arguments does `CountReads` take?" has three answers in the reference, and they are
different numbers. Milestone C has to reproduce one of them exactly; the other two are useful as
cross-checks and dangerous as sources.

| reading | how it is taken | `CountReads` | `CountVariants` | `PrintReads` |
|---|---|---:|---:|---:|
| the tool's own parser | `((CommandLineProgram) tool).getCommandLineParser()` | **70** | **72** | **70** |
| the usage text | what `-h` prints, which is what `tools/inventory` reads | 42 | 44 | 42 |
| a parser built from the instance | `new CommandLineArgumentParser(tool)` | 38 | 40 | 38 |

## The one a command line obeys is the first

`CommandLineProgram.getCommandLineParser()` builds the parser with the tool's plugin descriptors
and the standard argument collections. A parser constructed straight from the instance has neither,
so it is missing the four read-filter arguments (`--read-filter`, `--disable-read-filter`,
`--inverted-read-filter`, `--disable-tool-default-read-filters`) and the whole common and advanced
surface: 32 arguments per tool, every one of which a real command line may name.

The `tool-argument-declarations` suite holds both numbers side by side for exactly this reason. A
declaration generated from the shorter list would refuse command lines the reference accepts, and
nothing in the shorter list is wrong: it is a subset, which is the failure mode that reads as
correct until someone passes `--read-filter`.

## The usage text is the second reading, and it is the cross-check

`tools/inventory` is generated from the reference's own CLI, so it is an independent view of the
same tool: independent of reflection, and produced by a different code path. Every argument it
documents is in the declarations, with the same `required` and the same default, and
`tools/declarations/generate.py` fails if that stops being true.

Two renderings differ, and both are the same fact written twice:

  * an empty collection is `null` from `NamedArgumentDefinition.getDefaultValueAsString()` and `[]`
    in the usage text;
  * an empty string scalar is the empty string in one and `""` in the other.

The generator names those two and stops on anything else. That is how `--arguments_file` and
`--gcs-project-for-requester-pays` were found to differ in rendering rather than in fact.

## What this means for Milestone C

  * **#819 generates, it does not transcribe.** Seventy declarations a tool, over three hundred
    tools, is not a hand-writing job, and a hand-written one has no way to be checked.
  * **The generator's input is the golden**, which is the reference's answer, and its cross-check is
    the inventory, which is the reference's other answer. Neither is this repository's opinion.
  * **`preflight` runs `generate.py --check`**, so a regenerated golden that changes a declaration
    fails the build rather than drifting.
