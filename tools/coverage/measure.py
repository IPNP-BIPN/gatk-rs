#!/usr/bin/env python3
"""Turn `run_array.py`'s summary into a committable measurement.

The runner prints what it found; nothing until now kept the NUMBER, so the dashboard's
argument-coverage column has read `not measured` for every tool in this repository. This runs the
array against both sides and writes `tools/coverage/measured.json`, which the dashboard reads.

What is recorded, per tool:

* `rows`: the array's size at the declared strength;
* `rejected`: how many rows the REFERENCE refused. A rejected row is not a failed run: the
  reference refuses a block-compressed input whose output is not a `.tbi`, and reproducing that
  refusal is as much a part of the claim as reproducing an index;
* `distinct_outputs`: how many different outputs the accepted rows produced, which is the number
  that says whether the array is testing anything at all;
* `matched`: rows where the port answered exactly as the reference did, refusals included;
* `share`: `matched / rows`, which is what the dashboard prints.

What is NOT recorded is the outputs themselves. A Tribble index carries the last-modified time of
the file it was built from, and the corpus is rebuilt on every run, so two runs of this produce two
sets of digests. The comparison inside one run is exact; a golden of those digests would flake.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "tools" / "coverage" / "run_array.py"
MEASURED = REPO / "tools" / "coverage" / "measured.json"

SUMMARY = re.compile(
    r"tool=(?P<tool>\S+) t=(?P<t>\d+) rows=(?P<rows>\d+) rejected=(?P<rejected>\d+) "
    r"distinct_outputs=(?P<distinct>\d+) matched=(?P<matched>\d+) share=(?P<share>[\d.]+)"
)

COMMENT = [
    "Produced by tools/coverage/measure.py from run_array.py's own summary, in the pinned",
    "container on real x86-64. `share` is the fraction of the array's rows on which the port",
    "answered exactly as the reference did, refusals included.",
    "",
    "It is not a byte-identity claim over the argument surface: an array covers PAIRS of argument",
    "values, and `distinct_outputs` says how many of those pairs the corpus can actually observe.",
    "A tool whose accepted rows all produce one output has an array that covers its arguments",
    "without testing them.",
]


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tool", action="append", required=True)
    parser.add_argument("--port", required=True, help="the linux/amd64 port binary")
    parser.add_argument("--t", type=int, default=2)
    parser.add_argument(
        "--corpus-dir",
        help="where to write each tool's row-by-row corpus, which is what a divergence is read from",
    )
    options = parser.parse_args(argv)

    tools = {}
    for tool in options.tool:
        command = [
            sys.executable, str(RUNNER), "--tool", tool,
            "--t", str(options.t), "--port", options.port,
        ]
        if options.corpus_dir:
            corpus = Path(options.corpus_dir)
            corpus.mkdir(parents=True, exist_ok=True)
            command += ["--corpus", str(corpus / f"{tool}.t{options.t}.txt")]
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )
        print(result.stdout[-4000:])
        if result.returncode != 0:
            print(result.stderr[-2000:])
            raise SystemExit(f"the array run failed for {tool}")
        match = None
        for line in result.stdout.split("\n"):
            found = SUMMARY.match(line.strip())
            if found:
                match = found
        if match is None:
            raise SystemExit(f"no summary line for {tool}")
        tools[tool] = {
            "rows": int(match["rows"]),
            "accepted": int(match["rows"]) - int(match["rejected"]),
            "rejected": int(match["rejected"]),
            "distinct_outputs": int(match["distinct"]),
            "matched": int(match["matched"]),
            "t": int(match["t"]),
            "share": float(match["share"]),
        }

    MEASURED.write_text(
        json.dumps({"$comment": COMMENT, "tools": tools}, indent=2) + "\n"
    )
    print(f"wrote {MEASURED.relative_to(REPO)}: {len(tools)} tools")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
