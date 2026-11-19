#!/usr/bin/env python3
"""Time one tool's command lines inside the pinned container. Run by `bench.py`, not by hand.

The spec names three command lines over the same row: the reference cold (`java ... Main`), the
port, and the reference's steady-state loop (`Steady.java`). Every run is a child of this process,
so `wait4` returns exactly that child's CPU time and peak resident set, and each command is
`exec`ed so the child IS the tool rather than a shell around it.

The two sides alternate, oracle first on even repetitions and port first on odd ones, so a drift
in the host (a neighbour's load, a thermal step) lands on both rather than on whichever ran last.
Repetition zero warms the page cache for both and is reported but not summarised.
"""

import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

OUT = Path("/work/out")


def prepare(out=OUT):
    """An empty output directory and the two temporary ones, as `run_array.py` gives every row."""
    for path in (out, Path("/work/tmp"), Path("/work/tmp2")):
        if path.exists():
            for stale in path.iterdir():
                if stale.is_dir() and not stale.is_symlink():
                    shutil.rmtree(stale, ignore_errors=True)
                else:
                    stale.unlink(missing_ok=True)
        path.mkdir(parents=True, exist_ok=True)


def oom_kills():
    """How many processes the container's memory limit has killed so far (cgroup v2), or None."""
    try:
        for line in Path("/sys/fs/cgroup/memory.events").read_text().splitlines():
            if line.startswith("oom_kill "):
                return int(line.split()[1])
    except OSError:
        pass
    return None


def run(command, timeout, out=OUT, log=None):
    """One run: wall seconds, CPU seconds (user + system), peak RSS in KiB, exit code.

    A run the timeout kills says so (`timed_out`), and one the container's memory limit kills
    says so too (`oom_killed`): at scale those are results, not errors to retry.
    """
    prepare(out)
    kills = oom_kills()
    expired = threading.Event()
    start = time.perf_counter()
    sink = open(log, "w") if log else subprocess.DEVNULL
    child = subprocess.Popen(
        ["bash", "-c", f"exec {command}"],
        stdin=subprocess.DEVNULL,
        stdout=sink,
        stderr=sink,
        start_new_session=True,
    )
    def expire():
        expired.set()
        os.killpg(child.pid, signal.SIGKILL)

    timer = threading.Timer(timeout, expire)
    timer.start()
    _, status, usage = os.wait4(child.pid, 0)
    wall = time.perf_counter() - start
    timer.cancel()
    if log:
        sink.close()
    # Popen still holds the pid; tell it the child is gone so it does not wait again.
    child.returncode = os.waitstatus_to_exitcode(status)
    after = oom_kills()
    return {
        "wall": wall,
        "cpu": usage.ru_utime + usage.ru_stime,
        "rss_kib": usage.ru_maxrss,
        "exit": child.returncode,
        "timed_out": expired.is_set(),
        "oom_killed": kills is not None and after is not None and after > kills,
    }


def steady(command, iterations, timeout):
    """The reference's in-process loop, read back from the report `Steady.java` writes."""
    report = Path("/work/steady-report.txt")
    report.unlink(missing_ok=True)
    outer = run(f"{command.replace('@REPORT@', str(report))}", timeout)
    lines = report.read_text().splitlines() if report.exists() else []
    loop = []
    for line in lines:
        index, wall_ns, cpu_ns, status = line.split("\t")
        loop.append(
            {"wall": int(wall_ns) / 1e9, "cpu": int(cpu_ns) / 1e9, "status": status}
        )
    return {"process": outer, "iterations": loop, "requested": iterations}


def main(spec_path):
    spec = json.loads(Path(spec_path).read_text())
    reps, timeout = spec["reps"], spec["timeout"]
    runs = {"oracle": [], "port": []}
    if spec.get("scale"):
        # One run per side at scale, each into its own kept directory; the inputs are read once
        # first so neither side pays the page cache for the other.
        for path in spec.get("warm", []):
            subprocess.run(["cat", path], stdout=subprocess.DEVNULL, check=False)
        for side in ("oracle", "port"):
            out = Path(spec["out"][side])
            runs[side].append(run(spec[side], timeout, out, out.with_name(f"{side}.log")))
        reps = -1
    for rep in range(reps + 1):
        order = ("oracle", "port") if rep % 2 == 0 else ("port", "oracle")
        for side in order:
            runs[side].append(run(spec[side], timeout))
    result = {"cold": runs}
    if spec.get("steady"):
        result["steady"] = steady(spec["steady"], spec["iterations"], timeout)
    cpu = ""
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    result["host"] = {"cpu": cpu, "cores": os.cpu_count()}
    print(json.dumps(result))


if __name__ == "__main__":
    main(sys.argv[1])
