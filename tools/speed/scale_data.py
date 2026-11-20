#!/usr/bin/env python3
"""Derive the scale benchmark's data into DIR and check every file against its pinned sha256.

    python3 tools/speed/scale_data.py DIR          # derive what is missing, verify everything
    python3 tools/speed/scale_data.py DIR --pin    # derive, then write the hashes to scale.json

Milestone D, D1. `scale.json` names each file, the public source it comes from, the command that
derives it and the pinned image that command runs in. A file already in DIR with the pinned hash
is kept; any other is derived again. Verification fails on the first hash that does not match,
so a source that changed under its URL is caught before anything is timed on it.
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

SCALE = Path(__file__).resolve().parent / "scale.json"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def derive(entry, spec, data):
    command = entry["derive"].format(**{k: v["url"] for k, v in spec["sources"].items()})
    print(f"deriving {entry['name']}: {command}", flush=True)
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", "linux/amd64",
            "-v", f"{data.resolve()}:/data", "-w", "/data",
            spec["images"][entry["image"]], "sh", "-c", f"set -e; {command}",
        ],
    )
    if result.returncode != 0:
        (data / entry["name"]).unlink(missing_ok=True)
        raise SystemExit(f"could not derive {entry['name']}")


def main(argv):
    if not argv or argv[0].startswith("-"):
        raise SystemExit(__doc__)
    data, pin = Path(argv[0]), "--pin" in argv
    data.mkdir(parents=True, exist_ok=True)
    spec = json.loads(SCALE.read_text())
    for entry in spec["files"]:
        path = data / entry["name"]
        if not (path.exists() and entry.get("sha256") == sha256(path)):
            derive(entry, spec, data)
        digest = sha256(path)
        if pin:
            entry["sha256"], entry["bytes"] = digest, path.stat().st_size
        elif digest != entry.get("sha256"):
            raise SystemExit(f"{entry['name']}: sha256 {digest}, pinned {entry.get('sha256')}")
        print(f"ok {entry['name']} {digest} {path.stat().st_size}", flush=True)
    # htslib keeps the remote indexes it read (.crai, .tbi, .fai) beside the outputs; they are
    # the sources', not part of the data.
    listed = {entry["name"] for entry in spec["files"]}
    for stray in data.iterdir():
        if stray.name not in listed and stray.suffix in (".crai", ".tbi", ".fai", ".csi"):
            stray.unlink()
    if pin:
        # Read again, so an edit made to the steps while the data was derived is kept.
        current = json.loads(SCALE.read_text())
        current["files"] = spec["files"]
        SCALE.write_text(json.dumps(current, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
