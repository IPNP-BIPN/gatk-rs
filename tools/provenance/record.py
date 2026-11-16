#!/usr/bin/env python3
"""Write what a conformance run was made of, so a third party can tell whether theirs is the same.

A green run says the port and the reference agreed on some inputs, under some toolchain, in some
container. Without the record of which, the byte claim cannot be checked by anybody who was not
there, so CI writes this file and keeps it with the run's artefacts:

* the commit, and the pins of what it builds against, from `Cargo.lock`;
* the reference: the versions the inventory was generated from, and the oracle image's id;
* the toolchain: `rust-toolchain.toml`, `rustc -vV`, and the oracle's `java -version`;
* the compression provider, GKL or not, as each suite's Java properties set it, since a BAM's bytes
  depend on it;
* a SHA-256 of every input the comparison reads: the goldens, the covering arrays, the fixtures'
  generators and the workflow;
* on a device-path run (`GATK_RS_GPU` set, Milestone GPU.3), the device, its driver and the
  toolkit, because a kernel's result is a property of all three.

A value that could not be read is recorded as `null`: a record that left the key out could not be
told from an older format, and one that guessed would be worse than none.

    python3 tools/provenance/record.py --out provenance.json
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "tools" / "conformance" / "manifest.json"
INVENTORY = REPO / "tools" / "inventory" / "generated" / "inventory.json"

# Everything the comparison reads besides the build, relative to the repository.
INPUT_GLOBS = [
    "Cargo.lock",
    "rust-toolchain.toml",
    ".github/workflows/ci.yml",
    "tools/oracle/**/*",
    "tools/conformance/manifest.json",
    "tools/coverage/fixtures.json",
    "tools/coverage/MakeFixtures.java",
    "tools/coverage/arrays/**/*.json",
    "tools/fuzz/FuzzDriver.java",
    "tools/fuzz/run_fuzz.py",
]


def command(*argv):
    """The output of a command, or None when it cannot run here."""
    try:
        done = subprocess.run(argv, capture_output=True, text=True, check=False, cwd=REPO)
    except OSError:
        return None
    return (done.stdout + done.stderr).strip() if done.returncode == 0 else None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pins():
    """The git dependencies `Cargo.lock` resolved, each with its commit."""
    found = {}
    lock = (REPO / "Cargo.lock").read_text()
    for name, version, source in re.findall(
        r'name = "([^"]+)"\nversion = "([^"]+)"\nsource = "git\+([^"]+)"', lock
    ):
        found[name] = {"version": version, "source": source}
    return found


def compression(manifest):
    """Each distinct set of provider properties, with how many suites run under it."""
    counts = {}
    for suite in manifest["suites"]:
        props = suite.get("java_props", manifest.get("default_java_props", []))
        chosen = tuple(sorted(p for p in props if "samjdk." in p or "gkl" in p.lower()))
        counts[chosen] = counts.get(chosen, 0) + 1
    return [
        {"properties": list(props), "suites": count}
        for props, count in sorted(counts.items())
    ]


def inputs(manifest):
    """Path to SHA-256 of every input; a declared golden that is not there is recorded as null."""
    paths = set()
    for pattern in INPUT_GLOBS:
        paths.update(path for path in REPO.glob(pattern) if path.is_file())
    hashes = {str(path.relative_to(REPO)): sha256(path) for path in paths}
    for suite in manifest["suites"]:
        for case in suite["cases"]:
            golden = case.get("golden")
            if golden is not None:
                path = REPO / golden
                hashes[golden] = sha256(path) if path.is_file() else None
    return dict(sorted(hashes.items()))


def device():
    """The GPU a device-path run used, or None on a CPU run.

    Metal ships inside the operating system: its driver and the compiler that builds the kernels
    from source at run time are the OS build, and the bindings are the `metal` crate's version.
    """
    requested = os.environ.get("GATK_RS_GPU")
    if not requested or requested == "0":
        return None
    displays = command("system_profiler", "-json", "SPDisplaysDataType")
    gpus = json.loads(displays).get("SPDisplaysDataType", []) if displays else []
    lock = (REPO / "Cargo.lock").read_text()
    bindings = re.search(r'name = "metal"\nversion = "([^"]+)"', lock)
    return {
        "requested": requested,
        "devices": [
            {
                "model": gpu.get("sppci_model"),
                "cores": gpu.get("sppci_cores"),
                "metal": gpu.get("spdisplays_mtlgpufamilysupport"),
                "vendor": gpu.get("spdisplays_vendor"),
            }
            for gpu in gpus
        ] or None,
        "driver": command("sw_vers"),
        "toolkit": {"metal-crate": bindings.group(1) if bindings else None},
    }


def record():
    manifest = json.loads(MANIFEST.read_text())
    inventory = json.loads(INVENTORY.read_text())
    image = manifest["oracle"]["image"]
    run = {
        name: os.environ.get(name)
        for name in ["GITHUB_REPOSITORY", "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "RUNNER_OS", "RUNNER_ARCH"]
    }
    return {
        "format": 2,
        "commit": command("git", "rev-parse", "HEAD"),
        "run": run,
        "reference": inventory.get("reference"),
        "oracle": {
            "image": image,
            "platform": manifest["oracle"].get("platform"),
            "image_id": command("docker", "image", "inspect", "--format", "{{.Id}}", image),
            "java": command(
                "docker", "run", "--rm", "--platform", manifest["oracle"]["platform"],
                "--entrypoint", "java", image, "-version",
            ),
        },
        "toolchain": {
            "rustc": command("rustc", "-vV"),
            "cargo": command("cargo", "-V"),
        },
        "pins": pins(),
        "compression": compression(manifest),
        "inputs": inputs(manifest),
        "device": device(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    Path(args.out).write_text(json.dumps(record(), indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
