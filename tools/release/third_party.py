#!/usr/bin/env python3
"""Print the licence of every crate linked into one binary, for the release archive.

usage: third_party.py <root-package> <target-triple>

MIT, BSD and Apache-2.0 all ask that a binary distribution carry the licence text, so each
non-workspace crate the root reaches through normal (linked) dependencies contributes its LICENSE,
COPYING and NOTICE files, including those of the C library a `-sys` crate bundles one level down.
"""
import json
import pathlib
import subprocess
import sys

NAMES = ("LICENSE", "LICENCE", "COPYING", "NOTICE")


def licence_files(root):
    for path in sorted(root.glob("*")) + sorted(root.glob("*/*")):
        if path.is_file() and path.name.upper().startswith(NAMES):
            yield path


def main(root_name, target):
    meta = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--format-version", "1", "--locked", "--filter-platform", target]))
    packages = {p["id"]: p for p in meta["packages"]}
    nodes = {n["id"]: n for n in meta["resolve"]["nodes"]}
    members = set(meta["workspace_members"])
    root = next(i for i in members if packages[i]["name"] == root_name)

    seen, todo = set(), [root]
    while todo:
        node = nodes[todo.pop()]
        for dep in node["deps"]:
            linked = any(k["kind"] is None for k in dep["dep_kinds"])
            if linked and dep["pkg"] not in seen:
                seen.add(dep["pkg"])
                todo.append(dep["pkg"])

    for pid in sorted(seen - members, key=lambda i: (packages[i]["name"], packages[i]["version"])):
        p = packages[pid]
        print("=" * 78)
        print(f"{p['name']} {p['version']}  --  {p.get('license') or 'see below'}")
        if p.get("repository"):
            print(p["repository"])
        crate_dir = pathlib.Path(p["manifest_path"]).parent
        for f in licence_files(crate_dir):
            print(f"\n--- {f.relative_to(crate_dir)}\n")
            print(f.read_text(errors="replace").rstrip())
        print()


if __name__ == "__main__":
    main(*sys.argv[1:])
