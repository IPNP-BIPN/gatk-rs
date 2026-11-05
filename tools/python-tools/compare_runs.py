#!/usr/bin/env python3
"""Which files of two run directories have the same bytes, and how far apart the rest are.

Usage: compare_runs.py <first run> <second run> [more runs...]
Files are compared by sha256; for a pair of hashes that differ the size of the difference (bytes,
first offset) is printed. Nothing is read as HDF5 here, so this can run on the host.
"""
import hashlib
import os
import sys


def walk(root):
    found = {}
    for directory, _, names in os.walk(root):
        for name in names:
            path = os.path.join(directory, name)
            found[os.path.relpath(path, root)] = path
    return found


def main():
    runs = sys.argv[1:]
    trees = [walk(run) for run in runs]
    names = sorted(set().union(*trees))
    print("file\t" + "\t".join(os.path.basename(run.rstrip("/")) for run in runs) + "\tbytes differing vs first\tfirst offset")
    for name in names:
        hashes = []
        contents = []
        for tree in trees:
            if name in tree:
                with open(tree[name], "rb") as f:
                    content = f.read()
                contents.append(content)
                hashes.append(hashlib.sha256(content).hexdigest()[:10])
            else:
                contents.append(None)
                hashes.append("-")
        notes = []
        first = contents[0]
        for other in contents[1:]:
            if first is None or other is None:
                notes.append("-")
            elif first == other:
                notes.append("0")
            else:
                different = sum(a != b for a, b in zip(first, other)) + abs(len(first) - len(other))
                offset = next((i for i, (a, b) in enumerate(zip(first, other)) if a != b), min(len(first), len(other)))
                notes.append("%d@%d" % (different, offset))
        print(name + "\t" + "\t".join(hashes) + "\t" + ",".join(notes))


if __name__ == "__main__":
    main()
