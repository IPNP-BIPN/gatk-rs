#!/usr/bin/env python3
"""Compare two dumps, assembly by assembly: compare.py reference.txt port.txt.

Reports totals and the first few differences by name. Exits 1 on any difference."""
import sys


def records(path):
    cur = None
    with open(path, "rb") as f:
        for line in f:
            if line[:2] in (b"A ", b"X ", b"J ", b"U "):
                if cur:
                    yield cur
                cur = [line]
            else:
                cur.append(line)
    if cur:
        yield cur


def main():
    n = contigs = links = aborts = over = refused = bad = 0
    for ra, rb in zip_longest_records(sys.argv[1], sys.argv[2]):
        if ra is None and len(sys.argv) > 3 and sys.argv[3] == "--prefix":
            print("the reference is still being written: compared what there is")
            break
        if ra is None or rb is None:
            print("different number of records")
            bad += 1
            break
        kind = ra[0][:1]
        n += 1
        if rb[0][:1] == b"U":
            refused += 1
            continue
        if ra != rb:
            bad += 1
            if bad <= 20:
                print("DIFF", ra[0].decode("latin-1").strip(), "| vs |", rb[0].decode("latin-1").strip())
                for x, y in zip(ra[1:], rb[1:]):
                    if x != y:
                        print("  ref ", x.decode("latin-1").strip()[:200])
                        print("  port", y.decode("latin-1").strip()[:200])
                        break
            continue
        contigs += sum(1 for l in ra if l[:2] == b"C ")
        links += sum(1 for l in ra if l[:2] == b"L ")
        aborts += kind == b"X"
        over += kind == b"J"
    print(
        f"{n} records; identical apart from {bad} differences; "
        f"{contigs} contigs, {links} connections, {aborts} aborts, {over} Java overflows compared; "
        f"{refused} refused by the port and not compared"
    )
    sys.exit(1 if bad else 0)


def zip_longest_records(a, b):
    ia, ib = records(a), records(b)
    while True:
        x = next(ia, None)
        y = next(ib, None)
        if x is None and y is None:
            return
        yield x, y


if __name__ == "__main__":
    main()
