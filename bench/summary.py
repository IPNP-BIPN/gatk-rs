"""The bench's JSON files as one markdown table per runner and core count."""
import json, statistics as st, sys
for path in sys.argv[1:]:
    runs = json.load(open(path)); host = runs[0]["host"]
    print(f"\n### {path}: {host['cpu']}, {host['cores']} cores\n")
    print("| Case | Side | Wall | CPU | Peak RSS |\n|---|---|---:|---:|---:|")
    seen = {}
    for r in runs:
        for key, side in (("oracle", r["a"]), ("port", r["b"])):
            x = r["cold"][key][1:]
            seen.setdefault((r["case"], side), []).append((st.median(i["wall"] for i in x), st.median(i["cpu"] for i in x), st.median(i["rss_kib"] for i in x) / 1024))
    for (case, side), v in seen.items():
        w, c, m = (st.median(t[k] for t in v) for k in range(3))
        print(f"| {case} | {side} | {w:.2f} s | {c:.2f} s | {m:.0f} MiB |")
    for r in runs:
        if r["a"] == "old":
            ratio = st.median(b["cpu"] / a["cpu"] for a, b in zip(r["cold"]["oracle"][1:], r["cold"]["port"][1:]))
            print(f"\n{r['case']}: new/old CPU, paired median {ratio:.2f}")
