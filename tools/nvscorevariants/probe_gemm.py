import os, sys, torch, time
torch.set_num_threads(int(os.environ.get('T', '1')))
BIG = float(2**24)

def call(M, N, K, wcols, xcols, bias):
    """wcols/xcols: dict k -> value, uniform over rows and columns"""
    w = torch.zeros((M, K), dtype=torch.float32)
    x = torch.zeros((N, K), dtype=torch.float32)
    for k, v in wcols.items():
        w[:, k] = v
    for k, v in xcols.items():
        x[:, k] = v
    b = torch.full((M,), bias, dtype=torch.float32)
    return torch.addmm(b, x, w.t())  # (N, M)

def merged(M, N, K, big, j, l):
    """big = index or -1 for the bias. Returns (N,M) bool tensor: j,l merged before meeting big"""
    wc = {j: 1.0, l: 1.0}
    xc = {j: 1.0, l: 1.0}
    bias = 0.0
    if big < 0:
        bias = BIG
    else:
        wc[big] = BIG
        xc[big] = 1.0
    r = call(M, N, K, wc, xc, bias)
    return r != BIG

def classes(t):
    """summarise a (N,M) bool tensor as distinct row patterns"""
    return t

def adjacency(M, N, K):
    res = torch.zeros((K - 1, N, M), dtype=torch.bool)
    for j in range(K - 1):
        res[j] = merged(M, N, K, -1, j, j + 1)
    return res

def summarize(M, N, K):
    t0 = time.time()
    a = adjacency(M, N, K)
    print(f'M={M} N={N} K={K}: adjacency {time.time()-t0:.1f}s')
    # distinct patterns over (row,col)
    pats = {}
    for r in range(M):
        for c in range(N):
            breaks = tuple((~a[:, c, r]).nonzero().flatten().tolist())
            pats.setdefault(breaks, []).append((r, c))
    for br, who in pats.items():
        rows = sorted(set(x[0] for x in who))
        cols = sorted(set(x[1] for x in who))
        print(f'  breaks after k in {br[:40]}{"..." if len(br)>40 else ""} (n={len(br)}) rows {rows[0]}..{rows[-1]} ({len(rows)}) cols {cols[0]}..{cols[-1]} ({len(cols)}) positions {len(who)}')

which = sys.argv[1]
if which == 'adj':
    for shape in [(40, 64, 7), (40, 7, 7), (40, 64, 12200), (40, 7, 12200), (4, 64, 47), (4, 7, 47), (40, 1, 12200), (4, 1, 47), (40, 1, 7)]:
        summarize(*shape)
        sys.stdout.flush()
elif which == 'triple':
    M, N, K = [int(v) for v in sys.argv[2:5]]
    for r, c in [(0, 0), (M - 1, 0), (M - 9, 0), (0, N - 1), (M - 1, N - 1)]:
        print('row', r, 'col', c)
        for i in list(range(K)) + [-1]:
            out = []
            for j in range(K):
                for l in range(j + 1, K):
                    if j == i or l == i:
                        continue
                    if merged(M, N, K, i, j, l)[c, r]:
                        out.append((j, l))
            print('  big', i, 'merged before:', out)
        sys.stdout.flush()

def trees(M, N, K):
    keys = []
    results = []
    for i in list(range(K)) + [-1]:
        for j in range(K):
            for l in range(j + 1, K):
                if j == i or l == i:
                    continue
                keys.append((i, j, l))
                results.append(merged(M, N, K, i, j, l))
    stack = torch.stack(results)  # (T, N, M)
    sigs = {}
    for r in range(M):
        for c in range(N):
            sigs.setdefault(stack[:, c, r].numpy().tobytes(), []).append((r, c))
    print(f'M={M} N={N} K={K}: {len(sigs)} distinct trees')
    for sig, who in sigs.items():
        rows = sorted(set(x[0] for x in who))
        cols = sorted(set(x[1] for x in who))
        v = torch.frombuffer(bytearray(sig), dtype=torch.bool)
        by = {}
        for flag, key in zip(v.tolist(), keys):
            if flag:
                by.setdefault(key[0], []).append(key[1:])
        print(f'  rows {rows} cols {cols[0]}..{cols[-1]}({len(cols)}) positions {len(who)}')
        for i in list(range(K)) + [-1]:
            # present as partition of indices already merged: components
            pairs = by.get(i, [])
            parent = {}
            def find(a):
                parent.setdefault(a, a)
                while parent[a] != a:
                    a = parent[a]
                return a
            for a, b in pairs:
                parent[find(a)] = find(b)
            comp = {}
            for a in list(parent):
                comp.setdefault(find(a), []).append(a)
            print(f'    before {i:3}: ', sorted(sorted(x) for x in comp.values()))
    sys.stdout.flush()

if which == 'trees':
    for spec in sys.argv[2:]:
        trees(*[int(v) for v in spec.split(',')])

def subset(M, N, K, idx):
    keys = []
    results = []
    for i in idx + [-1]:
        for a, j in enumerate(idx):
            for l in idx[a + 1:]:
                if j == i or l == i:
                    continue
                keys.append((i, j, l))
                results.append(merged(M, N, K, i, j, l))
    stack = torch.stack(results)
    sigs = {}
    for r in range(M):
        for c in range(N):
            sigs.setdefault(stack[:, c, r].numpy().tobytes(), []).append((r, c))
    print(f'M={M} N={N} K={K} idx={idx}: {len(sigs)} distinct')
    for sig, who in sigs.items():
        rows = sorted(set(x[0] for x in who))
        cols = sorted(set(x[1] for x in who))
        v = torch.frombuffer(bytearray(sig), dtype=torch.bool)
        by = {}
        for flag, key in zip(v.tolist(), keys):
            if flag:
                by.setdefault(key[0], []).append(key[1:])
        print(f'  rows {rows[0]}..{rows[-1]}({len(rows)}) cols {cols} positions {len(who)}')
        for i in idx + [-1]:
            pairs = by.get(i, [])
            parent = {}
            def find(a):
                parent.setdefault(a, a)
                while parent[a] != a:
                    a = parent[a]
                return a
            for a, b in pairs:
                parent[find(a)] = find(b)
            comp = {}
            for a in list(parent):
                comp.setdefault(find(a), []).append(a)
            print(f'    before {i:5}: ', sorted(sorted(x) for x in comp.values()))
    sys.stdout.flush()

if which == 'subset':
    M, N, K = [int(v) for v in sys.argv[2:5]]
    subset(M, N, K, [int(v) for v in sys.argv[5].split(',')])

if which == 'classes':
    for M, K in [(40, 7), (4, 47), (64, 7), (4, 24), (40, 12200), (24, 15424)]:
        for N in list(range(1, 70)) + [127, 128, 129]:
            even = merged(M, N, K, 0, 1, 3)   # the odd chain before k=0
            seq = merged(M, N, K, 2, 0, 1)    # 0,1 together before k=2
            lines = {}
            for r in range(M):
                pattern = ''.join('e' if even[c, r] and not seq[c, r] else ('s' if seq[c, r] and not even[c, r] else '?') for c in range(N))
                lines.setdefault(pattern, []).append(r)
            desc = '; '.join(f'rows {v[0]}..{v[-1]}({len(v)}): {k}' for k, v in lines.items())
            print(f'M={M} K={K} N={N}: {desc}', flush=True)
