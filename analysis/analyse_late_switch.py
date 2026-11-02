import json, sys, numpy as np
from scipy.stats import ranksums
fl = lambda e: 0.0 if e <= 1e-8 else float(e)
for D in (30, 50):
    try: P = json.load(open(f"/home/claude/uni/posthoc_c14_d{D}.json"))
    except FileNotFoundError: continue
    B = json.load(open(f"/home/claude/heldout14/py_d{D}.json"))
    L = {}; X = {}
    for f, a, s, e in P:
        if a == "post_late20" and s < 10: L.setdefault(f, []).append(fl(e))
    for f, a, s, e in B:
        if a == "DRXO-U" and s < 10: X.setdefault(f, []).append(fl(e))
    c = [0, 0, 0]; det = []
    for f in range(1, 31):
        if len(L.get(f, [])) < 10: continue
        x, y = np.array(L[f]), np.array(X[f])
        if np.allclose(np.median(x), np.median(y), atol=1e-8) and np.allclose(x.mean(), y.mean(), rtol=1e-3, atol=1e-8): r = 1
        else:
            p = ranksums(x, y).pvalue
            r = 0 if (p < .05 and np.median(x) <= np.median(y) and x.mean() < y.mean()) else (2 if (p < .05 and np.median(x) >= np.median(y) and x.mean() > y.mean()) else 1)
        c[r] += 1
        if r != 1: det.append(("+" if r == 0 else "-") + f"F{f}")
    print(f"D{D} late20 vs DRXO (10 runs): {c}  {det}")
    for f in (1, 2, 3):
        if f in L: print(f"   F{f}: late mean {np.mean(L[f]):.3g} solved {sum(v==0 for v in L[f])}/10 | DRXO mean {np.mean(X[f]):.3g} solved {sum(v==0 for v in X[f])}/10")
