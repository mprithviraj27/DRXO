"""Holm-corrected pairwise comparisons on the held-out CEC 2014 suite (correction over the 30 functions of each
algorithm pair and dimension). A function counts as a win/loss only if the Holm-adjusted p < 0.05 and median and mean agree."""
import sys; sys.path.insert(0, "/home/claude/ms")
import numpy as np
from scipy.stats import ranksums
from build_tables import load14, ORDER14
def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for i, j in enumerate(o):
        run = max(run, min(1, (m - i) * p[j])); adj[j] = run
    return adj
L = ["| DRXO vs | D = 10 | D = 30 | D = 50 |", "|---|---|---|---|"]
for a in ORDER14[1:]:
    cells = []
    for D in (10, 30, 50):
        R = load14(D); ps, sign = [], []
        for f in range(1, 31):
            x, y = np.array(R["DRXO"][f]), np.array(R[a][f])
            if np.allclose(np.median(x), np.median(y), atol=1e-8) and np.allclose(x.mean(), y.mean(), rtol=1e-3, atol=1e-8):
                ps.append(1.0); sign.append(0); continue
            ps.append(ranksums(x, y).pvalue)
            sign.append(1 if (np.median(x) <= np.median(y) and x.mean() < y.mean()) else (-1 if (np.median(x) >= np.median(y) and x.mean() > y.mean()) else 0))
        adj = holm(ps); w = sum(1 for q, s in zip(adj, sign) if q < 0.05 and s == 1); l = sum(1 for q, s in zip(adj, sign) if q < 0.05 and s == -1)
        cells.append(f"{w} / {30 - w - l} / {l}")
    L.append(f"| {a} | " + " | ".join(cells) + " |")
open("/home/claude/ms/tables/holm.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
