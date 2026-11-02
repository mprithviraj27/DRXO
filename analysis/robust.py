"""Main-text robustness table for the held-out suite: uncorrected and Holm-corrected W/T/L, Vargha-Delaney large effects.
Also checks that the counts are unchanged when the tie-corrected Mann-Whitney U test replaces the rank-sum normal approximation."""
import sys; sys.path.insert(0, "/home/claude/ms")
import numpy as np
from scipy.stats import ranksums, mannwhitneyu
from build_tables import load14, ORDER14
from holm import holm
def A(x, y):
    x, y = np.asarray(x)[:, None], np.asarray(y)[None, :]; return float(((x < y) + 0.5 * (x == y)).mean())
def wtl(R, a, test, corr):
    ps, sg = [], []
    for f in range(1, 31):
        x, y = np.array(R["DRXO"][f]), np.array(R[a][f])
        if np.allclose(np.median(x), np.median(y), atol=1e-8) and np.allclose(x.mean(), y.mean(), rtol=1e-3, atol=1e-8):
            ps.append(1.0); sg.append(0); continue
        p = ranksums(x, y).pvalue if test == "rs" else (mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").pvalue if np.ptp(np.r_[x, y]) > 0 else 1.0)
        ps.append(p); sg.append(1 if (np.median(x) <= np.median(y) and x.mean() < y.mean()) else (-1 if (np.median(x) >= np.median(y) and x.mean() > y.mean()) else 0))
    q = holm(ps) if corr else np.array(ps)
    w = sum(1 for a_, s in zip(q, sg) if a_ < 0.05 and s == 1); l = sum(1 for a_, s in zip(q, sg) if a_ < 0.05 and s == -1)
    return w, 30 - w - l, l
algs = ORDER14[1:] + ["SR-skel"]
L = ["| DRXO vs | D | Wilcoxon | Holm-corrected | Large effects (DRXO / rival) | Median A |", "|---|---|---|---|---|---|"]
diff = []
for a in algs:
    for D in (10, 30, 50):
        R = load14(D)
        u = wtl(R, a, "rs", False); h = wtl(R, a, "rs", True); m = wtl(R, a, "mw", False)
        if m != u: diff.append((a, D, u, m))
        As = np.array([A(R["DRXO"][f], R[a][f]) for f in range(1, 31)])
        L.append(f"| {a} | {D} | {u[0]} / {u[1]} / {u[2]} | {h[0]} / {h[1]} / {h[2]} | {int((As >= 0.714).sum())} / {int((As <= 0.286).sum())} | {np.median(As):.2f} |")
open("/home/claude/ms/tables/robust.md", "w").write("\n".join(L) + "\n"); print("\n".join(L)); print("MW differs:", diff)
