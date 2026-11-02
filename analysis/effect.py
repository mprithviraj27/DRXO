"""Vargha-Delaney A (probability that a DRXO run is better than a rival run; ties count 1/2) on the held-out CEC 2014
suite. Thresholds (Vargha and Delaney 2000): |A - 0.5| >= 0.214 large, >= 0.138 medium, >= 0.064 small."""
import sys; sys.path.insert(0, "/home/claude/ms")
import numpy as np
from build_tables import load14, ORDER14
def A(x, y):
    x, y = np.asarray(x)[:, None], np.asarray(y)[None, :]
    return float(((x < y) + 0.5 * (x == y)).mean())
algs = [a for a in ORDER14[1:]] + (["SR-skel"] if "SR-skel" in load14(30) else [])
L = ["| DRXO vs | D | Large effect for DRXO | Large effect for rival | Medium (either way) | Small or negligible | Median A |", "|---|---|---|---|---|---|---|"]
for a in algs:
    for D in (10, 30, 50):
        R = load14(D)
        if a not in R or any(len(R[a].get(f, [])) < 25 for f in range(1, 31)): continue
        As = np.array([A(R["DRXO"][f], R[a][f]) for f in range(1, 31)])
        big_d = int((As >= 0.714).sum()); big_r = int((As <= 0.286).sum())
        med = int(((np.abs(As - 0.5) >= 0.138) & (np.abs(As - 0.5) < 0.214)).sum())
        L.append(f"| {a} | {D} | {big_d} | {big_r} | {med} | {30 - big_d - big_r - med} | {np.median(As):.2f} |")
open("/home/claude/ms/tables/effect.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
