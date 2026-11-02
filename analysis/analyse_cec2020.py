"""CEC 2020 RC15-RC33: DRXO, DRX-off and the official results of five competition entrants.
Each run is reduced to a key with Deb's rules (feasible f < any infeasible; infeasible ranked by mean violation).
Pairwise: Wilcoxon rank-sum on keys (p < 0.05; median and mean of keys agree). Ranking: Friedman on median keys."""
import os, json, glob
import numpy as np
from scipy.stats import ranksums, rankdata, friedmanchisquare
from rwc import FSTAR, C
O = "/home/claude/rw20/official"
SRC = {"SASS": ("GECCO COM109 codes and results/Submitted Result/SASS_RC{k}_F.txt", "GECCO COM109 codes and results/Submitted Result/SASS_RC{k}_CV.txt"),
       "sCMAgES": ("GECCO COM108 codes and results/Submitted Results/sCMAgES{k}_F.txt", "GECCO COM108 codes and results/Submitted Results/sCMAgES{k}_CV.txt"),
       "COLSHADE": ("CEC E-24586 codes and results/Submitted Results/COLSHADE_RC{k}_F.txt", "CEC E-24586 codes and results/Submitted Results/COLSHADE_RC{k}_CV.txt"),
       "EnMODE": ("CEC E-24363 codes and results/Submitted result/EnMODE_RC{k}__F.txt", "CEC E-24363 codes and results/Submitted result/EnMODE_RC{k}__CV.txt"),
       "BiPop-εMAgES": ("CEC E-24443 codes and results/Submitted Results/BiPopEpsMAgES_RC{k}_F-wbC.txt", "CEC E-24443 codes and results/Submitted Results/BiPopEpsMAgES_RC{k}_CV-wbC.txt")}
PROBS = list(range(15, 34)); ALGS = ["DRXO", "DRX-off"] + list(SRC)
def r8(x): return float(f"{x:.7e}")   # official results are reported with 8 significant digits
def keyof(f, v): return r8(f) if v <= 1e-12 else C + r8(v)
R = {}   # R[alg][k] = (keys, f, v)
for k, a, s, f, v, *_ in json.load(open("/home/claude/rw20/res_rw20.json")):
    R.setdefault(a, {}).setdefault(k, []).append((np.nan if f is None else f, v))
for a, (pf, pv) in SRC.items():
    for k in PROBS:
        F = np.loadtxt(os.path.join(O, pf.format(k=k)))[-1]; V = np.loadtxt(os.path.join(O, pv.format(k=k)))[-1]
        R.setdefault(a, {})[k] = list(zip(F, V))
def keys(a, k): return np.array([keyof(f, v) for f, v in R[a][k]])
def wtl(x, y):
    if np.allclose(np.median(x), np.median(y), rtol=1e-8, atol=1e-10) and np.allclose(x.mean(), y.mean(), rtol=1e-6, atol=1e-10): return 1
    p = ranksums(x, y).pvalue
    if p < 0.05 and np.median(x) <= np.median(y) and x.mean() < y.mean(): return 0
    if p < 0.05 and np.median(x) >= np.median(y) and x.mean() > y.mean(): return 2
    return 1
if __name__ == "__main__":
    ok = [a for a in ALGS if a in R and all(len(R[a].get(k, [])) >= 25 for k in PROBS)]
    print("complete:", ok)
    rows = []
    for k in PROBS:
        line = f"RC{k} f*={FSTAR[k]:.6g} | "
        for a in ok:
            K = keys(a, k); fr = np.mean([v <= 1e-12 for _, v in R[a][k]]) * 100
            feas = K[K < C]
            line += f"{a}: best {feas.min() if len(feas) else np.nan:.6g} med {np.median(K):.6g} FR {fr:.0f}% | "
        rows.append(line)
    print("\n".join(rows))
    med = np.array([[np.median(keys(a, k)) for a in ok] for k in PROBS]); rk = np.array([rankdata(m) for m in med])
    print("mean rank (median keys):", ", ".join(f"{a} {r:.2f}" for a, r in zip(ok, rk.mean(0))), " Friedman p", friedmanchisquare(*rk.T).pvalue)
    for b in ok[1:]:
        c = [0, 0, 0]; det = []
        for k in PROBS:
            r = wtl(keys("DRXO", k), keys(b, k)); c[r] += 1
            if r != 1: det.append(("+" if r == 0 else "-") + f"RC{k}")
        print(f"DRXO vs {b}: {c[0]}/{c[1]}/{c[2]}  {' '.join(det)}")
