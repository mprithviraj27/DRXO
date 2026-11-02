"""Held-out CEC 2014 analysis. Usage: python analyse14.py D
Sources: DRXO-U and DRX-off (our runs, 25 each), L-SRTDE and RDE (authors' C++ code, 25 runs, CEC 2014 functions),
jSO (official CEC 2017 code linked to CEC 2014, 25 runs), L-SHADE (official CEC 2014 competition results, 51 runs)."""
import sys, os, re, json
import numpy as np
from scipy.stats import ranksums, friedmanchisquare, rankdata
D = int(sys.argv[1]); FIDS = list(range(1, 31)); H = "/home/claude/heldout14"
fl = lambda e: 0.0 if e <= 1e-8 else float(e)
R = {}
p = f"{H}/py_d{D}.json"
if os.path.exists(p):
    for fid, a, s, e in json.load(open(p)): R.setdefault(a, {}).setdefault(fid, []).append(fl(e))
for name, folder in [("L-SRTDE", f"lsrtde_d{D}"), ("RDE", f"rde_d{D}")]:
    for f in FIDS:
        q = f"{H}/{folder}/FINAL_F{f}_D{D}.txt"
        if os.path.exists(q): R.setdefault(name, {})[f] = [fl(float(l.split()[0])) for l in open(q) if l.strip()]
for f in FIDS:
    q = f"{H}/jso_d{D}/F{f}.txt"
    if os.path.exists(q):
        R.setdefault("jSO", {})[f] = [fl(float(m)) for m in re.findall(r"error value = ([-+0-9.eE]+)", open(q).read())]
    q = f"/home/claude/cec14_official/L-SHADE/L-SHADE_{f}_{D}.txt"
    if os.path.exists(q):
        rows = [l.split() for l in open(q) if l.strip()]; R.setdefault("L-SHADE", {})[f] = [fl(float(v)) for v in rows[-1]]
ORDER = ["DRXO-U", "DRX-off", "L-SHADE", "jSO", "RDE", "L-SRTDE"]
need = {"L-SHADE": 51}
algs = [a for a in ORDER if a in R and all(len(R[a].get(f, [])) >= need.get(a, 25) for f in FIDS)]
miss = [a for a in ORDER if a not in algs]
lines = [f"CEC 2014 (held-out), D = {D}. Mean error (median in brackets).", "missing/incomplete: " + ", ".join(
    f"{a} ({sum(len(v) for v in R.get(a, {}).values())})" for a in miss)]
means = np.array([[np.mean(R[a][f]) for a in algs] for f in FIDS])
lines.append("F   " + "".join(f"{a:>13s}" for a in algs))
for i, f in enumerate(FIDS): lines.append(f"F{f:<3d}" + "".join(f"{means[i, j]:13.3e}" for j in range(len(algs))))
ranks = np.array([rankdata(m) for m in means]); lines.append("Rank" + "".join(f"{v:13.2f}" for v in ranks.mean(0)))
if len(algs) > 2: lines.append(f"Friedman p = {friedmanchisquare(*ranks.T).pvalue:.2e}")
def wtl(a, b):
    a, b = np.array(a), np.array(b)
    if np.allclose(np.median(a), np.median(b), atol=1e-8) and np.allclose(a.mean(), b.mean(), rtol=1e-3, atol=1e-8): return 1
    pv = ranksums(a, b).pvalue
    if pv < 0.05 and np.median(a) <= np.median(b) and a.mean() < b.mean(): return 0
    if pv < 0.05 and np.median(a) >= np.median(b) and a.mean() > b.mean(): return 2
    return 1
groups = {"unimodal F1-3": range(1, 4), "multimodal F4-16": range(4, 17), "hybrid F17-22": range(17, 23), "composition F23-30": range(23, 31)}
if "DRXO-U" in algs:
    lines.append("\nDRXO-U vs each (Wilcoxon rank-sum p < 0.05, median and mean agree): better / tie / worse")
    for a in algs[1:]:
        c = [0, 0, 0]; det = []; per = {g: [0, 0, 0] for g in groups}
        for f in FIDS:
            r = wtl(R["DRXO-U"][f], R[a][f]); c[r] += 1
            if r != 1: det.append(("+" if r == 0 else "-") + f"F{f}")
            for g, fs in groups.items():
                if f in fs: per[g][r] += 1
        lines.append(f"  vs {a:8s} {c[0]:2d} / {c[1]:2d} / {c[2]:2d}   {' '.join(det)}")
        lines.append("      by group: " + "; ".join(f"{g}: {v[0]}/{v[1]}/{v[2]}" for g, v in per.items()))
txt = "\n".join(lines); print(txt); open(f"{H}/result14_d{D}.txt", "w").write(txt)
