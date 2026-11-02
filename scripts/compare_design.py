"""Compare variants (and rivals) on design suites.
Usage: python cmp.py TARGET BASE1,BASE2,... [SUITE:D,...]
Rivals available as bases: L-SRTDE, RDE, L-SHADE (where files exist)."""
import sys, os, re, json, glob
import numpy as np
from scipy.stats import ranksums
sys.path.insert(0, '/home/claude')
fl = lambda e: 0.0 if e <= 1e-8 else float(e)
SETS = ["cec17:10", "cec22:10", "cec17:30", "cec22:20", "cec17:50"]

def fids(suite): return [k for k in range(1, 31) if k != 2] if suite == "cec17" else list(range(1, 13))

def load(suite, D):
    R = {}
    p = f"/home/claude/uni/res_{suite}_d{D}.json"
    if os.path.exists(p):
        for fid, v, s, e, *_ in json.load(open(p)): R.setdefault(v, {}).setdefault(fid, []).append(fl(e))
    def rd(name, path):
        for f in fids(suite):
            q = path.format(f=f, D=D)
            if os.path.exists(q): R.setdefault(name, {})[f] = [fl(float(l.split()[0])) for l in open(q) if l.strip()]
    if suite == "cec22":
        rd("L-SRTDE", "/home/claude/rivals_build/lsrtde/d{D}/FINAL_F{f}_D{D}.txt")
        from cec22off import BIAS
        q = f"/home/claude/official_lshade/off_lshade_d{D}.txt"; fid = None
        if os.path.exists(q):
            for line in open(q):
                m = re.search(r"Function = (\d+)", line)
                if m: fid = int(m.group(1))
                m = re.search(r"error value = ([-+0-9.eE]+)", line)
                if m: R.setdefault("L-SHADE", {}).setdefault(fid, []).append(fl(float(m.group(1)) - BIAS[fid]))
    else:
        rd("L-SRTDE", "/home/claude/drx_lsrtde/orig17/FINAL_F{f}_D{D}.txt" if D == 10 else "/home/claude/rivals_build/lsrtde17/d{D}/FINAL_F{f}_D{D}.txt")
        rd("RDE", "/home/claude/rivals_build/rde/d{D}/FINAL_F{f}_D{D}.txt")
        for f in fids(suite):   # jSO official CEC 2017 results (51 runs; last row = final error)
            q = f"/home/claude/rivals14/jso_src/jSO-SOURCE-RESULTS/rawResults/jSO_{f}_{D}.txt"
            if os.path.exists(q):
                rows = [l.split() for l in open(q) if l.strip()]
                R.setdefault("jSO", {})[f] = [fl(float(v)) for v in rows[-1]]
        paths = [f"/home/claude/rivals_build/lshade17/d{D}/out.txt"] + sorted(glob.glob(f"/home/claude/rivals_build/lshade17/d{D}/out_f*.txt"))
        fid = None
        for line in (l for q in paths if os.path.exists(q) for l in open(q)):
            m = re.search(r"Function = (\d+)", line)
            if m: fid = int(m.group(1))
            m = re.search(r"error value = ([-+0-9.eE]+)", line)
            if m: R.setdefault("L-SHADE", {}).setdefault(fid, []).append(fl(float(m.group(1)) - 100 * fid))
    return R

def wtl(a, b):
    a, b = np.array(a), np.array(b)
    if np.allclose(np.median(a), np.median(b), atol=1e-8) and np.allclose(a.mean(), b.mean(), rtol=1e-3, atol=1e-8): return 1
    pv = ranksums(a, b).pvalue
    if pv < 0.05 and np.median(a) <= np.median(b) and a.mean() < b.mean(): return 0
    if pv < 0.05 and np.median(a) >= np.median(b) and a.mean() > b.mean(): return 2
    return 1

def compare(T, bases, sets=SETS, quiet=False):
    tot = {b: [0, 0, 0] for b in bases}; lines = []
    for sd in sets:
        suite, D = sd.split(":"); D = int(D); R = load(suite, D)
        if T not in R: continue
        fs = [f for f in fids(suite) if f in R[T]]
        if len(fs) < len(fids(suite)): continue
        n = min(len(R[T][f]) for f in fs)
        row = f"{suite} D{D:<3d} (n={n:2d})"
        for b in bases:
            if b not in R or any(f not in R[b] for f in fs): row += f" | {b}: --"; continue
            c = [0, 0, 0]; det = []
            for f in fs:
                r = wtl(R[T][f], R[b][f]); c[r] += 1
                if r != 1: det.append(("+" if r == 0 else "-") + str(f))
            for i in range(3): tot[b][i] += c[i]
            row += f" | {b}: {c[0]}/{c[1]}/{c[2]}"
            if not quiet: row += " [" + " ".join(det) + "]"
        lines.append(row)
    lines.append("TOTAL " + " | ".join(f"{b}: {v[0]}/{v[1]}/{v[2]}" for b, v in tot.items()))
    return "\n".join(lines), tot

if __name__ == "__main__":
    T = sys.argv[1]; bases = sys.argv[2].split(",")
    sets = sys.argv[3].split(",") if len(sys.argv) > 3 else SETS
    print(f"== {T} vs {bases}"); print(compare(T, bases, sets)[0])
