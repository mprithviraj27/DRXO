"""Convergence traces of the frozen DRXO and its DRX-off skeleton on the held-out CEC 2014 suite (D = 50),
with the held-out seeds; the error is recorded at the 14 CEC 2014 checkpoints (0.01 ... 1.0 of the budget).
Each traced run must end on exactly the value recorded in the held-out test (checked)."""
import sys, os, json, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni')
from multiprocessing import Pool
import numpy as np
D = 50; FIDS = [int(x) for x in sys.argv[1].split(",")]; RUNS = 25
CP = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
OUT = "/home/claude/heldout14/traces_d50.json"
def job(a):
    from cec14off import make
    from drxo_u_final import drx_so, DRXO_U, SKELETON
    fid, alg, s = a; f0, b = make(fid, D); best = [np.inf]; hist = []
    marks = [int(c * 10000 * D) for c in CP]; nfe = [0]
    def f(x):
        v = f0(x); nfe[0] += 1
        if v < best[0]: best[0] = v
        while marks and nfe[0] >= marks[0]: hist.append(max(best[0] - b, 0.0)); marks.pop(0)
        return v
    cfg = DRXO_U if alg == "DRXO-U" else SKELETON
    r = drx_so(f, -100, 100, D, max_fe=10000 * D, seed=14000 + 100 * fid + s, **cfg)
    e = r[1] - b
    return fid, alg, s, [0.0 if h <= 1e-8 else h for h in hist], (0.0 if e <= 1e-8 else float(e))
if __name__ == "__main__":
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r[0], r[1], r[2]) for r in out}
    jobs = [(f, a, s) for s in range(RUNS) for f in FIDS for a in ("DRXO-U", "DRX-off") if (f, a, s) not in done]
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            out.append(list(r)); json.dump(out, open(OUT, "w"))
    ref = {(r[0], r[1], r[2]): r[3] for r in json.load(open("/home/claude/heldout14/py_d50.json"))}
    bad = [r[:3] for r in out if abs(r[4] - ref[(r[0], r[1], r[2])]) > 1e-9 * max(1, ref[(r[0], r[1], r[2])])]
    print("traced runs:", len(out), "mismatches with held-out finals:", len(bad), bad[:5]); print("FINISHED")
