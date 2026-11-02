"""Instrumented runs (identical search to the frozen code): nucleation activity, replacements and diversity.
DRXO on all CEC 2014 functions (2 runs, D = 10, 30, 50, held-out seeds); DRX-off and SR-skel on the same runs for diversity."""
import sys, json, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni')
from multiprocessing import Pool
import numpy as np
OUT = "/home/claude/uni/instr.json"
def job(a):
    from cec14off import make; from drx_so import drx_so; from variants import V
    D, fid, s, var = a; f, b = make(fid, D); st = {}
    r = drx_so(f, -100, 100, D, max_fe=10000 * D, seed=14000 + 100 * fid + s, stats=st, **V[var])
    L = st["log"]; n = len(L)
    pick = lambda q: L[min(n - 1, int(q * n))]
    return dict(D=D, fid=fid, s=s, var=var, err=max(r[1] - b, 0), passes=n, nuc=st.get("nuc_events", 0), grew=st.get("nuc_grew", 0),
                dead=st.get("nuc_dead", 0), acc_rate=float(np.mean([e["n_acc"] / e["N"] for e in L])),
                SR=float(np.mean([e["SR"] for e in L])), div25=pick(0.25)["div"], div50=pick(0.5)["div"], div75=pick(0.75)["div"])
if __name__ == "__main__":
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r["D"], r["fid"], r["s"], r["var"]) for r in out}
    jobs = [(D, f, s, v) for D in (10, 30, 50) for s in (0, 1) for f in range(1, 31) for v in ("C1sm", "skel2", "srskel") if (D, f, s, v) not in done]
    with Pool(1) as p:
        for r in p.imap_unordered(job, jobs):
            out.append(r); json.dump(out, open(OUT, "w"))
    print("FINISHED")
