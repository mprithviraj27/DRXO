"""Frozen DRXO and its DRX-off skeleton on CEC 2020 RC15-RC33 (25 runs, official MaxFEs). Resumable.
Rows: [problem, alg, run, best feasible f (nan if none), violation of best, fe]"""
import sys, os, json, time, warnings, hashlib; warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/claude/rw20"); sys.path.insert(0, "/home/claude/uni")
from multiprocessing import Pool
OUT = "/home/claude/rw20/res_rw20.json"; ALGS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["DRXO", "DRX-off"]; PROBS = list(range(15, 34)); RUNS = 25
def job(a):
    from rwc import make, maxfes, decode
    from drxo_u_final import drxo_u_final, drxo_u_skeleton
    k, alg, s = a; key, raw, lb, ub, D = make(k)
    if alg == "SR-skel":
        from drx_so import drx_so; from variants import V
        r = drx_so(key, lb, ub, D, max_fe=maxfes(D), seed=20000 + 100 * k + s, **V["srskel"])
    else:
        fn = drxo_u_final if alg == "DRXO" else drxo_u_skeleton
        r = fn(key, lb, ub, D, max_fe=maxfes(D), seed=20000 + 100 * k + s)
    f, v = decode(r[1]); return [k, alg, s, (None if f != f else float(f)), float(v), r[2]["fe"]]
if __name__ == "__main__":
    h = hashlib.sha256(open("/home/claude/uni/drxo_u_final.py", "rb").read()).hexdigest()
    assert h == open("/home/claude/uni/drxo_u_final.sha256").read().split()[0]
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r[0], r[1], r[2]) for r in out}
    jobs = [(k, a, s) for s in range(RUNS) for k in PROBS for a in ALGS if (k, a, s) not in done]
    print(len(jobs), "jobs", flush=True); t0 = time.time()
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            out.append(r)
            if i % 10 == 0 or i == len(jobs) - 1:
                json.dump(out, open(OUT + ".tmp", "w")); os.replace(OUT + ".tmp", OUT); print(f"{i+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    print("FINISHED")
