"""Post-hoc exploratory diagnostics (after the held-out test). Same seeds as the held-out runs.
Usage: python posthoc.py cec14 D variant[,variant]  |  python posthoc.py rw variant[,variant]"""
import sys, os, json, time, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni'); sys.path.insert(0, '/home/claude/rw20')
from multiprocessing import Pool
SUITE = sys.argv[1]
if SUITE == "cec14": D = int(sys.argv[2]); ALGS = sys.argv[3].split(","); RUNS = int(sys.argv[4]) if len(sys.argv) > 4 else 25; OUT = f"/home/claude/uni/posthoc_c14_d{D}.json"
else: ALGS = sys.argv[2].split(","); RUNS = 25; OUT = "/home/claude/uni/posthoc_rw.json"
def job(a):
    from drx_so import drx_so; from variants import V
    k, alg, s = a
    if SUITE == "cec14":
        from cec14off import make; f, b = make(k, D)
        e = drx_so(f, -100, 100, D, max_fe=10000 * D, seed=14000 + 100 * k + s, **V[alg])[1] - b
        return [k, alg, s, (0.0 if e <= 1e-8 else float(e))]
    from rwc import make, maxfes, decode
    key, raw, lb, ub, Dd = make(k)
    r = drx_so(key, lb, ub, Dd, max_fe=maxfes(Dd), seed=20000 + 100 * k + s, **V[alg])
    f, v = decode(r[1]); return [k, alg, s, (None if f != f else float(f)), float(v)]
if __name__ == "__main__":
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r[0], r[1], r[2]) for r in out}
    P = range(1, 31) if SUITE == "cec14" else range(15, 34)
    jobs = [(k, a, s) for s in range(RUNS) for a in ALGS for k in P if (k, a, s) not in done]
    print(len(jobs), "jobs", flush=True); t0 = time.time()
    with Pool(int(os.environ.get("NP", "1"))) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            out.append(r)
            if i % 10 == 0 or i == len(jobs) - 1:
                json.dump(out, open(OUT + ".tmp", "w")); os.replace(OUT + ".tmp", OUT); print(f"{i+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    print("FINISHED", flush=True)
