"""Run engine variants on CEC 2017 / CEC 2022 (design suites). Resumable.
Usage: python so_run.py SUITE D RUNS var1,var2,...      (SUITE = cec17 | cec22)
Results: res_SUITE_dD.json  rows [fid, variant, seed_index, error, fe]"""
import sys, os, json, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni')
from multiprocessing import Pool
SUITE, D, RUNS, VARS = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4].split(",")
assert SUITE in ("cec17", "cec22"), "held-out suites are not allowed here"
if SUITE == "cec17":
    from cec17off import make; FIDS = [k for k in range(1, 31) if k != 2]; FE = 10000 * D
    seed = lambda fid, s: 5000 + 100 * fid + s
else:
    from cec22off import make; FIDS = list(range(1, 13)); FE = {10: 200000, 20: 1000000}[D]
    seed = lambda fid, s: 1000 * fid + s
OUT = f"/home/claude/uni/res_{SUITE}_d{D}.json"

def job(a):
    from variants import V
    from drx_so import drx_so
    fid, var, s = a; f, b = make(fid, D)
    r = drx_so(f, -100, 100, D, max_fe=FE, seed=seed(fid, s), **V[var])
    e = r[1] - b
    return fid, var, s, (0.0 if e <= 1e-8 else float(e)), r[2]["fe"]

if __name__ == "__main__":
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r[0], r[1], r[2]) for r in out}
    jobs = [(fid, v, s) for s in range(RUNS) for v in VARS for fid in FIDS if (fid, v, s) not in done]
    jobs.sort(key=lambda j: (j[2], -j[0]))
    print(f"{SUITE} D{D}: {len(jobs)} jobs", flush=True); t0 = time.time()
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs, chunksize=1)):
            out.append(list(r))
            if i % 25 == 0 or i == len(jobs) - 1:
                json.dump(out, open(OUT + ".tmp", "w")); os.replace(OUT + ".tmp", OUT)
                print(f"{i+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    print("FINISHED", flush=True)
