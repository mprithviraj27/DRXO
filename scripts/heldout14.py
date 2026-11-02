"""HELD-OUT TEST on CEC 2014 (30 functions, 25 runs, 10000*D evaluations). Run only on FROZEN code.
Usage: python heldout14.py D alg1,alg2      (algs: DRXO-U = frozen engine, DRX-off = its skeleton)
Results: /home/claude/heldout14/py_d{D}.json rows [fid, alg, run, error]"""
import sys, os, json, time, warnings, hashlib
warnings.filterwarnings("ignore")
sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni')
from multiprocessing import Pool
D, ALGS = int(sys.argv[1]), sys.argv[2].split(",")
OUT = f"/home/claude/heldout14/py_d{D}.json"; os.makedirs("/home/claude/heldout14", exist_ok=True)
def job(a):
    from cec14off import make
    from drxo_u_final import drxo_u_final, drxo_u_skeleton
    fid, alg, s = a; f, b = make(fid, D)
    if alg == "SR-skel":   # L-SRTDE-type DRX-free baseline (added after the held-out test; DRXO code unchanged)
        from drx_so import drx_so; from variants import V
        e = drx_so(f, -100, 100, D, max_fe=10000 * D, seed=14000 + 100 * fid + s, **V["srskel"])[1] - b
    else:
        fn = drxo_u_final if alg == "DRXO-U" else drxo_u_skeleton
        e = fn(f, -100, 100, D, max_fe=10000 * D, seed=14000 + 100 * fid + s)[1] - b
    return fid, alg, s, (0.0 if e <= 1e-8 else float(e))
if __name__ == "__main__":
    h = hashlib.sha256(open("/home/claude/uni/drxo_u_final.py", "rb").read()).hexdigest()
    assert h == open("/home/claude/uni/drxo_u_final.sha256").read().split()[0], "code changed after freeze"
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    done = {(r[0], r[1], r[2]) for r in out}
    jobs = [(f, a, s) for s in range(25) for a in ALGS for f in range(1, 31) if (f, a, s) not in done]
    print(f"D{D}: {len(jobs)} jobs", flush=True); t0 = time.time()
    with Pool(2) as p:
        for i, r in enumerate(p.imap_unordered(job, jobs)):
            out.append(list(r))
            if i % 20 == 0 or i == len(jobs) - 1:
                json.dump(out, open(OUT + ".tmp", "w")); os.replace(OUT + ".tmp", OUT)
                print(f"{i+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    print("FINISHED", flush=True)
