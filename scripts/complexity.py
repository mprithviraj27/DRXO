"""Algorithm complexity, CEC 2014/2017 procedure: T0 = fixed arithmetic loop; T1 = 200,000 evaluations of F18;
T2 = mean time of 5 runs of the algorithm on F18 with 200,000 evaluations. Usage: python complexity.py"""
import sys, time, json, math; sys.path.insert(0, '/home/claude'); sys.path.insert(0, '/home/claude/uni')
import numpy as np
from cec17off import make
from drxo_u_final import drxo_u_final, drxo_u_skeleton
np.seterr(all='ignore'); t = time.perf_counter(); x = np.float64(0.55)
for i in range(1000000):
    x = x + x; x = x / 2; x = x * x; x = np.sqrt(x); x = np.log(x); x = np.exp(x); x = x / (x + 2)
T0 = time.perf_counter() - t; res = {"T0": T0}
for D in (10, 30, 50):
    f, b = make(18, D); X = np.random.default_rng(0).uniform(-100, 100, (200000, D))
    t = time.perf_counter()
    for xx in X: f(xx)
    T1 = time.perf_counter() - t
    for name, alg in (("DRXO", drxo_u_final), ("DRX-off", drxo_u_skeleton)):
        ts = []
        for s in range(5):
            t = time.perf_counter(); alg(f, -100, 100, D, max_fe=200000, seed=s); ts.append(time.perf_counter() - t)
        res[f"{name}_D{D}"] = dict(T1=T1, T2=float(np.mean(ts)), ratio=(np.mean(ts) - T1) / T0)
        print(name, D, res[f"{name}_D{D}"], flush=True)
json.dump(res, open("/home/claude/uni/complexity.json", "w"), indent=1); print("FINISHED")
