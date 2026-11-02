"""Post-hoc: rank-based stored energy (post_rankE) vs DRXO, DRX-off and SR-skel on CEC 2020 RC15-RC33."""
import sys, json, numpy as np
sys.path.insert(0, "/home/claude/rw20")
import analyse_rw as A
for k, a, s, f, v, *_ in json.load(open("/home/claude/uni/posthoc_rw.json")):
    A.R.setdefault("rankE", {}).setdefault(k, []).append((np.nan if f is None else f, v))
assert all(len(A.R["rankE"][k]) == 25 for k in A.PROBS)
for b in ["DRXO", "DRX-off", "SR-skel", "SASS", "COLSHADE", "EnMODE", "BiPop-εMAgES", "sCMAgES"]:
    c = [0, 0, 0]; det = []
    for k in A.PROBS:
        r = A.wtl(A.keys("rankE", k), A.keys(b, k)); c[r] += 1
        if r != 1: det.append(("+" if r == 0 else "-") + f"RC{k}")
    print(f"rankE vs {b}: {c[0]}/{c[1]}/{c[2]}  {' '.join(det)}")
fr = {k: np.mean([v <= 1e-12 for _, v in A.R["rankE"][k]]) for k in A.PROBS}; print("feasible rate min", min(fr.values()))
for k in (17, 22, 24, 25, 30):
    print(k, "median DRXO %.6g rankE %.6g off %.6g" % (np.median(A.keys("DRXO", k)), np.median(A.keys("rankE", k)), np.median(A.keys("DRX-off", k))))
