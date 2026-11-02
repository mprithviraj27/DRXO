import sys; sys.path.insert(0, "/home/claude/rw20")
from analyse_rw import *
NAMES = {15: "Speed reducer", 16: "Industrial refrigeration system", 17: "Tension/compression spring, case 1", 18: "Pressure vessel",
         19: "Welded beam", 20: "Three-bar truss", 21: "Multiple-disk clutch brake", 22: "Planetary gear train", 23: "Step-cone pulley",
         24: "Robot gripper", 25: "Hydrostatic thrust bearing", 26: "Four-stage gearbox", 27: "10-bar truss", 28: "Rolling-element bearing",
         29: "Gas transmission compressor", 30: "Tension/compression spring, case 2", 31: "Gear train", 32: "Himmelblau's function", 33: "Topology optimization"}
from rwc import D_, GN, HN
def g(x):
    if x >= C: return "infeasible"
    return f"{x:.6g}"
L = ["| Problem | D | f* (reported) | DRXO best | DRXO median | DRX-off median | Best entrant median (algorithm) |", "|---|---|---|---|---|---|---|"]
for k in PROBS:
    kd = keys("DRXO", k); ko = keys("DRX-off", k)
    ent = {a: np.median(keys(a, k)) for a in SRC}; ba = min(ent, key=ent.get)
    L.append(f"| RC{k} {NAMES[k]} | {D_[k-1]} | {FSTAR[k]:.6g} | {g(kd.min())} | {g(np.median(kd))} | {g(np.median(ko))} | {g(ent[ba])} ({ba}) |")
open("/home/claude/ms/tables/rw.md", "w").write("\n".join(L) + "\n")
ok = ALGS
med = np.array([[np.median(keys(a, k)) for a in ok] for k in PROBS]); rk = np.array([rankdata(m) for m in med]).mean(0)
L2 = ["| Measure | DRX-off | SASS | sCMAgES | COLSHADE | EnMODE | BiPop-εMAgES |", "|---|---|---|---|---|---|---|"]
row = []
for b in ok[1:]:
    c = [0, 0, 0]
    for k in PROBS: c[wtl(keys("DRXO", k), keys(b, k))] += 1
    row.append(f"{c[0]} / {c[1]} / {c[2]}")
L2.append("| DRXO vs (better / tie / worse) | " + " | ".join(row) + " |")
L2.append("| Mean rank (DRXO " + f"{rk[0]:.2f}) | " + " | ".join(f"{r:.2f}" for r in rk[1:]) + " |")
open("/home/claude/ms/tables/rw_summary.md", "w").write("\n".join(L2) + "\n")
print(open("/home/claude/ms/tables/rw_summary.md").read())
