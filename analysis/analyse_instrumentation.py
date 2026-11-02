import json, numpy as np
d = json.load(open("/home/claude/uni/instr.json"))
print("n rows", len(d))
for D in (10, 30, 50):
    for v in ("C1sm", "skel2", "srskel"):
        R = [r for r in d if r["D"] == D and r["var"] == v]
        if not R: continue
        nuc = np.array([r["nuc"] for r in R]); grew = sum(r["grew"] for r in R); dead = sum(r["dead"] for r in R)
        print(f"D{D} {v:7s} n={len(R):3d} runs_with_nuc={np.mean(nuc>0):.2f} fns_with_nuc={len({r['fid'] for r in R if r['nuc']>0})} nuc/run mean={nuc.mean():.2f} max={nuc.max()} "
              f"grew={grew} dead={dead} acc={np.mean([r['acc_rate'] for r in R]):.3f} SR={np.mean([r['SR'] for r in R]):.3f} "
              f"div25={np.median([r['div25'] for r in R]):.3g} div50={np.median([r['div50'] for r in R]):.3g} div75={np.median([r['div75'] for r in R]):.3g}")
    R = [r for r in d if r["D"] == D and r["var"] == "C1sm" and r["nuc"] > 0]
    print("   functions with nucleation:", sorted({r["fid"] for r in R}), " events per pass x1e3:", [round(1e3 * r["nuc"] / r["passes"], 2) for r in R][:12])
