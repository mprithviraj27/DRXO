"""Generate all manuscript tables (markdown) from the raw results. Usage: python build_tables.py"""
import sys, os, re, json, glob
import numpy as np
from scipy.stats import ranksums, friedmanchisquare, rankdata
sys.path.insert(0, "/home/claude/uni"); sys.path.insert(0, "/home/claude")
T = "/home/claude/ms/tables"; H = "/home/claude/heldout14"
fl = lambda e: 0.0 if e <= 1e-8 else float(e)

def wtl(a, b):
    a, b = np.array(a), np.array(b)
    if np.allclose(np.median(a), np.median(b), atol=1e-8) and np.allclose(a.mean(), b.mean(), rtol=1e-3, atol=1e-8): return 1
    pv = ranksums(a, b).pvalue
    if pv < 0.05 and np.median(a) <= np.median(b) and a.mean() < b.mean(): return 0
    if pv < 0.05 and np.median(a) >= np.median(b) and a.mean() > b.mean(): return 2
    return 1

def fmt(x):
    if x == 0: return "0"
    return f"{x:.2e}".replace("e+0", "e+").replace("e-0", "e−").replace("e+", "e").replace("-", "−")

# ------------------------------------------------------------------ held-out CEC 2014
def load14(D):
    R = {}
    p = f"{H}/py_d{D}.json"
    if os.path.exists(p):
        for fid, a, s, e in json.load(open(p)): R.setdefault({"DRXO-U": "DRXO", "DRX-off": "DRX-off", "SR-skel": "SR-skel"}[a], {}).setdefault(fid, []).append(fl(e))
    for name, folder in [("L-SRTDE", f"lsrtde_d{D}"), ("RDE", f"rde_d{D}")]:
        for f in range(1, 31):
            q = f"{H}/{folder}/FINAL_F{f}_D{D}.txt"
            if os.path.exists(q): R.setdefault(name, {})[f] = [fl(float(l.split()[0])) for l in open(q) if l.strip()]
    for f in range(1, 31):
        q = f"{H}/jso_d{D}/F{f}.txt"
        if os.path.exists(q): R.setdefault("jSO", {})[f] = [fl(float(m)) for m in re.findall(r"error value = ([-+0-9.eE]+)", open(q).read())]
        q = f"/home/claude/cec14_official/L-SHADE/L-SHADE_{f}_{D}.txt"
        if os.path.exists(q): R.setdefault("L-SHADE", {})[f] = [fl(float(v)) for v in [l.split() for l in open(q) if l.strip()][-1]]
    return R

ORDER14 = ["DRXO", "DRX-off", "L-SHADE", "jSO", "RDE", "L-SRTDE"]
GROUPS14 = [("Unimodal (3)", range(1, 4)), ("Multimodal (13)", range(4, 17)), ("Hybrid (6)", range(17, 23)), ("Composition (8)", range(23, 31))]

def complete(R, a, fids, n=25): return a in R and all(len(R[a].get(f, [])) >= n for f in fids)

def heldout_tables():
    dims = [D for D in (10, 30, 50) if all(complete(load14(D), a, range(1, 31)) for a in ORDER14)]
    ranks, wt = {}, {}
    for D in dims:
        R = load14(D); F = range(1, 31)
        means = np.array([[np.mean(R[a][f]) for a in ORDER14] for f in F])
        rk = np.array([rankdata(m) for m in means]); ranks[D] = (rk.mean(0), friedmanchisquare(*rk.T).pvalue)
        for a in ORDER14[1:] + (["SR-skel"] if complete(R, "SR-skel", range(1, 31)) else []):
            row = []
            for g, fs in GROUPS14:
                c = [0, 0, 0]
                for f in fs: c[wtl(R["DRXO"][f], R[a][f])] += 1
                row.append(c)
            tot = [sum(x[i] for x in row) for i in range(3)]
            wt[(a, D)] = row + [tot]
        # appendix: per-function mean (std)
        L = [f"**Table A{dims.index(D)+1}.** Mean (standard deviation) of the final error on CEC 2014, D = {D} (25 runs; L-SHADE: 51 official runs). Bold marks the lowest mean.", "",
             "| F | " + " | ".join(ORDER14) + " |", "|---|" + "---|" * len(ORDER14)]
        for f in F:
            m = [np.mean(R[a][f]) for a in ORDER14]; best = min(m)
            cells = []
            for a, mm in zip(ORDER14, m):
                c = f"{fmt(mm)} ({fmt(np.std(R[a][f]))})"
                cells.append(f"**{c}**" if np.isclose(mm, best, rtol=1e-9, atol=1e-12) else c)
            L.append(f"| F{f} | " + " | ".join(cells) + " |")
        open(f"{T}/A_cec14_d{D}.md", "w").write("\n".join(L) + "\n")
    # Table: mean ranks
    L = ["| Algorithm | " + " | ".join(f"D = {D}" for D in dims) + " |", "|---|" + "---|" * len(dims)]
    for j, a in enumerate(ORDER14):
        cells = []
        for D in dims:
            v = ranks[D][0][j]; best = ranks[D][0].min()
            cells.append(f"**{v:.2f}**" if np.isclose(v, best) else f"{v:.2f}")
        L.append(f"| {a} | " + " | ".join(cells) + " |")
    L.append("| Friedman p | " + " | ".join(f"{ranks[D][1]:.1e}".replace("e-0", " × 10^−").replace("e-", " × 10^−") + "^" for D in dims) + " |")
    open(f"{T}/heldout_ranks.md", "w").write("\n".join(L) + "\n")
    L = ["| Comparison | D | " + " | ".join(g for g, _ in GROUPS14) + " | Total (30) |", "|---|---|" + "---|" * 5]
    for a in ORDER14[1:] + ["SR-skel"]:
        for D in dims:
            if (a, D) not in wt: continue
            L.append(f"| vs {a} | {D} | " + " | ".join(f"{c[0]} / {c[1]} / {c[2]}" for c in wt[(a, D)]) + " |")
    open(f"{T}/heldout_wtl.md", "w").write("\n".join(L) + "\n")
    json.dump({"dims": dims, "ranks": {D: list(ranks[D][0]) for D in dims}, "p": {D: ranks[D][1] for D in dims},
               "wtl": {f"{a}|{D}": wt[(a, D)][-1] for (a, D) in wt}}, open(f"{T}/heldout_summary.json", "w"), indent=1)
    print("held-out dims:", dims)

# ------------------------------------------------------------------ design suites
from cmp import load as load_design
F17 = [k for k in range(1, 31) if k != 2]
def runs_ok(R, a, fids, n): return a in R and all(len(R[a].get(f, [])) >= n for f in fids)

def cec17_table():
    L = ["| D | vs DRX-off | vs L-SHADE | vs jSO | vs RDE | vs L-SRTDE | Algorithms ranked | Mean rank of DRXO | Best mean rank |", "|---|---|---|---|---|---|---|---|---|"]
    info = {}
    for D in (10, 30, 50):
        R = load_design("cec17", D); R["DRXO"] = R.get("C1sm", {}); R["DRX-off"] = R.get("skel2", {})
        n = 25 if D != 10 else 25
        if not runs_ok(R, "DRXO", F17, n): print("cec17 D", D, "DRXO incomplete"); continue
        cells = []; algs = ["DRXO"]
        for a in ["DRX-off", "L-SHADE", "jSO", "RDE", "L-SRTDE"]:
            if not runs_ok(R, a, F17, 25): cells.append("—"); continue
            c = [0, 0, 0]
            for f in F17: c[wtl(R["DRXO"][f], R[a][f])] += 1
            cells.append(f"{c[0]} / {c[1]} / {c[2]}"); algs.append(a); info[(D, a)] = c
        means = np.array([[np.mean(R[a][f]) for a in algs] for f in F17]); rk = np.array([rankdata(m) for m in means]).mean(0)
        L.append(f"| {D} | " + " | ".join(cells) + f" | {len(algs)} | {rk[0]:.2f} | {rk.min():.2f} ({algs[int(np.argmin(rk))]}) |")
    open(f"{T}/cec17.md", "w").write("\n".join(L) + "\n"); return info

def cec22_table():
    import openpyxl, re as _re
    BIAS = {1: 300, 2: 400, 3: 600, 4: 800, 5: 900, 6: 1800, 7: 2000, 8: 2200, 9: 2300, 10: 2400, 11: 2600, 12: 2700}
    COMP = {"EA4eig": "EA4eigN100_10_{D}.xlsx", "NL-SHADE-LBC": "NL-SHADE-LBC_{D}.xlsx", "NL-SHADE-RSP-MID": "NL-SHADE-RSP-MID_{D}.xlsx",
            "S-LSHADE-DP": "S_LSHADE_DP_{D}.xlsx", "jSObinexpEig": "jSObinexpEig_{D}.xlsx"}
    res = {}
    for D in (10, 20):
        R = load_design("cec22", D); R["DRXO"] = R.get("C1sm", {})
        if not runs_ok(R, "DRXO", range(1, 13), 30): print("cec22 D", D, "incomplete"); continue
        for a, fn in COMP.items():
            ws = openpyxl.load_workbook(f"/home/claude/cec22_comp/{fn.format(D=D)}", read_only=True)["Sheet1"]
            rows = [r[:12] for r in ws.iter_rows(values_only=True)]
            R[a] = {f: [fl(r[f - 1]) for r in rows if r[f - 1] is not None] for f in range(1, 13)}
        algs = ["DRXO", "L-SHADE", "L-SRTDE"] + list(COMP)
        means = np.array([[np.mean(R[a][f]) for a in algs] for f in range(1, 13)]); rk = np.array([rankdata(m) for m in means])
        wt = {}
        for a in algs[1:]:
            c = [0, 0, 0]
            for f in range(1, 13): c[wtl(R["DRXO"][f], R[a][f])] += 1
            wt[a] = c
        res[D] = (algs, rk.mean(0), friedmanchisquare(*rk.T).pvalue, wt)
    if not res: return None
    Ds = sorted(res); algs = res[Ds[0]][0]
    L = ["| Algorithm | " + " | ".join(f"Mean rank, D = {D} | DRXO better / tie / worse, D = {D}" for D in Ds) + " |", "|---|" + "---|---|" * len(Ds)]
    for j, a in enumerate(algs):
        cells = []
        for D in Ds:
            rk = res[D][1]; v = rk[j]
            cells.append(f"**{v:.2f}**" if np.isclose(v, rk.min()) else f"{v:.2f}")
            cells.append("—" if a == "DRXO" else "{} / {} / {}".format(*res[D][3][a]))
        L.append(f"| {a} | " + " | ".join(cells) + " |")
    open(f"{T}/cec22_friedman.json", "w").write(json.dumps({D: res[D][2] for D in Ds}))
    open(f"{T}/cec22.md", "w").write("\n".join(L) + "\n"); return res

ABL = [("abl_notemp", "Thermal state and Arrhenius mobility (SHADE memory instead)"), ("abl_parent", "Grain consumption (parent replacement instead)"),
       ("abl_nodrive", "Stored-energy driving force"), ("abl_nonuc", "Strain-triggered nucleation"), ("abl_noZ", "Zener–Hollomon coupling of nucleation")]
SETS = [("cec17", 10, 25), ("cec22", 10, 30), ("cec17", 30, 10), ("cec17", 50, 10)]
def ablation_table():
    hdr = ["CEC 2017 D10", "CEC 2022 D10", "CEC 2017 D30", "CEC 2017 D50", "Total"]
    L = ["| Mechanism removed | " + " | ".join(hdr) + " |", "|---|" + "---|" * len(hdr)]; out = {}
    for v, name in ABL + [("skel2", "All DRX mechanisms (DRX-off skeleton)")]:
        cells = []; tot = [0, 0]; det = []
        for suite, D, n in SETS:
            R = load_design(suite, D); F = F17 if suite == "cec17" else list(range(1, 13))
            if not (runs_ok(R, v, F, n) and runs_ok(R, "C1sm", F, n)): cells.append("—"); continue
            worse = better = 0
            for f in F:
                a = R[v][f][:n] if v != "skel2" else R[v][f][:n]; b = R["C1sm"][f][:n]
                r = wtl(a, b)
                if r == 2: worse += 1
                elif r == 0: better += 1
            cells.append(f"{worse} / {better}"); tot[0] += worse; tot[1] += better
        cells.append(f"**{tot[0]} / {tot[1]}**"); out[v] = tot
        L.append(f"| {name} | " + " | ".join(cells) + " |")
    open(f"{T}/ablation.md", "w").write("\n".join(L) + "\n"); return out

SENS = [("sens_eps3", "Critical strain ε~c0~ = 3"), ("sens_eps9", "Critical strain ε~c0~ = 9"), ("sens_mu055", "Mobility μ~max~ = 0.55"),
        ("sens_mu075", "Mobility μ~max~ = 0.75"), ("sens_heat03", "Thermal inertia c~h~ = 0.3")]
def sensitivity_table():
    L = ["| Setting | D = 10 | D = 30 |", "|---|---|---|"]; out = {}
    for v, name in SENS:
        cells = []
        for D, n in ((10, 25), (30, 10)):
            R = load_design("cec17", D)
            if not (runs_ok(R, v, F17, n) and runs_ok(R, "C1sm", F17, n)): cells.append("—"); continue
            c = [0, 0, 0]
            for f in F17: c[wtl(R[v][f][:n], R["C1sm"][f][:n])] += 1
            cells.append(f"{c[0]} / {c[1]} / {c[2]}"); out[(v, D)] = c
        L.append(f"| {name} | " + " | ".join(cells) + " |")
    open(f"{T}/sensitivity.md", "w").write("\n".join(L) + "\n"); return out

def complexity_table():
    p = "/home/claude/uni/complexity.json"
    if not os.path.exists(p): return None
    c = json.load(open(p))
    L = ["| Algorithm | D | T0 | T1 | T2 | (T2 − T1)/T0 |", "|---|---|---|---|---|---|"]
    for name in ("DRXO", "DRX-off"):
        for D in (10, 30, 50):
            r = c[f"{name}_D{D}"]; L.append(f"| {name} | {D} | {c['T0']:.3f} | {r['T1']:.2f} | {r['T2']:.2f} | {r['ratio']:.1f} |")
    open(f"{T}/complexity.md", "w").write("\n".join(L) + "\n"); return c

def convergence_figure():
    p = f"{H}/traces_d50.json"
    if not os.path.exists(p): return None
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    tr = json.load(open(p)); CP = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    fids = sorted({r[0] for r in tr}); os.makedirs("/home/claude/ms/figs", exist_ok=True)
    fig, axs = plt.subplots(2, 3, figsize=(10, 6.2)); axs = axs.ravel()
    for ax, f in zip(axs, fids):
        for alg, col, ls in (("DRXO-U", "#1f4e79", "-"), ("DRX-off", "#9e9e9e", "--")):
            H_ = np.array([r[3] for r in tr if r[0] == f and r[1] == alg]); med = np.median(H_, 0)
            ax.plot(np.array(CP) * 500000, np.maximum(med, 1e-8), ls, color=col, lw=1.8, label="DRXO" if alg == "DRXO-U" else "DRX-off")
        q = f"/home/claude/cec14_official/L-SHADE/L-SHADE_{f}_50.txt"
        if os.path.exists(q):
            rows = np.array([[float(v) for v in l.split()] for l in open(q) if l.strip()])
            ax.plot(np.array(CP) * 500000, np.maximum(np.median(rows, 1), 1e-8), ":", color="#c55a11", lw=1.8, label="L-SHADE")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"F{f}", fontsize=10); ax.tick_params(labelsize=8); ax.grid(alpha=0.3)
        ax.set_xlabel("Function evaluations", fontsize=8); ax.set_ylabel("Median error", fontsize=8)
    axs[0].legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig("/home/claude/ms/figs/fig1_convergence.png", dpi=200); plt.close(fig)
    return fids

def appendix():
    parts = ["## Appendix A. Per-function results on the held-out CEC 2014 suite"]
    for D in (10, 30, 50):
        q = f"{T}/A_cec14_d{D}.md"
        if os.path.exists(q): parts.append(open(q).read())
    open(f"{T}/appendix.md", "w").write("\n\n".join(parts) + "\n")

if __name__ == "__main__":
    heldout_tables()
    print("cec17", cec17_table()); print("cec22", None if cec22_table() is None else "ok")
    print("ablation", ablation_table()); print("sens", sensitivity_table()); print("cplx", complexity_table() is not None)
    print("fig", convergence_figure()); appendix()
