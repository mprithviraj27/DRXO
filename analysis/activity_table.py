import json, numpy as np
def sci(v):
    m, e = f'{v:.1e}'.split('e'); e = int(e)
    return m if e == 0 else f'{m} × 10^{"−" if e < 0 else ""}{abs(e)}^'
d = json.load(open("/home/claude/uni/instr.json"))
NM = {"C1sm": "DRXO", "skel2": "DRX-off", "srskel": "SR-skel"}
L = ["| D | Algorithm | Runs with nucleation | Nuclei per 1,000 passes | Nuclei that grew | Mean success rate | Mean acceptance rate | Diversity at 25 % / 50 % of budget |", "|---|---|---|---|---|---|---|---|"]
for D in (10, 30, 50):
    for v in ("C1sm", "srskel", "skel2"):
        R = [r for r in d if r["D"] == D and r["var"] == v]
        nuc = sum(r["nuc"] for r in R); ps = sum(r["passes"] for r in R); g = sum(r["grew"] for r in R); de = sum(r["dead"] for r in R)
        a = f"{np.mean([r['nuc'] > 0 for r in R]) * 100:.0f} %" if v == "C1sm" else "—"
        b = f"{1000 * nuc / ps:.1f}" if v == "C1sm" else "—"
        c = f"{100 * g / (g + de):.0f} %" if v == "C1sm" else "—"
        L.append(f"| {D} | {NM[v]} | {a} | {b} | {c} | {np.mean([r['SR'] for r in R]):.2f} | {np.mean([r['acc_rate'] for r in R]):.2f} | {sci(np.median([r['div25'] for r in R]))} / {sci(np.median([r['div50'] for r in R]))} |")
open("tables/activity.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))
