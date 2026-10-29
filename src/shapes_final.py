"""Load-shape extension: table, macros. Uses stored OLTC base-case headroom series."""
import numpy as np, pandas as pd
from shapes import load_shapes, energy_hc, hour_hc, hours_curtailed, GRIDS, RES

SPEC = [("Flat", "Flat (base case)"), ("DC20@15", "$\\pm$10\\%, peak 15:00"), ("DC20@19", "$\\pm$10\\%, peak 19:00"),
        ("DC20@3", "$\\pm$10\\%, peak 03:00"), ("DC40@15", "$\\pm$20\\%, peak 15:00"), ("DC40@3", "$\\pm$20\\%, peak 03:00")]
shapes, t = load_shapes([s for s, _ in SPEC[1:]])
rows = []
track = []
for g in GRIDS:
    PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    for i in range(PM.shape[0]):
        pb = PM[i]
        track.append(pb.mean() / pb.min())
        for s, _ in SPEC:
            f = shapes[s]
            r = pb / f
            H0 = r.min()
            e1 = energy_hc(pb, f, 0.01)
            rows.append(dict(g=g, i=i, s=s, lf=f.mean(), firm=H0, e1=e1, h100=hour_hc(r, 100), hrs1=hours_curtailed(pb, f, e1)))
df = pd.DataFrame(rows)
flat = df[df.s == "Flat"].set_index(["g", "i"])
DIG = {str(i): w for i, w in enumerate(["Z", "O", "T", "Th", "Fo", "F", "Si", "Se", "E", "N"])}
mac = {}
lines = []
for s, label in SPEC:
    d = df[df.s == s].set_index(["g", "i"])
    lf = d.lf.iloc[0]
    fp = np.median(d.firm / flat.firm - 1) * 100
    fe = np.median(d.firm * lf / flat.firm - 1) * 100
    g1 = np.median(d.e1 / d.firm - 1) * 100
    hr = d.hrs1.median()
    sg = lambda x: (f"$-${abs(x):.1f}" if x < -0.05 else f"+{x:.1f}" if x > 0.05 else "0.0")
    lines.append(f"{label} & {lf:.2f} & {sg(fp)} & {sg(fe)} & {g1:.1f} & {hr:.0f} \\\\")
    tag = "".join(DIG.get(c, c) for c in s.replace("@", "at"))
    neg = lambda x, f="{:.1f}": ("$-$" + f.format(abs(x))) if x < 0 else f.format(x)
    mac[f"shp{tag}fp"] = neg(fp); mac[f"shp{tag}fe"] = neg(fe); mac[f"shp{tag}gO"] = f"{g1:.1f}"; mac[f"shp{tag}hrs"] = f"{hr:.0f}"
    print(label, round(lf, 2), round(fp, 1), round(fe, 1), round(g1, 1), hr)
tr = (np.median(track) - 1) * 100
mac["shpTrack"] = f"{tr:.1f}"
print("tracking bound median energy gain", tr)
with open("../tex/tab_shapes.tex", "w") as f:
    f.write("\\begin{table}[t]\n\\centering\n\\caption{Effect of a Daily Swing in the New Load (Pooled Medians)}\n\\label{tab:shapes}\n\\setlength{\\tabcolsep}{2.6pt}\n")
    f.write("\\begin{tabular}{lccccc}\n\\toprule\n")
    f.write("Daily swing & Load & \\multicolumn{2}{c}{Firm HC vs flat (\\%)} & Gain & Hours \\\\\n\\cmidrule(lr){3-4}\n")
    f.write(" & factor & Peak & Energy & at 1\\% (\\%) & at 1\\% \\\\\n\\midrule\n")
    f.write("\n".join(lines) + "\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\vspace{2pt}\n\\begin{minipage}{\\columnwidth}\\footnotesize Peak: firm nameplate $H_k^{\\mathrm{firm}}$. Energy: firm nameplate times load factor, i.e., annual energy that can be served without curtailment. Gain and hours refer to $\\alpha=1\\%$ of the load's own energy.\\end{minipage}\n\\end{table}\n")
with open("../tex/numbers_shapes.tex", "w") as f:
    for k, v in sorted(mac.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
