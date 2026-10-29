"""Flexible-connection contract parameters implied by the results (pooled medians over all MV buses, OLTC base case).

Columns: energy budget alpha = 0.25 % and 1 %, hour budget 100 h/yr. Rows: gain over firm HC, curtailed hours,
events, median and longest event, deepest curtailment, curtailed energy share. Output: tex/numbers_contract.tex
"""
import json
import numpy as np
from core import flexible_hc, curtailment_stats

RES = "../results"
import os
rows = []
HAS_N1 = all(os.path.exists(f"{RES}/{g}_oltc_N1.npz") for g in ["rural", "semiurb", "urban", "comm"])
for g in ["rural", "semiurb", "urban", "comm"]:
    PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    PN = np.load(f"{RES}/{g}_oltc_N1.npz")["PM"].astype(float) if HAS_N1 else None
    for i, pm in enumerate(PM):
        firm = pm.min(); s = np.sort(pm)
        cases = [("A", pm, flexible_hc(pm, 0.0025)), ("B", pm, flexible_hc(pm, 0.01)), ("H", pm, s[400])]
        if HAS_N1:
            cases.append(("N", PN[i], flexible_hc(PN[i], 0.01)))
        for tag, pm, P in cases:
            firm = pm.min()
            st = curtailment_stats(pm, P)
            rows.append(dict(tag=tag, firm=firm, add=P - firm, gain=(P / firm - 1) * 100, hrs=st["hours"], nev=st["n_events"],
                             med=st["med_event_h"], mx=st["max_event_h"], dep=st["max_depth"] * 100,
                             en=st["energy_frac"] * 100))
mac = {}
for tag in ("ABHN" if HAS_N1 else "ABH"):
    r = [x for x in rows if x["tag"] == tag]
    med = {k: float(np.median([x[k] for x in r])) for k in r[0] if k not in ("tag",)}
    mac.update({f"ct{tag}add": f"{med['add']:.2f}", f"ct{tag}firm": f"{float(np.median([x['firm'] for x in r])):.2f}", f"ct{tag}gain": f"{med['gain']:.1f}", f"ct{tag}hrs": f"{med['hrs']:.0f}", f"ct{tag}nev": f"{med['nev']:.0f}",
                f"ct{tag}med": f"{med['med'] * 60:.0f}", f"ct{tag}mx": f"{med['mx'] + 1e-9:.1f}", f"ct{tag}dep": f"{med['dep']:.0f}",
                f"ct{tag}en": f"{med['en']:.2f}"})
# hour budget matching the bulk-system study's curtailment at alpha = 0.25 % (about 85 h/yr)
g85 = []
for g in ["rural", "semiurb", "urban", "comm"]:
    for pm in np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float):
        g85.append(np.sort(pm)[340] / pm.min() - 1)
mac["gEightyFive"] = f"{np.median(g85) * 100:.1f}"
print(json.dumps(mac, indent=1))
with open("../tex/numbers_contract.tex", "w") as f:
    for k, v in sorted(mac.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

fc = {}
try:
    for line in open("../tex/numbers_fc.tex"):
        k, v = line.strip()[len("\\newcommand{\\"):].split("}{")
        fc[k] = v.rstrip("}")
except FileNotFoundError:
    pass
R = [("Firm HC (MW)", "firm"), ("Added capacity (MW)", "add"), ("Gain over firm HC (\\%)", "gain"), ("Curtailed hours (h/yr)", "hrs"), ("Events (per year)", "nev"),
     ("Median event length (min)", "med"), ("Longest event (h)", "mx"), ("Deepest curtailment (\\% of $P$)", "dep"),
     ("Curtailed energy (\\%)", "en")]
with open("../tex/tab_contract.tex", "w") as f:
    f.write("\\begin{table}[t]\n\\centering\n\\caption{Contract Parameters Implied by Each Budget (Pooled Medians, Perfect Foresight)}\n\\label{tab:contract}\n")
    cols = "ABHN" if HAS_N1 else "ABH"
    f.write("\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{l" + "c" * len(cols) + "}\n\\toprule\n")
    if HAS_N1:
        f.write(" & \\multicolumn{3}{c}{N-0} & N-1 \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-5}\n")
    f.write("Parameter & $\\alpha=0.25\\%$ & $\\alpha=1\\%$ & 100 h" + (" & $\\alpha=1\\%$" if HAS_N1 else "") + " \\\\\n\\midrule\n")
    for lab, k in R:
        f.write(lab + "".join(f" & {mac['ct' + t + k]}" for t in cols) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\vspace{2pt}\n\\begin{minipage}{\\columnwidth}\\footnotesize ")
    f.write("Medians over buses of each statistic, so they need not multiply. N-1: planning band after outages. All values assume the headroom is known in each 15-min interval (forecasts: Section~\\ref{sec:sens}).")
    f.write("\\end{minipage}\n\\end{table}\n")
