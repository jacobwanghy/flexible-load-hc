"""Pooled metrics for the three network configurations -> tab_sens.tex and macros."""
import os, json
import numpy as np
import pandas as pd

CONF = [("_oltc", "Base (ideal OLTC)", "Base"), ("_fixed", "Fixed taps", "Fix"),
        ("_oltc_N1", "N-1, planning band", "Nmo"), ("_oltc_N1w", "N-1, 0.90--1.10 p.u.", "Nmow"),
        ("_oltc_S2", "Scenario 2 (future)", "Stwo")]
DIG = {str(i): w for i, w in enumerate(["Z", "O", "T", "Th", "Fo", "F", "Si", "Se", "E", "N"])}
rows, macros = [], {}
for suf, label, tag in CONF:
    os.environ["HC_SUFFIX"] = suf
    import importlib, analyze
    importlib.reload(analyze)
    avail, PB, _ = analyze.main()
    if len(avail) < 4:
        print("missing", suf, avail)
        continue
    allr = [r for g in avail for r in PB[g]]
    firm = np.array([r["firm"] for r in allr])
    nzero = int((firm <= 0.005).sum())
    if nzero:                                  # gains are undefined where firm HC is zero; report them separately
        print(tag, "buses with zero firm HC:", nzero)
        keepz = firm > 0.005
        allr = [r for r, k in zip(allr, keepz) if k]; firm = firm[keepz]
    g = lambda k: np.median((np.array([r[k] for r in allr]) / firm - 1) * 100)
    bt = np.array([r["bind_firm"] for r in allr])
    # The search labels the binding class at the lower (feasible) bracket. Next to a pre-existing overvoltage this
    # label can read "vhigh" although added load cannot worsen an overvoltage; for such buses use the class that the
    # boundary audit finds violated just above the headroom (dominant over the year).
    if (bt == "vhigh").any() and all(os.path.exists(f"../results/audit/boundary_{gname}{suf}.json") for gname in avail):
        k = 0
        for gname in avail:
            aud = json.load(open(f"../results/audit/boundary_{gname}{suf}.json"))["per_bus"]
            for j, r in enumerate(PB[gname]):
                if r["bind_firm"] == "vhigh":
                    ra = aud[j]["reason_above"]
                    bt[k + j] = max((c for c in ["vlow", "line", "trafo"]), key=lambda c: ra[c])
            k += len(PB[gname])
    stats_relabelled = int((np.array([r["bind_firm"] for r in allr]) != bt).sum())
    print(tag, "relabelled", stats_relabelled)
    g1 = (np.array([r["flex_0.01"] for r in allr]) / firm - 1) * 100
    d = dict(label=label, nbus=len(allr), firm=np.median(firm), g025=g("flex_0.0025"), g1=g("flex_0.01"), h100=g("hcH_100"),
             hrs1=np.median([r["hrs_0.01"] for r in allr]),
             nv=int((bt == "vlow").sum()), nl=int((bt == "line").sum()), nt=int((bt == "trafo").sum()),
             gv=np.median(g1[bt == "vlow"]) if (bt == "vlow").any() else float("nan"),
             gl=np.median(g1[bt == "line"]) if (bt == "line").any() else float("nan"))
    rows.append(d)
    for k, v in d.items():
        if k != "label":
            name = "".join(DIG.get(c, c) for c in f"sens{tag}{k}")
            macros[name] = f"{v:d}" if isinstance(v, int) else f"{v:.1f}"
    print(tag, {k: (round(v, 2) if isinstance(v, float) else v) for k, v in d.items()})

if rows:
    g1s = [d["g1"] for d in rows]; hs = [d["h100"] for d in rows]
    macros["sensRangeLo"] = f"{np.floor(min(g1s)):.0f}"; macros["sensRangeHi"] = f"{np.ceil(max(g1s)):.0f}"
    macros["sensHRangeLo"] = f"{np.floor(min(hs)):.0f}"; macros["sensHRangeHi"] = f"{np.ceil(max(hs)):.0f}"
with open("../tex/tab_sens.tex", "w") as f:
    f.write("\\begin{table}[t]\n\\centering\n\\caption{Sensitivity to Network Assumptions (Pooled Medians, All Buses)}\n\\label{tab:sens}\n\\setlength{\\tabcolsep}{2.4pt}\n")
    f.write("\\begin{tabular}{lcccccc}\n\\toprule\n")
    f.write("Configuration & Firm & \\multicolumn{2}{c}{Gain (\\%)} & Gain & Hours & Buses \\\\\n\\cmidrule(lr){3-4}\n")
    f.write(" & (MW) & 0.25\\% & 1\\% & 100 h (\\%) & at 1\\% & V/L/T \\\\\n\\midrule\n")
    for d in rows:
        f.write(f"{d['label']} & {d['firm']:.2f} & {d['g025']:.1f} & {d['g1']:.1f} & {d['h100']:.1f} & {d['hrs1']:.0f} & {d['nv']}/{d['nl']}/{d['nt']} \\\\\n")
        if d["label"].startswith("Fixed taps"):
            # non-flat new load (from shapes.py): daily swing of +-20 % around the mean, peak at 15:00
            sh = pd.read_csv("../results/shapes.csv"); x = sh[sh["shape"] == "DC40@15"]; fm = x.firm
            sv = dict(firm=fm.median(), g025=((x.e025 / fm - 1) * 100).median(), g1=((x.e1 / fm - 1) * 100).median(),
                      h100=((x.h100 / fm - 1) * 100).median(), hrs1=x.hrs1.median())
            f.write(f"Daily swing $\\pm$20\\% & {sv['firm']:.2f} & {sv['g025']:.1f} & {sv['g1']:.1f} & {sv['h100']:.1f} & {sv['hrs1']:.0f} & -- \\\\\n")
            macros.update(shpRowFirm=f"{sv['firm']:.2f}", shpRowGone=f"{sv['g1']:.1f}", shpRowHrs=f"{sv['hrs1']:.0f}")
    f.write("\\bottomrule\n\\end{tabular}\n\\vspace{2pt}\n\\begin{minipage}{\\columnwidth}\\footnotesize Daily swing: load $Pf(t)$ peaking at 15:00, mean $0.8P$; firm HC is its peak $P$; budgets relative to its own energy. Scenario 2 has 460 MV buses. V/L/T: buses whose firm HC is set by undervoltage, line, or transformer limits; one scenario-2 bus next to a pre-existing overvoltage is counted by the limit violated just above its headroom.\\end{minipage}\n\\end{table}\n")
with open("../tex/numbers_sens.tex", "w") as f:
    for k, v in sorted(macros.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
