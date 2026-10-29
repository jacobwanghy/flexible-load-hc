"""Figures, tables and LaTeX number macros from the stored headroom series."""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze import main, KEYS, ALPHA_GRID, LABEL, COL, FIG, RES, hours_hc, SUFFIX
from core import flexible_hc, curtailment_stats, DT_H

os.makedirs(FIG, exist_ok=True)
MK = {"rural": "o", "semiurb": "s", "urban": "^", "comm": "D"}
LS = {"rural": "-", "semiurb": "--", "urban": "-.", "comm": ":"}
avail, PB, SUM = main()
PK = json.load(open(f"{RES}/peaks.json"))
VAL = json.load(open(f"{RES}/validation2.json"))
PMs = {g: np.load(f"{RES}/{g}{SUFFIX}.npz")["PM"].astype(float) for g in avail}
macros = {}


DIG = {str(i): w for i, w in enumerate(["Z", "O", "T", "Th", "Fo", "F", "Si", "Se", "E", "N"])}


def mac(name, val, fmt="{:.1f}"):
    name = "".join(DIG.get(c, c) for c in name.replace("_", "").replace(".", ""))
    macros[name] = fmt.format(val) if not isinstance(val, str) else val


def gtag(g):
    return {"rural": "Rur", "semiurb": "Semi", "urban": "Urb", "comm": "Com"}[g]


# ---------------- Fig. 1: HC gain versus curtailed hours per year (pooled over all buses)
# For any nameplate P, the curtailed intervals are those with Pbar(t) < P, so the hour count depends on P alone:
# the hour-budget curve is the pooled sorted headroom, and each energy budget is one point on it.
plt.rcParams["font.family"] = "serif"
allPM = np.vstack([PMs[g] for g in avail])
srt = np.sort(allPM, axis=1)
firm_v = srt[:, 0]
Ns = np.unique(np.clip(np.round(np.logspace(np.log10(0.25), np.log10(8784), 140) / DT_H).astype(int), 1, allPM.shape[1] - 1))
G = srt[:, Ns] / firm_v[:, None] - 1
fig, ax = plt.subplots(figsize=(3.45, 2.15))
ax.fill_between(Ns * DT_H, np.percentile(G, 25, axis=0) * 100, np.percentile(G, 75, axis=0) * 100,
                color="#9a9a9a", alpha=0.25, lw=0, label="IQR, all buses")
ax.plot(Ns * DT_H, np.median(G, axis=0) * 100, "-", color="#1b6ca8", label="Median, hour budget $N$")
allrows = [r for g in avail for r in PB[g]]
for a, mk in [(0.0025, "o"), (0.005, "s"), (0.01, "D"), (0.05, "^")]:
    hh = np.median([r[f"hrs_{a}"] for r in allrows])
    gg = np.median([r[f"flex_{a}"] / r["firm"] - 1 for r in allrows]) * 100
    ax.plot(hh, gg, mk, color="#d1495b", ms=4.5, mfc="white", mew=1.1)
    ax.annotate(f"$\\alpha$={a*100:g}%", (hh, gg), xytext=(-4, 5), textcoords="offset points", ha="right", fontsize=6.5, color="#d1495b")
ax.plot([], [], "D", color="#d1495b", ms=4, mfc="white", label="Median, energy budget $\\alpha$")
ax.set_xscale("log")
ax.set_xlim(0.25, 8784)
ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1, 10, 100, 1000, 8784]))
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(0, 32)
ax.set_xlabel("Curtailed hours per year (h)")
ax.set_ylabel("HC gain over firm HC (%)")
ax.grid(True, which="major", color="#dddddd", lw=0.4)
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
for s_ in ["top", "right"]:
    ax.spines[s_].set_visible(False)
fig.tight_layout(pad=0.2)
fig.savefig(f"{FIG}/fig_duration.pdf")
plt.close(fig)

# ---------------- Fig. 2: flexible HC gain vs annual curtailment budget
fig, ax = plt.subplots(figsize=(3.45, 2.1))
x = ALPHA_GRID * 100
Call = np.array([r["curve"] for g in avail for r in PB[g]]) - 1
ax.fill_between(x, np.percentile(Call, 25, axis=0) * 100, np.percentile(Call, 75, axis=0) * 100,
                color="#9a9a9a", alpha=0.25, lw=0, label="IQR, all buses")
for g in avail:
    C = np.array([r["curve"] for r in PB[g]]) - 1
    med = np.median(C, axis=0) * 100
    ax.plot(x, med, LS[g], color=COL[g], marker=MK[g], ms=3, label=f"{LABEL[g]} median")
ax.set_xlabel(r"Annual curtailment budget $\alpha$ (% of energy)")
ax.set_ylabel(r"$H_k(\alpha)/H_k^{\mathrm{firm}}-1$ (%)")
ax.set_xlim(0, 5)
ax.set_ylim(bottom=0)
ax.grid(True, color="#dddddd", lw=0.4)
ax.set_ylim(0, 38)
ax.legend(frameon=False, loc="upper left", fontsize=6.5, ncol=2)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=0.2)
fig.savefig(f"{FIG}/fig_gain.pdf")
plt.close(fig)

# ---------------- Fig. 3: gain at 1% vs firm HC, by binding constraint
fig, ax = plt.subplots(figsize=(3.45, 2.1))
BM = {"vlow": ("o", "Undervoltage"), "line": ("s", "Line thermal"), "trafo": ("^", "Transformer"), "vhigh": ("x", "Overvoltage")}
BC = {"vlow": "#1b6ca8", "line": "#d1495b", "trafo": "#2a9d8f", "vhigh": "#7a7a7a"}
for bt in ["vlow", "line", "trafo"]:
    xs, ys = [], []
    for g in avail:
        for r in PB[g]:
            if r["bind_firm"] == bt:
                xs.append(r["firm"]); ys.append((r["flex_0.01"] / r["firm"] - 1) * 100)
    if xs:
        ax.scatter(xs, ys, s=7, marker=BM[bt][0], facecolors="none" if bt != "trafo" else BC[bt],
                   edgecolors=BC[bt], linewidths=0.6, label=f"{BM[bt][1]} ({len(xs)})")
ax.set_xscale("log")
from matplotlib.ticker import LogLocator, ScalarFormatter, NullFormatter
ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1, 2, 5, 10, 20, 50]))
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r"Firm hosting capacity $H_k^{\mathrm{firm}}$ (MW)")
ax.set_ylabel(r"Gain at $\alpha=1\%$ (%)")
ax.grid(True, color="#dddddd", lw=0.4)
ax.set_ylim(-2, 58)
ax.set_yticks([0, 10, 20, 30, 40])
leg = ax.legend(frameon=True, framealpha=0.92, edgecolor="none", loc="upper center", ncol=3, title="Binding constraint at worst interval", title_fontsize=7, fontsize=6.5)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout(pad=0.2)
fig.savefig(f"{FIG}/fig_binding.pdf")
plt.close(fig)

# ---------------- numbers for text and tables
allrows = [r for g in avail for r in PB[g]]
mac("nBus", len(allrows), "{:d}")
mac("nGrids", len(avail), "{:d}")
def sci(x):
    e = int(np.floor(np.log10(x)))
    return f"$10^{{{e+1}}}$"
macros["maxdV"] = sci(max(v["max_dV_pu"] for v in VAL.values()))
macros["maxdL"] = sci(max(max(v["max_dline_pct"], v["max_dtrafo_pct"]) for v in VAL.values()))
tab = []
for g in avail:
    rows = PB[g]; PM = PMs[g]; t = gtag(g)
    firm = np.array([r["firm"] for r in rows])
    snap = PM[:, PK[g]["t_peak_net"]] / firm - 1
    snapL = PM[:, PK[g]["t_peak_load"]] / firm - 1
    gains = {a: np.array([r[f"flex_{a}"] / r["firm"] - 1 for r in rows]) * 100 for a in [0.0025, 0.005, 0.01, 0.05]}
    hrs = {a: np.array([r[f"hrs_{a}"] for r in rows]) for a in [0.0025, 0.01]}
    ev = np.array([r["maxev_0.01"] for r in rows]); dep = np.array([r["depth_0.01"] for r in rows]) * 100
    hg = {H: np.array([r[f"hcH_{H}"] / r["firm"] - 1 for r in rows]) * 100 for H in [50, 100, 200]}
    bind = {k: sum(r["bind_firm"] == k for r in rows) for k in KEYS}
    tab.append(dict(g=g, n=len(rows), firm_med=np.median(firm), firm_min=firm.min(), firm_max=firm.max(),
                    snap_med=np.median(snap) * 100, snap_max=snap.max() * 100, snapL_max=snapL.max() * 100,
                    g025=np.median(gains[0.0025]), g05=np.median(gains[0.005]), g1=np.median(gains[0.01]),
                    g1q1=np.percentile(gains[0.01], 25), g1q3=np.percentile(gains[0.01], 75), g1max=gains[0.01].max(),
                    g5=np.median(gains[0.05]),
                    h025=np.median(hrs[0.0025]), h1=np.median(hrs[0.01]), ev1=np.median(ev), dep1=np.median(dep),
                    hg100=np.median(hg[100]), hg50=np.median(hg[50]), hg200=np.median(hg[200]), **{f"b_{k}": v for k, v in bind.items()}))
    for k, v in tab[-1].items():
        if k != "g":
            mac(f"{k}{t}", v,
                "{:d}" if isinstance(v, (int, np.integer)) else "{:.1f}")
df = pd.DataFrame(tab)
df.to_csv(f"{RES}/summary.csv", index=False)

# pooled statistics
def pooled(key):
    return np.array([r[key] for r in allrows])
firm_all = pooled("firm")
for a, nm in [(0.0025, "Ga"), (0.005, "Gb"), (0.01, "Gc"), (0.05, "Gd")]:
    gg = (pooled(f"flex_{a}") / firm_all - 1) * 100
    mac(f"pool{nm}med", np.median(gg)); mac(f"pool{nm}qa", np.percentile(gg, 25)); mac(f"pool{nm}qc", np.percentile(gg, 75)); mac(f"pool{nm}max", gg.max())
for a, nm in [(0.0025, "Ga"), (0.01, "Gc")]:
    mac(f"pool{nm}hrs", np.median(pooled(f"hrs_{a}")), "{:.0f}")
    mac(f"pool{nm}ev", np.median(pooled(f"maxev_{a}")), "{:.1f}")
    mac(f"pool{nm}dep", np.median(pooled(f"depth_{a}")) * 100, "{:.0f}")
    mac(f"pool{nm}nev", np.median(pooled(f"nev_{a}")), "{:.0f}")
    mac(f"pool{nm}depmax", np.max(pooled(f"depth_{a}")) * 100, "{:.0f}")
for a, nm in [(0.01, "Gc")]:
    mac(f"pool{nm}medev", np.median(pooled(f"medev_{a}")), "{:.1f}")
    mac(f"pool{nm}medevmin", np.median(pooled(f"medev_{a}")) * 60, "{:.0f}")
    mac(f"pool{nm}single", np.median(pooled(f"single_{a}")) * 100, "{:.0f}")
zz = pooled("zth")
mac("corrZth", np.corrcoef(zz, (pooled("flex_0.01") / firm_all - 1))[0, 1], "{:.2f}")
mac("corrFirm", np.corrcoef(np.log(firm_all), (pooled("flex_0.01") / firm_all - 1))[0, 1], "{:.2f}")
for H in [1, 10, 50, 100, 200]:
    gg = (pooled(f"hcH_{H}") / firm_all - 1) * 100
    mac(f"poolH{H}med", np.median(gg))
    mac(f"poolH{H}en", np.median(pooled(f"enH_{H}")) * 100, "{:.3f}")
dd = pooled("depth_0.01") * 100
mac("depthShareOver", (dd > 25).mean() * 100)
mac("depthPninety", np.percentile(dd, 90), "{:.0f}")
mac("firmPoolMed", np.median(firm_all), "{:.1f}")
bt = np.array([r["bind_firm"] for r in allrows])
g1 = (pooled("flex_0.01") / firm_all - 1) * 100
for k, nm in [("vlow", "Vlow"), ("line", "Line"), ("trafo", "Trafo")]:
    sel = bt == k
    mac(f"n{nm}", int(sel.sum()), "{:d}")
    if sel.any():
        mac(f"g{nm}", np.median(g1[sel]))
spread = (pooled("pmed") / firm_all - 1) * 100
mac("spreadMed", np.median(spread))
mac("corrSpread", np.corrcoef(g1, spread)[0, 1], "{:.2f}")
# timing of curtailment at alpha = 1%, pooled
cnt = np.zeros(PMs[avail[0]].shape[1])
for g in avail:
    for i, r in enumerate(PB[g]):
        cnt += PMs[g][i] < r["flex_0.01"]
# SimBench timestamps are local clock time (with daylight saving), used for time-of-day statistics
_tcol = pd.read_csv("/usr/local/lib/python3.13/site-packages/simbench/networks/1-complete_data-mixed-all-0-sw/LoadProfile.csv", sep=";", usecols=["time"]).time
idx = pd.DatetimeIndex(pd.to_datetime(_tcol, format="%d.%m.%Y %H:%M"))
s = pd.Series(cnt, idx)
ev_share = s[(s.index.hour >= 16) & (s.index.hour < 22)].sum() / s.sum() * 100
mac("eveShare", ev_share, "{:.0f}")
win = s[s.index.month.isin([11, 12, 1, 2])].sum() / s.sum() * 100
mac("winShare", win, "{:.0f}")
snap_pool = np.concatenate([PMs[g][:, PK[g]["t_peak_net"]] / np.array([r["firm"] for r in PB[g]]) - 1 for g in avail]) * 100
mac("snapOverShare", (snap_pool > 1).mean() * 100, "{:.0f}")
mac("snapOverMax", snap_pool.max())

with open("../tex/numbers.tex", "w") as f:
    for k, v in sorted(macros.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

# grid table data
gt = []
for g in avail:
    m = json.load(open(f"{RES}/{g}{SUFFIX}_meta.json"))
    gt.append(dict(g=LABEL[g], vn=m["vn_kv"], buses=m["n_bus"], lines=m["n_line"], trafo=m["trafo_mva"],
                   peak=PK[g]["peak_load"], peaknet=PK[g]["peak_net"], gen=m["sgen_mw"], viol=m["base_viol"]))
json.dump(gt, open(f"{RES}/gridtable.json", "w"), indent=1, default=float)
print(df.round(1).to_string())
print(json.dumps({k: macros[k] for k in sorted(macros) if k.startswith(("pool", "g", "n", "corr", "eve", "win", "snap", "max"))}, indent=0))

# ---------------- LaTeX tables
with open("../tex/tab_grids.tex", "w") as f:
    f.write("\\begin{table}[t]\n\\centering\n\\caption{SimBench MV Grids (Scenario 0)}\n\\label{tab:grids}\n\\setlength{\\tabcolsep}{2.8pt}\n")
    f.write("\\begin{tabular}{lcccccc}\n\\toprule\n")
    f.write("Grid & $V_n$ & Buses & HV/MV & Peak & Peak net & DG \\\\\n & (kV) & (MV) & (MVA) & load (MW) & load (MW) & (MW) \\\\\n\\midrule\n")
    for g in avail:
        m = json.load(open(f"{RES}/{g}{SUFFIX}_meta.json"))
        f.write(f"{LABEL[g]} & {m['vn_kv']:.0f} & {len(PB[g])} & $2\\times{m['trafo_mva']/2:.0f}$ & {PK[g]['peak_load']:.1f} & {PK[g]['peak_net']:.1f} & {m['sgen_mw']:.1f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\vspace{2pt}\n\\begin{minipage}{\\columnwidth}\\footnotesize MV buses joined by closed bus-bus switches are counted once.\\end{minipage}\n\\end{table}\n")

with open("../tex/tab_results.tex", "w") as f:
    f.write("\\begin{table}[t]\n\\centering\n\\caption{Firm HC and Median Flexibility Gains by Grid}\n\\label{tab:res}\n\\setlength{\\tabcolsep}{2.5pt}\n")
    f.write("\\begin{tabular}{lccccccc}\n\\toprule\n")
    f.write(" & Firm HC & \\multicolumn{3}{c}{Gain, energy budget (\\%)} & Gain & \\multicolumn{2}{c}{At $\\alpha=1\\%$} \\\\\n\\cmidrule(lr){3-5}\\cmidrule(lr){7-8}\n")
    f.write("Grid & median & 0.25\\% & 0.5\\% & 1\\% & 100 h & Hours & Longest \\\\\n & (MW) & & & & (\\%) & (h/yr) & event (h) \\\\\n\\midrule\n")
    for r in tab:
        f.write(f"{LABEL[r['g']]} & {r['firm_med']:.2f} & {r['g025']:.1f} & {r['g05']:.1f} & {r['g1']:.1f} & {r['hg100']:.1f} & {r['h1']:.0f} & {r['ev1']:.1f} \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n\\vspace{2pt}\n\\begin{minipage}{\\columnwidth}\\footnotesize  Gains are medians over buses of $H_k/H_k^{\\mathrm{firm}}-1$. The 100~h column allows curtailment in at most 100 hours per year.\\end{minipage}\n\\end{table}\n")
print("tables written")
