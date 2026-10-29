"""Post-process per-bus headroom series into hosting-capacity metrics, tables and figures."""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from core import flexible_hc, curtailment_stats, DT_H

GRIDS = [g for g in ["rural", "semiurb", "urban", "comm"]]
LABEL = {"rural": "Rural", "semiurb": "Semi-urban", "urban": "Urban", "comm": "Commercial"}
KEYS = ["vlow", "vhigh", "line", "trafo"]
ALPHAS = [0.0025, 0.005, 0.01, 0.02, 0.05]
ALPHA_GRID = np.array([0, 0.0005, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05])
HOURS = [1, 10, 50, 100, 200]
RES = "../results"
import os as _os
SUFFIX = _os.environ.get("HC_SUFFIX", "_oltc")
FIG = "../tex/figs"

plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "Times New Roman", "Times", "DejaVu Serif"], "mathtext.fontset": "stix",
                     "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.0,
                     "savefig.dpi": 300, "pdf.fonttype": 42})
COL = {"rural": "#1b6ca8", "semiurb": "#d1495b", "urban": "#2a9d8f", "comm": "#b8770b"}


def load(g):
    d = np.load(f"{RES}/{g}{SUFFIX}.npz")
    m = json.load(open(f"{RES}/{g}{SUFFIX}_meta.json"))
    return d, m


def hours_hc(pm, H):
    k = int(round(H / DT_H))
    s = np.sort(pm)
    return s[min(k, len(s) - 1)]


def main():
    avail = [g for g in GRIDS if __import__("os").path.exists(f"{RES}/{g}{SUFFIX}_meta.json")]
    summary, per_bus = {}, {}
    for g in avail:
        d, m = load(g)
        PM, BIND = d["PM"].astype(float), d["BIND"]
        T = PM.shape[1]
        # coincident-peak snapshot: step with the highest total existing net load
        # (recomputed here from the stored headroom is not possible; use argmin of feeder-wide mean headroom)
        rows = []
        for i, meta in enumerate(m["meta"]):
            pm = PM[i]
            firm = pm.min()
            r = dict(bus=meta["bus"], zth=meta["zth_pu"], firm=firm,
                     bind_firm=KEYS[BIND[i, pm.argmin()]],
                     cv=pm.std() / pm.mean(), pmed=np.median(pm), pmax_=pm.max(),
                     )
            for a in ALPHAS:
                P = flexible_hc(pm, a)
                st = curtailment_stats(pm, P)
                r[f"flex_{a}"] = P
                r[f"hrs_{a}"] = st["hours"]
                r[f"nev_{a}"] = st["n_events"]
                r[f"maxev_{a}"] = st["max_event_h"]
                r[f"depth_{a}"] = st["max_depth"]
                r[f"medev_{a}"] = st["med_event_h"]
                r[f"single_{a}"] = st["single_share"]
            for H in HOURS:
                P = hours_hc(pm, H)
                r[f"hcH_{H}"] = P
                r[f"enH_{H}"] = curtailment_stats(pm, P)["energy_frac"]
            r["curve"] = [flexible_hc(pm, a) / firm for a in ALPHA_GRID]
            # share of time each constraint binds
            r["bind_share"] = {k: float((BIND[i] == j).mean()) for j, k in enumerate(KEYS)}
            rows.append(r)
        per_bus[g] = rows
        summary[g] = dict(meta={k: v for k, v in m.items() if k != "meta"}, n=len(rows))
    return avail, per_bus, summary


if __name__ == "__main__":
    avail, per_bus, summary = main()
    for g in avail:
        rows = per_bus[g]
        firm = np.array([r["firm"] for r in rows])
        print(g, summary[g]["n"], "firm median", np.median(firm).round(2), "range", firm.min().round(2), firm.max().round(2))
        for a in ALPHAS:
            gain = np.array([r[f"flex_{a}"] / r["firm"] - 1 for r in rows]) * 100
            hrs = np.array([r[f"hrs_{a}"] for r in rows])
            print(f"  alpha {a*100:.2f}%  gain median {np.median(gain):.1f}% IQR [{np.percentile(gain,25):.1f},{np.percentile(gain,75):.1f}] max {gain.max():.1f}  hours median {np.median(hrs):.0f}")
        for H in HOURS:
            gain = np.array([r[f"hcH_{H}"] / r["firm"] - 1 for r in rows]) * 100
            en = np.array([r[f"enH_{H}"] for r in rows]) * 100
            print(f"  H={H}h gain median {np.median(gain):.1f}%  energy median {np.median(en):.3f}%")
        b = {}
        for r in rows:
            b[r["bind_firm"]] = b.get(r["bind_firm"], 0) + 1
        print("  binding at firm:", b, " base viol", summary[g]["meta"]["base_viol"])
