"""Yardstick: on-site battery needed to firm the flexible nameplate (brief reviewer concern).

For a load of nameplate P = H_k(alpha) at bus k, a battery must supply P - Pbar(t) whenever positive and may
recharge from the grid only within the remaining headroom Pbar(t) - P. Required energy = max of the deficit
recursion d(t) = max(0, d(t-1) + (P - Pbar(t)) dt) (perfect foresight, lossless, starting full).
Output: median battery power and energy per MW of added capacity, and in hours of added capacity.
"""
import numpy as np, pandas as pd
from numba import njit
from core import flexible_hc

@njit(cache=True)
def deficit(pb, P, dt):
    d = 0.0; dmax = 0.0
    for v in pb:
        d = max(0.0, d + (P - v) * dt)
        if d > dmax:
            dmax = d
    return dmax

rows = []
for g in ["rural", "semiurb", "urban", "comm"]:
    PM = np.load(f"../results/{g}_oltc.npz")["PM"].astype(float)
    for i, pb in enumerate(PM):
        firm = pb.min()
        for a in [0.0025, 0.01]:
            P = flexible_hc(pb, a)
            add = P - firm
            if add <= 1e-6:
                continue
            Pb = float(np.max(P - pb)); Eb = deficit(pb, P, 0.25)
            rows.append(dict(grid=g, bus=i, alpha=a, add_MW=add, batt_MW=Pb, batt_MWh=Eb, batt_h_of_add=Eb / add, batt_MW_per_add=Pb / add))
d = pd.DataFrame(rows); d.to_csv("../results/battery.csv", index=False)
s = d.groupby("alpha")[["add_MW", "batt_MW", "batt_MWh", "batt_h_of_add", "batt_MW_per_add"]].median()
print(s.round(2))
with open("../tex/numbers_batt.tex", "w") as f:
    for a, tag in [(0.01, "O"), (0.0025, "Q")]:
        f.write(f"\\newcommand{{\\batt{tag}MWh}}{{{s.loc[a,'batt_MWh']:.1f}}}\n\\newcommand{{\\batt{tag}h}}{{{s.loc[a,'batt_h_of_add']:.0f}}}\n\\newcommand{{\\batt{tag}MW}}{{{s.loc[a,'batt_MW']:.2f}}}\n\\newcommand{{\\batt{tag}add}}{{{s.loc[a,'add_MW']:.2f}}}\n")
