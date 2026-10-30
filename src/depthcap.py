"""Hosting capacity with curtailment depth capped at a share d of nameplate (d = 25 %).

The load may be curtailed by at most d*P in any interval, so P - Pbar(t) <= d P for all t, i.e. P <= firm/(1-d).
With an energy budget alpha: H_cap = min(H(alpha), firm/(1-d)). Reported for N-0 and N-1 (planning band), alpha = 1 %.
Output: tex/numbers_depth.tex
"""
import numpy as np
from core import flexible_hc

RES = "../results"
D = 0.25
mac = {}
for tag, suf in [("NZ", "_oltc"), ("NO", "_oltc_N1")]:
    g, gc, cap = [], [], []
    for grid in ["rural", "semiurb", "urban", "comm"]:
        for pm in np.load(f"{RES}/{grid}{suf}.npz")["PM"].astype(float):
            f = pm.min(); H = flexible_hc(pm, 0.01); Hc = min(H, f / (1 - D))
            g.append((H / f - 1) * 100); gc.append((Hc / f - 1) * 100); cap.append(H > f / (1 - D) + 1e-9)
    mac[f"dep{tag}gain"] = f"{np.median(gc):.1f}"
    mac[f"dep{tag}share"] = f"{np.mean(cap) * 100:.0f}"
    print(suf, "uncapped median %.1f, capped median %.1f, share of buses capped %.0f%%" % (np.median(g), np.median(gc), np.mean(cap) * 100))
mac["depCap"] = f"{D * 100:.0f}"
mac["depMaxGain"] = f"{D / (1 - D) * 100:.0f}"
with open("../tex/numbers_depth.tex", "w") as fh:
    for k, v in sorted(mac.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
