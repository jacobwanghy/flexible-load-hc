"""Bounded corrective N-1 case (no outage statistics).

Preventive: the load is curtailed whenever an outage WOULD cause a violation, i.e., it follows the N-1 headroom.
Corrective: in normal operation only the intact-network (N-0) headroom limits the load; when an outage occurs, the load
must drop to the post-outage headroom until repair. Without outage rates and durations, curtailed energy cannot be
estimated, so this reports the bound: the nameplate allowed by the N-0 headroom at a 1 % budget, and the depth the load
must be able to shed after the worst outage, 1 - min_t PbarN1(t)/P (as % of P), per bus.
Output: tex/numbers_corr.tex
"""
import numpy as np
from core import flexible_hc

RES = "../results"
pre, cor, depth, extra = [], [], [], []
for g in ["rural", "semiurb", "urban", "comm"]:
    P0 = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    P1 = np.load(f"{RES}/{g}_oltc_N1.npz")["PM"].astype(float)
    for a, b in zip(P0, P1):
        Hp = flexible_hc(b, 0.01); Hc = flexible_hc(a, 0.01)
        pre.append(Hp); cor.append(Hc); extra.append(Hc - Hp)
        depth.append((1 - b.min() / Hc) * 100)
mac = dict(corPreMW=f"{np.median(pre):.2f}", corCorMW=f"{np.median(cor):.2f}", corExtraMW=f"{np.median(extra):.2f}",
           corDepth=f"{np.median(depth):.0f}", corDepthQa=f"{np.percentile(depth, 25):.0f}",
           corDepthQc=f"{np.percentile(depth, 75):.0f}")
print(mac)
with open("../tex/numbers_corr.tex", "w") as f:
    for k, v in sorted(mac.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
