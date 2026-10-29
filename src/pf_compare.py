"""Unity power-factor sensitivity on a random 25-bus subset per grid (same buses as the base case)."""
import json, numpy as np
from core import flexible_hc
rows = []
for g in ["rural", "semiurb", "urban", "comm"]:
    b1 = json.load(open(f"../results/{g}_oltc_pf100_meta.json"))["meta"]; P1 = np.load(f"../results/{g}_oltc_pf100.npz")["PM"].astype(float)
    b0 = {m["bus"]: i for i, m in enumerate(json.load(open(f"../results/{g}_oltc_meta.json"))["meta"])}
    P0 = np.load(f"../results/{g}_oltc.npz")["PM"].astype(float)
    for i, m in enumerate(b1):
        j = b0[m["bus"]]
        f0, f1 = P0[j].min(), P1[i].min()
        rows.append((f1 / f0 - 1, flexible_hc(P0[j], 0.01) / f0 - 1, flexible_hc(P1[i], 0.01) / f1 - 1))
r = np.array(rows) * 100
print("n", len(r), "firm change median %.1f%%, gain 0.95: %.1f%%, gain 1.0: %.1f%%" % tuple(np.median(r, axis=0)))
open("../tex/numbers_pf.tex", "w").write("\\newcommand{\\pfN}{%d}\n\\newcommand{\\pfFirm}{%.1f}\n\\newcommand{\\pfGainA}{%.1f}\n\\newcommand{\\pfGainB}{%.1f}\n" % ((len(r),) + tuple(np.median(r, axis=0))))
