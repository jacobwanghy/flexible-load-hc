"""Small derived results used in the text: copperplate reconciliation with Gu et al., rank correlations,
deepest curtailment at the 100-h budget. Writes numbers_cp.tex and numbers_rank.tex."""
import json
import numpy as np, pandas as pd, simbench as sb
from scipy.stats import spearmanr
from core import GRIDS, flexible_hc

# (1) Gu et al. copperplate setup on SimBench profiles: transformer 10% above peak gross load, no DG,
#     flat new load, hour budget of 1% of intervals.
cp = {}
for g, c in GRIDS.items():
    p = sb.get_absolute_values(sb.get_simbench_net(c), profiles_instead_of_study_cases=True)
    L = p[("load", "p_mw")].values.sum(1)
    H = 1.1 * L.max() - L; s = np.sort(H); k = int(0.01 * len(H))
    cp[g] = float((s[k] / s[0] - 1) * 100)
json.dump(cp, open("../results/audit/copperplate_reconcile.json", "w"), indent=1)
open("../tex/numbers_cp.tex", "w").write("\\newcommand{\\cpLo}{%d}\n\\newcommand{\\cpHi}{%d}\n" % (int(round(min(cp.values()))), int(round(max(cp.values())))))

# (2) rank correlations and (3) depth at the 100-h budget
rows = []
for g in GRIDS:
    PM = np.load(f"../results/{g}_oltc.npz")["PM"].astype(float)
    meta = json.load(open(f"../results/{g}_oltc_meta.json"))["meta"]
    for i, pb in enumerate(PM):
        firm = pb.min(); P = flexible_hc(pb, 0.01); P100 = np.sort(pb)[400]
        cur = pb[pb < P100]
        rows.append((g, P / firm - 1, meta[i]["zth_pu"], firm, (P100 - cur).max() / P100 if len(cur) else 0.0))
d = pd.DataFrame(rows, columns=["g", "gain", "zth", "firm", "depth100"])
per = {g: (spearmanr(x.gain, x.zth)[0], spearmanr(x.gain, x.firm)[0]) for g, x in d.groupby("g")}
print("copperplate", cp); print("per-grid spearman (zth, firm)", per)
neg = lambda x: ("$-$%.2f" % abs(x)) if x < 0 else "%.2f" % x
open("../tex/numbers_rank.tex", "w").write("\\newcommand{\\spZth}{%s}\n\\newcommand{\\spFirm}{%s}\n\\newcommand{\\depthHundred}{%.0f}\n" % (
    neg(spearmanr(d.gain, d.zth)[0]), neg(spearmanr(d.gain, d.firm)[0]), d.depth100.median() * 100))
