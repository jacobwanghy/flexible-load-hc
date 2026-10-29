"""Uncertainty of the headline medians (pooled over buses).

Buses on one feeder share a binding constraint and are not independent. Confidence intervals therefore use a
cluster bootstrap: within each grid (stratum), feeder groups (ring groups: connected MV network once the MV busbars are
removed and all normally open ties closed; a candidate on a busbar is its own cluster) are resampled with replacement
and all their buses kept; the pooled median is recomputed (B = 5000, percentile 95 % interval, seed 7).
Leave-one-grid-out: the pooled median of the 1 % gain without each grid in turn.
Output: results/stats.json, tex/numbers_stats.tex
"""
import json, copy
import numpy as np, simbench as sb
from core import GRIDS, flexible_hc, curtailment_stats
from n1 import ring_groups

RES = "../results"
B, SEED = 5000, 7


def per_bus():
    rows = []
    for g, code in GRIDS.items():
        PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
        meta = json.load(open(f"{RES}/{g}_oltc_meta.json"))["meta"]
        comps, _ = ring_groups(sb.get_simbench_net(code))
        for i, m in enumerate(meta):
            pm = PM[i]; firm = pm.min()
            gi = [j for j, c in enumerate(comps) if m["bus"] in c]
            P1 = flexible_hc(pm, 0.01); st = curtailment_stats(pm, P1)
            s = np.sort(pm)
            rows.append(dict(grid=g, cluster=f"{g}-{gi[0]}" if gi else f"{g}-bb{m['bus']}", firm=firm,
                             g025=(flexible_hc(pm, 0.0025) / firm - 1) * 100, g1=(P1 / firm - 1) * 100,
                             h100=(s[400] / firm - 1) * 100, add1=P1 - firm, hrs1=st["hours"],
                             dep1=st["max_depth"] * 100))
    return rows


if __name__ == "__main__":
    rows = per_bus()
    keys = ["firm", "g025", "g1", "h100", "add1", "hrs1", "dep1"]
    X = {k: np.array([r[k] for r in rows]) for k in keys}
    grid = np.array([r["grid"] for r in rows]); clus = np.array([r["cluster"] for r in rows])
    strata = {g: {c: np.flatnonzero(clus == c) for c in np.unique(clus[grid == g])} for g in GRIDS}
    rng = np.random.default_rng(SEED)
    boot = {k: [] for k in keys}
    for _ in range(B):
        idx = []
        for g, cl in strata.items():
            names = list(cl)
            for c in rng.choice(len(names), len(names), replace=True):
                idx.extend(cl[names[c]])
        idx = np.array(idx)
        for k in keys:
            boot[k].append(np.median(X[k][idx]))
    out = {k: dict(median=float(np.median(X[k])), lo=float(np.percentile(boot[k], 2.5)),
                   hi=float(np.percentile(boot[k], 97.5))) for k in keys}
    out["clusters"] = {g: len(c) for g, c in strata.items()}
    out["logo_g1"] = {g: float(np.median(X["g1"][grid != g])) for g in GRIDS}
    out["per_grid_g1"] = {g: float(np.median(X["g1"][grid == g])) for g in GRIDS}
    json.dump(out, open(f"{RES}/stats.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
    lo = min(out["logo_g1"].values()); hi = max(out["logo_g1"].values())
    mac = dict(ciGoneLo=f"{out['g1']['lo']:.1f}", ciGoneHi=f"{out['g1']['hi']:.1f}",
               ciGqLo=f"{out['g025']['lo']:.1f}", ciGqHi=f"{out['g025']['hi']:.1f}",
               ciHLo=f"{out['h100']['lo']:.1f}", ciHHi=f"{out['h100']['hi']:.1f}",
               ciFirmLo=f"{out['firm']['lo']:.2f}", ciFirmHi=f"{out['firm']['hi']:.2f}",
               ciHrsLo=f"{out['hrs1']['lo']:.0f}", ciHrsHi=f"{out['hrs1']['hi']:.0f}",
               ciAddLo=f"{out['add1']['lo']:.2f}", ciAddHi=f"{out['add1']['hi']:.2f}",
               logoLo=f"{lo:.1f}", logoHi=f"{hi:.1f}", nClusters=str(sum(out["clusters"].values())),
               depOne=f"{out['dep1']['median']:.0f}", hrsOneShare=f"{out['hrs1']['median'] / 8784 * 100:.0f}")
    with open("../tex/numbers_stats.tex", "w") as f:
        for k, v in sorted(mac.items()):
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
