"""Full-population boundary audit with the vectorized solver (brief 3.5, 3.6).

For every bus and every interval of the base case (ideal OLTC):
  * solve at P = Pbar(t) and at P = Pbar(t) + DELTA;
  * feasible at Pbar: every constraint class satisfied (pre-existing violations not worsened);
  * infeasible at Pbar + DELTA, and why (constraint class or non-convergence);
  * element-level check that no bus overvoltage is worsened (class-level relaxation cannot hide one);
  * count non-converged columns.
Usage: python3 audit_boundary.py <grid> [mode] [line_limit]
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from core import GridModel, GRIDS, TAN_PHI

DELTA = 0.006           # MW; search tolerance is 0.005 MW
STORE = 1e-5            # MW; headroom is stored as float32, so test feasibility 10 W below the stored value
TOL = 1e-7              # same rounding allowance as the search
grid = sys.argv[1]
mode = sys.argv[2] if len(sys.argv) > 2 else "oltc"
ll = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
scen = int(sys.argv[4]) if len(sys.argv) > 4 else 0
tag = f"{grid}_{mode}" + ("" if ll == 1.0 else f"_L{int(ll*100)}") + ("" if scen == 0 else f"_S{scen}")
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
PM = np.load(f"{RES}/{tag}.npz")["PM"].astype(float)
meta = json.load(open(f"{RES}/{tag}_meta.json"))["meta"]
gm = GridModel(GRIDS[grid].replace("--0-", f"--{scen}-"), mode=mode, line_limit=ll)
keys = ["vlow", "vhigh", "line", "trafo"]
lim = gm.vmin > 0

V0, c0 = gm.solve(gm.S0)
mg0 = gm.margins(V0, c0, gm.S0)
base = np.vstack([mg0[k] for k in keys])
relax = np.minimum(base, 0)
over0 = np.maximum(np.abs(V0[lim]) - gm.vmax[lim, None], 0)   # element-level overvoltage before connection

out = []
t0 = time.time()
for i, m in enumerate(meta):
    p = m["ppc"]
    e = np.zeros((gm.nb, 1)); e[p] = 1.0
    unit = -(1 + 1j * TAN_PHI) / gm.base * e
    res = {"bus": m["bus"]}
    for name, P in [("at", PM[i] - STORE), ("above", PM[i] + DELTA)]:
        S = gm.S0 + unit * P[None, :]
        V, c = gm.solve(S, V0=V0)
        mg = gm.margins(V, c, S)
        M = np.vstack([mg[k] for k in keys]) - relax + np.where(relax < 0, TOL, 0.0)
        feas = (M.min(axis=0) >= 0) & c
        res[f"nonconv_{name}"] = int((~c).sum())
        if name == "at":
            res["infeasible_at"] = int((~feas).sum())
            over = np.maximum(np.abs(V[lim]) - gm.vmax[lim, None], 0)
            res["overvoltage_worsened"] = int(((over - over0) > 1e-6).any(axis=0).sum())
            res["minmargin_at"] = float(M.min())
        else:
            res["feasible_above"] = int(feas.sum())
            viol = ~feas
            why = np.where(~c, "nonconv", np.array(keys)[M.argmin(axis=0)])
            res["reason_above"] = {k: int(((why == k) & viol).sum()) for k in keys + ["nonconv"]}
    out.append(res)
    if (i + 1) % 10 == 0:
        print(tag, i + 1, len(meta), f"{time.time()-t0:.0f}s", flush=True)

summ = dict(tag=tag, buses=len(out), intervals=len(out) * PM.shape[1],
            infeasible_at=sum(r["infeasible_at"] for r in out),
            feasible_above=sum(r["feasible_above"] for r in out),
            nonconv_at=sum(r["nonconv_at"] for r in out),
            nonconv_above=sum(r["nonconv_above"] for r in out),
            overvoltage_worsened=sum(r["overvoltage_worsened"] for r in out),
            reason_above={k: sum(r["reason_above"][k] for r in out) for k in keys + ["nonconv"]},
            base_intervals_with_preexisting_violation={k: int((base[j] < 0).sum()) for j, k in enumerate(keys)})
json.dump(dict(summary=summ, per_bus=out), open(f"{RES}/audit/boundary_{tag}.json", "w"), indent=1)
print("SUMMARY", json.dumps(summ))
