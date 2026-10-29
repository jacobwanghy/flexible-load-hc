"""Monotonicity of feasibility in added load (brief 3.5).

For every bus: a stratified set of intervals (its 20 lowest-headroom intervals, up to 20 intervals with
a pre-existing overvoltage, 60 random intervals), each swept over 61 load levels from 0 to 1.5 * Pbar(t).
Checks that the feasible set is exactly the interval [0, Pbar(t)] (one switch from feasible to infeasible,
at Pbar within the 5 kW search tolerance), and that each constraint-class margin moves monotonically
(undervoltage, line and transformer margins non-increasing; overvoltage margin non-decreasing).
Also writes margin-vs-load curves for one bus per grid and binding class (for the audit report figure).
Usage: python3 audit_monotone.py <grid>
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from core import GridModel, GRIDS, TAN_PHI

grid = sys.argv[1]
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
d = np.load(f"{RES}/{grid}_oltc.npz"); PM = d["PM"].astype(float); BIND = d["BIND"]
meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
gm = GridModel(GRIDS[grid], mode="oltc")
keys = ["vlow", "vhigh", "line", "trafo"]
V0, c0 = gm.solve(gm.S0)
mg0 = gm.margins(V0, c0, gm.S0)
base = np.vstack([mg0[k] for k in keys])
relax = np.minimum(base, 0)
ovt = np.where(base[1] < 0)[0]
rng = np.random.default_rng(11)
NP = 61
fr = np.linspace(0, 1.5, NP)
stats = dict(grid=grid, cases=0, noncontiguous=0, boundary_err_max_kW=0.0, nonmonotone={k: 0 for k in keys}, nonconv=0)
curves = {}
for i, m in enumerate(meta):
    pb = PM[i]
    ts = np.unique(np.concatenate([np.argsort(pb)[:20], rng.choice(ovt, min(20, len(ovt)), replace=False) if len(ovt) else [],
                                   rng.choice(len(pb), 60, replace=False)]).astype(int))
    tt = np.repeat(ts, NP)
    P = (pb[ts][:, None] * fr[None, :]).ravel()
    e = np.zeros((gm.nb, 1)); e[m["ppc"]] = 1.0
    S = gm.S0[:, tt] - (1 + 1j * TAN_PHI) / gm.base * e * P[None, :]
    V, c = gm.solve(S, V0=V0[:, tt])
    mg = gm.margins(V, c, S)
    M = np.vstack([mg[k] for k in keys]) - relax[:, tt] + np.where(relax[:, tt] < 0, 1e-7, 0.0)
    feas = ((M.min(axis=0) >= 0) & c).reshape(len(ts), NP)
    stats["nonconv"] += int((~c).sum())
    stats["cases"] += len(ts)
    for j in range(len(ts)):
        f = feas[j]
        sw = np.flatnonzero(np.diff(f.astype(int)))
        if not (f[0] and len(sw) == 1 and not f[-1]):
            stats["noncontiguous"] += 1
        else:
            last = fr[sw[0]] * pb[ts[j]]; first_bad = fr[sw[0] + 1] * pb[ts[j]]
            # Pbar must lie between the last feasible and first infeasible sweep point (5 kW tolerance)
            err = max(0.0, last - pb[ts[j]] - 0.005, pb[ts[j]] - first_bad - 0.005) * 1000
            stats["boundary_err_max_kW"] = max(stats["boundary_err_max_kW"], err)
    for q, k in enumerate(keys):
        raw = np.vstack([mg[k]]).reshape(len(ts), NP)
        dd = np.diff(raw, axis=1)
        ok = c.reshape(len(ts), NP)[:, 1:] & c.reshape(len(ts), NP)[:, :-1]
        bad = (dd < -1e-9) & ok if k == "vhigh" else (dd > 1e-9) & ok
        stats["nonmonotone"][k] += int(bad.any(axis=1).sum())
    bc = keys[BIND[i, pb.argmin()]]
    if bc not in curves:
        t0 = int(pb.argmin())
        j = int(np.flatnonzero(ts == t0)[0])
        curves[bc] = dict(bus=m["bus"], interval=t0, pbar=float(pb[t0]), P=(fr * pb[t0]).tolist(),
                          margins={k: np.asarray(mg[k]).reshape(len(ts), NP)[j].tolist() for k in keys})
json.dump(dict(stats=stats, curves=curves), open(f"{RES}/audit/monotone_{grid}.json", "w"))
print("SUMMARY", json.dumps(stats))
