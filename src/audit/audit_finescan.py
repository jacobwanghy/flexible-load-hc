"""Fine-resolution monotonicity scan (brief 2.2), own solver, element-level feasibility rule.

Stratified bus-intervals per grid: every bus's worst interval; for each bus, the interval with the largest reverse
power flow (lowest aggregate net load), an interval with pre-existing overvoltage if any, and two random intervals.
Each case is scanned from 0 to 1.3 * Pbar(t) in 5 kW steps (search tolerance). Feasibility uses the ELEMENT-LEVEL rule:
every element's margin m_e(P) >= min(0, m_e(0)) - 1e-7 (no element may be worsened beyond its own base violation,
no new violation anywhere). Reports: cases, any feasible point above the first infeasible point (non-monotone),
first-infeasible step vs stored Pbar (should be within 5 kW), non-convergence.
Usage: python3 audit_finescan.py <grid>
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from core import GridModel, GRIDS, TAN_PHI

grid = sys.argv[1]
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
PM = np.load(f"{RES}/{grid}_oltc.npz")["PM"].astype(float)
meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
gm = GridModel(GRIDS[grid], mode="oltc")
lim = gm.vmin > 0
STEP, TOLE = 0.005, 1e-7


def element_margins(V, S):
    Vm = np.abs(V)
    a, b = gm.line_rng
    Il = np.maximum(np.abs(gm.Yf[a:b] @ V), np.abs(gm.Yt[a:b] @ V)) / gm.imax_line_pu[:, None]
    return np.vstack([(Vm[lim] - gm.vmin[lim, None]) / 0.05, (gm.vmax[lim, None] - Vm[lim]) / 0.05,
                      1 - Il, 1 - gm.trafo_loading(V, S)])


V0, c0 = gm.solve(gm.S0)
E0 = element_margins(V0, gm.S0)
ovt = np.where((E0[lim.sum():2 * lim.sum()] < 0).any(axis=0))[0]
netload = -gm.S0.real.sum(axis=0)
t_rev = int(np.argmin(netload))
rng = np.random.default_rng(21)
cases = []
for i in range(len(meta)):
    ts = {int(PM[i].argmin()), t_rev, *map(int, rng.choice(PM.shape[1], 2, replace=False))}
    if len(ovt):
        ts.add(int(rng.choice(ovt)))
    cases += [(i, t) for t in ts]
stats = dict(grid=grid, cases=len(cases), nonmonotone=0, first_infeasible_err_kW_max=0.0, nonconv=0,
             steps_total=0)
bad = []
t0 = time.time()
for (i, t) in cases:
    pb = PM[i, t]
    P = np.arange(0, 1.3 * pb + STEP, STEP)
    e = np.zeros((gm.nb, 1)); e[meta[i]["ppc"]] = 1.0
    S = gm.S0[:, [t]] - (1 + 1j * TAN_PHI) / gm.base * e * P[None, :]
    V, c = gm.solve(S, V0=np.repeat(V0[:, [t]], len(P), axis=1))
    E = element_margins(V, S)
    feas = np.all(E >= np.minimum(0, E0[:, [t]]) - TOLE, axis=0) & c
    stats["nonconv"] += int((~c).sum()); stats["steps_total"] += len(P)
    first_bad = np.flatnonzero(~feas)
    if len(first_bad) == 0:
        bad.append(dict(bus=meta[i]["bus"], t=t, issue="no infeasible point within 1.3*Pbar")); continue
    fb = first_bad[0]
    if feas[fb:].any():
        stats["nonmonotone"] += 1
        bad.append(dict(bus=meta[i]["bus"], t=t, issue="feasible again after first infeasible point"))
    P_first_bad = P[fb]
    err = max(0.0, pb - P_first_bad, P_first_bad - pb - 2 * STEP) * 1000      # Pbar must lie in [P_first_bad - step, P_first_bad]
    stats["first_infeasible_err_kW_max"] = max(stats["first_infeasible_err_kW_max"], err)
    if err > 0:
        bad.append(dict(bus=meta[i]["bus"], t=t, issue=f"boundary mismatch {err:.2f} kW", pbar=pb, first_bad=P_first_bad))
stats["seconds"] = round(time.time() - t0)
json.dump(dict(stats=stats, issues=bad[:200]), open(f"{RES}/audit/finescan_{grid}.json", "w"), indent=1)
print("SUMMARY", json.dumps(stats))
