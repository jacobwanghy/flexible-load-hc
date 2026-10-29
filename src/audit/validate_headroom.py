"""Stratified validation of the computed headroom against an independent pandapower Newton-Raphson search
(brief section 5).

Sample (per grid): for each constraint class that sets firm HC (undervoltage, line, transformer), up to 8 buses;
for each bus its worst and second-worst intervals, two random intervals, and one interval with a
pre-existing overvoltage when the grid has any. For each case pandapower finds the headroom independently by
bisection to 1 kW, using its own network model with the ideal OLTC represented by external grids at the MV
busbars. Reports the headroom error (kW), false-feasible / false-infeasible counts at the computed boundary,
and voltage / line / transformer loading differences at the boundary operating point.
Usage: python3 validate_headroom.py <grid>
"""
import sys, os, json, time, copy
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandapower as pp, simbench as sb
from core import GridModel, GRIDS, TAN_PHI

grid = sys.argv[1]
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
d = np.load(f"{RES}/{grid}_oltc.npz"); PM = d["PM"].astype(float); BIND = d["BIND"]
meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
keys = ["vlow", "vhigh", "line", "trafo"]
gm = GridModel(GRIDS[grid], mode="oltc")
V0all, c0 = gm.solve(gm.S0)
mg0 = gm.margins(V0all, c0, gm.S0)
ovt = np.where(mg0["vhigh"] < 0)[0]

# ---------------- independent pandapower model (ideal OLTC: busbars held by external grids)
net = sb.get_simbench_net(GRIDS[grid])
prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
setp = float(net.trafo.autoTapSetp.iloc[0])
lvb = net.trafo.lv_bus.values
smax = float((net.trafo.sn_mva * net.trafo.max_loading_percent / 100).iloc[0])
net.trafo.in_service = False
net.ext_grid.in_service = False
eg_idx = [pp.create_ext_grid(net, b, vm_pu=setp) for b in lvb]
new = pp.create_load(net, net.bus.index[0], p_mw=0.0, q_mvar=0.0, name="new")
mv = (net.bus.vn_kv < 50).values
LP, LQ, GP = prof[("load", "p_mw")].values, prof[("load", "q_mvar")].values, prof[("sgen", "p_mw")].values
nload = len(LP[0])


def run(bus, t, P):
    net.load.loc[net.load.index[:nload], "p_mw"] = LP[t]
    net.load.loc[net.load.index[:nload], "q_mvar"] = LQ[t]
    net.sgen.p_mw = GP[t]
    net.load.at[new, "bus"] = bus
    net.load.at[new, "p_mw"] = P
    net.load.at[new, "q_mvar"] = P * TAN_PHI
    try:
        pp.runpp(net, tolerance_mva=1e-9)
    except Exception:
        return None
    vm = net.res_bus.vm_pu.values[mv]
    eg = net.res_ext_grid.loc[eg_idx]
    grp = gm.lookup[lvb]
    I = []
    for L in np.unique(grp):
        sel = grp == L
        I += [abs((eg.p_mw.values[sel] + 1j * eg.q_mvar.values[sel]).sum()) / setp / sel.sum()] * int(sel.sum())
    tl = np.array(I) / smax
    m = dict(vlow=((vm - net.bus.min_vm_pu.values[mv]) / 0.05).min(),
             vhigh=((net.bus.max_vm_pu.values[mv] - vm) / 0.05).min(),
             line=(1 - net.res_line.loading_percent.values / net.line.max_loading_percent.values).min(),
             trafo=(1 - tl).min())
    return m, vm, net.res_line.loading_percent.values.copy(), tl * 100


def feasible(m, relax):
    return m is not None and all(m[k] - relax[k] + (1e-7 if relax[k] < 0 else 0) >= 0 for k in keys)


rng = np.random.default_rng(5)
cases = []
for cls in ["vlow", "line", "trafo"]:
    cand = [i for i, mm in enumerate(meta) if keys[BIND[i, PM[i].argmin()]] == cls]
    for i in rng.choice(cand, min(8, len(cand)), replace=False) if cand else []:
        ts = [int(PM[i].argmin()), int(np.argsort(PM[i])[1]), int(rng.integers(PM.shape[1])), int(rng.integers(PM.shape[1]))]
        if len(ovt):
            ts.append(int(rng.choice(ovt)))
        for t in ts:
            cases.append((int(i), t, cls))

out = []
t0 = time.time()
for (i, t, cls) in cases:
    bus = meta[i]["bus"]; pb = PM[i, t]
    b0 = run(bus, t, 0.0)
    relax = {k: min(b0[0][k], 0) for k in keys}
    # independent bisection: bracket [0, hi] with hi infeasible
    lo, hi = 0.0, max(2 * pb, 1.0)
    while True:
        r = run(bus, t, hi)
        if not feasible(r[0] if r else None, relax):
            break
        lo, hi = hi, hi * 2
    while hi - lo > 0.001:
        mid = 0.5 * (lo + hi)
        r = run(bus, t, mid)
        if feasible(r[0] if r else None, relax):
            lo = mid
        else:
            hi = mid
    r_at = run(bus, t, pb - 1e-5)
    r_up = run(bus, t, pb + 0.006)
    # our solver at the same operating point, for V / loading agreement
    e = np.zeros((gm.nb, 1)); e[meta[i]["ppc"]] = 1.0
    S = gm.S0[:, [t]] - (1 + 1j * TAN_PHI) / gm.base * e * (pb - 1e-5)
    V, c = gm.solve(S)
    vm_me = np.abs(V[gm.lookup[net.bus.index[mv]], 0])
    a, b = gm.line_rng
    Il = np.maximum(np.abs(gm.Yf[a:b] @ V[:, 0]), np.abs(gm.Yt[a:b] @ V[:, 0])) / gm.imax_line_pu * 100
    tl_me = gm.trafo_loading(V, S)[:, 0] * 100
    out.append(dict(grid=grid, bus=bus, interval=t, firm_class=cls, preexisting_overvoltage=bool(relax["vhigh"] < 0),
                    pbar_ours=pb, pbar_pandapower=0.5 * (lo + hi), err_kW=(pb - 0.5 * (lo + hi)) * 1000,
                    feasible_at_ours=feasible(r_at[0] if r_at else None, relax),
                    infeasible_above_ours=not feasible(r_up[0] if r_up else None, relax),
                    dV_pu=float(np.abs(r_at[1] - vm_me).max()), dline_pct=float(np.abs(r_at[2] - Il).max()),
                    dtrafo_pct=float(np.abs(np.sort(r_at[3]) - np.sort(tl_me)).max())))
    print(grid, len(out), len(cases), f"{time.time()-t0:.0f}s", f"err {out[-1]['err_kW']:+.2f} kW", flush=True)
json.dump(out, open(f"{RES}/audit/validate_headroom_{grid}.json", "w"), indent=1)
print("DONE", grid)
