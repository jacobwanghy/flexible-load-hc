"""Independent pandapower checks on the stratified validation cases (brief 2.2, 3.2).

(a) For all cases in validate_headroom_all.csv: the limiting constraint class just above the computed headroom
    (Pbar + 6 kW) according to pandapower vs according to our solver.
(b) For up to 15 cases per grid: an independent pandapower scan from 0 to 1.3 * Pbar in 60 steps, checking that
    feasibility switches exactly once and that the switch brackets the computed headroom.
Element-level feasibility rule in both solvers: no element worse than min(0, its own base margin) - 1e-7.
Usage: python3 pp_scan.py <grid>
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, pandapower as pp, simbench as sb
from core import GridModel, GRIDS, TAN_PHI

grid = sys.argv[1]
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
cases = pd.read_csv(f"{RES}/audit/validate_headroom_all.csv")
cases = cases[cases.grid == grid].reset_index(drop=True)
meta = {m["bus"]: m for m in json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]}
gm = GridModel(GRIDS[grid], mode="oltc")
lim = gm.vmin > 0
keys = ["vlow", "vhigh", "line", "trafo"]

net = sb.get_simbench_net(GRIDS[grid])
prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
setp = float(net.trafo.autoTapSetp.iloc[0]); lvb = net.trafo.lv_bus.values
smax = float((net.trafo.sn_mva * net.trafo.max_loading_percent / 100).iloc[0])
net.trafo.in_service = False; net.ext_grid.in_service = False
eg_idx = [pp.create_ext_grid(net, b, vm_pu=setp) for b in lvb]
new = pp.create_load(net, net.bus.index[0], p_mw=0.0, q_mvar=0.0)
mv = (net.bus.vn_kv < 50).values
LP, LQ, GP = prof[("load", "p_mw")].values, prof[("load", "q_mvar")].values, prof[("sgen", "p_mw")].values
nload = len(LP[0]); grp = gm.lookup[lvb]


def pp_margins(bus, t, P):
    net.load.loc[net.load.index[:nload], "p_mw"] = LP[t]; net.load.loc[net.load.index[:nload], "q_mvar"] = LQ[t]
    net.sgen.p_mw = GP[t]; net.load.at[new, "bus"] = bus; net.load.at[new, "p_mw"] = P; net.load.at[new, "q_mvar"] = P * TAN_PHI
    try:
        pp.runpp(net, tolerance_mva=1e-9)
    except Exception:
        return None
    vm = net.res_bus.vm_pu.values[mv]
    eg = net.res_ext_grid.loc[eg_idx]
    I = []
    for L in np.unique(grp):
        sel = grp == L
        I += [abs((eg.p_mw.values[sel] + 1j * eg.q_mvar.values[sel]).sum()) / setp / sel.sum()] * int(sel.sum())
    return [(vm - net.bus.min_vm_pu.values[mv]) / 0.05, (net.bus.max_vm_pu.values[mv] - vm) / 0.05,
            1 - net.res_line.loading_percent.values / net.line.max_loading_percent.values, 1 - np.array(I) / smax]


def feasible_and_class(M, M0):
    worst, cls = np.inf, None
    ok = True
    for k, (m, m0) in enumerate(zip(M, M0)):
        slack = m - np.minimum(0, m0) + 1e-7
        if slack.min() < 0:
            ok = False
        if slack.min() < worst:
            worst, cls = slack.min(), keys[k]
    return ok, cls


def our_class(bus, t, P):
    e = np.zeros((gm.nb, 1)); e[meta[bus]["ppc"]] = 1.0
    S0 = gm.S0[:, [t]]; S = S0 - (1 + 1j * TAN_PHI) / gm.base * e * P
    out = []
    for SS in (S0, S):
        V, c = gm.solve(SS)
        Vm = np.abs(V[:, 0]); a, b = gm.line_rng
        Il = np.maximum(np.abs(gm.Yf[a:b] @ V[:, 0]), np.abs(gm.Yt[a:b] @ V[:, 0])) / gm.imax_line_pu
        out.append([(Vm[lim] - gm.vmin[lim]) / 0.05, (gm.vmax[lim] - Vm[lim]) / 0.05, 1 - Il, 1 - gm.trafo_loading(V, SS)[:, 0]])
    return feasible_and_class(out[1], out[0])[1]


res_a, res_b = [], []
t0 = time.time()
for j, r in cases.iterrows():
    bus, t, pb = int(r.bus), int(r.interval), float(r.pbar_ours)
    M0 = pp_margins(bus, t, 0.0)
    Mup = pp_margins(bus, t, pb + 0.006)
    ok_up, cls_pp = feasible_and_class(Mup, M0)
    res_a.append(dict(grid=grid, bus=bus, interval=t, class_pandapower=cls_pp, class_ours=our_class(bus, t, pb + 0.006), pp_feasible_above=ok_up))
    if j < 15:
        P = np.linspace(0, 1.3 * pb, 61)
        f = np.array([feasible_and_class(pp_margins(bus, t, p), M0)[0] for p in P])
        sw = np.flatnonzero(np.diff(f.astype(int)))
        contiguous = bool(f[0] and len(sw) == 1 and not f[-1])
        brackets = contiguous and P[sw[0]] <= pb + 1e-9 and pb <= P[sw[0] + 1] + 1e-9
        res_b.append(dict(grid=grid, bus=bus, interval=t, contiguous=contiguous, brackets_pbar=bool(brackets)))
    print(grid, j + 1, len(cases), f"{time.time()-t0:.0f}s", flush=True)
json.dump(dict(class_agreement=res_a, scans=res_b), open(f"{RES}/audit/ppscan_{grid}.json", "w"), indent=1)
a = pd.DataFrame(res_a); b = pd.DataFrame(res_b)
print("SUMMARY", grid, "class agree", int((a.class_pandapower == a.class_ours).sum()), "/", len(a),
      "pp feasible above", int(a.pp_feasible_above.sum()), "scans contiguous", int(b.contiguous.sum()), "/", len(b),
      "bracket", int(b.brackets_pbar.sum()))
