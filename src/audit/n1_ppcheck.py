"""Independent pandapower check of the N-1 headroom.

For random bus-intervals whose N-1 headroom is set by an outage, rebuild that post-restoration topology in pandapower
(outaged line's switch opened, chosen tie closed, transformers replaced by external grids at the OLTC set point as in
the base-case audits), and check with Newton-Raphson that the stored headroom is feasible and 6 kW above it is not,
using the element-level rule of Eq. (2) relative to the same topology without the new load.
Usage: python3 n1_ppcheck.py <grid> [n_cases]  -> results/audit/n1_ppcheck_<grid>.json
"""
import sys, os, json, copy
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandapower as pp, simbench as sb
from core import GRIDS, TAN_PHI

grid = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 15
RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
z = np.load(f"{RES}/{grid}_oltc_N1.npz"); PM, W = z["PM"].astype(float), z["WORST"]
meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
log = {x["line"]: x for x in json.load(open(f"{RES}/{grid}_oltc_N1_run.json"))["log"]}
net0 = sb.get_simbench_net(GRIDS[grid])
prof = sb.get_absolute_values(net0, profiles_instead_of_study_cases=True)
LP, LQ, GP = prof[("load", "p_mw")].values, prof[("load", "q_mvar")].values, prof[("sgen", "p_mw")].values
SP = prof[("storage", "p_mw")].values if ("storage", "p_mw") in prof and len(net0.storage) else None
setp = float(net0.trafo.autoTapSetp.iloc[0])
smax = (net0.trafo.sn_mva * net0.trafo.max_loading_percent / 100).values
tie_sw = net0.switch[(~net0.switch.closed) & (net0.switch.et == "l")]
rng = np.random.default_rng(5)
cand = np.argwhere(W >= 0)
pick = cand[rng.choice(len(cand), min(N, len(cand)), replace=False)]


def build(line):
    net = copy.deepcopy(net0)
    c = log[line]
    lsw = net.switch.index[(net.switch.et == "l") & (net.switch.element == line)]
    if not c["lost"]:
        net.switch.loc[lsw[0], "closed"] = False
        if c["tie"] is not None:
            net.switch.loc[tie_sw.index[tie_sw.element == c["tie"]], "closed"] = True
    lvb = net.trafo.lv_bus.values
    net.trafo.in_service = False; net.ext_grid.in_service = False
    eg = [pp.create_ext_grid(net, b, vm_pu=setp) for b in lvb]
    new = pp.create_load(net, net.bus.index[0], p_mw=0.0)
    return net, eg, new, set(c["lost"])


def margins(net, eg, new, lost, t, bus, P):
    nl = len(LP[0])
    net.load.loc[net.load.index[:nl], "p_mw"] = LP[t]; net.load.loc[net.load.index[:nl], "q_mvar"] = LQ[t]
    net.sgen.p_mw = GP[t]
    if SP is not None:
        net.storage.p_mw = SP[t]
    for b in lost:                                   # unsupplied buses: no load, no DG
        net.load.loc[net.load.bus == b, ["p_mw", "q_mvar"]] = 0.0; net.sgen.loc[net.sgen.bus == b, "p_mw"] = 0.0
    net.load.at[new, "bus"] = bus; net.load.at[new, "p_mw"] = P; net.load.at[new, "q_mvar"] = P * TAN_PHI
    pp.runpp(net, tolerance_mva=1e-9)
    mv = (net.bus.vn_kv < 50).values & ~net.bus.index.isin(list(lost))
    vm = net.res_bus.vm_pu.values[mv]
    m = [(vm - net.bus.min_vm_pu.values[mv]) / 0.05, (net.bus.max_vm_pu.values[mv] - vm) / 0.05]
    ok = net.line.in_service.values & ~net.line.from_bus.isin(list(lost)).values & ~net.line.to_bus.isin(list(lost)).values
    ld = net.res_line.loading_percent.values[ok] / net.line.max_loading_percent.values[ok]
    m.append(1 - np.nan_to_num(ld))
    r = net.res_ext_grid.loc[eg]
    lvb = net0.trafo.lv_bus.values
    # transformer current: busbar injection shared by the units feeding a busbar group (as in the solver)
    S = np.abs(r.p_mw.values + 1j * r.q_mvar.values)
    m.append(1 - S / setp / smax)
    return m


out = []
for i, t in pick:
    pbar = PM[i, t]; line = int(W[i, t]); bus = meta[i]["bus"]
    net, eg, new, lost = build(line)
    M0 = margins(net, eg, new, lost, t, bus, 0.0)
    feas = lambda M: all((m - np.minimum(0, m0) + 1e-7).min() >= 0 for m, m0 in zip(M, M0))
    a = feas(margins(net, eg, new, lost, t, bus, max(pbar - 1e-5, 0)))
    b = feas(margins(net, eg, new, lost, t, bus, pbar + 0.006))
    out.append(dict(bus=bus, interval=int(t), line=line, pbar=pbar, feasible_at=bool(a), infeasible_above=not b))
    print(grid, bus, t, line, f"{pbar:.3f}", a, not b, flush=True)
json.dump(out, open(f"{RES}/audit/n1_ppcheck_{grid}.json", "w"), indent=1)
print("SUMMARY", grid, "cases", len(out), "feasible at", sum(x["feasible_at"] for x in out),
      "infeasible 6 kW above", sum(x["infeasible_above"] for x in out))
