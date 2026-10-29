"""Independent check of stored headroom values with pandapower Newton-Raphson.

For random (bus, interval) pairs: at P = Pbar the network must satisfy all limits
(except limits already violated before connection, which must not worsen), and at
P = Pbar + 20 kW at least one limit must be violated or worsened.
"""
import json, sys
import numpy as np, pandapower as pp, simbench as sb
from core import GRIDS, TAN_PHI

rng = np.random.default_rng(7)
TOL = 1e-6
report = {}
for g, code in GRIDS.items():
    d = np.load(f"../results/{g}_oltc.npz"); PM = d["PM"].astype(float)
    meta = json.load(open(f"../results/{g}_oltc_meta.json"))["meta"]
    base = sb.get_simbench_net(code)
    prof = sb.get_absolute_values(base, profiles_instead_of_study_cases=True)
    ok_at, viol_above, n = 0, 0, 30
    worst = []
    for _ in range(n):
        i = rng.integers(len(meta)); t = rng.integers(PM.shape[1])
        bus = meta[i]["bus"]; P = PM[i, t]

        def run(Padd):
            net = sb.get_simbench_net(code)
            net.load.p_mw = prof[("load", "p_mw")].iloc[t].values
            net.load.q_mvar = prof[("load", "q_mvar")].iloc[t].values
            net.sgen.p_mw = prof[("sgen", "p_mw")].iloc[t].values
            if Padd > 0:
                pp.create_load(net, bus, p_mw=Padd, q_mvar=Padd * TAN_PHI)
            # ideal OLTC: MV busbars held at the tap set point; transformer flows from the busbar injections
            setp = float(net.trafo.autoTapSetp.iloc[0])
            lvb = net.trafo.lv_bus.values
            smax = (net.trafo.sn_mva * net.trafo.max_loading_percent / 100).values
            net.trafo.in_service = False; net.ext_grid.in_service = False
            for b in lvb:
                pp.create_ext_grid(net, b, vm_pu=setp)
            try:
                pp.runpp(net)
            except Exception:
                return None
            mv = (net.bus.vn_kv < 50).values
            vm = net.res_bus.vm_pu.values[mv]
            eg = net.res_ext_grid.iloc[1:]
            grp = net._pd2ppc_lookups["bus"][net.ext_grid.bus.values[1:]]   # fused busbars share one group
            I_per = []
            for L in np.unique(grp):
                sel = grp == L
                Sg = np.abs((eg.p_mw.values[sel] + 1j * eg.q_mvar.values[sel]).sum())
                I_per += [Sg / setp / sel.sum()] * sel.sum()   # equal sharing between identical transformers
            I_per = np.array(I_per)
            m = np.concatenate([vm - net.bus.min_vm_pu.values[mv], net.bus.max_vm_pu.values[mv] - vm,
                                (net.line.max_loading_percent.values - net.res_line.loading_percent.values) / 100 * 0.05,
                                (1 - I_per / smax[0]) * 0.05])
            return m

        m0 = run(0.0)
        relax = np.minimum(m0, 0)
        mP = run(P)
        mU = run(P + 0.02)
        a = mP is not None and np.all(mP - relax >= -TOL)
        b = mU is None or np.any(mU - relax < -TOL)
        ok_at += a; viol_above += b
        if not (a and b):
            worst.append(dict(bus=int(bus), t=int(t), P=float(P), feasible_at=bool(a), violated_above=bool(b)))
    report[g] = dict(samples=n, feasible_at_Pbar=int(ok_at), violated_at_Pbar_plus_20kW=int(viol_above), failures=worst)
    print(g, report[g], flush=True)
json.dump(report, open("../results/boundary_check.json", "w"), indent=1)
