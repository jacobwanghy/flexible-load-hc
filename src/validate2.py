"""Validate both solver modes against pandapower Newton-Raphson, including transformer loading."""
import json, time
import numpy as np, pandapower as pp, simbench as sb
from core import GridModel, GRIDS, TAN_PHI

rng = np.random.default_rng(1)
out = {}
for mode in ["fixed", "oltc"]:
    if mode == "fixed" and False: continue
    for name, code in GRIDS.items():
        gm = GridModel(code, mode=mode)
        prof = sb.get_absolute_values(sb.get_simbench_net(code), profiles_instead_of_study_cases=True)
        net0 = sb.get_simbench_net(code)
        mv = net0.bus.index[net0.bus.vn_kv < 50].values
        ts = rng.choice(gm.T, 25, replace=False)
        ks = rng.choice(mv, 25)
        P = rng.uniform(0, 4, 25)
        S = gm.S0[:, ts].copy()
        for j, k in enumerate(ks):
            S[gm.lookup[k], j] -= (1 + 1j * TAN_PHI) * P[j] / gm.base
        V, c = gm.solve(S)
        TL = gm.trafo_loading(V, S)
        ev, el, et = [], [], []
        for j, t in enumerate(ts):
            net = sb.get_simbench_net(code)
            net.load.p_mw = prof[("load", "p_mw")].iloc[t].values
            net.load.q_mvar = prof[("load", "q_mvar")].iloc[t].values
            net.sgen.p_mw = prof[("sgen", "p_mw")].iloc[t].values
            pp.create_load(net, ks[j], p_mw=P[j], q_mvar=P[j] * TAN_PHI)
            if mode == "oltc":
                # equivalent reference: MV busbars as slack at the tap set point, transformers removed
                setp = float(net.trafo.autoTapSetp.iloc[0])
                lvb = net.trafo.lv_bus.values
                net.trafo.in_service = False
                net.ext_grid.in_service = False
                for b in np.unique(lvb):
                    pp.create_ext_grid(net, b, vm_pu=setp)
            pp.runpp(net, tolerance_mva=1e-10)
            mvb = net.bus.vn_kv < 50
            vm_pp = net.res_bus.vm_pu[mvb].values
            vm_me = np.abs(V[gm.lookup[net.bus.index[mvb]], j])
            ev.append(np.nanmax(np.abs(vm_pp - vm_me)))
            a, b = gm.line_rng
            Il = np.maximum(np.abs(gm.Yf[a:b] @ V[:, j]), np.abs(gm.Yt[a:b] @ V[:, j])) / gm.imax_line_pu * 100
            el.append(np.abs(net.res_line.loading_percent.values - Il).max())
            if mode == "fixed":
                tl_pp = net.res_trafo.loading_percent.values / net.trafo.max_loading_percent.values
            else:
                eg = net.res_ext_grid.iloc[1:]
                grp = gm.lookup[net.ext_grid.bus.values[1:]]
                smax = (net.trafo.sn_mva * net.trafo.max_loading_percent / 100).values[0]
                tl_pp = []
                for L in np.unique(grp):
                    sel = grp == L
                    Stot = (eg.p_mw.values[sel] + 1j * eg.q_mvar.values[sel]).sum()
                    vb = abs(V[L, j])
                    n = gm.n_trafo_at[int(L)]
                    tl_pp += [abs(Stot) / vb / n / smax] * n
                tl_pp = np.array(tl_pp)
            tl_me = np.sort(TL[:, j]); tl_pp = np.sort(tl_pp)
            et.append(np.abs(tl_me - tl_pp).max() * 100)
        out[f"{name}_{mode}"] = dict(max_dV_pu=float(max(ev)), max_dline_pct=float(max(el)), max_dtrafo_pct=float(max(et)), conv=bool(c.all()))
        print(name, mode, out[f"{name}_{mode}"], flush=True)
json.dump(out, open("../results/validation2.json", "w"), indent=1)
