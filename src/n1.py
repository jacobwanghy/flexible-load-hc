"""N-1 line-outage sensitivity (preventive criterion with restoration switching).

For every in-service MV line l of a grid (contingency), the line is taken out of service and the normally open tie
line that reconnects every bus it cut off is closed (if several ties can do so, the one with the largest worst-case
margin over the year without new load is used). If no tie can, the cut-off buses are treated as unsupplied: their
loads and DG are set to zero and their own limits are not checked. The intact network is included.
For a candidate bus k, the N-1 headroom is
    Pbar_k^{N-1}(t) = min over the intact network and all contingencies of Pbar_k^c(t),
where Pbar_k^c is Eq. (2) evaluated in the post-restoration topology c (pre-existing violations in topology c may not
worsen). Limits are the same as in the base case (voltage band, 100 % line rating, SimBench transformer limit).

Speed: with the ideal OLTC the MV busbars are slack buses, so a contingency changes voltages and currents only in the
ring group that contains the outaged line. The power flow is solved on that group only; the other feeders' currents
into the shared busbar are taken from a full solve of topology c and added for the transformer limit (exact, since they
do not depend on the new load). A contingency outside bus k's group changes k's results only through busbar losses,
which are ignored. Each bus's headroom is therefore the minimum over the intact case and the contingencies of its group.
Usage: python3 n1.py <grid>   ->  results/<grid>_oltc_N1.npz (+ _meta.json)
"""
import sys, os, json, time, copy
import numpy as np, scipy.sparse as spr
import simbench as sb, pandapower as pp
from core import GridModel, GRIDS, pmax_series, TAN_PHI
from n1_topo import contingencies

KEYS = ["vlow", "vhigh", "line", "trafo"]


class LocalModel:
    """Power flow and margins restricted to one ring group of topology gm (gm: GridModel without profiles)."""

    def __init__(self, gm, S0full, group_ppc, lost_ppc=(), V_full=None):
        self.base = gm.base
        lost = set(int(x) for x in lost_ppc)
        S0full = S0full.copy()
        if lost:
            S0full[list(lost)] = 0.0                      # unsupplied buses: no load, no DG
        slack = [int(s) for s in gm.slack]
        g = [int(x) for x in group_ppc if int(x) not in set(slack)]
        self.loc = np.array(g + slack)
        self.ng = len(g)
        self.pos = {int(p): i for i, p in enumerate(self.loc)}
        self.nb = len(self.loc)
        Y = gm.Yprime
        Ygg = Y[np.ix_(g, g)]
        self.Z = np.linalg.inv(Ygg)
        self.vS = gm.vS.copy()
        self.w = -(self.Z @ (Y[np.ix_(g, slack)] @ self.vS))[:, None]
        self.S0 = S0full[self.loc]
        # voltage limits of supplied group buses
        self.vmin = np.array([gm.vmin[p] if p not in lost else 0.0 for p in g])
        self.vmax = np.array([gm.vmax[p] for p in g])
        self.lim = self.vmin > 0
        # lines touching the group (both ends in group or slack), excluding lines into lost buses
        a, b = gm.line_rng
        Yf = gm.Yf[a:b].tocsr(); Yt = gm.Yt[a:b].tocsr()
        locset = set(self.loc.tolist()); gset = set(g)
        fb = gm.line_fbus
        tb = gm.line_tbus
        keep = []
        for i in range(b - a):
            ends = {int(fb[i]), int(tb[i])}
            if ends <= locset and ends & gset and not ends & lost and Yf[i].nnz:
                keep.append(i)
        cols = self.loc
        self.Yf = Yf[keep][:, cols].toarray()          # small dense blocks: BLAS is much faster than sparse @ dense
        self.Yt = Yt[keep][:, cols].toarray()
        self.imax = gm.imax_line_pu[keep]
        self.line_idx = np.array(keep)
        # transformer current at busbars: group part + rest of the network (fixed, from a full solve of topology c)
        if V_full is None:
            V_full, _ = gm.solve(S0full)
        self.tr = []
        ta, tb_ = gm.trafo_rng
        for j in range(tb_ - ta):
            L = int(gm.trafo_tbus[j])
            rest = np.setdiff1d(np.arange(gm.nb), self.loc)
            I_rest = Y[L, rest] @ V_full[rest] if len(rest) else 0.0
            self.tr.append((self.pos[L], Y[L, self.loc], I_rest, gm.n_trafo_at[L], gm.smax_trafo_pu[j]))

    def subset(self, cols):
        """Same model restricted to a subset of intervals."""
        c = copy.copy(self)
        if hasattr(c, "_base"):
            del c._base
        c.S0 = self.S0[:, cols]
        c.tr = [(pL, yrow, (I_rest[cols] if np.ndim(I_rest) else I_rest), ntr, smax) for (pL, yrow, I_rest, ntr, smax) in self.tr]
        return c

    def feasible_at(self, k, P, tan_phi):
        """Eq. (2) feasibility of load P(t) at local bus k, per interval (same rule and tolerance as pmax_series)."""
        unit = np.zeros((self.nb, 1), dtype=complex); unit[k] = -(1 + 1j * tan_phi) / self.base
        if not hasattr(self, "_base"):                 # base state of this topology, shared by all candidate buses
            V0, c0 = self.solve(self.S0)
            b = self.margins(V0, c0, self.S0)
            self._base = (V0, np.vstack([np.minimum(b[q], 0) for q in KEYS]))
        V0, relax = self._base
        S = self.S0 + unit * P[None, :]
        V, c = self.solve(S, V0=V0)
        m = self.margins(V, c, S)
        M = np.vstack([m[q] for q in KEYS]) - relax + np.where(relax < 0, 1e-7, 0.0)
        return M.min(axis=0) >= 0

    def solve(self, S, vs=None, V0=None, tol=1e-9, maxit=60):
        Sn = S[:self.ng]
        Vn = np.ones_like(Sn) * self.vS[0] if V0 is None else V0[:self.ng].copy()
        for it in range(maxit):
            Vnew = self.w + self.Z @ np.conj(Sn / Vn)
            d = np.abs(Vnew - Vn).max(axis=0)
            Vn = Vnew
            if np.all(d < tol):
                break
        conv = (d < 1e-6) & np.isfinite(Vn).all(axis=0) & (np.abs(Vn).min(axis=0) >= 0.5)
        V = np.empty((self.nb, S.shape[1]), dtype=complex)
        V[:self.ng] = Vn
        V[self.ng:] = self.vS[:, None]
        return V, conv

    def margins(self, V, conv, S):
        Vm = np.abs(V[:self.ng])
        if self.lim.any():
            vl = ((Vm[self.lim] - self.vmin[self.lim, None]) / 0.05).min(axis=0)
            vh = ((self.vmax[self.lim, None] - Vm[self.lim]) / 0.05).min(axis=0)
        else:
            vl = vh = np.ones(V.shape[1])
        if len(self.imax):
            I = np.maximum(np.abs(self.Yf @ V), np.abs(self.Yt @ V)) / self.imax[:, None]
            ln = (1 - I).min(axis=0)
        else:
            ln = np.ones(V.shape[1])
        tm = np.ones(V.shape[1])
        for (pL, yrow, I_rest, ntr, smax) in self.tr:
            inj = yrow @ V + I_rest - np.conj(S[pL] / V[pL])
            tm = np.minimum(tm, 1 - np.abs(inj) / ntr / smax)
        out = dict(vlow=vl, vhigh=vh, line=ln, trafo=tm)
        return {k: np.where(conv, v, -1.0) for k, v in out.items()}


def injections_by_pp_bus(net):
    prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
    T = prof[("load", "p_mw")].shape[0]
    S = np.zeros((net.bus.index.max() + 1, T), dtype=complex)
    np.add.at(S, net.load.bus.values, -(prof[("load", "p_mw")].values.T + 1j * prof[("load", "q_mvar")].values.T))
    np.add.at(S, net.sgen.bus.values, prof[("sgen", "p_mw")].values.T)
    if ("storage", "p_mw") in prof and prof[("storage", "p_mw")].shape[1] and len(net.storage):
        np.add.at(S, net.storage.bus.values, -prof[("storage", "p_mw")].values.T)
    return S


def to_ppc(gm, Spp):
    S = np.zeros((gm.nb, Spp.shape[1]), dtype=complex)
    idx = gm.net.bus.index.values
    np.add.at(S, gm.lookup[idx], Spp[idx])
    return S / gm.base


def ring_groups(net):
    """pandapower MV buses per ring group (all ties closed, busbars removed)."""
    import pandapower.topology as top, networkx as nx
    pp.runpp(net, numba=False)
    lk = net._pd2ppc_lookups["bus"]
    bb = set(lk[net.trafo.lv_bus.values])
    G = top.create_nxgraph(net, respect_switches=False, include_trafos=False)
    mv = [b for b in G.nodes if net.bus.vn_kv[b] < 50 and lk[b] not in bb]
    comps = list(nx.connected_components(G.subgraph(mv)))
    line_group = {}
    for l in net.line.index:
        for i, c in enumerate(comps):
            if net.line.from_bus[l] in c or net.line.to_bus[l] in c:
                line_group[int(l)] = i
    return comps, line_group


def component(gm, Yb, p):
    """Non-slack buses electrically connected to ppc bus p in topology gm (feeders meet only at slack busbars)."""
    stop = set(gm.slack.tolist()) | set(gm.hv.tolist())
    seen = {int(p)}; st = [int(p)]
    while st:
        u = st.pop()
        for v in Yb[u].indices:
            if v not in seen and v not in stop:
                seen.add(int(v)); st.append(int(v))
    return frozenset(seen)


if __name__ == "__main__":
    grid = sys.argv[1]
    # post-outage voltage band: "plan" = the base-case planning band; "wide" = 0.90-1.10 p.u. (EN 50160 supply range)
    band = sys.argv[2] if len(sys.argv) > 2 else "plan"
    RES = "../results"
    code = GRIDS[grid]
    net0 = sb.get_simbench_net(code)
    Spp = injections_by_pp_bus(net0)
    meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
    PM0 = np.load(f"{RES}/{grid}_oltc.npz")["PM"].astype(float)
    comps, line_group = ring_groups(copy.deepcopy(net0))
    cs, ties = contingencies(copy.deepcopy(net0))
    tie_sw = net0.switch[(~net0.switch.closed) & (net0.switch.et == "l")]
    bus_group = {int(b): i for i, c in enumerate(comps) for b in c}
    gm0 = GridModel(code, profiles=False, mode="oltc", net=copy.deepcopy(net0))
    S00 = to_ppc(gm0, Spp)
    V00, _ = gm0.solve(S00)
    TL0 = gm0.trafo_loading(V00, S00)
    trafo_lv_pp = net0.trafo.lv_bus.values
    PM1 = PM0.copy()
    WORST = np.full(PM0.shape, -1, dtype=np.int16)        # contingency (line id) setting the N-1 headroom, -1 = intact
    BIND1 = np.load(f"{RES}/{grid}_oltc.npz")["BIND"].copy()
    log = []
    t0 = time.time()
    tag = f"{grid}_oltc_N1" + ("" if band == "plan" else "w")
    ckpt = f"{RES}/{tag}_checkpoint.npz"
    done = 0
    if os.path.exists(ckpt):                              # resume an interrupted run
        z = np.load(ckpt, allow_pickle=True)
        PM1, WORST, BIND1 = z["PM"].astype(float), z["WORST"], z["BIND"]
        log = json.loads(str(z["log"])); done = len(log)
        print("resuming after", done, "contingencies", flush=True)
    for ci, c in enumerate(cs):
        if ci < done:
            continue
        gi = line_group[c["line"]]
        net = copy.deepcopy(net0)
        # outage: open the line's switch at one end (the line stays in pandapower's model as a dangling, unloaded line,
        # exactly as the normally open tie lines are modelled in the base case, so branch indices are unchanged)
        lsw = net.switch.index[(net.switch.et == "l") & (net.switch.element == c["line"])]
        net.switch.loc[lsw[0], "closed"] = False
        options = c["ties"] if c["lost"] and c["ties"] else [None]
        best = None
        for tl in options:
            n2 = copy.deepcopy(net)
            if tl is not None:
                n2.switch.loc[tie_sw.index[tie_sw.element == tl], "closed"] = True
            lost = c["lost"] if tl is None else []
            if lost:                                  # unrestorable: keep the line closed (no isolated buses), zero injection beyond it
                n2.switch.loc[lsw[0], "closed"] = True
            gm = GridModel(code, profiles=False, mode="oltc", net=n2)
            if band == "wide":
                gm.vmin = np.where(gm.vmin > 0, 0.90, 0.0); gm.vmax = np.where(gm.vmax < 9, 1.10, 9.0)
            S0 = to_ppc(gm, Spp)
            lost_ppc = sorted(set(int(gm.lookup[b]) for b in lost))
            if lost_ppc:
                S0[lost_ppc] = 0.0
            Vf, cf = gm.solve(S0)
            mg = gm.margins(Vf, cf, S0)
            if lost_ppc:                              # unsupplied buses are not checked
                keep_v = np.ones(gm.nb, bool); keep_v[lost_ppc] = False
            score = min(v.min() for v in mg.values())
            if best is None or score > best[0]:
                best = (score, tl, gm, S0, Vf, lost, lost_ppc, {k: float(v.min()) for k, v in mg.items()})
        score, tl, gm, S0, Vf, lost, lost_ppc, base_min = best
        # candidate buses: the outaged line's ring group, plus buses fed from a busbar whose transformer loading changed
        TL = gm.trafo_loading(Vf, S0)
        chg_lv = set(int(trafo_lv_pp[j]) for j in range(len(TL)) if np.abs(TL[j] - TL0[j]).max() > 0.005)
        chg_ppc = set(int(gm.lookup[b]) for b in chg_lv)
        Yb = spr.csr_matrix(gm.Yprime != 0)
        cache = {}
        n_eval = 0; n_int = 0
        for i, m in enumerate(meta):
            b = m["bus"]
            if b in lost:
                continue
            p = int(gm.lookup[b])
            if p in set(gm.slack.tolist()):
                if p not in chg_ppc:
                    continue
                # candidate on a busbar (transformer-limited): full-network search, all intervals
                pm, bd, _ = pmax_series(gm, p, S0=S0)
                upd = pm < PM1[i]
                PM1[i, upd] = pm[upd]; WORST[i, upd] = c["line"]
                BIND1[i, upd] = np.array([KEYS.index(x) for x in bd[upd]], dtype=np.int8)
                n_eval += 1; n_int += gm.T if hasattr(gm, "T") else PM1.shape[1]
                continue
            else:
                comp = component(gm, Yb, p)
                feeds = set(int(v) for u in comp for v in Yb[u].indices if v in set(gm.slack.tolist()))
                if bus_group.get(b) != gi and not (feeds & chg_ppc):
                    continue
            if comp not in cache:
                cache[comp] = LocalModel(gm, S0, sorted(comp), lost_ppc, V_full=Vf)
            lm = cache[comp]
            # only intervals where the current N-1 headroom is infeasible in topology c can lower it
            ok = lm.feasible_at(lm.pos[p], PM1[i], TAN_PHI)
            cols = np.flatnonzero(~ok)
            if len(cols):
                pm, bd, _ = pmax_series(lm.subset(cols), lm.pos[p])
                upd = pm < PM1[i, cols]
                PM1[i, cols[upd]] = pm[upd]; WORST[i, cols[upd]] = c["line"]
                BIND1[i, cols[upd]] = np.array([KEYS.index(x) for x in bd[upd]], dtype=np.int8)
                n_int += len(cols)
            n_eval += 1
        log.append(dict(line=c["line"], group=gi, tie=tl, lost=lost, n_options=len(options), base_min_margin=base_min,
                        buses=n_eval, intervals_searched=n_int, busbars_changed=sorted(chg_lv)))
        print(grid, ci + 1, len(cs), "group", gi, "tie", tl, "buses", n_eval, "intervals", n_int,
              f"basemin={score:.3f}", f"{time.time()-t0:.0f}s", flush=True)
        if (ci + 1) % 5 == 0:
            np.savez(ckpt, PM=PM1, WORST=WORST, BIND=BIND1, log=json.dumps(log))
    np.savez_compressed(f"{RES}/{tag}.npz", PM=PM1.astype(np.float32), WORST=WORST, BIND=BIND1)
    json.dump(dict(grid=grid, log=log, seconds=time.time() - t0), open(f"{RES}/{tag}_run.json", "w"), indent=1)
    if os.path.exists(ckpt):
        os.remove(ckpt)
    print("DONE", tag, time.time() - t0)
