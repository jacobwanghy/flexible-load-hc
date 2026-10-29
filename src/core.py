"""Vectorized time-series hosting-capacity engine for SimBench MV grids.

The network model (Ybus, branch admittances, limits) is taken from pandapower's
internal case so that the solver uses exactly the same physics as pandapower.
Power flows for all 15-min steps of a year are solved at once with the implicit
Z-bus fixed-point method:  V_n = Z_nn (conj(S_n / V_n) - Y_ns V_s).
"""
import numpy as np
import pandapower as pp
import simbench as sb

GRIDS = {
    "rural": "1-MV-rural--0-sw",
    "semiurb": "1-MV-semiurb--0-sw",
    "urban": "1-MV-urban--0-sw",
    "comm": "1-MV-comm--0-sw",
}

PF_NEW = 0.95                      # power factor of the new large load (lagging)
TAN_PHI = np.tan(np.arccos(PF_NEW))
DT_H = 0.25                        # 15-min resolution


class GridModel:
    def __init__(self, code, profiles=True, mode="oltc", line_limit=1.0, net=None):
        """mode: 'oltc' = ideal on-load tap changer holding each HV/MV transformer's MV
        busbar at its SimBench set point; 'fixed' = taps held at their given positions."""
        self.mode = mode
        net = sb.get_simbench_net(code) if net is None else net   # a modified net may be passed (N-1 study)
        self.code = code
        self.net = net
        pp.runpp(net, numba=False)
        ppci = net._ppc["internal"]
        self.base = ppci["baseMVA"]
        Y = ppci["Ybus"].toarray()
        bus_types = ppci["bus"][:, 1].astype(int)
        self.Yf = ppci["Yf"]
        self.Yt = ppci["Yt"]
        self.nb = Y.shape[0]
        self.lookup = net._pd2ppc_lookups["bus"]
        br = net._pd2ppc_lookups["branch"]
        self.line_rng = br["line"]
        self.trafo_rng = br["trafo"]
        branch = ppci["branch"]
        fbus = branch[:, 0].real.astype(int)
        tbus = branch[:, 1].real.astype(int)
        ta, tb = self.trafo_rng
        if mode == "fixed":
            self.slack = np.where(bus_types == 3)[0]
            assert len(self.slack) == 1
            self.vS = np.full(1, net.ext_grid.vm_pu.values[0] * np.exp(1j * np.deg2rad(net.ext_grid.va_degree.values[0])))
            self.hv = np.array([], dtype=int)
        else:
            # remove transformer branches; MV busbars become slack buses at the tap set point
            nbr = branch.shape[0]
            import scipy.sparse as spr
            rows = np.arange(ta, tb)
            Cf = spr.csr_matrix((np.ones(tb - ta), (np.arange(tb - ta), fbus[ta:tb])), shape=(tb - ta, self.nb))
            Ct = spr.csr_matrix((np.ones(tb - ta), (np.arange(tb - ta), tbus[ta:tb])), shape=(tb - ta, self.nb))
            Ytr = (Cf.T @ self.Yf[ta:tb] + Ct.T @ self.Yt[ta:tb]).toarray()
            Y = Y - Ytr
            self.hv = np.unique(fbus[ta:tb])
            self.slack = np.unique(tbus[ta:tb])
            setp = float(net.trafo.autoTapSetp.iloc[0]) if "autoTapSetp" in net.trafo else 1.025
            if not np.isfinite(setp):          # some SimBench scenario files omit the set point; use scenario 0's 1.025 p.u.
                setp = 1.025
            self.vS = np.full(len(self.slack), setp + 0j)
            self.n_trafo_at = {int(b): int((tbus[ta:tb] == b).sum()) for b in self.slack}
        self.Yprime = Y
        excl = set(self.slack.tolist()) | set(self.hv.tolist())
        self.n_idx = np.array([i for i in range(self.nb) if i not in excl])
        Ynn = Y[np.ix_(self.n_idx, self.n_idx)]
        self.YnS = Y[np.ix_(self.n_idx, self.slack)]
        self.Znn = np.linalg.inv(Ynn)
        # line current limits in pu of the from-bus base
        vn_from = ppci["bus"][fbus, 9].real                      # BASE_KV
        ln = net.line
        imax_ka = (ln.max_i_ka * ln.parallel * ln.df).values * ln.max_loading_percent.values / 100 * line_limit
        a, b = self.line_rng
        self.line_fbus = fbus[a:b]
        self.line_tbus = tbus[a:b]
        self.imax_line_pu = imax_ka * np.sqrt(3) * vn_from[a:b] / self.base
        ta, tb = self.trafo_rng
        tr = net.trafo
        self.smax_trafo_pu = (tr.sn_mva * tr.parallel * tr.df).values * tr.max_loading_percent.values / 100 / self.base
        self.trafo_fbus = fbus[ta:tb]
        self.trafo_tbus = tbus[ta:tb]
        # voltage limits per ppc bus (from SimBench bus data)
        self.vmin = np.full(self.nb, 0.0)
        self.vmax = np.full(self.nb, 9.0)
        for pb, row in net.bus.iterrows():
            i = self.lookup[pb]
            self.vmin[i] = max(self.vmin[i], row.min_vm_pu)
            self.vmax[i] = min(self.vmax[i], row.max_vm_pu)
        # HV buses are outside the model in OLTC mode
        self.vmin[self.hv] = 0.0
        self.vmax[self.hv] = 9.0
        if profiles:
            self._load_profiles()

    def _load_profiles(self):
        net = self.net
        prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
        self.T = prof[("load", "p_mw")].shape[0]
        S = np.zeros((self.nb, self.T), dtype=complex)
        lb = self.lookup[net.load.bus.values]
        P = prof[("load", "p_mw")].values.T
        Q = prof[("load", "q_mvar")].values.T
        np.add.at(S, lb, -(P + 1j * Q) / self.base)
        gb = self.lookup[net.sgen.bus.values]
        G = prof[("sgen", "p_mw")].values.T
        np.add.at(S, gb, G / self.base)
        if ("storage", "p_mw") in prof and prof[("storage", "p_mw")].shape[1] and len(net.storage):
            # pandapower convention: storage p_mw > 0 is charging (consumption)
            sbus = self.lookup[net.storage.bus.values]
            np.add.at(S, sbus, -prof[("storage", "p_mw")].values.T / self.base)
        self.S0 = S                                   # net injection, pu
        self.index = prof[("load", "p_mw")].index

    def study_case(self, case):
        """Injection matrix (nb x 1) for a SimBench study case, plus slack voltage."""
        net = self.net
        lc = net.loadcases.loc[case]
        S = np.zeros((self.nb, 1), dtype=complex)
        lb = self.lookup[net.load.bus.values]
        np.add.at(S[:, 0], lb, -(net.load.p_mw.values * lc.pload + 1j * net.load.q_mvar.values * lc.qload) / self.base)
        gb = self.lookup[net.sgen.bus.values]
        fac = np.where(net.sgen.type.str.contains("PV"), lc.PV_p,
                       np.where(net.sgen.type.str.contains("Wind"), lc.Wind_p, lc.RES_p))
        np.add.at(S[:, 0], gb, net.sgen.p_mw.values * fac / self.base)
        return S, lc.Slack_vm

    def solve(self, S, vs=None, V0=None, tol=1e-9, maxit=60):
        """Solve power flow for injection matrix S (nb x T). Returns V (nb x T), converged mask.
        vs (optional) overrides the slack voltage magnitude (fixed-tap mode only)."""
        vS = self.vS if vs is None else np.full(len(self.slack), vs + 0j)
        Sn = S[self.n_idx]
        Vn = np.ones_like(Sn) * vS[0] if V0 is None else V0[self.n_idx].copy()
        w = -(self.Znn @ (self.YnS @ vS))[:, None]
        conv = np.zeros(S.shape[1], dtype=bool)
        for it in range(maxit):
            Vnew = w + self.Znn @ np.conj(Sn / Vn)
            d = np.abs(Vnew - Vn).max(axis=0)
            Vn = Vnew
            if np.all(d < tol):
                break
        conv = d < 1e-6
        bad = ~np.isfinite(Vn).all(axis=0) | (np.abs(Vn).min(axis=0) < 0.5)
        conv &= ~bad
        V = np.empty((self.nb, S.shape[1]), dtype=complex)
        V[self.slack] = vS[:, None]
        V[self.hv] = 1.0
        V[self.n_idx] = Vn
        return V, conv

    def trafo_loading(self, V, S):
        """Transformer loading (fraction of limit), current-based as in pandapower. (ntrafo x T)"""
        ta, tb = self.trafo_rng
        if self.mode == "fixed":
            If = np.abs(self.Yf[ta:tb] @ V)
            It = np.abs(self.Yt[ta:tb] @ V)
            I = np.maximum(If, It)
        else:
            I = np.empty((tb - ta, V.shape[1]))
            for j in range(tb - ta):
                L = self.trafo_tbus[j]
                inj = self.Yprime[L] @ V - np.conj(S[L] / V[L])
                I[j] = np.abs(inj) / self.n_trafo_at[int(L)]
        return I / self.smax_trafo_pu[:, None]

    def margins(self, V, conv, S=None):
        """Per-constraint-class margins (>=0 feasible). Returns dict of (T,) arrays and argmins."""
        Vm = np.abs(V)
        lim = self.vmin > 0
        m_vlow = ((Vm[lim] - self.vmin[lim, None]) / 0.05).min(axis=0)
        m_vhigh = ((self.vmax[lim, None] - Vm[lim]) / 0.05).min(axis=0)
        If = self.Yf @ V
        It = self.Yt @ V
        a, b = self.line_rng
        Iline = np.maximum(np.abs(If[a:b]), np.abs(It[a:b]))
        m_line = (1 - Iline / self.imax_line_pu[:, None]).min(axis=0)
        m_trafo = (1 - self.trafo_loading(V, S)).min(axis=0)
        out = dict(vlow=m_vlow, vhigh=m_vhigh, line=m_line, trafo=m_trafo)
        for k in out:
            out[k] = np.where(conv, out[k], -1.0)
        return out


def combined_margin(mg, relax=None):
    keys = ["vlow", "vhigh", "line", "trafo"]
    M = np.vstack([mg[k] for k in keys])
    if relax is not None:
        # constraints violated in the base case may not be worsened beyond their base value
        M = M - np.minimum(relax, 0)[:, None] if relax.ndim == 1 else M - np.minimum(relax, 0)
    return M.min(axis=0), np.array(keys)[M.argmin(axis=0)]


def pmax_series(gm, bus_ppc, S0=None, vs=None, p_hi0=None, tol_mw=0.005, maxit=40, tan_phi=None):
    """Largest additional load P(t) [MW] at ppc bus with all constraints met, per time step.

    Uses a bracketed Illinois (modified regula-falsi) search on the combined margin.
    Constraints already violated in the base case are relaxed to their base-case value.
    """
    S0 = gm.S0 if S0 is None else S0
    T = S0.shape[1]
    e = np.zeros((gm.nb, 1))
    e[bus_ppc] = 1.0
    tp = TAN_PHI if tan_phi is None else tan_phi
    unit = -(1 + 1j * tp) / gm.base * e               # injection per MW of new load

    V0, c0 = gm.solve(S0, vs=vs)
    mg0 = gm.margins(V0, c0, S0)
    keys = ["vlow", "vhigh", "line", "trafo"]
    base = np.vstack([mg0[k] for k in keys])           # (4, T)
    relax = np.minimum(base, 0)                         # negative where base case violated

    def evalm(P, Vinit):
        Sx = S0 + unit * P[None, :]
        V, c = gm.solve(Sx, vs=vs, V0=Vinit)
        mg = gm.margins(V, c, Sx)
        M = np.vstack([mg[k] for k in keys]) - relax
        # a constraint violated before connection only needs to not get worse; allow for rounding
        # (with the MV busbar voltage held fixed, feeders other than the one with the new load are
        # electrically unaffected, so their margin is unchanged up to floating-point noise)
        M = M + np.where(relax < 0, 1e-7, 0.0)
        return M.min(axis=0), M.argmin(axis=0), V

    lo = np.zeros(T)
    f_lo, _, Vlo = evalm(lo, V0)
    hi = np.full(T, p_hi0 if p_hi0 else 5.0)
    f_hi, _, Vhi = evalm(hi, V0)
    # expand hi until infeasible everywhere
    for _ in range(12):
        need = f_hi >= 0
        if not need.any():
            break
        lo = np.where(need, hi, lo)
        f_lo = np.where(need, f_hi, f_lo)
        hi = np.where(need, hi * 2, hi)
        f_hi2, _, Vh2 = evalm(hi, Vhi)
        f_hi = np.where(need, f_hi2, f_hi)
        Vhi = np.where(need[None, :], Vh2, Vhi)
        Vlo = np.where(need[None, :], Vhi, Vlo)
    side = np.zeros(T, dtype=int)
    for it in range(maxit):
        width = hi - lo
        if width.max() < tol_mw:
            break
        # Illinois step, safeguarded with bisection when the bracket is wide or f is a divergence flag
        denom = f_lo - f_hi
        x = np.where(denom > 0, (lo * (-f_hi) + hi * f_lo) / np.where(denom > 0, denom, 1), 0.5 * (lo + hi))
        use_bis = (f_hi <= -0.999) | (it % 3 == 2)
        x = np.where(use_bis, 0.5 * (lo + hi), x)
        x = np.clip(x, lo + 0.05 * width, hi - 0.05 * width)
        fx, _, Vx = evalm(x, Vlo)
        pos = fx >= 0
        lo = np.where(pos, x, lo)
        f_lo_new = np.where(pos, fx, f_lo)
        f_hi_new = np.where(~pos, fx, f_hi)
        # Illinois modification: halve the retained endpoint value when the same side is kept twice
        f_hi_new = np.where(pos & (side == 1), f_hi_new * 0.5, f_hi_new)
        f_lo_new = np.where(~pos & (side == -1), f_lo_new * 0.5, f_lo_new)
        side = np.where(pos, 1, -1)
        hi = np.where(~pos, x, hi)
        f_lo, f_hi = f_lo_new, f_hi_new
        Vlo = np.where(pos[None, :], Vx, Vlo)
    _, arg, _ = evalm(lo, Vlo)
    binding = np.array(keys)[arg]
    return lo, binding, base


def flexible_hc(pmax, alpha):
    """Largest nameplate P <= max_t pmax(t) with curtailed energy <= alpha * P * T (all steps equal length)."""
    pmax = np.sort(pmax)
    if alpha <= 0:
        return pmax[0]
    T = len(pmax)
    lo, hi = pmax[0], pmax[-1]
    def frac(P):
        return np.maximum(P - pmax, 0).sum() / (P * T)
    if frac(hi) <= alpha:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if frac(mid) <= alpha:
            lo = mid
        else:
            hi = mid
    return lo


def curtailment_stats(pmax, P):
    """Event statistics for a load of nameplate P curtailed down to pmax(t)."""
    short = np.maximum(P - pmax, 0)
    on = short > 1e-9
    hours = on.sum() * DT_H
    # contiguous events
    d = np.diff(np.concatenate([[0], on.astype(int), [0]]))
    starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
    lens = (ends - starts) * DT_H
    return dict(
        energy_frac=short.sum() / (P * len(pmax)),
        hours=hours,
        n_events=len(lens),
        max_event_h=lens.max() if len(lens) else 0.0,
        med_event_h=float(np.median(lens)) if len(lens) else 0.0,
        single_share=float((lens <= DT_H).mean()) if len(lens) else 0.0,
        max_depth=(short.max() / P) if P > 0 else 0.0,
    )
