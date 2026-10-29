"""Net-load-based headroom forecasts (imperfect-foresight baselines that use only past information).

The existing loads, DG and storage injections s(t) of every bus are forecast, and the headroom of Eq. (2) is then
computed from the forecast injections with the same power flow and search as the true headroom:
  day-ahead  ("NL-DA"):  s_hat(t) = mean of s at the same time of day on the last 4 days of the same day type
                         (Mon-Fri, Saturday, Sunday) before the day of t (information up to the previous midnight);
  15-min     ("NL-15"):  s_hat(t) = s(t-1) + [s_bar(t) - s_bar(t-1)], i.e. the last measured injections plus the
                         change between the two times in the day-ahead profile s_bar above.
Where fewer than 4 earlier days of the type exist, the available ones are used; on the first day of each type the
forecast equals the actual value (optimistic, at most 3 days of the year).
Usage: python3 forecast_headroom.py <grid>  -> results/<grid>_oltc_fcNL.npz with PM_DA, PM_15 (MW, buses x T)
"""
import sys, json, time, copy
import numpy as np, pandas as pd, scipy.sparse as spr, simbench as sb
from core import GridModel, GRIDS, pmax_series
from n1 import LocalModel, component

T_DAY = 96


def day_types(T):
    idx = pd.date_range("2016-01-01", periods=T, freq="15min")
    dow = idx.dayofweek.values[::T_DAY]
    return np.where(dow < 5, 0, np.where(dow == 5, 1, 2))


def dayahead_profile(S, k=4):
    """Mean injections over the last k earlier days of the same type (columns = intervals)."""
    nb, T = S.shape
    nd = T // T_DAY
    typ = day_types(T)
    D = S.reshape(nb, nd, T_DAY)
    out = np.empty_like(D)
    for d in range(nd):
        prev = [e for e in range(d - 1, -1, -1) if typ[e] == typ[d]][:k]
        out[:, d] = D[:, prev].mean(axis=1) if prev else D[:, d]
    return out.reshape(nb, T)


if __name__ == "__main__":
    grid = sys.argv[1]
    RES = "../results"
    gm = GridModel(GRIDS[grid], mode="oltc")
    S0 = gm.S0
    Sbar = dayahead_profile(S0)
    S15 = np.concatenate([S0[:, :1], S0[:, :-1]], axis=1) + Sbar - np.concatenate([Sbar[:, :1], Sbar[:, :-1]], axis=1)
    meta = json.load(open(f"{RES}/{grid}_oltc_meta.json"))["meta"]
    Yb = spr.csr_matrix(gm.Yprime != 0)
    stop = set(gm.slack.tolist()) | set(gm.hv.tolist())
    out = {}
    t0 = time.time()
    for name, Sh in [("DA", Sbar), ("15", S15)]:
        PM = np.zeros((len(meta), gm.T), dtype=np.float32)
        cache = {}
        Vf, _ = gm.solve(Sh)
        for i, m in enumerate(meta):
            b, p = m["bus"], m["ppc"]
            if p in stop:                              # candidate on an MV busbar: full model
                pm, _, _ = pmax_series(gm, p, S0=Sh)
            else:                                      # its feeder only (exact: feeders meet at slack busbars)
                comp = component(gm, Yb, p)
                if comp not in cache:
                    cache[comp] = LocalModel(gm, Sh, sorted(comp), V_full=Vf)
                lm = cache[comp]
                pm, _, _ = pmax_series(lm, lm.pos[p])
            PM[i] = pm
            print(grid, name, i + 1, len(meta), f"{time.time()-t0:.0f}s", flush=True)
        out[f"PM_{name}"] = PM
    np.savez_compressed(f"{RES}/{grid}_oltc_fcNL.npz", **out)
    print("DONE", grid, time.time() - t0)
