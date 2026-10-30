"""Imperfect foresight with stronger non-oracle forecasts, alternating-weeks protocol.

Forecast F(t) of the headroom Pbar(t), using only information available at the decision time:
  15-min horizon (information up to t-1):
    P15   persistence                       Pbar(t-1)
    SA15  seasonally adjusted persistence   Pbar(t-1) + [B(t) - B(t-1)], B = day-type profile of Pbar (below)
    NL15  net-load based                    headroom computed from forecast injections s(t-1) + [s_bar(t) - s_bar(t-1)]
  day-ahead (information up to the previous midnight):
    PDA   persistence                       Pbar(t-96)
    SNW   seasonal naive, previous week     Pbar(t-672)
    SNDT  seasonal naive, day type          Pbar at the same time on the latest earlier day of the same type
                                            (Mon-Fri, Sat, Sun)
    NLDA  net-load based                    headroom from s_bar(t): mean injections at the same time on the last 4
                                            earlier days of the same type
B(t) is the mean of Pbar at the same time on the last 4 earlier days of the same type.
Control law and selection as in robustness.py: L(t) = min(P, (1-m) F(t)); P and m (0..30 %) are chosen on the training
weeks to maximize P with no violation (L <= Pbar) and curtailed energy <= alpha P |T_train| dt; they are then applied to
the test weeks. Both directions (even -> odd weeks, odd -> even). A case is one bus and one direction.
Blocked check: the same with training on January-June and testing on July-December, and the reverse.
Where no earlier data exist (start of the year) F = Pbar (optimistic, at most the first week).
Output: results/robustness_fc.csv, tex/numbers_fc.tex
"""
import numpy as np, pandas as pd
from robustness import best_P, ALPHA
from forecast_headroom import day_types, dayahead_profile, T_DAY

RES = "../results"
GRIDS = ["rural", "semiurb", "urban", "comm"]


def lag(x, k):
    return np.concatenate([x[:k], x[:-k]])


def same_type_prev(pb, T):
    typ = day_types(T); nd = T // T_DAY
    D = pb.reshape(nd, T_DAY); out = D.copy()
    for d in range(nd):
        prev = [e for e in range(d - 1, -1, -1) if typ[e] == typ[d]]
        if prev:
            out[d] = D[prev[0]]
    return out.reshape(T)


def forecasts(pb, nl):
    T = len(pb)
    B = dayahead_profile(pb[None, :])[0]
    F = dict(P15=lag(pb, 1), SA15=lag(pb, 1) + B - lag(B, 1), PDA=lag(pb, 96), SNW=lag(pb, 672),
             SNDT=same_type_prev(pb, T))
    if nl is not None:
        F["NL15"], F["NLDA"] = nl
    return F


def run():
    rows = []
    for g in GRIDS:
        PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
        try:
            z = np.load(f"{RES}/{g}_oltc_fcNL.npz"); NL = (z["PM_15"].astype(float), z["PM_DA"].astype(float))
        except FileNotFoundError:
            NL = None
        T = PM.shape[1]
        wk = (np.arange(T) // (96 * 7)) % 2 == 0
        h1 = np.arange(T) < T // 2                     # blocked split: Jan-Jun vs Jul-Dec
        for i, pb in enumerate(PM):
            firm = pb.min()
            Fs = forecasts(pb, None if NL is None else (NL[0][i], NL[1][i]))
            for name, F in Fs.items():
                for split, (tr, te) in {"EvenWkToOdd": (wk, ~wk), "OddWkToEven": (~wk, wk),
                                          "H1toH2": (h1, ~h1), "H2toH1": (~h1, h1)}.items():
                    P, m = best_P(pb[tr], F[tr], ALPHA)
                    L = np.minimum(P, (1 - m) * F[te])
                    viol = L > pb[te] + 1e-9
                    rows.append(dict(grid=g, bus=i, forecast=name, split=split, m=m * 100,
                                     gain=(P / firm - 1) * 100, viol_any=bool(viol.any()),
                                     viol_intervals=int(viol.sum()),
                                     max_excess_pct=float(((L - pb[te]) / P).max() * 100) if viol.any() else 0.0,
                                     test_energy_share=float((P - L).sum() / (P * L.size) * 100),
                                     mape=float(np.mean(np.abs(F - pb) / pb) * 100)))
        print(g, "done", "NL" if NL is not None else "no NL", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(f"{RES}/robustness_fc.csv", index=False)
    return d


if __name__ == "__main__":
    import sys
    dall = pd.read_csv(f"{RES}/robustness_fc.csv") if "--summary" in sys.argv else run()   # --summary: macros only
    blocked = dall[dall.split.str.startswith("H")]
    d = dall[~dall.split.str.startswith("H")]          # main protocol: alternating weeks
    s = d.groupby("forecast").agg(gain=("gain", "median"), gain_q1=("gain", lambda x: x.quantile(.25)),
                                  gain_q3=("gain", lambda x: x.quantile(.75)), viol=("viol_any", "mean"),
                                  vint=("viol_intervals", "median"), m=("m", "median"), mape=("mape", "median"),
                                  energy=("test_energy_share", "median"), n=("gain", "size"))
    s["viol"] *= 100
    v = d[d.viol_any]
    mac_extra = dict(fcMapeLo=f"{s.mape.min():.1f}", fcMapeHi=f"{s.mape.max():.1f}",
                     fcExcLo=f"{v.groupby('forecast').max_excess_pct.median().min():.1f}",
                     fcExcHi=f"{v.groupby('forecast').max_excess_pct.median().max():.1f}",
                     fcVintLo=f"{v.groupby('forecast').viol_intervals.median().min():.0f}",
                     fcVintHi=f"{v.groupby('forecast').viol_intervals.median().max():.0f}")
    print(s.round(2))
    mac = dict(mac_extra)
    for f in s.index:
        t = f.replace("15", "Fifteen").replace("DA", "Da").replace("DT", "Dt")
        mac[f"fc{t}gain"] = f"{s.loc[f, 'gain']:.1f}"
        mac[f"fc{t}viol"] = f"{s.loc[f, 'viol']:.0f}"
        mac[f"fc{t}m"] = f"{s.loc[f, 'm']:.0f}"
        mac[f"fc{t}mape"] = f"{s.loc[f, 'mape']:.1f}"
    H15, HDA = ["P15", "SA15", "NL15"], ["PDA", "SNW", "SNDT", "NLDA"]
    for hz, lst in [("Fifteen", H15), ("Da", HDA)]:
        lst = [f for f in lst if f in s.index]
        b = s.loc[lst, "gain"].idxmax()
        mac[f"fcBest{hz}gain"] = f"{s.loc[b, 'gain']:.1f}"
        mac[f"fcBest{hz}viol"] = f"{s.loc[b, 'viol']:.0f}"
        mac[f"fcRange{hz}gainLo"] = f"{s.loc[lst, 'gain'].min():.1f}"
        mac[f"fcRange{hz}gainHi"] = f"{s.loc[lst, 'gain'].max():.1f}"
        mac[f"fcRange{hz}violLo"] = f"{s.loc[lst, 'viol'].min():.0f}"
        mac[f"fcRange{hz}violHi"] = f"{s.loc[lst, 'viol'].max():.0f}"
        print(hz, "best", b)
    # blocked (half-year) split: same statistics
    sb_ = blocked.groupby("forecast").agg(gain=("gain", "median"), viol=("viol_any", "mean"))
    sb_["viol"] *= 100
    print("blocked half-year split\n", sb_.round(2))
    for hz, lst in [("Fifteen", H15), ("Da", HDA)]:
        lst = [f for f in lst if f in sb_.index]
        mac[f"fcBlk{hz}gainLo"] = f"{sb_.loc[lst, 'gain'].min():.1f}"
        mac[f"fcBlk{hz}gainHi"] = f"{sb_.loc[lst, 'gain'].max():.1f}"
        mac[f"fcBlk{hz}violLo"] = f"{sb_.loc[lst, 'viol'].min():.0f}"
        mac[f"fcBlk{hz}violHi"] = f"{sb_.loc[lst, 'viol'].max():.0f}"
    sb_.round(3).to_csv(f"{RES}/robustness_fc_blocked_summary.csv")
    mac["fcDaMaxAll"] = f"{max(s.loc[[x for x in HDA if x in s.index], 'gain'].max(), sb_.loc[[x for x in HDA if x in sb_.index], 'gain'].max()):.1f}"
    with open("../tex/numbers_fc.tex", "w") as fh:
        for k, v in sorted(mac.items()):
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
