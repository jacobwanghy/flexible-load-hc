"""Out-of-sample check of the imperfect-foresight result.

For each bus and forecast (15-min and day-ahead persistence), the nameplate P and margin m are chosen on one
half-year (training) exactly as in robustness.py, then applied unchanged to the other half (test).
Reports the share of buses with any limit violation in the test half, the test-half curtailed-energy share,
and the gain. Both splits are run (Jan-Jun -> Jul-Dec and Jul-Dec -> Jan-Jun).
"""
import numpy as np, pandas as pd
from robustness import best_P, ALPHA

RES = "../results"
GRIDS = ["rural", "semiurb", "urban", "comm"]
rows = []
for g in GRIDS:
    PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    T = PM.shape[1]; half = T // 2
    for i, pb in enumerate(PM):
        firm_year = pb.min()
        for name, lag in [("15min", 1), ("dayahead", 96)]:
            F = np.concatenate([pb[:lag], pb[:-lag]])
            wk = (np.arange(T) // (96 * 7)) % 2 == 0
            for split, (tr, te) in {"H1toH2": (slice(0, half), slice(half, T)), "H2toH1": (slice(half, T), slice(0, half)),
                                    "EvenWkToOdd": (wk, ~wk), "OddWkToEven": (~wk, wk)}.items():
                P, m = best_P(pb[tr], F[tr], ALPHA)
                L = np.minimum(P, (1 - m) * F[te])
                viol = (L > pb[te] + 1e-9)
                rows.append(dict(grid=g, bus=i, forecast=name, split=split, m=m * 100,
                                 gain_vs_year_firm=(P / firm_year - 1) * 100,
                                 test_violation_intervals=int(viol.sum()), test_any_violation=bool(viol.any()),
                                 test_max_excess_MW=float(np.max(L - pb[te])) if viol.any() else 0.0,
                                 test_energy_share=float((P - L).sum() / (P * L.size)) * 100))
    print(g, "done", flush=True)
d = pd.DataFrame(rows)
d.to_csv(f"{RES}/robustness_oos.csv", index=False)
d["design"] = np.where(d.split.str.contains("Wk"), "weeks", "halves")
print(d.groupby(["design", "forecast"]).agg(viol=("test_any_violation", "mean"), gain=("gain_vs_year_firm", "median"),
      vint=("test_violation_intervals", "median")).round(3))
s = d[d.design == "weeks"].groupby("forecast").agg(share_buses_violating=("test_any_violation", "mean"),
                              median_violation_intervals=("test_violation_intervals", "median"),
                              median_gain=("gain_vs_year_firm", "median"),
                              median_energy_share=("test_energy_share", "median"))
s["share_buses_violating"] *= 100
print(s.round(2))
mac = {}
for f, tag in [("15min", "Fifteen"), ("dayahead", "Day")]:
    mac[f"oos{tag}viol"] = f"{s.loc[f, 'share_buses_violating']:.0f}"
    mac[f"oos{tag}gain"] = f"{s.loc[f, 'median_gain']:.1f}"
    h = d[(d.design == "halves") & (d.forecast == f)]
    mac[f"oosH{tag}viol"] = f"{h.test_any_violation.mean() * 100:.0f}"
with open("../tex/numbers_oos.tex", "w") as fh:
    for k, v in sorted(mac.items()):
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
