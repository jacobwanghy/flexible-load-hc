"""Static time-of-day baseline (no real-time control).

The load runs at nameplate P outside a fixed window W and is held at its firm HC inside W, every year, with no
measurement or forecast. It is feasible if P <= Pbar(t) for every t outside W, so P = min_{t not in W} Pbar(t).
Windows: 'eve' 16:00-22:00 every day; 'wintereve' 16:00-22:00 in November-February; 'winter' all hours of
November-February. Reported per bus: gain P/firm - 1, restricted hours |W| (counted only if P > firm), and the energy
not served, (P - firm)|W| / (P |T|). Out of sample: P chosen on January-June, applied to July-December and the reverse;
a case violates if P exceeds the headroom in a test interval outside W.
Output: results/tod.csv, tex/numbers_tod.tex
"""
import numpy as np, pandas as pd

RES = "../results"
GRIDS = ["rural", "semiurb", "urban", "comm"]
T = 35136
idx = pd.date_range("2016-01-01", periods=T, freq="15min")
hr, mon = idx.hour.values, idx.month.values
eve = (hr >= 16) & (hr < 22)
win = np.isin(mon, [11, 12, 1, 2])
WINDOWS = {"eve": eve, "wintereve": eve & win, "winter": win}
h1 = np.arange(T) < T // 2

rows = []
for g in GRIDS:
    PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    for i, pb in enumerate(PM):
        firm = pb.min()
        for w, W in WINDOWS.items():
            P = max(pb[~W].min(), firm)
            en = (P - firm) * W.sum() / (P * T) * 100
            r = dict(grid=g, bus=i, window=w, gain=(P / firm - 1) * 100, hours=W.sum() * 0.25 if P > firm else 0.0,
                     energy=en)
            viol = []
            for tr, te in [(h1, ~h1), (~h1, h1)]:
                Ptr = max(pb[tr & ~W].min(), pb[tr].min())
                viol.append(bool((Ptr > pb[te & ~W] + 1e-9).any()))
            r["oos_viol_share"] = np.mean(viol)
            r["viol_h1_to_h2"], r["viol_h2_to_h1"] = viol
            rows.append(r)
d = pd.DataFrame(rows)
d.to_csv(f"{RES}/tod.csv", index=False)
s = d.groupby("window").agg(gain=("gain", "median"), hours=("hours", "median"), energy=("energy", "median"),
                            viol=("oos_viol_share", "mean"), v12=("viol_h1_to_h2", "mean"),
                            v21=("viol_h2_to_h1", "mean"))
s[["viol", "v12", "v21"]] *= 100
print(s.round(2))
mac = {}
for w, tag in [("eve", "Eve"), ("wintereve", "WinEve"), ("winter", "Win")]:
    mac[f"tod{tag}gain"] = f"{s.loc[w, 'gain']:.1f}"
    mac[f"tod{tag}hrs"] = f"{s.loc[w, 'hours']:.0f}"
    mac[f"tod{tag}en"] = f"{s.loc[w, 'energy']:.1f}"
    mac[f"tod{tag}viol"] = f"{s.loc[w, 'viol']:.0f}"
    mac[f"tod{tag}vOT"] = f"{s.loc[w, 'v12']:.0f}"
    mac[f"tod{tag}vTO"] = f"{s.loc[w, 'v21']:.0f}"
with open("../tex/numbers_tod.tex", "w") as f:
    for k, v in sorted(mac.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
