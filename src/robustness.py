"""Robustness of the flexible-HC gain to imperfect foresight (brief section 8).

Instead of following the true headroom Pbar(t), the flexible load sets its power in interval t from a
forecast F(t) with a safety margin m:   L(t) = min(P, (1 - m) F(t)).
Forecasts: 15-min persistence F(t) = Pbar(t-1) (curtailment based on the last measured interval) and
day-ahead persistence F(t) = Pbar(t-96) (schedule fixed from the same interval of the previous day).
A nameplate P is admissible if L(t) <= Pbar(t) in every interval (no limit violated) and the curtailed
energy sum(P - L(t)) stays within alpha of the load's annual energy. For each bus and forecast the best
margin m in 0..30% (1% steps) is chosen. Output: robustness table + macros.
"""
import json
import numpy as np, pandas as pd

RES = "../results"
GRIDS = ["rural", "semiurb", "urban", "comm"]
ALPHA = 0.01
MARGINS = np.round(np.arange(0, 0.31, 0.01), 2)


def best_P(pb, F, alpha):
    T = len(pb)
    firm = pb.min()
    best = (firm, 0.0)
    for m in MARGINS:
        cap = (1 - m) * F

        def ok(P):
            L = np.minimum(P, cap)
            return np.all(L <= pb + 1e-12) and (P - L).sum() <= alpha * P * T

        lo, hi = firm, pb.max() * 1.5
        if not ok(lo):
            continue
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if ok(mid) else (lo, mid)
        if lo > best[0]:
            best = (lo, m)
    return best


if __name__ == "__main__":
    rows = []
    for g in GRIDS:
        PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
        for i, pb in enumerate(PM):
            firm = pb.min()
            r = dict(grid=g, bus=i, firm=firm)
            for name, lag in [("perfect", 0), ("15min", 1), ("dayahead", 96)]:
                F = pb if lag == 0 else np.concatenate([pb[:lag], pb[:-lag]])   # first lag intervals: assume known
                P, m = best_P(pb, F, ALPHA) if lag else best_P(pb, pb, ALPHA)
                L = np.minimum(P, (1 - m) * F)
                r[f"gain_{name}"] = (P / firm - 1) * 100
                r[f"addMW_{name}"] = P - firm
                r[f"m_{name}"] = m * 100
                r[f"hrs_{name}"] = ((P - L) > 1e-9).sum() * 0.25
            rows.append(r)
        print(g, "done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/robustness.csv", index=False)
    mac = {}
    for name, tag in [("perfect", "Perf"), ("15min", "Fifteen"), ("dayahead", "Day")]:
        mac[f"rob{tag}gain"] = f"{df[f'gain_{name}'].median():.1f}"
        mac[f"rob{tag}addMW"] = f"{df[f'addMW_{name}'].median():.2f}"
        mac[f"rob{tag}m"] = f"{df[f'm_{name}'].median():.0f}"
        mac[f"rob{tag}hrs"] = f"{df[f'hrs_{name}'].median():.0f}"
        mac[f"rob{tag}qa"] = f"{df[f'gain_{name}'].quantile(0.25):.1f}"
        mac[f"rob{tag}qc"] = f"{df[f'gain_{name}'].quantile(0.75):.1f}"
        print(name, {k: v for k, v in mac.items() if k.startswith(f"rob{tag}")})
    with open("../tex/numbers_rob.tex", "w") as f:
        for k, v in sorted(mac.items()):
            f.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
