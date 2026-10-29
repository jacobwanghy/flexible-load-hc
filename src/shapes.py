"""Load-shape extension: hosting capacity for new loads that follow a daily/weekly profile f(t).

The stored headroom series Pbar_k(t) does not depend on the new load's shape (power factor fixed),
so a load P*f(t) (max f = 1, P = nameplate peak) is feasible in interval t iff P*f(t) <= Pbar_k(t).
"""
import json, os, sys
import numpy as np, pandas as pd

LIB = "/usr/local/lib/python3.13/site-packages/simbench/networks/1-complete_data-mixed-all-0-sw/LoadProfile.csv"
RES = "../results"
GRIDS = ["rural", "semiurb", "urban", "comm"]
DT = 0.25


def load_shapes(names, smooth=True):
    """Typical-day shapes: SimBench profiles averaged over (month, day type, 15-min slot), scaled to max 1."""
    p = pd.read_csv(LIB, sep=";")
    t = pd.to_datetime(p.time, format="%d.%m.%Y %H:%M")
    dtype = np.select([t.dt.dayofweek < 5, t.dt.dayofweek == 5], [0, 1], 2)
    key = [t.dt.month.values, dtype, (t.dt.hour * 4 + t.dt.minute // 15).values]
    out = {"Flat": np.ones(len(p))}
    h = (t.dt.hour + t.dt.minute / 60).values
    for n in names:
        if n.startswith("DC"):          # DC<amp%>@<peak hour>: smooth daily swing, max 1 at the peak hour
            amp, hp = n[2:].split("@")
            a = float(amp) / 100
            out[n] = 1 - a / 2 * (1 - np.cos(2 * np.pi * (h - float(hp)) / 24))
            continue
        x = pd.Series(p[f"{n}_pload"].values.astype(float))
        if smooth:
            x = x.groupby(key).transform("mean")
        x = x.values
        out[n] = x / x.max()
    return out, t


def firm(r):
    return r.min()


def hour_hc(r, H):
    s = np.sort(r)
    return s[min(int(round(H / DT)), len(s) - 1)]


def energy_hc(pbar, f, alpha):
    m = f > 0
    r = pbar[m] / f[m]
    lo, hi = r.min(), r.max()
    F = f.sum()
    frac = lambda P: np.maximum(P * f - pbar, 0).sum() / (P * F)
    if frac(hi) <= alpha:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if frac(mid) <= alpha else (lo, mid)
    return lo


def hours_curtailed(pbar, f, P):
    return ((P * f - pbar) > 1e-9).sum() * DT


if __name__ == "__main__":
    names = sys.argv[1:] or ["G1-A", "G2-A", "G3-A"]
    shapes, t = load_shapes(names)
    rows = []
    for g in GRIDS:
        PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
        for i in range(PM.shape[0]):
            pb = PM[i]
            for n, f in shapes.items():
                m = f > 0
                r = pb[m] / f[m]
                H0 = firm(r)
                e1 = energy_hc(pb, f, 0.01)
                rows.append(dict(grid=g, bus=i, shape=n, lf=f.mean(), firm=H0, e025=energy_hc(pb, f, 0.0025), e1=e1,
                                 h100=hour_hc(r, 100), hrs1=hours_curtailed(pb, f, e1)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/shapes.csv", index=False)
    flat = df[df["shape"] == "Flat"].set_index(["grid", "bus"])
    for n in shapes:
        d = df[df["shape"] == n].set_index(["grid", "bus"])
        print(f"{n:6s} LF {d.lf.iloc[0]:.2f} | firm peak med {d.firm.median():.2f} MW, vs flat {np.median(d.firm/flat.firm)*100-100:+.1f}% | "
              f"firm energy (P*LF) vs flat {np.median(d.firm*d.lf/flat.firm)*100-100:+.1f}% | gain 0.25% {np.median(d.e025/d.firm-1)*100:.1f} "
              f"1% {np.median(d.e1/d.firm-1)*100:.1f} 100h {np.median(d.h100/d.firm-1)*100:.1f} | hrs@1% {d.hrs1.median():.0f}")
