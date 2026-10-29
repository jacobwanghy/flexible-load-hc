"""Unit tests for the hosting-capacity metric definitions (brief 3.1-3.4).

Run: python3 test_metrics.py   (prints PASS/FAIL lines and writes ../../results/audit/metric_tests.json)
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from core import flexible_hc, curtailment_stats, DT_H
from analyze import hours_hc
from shapes import energy_hc, hour_hc, hours_curtailed

results = []


def check(name, cond, detail=""):
    results.append(dict(test=name, passed=bool(cond), detail=str(detail)))
    print("PASS" if cond else "FAIL", name, detail)


rng = np.random.default_rng(0)
T = 35136
pbar = 5 + rng.random(T)                     # continuous values, no ties

# --- hour budget, flat load: at most N hours curtailed, and maximal
for N in [0, 0.25, 1, 50, 100, 200]:
    P = hours_hc(pbar, N)
    st = curtailment_stats(pbar, P)
    k = int(round(N / DT_H))
    s = np.sort(pbar)
    # maximality: any P' > P curtails more than N hours (the next order statistic becomes curtailed)
    Pnext = s[k + 1] if k + 1 < T else np.inf
    st2 = curtailment_stats(pbar, (P + Pnext) / 2) if np.isfinite(Pnext) else None
    check(f"hour budget N={N}: curtailed hours <= N", st["hours"] <= N + 1e-9, f"hours={st['hours']}")
    check(f"hour budget N={N}: exactly N hours (no ties)", abs(st["hours"] - N) < 1e-9, f"hours={st['hours']}")
    if st2 is not None:
        check(f"hour budget N={N}: maximal (slightly larger P exceeds N)", st2["hours"] > N, f"hours={st2['hours']}")
check("hour budget N=0 equals firm HC", hours_hc(pbar, 0) == pbar.min())

# --- ties: equality means not curtailed
tied = np.array([1.0, 1.0, 1.0, 2.0, 3.0])
st = curtailment_stats(tied, 1.0)
check("tie: P equal to headroom is not curtailed", st["hours"] == 0, st)
P = hours_hc(tied, 0.25)                   # one interval allowed -> 2nd smallest = 1.0 (ties keep it at 1.0)
check("tie: hour budget with tied minimum stays at tied value", P == 1.0 and curtailment_stats(tied, P)["hours"] <= 0.25, P)

# --- energy budget, flat load: definition E(P) <= alpha * P * T * dt, and maximality
for a in [0.0025, 0.005, 0.01, 0.05]:
    P = flexible_hc(pbar, a)
    frac = np.maximum(P - pbar, 0).sum() / (P * T)
    frac_up = np.maximum(P * 1.0001 - pbar, 0).sum() / (P * 1.0001 * T)
    check(f"energy budget a={a}: within budget", frac <= a + 1e-12, f"frac={frac:.6f}")
    check(f"energy budget a={a}: maximal", frac_up > a or P >= pbar.max(), f"frac_up={frac_up:.6f}")

# --- monotonicity of E(P)/P (used for bisection)
Ps = np.linspace(pbar.min(), pbar.max() * 1.5, 200)
r = [np.maximum(P - pbar, 0).sum() / P for P in Ps]
check("E(P)/P non-decreasing in P", np.all(np.diff(r) >= -1e-12))

# --- shaped loads: denominator is the load's own energy alpha * P * sum f * dt
h = (np.arange(T) % 96) / 4
f = 1 - 0.2 / 2 * (1 - np.cos(2 * np.pi * (h - 3) / 24))
for a in [0.0025, 0.01]:
    P = energy_hc(pbar, f, a)
    E = DT_H * np.maximum(P * f - pbar, 0).sum()
    own = a * P * DT_H * f.sum()
    flatdenom = a * P * DT_H * T
    check(f"shaped a={a}: curtailed energy <= alpha * own energy", E <= own * (1 + 1e-9), f"E={E:.4f} own={own:.4f} (flat-style denominator would be {flatdenom:.4f})")
    Pup = P * 1.0001
    Eup = DT_H * np.maximum(Pup * f - pbar, 0).sum()
    check(f"shaped a={a}: maximal", Eup > a * Pup * DT_H * f.sum(), "")
check("flat shape reproduces flat energy budget", abs(energy_hc(pbar, np.ones(T), 0.01) - flexible_hc(pbar, 0.01)) < 1e-6)

# --- shaped hour budget uses transformed series pbar/f
ratio = pbar / f
for N in [0, 1, 100]:
    P = hour_hc(ratio, N)
    check(f"shaped hour budget N={N}: curtailed hours <= N", hours_curtailed(pbar, f, P) <= N + 1e-9, hours_curtailed(pbar, f, P))
check("shaped firm = min pbar/f is feasible everywhere", np.all(ratio.min() * f <= pbar + 1e-12))
check("shaped profile has no zero intervals (f>0)", f.min() > 0, f"min f={f.min():.2f}")

# --- event statistics consistency
x = np.array([5, 1, 1, 5, 5, 1, 5, 1, 1, 1], dtype=float)
st = curtailment_stats(x, 2.0)
check("events: 3 events, 6 curtailed intervals = 1.5 h, longest 0.75 h", st["n_events"] == 3 and abs(st["hours"] - 1.5) < 1e-9 and abs(st["max_event_h"] - 0.75) < 1e-9, st)
check("events: deepest curtailment (2-1)/2 = 50%", abs(st["max_depth"] - 0.5) < 1e-12, st["max_depth"])

json.dump(results, open(os.path.join(os.path.dirname(__file__), "..", "..", "results", "audit", "metric_tests.json"), "w"), indent=1)
print(f"{sum(r['passed'] for r in results)}/{len(results)} passed")
