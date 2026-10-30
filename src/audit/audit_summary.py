"""Collect audit results into LaTeX macros (numbers_audit.tex) and a JSON summary."""
import json, glob, os
import numpy as np, pandas as pd

RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
A = f"{RES}/audit"
b = [json.load(open(f))["summary"] for f in sorted(glob.glob(f"{A}/boundary_*_oltc.json"))]
m = [json.load(open(f))["stats"] for f in sorted(glob.glob(f"{A}/monotone_*.json"))]
v = pd.read_csv(f"{A}/validate_headroom_all.csv") if os.path.exists(f"{A}/validate_headroom_all.csv") else \
    pd.DataFrame(sum([json.load(open(f)) for f in sorted(glob.glob(f"{A}/validate_headroom_*.json"))], []))
tests = json.load(open(f"{A}/metric_tests.json"))
S = dict(
    boundary_grids=len(b), boundary_intervals=sum(x["intervals"] for x in b),
    boundary_infeasible_at=sum(x["infeasible_at"] for x in b), boundary_feasible_above=sum(x["feasible_above"] for x in b),
    nonconv=sum(x["nonconv_at"] + x["nonconv_above"] for x in b), overvoltage_worsened=sum(x["overvoltage_worsened"] for x in b),
    monotone_grids=len(m), monotone_cases=sum(x["cases"] for x in m), noncontiguous=sum(x["noncontiguous"] for x in m),
    monotone_nonconv=sum(x["nonconv"] for x in m),
    nonmonotone_line=sum(x["nonmonotone"]["line"] for x in m), nonmonotone_trafo=sum(x["nonmonotone"]["trafo"] for x in m),
    val_cases=len(v), val_ovt=int(v.preexisting_overvoltage.sum()), val_err_max=float(v.err_kW.abs().max()),
    val_err_med=float(v.err_kW.abs().median()), val_err_p95=float(v.err_kW.abs().quantile(0.95)),
    val_false_feasible=int((~v.feasible_at_ours.astype(bool)).sum()), val_false_infeasible=int((~v.infeasible_above_ours.astype(bool)).sum()),
    tests_passed=sum(t["passed"] for t in tests), tests_total=len(tests))
json.dump(S, open(f"{A}/audit_summary.json", "w"), indent=1)
mac = dict(audBoundaryM=f"{S['boundary_intervals']/1e6:.1f}", audMonoK=f"{S['monotone_cases']/1e3:.0f}",
           audValN=str(S["val_cases"]), audValOvt=str(S["val_ovt"]), audValMax=f"{S['val_err_max']:.1f}",
           audValMed=f"{S['val_err_med']:.1f}")
# like-for-like comparison with Gu et al.: hour budget of 1% of intervals (351 of 35136)
g1pi = []
for g in ["rural", "semiurb", "urban", "comm"]:
    PM = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    k = int(0.01 * PM.shape[1]); srt = np.sort(PM, axis=1)
    g1pi += list(srt[:, k] / srt[:, 0] - 1)
n1 = sum([json.load(open(f)) for f in sorted(glob.glob(f"{A}/n1_ppcheck_*.json"))], [])
mac["nOnePPN"] = str(len(n1)); mac["nOnePPOK"] = str(sum(x["feasible_at"] and x["infeasible_above"] for x in n1))
mac["audTestsN"] = str(S["tests_total"])
mac["gonepctint"] = f"{np.median(g1pi) * 100:.1f}"
with open(os.path.join(os.path.dirname(__file__), "..", "..", "tex", "numbers_audit.tex"), "w") as f:
    for k, x in sorted(mac.items()):
        f.write(f"\\newcommand{{\\{k}}}{{{x}}}\n")
print(json.dumps(S, indent=1))
