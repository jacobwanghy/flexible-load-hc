"""Collect the fine monotonicity scans into LaTeX macros (numbers_fs.tex, numbers_fsx.tex)."""
import json, glob, os
A = os.path.join(os.path.dirname(__file__), "..", "..", "results", "audit")
T = os.path.join(os.path.dirname(__file__), "..", "..", "tex")
fs = [json.load(open(f))["stats"] for f in sorted(glob.glob(f"{A}/finescan_*.json")) if "_x3" not in f]
fx = [json.load(open(f))["stats"] for f in sorted(glob.glob(f"{A}/finescan_*_x3.json"))]
print("1.3x:", sum(s["cases"] for s in fs), "cases,", sum(s["nonmonotone"] for s in fs), "non-monotone,",
      sum(s["nonconv"] for s in fs), "non-converged")
print("3x:", sum(s["cases"] for s in fx), "cases,", sum(s.get("nonmonotone", 0) for s in fx), "non-monotone")
open(f"{T}/numbers_fs.tex", "w").write("\\newcommand{\\fsCases}{%d}\n\\newcommand{\\fsSteps}{%.1f}\n" % (
    sum(s["cases"] for s in fs), sum(s["steps_total"] for s in fs) / 1e6))
open(f"{T}/numbers_fsx.tex", "w").write("\\newcommand{\\fsxCases}{%d}\n" % sum(s["cases"] for s in fx))
