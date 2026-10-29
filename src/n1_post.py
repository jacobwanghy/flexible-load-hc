"""Write analysis metadata for the N-1 results (same bus list as the base case) and N-1 summary macros."""
import json, os
import numpy as np

RES = "../results"
GR = ["rural", "semiurb", "urban", "comm"]
for g, tag in [(g, t) for t in ["_oltc_N1", "_oltc_N1w"] for g in GR]:
    f = f"{RES}/{g}{tag}.npz"
    if not os.path.exists(f):
        print("missing", g); continue
    base = json.load(open(f"{RES}/{g}_oltc_meta.json"))
    run = json.load(open(f"{RES}/{g}{tag}_run.json")) if os.path.exists(f"{RES}/{g}{tag}_run.json") else \
        json.load(open(f"{RES}/{g}{tag}_meta.json"))
    if "log" in run:                                    # first call: keep the run log separately
        json.dump(run, open(f"{RES}/{g}{tag}_run.json", "w"), indent=1)
    base["contingency_mode"] = "N-1 line outages with restoration switching (n1.py)"
    base["n_contingencies"] = len(run["log"])
    base["n_unrestorable"] = sum(1 for x in run["log"] if x["lost"])
    json.dump(base, open(f"{RES}/{g}{tag}_meta.json", "w"), indent=1)
    PM1 = np.load(f)["PM"].astype(float); PM0 = np.load(f"{RES}/{g}_oltc.npz")["PM"].astype(float)
    W = np.load(f)["WORST"]
    print(g, tag, "contingencies", base["n_contingencies"], "unrestorable", base["n_unrestorable"],
          "firm N0 median %.2f  N1 median %.2f" % (np.median(PM0.min(1)), np.median(PM1.min(1))),
          "buses with N1 firm < N0 firm:", int((PM1.min(1) < PM0.min(1) - 0.005).sum()), "/", len(PM1),
          "zero firm:", int((PM1.min(1) <= 0.005).sum()))

# macros: number of contingencies and unrestorable (spur) outages over the four grids
runs = [json.load(open(f"{RES}/{g}_oltc_N1_run.json")) for g in GR if os.path.exists(f"{RES}/{g}_oltc_N1_run.json")]
if len(runs) == 4:
    n = sum(len(r["log"]) for r in runs); u = sum(sum(1 for x in r["log"] if x["lost"]) for r in runs)
    multi = sum(sum(1 for x in r["log"] if x["n_options"] > 1) for r in runs)
    sp = []
    for g in GR:
        P = np.load(f"{RES}/{g}_oltc_N1.npz")["PM"].astype(float); sp += list(np.median(P, 1) / P.min(1) - 1)
    open("../tex/numbers_nmo.tex", "w").write("\\newcommand{\\nmoN}{%d}\n\\newcommand{\\nmoUnr}{%d}\n\\newcommand{\\nmoMulti}{%d}\n\\newcommand{\\nmoSpread}{%.0f}\n" % (n, u, multi, np.median(sp) * 100))
    print("contingencies", n, "unrestorable", u, "several restoration options", multi)
