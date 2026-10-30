"""Rough check behind the sentence on transformer N-1 limits (Section IV-D).

For the rural, semi-urban and commercial grids (transformers at 100 % in SimBench), an N-1 transformer limit would be
half the installed capacity. At the grid's peak net load the remaining active-power room for new load is about
0.5 * capacity - peak net load (reactive power and losses ignored). The script counts buses whose firm HC (N-0 and N-1)
exceeds that bound minus 1 MW, i.e., where the halved limit could bind. Output: results/audit/trafo_n1_bound.json
"""
import json, os
import numpy as np

RES = os.path.join(os.path.dirname(__file__), "..", "..", "results")
out = {}
for g in ["rural", "semiurb", "comm"]:
    m = json.load(open(f"{RES}/{g}_oltc_meta.json"))
    peaks = json.load(open(f"{RES}/peaks.json"))[g]
    bound = 0.5 * m["trafo_mva"] - peaks["peak_net"]
    r = dict(trafo_mva=m["trafo_mva"], peak_net_MW=peaks["peak_net"], bound_MW=bound)
    for tag in ["_oltc", "_oltc_N1"]:
        P = np.load(f"{RES}/{g}{tag}.npz")["PM"].astype(float)
        r[f"buses_near_bound{tag}"] = int((P.min(1) > bound - 1).sum())
    r["n_buses"] = len(m["meta"])
    out[g] = r
    print(g, r)
json.dump(out, open(f"{RES}/audit/trafo_n1_bound.json", "w"), indent=1)
