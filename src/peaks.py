"""Aggregate peak-load and peak-net-load intervals per grid (snapshot planning interval)."""
import json
import simbench as sb
from core import GRIDS

out = {}
for n, c in GRIDS.items():
    net = sb.get_simbench_net(c)
    prof = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
    L = prof[("load", "p_mw")].values.sum(1)
    G = prof[("sgen", "p_mw")].values.sum(1)
    out[n] = dict(t_peak_load=int(L.argmax()), t_peak_net=int((L - G).argmax()), peak_load=float(L.max()),
                  peak_net=float((L - G).max()), min_net=float((L - G).min()),
                  energy_load=float(L.sum() * 0.25), energy_gen=float(G.sum() * 0.25))
    print(n, out[n])
json.dump(out, open("../results/peaks.json", "w"), indent=1)
