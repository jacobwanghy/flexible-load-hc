"""Runtime benchmark (brief section 4): model setup, one full-year power flow, and the headroom search per bus."""
import sys, os, time, json, platform
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np
from core import GridModel, GRIDS, pmax_series
out = {"cpu": platform.processor() or platform.machine(), "threads": os.environ.get("OMP_NUM_THREADS")}
try:
    out["cpu_model"] = [l.split(":")[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")][0]
    out["cores"] = os.cpu_count()
except Exception:
    pass
for g in ["rural", "urban"]:
    t0 = time.time(); gm = GridModel(GRIDS[g], mode="oltc"); t_setup = time.time() - t0
    t0 = time.time(); V, c = gm.solve(gm.S0); t_solve = time.time() - t0
    mv = gm.net.bus.index[gm.net.bus.vn_kv < 50][:5]
    ts = []
    for b in mv:
        t0 = time.time(); pmax_series(gm, gm.lookup[b]); ts.append(time.time() - t0)
    out[g] = dict(setup_s=round(t_setup, 1), full_year_solve_s=round(t_solve, 2), search_per_bus_s=round(float(np.mean(ts)), 1),
                  note="setup includes SimBench load and profile assembly; search excludes I/O and post-processing")
    print(g, out[g], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "..", "results", "audit", "benchmark.json"), "w"), indent=1)
print(out)
