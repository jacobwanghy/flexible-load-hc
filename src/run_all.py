"""Compute per-bus, per-step additional-load headroom P_max(t) for a SimBench MV grid."""
import sys, time, json
import numpy as np
from core import GridModel, GRIDS, pmax_series

name = sys.argv[1]
mode = sys.argv[2] if len(sys.argv) > 2 else "oltc"
ll = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
scen = int(sys.argv[4]) if len(sys.argv) > 4 else 0
pf = float(sys.argv[5]) if len(sys.argv) > 5 else 0.95
bus_subset = int(sys.argv[6]) if len(sys.argv) > 6 else 0      # >0: random subset of this many buses
tag = f"{name}_{mode}" + ("" if ll == 1.0 else f"_L{int(ll*100)}") + ("" if scen == 0 else f"_S{scen}") + ("" if pf == 0.95 else f"_pf{int(pf*100)}")
code = GRIDS[name].replace("--0-", f"--{scen}-")
gm = GridModel(code, mode=mode, line_limit=ll)
tanphi = float(np.tan(np.arccos(pf)))
net = gm.net
mv = net.bus.index[net.bus.vn_kv < 50]
seen, buses = set(), []
for b in mv:
    p = int(gm.lookup[b])
    if p not in seen:
        seen.add(p); buses.append((int(b), p))

# base-case statistics
V0, c0 = gm.solve(gm.S0)
mg0 = gm.margins(V0, c0, gm.S0)
base_viol = {k: int((v < 0).sum()) for k, v in mg0.items()}


zth = np.abs(np.diag(gm.Znn))
nidx = {int(p): i for i, p in enumerate(gm.n_idx)}

PM = np.zeros((len(buses), gm.T), dtype=np.float32)
BIND = np.zeros((len(buses), gm.T), dtype=np.int8)
snap, snap_bind, meta = [], [], []
keys = ["vlow", "vhigh", "line", "trafo"]
t0 = time.time()
if bus_subset > 0:
    rs = np.random.default_rng(123)
    buses = [buses[j] for j in sorted(rs.choice(len(buses), min(bus_subset, len(buses)), replace=False))]
PM = np.zeros((len(buses), gm.T), dtype=np.float32)
BIND = np.zeros((len(buses), gm.T), dtype=np.int8)
for i, (b, p) in enumerate(buses):
    pm, bind, _ = pmax_series(gm, p, tan_phi=tanphi)
    PM[i] = pm
    BIND[i] = np.searchsorted(keys, bind) if False else np.array([keys.index(x) for x in bind], dtype=np.int8)
    meta.append(dict(bus=b, ppc=p, name=str(net.bus.name[b]), zth_pu=float(zth[nidx[p]]) if p in nidx else 0.0,
                     subnet=str(net.bus.subnet[b]) if 'subnet' in net.bus else ''))
    print(name, i + 1, len(buses), f"{time.time()-t0:.0f}s", f"firm={pm.min():.3f}", flush=True)

np.savez_compressed(f"../results/{tag}.npz", PM=PM, BIND=BIND)
json.dump(dict(grid=name, code=code, scenario=scen, pf=pf, mode=mode, line_limit=ll, meta=meta, base_viol=base_viol, T=gm.T,
               n_bus=len(net.bus), n_line=len(net.line), load_mw=float(net.load.p_mw.sum()),
               sgen_mw=float(net.sgen.p_mw.sum()), trafo_mva=float((net.trafo.sn_mva).sum()),
               peak_load_prof=float((-gm.S0.real).clip(min=0).sum(0).max()),
               vn_kv=float(net.bus.vn_kv[mv].iloc[0])),
          open(f"../results/{tag}_meta.json", "w"), indent=1)
print("DONE", tag, time.time() - t0)
