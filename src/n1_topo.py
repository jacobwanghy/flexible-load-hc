"""Contingency enumeration for the N-1 line-outage sensitivity.

For each SimBench MV grid: ring groups (connected parts of the MV network once all normally open line switches are
closed and the MV busbars are removed); for every in-service line, the buses it disconnects when it fails and the
normally open tie line(s) whose closing restores all of them.
"""
import numpy as np, networkx as nx, pandapower.topology as top, simbench as sb, pandapower as pp


def contingencies(net):
    sw = net.switch
    tie_sw = sw[(~sw.closed) & (sw.et == "l")]
    ties = sorted(set(tie_sw.element.astype(int)))
    pp.runpp(net, numba=False)
    lk = net._pd2ppc_lookups["bus"]
    bb_ppc = set(lk[net.trafo.lv_bus.values])
    mv = [b for b in net.bus.index if net.bus.vn_kv[b] < 50]
    busbars = [b for b in mv if lk[b] in bb_ppc]
    out = []
    for l in net.line.index:
        if l in ties or not net.line.in_service[l]:
            continue
        net.line.at[l, "in_service"] = False
        G = top.create_nxgraph(net, respect_switches=True, include_trafos=False)
        reach = set()
        for b in busbars:
            reach |= nx.node_connected_component(G, b)
        D = sorted(set(mv) - reach)
        restore = []
        if D:
            for tl in ties:
                idx = tie_sw.index[tie_sw.element == tl]
                net.switch.loc[idx, "closed"] = True
                G2 = top.create_nxgraph(net, respect_switches=True, include_trafos=False)
                r2 = set()
                for b in busbars:
                    r2 |= nx.node_connected_component(G2, b)
                if set(D) <= r2:
                    restore.append(tl)
                net.switch.loc[idx, "closed"] = False
        net.line.at[l, "in_service"] = True
        out.append(dict(line=int(l), lost=[int(x) for x in D], ties=restore))
    return out, ties


if __name__ == "__main__":
    import json, sys
    from core import GRIDS
    res = {}
    for g, c in GRIDS.items():
        net = sb.get_simbench_net(c)
        cs, ties = contingencies(net)
        n_none = sum(1 for x in cs if not x["lost"])
        n_rest = sum(1 for x in cs if x["lost"] and x["ties"])
        n_unr = sum(1 for x in cs if x["lost"] and not x["ties"])
        multi = sum(1 for x in cs if len(x["ties"]) > 1)
        print(g, "lines", len(cs), "no-loss", n_none, "restorable", n_rest, "unrestorable", n_unr, "multi-tie", multi,
              "lost buses (unrestorable, max)", max([len(x["lost"]) for x in cs if not x["ties"]] or [0]))
        res[g] = dict(contingencies=cs, ties=ties)
    json.dump(res, open("../results/n1_contingencies.json", "w"), indent=1)
