# Flexible-connection hosting capacity of large loads in SimBench MV grids

Code and results for the paper *Hosting Capacity Gains and Curtailment Burden of Flexible Connections for Large Medium-Voltage Loads*
(H. Y. J. Wang, submitted to IEEE PES General Meeting 2027).

For every medium-voltage (MV) bus of the four SimBench MV grids (scenario 0, switching variant), the code computes the
**headroom series**: the largest additional load (flat, cos phi = 0.95 lagging) that can be added in each 15-min
interval of the 2016 profile year without creating or worsening a violation of voltage, line, or transformer limits
(Eq. (2) of the paper). Firm, snapshot, energy-budget, and hour-budget hosting capacity, curtailment statistics,
sensitivities (tap model, line reserve, future scenario, power factor, N-1 line outages), imperfect-foresight tests,
and all tables and figures are derived from these series.

## Environment

Tested with Python 3.13.18 and the exact versions in `requirements.txt`:

```
python3.13 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

`pdflatex` and `bibtex` (TeX Live) are needed only to build the PDF. A clean reinstall in a fresh virtual
environment reproduced the stored headroom of 12 random buses within the 5-kW search tolerance
(`results/audit/clean_repro.json`, environment listings in `results/audit/environment_*.txt`).

## Result data

The headroom series (`results/*.npz`, about 420 MB) are attached to the repository's v1.0 release as
`results_data_part1_base_configs.zip` and `results_data_part2_N1_forecasts.zip`. Download both and unzip them in the
repository root; they extract into `results/`.

## Reproducing the paper

```
cd src
./reproduce.sh tables   # all numbers, tables and figures from the stored headroom series (minutes)
./reproduce.sh full     # recompute everything from SimBench first (several days on two cores; the N-1 runs
                        # checkpoint every 5 outages and resume when restarted)
```

Every number in the paper is a LaTeX macro written by one of the scripts below into `tex/numbers_*.tex`; tables are
written to `tex/tab_*.tex` and figures to `tex/figs/`.

## Layout

| Path | Content |
|---|---|
| `src/core.py` | Grid model from pandapower's internal case; vectorized Z-bus power flow over all 35,136 intervals (ideal-OLTC and fixed-tap modes); headroom search (bracketed Illinois, 5 kW); HC metrics and curtailment statistics |
| `src/run_all.py` | `python3 run_all.py <grid> <oltc\|fixed> [line_limit] [scenario] [pf] [n_random_buses]`: headroom series of all MV buses -> `results/<tag>.npz` |
| `src/n1_topo.py`, `src/n1.py` | N-1 line-outage sensitivity: contingency list with restoration switching; headroom as the minimum over the intact network and all contingencies |
| `src/n1_post.py` | Metadata for the N-1 results (used by `sensitivity.py`) |
| `src/forecast_headroom.py` | Net-load-based headroom forecasts (day-ahead and 15-min) |
| `src/robustness.py`, `src/robustness_oos.py`, `src/robustness_fc.py` | Imperfect foresight: persistence (in sample, alternating weeks, half-years) and stronger forecasts (seasonal-naive, seasonally adjusted persistence, net-load based; alternating weeks and blocked half-years) |
| `src/tod.py` | Fixed-schedule baseline (fixed curtailment window, no measurement) |
| `src/depthcap.py` | Hosting capacity with curtailment depth capped at 25 % of nameplate |
| `src/corrective.py` | Bounded corrective N-1 case (curtailment only after an outage; no outage statistics) |
| `src/figures.py`, `src/analyze.py` | Fig. 1, Fig. 2, Tables I and II, main number macros |
| `src/sensitivity.py` | Table IV (network assumptions incl. N-1) |
| `src/contract.py` | Contract-parameter table (Table III) |
| `src/stats.py` | Cluster-bootstrap confidence intervals (feeder groups within grids) and leave-one-grid-out check |
| `src/shapes.py`, `src/shapes_final.py` | Loads with a daily swing |
| `src/battery.py`, `src/pf_compare.py`, `src/extras.py`, `src/peaks.py` | Battery yardstick, unity power factor, copperplate reconciliation with Gu et al., rank correlations, peak intervals |
| `src/validate2.py`, `src/verify_boundary.py` | Solver vs pandapower Newton-Raphson; boundary spot check |
| `src/audit/` | Unit tests of the metric definitions; full-population boundary audit; pandapower check of N-1 headroom (`n1_ppcheck.py`); monotonicity scans (5 kW to 1.3x, 20 kW to 3x); independent pandapower bisection and scans; benchmark; macro collection |
| `results/` | Headroom series (`<grid>_<config>.npz`, float32 MW, buses x intervals) with metadata (`*_meta.json`), derived CSV tables, audit logs (`results/audit/`) |
| `tex/` | LaTeX source of the paper |

Grid names: `rural`, `semiurb`, `urban`, `comm` (SimBench codes `1-MV-<grid>--0-sw`). Configuration tags: `oltc`
(base case), `fixed` (fixed taps), `oltc_L60` (lines at 60 %), `oltc_S2` (scenario 2), `oltc_pf100` (unity power
factor, random subset), `oltc_N1` (N-1 line outages, planning voltage band), `oltc_N1w` (N-1, 0.90-1.10 p.u. after outages), `oltc_fcNL` (net-load-based forecasts).

## Modelling assumptions (see the paper for details)

Balanced positive-sequence model from pandapower; ideal OLTC holds each MV busbar at 1.025 p.u. (transformer impedance
and tap range not modelled); SimBench per-bus voltage bands (0.965-1.055 p.u. at MV); lines at 100 % of rating (N-0)
unless stated; transformers at the SimBench maximum loading; DG at unity power factor without voltage control; one new
connection at a time.

## License

Code: MIT (see `LICENSE`). Results and figures: CC BY 4.0. SimBench data are downloaded by the `simbench` package
under their own license and are not redistributed here.

## Citation

See `CITATION.cff`.

## AI use

The code was written by Claude (Anthropic) under the author's direction; see the paper's acknowledgment.
