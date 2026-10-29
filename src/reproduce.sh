#!/usr/bin/env bash
# Regenerate the paper's numbers, tables and figures.
#
#   ./reproduce.sh tables   from the stored headroom series in ../results (minutes):
#                           every LaTeX number macro, every table, every figure, then the PDF if pdflatex is present
#   ./reproduce.sh full     recompute everything from SimBench (about 2-3 days on two cores), then "tables"
#
# Run from this directory. Requires the pinned environment in requirements.txt (Python 3.13).
set -euo pipefail
cd "$(dirname "$0")"
GRIDS="rural semiurb urban comm"
MODE="${1:-tables}"

if [ "$MODE" = "full" ]; then
  python3 validate2.py                       # solver vs pandapower Newton-Raphson
  python3 peaks.py
  for g in $GRIDS; do
    python3 run_all.py $g oltc               # base case (ideal OLTC)
    python3 run_all.py $g fixed              # fixed taps
    python3 run_all.py $g oltc 0.6           # lines limited to 60 %
    python3 run_all.py $g oltc 1.0 2         # SimBench scenario 2
    python3 run_all.py $g oltc 1.0 0 1.0 25  # unity power factor, 25 random buses per grid
  done
  python3 n1_topo.py                         # contingency list
  for g in $GRIDS; do
    python3 n1.py $g plan                    # N-1 line outages with restoration switching
    python3 n1.py $g wide                    # same, 0.90-1.10 p.u. after outages
    python3 forecast_headroom.py $g          # net-load-based headroom forecasts
  done
  ( cd audit
    python3 test_metrics.py
    for g in $GRIDS; do
      python3 audit_boundary.py $g; python3 audit_boundary.py $g oltc 1.0 2
      python3 audit_monotone.py $g; python3 validate_headroom.py $g
      python3 audit_finescan.py $g; python3 audit_finescan_x.py $g 3.0 0.02; python3 pp_scan.py $g
      python3 n1_ppcheck.py $g 15
    done
    python3 benchmark.py )
fi

# ---- numbers, tables and figures from stored results ----
python3 figures.py            # Fig. 1, Fig. 2, Table I, Table II, numbers.tex
python3 n1_post.py            # N-1 metadata for the sensitivity table
python3 sensitivity.py        # Table IV, numbers_sens.tex
python3 robustness.py         # in-sample foresight, numbers_rob.tex
python3 robustness_oos.py     # persistence, alternating weeks and half-years, numbers_oos.tex
python3 robustness_fc.py      # stronger forecasts, alternating weeks, numbers_fc.tex
python3 shapes_final.py       # daily-shape cases, numbers_shapes.tex
python3 battery.py            # numbers_batt.tex
python3 pf_compare.py         # numbers_pf.tex
python3 extras.py             # copperplate reconciliation, rank correlations, numbers_cp.tex, numbers_rank.tex
python3 stats.py              # cluster-bootstrap intervals, leave-one-grid-out, numbers_stats.tex
python3 contract.py           # Table III (contract parameters), numbers_contract.tex
( cd audit && python3 audit_summary.py && python3 finescan_summary.py )

if command -v pdflatex >/dev/null; then
  ( cd ../tex && pdflatex -interaction=nonstopmode main >/dev/null && bibtex main >/dev/null \
      && pdflatex -interaction=nonstopmode main >/dev/null && pdflatex -interaction=nonstopmode main >/dev/null \
      && echo "PDF: tex/main.pdf" )
fi
echo "done ($MODE)"
