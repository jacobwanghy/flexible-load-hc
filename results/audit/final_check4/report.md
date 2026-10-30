# Final check 4: independent recomputation of main.pdf

Sources: everything below was recomputed from `results/<grid>_<config>.npz`, `robustness_fc.csv`, `tod.csv`, `peaks.json`, `*_meta.json` and the audit JSONs, using `core.flexible_hc` and `curtailment_stats` only. No macro file was trusted. Scripts are in the session scratchpad (rc.py, rc2.py, rc3.py).

## Checked numbers

| Location | Paper | Recomputed | Status |
|---|---|---|---|
| Abstract / IV-B: 1% gain | 15.0% | 14.95 | MATCH |
| Abstract: hours at 1% | ~3350 h/yr | 3350.25 | MATCH |
| Abstract / IV-D: static winter gain | 3.1% | 3.10 | MATCH |
| Abstract: forecasts keep < half | best 6.7 (alt), 6.5 (blocked) vs 15.0 | 45% / 43% | MATCH |
| Abstract / IV-D / V: N-1 1% gain | 50.3% | 50.31 | MATCH |
| Abstract / Intro: 449 buses | 449 | 94+114+138+103 | MATCH |
| Table I peak load (MW) | 7.5/13.4/21.0/14.9 | 7.51/13.40/20.98/14.94 | MATCH |
| Table I peak net load (MW) | 6.3/12.1/19.8/13.4 | 6.30/12.08/19.82/13.37 | MATCH |
| Table I DG (MW) | 25.6/23.8/13.6/16.6 | 25.57/23.80/13.56/16.63 | MATCH |
| Table I HV/MV (MVA), Vn | 2x25, 2x40, 2x63, 2x40; 20/20/10/20 | 50/80/126/80; 20/20/10/20 | MATCH |
| IV-A firm median | 6.1 MW | 6.134 | MATCH |
| IV-A firm 95% interval | 5.85-6.61 | stats.json 5.854-6.609 (bootstrap not re-run) | MATCH (vs stats.json) |
| IV-A rural firm range | 1.5-42.8 MW | 1.484-42.78 | MATCH |
| IV-A median headroom over firm | 16.0% | 15.99 | MATCH |
| IV-A snapshot >1% overstated | 42% of buses, up to 12.5% | 41.6%, 12.53 | MATCH |
| IV-B 0.25% / 5% | 11.5 / 21.7 | 11.48 / 21.70 | MATCH |
| IV-B IQR at 1% | 10.9-19.9 | 10.90-19.91 | MATCH |
| IV-B fixed taps at 1% (new) | 16.7% | 16.70 | MATCH |
| IV-B 17 feeder groups, 3-5 per grid | 17 | 5+3+5+4 | MATCH |
| IV-B bootstrap 1% interval | 12.9-16.1 | stats.json 12.91-16.14 (not re-run) | MATCH (vs stats.json) |
| IV-B leave-one-grid-out | 14.7-15.3 | 14.67/14.95/15.29/15.32 | MATCH |
| IV-B eve share of curtailed bus-intervals | 38% (25% of intervals) | 38.0% (25.0%) | MATCH |
| IV-B Nov-Feb share | 43% (33%) | 43.5% (33.1%) | MATCH |
| IV-B added MW at 1% | 0.82 (0.69-0.96) | 0.825 (stats.json 0.686-0.963) | MATCH |
| IV-B 3350 h = 38% of year | 38% | 3350/8784 = 38.1% | MATCH |
| IV-B events / median length / depth | 1263 / 45 min / 13% | 1263 / 45 / 13.0 | MATCH |
| IV-B 100 h gain / depth | 6.8% / 6% | 6.82 / 6.38 | MATCH |
| IV-B 85 h gain | 6.6% | 6.60 | MATCH |
| IV-B 1% of intervals (Gu) | 6.7% | 6.66 | MATCH |
| IV-B copperplate | 185-258% | 184.7-258.0 (copperplate_reconcile.json) | MATCH |
| Table II firm (MW) | 5.62/6.18/5.77/7.22 | 5.62/6.18/5.77/7.22 | MATCH |
| Table II 0.25% | 9.8/11.0/11.6/13.1 | 9.85/10.96/11.60/13.06 | MATCH |
| Table II 0.5% | 11.4/12.4/13.2/14.7 | 11.43/12.44/13.20/14.74 | MATCH |
| Table II 1% | 13.3/14.2/15.1/16.8 | 13.26/14.21/15.09/16.84 | MATCH |
| Table II 100 h | 5.4/6.7/7.0/8.5 | 5.36/6.69/7.03/8.53 | MATCH |
| Table II hours | 3072/3518/3472/2945 | 3072/3518/3472.4/2945 | MATCH |
| Table II longest event | 19.0/20.8/18.5/18.8 | 19.0/20.75/18.5/18.75 | MATCH |
| Table III a=0.25%: firm, add, gain, h, ev, len, longest, depth, energy | 6.13, 0.62, 11.5, 1213, 926, 30, 15.0, 10, 0.25 | 6.134, 0.620, 11.48, 1212.75, 926, 30, 15.0, 10.3, 0.25 | MATCH |
| Table III a=1% | 6.13, 0.82, 15.0, 3350, 1263, 45, 19.0, 13, 1.00 | 6.134, 0.825, 14.95, 3350.25, 1263, 45, 19.0, 13.0, 1.00 | MATCH |
| Table III 100 h | 6.13, 0.34, 6.8, 100, 224, 15, 3.3, 6, 0.01 | 6.134, 0.337, 6.82, 100, 224, 15, 3.25, 6.4, 0.013 | MATCH (3.25 rounded half-up) |
| Table III N-1 a=1% | 2.88, 1.32, 50.3, 1536, 1060, 30, 15.8, 33, 1.00 | 2.877, 1.323, 50.31, 1535.5, 1060, 30, 15.75, 33.5, 1.00 | MATCH |
| Table IV base | 6.13, 11.5, 15.0, 6.8, 3350, 79/365/5 | 6.134, 11.48, 14.95, 6.82, 3350, 79/365/5 | MATCH |
| Table IV fixed taps | 5.89, 13.2, 16.7, 8.2, 3190, 106/341/2 | 5.894, 13.17, 16.70, 8.20, 3189.75, 106/341/2 | MATCH |
| Table IV N-1 planning band | 2.88, 38.7, 50.3, 29.4, 1536, 241/203/5 | 2.877, 38.66, 50.31, 29.44, 1535.5, 241/203/5 | MATCH |
| Table IV N-1 0.90-1.10 | 3.60, 29.2, 38.8, 20.2, 1781, 116/328/5 | 3.597, 29.19, 38.82, 20.23, 1780.75, 116/328/5 | MATCH |
| Table IV scenario 2 | 5.92, 14.1, 18.8, 8.1, 2493, 82/373/5 | 5.916, 14.07, 18.83, 8.11, 2492.75, raw BIND 82 vlow/372 line/5 trafo/1 vhigh | MATCH (one vhigh bus reassigned to line per footnote; reassignment itself NOT-CHECKABLE without power flow) |
| IV-C binding-group gains, OLTC | 15.7 (79) / 14.7 (365) | 15.65 / 14.69 | MATCH |
| IV-C fixed taps voltage vs line | 28.5 / 15.3 | 28.48 / 15.31 | MATCH |
| IV-C rank corr. Zth / firm | 0.01 / -0.18, signs differ by grid | 0.011 / -0.177; per grid Zth +,+,-,+; firm -,-,+,- | MATCH |
| IV-C ratio correlation | r = 0.98 | Pearson 0.982 | MATCH |
| IV-D S2 buses / gain | 460 / 18.8 | 460 / 18.83 | MATCH |
| IV-D unity pf (100 buses) | firm +5.7%, gain 15.1 -> 14.6 | n=100, +5.74%, 15.09 -> 14.57 | MATCH |
| IV-D N-1 outages / unsupplied | 444 / 20 | 444 / 20 (meta) | MATCH |
| IV-D N-1 firm 6.1 -> 2.9; 241 vs 79; 3.6 MW; 38.8%; 29.4% for 100 h; 1.32 vs 0.82 MW | all | 2.88, 241/79, 3.60, 38.82, 29.44, 1.323/0.825 | MATCH |
| IV-D "fewer hours, deeper curtailment" | - | 1536 < 3350 h; 33% > 13% | MATCH |
| IV-D halving non-urban transformer limits binds ~1 busbar bus per grid | - | no script/result found | NOT-CHECKABLE |
| IV-D alt weeks P15 / SA15 / NL15 | 6.7 / 5.7 / 5.7 | 6.69 / 5.69 / 5.69 | MATCH |
| IV-D alt day-ahead range, NL best | 2.3-3.2, NL best | 2.29 (SNW) - 3.23 (NLDA) | MATCH |
| IV-D "about half the 15-min gain" | - | 3.23/6.69 = 0.48 | MATCH |
| IV-D MAPE | 0.6-1.7% | 0.625-1.716 | MATCH |
| IV-D violation share 15-min / DA | 35-41% / 48-53% | 35.0-41.2 / 48.1-53.2 | MATCH |
| IV-D violating intervals | median 1-4 | 1-4 | MATCH |
| IV-D excess | 1.0-1.6% of nameplate | 1.035-1.564 | MATCH |
| IV-D blocked 15-min gain / viol (new) | 5.8-6.5 / 33-43% | 5.84-6.52 / 32.7-42.9 | MATCH |
| IV-D blocked day-ahead gain / viol (new) | 2.2-3.3 / 49-54% | 2.16-3.28 / 49.3-53.7 | MATCH |
| IV-D static winter gain / energy (new) | 3.1% / 1.0% | 3.10 / 0.99 | MATCH |
| IV-D static eve gain (new) | 0.2% | 0.22 | MATCH |
| IV-D static OOS violation (new) | 50% | 49.9% (448 buses 0.5, 1 bus 0.0) | MATCH (see problem 3) |
| V: N-1 50.3% (1.32 MW) | - | - | MATCH |
| V: "none of the forecasts ... kept more than 3.2%" | 3.2 | alt max DA 3.23, blocked max DA 3.28 | MISMATCH (see problem 1) |
| II-D verification: 15.8 M bus-intervals | 15.8 M | 15,776,064 | MATCH |
| II-D 295 cases / 55 overvoltage / 3.4 kW / 0.4 kW | - | 295 / 55 / 3.375 / 0.367 | MATCH |
| II-D 2103 bus-intervals / 4.1 M power flows | - | 470+570+552+511 = 2103; 4,128,517 steps | MATCH |
| II-D 60 pandapower scans; 60/60 N-1 cases | - | 60; 60/60 | MATCH |
| II-B |T| = 35136 (8784 h) | - | - | MATCH |
| Page count | <= 5 | 5 | MATCH |

## Problems found

1. **Discussion contradicts the blocked half-year result (Section V, 1st paragraph, last sentence).** It says "with day-ahead notice ... none of the forecasts tested here kept more than 3.2%". Section IV-D (blocked half-years) reports day-ahead gains of 2.2-3.3%: NLDA H1/H2 gives 3.28%. Fix by writing "about 3%" or "more than 3.3%", or by scoping the sentence to alternating weeks.
2. **"so alternating weeks do not flatter the forecasts" (IV-D, Imperfect foresight) is slightly overstated.** The best 15-min gain is 6.7% with alternating weeks and 6.5% blocked, so alternating weeks flatter the best forecast slightly. P15 violations are lower when blocked (32.7% vs 35.0%), and SA/NL gains are higher. "do not materially flatter" or "do not systematically flatter" would be accurate.
3. **Static baseline OOS violation of 50% is largely true by construction (IV-D, Static baseline).** For each bus, tod.py trains on one half and tests on the other. The direction that trains on the half with the higher out-of-window minimum always violates, so every bus without an exact tie scores 0.5 (448 of 449). The number matches. However, "it violated limits in the other in 50% of cases" reads like an empirical rate, when it is an artifact of the two-direction design. Consider "in one of the two directions at every bus but one" or similar.
4. **Rounding convention (Table III, 100 h, longest event): minor.** The exact median is 3.25, shown as 3.3 (half-up). Python's `%.1f` would give 3.2, so the table was rounded by a different rule than a plain format call. Not an error, noted for consistency only. 15.75 -> 15.8, 20.75 -> 20.8 and 18.75 -> 18.8 are consistent with half-up.
5. **Half-year split boundary: trivial.** Both robustness_fc.py and tod.py split at T//2 = 17568 intervals = 183 days, so "H1" runs Jan 1 to Jul 1 inclusive. The text says "January-June". This is a one-day difference with no effect on the numbers.
6. **NOT-CHECKABLE.** The claim "halving the others' limits would bind at about one busbar bus per grid" (end of IV-D) has no script or result file in results/ or results/audit/. The same applies to the urban 50% transformer maximum loading, which is not in the meta files.

## LaTeX / references
- main.log: no undefined references or citations, no LaTeX warnings, and no rerun request. The rendered text contains no "??". There are no TODO or placeholder strings in main.tex, ack.tex or tab_*.tex.
- Cite keys: 32 cited, 32 in refs.bib. None are missing and none are uncited.
- All \ref and \eqref targets resolve. Labels eq:zbus and sec:val are defined but never referenced (harmless). tab_shapes.tex (label tab:shapes) exists but is not \input. That is fine because the table is not in the paper.
- main.pdf (00:36:46) is newer than main.tex, ack.tex, the tables and the numbers files, so it is up to date.
- Page count: 5.

## Acknowledgment (tex/ack.tex)
It mentions only "the N-1 and forecasting analyses". Both remain in the paper (IV-D). It does not mention any removed analysis (battery, shapes, L60, etc.). OK.
