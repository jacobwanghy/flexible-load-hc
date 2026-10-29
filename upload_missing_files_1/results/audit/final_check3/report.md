# Final check 3: independent number and consistency audit

Basis: main.tex and main.pdf as of 2026-10-29 22:46 (the version with the battery sentence removed and fani2024, zakaria2025 and umoh2023 restored). The text was re-extracted with `pdftotext -layout` after the update.
How the values were recomputed: from the raw `results/<grid>_<config>.npz` arrays, using my own exact piecewise-linear solver for H(alpha). It agrees with `core.flexible_hc` to 1e-15. The hour-budget value is sorted[4N]. Curtailment statistics count intervals with PM < P. All medians are pooled over buses with plain numpy medians. The macro files were not used as inputs.

## 1. Numbers checked

| Claim (location) | Paper | Recomputed | Source | Status |
|---|---|---|---|---|
| Buses (abstract, Sec. I, III, Fig. 1) | 449 | 94+114+138+103 = 449 npz rows | \nBus | MATCH |
| Pooled 1% gain (abstract, IV-B, IV-D, Tab III/IV) | 15.0% | 14.953 | \poolGcmed | MATCH |
| Curtailed hours at 1% (abstract, IV-B, Tab III/IV) | 3350 | 3350.25 | \poolGchrs | MATCH |
| 100 h gain (abstract, IV-B, Tab III/IV) | 6.8% | 6.820 | \poolHOZZmed | MATCH |
| Best 15-min forecast (abstract, Tab III note, IV-D) | 6.7% | P15 median 6.691 | robustness_fc.csv | MATCH |
| Best day-ahead forecast (abstract, Tab III note, V) | 3.2% | NLDA median 3.233 | robustness_fc.csv | MATCH |
| Firm HC base to N-1 (abstract, IV-D) | 6.1 to 2.9 MW | 6.134 / 2.877 | sens macros | MATCH |
| N-1 1% gain (abstract, IV-D, Tab III/IV) | 50.3% | 50.313 | \sensNmogO | MATCH |
| 0.25% / 1% / 5% hours (Sec. II-B: 35 136 intervals, 8784 h) | 35136 / 8784 | 35136 columns; 366 x 24 | - | MATCH |
| Pooled firm HC (IV-A, IV-B) | 6.1 MW | 6.134 | \firmPoolMed | MATCH |
| Firm HC 95% CI | 5.85–6.61 | 5.854–6.609 | stats.json | MATCH |
| Rural firm range | 1.5–42.8 MW | 1.484–42.780 | npz | MATCH |
| Median headroom over firm | 16.0% | 15.987 | \spreadMed | MATCH |
| Snapshot overstatement >1% | 42% of buses | 41.65% (t* = peaks.json t_peak_net) | npz + peaks.json | MATCH |
| Snapshot max overstatement | 12.5% | 12.532 | npz | MATCH |
| Gain at 0.25% | 11.5% | 11.476 | \poolGamed | MATCH |
| IQR at 1% | 10.9–19.9% | 10.896–19.914 | npz | MATCH |
| Gain at 5% | 21.7% | 21.701 | \poolGdmed | MATCH |
| Feeder groups | 17 (3 to 5 per grid) | 5+3+5+4 | stats.json | MATCH |
| CI at 0.25% | 9.6–12.7 | 9.637–12.703 | stats.json | MATCH |
| CI at 1% | 12.9–16.1 | 12.908–16.143 | stats.json | MATCH |
| CI at 100 h | 5.4–7.8 | 5.444–7.804 | stats.json | MATCH |
| Leave-one-grid-out 1% | 14.7–15.3 | 14.669–15.317 | stats.json | MATCH |
| Evening share of curtailed bus-intervals | 38% (25% of intervals) | 38.01% (0.250) | npz | MATCH |
| Nov–Feb share | 43% (33%) | 43.47% (33.06%) | npz | MATCH |
| Added MW at 1% | 0.82 | 0.8248 | npz | MATCH |
| Added MW CI | 0.69–0.96 | 0.686–0.963 | stats.json | MATCH |
| Share of year | 38% | 3350.25/8784 = 38.1% | - | MATCH |
| Events at 1% | 1263 | 1263 | npz | MATCH |
| Median event length | 45 min | 45 | npz | MATCH |
| Depth at 1% | 13% | 13.01 | npz | MATCH |
| Depth at 100 h | 6% | 6.38 | npz | MATCH |
| Commitment of about 85 h buys | 5.9–6.8% | 50 h: 5.94; 100 h: 6.82; **85 h: 6.60** | \poolHFZmed, \poolHOZZmed | MATCH to macros, see Issue 2 |
| Gu budget (1% of intervals) | 6.7% | order statistic 351/352: 6.660/6.667 | \gonepctint | MATCH |
| Copperplate setup | 185–258% | 184.67–258.03 (per grid) | audit/copperplate_reconcile.json | MATCH to JSON (copperplate not re-simulated) |
| Undervoltage / line buses | 79 / 365 | 79 / 365 (vhigh 0) | BIND at argmin | MATCH |
| Gain, undervoltage / line | 15.7 / 14.7 | 15.652 / 14.691 | npz | MATCH |
| Fixed taps, voltage / line | 28.5 / 15.3 | 28.476 / 15.311 | npz | MATCH |
| Spearman with Zth | 0.01 | 0.0110 | npz + meta zth_pu | MATCH |
| Spearman with firm HC | −0.18 | −0.1767 | npz | MATCH |
| Signs differ between grids | claim | Zth: urban −0.29, others + ; firm: urban +0.41, others − | npz | MATCH |
| r with spread | 0.98 | Pearson 0.982 | npz | MATCH |
| Fixed taps 1% gain (IV-D) | 16.7 | 16.695 | npz | MATCH |
| Scenario 2 buses | 460 | 96+119+138+107 = 460 | npz | MATCH |
| Scenario 2 1% gain | 18.8 | 18.834 | npz | MATCH |
| pf=1 subset size | 100 | 4 x 25 | npz | MATCH |
| pf=1 firm rise | 5.7% | 5.745 (bus-matched) | npz | MATCH |
| pf gain change | 15.1 to 14.6 | 15.091 (base, same 100 buses) to 14.567 | npz | MATCH |
| In-service line outages | 444 | 93+113+136+102 | n1_contingencies.json / N1_run logs | MATCH |
| Unrestorable spur outages | 20 | 6+5+1+8 entries with non-empty `lost` | *_oltc_N1_run.json | MATCH |
| N-1 undervoltage buses | 241 vs 79 | 241 / 79 | npz | MATCH |
| N-1 wide-band firm HC | 3.6 MW | 3.597 | npz | MATCH |
| N-1 spread | 64% vs 16.0% | 64.15 / 15.99 | npz | MATCH |
| N-1 wide-band 1% gain | 38.8 | 38.815 | npz | MATCH |
| N-1 100 h gain | 29.4 | 29.438 | npz | MATCH |
| N-1 added MW vs base | 1.32 vs 0.82 | 1.3232 / 0.8248 | npz | MATCH |
| Forecasts: SA15 / NL15 | 5.7 / 5.7 | 5.686 / 5.686 (not identical rows; r = 0.999) | csv | MATCH |
| Day-ahead range | 2.3–3.2 | SNW 2.288, SNDT 2.314, PDA 2.670, NLDA 3.233 | csv | MATCH |
| "About half the 15-min gain" | - | 3.233/6.691 = 0.48 | csv | MATCH |
| MAPE | 0.6–1.7% | 0.625–1.716 | csv | MATCH |
| Violation share 15-min | 35–41% | 34.97, 40.09, 41.20 | mean viol_any | MATCH |
| Violation share day-ahead | 48–53% | 48.11, 49.11, 52.45, 53.23 | mean viol_any | MATCH |
| Violating intervals | median 1–4 | 1–4 (medians over violating rows) | csv | MATCH (only if conditioned on viol_any; over all rows it is 0–1) |
| Excess | 1.0–1.6% | 1.035–1.564 (violating rows) | csv | MATCH (same conditioning) |
| Shape, ±10% at 03:00 firm | +9.8 | 9.785 | shapes.csv | MATCH |
| Shape, ±20% at 03:00 firm | +13.7 | 13.747 | shapes.csv | MATCH |
| Shapes at 15:00/19:00 under 2% | <2 | 1.29, 0.53, 1.73 (11:00: 1.23) | shapes.csv | MATCH |
| Energy change range | −1.2 to −18.6 | −1.19 to −18.62 (DC shapes) | shapes.csv | MATCH |
| Validation cases / overvoltage | 295 / 55 | 295 / 55 | audit_summary.json | MATCH |
| Validation error max / median | 3.4 / 0.4 kW | 3.375 / 0.367 | audit_summary.json | MATCH |
| dV / dL | <1e-10 p.u. / <1e-7 pp | 5.0e-11 / 9.8e-8 | val_runtime_check.json, validation.json | MATCH |
| Bus-intervals in boundary check | 15.8 million | 15,776,064; 0 failures | audit_summary.json | MATCH |
| Monotonicity scans | 2103 bus-intervals, 4.1 M PF | 2103; 4,128,517 steps; 0 nonmonotone | finescan_*.json | MATCH |
| N-1 pandapower check | 60 of 60 | 60/60 feasible_at and infeasible_above | n1_ppcheck_*.json | MATCH |
| Runtime | ~1 s; 10–35 s/bus | 0.89–1.31 s; 10.4–31.2 s | benchmark.json, val_runtime_check.json | MATCH |
| Solver constants | 1e-9, 60 it, 1e-6, 5 kW, 1e-7, 0.05 p.u. | identical in core.py | code | MATCH |

### Table I (gridtable.json, peaks.json)
Vn, transformers (2x25, 2x40, 2x63, 2x40), peak load (7.5/13.4/21.0/14.9), peak net load (6.3/12.1/19.8/13.4) and DG (25.6/23.8/13.6/16.6): all MATCH. MV buses 94/114/138/103 match the npz rows. gridtable.json gives 97/117/144/107 before merging switch-joined buses, which the table note covers.

### Table II (per grid, recomputed)
| Grid | Firm | 0.25% | 0.5% | 1% | 100 h | Hours | Longest |
|---|---|---|---|---|---|---|---|
| Rural | 5.6198 | 9.848 | 11.426 | 13.258 | 5.365 | 3072.00 | 19.00 |
| Semi-urban | 6.1844 | 10.959 | 12.438 | 14.214 | 6.686 | 3518.00 | 20.75 |
| Urban | 5.7724 | 11.604 | 13.199 | 15.088 | 7.032 | 3472.38 | 18.50 |
| Commercial | 7.2247 | 13.059 | 14.741 | 16.845 | 8.530 | 2945.00 | 18.75 |

All match the paper. Rounding of the .x5 longest-event values is discussed under Issue 1.

### Table III (recomputed: 0.25% / 1% / 100 h / N-1 1%)
- Firm HC: 6.134 / 6.134 / 6.134 / 2.877
- Added capacity: 0.6200 / 0.8248 / 0.3372 / 1.3232
- Gain: 11.476 / 14.953 / 6.820 / 50.313
- Hours: 1212.75 / 3350.25 / 100 / 1535.5
- Events: 926 / 1263 / 224 / 1060
- Median event length: 30 / 45 / 15 / 30 min
- Longest event: 15.0 / 19.0 / **3.25** / 15.75 h
- Depth: 10.29 / 13.01 / 6.38 / 33.47 %
- Curtailed energy: 0.250 / 1.000 / 0.0134 / 1.000 %

Every cell matches. The exception is rounding: 3.25 is printed as 3.2 and 15.75 as 15.8 (Issue 1).

### Table IV (recomputed)
| Configuration | Firm | 0.25% | 1% | 100 h | Hours | V/L/T |
|---|---|---|---|---|---|---|
| Base | 6.134 | 11.476 | 14.953 | 6.820 | 3350.25 | 79/365/5 |
| Fixed | 5.894 | 13.170 | 16.695 | 8.201 | 3189.75 | 106/341/2 |
| L60 | 3.214 | 17.139 | 23.214 | 10.620 | 2363.25 | 29/415/5 |
| N-1 | 2.877 | 38.656 | 50.313 | 29.438 | 1535.5 | 241/203/5 |
| N-1 wide | 3.597 | 29.188 | 38.815 | 20.229 | 1780.75 | 116/328/5 |
| Scenario 2 | 5.916 | 14.073 | 18.834 | 8.107 | 2492.75 | 82/**1 vhigh**/372/5 |

All match. The one difference is in Scenario 2: the raw BIND labels one bus "vhigh" (rural S2 bus 24). sensitivity.py relabels it to "line" from audit/boundary_rural_oltc_S2.json, where the reason above the headroom is line in all 35136 intervals. That gives 82/373/5, as printed. The value is justified, but the relabel is not mentioned in the paper.

## 2. Mismatches and issues

No numerical mismatch exceeds the last displayed digit. The minor issues are:

1. **Rounding convention (Table III, Table II).** The 100 h longest event is 3.25 h exactly but printed "3.2". At the same time 20.75 prints as "20.8" (Table II, semi-urban), 18.75 as "18.8" (Table II, commercial) and 15.75 as "15.8" (Table III, N-1). All four follow Python round-half-even and are internally consistent. A reader who assumes half-up rounding would expect 3.3. To avoid the question, print 3.25, or apply half-up rounding throughout.
2. **"a commitment of that size buys only about 5.9–6.8%" (IV-B, second paragraph after Table III).** "That size" refers to about 85 h/yr. The lower bound 5.9% is the 50 h value (\poolHFZmed). The hour-budget gain at 85 h is 6.60%, and 6.8% is the 100 h value. The range is a 50–100 h bracket, but it reads as the 85 h value. Suggest "about 6.6%" or "6.6% (6.8% at 100 h)".
3. **"(Table II)" after the pooled medians 11.5/15.0/21.7% (IV-B, first sentence).** Table II gives only per-grid values. The pooled 11.5 and 15.0 appear in Tables III and IV, and the 5% gain (21.7%) appears in no table. The cross-reference is valid but points to a table that does not contain the quoted numbers.
4. **"20 spur outages cannot be restored, and their buses are dropped" (IV-D, N-1).** All 449 buses remain in the N-1 results: the npz has 449 rows, no NaN, and Table IV counts 241+203+5 = 449. The disconnected buses are dropped only from those contingencies. As worded, a reader could take it to mean the buses are removed from the study. Suggest "their disconnected buses are excluded for that outage".
5. **Forecast violation statistics (IV-D).** "in a median of 1–4 intervals and by about 1.0–1.6% of nameplate" holds only for medians over the test cases that had a violation. Over all rows the medians are 0–1 intervals and 0–0.13%. The sentence ("had at least one violation, in a median of ...") implies the conditioning, so it is acceptable. It could be made explicit with "among those".
6. **Scenario 2 V/L/T relabel.** See the Table IV note above. One bus's class comes from the boundary audit rather than from the search label. This is defensible but undocumented.

## 3. Cross-references, formatting, build

- Page count: 5 (pdfinfo). OK.
- main.log: no undefined references or citations and no multiply-defined labels. There is one harmless Underfull hbox in the bibliography (lines 174–176, the Sarker reference). No "??" in the PDF.
- Placeholders: only "[REPOSITORY URL]" (intended). The keyword grep matched nothing else; "replaces" in IV-D was a false positive.
- Cross-references resolve correctly:
  - Table I is grids, Table II per-grid results, Table III contract terms, Table IV sensitivity.
  - Fig. 1 is the duration figure and Fig. 2 the binding-constraint figure.
  - "Section IV-D" points to Robustness, which holds both the N-1 case and the forecasts.
  - Eq. (2) is \eqref{eq:head} and Eq. (3) is the energy-budget equation.
  - "(1)" is hard-coded in "reduces to the limits in (1)" (Sec. II-A). It is correct, since the first equation is the unlabelled (1), but it is fragile.
- Tables and figures are first referenced in order.
- Same-quantity consistency: firm 6.1 (text) vs 6.13 (tables), N-1 2.9 vs 2.88, 15.0, 3350, 6.8, 50.3, 6.7 and 3.2 are identical in the abstract, body, Table III note and Section V.
- Every "from A to B" claim matches the tables or the raw data: 6.1 to 2.9 MW; 79 to 241 buses; 16.0% to 64%; 0.82 to 1.32 MW; 15.1% to 14.6%; fewer hours (3350 to 1536) and deeper curtailment (13 to 33%).
- tab_shapes.tex is not included in the paper, which is fine. Its values, as quoted in the text, were checked against shapes.csv.
- Unused macro files remain \input but do not affect the output: numbers_batt.tex (the battery sentence is now removed), numbers_oos.tex, and others. The numbers.tex macros \corrFirm = −0.17 and \corrZth = 0.08 differ from the printed \spFirm/\spZth values, but they are unused.

## 4. References
- 32 distinct \cite keys and 32 refs.bib entries. Every cited key exists, and every bib entry is cited. main.bbl has 32 \bibitem.
- The restored keys fani2024 [11], zakaria2025 [12] and umoh2023 [22] are cited and render.

## 5. Possibly unsupported or soft claims
- Intro figures (4.4%, 6.7–12%, 22 BAs, 76–126 GW, about 85 h/yr) and "Gu et al. reported a 78% gain" are external and NOT-CHECKABLE from the results.
- "so the difference stems from how much spare capacity the network has, not from the load data": supported only by the per-grid copperplate JSON (185–258%). The copperplate runs were not re-simulated here.
- "Most of the value requires near-real-time headroom estimation": consistent with the forecast results (6.7 vs 3.2 vs 15.0).
- "Curtailment ... by up to 13% of nameplate" is the median over buses of the per-bus maximum depth, not a maximum across buses (the pooled maximum is about 28%). The Table III note says "Medians over buses of each statistic", but the body sentence reads like an upper bound.
