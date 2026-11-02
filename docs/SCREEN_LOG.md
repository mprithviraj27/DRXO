# Unified DRX engine — single-objective screening log

Design (screening) sets: CEC 2017 D10 (25 runs), CEC 2022 D10 (30 runs), CEC 2017 D30 (10 runs).
Statistic: Wilcoxon rank-sum p < 0.05 with median and mean agreeing; better / tie / worse.
Held-out CEC 2014 untouched. All variants defined in `variants.py`; engine `drx_so.py`
(all new switches off = frozen v0.7, checked bit-for-bit by `test_equiv.py`).

## Baselines (D30, n = 10)
- v0.7 vs skeleton 15/14/0; v0.7 vs L-SRTDE 1/8/20; v0.7 vs L-SHADE 9/14/6.

## Round 1 — one module on top of v0.7 (D30 vs v0.7)
| module | D10 vs v0.7 | D30 vs v0.7 | decision |
|---|---|---|---|
| adaptive temperature → mobility (at_mob) | 8/17/4 | 5/19/5 (gains F5, F7, F8, F10; losses F13–F15) | keep |
| adaptive temperature → growth-front fraction | 2/26/1 | 3/21/5 | drop |
| adaptive temperature → nucleation only | 4/25/0 | 1/27/1 | (included via temperature) |
| Boltzmann partner selection | 2/24/3 | 2/26/1 | drop |
| texture crossover (p = 0.5) | 3/25/1 | 4/20/5 | retest later |
| SRX polishing (CMA-ES, last 5 %) | 0/27/2 | 0/15/14 | drop (harmful) |
| CDRX dual mode (θc = 0.1) | 0/29/0 | 0/29/0 | never triggered; drop |
| reheating | 0/27/2 | 0/28/1 | drop |

## Diagnosis of the hybrid-function gap (F11–F20) to L-SRTDE
Crossover-memory variants and bound handling: no consistent effect.
Replacement rule is the key: a successful trial replacing the *oldest* grain (grain consumption)
instead of its parent: F10 134 → 12, F13 61 → 18, F15 15 → 5, F20 25 → 7 (6 runs).

## Round 2 — grain consumption (A1 = adaptive mobility + age consumption)
D30: A1 vs v0.7 8/20/1; A1 + texture 11/17/1. Still 2/12/15 vs L-SRTDE.
Diagnosis: the prior-boundary archive duplicates the surviving old grains; removing it gives
F5 7 → 1.4, F8 8 → 1.8, F14 24 → 17 (5 runs).

## Round 3 — no archive (B0 = A1 + texture, r_arc = 0)
D30: B0 vs v0.7 19/10/0, vs L-SHADE 19/10/0, vs RDE 13/11/5, vs L-SRTDE 4/16/9.
Ablation at D30: nucleation 5/21/3; texture 0/29/0; ε_c0 = 6 vs 3: 4/25/0.
D10: B0 vs v0.7 8/13/8, vs L-SRTDE 2/22/5. CEC 2022 D10: mostly ties (saturated).
Decision: drop texture (no effect), raise critical strain to ε_c0 = 6 → candidate C1.

## Round 4 — C1 (texture removed, eps_c0 = 6) and realised-bulge (crossover) memory
C1sm vs C1 over all screening sets: 6/62/2 -> C1sm adopted. C1sm at D50 (10 runs): vs skel2 23/4/2, jSO 22/5/2,
L-SHADE 23/4/2, RDE 16/6/7, L-SRTDE 5/13/11; losses on F1, F3 (final precision).

## Round 5 — consumption only of higher-energy grains (C2 = C1sm + accept='age_worse')
C2 vs C1sm: D10 3/23/3, CEC22 D10 3/8/1, D30 0/16/13, D50 4/10/15 -> rejected.

## Freeze
C1sm frozen as drxo_u_final.py (SHA-256 bb924c68...aa36bb1) before any run on CEC 2014 or CEC 2020 RC15-RC33.
Held-out dimensions fixed beforehand: D = 30, 50 (D = 10 added later, code unchanged).
