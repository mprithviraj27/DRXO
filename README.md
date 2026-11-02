# DRXO — Dynamic Recrystallization Optimizer

Code and raw results for the article *"A dynamic-recrystallization model for adaptive differential evolution: thermal mobility, grain consumption and stored-energy migration"* (M. Prithvi Raj).

DRXO is a differential evolution algorithm derived from a coupled model of dynamic recrystallization: a thermal state set by the success rate controls the scale factor (Arrhenius mobility) and the nucleation parameters (Zener–Hollomon factor); accepted trials replace the oldest grains; a stored-energy driving force scales the attraction toward the growth front; strain triggers nucleation.

## Frozen algorithm

`drxo/drxo.py` is the exact file that was frozen before the held-out test.

```
SHA-256  bb924c68b0c4db67ee95a9984cb232729759547417719c4d5d3649f5eaa36bb1  drxo/drxo.py
```

Check with `sha256sum drxo/drxo.py`. The file needs only NumPy.

```python
from drxo import drxo_u_final as drxo
x_best, f_best, info = drxo(f, lb, ub, dim, max_fe=10000 * dim, seed=1)
```

`drxo_u_skeleton` in the same file is the DRX-off skeleton used as the ablation reference.

`drxo/drx_engine.py` is the development engine with all switches (ablation, factorial and post-hoc variants); `drxo/variants.py` names every configuration used in the paper (`C1sm` = DRXO, `skel2` = DRX-off skeleton, `srskel` = SR-skel, `abl_*`, `sens_*`, `F_*` factorial variants, `post_*` post-hoc diagnostics). With the `C1sm` configuration the engine reproduces the frozen file bit for bit.

## Benchmarks

`benchmarks/src` holds the organisers' C/C++ codes. Build the shared libraries, for example:

```
g++ -O2 -fPIC -shared benchmarks/src/cec14.cpp -o benchmarks/libcec14.so
g++ -O2 -fPIC -shared benchmarks/src/cec17.cpp -o benchmarks/libcec17.so
g++ -O2 -fPIC -shared benchmarks/src/cec22.cpp -o benchmarks/libcec22.so
g++ -O2 -fPIC -shared benchmarks/src/cec2020rwc.cpp -o benchmarks/cec2020rwc.so
```

The `input_data` folders of the official CEC 2014, 2017 and 2022 distributions must be in the working directory. The CEC 2014 source here contains one correction: input data are read with `%lf` instead of `%Lf`, which corrupts the data on Linux. Set `CEC_LIB_DIR` if the libraries are elsewhere. For CEC 2020, `rwc.py` passes a copy of each solution, because several problems round the input vector in place.

## Seeds and protocol

| Suite | Runs | Budget | Seed of run r on function/problem k |
|---|---|---|---|
| CEC 2014 (held-out) | 25 | 10,000·D | 14000 + 100k + r |
| CEC 2020 RC15–RC33 (held-out) | 25 | official MaxFEs | 20000 + 100k + r |
| CEC 2017 / CEC 2022 (design) | 25 / 30 (10 in screening) | 10,000·D / 200k, 1M | see `scripts/design_runs.py` |

Errors ≤ 1e-8 are set to zero. `docs/SCREEN_LOG.md` is the screening log written during development.

## Results

* `results/cec14/` — final errors on the held-out suite: `py_d{D}.json` rows `[function, algorithm, run, error]` for DRXO-U (DRXO), DRX-off and SR-skel; `lsrtde_d{D}.json`, `rde_d{D}.json`, `jso_d{D}.json` rows `[function, run, error]`; `traces_d50.json` convergence traces. L-SHADE values are the official CEC 2014 results and are not redistributed.
* `results/cec2020/res_rw20.json` — rows `[problem, algorithm, run, best feasible f, violation, evaluations]`. Official competition results of the CEC 2020 entrants are available from the organisers.
* `results/design/` — CEC 2017 and CEC 2022 runs of every screened, ablation, sensitivity and factorial variant, rows `[function, variant, run, error, ...]`.
* `results/diagnostics/` — instrumented runs (nucleation counts, acceptance, diversity), the F29 failure trace, the complexity measurement and the post-hoc late-switch runs (`posthoc_c14_d*.json`). The post-hoc rank-energy runs on CEC 2020 are in `results/cec2020/posthoc_rw.json`. The factorial variants (`F_sr_*`, `F_dx_*`, `srskel`) are in `results/design/res_cec17_d30.json` and `res_cec17_d50.json`.

## Scripts

`scripts/` contains the run scripts and `analysis/` the scripts that build every table of the paper (including the factorial decomposition at D = 30 and 50, the instrumented runs, the late-switch and rank-energy diagnostics), as used in the original environment; their absolute paths (`/home/claude/...`) must be adapted. `examples/run_example.py` is a self-contained starting point. `rivals/run_rivals14.sh` documents how L-SRTDE, RDE and jSO (authors' codes, obtained from the authors' repositories and not redistributed here) were linked to the CEC 2014 functions; only the benchmark call, seeds and an output line were changed.

## Citation

Please cite the article (details will be added on publication).

## License

MIT for the DRXO code; benchmark codes remain under the organisers' terms.
