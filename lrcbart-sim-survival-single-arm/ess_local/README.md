# Survival single-arm ESS workflow

<!-- LRC-BART MODIFICATION START -->

ESS curve construction runs inside each `lrcBART.R` replicate. It is not a
separate simulation pass.

For each dataset, `ess_cal.R` fits a censored Stage-1 AFT-BART to the 300 RWD
controls and predicts the counterfactual control log-time surface at all
treated RCT profiles. It uses the resulting Monte Carlo variance, the RWD
residual variance, and draws from the discrepancy-tree prior to construct
`ESS_tau0(s0^2)` over

```r
10^seq(-6, 1, length.out = 141)
```

The target table contains `100`, `f90`, and `s0min`. `100` is capped at the
replicate's attainable ESS ceiling when necessary; `f90` is 90% of that
ceiling; `s0min` is the fixed near-complete-pooling reference with
`s0_sq = 1e-6`. Each replicate's checkpoint is saved under `res/ess/`
and reused by the default and sensitivity fits.

The default LRC-BART configuration fits `100` and `f90` with fixed `w = 1`.
The `w0.9` sensitivity fits the same two targets with fixed `w = 0.9`, and the
`s0min` sensitivity uses fixed `w = 1`. With no RCT controls, `g(X)` is drawn
from its prior and the control prediction is always `f(X) + g(X)`. The sampler
uses the RWD residual variance for both control variance slots.

`run_all.R` creates one ESS figure after the complete analysis sweep, using
the last requested replicate for each sample-size/scenario combination. The
plot scripts only read checkpoints; they do not refit a model. To plot a
specific checkpoint manually:

```bash
cd /Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-survival-single-arm/ess_local
Rscript ess_plot.R --res /absolute/path/to/checkpoint.RData
```

No C++ source is copied into this folder. `lrcBART.R` loads the active
standalone implementations from `/Users/oliviazhang/Desktop/lrcBART/`.


ESS checkpoint files live in the subproject's `res/ess/` folder. ESS plot
scripts save figures to the subproject's `inserts/` folder by default; `--out`
can select a different figure path. The R scripts remain in `ess_local/`.

<!-- LRC-BART MODIFICATION END -->
