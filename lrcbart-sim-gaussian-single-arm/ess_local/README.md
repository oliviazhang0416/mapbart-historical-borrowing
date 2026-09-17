# Gaussian single-arm ESS workflow

<!-- LRC-BART MODIFICATION START -->

ESS calibration runs inside each `lrcBART.R` replicate. It is not a separate
simulation step.

For each dataset, `ess_cal.R` fits the Stage-1 RWD control BART, predicts the
counterfactual control mean at the RCT treatment profiles, estimates the
single-arm ESS ceiling, and constructs the `100`, `f90`, and `s0min` settings.
`s0min` uses the smallest grid value, `s0_sq = 1e-6`. The checkpoint is written
to `res/ess/` and reused by every LRC-BART configuration for that
replicate.

The `default` configuration fits the `100` and `f90` targets with fixed
`w = 1`. The same calibration checkpoint supports the `w0.9` configuration at
both targets and the near-complete-pooling `s0min` sensitivity fit. Because this
is a single-arm design, there are no RCT controls:
the discrepancy forest `g` is prior-only and `tau0` is not updated from a
missing RCT-control likelihood.

The optional plot scripts read a completed checkpoint and never refit a model:

```bash
cd /Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-gaussian-single-arm/ess_local
Rscript ess_plot.R --scenario sc1
Rscript ess_plot_simple.R --scenario sc1
```

Pass `--res /absolute/path/to/checkpoint.RData` to plot an exact replicate.
No C++ source is copied into this folder; calibration uses the active standalone
`/Users/oliviazhang/Desktop/lrcBART/cess.cpp` implementation loaded by
`lrcBART.R`.


ESS checkpoint files live in the subproject's `res/ess/` folder. ESS plot
scripts save figures to the subproject's `inserts/` folder by default; `--out`
can select a different figure path. The R scripts remain in `ess_local/`.

<!-- LRC-BART MODIFICATION END -->
