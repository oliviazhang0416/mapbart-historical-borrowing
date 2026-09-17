# Gaussian two-arm ESS workflow

<!-- LRC-BART MODIFICATION START -->

ESS calibration is part of each method's replicate loop; it is not a separate
project-wide run.

- `lrcBART.R` sources `ess_cal.R`. For each dataset, it fits the Stage-1 RWD
  BART, estimates the RCT-control residual variance, calibrates all standard
  LRC-BART targets, and saves a checkpoint in `res/ess/`.
- The standard targets `100`, `75`, `50`, `f90`, `f50`, and `f25` reuse that
  one checkpoint. The `w0` and `w1` calls reuse its target-100 calibration;
  `Hg10` receives a separate `H_g=10` checkpoint.
- `MAP.R` sources `ess_cal_map.R` and checkpoints the inherited RBesT/gMAP
  calibration for targets 100, 75, and 50 for each dataset.
- No C++ source is copied here. `lrcBART.R` loads the active standalone files
  `/Users/oliviazhang/Desktop/lrcBART/clrcbart.cpp` and
  `/Users/oliviazhang/Desktop/lrcBART/cess.cpp` before sourcing `ess_cal.R`.

The optional plot scripts read completed checkpoints and never refit a model:

```bash
cd /Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-gaussian/ess_local
Rscript ess_plot.R --scenario sc1
Rscript ess_plot_simple.R --scenario sc1
```

Pass `--res /absolute/path/to/checkpoint.RData` to select an exact replicate.


ESS checkpoint files live in the subproject's `res/ess/` folder. ESS plot
scripts save figures to the subproject's `inserts/` folder by default; `--out`
can select a different figure path. The R scripts remain in `ess_local/`.

<!-- LRC-BART MODIFICATION END -->
