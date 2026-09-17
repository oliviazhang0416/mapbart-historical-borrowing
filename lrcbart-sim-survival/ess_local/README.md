# Local ESS workflow (survival two-arm)

`ess_cal.R` is sourced by `lrcBART.R` inside its simulation-replicate loop.
For the current dataset it fits the censored-survival calibration models,
selects the `s_0^2` values for 100, 75, 50, f90, f50, and f25, and saves a
resume checkpoint under `res/ess/`.

The calculation uses the active standalone sources directly:

- `/Users/oliviazhang/Desktop/lrcBART/clrcbart.cpp`
- `/Users/oliviazhang/Desktop/lrcBART/cess.cpp`

No C++ source or include tree is copied into this subproject. Censored rows
have status `0`; their centered observed log time is passed as the lower bound
for latent-time augmentation.

`ess_plot.R` and `ess_plot_simple.R` only plot checkpoints already produced by
`lrcBART.R`; they do not run a separate simulation loop.

ESS checkpoint files live in the subproject's `res/ess/` folder. ESS plot
scripts save figures to the subproject's `inserts/` folder by default; `--out`
can select a different figure path. The R scripts remain in `ess_local/`.
