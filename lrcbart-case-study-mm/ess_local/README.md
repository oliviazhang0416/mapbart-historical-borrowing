# Case-study ESS curve

`ess_cal.R` is sourced by `lrcBART.R`. For each outcome and requested
`H_f`, it fits the UCMM Stage-1 censored AFT-BART model, constructs the
Monte Carlo `ESS_tau0` curve over `s_0^2 = 10^seq(-6, 1, length.out = 141)`,
and saves a reusable checkpoint under `res/ess/`.

The selected targets are `100` (100% of the 30-patient EloKRd trial;
requested ESS 30), `75` (requested ESS 22.5), `50` (requested ESS 15), and
90%, 50%, and 25% of the feasible ESS ceiling. `ess_plot.R` overlays the `H_f=10`
and `H_f=50` curves after the analysis finishes. The default `H_g=5` is not
included in checkpoint or plot-name suffixes.

The active compiled ESS implementation is
`/Users/oliviazhang/Desktop/lrcBART/cess.cpp`; it is loaded directly and is
not copied into this application folder.

ESS checkpoint files live in the subproject's `res/ess/` folder. ESS plot
scripts save figures to the subproject's `inserts/` folder by default; `--out`
can select a different figure path. The R scripts remain in `ess_local/`.
