# LRC-BART joint model

> ### Reviewers start here — [HANDOFF-REVISION-LRC-BART-20260925.md](HANDOFF-REVISION-LRC-BART-20260925.md)
>
> This branch migrates PR#2's joint-treatment formulation into the code. Three findings, in brief:
>
> 1. **Joint LRC-BART performs about the same as the non-joint version** on `revision/lrc-bart-20260917`, which does not model the treatment effect jointly. Parity with the simpler model argues against the joint layer, not for it.
> 2. **The claimed advantage over the BART counterparts is unclear.** In PR#2's own validation only 1 of 5 scenarios favours joint LRC-BART; it is *worse* than the separated baseline in two, and the remaining gaps sit inside Monte Carlo error.
> 3. **PR#2's evidence covers two-arm Gaussian only** — no survival, no single-arm, no application — yet it makes joint modeling primary for the whole paper.
>
> The note carries the numbers, the MCSE caveats and what would settle each point.

<!-- JOINT LRC-BART MODIFICATION START -->
PR#2 joint-model migration from `lrcbart-historical-borrowing`, following the
folder's migration plans. The shared C++ model and Gaussian/survival two-arm
fitting workflows are implemented. Gaussian single-arm generation, fitting,
reporting and plotting are also migrated using PR#2’s separate-treatment
formulation, as are the survival single-arm workflows. The case-study scripts are also migrated using the separate-treatment
formulation; manuscript migration remains deferred.

Both migrated generators retain 100 replicates per configuration. Their 15
configurations use Sc1, Sc1a, Sc1b, Sc1c, Sc1c-i, Sc1c-ii, Sc2, and Sc3.
Joint LRC fits use one chain, 1,000 warm-up iterations, 1,000 retained
draws, no thinning, and one discrepancy sweep per iteration. Preliminary ESS calibration stays
at 1,000/1,000. These sampling defaults match `lrcbart-historical-borrowing`; migration checks use temporary
outputs and explicitly reduced settings.

The Gaussian and survival two-arm `run_all.R` scripts retain their full method
lists and run all configured scenarios by default. Pass
`--scenarios sc1c,c-i,c-ii` to select only those three configurations. Balance, result and ESS plots are enabled. Add `--plots-only` to regenerate
figures from saved inputs without generating data or fitting models. Joint results use `LRC-BART-joint_` filenames and save compact
reporting rows after each replicate, with run settings and partial/completed
status (`planned_replicates` and `complete` attributes). Posterior draws,
chain summaries, and diagnostic fields are not saved. All configured chains
still contribute to the calculations in memory. There is no chain-level resume;
a new invocation refits its selected work. Compact ESS files retain the curves
and calibration quantities needed for fitting and plots. No old data, results, inserts,
or patient datasets are copied into this folder.

Gaussian single-arm retains 200 treated trial participants, 300 historical
controls, and 100 replicates. Its 12 configurations cover Sc1, Sc1a, Sc1b,
Sc1c, Sc1c-i, Sc1c-ii and Sc2. Sc1c-i is heterogeneous and Sc1c-ii is null.
Treatment is fitted independently; the historical-plus-discrepancy model
predicts hypothetical trial controls. The discrepancy remains prior-driven,
with equal trial-control and historical-control residual variances for that
prediction. The independent treatment model has its own residual variance.
Five LRC fits are run: fixed w=1 with ESS 100/f90, fixed w=0.9 with ESS
100/f90, and fixed w=1 with s0min=1e-6. Each fixes tau0_sq at
min(s0_sq, tau1_sq/2). Fits and preliminary calibration use 1,000 burn-in
and 1,000 retained draws. LMv2 (historical weight 1), BARTv2 and hierLM
(prior scale 0.05) retain their original models and sampling settings.
Single-arm result names retain `LRC-BART_`, with `joint_model=FALSE` in
reporting rows. Treatment predictions and control forests stay in memory;
no treatment-fit cache, chain files or diagnostic summaries are saved.
The single-arm runner supports the same `--scenarios` and `--plots-only`
arguments as the two-arm runners.

Survival single-arm uses the same separate-treatment formulation and five LRC
configurations. Both treated sample sizes (30 and 200) run by default with 300
historical controls and 100 replicates in each of the 12 scenario configurations.
The original dropout and administrative-censoring distributions are retained;
scenario shifts are applied to latent times before recomputing observed times
and event indicators. Heterogeneous truth is evaluated on each trial's profiles.
Median-survival-ratio and five-year RMST-ratio summaries retain the original
separate-arm prediction and numerical integration conventions. AFTv2 (historical
weight 1), BARTv2 and hierAFT (prior scale 0.05) retain their statistical models.
LRC fits and preliminary calibration use 1,000 burn-in and 1,000 retained draws.
Results remain in `res/n30` and `res/n200` with sample-size-tagged filenames;
ESS files remain in `res/ess`. The runner supports `--scenarios` and
`--plots-only`, and enables result, balance and ESS plots. Only compact reporting
and ESS files are saved; no treatment-fit cache, posterior chain files or
chain summaries are written.
<!-- JOINT LRC-BART MODIFICATION END -->

## Layout

| Path | Contents |
|---|---|
| `lrcBART/` | LRC-BART C++ engine (`clrcbart.cpp`, `cess.cpp`, `include.lrc/`) |
| `wBART/`, `aBART/` | Gaussian and AFT comparator engines |
| `psrwe/` | Vendored psrwe (Wang et al., GPL >= 3) for the PSCL comparator |
| `bartModelMatrix.R` | Shared design-matrix helper |
| `lrcbart-sim-gaussian{,-single-arm}/` | Gaussian simulations |
| `lrcbart-sim-survival{,-single-arm}/` | Survival simulations |
| `lrcbart-case-study-mm/` | Myeloma application (EloKRd vs UCMM) |
| `manuscript/` | Deferred; not copied during this step |

## Running

Each project is driven by its own `run_all.R`:

```bash
cd lrcbart-sim-gaussian
Rscript run_all.R
```

Scripts locate the repository root by walking up from the script's own
location (and, failing that, the working directory) until they find the
`.lrcbart-root` marker. You can therefore invoke them from any directory:

```bash
Rscript /path/to/lrcbart-joint-model/lrcbart-sim-gaussian/run_all.R
```

No paths need editing.

Requires R with `Rcpp` and `RcppEigen`; joint result provenance also uses `digest`. The C++ engines compile on first use
via `Rcpp::sourceCpp`. `PSCL.R` additionally needs `miceadds` and `dplyr`, and
`MAP.R` needs `RBesT`.

### Comparator provenance

`PSCL.R` sources `psrwe/R/` directly from the vendored copy in this
repository, so it runs pinned local source rather than an installed build.
`PSRWE_DIR` overrides the location if you need a different checkout.

## Not in this repository

`res/`, `inserts/`, and the simulation `data/` directories hold regenerable
output and are excluded from this copy. Rebuild simulated data with each
project's `data_gen_p10.R`, then run `run_all.R`.

`MAP.R` and `ess_local/ess_cal_map.R` use the **installed** `RBesT` package
(`library(RBesT)`), not a local copy. Results depend on the installed version,
so record it when reporting: `packageVersion("RBesT")`.

Every other comparator is included in this repository.

Patient-level inputs are excluded. The case-study runner reads the existing
merged cohort from the original historical-borrowing folder without copying it.
Set `MERGED_FILE` to override that input path. Case-study fits retain four
chains with 2,000 burn-in and 2,000 retained draws; preliminary ESS calibration
retains 4,000 draws after 2,000 burn-in. Only compact reporting and ESS outputs
are saved. `--plots-only` regenerates ESS plots and prints saved summaries.
Case-study execution and validation were skipped at the user’s request.
