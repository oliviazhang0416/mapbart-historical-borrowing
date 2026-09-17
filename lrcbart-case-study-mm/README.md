# LRC-BART myeloma case study

This application compares the 30-patient EloKRd trial treatment cohort with
253 UCMM external controls for progression-free survival (PFS) and overall
survival (OS). All adjusted methods use one harmonized seven-column covariate
matrix:

```text
age, male, race_Black, race_Other, hispanic, high_risk_cyto, asct
```

The application does not construct or run an original-coding analysis.

## Run the full application

```bash
cd /Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-case-study-mm
Rscript run_all.R
```

`run_all.R` cleans and merges the private data, runs PFS and OS, constructs
the ESS plots after all model fits finish, and writes the result table.

## Run settings

All configured methods run by default. `SCRIPTS` can select case-study
analysis methods while leaving data preparation, available ESS plotting,
and the result table enabled.

`OUTCOMES` and `LRC_CONFIGS` select both analyses and summary rows. For example,
`OUTCOMES=PFS LRC_CONFIGS=Hf50 Rscript run_all.R` runs PFS analyses with the
selected LRC configuration, plots its H_f=50 ESS curve, and summarizes that
configuration alongside available PFS comparator results.

Independent stages continue after a failure, but failed preparation and model
outputs are excluded from dependent work. The final status summary raises an
error if any requested work failed or was blocked.

## Production configurations

- KM: unadjusted five-year RMST ratio, once per outcome.
- AFT-CP: harmonized covariates with UCMM likelihood weight `w=1`.
- HierAFT: harmonized covariates with discrepancy prior values `0.05` and
  `0.5`; one chain, 1,600 total iterations, and 100 warm-up iterations.
- Standard BART: harmonized covariates, 50 trees per arm, 2,000 burn-in and
  2,000 retained draws.
- LRC-BART: harmonized covariates; `H_g=5`; four chains with 2,000 burn-in and
  2,000 retained draws. The four configurations are `default` (`H_f=10,
  w=1`), `Hf50`, `w0.9`, and `Hf50_w0.9`. Each runs targets `100`
  (100% of the 30-patient EloKRd trial, requested ESS 30), `75` (requested
  ESS 22.5), `50` (requested ESS 15), `f90`, `f50`, and `f25`.

The ESS Stage-1 fit uses 2,000 burn-in and 4,000 retained draws. The supplied
calibration range is `s_0^2 = 10^seq(-6, 1, length.out=141)`.

## Control estimates and RMST reporting

For PFS and OS, every model reports treatment, hypothetical-control, and
UCMM-control RMST up to five years. All three use the **EloKRd patients'
covariates**, including the same centering used for model fitting. Each
posterior draw averages individual patients' RMSTs over EloKRd; predictions
are not evaluated at a single average covariate profile.

The UCMM-control estimate is therefore the external-control model standardized
to EloKRd, rather than an average over the observed UCMM population:

- LRC-BART uses `f_test + y_center` with `sqrt(sigma2_sq)` for UCMM and
  `f_test + g_test + y_center` with `sqrt(sigma1_sq)` for hypothetical control.
- HierAFT uses the external `alpha2`/`beta2` layer for UCMM and the concurrent
  `alpha1`/`beta1` layer for hypothetical control, with the same residual variance.
- AFT-CP and Standard BART have no discrepancy layer, so their UCMM and
  hypothetical-control estimates are identical.

Model result files save `rmst_ucmm` posterior draws and `rmst_ucmm_est`
(`est`, `lo`, `hi`: posterior median and 95% equal-tailed credible interval),
with `rmst_ucmm_population = "EloKRd"`. Existing `rmst_ctrl` and `rmst_ctrl_est`
continue to represent hypothetical control. The treatment/control ratio and
difference continue to use that hypothetical control.

KM is unadjusted: its `rmst_ucmm_est` describes the observed UCMM cohort,
with `rmst_ucmm_population = "UCMM"` and a 95% confidence interval. Its horizon
is capped at the shorter arm's maximum follow-up, and its treatment contrast
uses observed UCMM control. KM has no hypothetical-control estimate.

`res/results_table.csv` and `res/results_table.RData` include numeric
`rmst_trt_*`, `rmst_hyp_ctrl_*`, and `rmst_ucmm_*` estimate/lower/upper columns,
plus the UCMM population, RMST horizon, interval type, and contrast-control
label. The printed summary also shows the three arm estimates. These additions
run automatically through `run_all.R` for all selected outcomes/configurations.
Older LRC-BART and HierAFT files without the new UCMM summary show missing UCMM
values until their analysis scripts are rerun.

## Private-data preparation

The source data and cleaning scripts are under `private_data/`:

```bash
Rscript private_data/data_cleaning_ucmm.R
Rscript private_data/data_cleaning_elokrd.R
Rscript private_data/data_merge.R
```

The last command writes the de-identified merged analysis files and the seven
harmonized model columns under `data_cleaned/`. Model results go to `res/`;
KM plots go to `res/`; ESS checkpoints go to `res/ess/`, and ESS plots
go to `inserts/`.

The active model sources are loaded directly from:

```text
/Users/oliviazhang/Desktop/lrcBART/clrcbart.cpp
/Users/oliviazhang/Desktop/lrcBART/cess.cpp
```
