# Joint-validation reproducibility source pack

This is a portable source-and-summary companion for the 200-dataset joint Gaussian validation. It supports source inspection, package rebuilding and no-fit recomputation of the published method-condition summaries. It is not a new fit, a replacement for the original numerical archive or a publication claim.

The bundle preserves the source paths used by the validation: `01-code/R/dgp.R`, `08-interaction-pilot/20260918/code/scenarios.R`, `09-method-exploration/state.json`, the validation scripts and ledgers under `09-method-exploration/validation/`, and the source-only `09-method-exploration/joint-model/lrcbart-joint/` package. The selected result set contains `canonical.csv` (about 5 MB), `summary.csv`, `convergence.csv`, `convergence_summary.csv`, `jobs_selected.csv`, the run report/configuration, and the three reviewed paired-comparison summary directories. The current table values are retained in `summary.csv` and are recomputable from `canonical.csv`; no alternate fitted table is generated. The 1,800 raw chain files and data/calibration caches are intentionally excluded; the original raw archive remains in the parent revision workspace under `09-method-exploration/validation/results/`.

The frozen preselection timestamp is `2026-09-19T22:32:14Z`, taken from `09-method-exploration/state.json` (`candidate_frozen_at_utc`). The full numerical design is 1,800 jobs (200 datasets, three methods, three chains), with 30,000 burn-in iterations and 12,000 retained draws per job. Full execution is intentionally never automatic and is expected to require substantial compute. The supplied generators can recreate the simulated data and calibration caches in a fresh run.

## Dependencies and isolated package build

The package `DESCRIPTION` declares `Rcpp` (and `stats`) for the package build, with `testthat`, `BART` and `survival` in `Suggests`/runtime support. The validation scripts additionally use the installed R packages `digest`, `jsonlite` and `posterior`; `parallel`, `stats` and `tools` are base/recommended R packages. Use a private library inside the bundle:

```sh
cd 05-writing/journal-revision/reproducibility/joint-validation
mkdir -p 09-method-exploration/joint-model/private-library
R CMD INSTALL --library=09-method-exploration/joint-model/private-library \
  09-method-exploration/joint-model/lrcbart-joint
```

The original canonical rows retain code/library hashes that encode the historical absolute paths and compiled library provenance. The source-only pack is therefore portable for inspection and rebuilding, but it does not promise bitwise identity across platforms or toolchain versions.

## No-fit verification

The standard-library Python check reads only `09-method-exploration/validation/results/canonical.csv` and `summary.csv`. It recomputes per-group RMSE, the delta-method RMSE MCSE, 95% coverage counts and ATE flags (`R-hat > 1.05` or bulk ESS `< 400`):

```sh
python3 verify_joint_validation.py
```

No R fit, MCMC job, data cache or raw chain is used by this check. The exact successful output is retained in `VERIFY-OUTPUT.txt`, and `SHA256SUMS.json` lists SHA-256 hashes for every bundled file except the manifest itself.

## Optional rerun commands

The source scripts are retained for a deliberate, separately reviewed rerun. Work in a copy of this bundle. After the private package build and dependency setup, move the supplied summary directory aside to preserve it, then run from the copied bundle root:

```sh
mv 09-method-exploration/validation/results 09-method-exploration/validation/results-supplied
export VALIDATION_PRESELECTION_UTC=2026-09-19T22:32:14Z
export VALIDATION_PRIMARY_CANDIDATE_ID=JOINT-LRC-01
export VALIDATION_PRIMARY_BASELINE_ID=SEPARATED-SOURCE-01
R_LIBS_USER="$PWD/09-method-exploration/joint-model/private-library" \
  Rscript 09-method-exploration/validation/run_validation.R
R_LIBS_USER="$PWD/09-method-exploration/joint-model/private-library" \
  Rscript 09-method-exploration/validation/summarize_validation.R
```

Those commands target the heavy 1,800-job, 30,000/12,000 configuration and are documented for reproducibility only; this bundle does not autorun them. The portable pack contains no data cache; the supplied generators and calibration code create fresh caches in the copied bundle. Do not mix freshly compiled fits with caches from the original binary, because the runner checks code and library hashes. `checks.R` is a small design/API check and does not replace the full validation.

## Saved-chain sensitivity

The `chain-sensitivity/` companion reuses the reviewed saved-trace summary `INPUT-constant_chain_rmse.csv`. It reports RMSE for individual chains 1--3 across all 120 constant-effect method--dataset groups and computes LRC percentage reductions against both comparators. It does not filter convergence flags, refit models or establish convergence. Re-run the derivation with:

```sh
python3 chain-sensitivity/derive_chain_sensitivity.py \
  --input chain-sensitivity/INPUT-constant_chain_rmse.csv \
  --out-dir chain-sensitivity
```
