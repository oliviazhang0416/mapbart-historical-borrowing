# Fresh validation design

This folder implements the frozen validation queue in the parent `VALIDATION-DESIGN.md` and `WORKFLOW.md`. It prepares a resumable run but does not launch fresh validation fits. The primary candidate is `JOINT-LRC-01`; `JOINT-SOURCE-01` is its matched joint source-indicator comparator, and `SEPARATED-SOURCE-01` is the variance-matched separated baseline.

The queue contains 200 datasets, three chains per method, and three methods: 1,800 posterior jobs. Every dataset has 300 randomized trial patients, 300 external controls, ten covariates, trial residual SD 1.5, and external residual SD 1.8. The joint methods use the private `joint-model/private-library` package and the signatures `joint_lrc_bart(y, x, source, treatment, ...)` and `joint_source_bart(y, x, source, treatment, ...)`. The treatment coefficient has prior \(\tau\sim N(0,100)\). `JOINT-LRC-01` uses 50 f trees, 5 g trees, the original \((\alpha_g,\beta_g)=(0.5,3)\) prior, five g sweeps, and the original training-only calibration object. `JOINT-SOURCE-01` uses 50 f trees, appends the binary source indicator, and sets \(H_g=0\). Both retain source-specific residual variances.

`SEPARATED-SOURCE-01` fits a separate 50-tree treatment BART and a 50-tree control BART with the binary source indicator and original source codes. Its two control variances remain source-specific. The joint models also share one residual variance across trial arms, whereas the separated treatment fit estimates its own variance. The comparison therefore evaluates the joint formulation as a whole; it does not isolate mean sharing from trial-arm variance pooling.

The constant-effect conditions use common seeds 92026001--92026040: `compatible` has no external shift, `one` has the known upward shift 2 outside ONE50, and `two` has the known upward shift 2 outside XOR50. The generator is called with the same seed for the three conditions, and assertions require identical trial covariates, assignments, and outcomes. The heterogeneous condition uses XOR50 with fresh seeds 92126001--92126040 and adds

\[
\tau_i=1+0.5\{X_{5i}-\overline X_{5,\mathrm{trial}}\}
\]

to treated outcomes, so its trial-profile ATE is exactly 1. The null condition uses XOR50 and seeds 92226001--92226040, subtracting 1 from treated outcomes so its true ATE is 0. The null gate is two-sided 95% exclusion of zero, with binomial uncertainty.

For each dataset, the calibration routine reproduces the existing setup-generation procedure: a 1,000-burn/1,000-draw 50-tree BART on RCT controls, a matching stage fit on external controls, and `lrc_ess_calibrate()` with `N_target=100`, \(H_g=5\), \((0.5,3)\), and 1,000/1,000 calibration settings. This is a practical prior-scale procedure; the old ESS=100 interpretation is not claimed for the joint model.

Each raw result stores posterior traces, ATE and control-surface draws, discrepancy summaries, source counts for both the separated controls and the full 300+300 joint data, convergence inputs, DGP/data hash, code and private-library hashes, seed lineage, and prior/source-variance configuration. Forest serializations are stripped after prediction. Writes are atomic and cache reuse requires an exact method, condition, replicate, chain, seed, burn-in, and draw-count match. The frozen full ledger is [jobs.csv](jobs.csv); filters are applied only to a copy written as `jobs_selected.csv`.

The preparation code is [validation_core.R](validation_core.R), the runner is [run_validation.R](run_validation.R), and the canonical summarizer is [summarize_validation.R](summarize_validation.R). The summarizer writes `canonical.csv` with the exact fields required by `evaluate.R`: `method`, `scenario`, `rep`, `ate_mean`, `truth`, `lower`, `upper`, `max_rhat`, and `min_bulk_ess`. It keeps g-contrast convergence in a separate `convergence.csv` record.

The prior pilot used 36 joint fits in about six minutes with three workers. The full batch has 1,800 jobs, with an extra treatment fit within each separated-baseline job and 200 dataset calibrations. The reviewed driver uses six fit workers; calibration is serial. Allow roughly three to four hours as a planning estimate, with runtime depending on tree sizes and filesystem overhead. This is not a completed runtime benchmark.
