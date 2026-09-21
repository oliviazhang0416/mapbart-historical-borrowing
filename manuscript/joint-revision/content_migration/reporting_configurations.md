# Simulation and application reporting configurations

Status: settled with the user on 2026-09-16.

Subsequent drafting update: the approved detailed structure is now implemented in the three manuscript documents. Statements below about table layout being pending describe the earlier configuration-settlement checkpoint; see the collaborator handoff for the current draft and outstanding scientific issues.


Purpose: authoritative reporting set for the four simulation folders and the myeloma application. Report all completed configurations listed here, using the stated default setting in the main comparison and the remaining settings in sensitivity tables. This settles configuration inclusion and the primary references; exact table layouts, performance columns, scenario allocation between main text and appendix, and interpretation of findings remain to be settled.

Saving this document does not change analyses, rerun models, or resolve the implementation qualifications recorded in [implementation notes](implementation_notes.md).

## 1. Gaussian two-arm

Folder: `lrcbart-sim-gaussian`.

LRC-BART settings: H_f = 50, H_g = 5; update w under a Beta(1, 1) prior.

| Role | Requested ESS setting |
| --- | --- |
| Main comparison | 100 |
| Sensitivity analyses | 75, 50, f90, f50, f25 |

Comparators:

- LM: trial-only, complete pooling, and source-indicator adjustment (`LMv1`, `LMv2`, `LMv3`).
- BART: trial-only, complete pooling, and source-indicator adjustment (`BARTv1`, `BARTv2`, `BARTv3`).
- Hierarchical LM: prior setting 0.25.
- PSCL.
- MAP: requested target 100 in the main comparison; targets 75 and 50 in the target-sensitivity comparison.

There are six LRC-BART configurations and 17 method/configuration combinations per scenario setting in the saved reporting set. Existing results cover Sc1-Sc5 through 12 scenario settings; their exact placement in main and supplementary tables remains open.

Source references: [run_all.R](../../lrcbart-sim-gaussian/run_all.R), [lrcBART.R](../../lrcbart-sim-gaussian/lrcBART.R), and [ESS calibration](../../lrcbart-sim-gaussian/ess_local/ess_cal.R).

## 2. Survival two-arm

Folder: `lrcbart-sim-survival`.

LRC-BART settings: H_f = 50, H_g = 5; update w under a Beta(1, 1) prior.

| Role | Requested ESS setting |
| --- | --- |
| Main comparison | 100 |
| Sensitivity analyses | 75, 50, f90, f50, f25 |

Comparators:

- AFT: trial-only, complete pooling, and source-indicator adjustment (`AFTv1`, `AFTv2`, `AFTv3`).
- BART: trial-only, complete pooling, and source-indicator adjustment (`BARTv1`, `BARTv2`, `BARTv3`).
- Hierarchical AFT: prior setting 0.25.

There are six LRC-BART configurations and 13 method/configuration combinations per scenario setting. Existing results cover the same 12 Sc1-Sc5 settings as the Gaussian two-arm study. MAP and PSCL are not among this folder's saved methods.

Median-survival-ratio and RMST-ratio summaries are available. Their reporting priority remains to be settled, retaining the recorded RMST implementation qualification.

Source references: [run_all.R](../../lrcbart-sim-survival/run_all.R) and [lrcBART.R](../../lrcbart-sim-survival/lrcBART.R).

## 3. Gaussian single-arm

Folder: `lrcbart-sim-gaussian-single-arm`.

LRC-BART settings: H_f = 50 and H_g = 5. Report these five configurations:

| Role | Configuration | Fixed w | ESS/prior setting |
| --- | --- | --- | --- |
| Main comparison | default | 1 | Target 100 |
| Sensitivity | default | 1 | f90 |
| Sensitivity | w0.9 | 0.9 | Target 100 |
| Sensitivity | w0.9 | 0.9 | f90 |
| Sensitivity | s0min | 1 | s_0^2 = 10^-6 |

Comparators:

- LM with full RWD likelihood weight (`LMv2`, RWD weight 1).
- Standard BART (`BARTv2`).
- Hierarchical LM with prior setting 0.05.

Report Sc1-Sc2 with 200 treated trial patients and 300 historical controls. There are nine method/configuration combinations per scenario. Targets 75, 50, f50, and f25 are not part of the saved single-arm reporting set.

Source references: [run_all.R](../../lrcbart-sim-gaussian-single-arm/run_all.R), [lrcBART.R](../../lrcbart-sim-gaussian-single-arm/lrcBART.R), and [data generation](../../lrcbart-sim-gaussian-single-arm/data_gen_p10.R).

## 4. Survival single-arm

Folder: `lrcbart-sim-survival-single-arm`.

Report the same five LRC-BART configurations as Gaussian single-arm, with H_f = 50 and H_g = 5:

| Role | Configuration | Fixed w | ESS/prior setting |
| --- | --- | --- | --- |
| Main comparison | default | 1 | Target 100 |
| Sensitivity | default | 1 | f90 |
| Sensitivity | w0.9 | 0.9 | Target 100 |
| Sensitivity | w0.9 | 0.9 | f90 |
| Sensitivity | s0min | 1 | s_0^2 = 10^-6 |

Report them separately for n_1 = 30 and n_1 = 200, both with 300 historical controls and Sc1-Sc2. The requested target 100 is an absolute ESS count in both sample-size settings; it is not 100% of n_1.

Comparators:

- AFT with full RWD likelihood weight (`AFTv2`, RWD weight 1).
- Standard BART (`BARTv2`).
- Hierarchical AFT with prior setting 0.05.

There are nine method/configuration combinations per scenario and trial sample size. Targets 75, 50, f50, and f25 are not part of the saved single-arm reporting set.

Source references: [run_all.R](../../lrcbart-sim-survival-single-arm/run_all.R) and [lrcBART.R](../../lrcbart-sim-survival-single-arm/lrcBART.R).

## 5. Myeloma application

Folder: `lrcbart-case-study-mm`.

Use harmonized covariates, 30 EloKRd patients, and 200 verified triplet-treated UCMM controls.

Report the complete LRC-BART 2-by-2 configuration set:

| Configuration | H_f | H_g | Fixed w |
| --- | --- | --- | --- |
| default | 10 | 5 | 1 |
| Hf50 | 50 | 5 | 1 |
| w0.9 | 10 | 5 | 0.9 |
| Hf50_w0.9 | 50 | 5 | 0.9 |

Each configuration has six ESS target settings:

| Saved label | Requested ESS |
| --- | --- |
| 100 | 30, or 100% of the trial sample size |
| 75 | 22.5, or 75% of the trial sample size |
| 50 | 15, or 50% of the trial sample size |
| f90 | 90% of the calibration ceiling |
| f50 | 50% of the calibration ceiling |
| f25 | 25% of the calibration ceiling |

There are 24 LRC-BART configurations per endpoint. Use **default, requested ESS 30** as the main reference. Report the remaining 23 configurations per endpoint in sensitivity tables.

Comparators:

- KM, unadjusted.
- AFT-CP, with full RWD likelihood weight.
- Standard BART.
- Hierarchical AFT with both prior settings 0.05 and 0.5.

Both PFS and OS are included in the reporting set. There are 29 method/configuration rows per endpoint, 58 in total. Five-year RMST estimates by arm, ratios, differences, and intervals are saved. Table 7 reports PFS patient-profile prior ESS and the UCMM count matching each displayed age-group and categorical-covariate combination; OS remains in the application outcome tables.

Source references: [application README](../../lrcbart-case-study-mm/README.md), [lrcBART.R](../../lrcbart-case-study-mm/lrcBART.R), and [summarize.R](../../lrcbart-case-study-mm/summarize.R).

## Shared naming and interpretation rules

- Use `default` for the default configuration, with no added "config" label. Older saved Gaussian two-arm metadata may say `standard`; its reporting label is `default`.
- In simulations, 100, 75, and 50 denote requested ESS counts. In the application, those saved labels denote percentages of the 30-person trial, giving requested ESS 30, 22.5, and 15. Tables must make the actual requested ESS explicit.
- f90, f50, and f25 denote 90%, 50%, and 25% of the calibration ceiling, not percentages of the historical sample size. The ceiling can vary by simulated dataset or application calibration setting.
- Requested targets can be capped at the calibration ceiling. Do not describe requested targets as achieved or realized borrowing without distinguishing the relevant saved quantity.
- s0min uses the minimum supplied s_0^2 grid value, 10^-6, and bypasses ESS-target selection. It is a near-pooling prior benchmark, not exact complete pooling.
- Keep the LRC-BART mixture probability w distinct from a comparator's RWD likelihood weight. Both single-arm LRC-BART w values are fixed; the default two-arm w is updated.
- The two-arm w0, w1, and Hg10 sensitivities are outside the current reporting set because those configurations have no saved result files in the current folders.
- Do not merge prior ESS measures with MCMC effective sample sizes or present one as the other.

## Verified result coverage and remaining work

The inspected simulation reporting set comprises 204 Gaussian two-arm, 156 survival two-arm, 18 Gaussian single-arm, and 36 survival single-arm result files. All 414 contain complete distinct replicate IDs 1-100. The saved simulation files contain alternative-hypothesis runs only; Type I error requires additional null results.

The remaining table decisions include exact row/column layouts, main-text scenario subsets and full supplementary coverage, performance measures, survival-estimand priority, and application endpoint placement. No numerical findings or completed manuscript tables are settled by this configuration decision. Previously recorded model, ESS, and RMST implementation differences remain open.
