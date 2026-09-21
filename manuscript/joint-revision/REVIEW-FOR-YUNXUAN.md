# Joint LRC-BART revision for Yunxuan's review

Yunxuan, we have revised the manuscript around the joint treatment model while retaining your compact main-paper organization, separate tables and detailed appendix. We would like your scientific review before deciding which additional analyses to run for submission.

## Reading files

| Document | Clean | Tracked |
|---|---|---|
| Main paper, 13 pages | [PDF](main_critical_final.pdf) · [LaTeX](main_critical_final.tex) | [PDF](main_tracked.pdf) · [full source patch](main_tracked.patch) |
| Main tables, 4 pages | [PDF](main_tables_figures.pdf) · [LaTeX](main_tables_figures.tex) | [PDF](main_tables_figures_tracked.pdf) · [full source patch](main_tables_figures_tracked.patch) |
| Appendix, 76 pages | [PDF](appendix.pdf) · [LaTeX](appendix.tex) | [PDF](appendix_tracked.pdf) · [full source patch](appendix_tracked.patch) |

The tracked files compare against your commit `51b269deff1bd3dccf79b655092feaa242e8078c`. Table and appendix redlines use additions-only color markup because deletion markup broke table alignment; their full source patches retain all deletions and mathematical changes. The parent manuscript files remain available as the baseline.

## What changed

The primary Gaussian model is

`Y_i = f(x_i) + S_i g(x_i) + A_i psi + epsilon_i`.

Here `S=1` for both randomized arms and `A=1` only for treated trial participants. We use a common treatment effect and shared trial residual variance. Both trial arms inform the trial response surface through treatment-adjusted residuals. The main text and appendix now follow that likelihood in the estimand, priors, posterior updates and interpretation.

The primary numerical comparison is the joint Gaussian experiment with 200 datasets: 40 each under compatible controls, a one-variable source shift, a two-variable source shift, heterogeneous effects and the null. The earlier 100-replicate studies and the n230 single-arm survival application remain separate-treatment variants. We have not fitted a joint survival model or a heterogeneous treatment-effect forest in this revision.

## Review priorities

1. **Joint model and identification — Section 2 and Appendix A.** Are the shared trial surface, common-effect assumption, source variances and treatment-adjusted updates stated correctly? Without concurrent controls, a constant can move between `g` and `psi` without changing the likelihood; please check that boundary and the single-arm interpretation.
2. **Theory — Section 4 and Appendix C.** Please examine the fixed-region joint-learning theorem and its proof. The result conditions on fixed priors and known variances; it does not prove consistency of learned BART partitions or a treatment-risk improvement. Your control-risk result is retained as a separate benchmark.
3. **ESS — Section 3 and Appendix B.** We retain the control-mean prior information reference and distinguish it from treatment-effect information under the joint model. Please assess the Schur-complement explanation, the marginal ELIR discussion, and how prominently we should address differences between the calibration and fitting priors.
4. **Empirical evidence — Section 5 and Appendix G.** Please assess whether the current conclusions match the results. The pooled RMSE reductions are 4.64% against separated source-BART and 6.9% against joint source-BART across Compatible/One-variable/Two-variable, with no general advantage across all conditions. Joint LRC has ATE diagnostic flags in 18/200 method–dataset groups; discrepancy diagnostics are weaker. Saved-chain sensitivity preserves all datasets but does not establish convergence. The comparisons also differ in trial variance sharing and prior scales.
5. **Paper scope and next analyses — Sections 6–8.** Should the separate-treatment survival illustration remain in the main paper or move to the supplement? For a top statistics journal, our priorities are resolving flagged fits, matched variance/prior comparisons, calibration alignment, stronger null/heterogeneity assessment and an application that directly evaluates the joint model. Please help prioritize these and assess the contribution relative to prior work.

## Checks and reproducibility

The existing numerical check recomputes all 15 method-condition RMSE and coverage rows, Monte Carlo errors and ATE flag counts from the 600 canonical rows. The supplied source-only bundle and checksums are in [reproducibility/joint-validation](reproducibility/joint-validation/README.md). No new MCMC fit was run during this manuscript rewrite. The original joint source/output freeze was verified unchanged. Application results were reconciled to the n230 saved summaries; missing original posterior files prevent an independent application refit.

The clean and tracked PDFs were compiled and visually checked. The package contains internal LLM-assisted editing audits for transparency; they are not peer review and do not close the remaining scientific questions. Author fields are still blank.

Build the clean package with `python3 source/scripts/build_review_package.py` from this directory. Run the numerical summary check with `python3 reproducibility/joint-validation/verify_joint_validation.py`.
