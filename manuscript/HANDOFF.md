# Handoff: changes from the revision manuscript

18 September 2026

This note summarizes the changes from `mapbart-historical-borrowing-pr1/revision/05-writing/main_draft.tex` and its included theory and results sections to the current manuscript in `lrcbart-historical-borrowing/manuscript/`. It distinguishes changes in exposition from changes in the reported analyses. The revision manuscript remains unchanged.

## 1. Organization and writing

The draft has been rewritten in the style of a statistical methods submission, with standard statistical and BART terminology. The main text follows Introduction, Methods, Simulation studies, Multiple-myeloma application, and Discussion. Detailed calculations, proofs, implementation qualifications and complete results are organized into Appendix A–F.

The revision assembled the main article, displays and web appendices into one document. The current version separates them into three documents, each with its editable `.tex` file alongside the PDF:

- `main`: the article, including the main equations and three conditional/calibration theorems.
- `main_tables_figures`: seven tables and Figures 1--2.
- `appendix`: computation, ESS derivation, conditional theory, survival calculations, and complete simulation/application results.

The approved explanatory material is retained in `content_migration/`. These notes are the basis of the redraft, rather than additional sections to insert verbatim. The handoff is a description of manuscript changes, not a study-running guide.

As of 19 September 2026, `main.pdf` is 12 pages and `appendix.pdf` is 56 pages. Both build without undefined references or overfull boxes. The appendix retains one known duplicate PDF-destination warning caused by splitting Table S2 across two `longtable` environments; the visible table is complete. The ESS comparison HTML and the model, computation, ESS and application walkthroughs in `content_migration/` have been synchronized with the manuscript notation.

## 2. Model specification: clarify the data groups and conditioning

The underlying historical-plus-discrepancy construction is retained. The revision introduced the control model as `y_i = f(x_i) + D_i g(x_i) + error`. The current main text writes all three outcome models explicitly and uses group labels `trt`, `1` and `2` for randomized treatment, randomized control and RWD control. Their conditional means are `f_trt(x)`, `f(x)+g(x)` and `f(x)`, with residual variances `sigma_trt^2`, `sigma_1^2` and `sigma_2^2`. Trial-standardized estimands average over `\mathbf{X}_{trt}` and `\mathbf{X}_1`, with one contribution per participant; in a single-arm trial, `n_1=0`. This makes clear that the control-only equation is not the complete model for all trial participants.

The revised explanation distinguishes the following:

- Both control sources inform `f` in the final two-arm fit; only trial-control partial residuals directly inform `g`.
- The preliminary historical fit used for ESS is different: its likelihood contains historical outcomes only.
- The latent `z_{h^g,l}` is a discrepancy terminal-node mixture indicator, not an observed source indicator. The shared `w` is a prior mixture probability, not a historical-patient likelihood weight.
- Tree indices are `h^trt`, `h^f` and `h^g` for the `f_trt`, `f` and `g` ensembles. Sharing candidate variables or cutpoints does not make the discrepancy trees inherit the prognostic trees.

The single-arm explanation is brief in the main text. Its target control surface remains `f+g`; without trial controls, the discrepancy receives no outcome-likelihood information. This property was already present in the revision. What is now made explicit is the current implementation: the reported single-arm/application fits fix `w` and fix `tau_0^2 = min(s_0^2, tau_1^2/2)`, whereas the revision's description treated the discrepancy scale as a prior draw. Calibration still averages over an untruncated scale law. Those specifications must not be described as identical.

## 3. ESS: rebuild the argument and distinguish the averaging operations

The ESS section now begins with the information-matching principle, using ELIR as motivation. Prior ESS expresses prior information in units of observations from a reference likelihood. For the proposed definition, target-prior information is the inverse conditional marginal variance of the RWD-informed prior for the current-control mean, while one normal randomized-control reference observation supplies Fisher information `1/sigma_1^2`. Equating these quantities defines a conditional equivalent count `m`; averaging those counts over the spike-variance calibration law defines prior ESS.

Definitions use `\coloneqq`, whereas derived equalities retain `=`. Bold `\mathbf{X}_a` and `\mathbf{Y}_a` distinguish group data collections from scalar simulation predictors. Tree structures and terminal-node parameters use `\mathcal T_{h^f}^f, \theta^f_{h^f\ell}` and `\mathcal T_{h^g}^g, \theta^g_{h^g\ell}`; the separate treatment surface is `f_trt(x)` with tree index `h^trt`. `V_f` is the posterior variance of the RWD-control mean standardized to `\mathbf{X}_1`, `V_g` is the corresponding discrepancy-prior variance, and `V_pi=V_f+V_g` is the RWD-informed prior variance. The preliminary historical fit conditions on `\mathbf{Y}_2,\mathbf{X}_2` and is evaluated at `\mathbf{X}_1`. The derivation explains the profile weights `a_{h^g\ell}` and retains covariance among predictions when calculating an overall mean.

Appendix B now places the proposed definition beside Neuenschwander's variance calibration, Morita--Thall--M\"uller's curvature match and ELIR. It uses `p(Y\mid\theta_\star)` for a generic sampling density so that generic likelihood notation is not confused with the manuscript's prognostic surface `f`. Neuenschwander's complete-pooling row is presented as a calibration anchor rather than as a hypothetical-likelihood reference. Morita's reference is the expected curvature of a posterior formed from a low-information prior and hypothetical data; ELIR uses one-observation Fisher information. The proposed method uses the same one-observation reference-information concept as ELIR but measures RWD-informed-prior information by inverse conditional marginal variance.

There is an important distinction from the revision's general definition. Its Appendix D first defined a working ESS by averaging reciprocal conditional variances over leaf indicators and tree structures. The new general definition first integrates the indicators and trees **inside the discrepancy variance**, and then averages the matched counts over the spike variance **outside the reciprocal**:

\[
V_g(w,\tau_0^2)
=\{w\tau_0^2+(1-w)\tau_1^2\}
 E_{\mathcal T^g}\!\left[\sum_{h^g,\ell}a_{h^g\ell}^2\right],
\qquad
\mathrm{ESS}(w,s_0^2)
=E_{\tau_0^2}\!\left[\frac{\sigma_1^2}{V_f+V_g(w,\tau_0^2)}\right].
\]

These orders of averaging are generally different. For the all-spike reference, the resulting formula retains the revision's implemented expected-tree-weight coefficient before inversion; that numerical convention is not newly introduced. The new walkthrough derives it through the variance definition rather than leaving the relationship implicit.

The ceiling `sigma_1^2/V_f` is retained and derived as the zero-spike-scale limit. The text also distinguishes all-spike calibration from allowing slab assignments in the final fit. It does not claim a global maximum-borrowing property over an untruncated scale law that can exceed the slab variance.

Several qualifications already appeared in the revision and are retained: this is not exact marginal ELIR; patient-profile ESS values do not sum or average to overall ESS; and survival reference units are uncensored normal log-time observations. The redraft develops their interpretation more directly. It also separates the full-variance mathematical ESS from the saved block-average selector, finite grid and ceiling override. The current numerical selector has not been replaced as part of manuscript drafting.

## 4. Posterior computation and theoretical scope

The computation walkthrough has been shortened and reorganized around partial residuals, marginal leaf factors, the Metropolis–Hastings ratio, and conditional parameter updates. Both `M_f` and `M_g` now use the same precision-sum notation `A` and `B`. The density-ratio representation for `M_g` and the auxiliary residual-mean variance notation have been removed. This is an algebraic and presentational rewrite, not a change to the tree acceptance calculation.

The text states explicitly that discrepancy tree proposals integrate out both the terminal-node parameter and `z_{h^g\ell}`, conditional on the shared hyperparameters. After a tree is retained, including after rejection of a proposed move, its indicators and node parameters are sampled conditionally. Appendix A provides the normal integral, sampling steps and acceptance probability. The posterior spike probability is denoted `gamma(rbar)` and the additional mixture shrinkage by `Delta_mix(rbar)`, avoiding collisions with the RWD-informed prior `pi` and information notation `I`.

The theoretical presentation is more compact and narrower. Calibration is now Theorem 1; bounded additional shrinkage relative to the slab-only leaf estimate is Theorem 2; conditional learning of a leaf discrepancy is Theorem 3. The revision's extended fixed-partition ensemble results, detailed detection-threshold discussion and block-rule asymptotics have not all been carried into the new theorem statements. This is a reduction in scope, not a new proof of full-model robustness. The revision already qualified its theory; the redraft preserves the distinctions between leaf inference, latent-component selection, and inference with unknown tree structures.

## 5. Simulation reporting: use the current study results

The numerical results have been replaced by summaries from the four current simulation folders, rather than copied from the revision's tables. The revision's main two-arm study reported 500 replicates; the current completed results contain 100 replicates per setting. After excluding the discontinued single-arm hierarchical setting $s_\tau^2=0.5$, the current reporting set comprises 408 result files and 596 performance rows, including supplementary survival RMST estimands.

The main differences collaborators should track are:

| Item | Revision manuscript | Current manuscript |
|---|---|---|
| Scenario numbering | Sc4 regional; Sc5 global | Sc4 global; Sc5 regional |
| Two-arm reporting | Included additional nested-model and tuning comparisons | Primary target 100; targets 75, 50, f90, f50 and f25 in sensitivity results |
| Extra experiments | CAHB transfer, development/oracle analyses, and `H_g=10`, `w=0`, `w=1` comparisons | Not carried into the current empirical comparison |
| Single-arm reporting | Single-arm results and sensitivities in the revision study | Gaussian `n_trt=200`; survival `n_trt=30,200`; Sc1/Sc2; five lrcBART configurations plus three comparators |
| Minimum-scale label | Earlier full-borrowing terminology | `s0min`, explicitly a small-scale prior benchmark |
| Decision frequency | Described through interval exclusion of the null | Reports the actual saved, uncalibrated decision rules; Gaussian benefit threshold 0.5, survival ratio threshold 1, with method-specific rules documented |

The single-arm lrcBART configurations are fixed `w=1` with target 100, f90 and s0min, plus fixed `w=0.9` with target 100 and f90. The two-arm primary configuration retains `H_f=50`, `H_g=5` and an updated `w` under a Beta(1,1) prior. MAP and PSCL are reported for Gaussian two-arm simulations, not added to the survival comparison.

The new tables include Monte Carlo standard errors in the appendix and distinguish requested ESS, capped targets, evaluated curve ESS and ceilings. The available runs are alternative-hypothesis simulations; the draft does not claim to provide null Type I error or calibrated-power results. Design-specific generators, comparator implementations and survival numerical conventions are documented from the current code. Differences from the revision's numerical estimates must not be interpreted as a paired demonstration of improved performance or as a change caused solely by rewriting the manuscript.

## 6. Application: retain the cohorts, change the reporting and local ESS target

The application retains 30 EloKRd patients and 200 verified triplet-treated UCMM historical controls. Harmonization was already part of the revision. The current analysis reports **harmonized coding only**, removing the original-coding sensitivity comparison. PFS and OS five-year RMST results remain, with arm summaries, ratios, differences and intervals. The main lrcBART setting remains `H_f=10`, `H_g=5`, fixed `w=1`, and requested ESS 30. The appendix reports all 58 method/endpoint rows, including the `H_f=50` and `w=0.9` sensitivities. Application labels 100/75/50 correspond to counts 30/22.5/15; they are absolute counts in the simulations.

Data-construction reporting is more explicit: the UCMM selection flow, harmonized category values, differing age reference dates, endpoint definitions, and the high-risk-cytogenetics fallback are described. The postinduction nature of ASCT and descriptive interpretation of the comparison are retained qualifications, not newly discovered features. The draft also avoids equating the available 30-person extract with the full published trial enrollment without clinical confirmation.

The application results were refreshed from the rerun case-study outputs. Table 6 and Table S10 report hypothetical EloKRd control RMST for adjusted methods and observed UCMM RMST only for KM. The redundant historical-model RMST standardized to EloKRd covariates was removed from adjusted-method outputs and tables. Updated contrasts, sensitivity ranges and diagnostics replace the earlier numerical values. Simulation results and the preliminary PFS calibration were retained.

### Figure 2: preserve the revision's visual style, update its inputs

Figure 2 retains the revision's two panels against `X_5`: the probability of a small discrepancy and posterior mean discrepancy, with a shared scenario legend, distinct point shapes, explicit region labels and dashed truth references. It now uses the first 40 replicates of the current Gaussian primary analyses for three regional shifts (0.5, 1 and 2) in current Sc5. The global-shift curve is omitted. Profile-level summaries were recovered with the saved datasets, scales and seeds and checked against the saved aggregate results. The old figure's numerical curves were not copied.

The original MAPBART manuscript's Figure 3 DAGs are migrated as Figure 1, illustrating the three simulation source-selection mechanisms. The DAG precedes the regional plot in order of first citation. Its structure is retained, with notation defined locally in the caption and a reference in the simulation-design text.

### Table 7: patient-profile PFS prior ESS

Table 7 reports one PFS prior ESS value for each of the 30 EloKRd patient covariate profiles. Its columns are age group, age, sex, race, Hispanic/Latino status, cytogenetic risk, ASCT, UCMM cohort count and prior ESS. Rows are ordered by age group and exact age, and every row displays its age group and UCMM count. Age is rounded to one decimal year for display; exact age is used in the calculation. A row's UCMM count is the number of triplet-treated historical controls matching its age group and five displayed categorical covariates.

For profile $x_i$, the target is the conditional mean control log-PFS $\mu_i=f(x_i)+g(x_i)$. The reference experiment uses the same parameter, $Y_{ij}^{\mathrm{ref}}\mid x_i\sim N(\mu_i,\sigma_1^2)$. The preliminary $H_f=10$ historical model supplies $V_f(x_i)$, and the all-spike discrepancy reference supplies $H_g\tau_0^2$. The reported value is $E_{\tau_0^2}[\sigma_1^2/\{V_f(x_i)+H_g\tau_0^2\}]$ under the globally selected PFS scale.

Every row is informed by the full 200-patient triplet UCMM cohort and uses the same residual variance and scale law. Values range from 3.8 to 24.5, are expressed in uncensored normal-reference units on the log-PFS scale, and do not measure compatibility with unobserved EloKRd controls. The pointwise values do not sum or average to the overall ESS. OS remains in the outcome tables.

Hierarchical comparator labels are `hier. LM` and `hier. AFT`; discrepancy-variance prior scales appear in table footnotes rather than method-name parentheses. Single-arm simulations and main Table 6 report $s_\tau^2=0.05$. The application sensitivity Table S10 retains both $s_\tau^2=0.05$ and $s_\tau^2=0.5$.

Comparator display labels use -NP for trial-only controls, -CP for complete pooling, and -PP for source-indicator-adjusted pooling. In single-arm settings, -CP denotes use of all historical controls without a discrepancy term. These are display changes; fitted methods and numerical results are unchanged.

## 7. Remaining limitations

The implemented prior-scale calibration does not fully coincide with the theoretical prior-ESS definition. The reported selector uses block-specific estimates of historical uncertainty, a finite scale grid and implementation-specific spike-variance and tree-generation rules. Consequently, the reported prior ESS should be interpreted as the information summary under the stated calibration procedure rather than as an exact description of the prior used in every final fit. Full alignment would require a methodological revision and rerunning the affected analyses.

The supplementary simulation RMST summaries use each arm's own fitted residual standard deviation in both the two-arm and single-arm workflows, together with a 0.05 lower bound on control RMST denominators. These conventions affect only the supplementary RMST summaries; the primary median-survival analyses and the application's analytic RMST summaries do not use them. The saved two-arm results predate the own-scale setting and must be regenerated before Table S5 and the manuscript PDFs are refreshed.

The simulations contain 100 replicates per setting and evaluate alternative-hypothesis scenarios only. They therefore provide limited precision for comparisons between methods and do not establish Type I error control or decision-rule calibration. The conditional theoretical results apply within a fixed discrepancy-tree leaf and do not establish robustness or posterior consistency for the full sum-of-trees model. In the single-arm setting, trial outcomes do not identify the untreated-control discrepancy, so conclusions remain dependent on the discrepancy prior, covariate support and equal-residual-variance assumption.

## 8. Outstanding submission items

Author information, funding, conflicts of interest and the code/data-access statement remain to be supplied before submission.
