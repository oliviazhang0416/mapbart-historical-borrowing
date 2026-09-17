# Handoff: changes from the revision manuscript

17 September 2026

This note summarizes the changes from `mapbart-historical-borrowing-pr1/revision/05-writing/main_draft.tex` and its included theory and results sections to the current manuscript in `lrcbart-historical-borrowing/manuscript/`. It distinguishes changes in exposition from changes in the reported analyses. The revision manuscript remains unchanged.

## 1. Organization and writing

The draft has been rewritten in the style of a statistical methods submission, with standard statistical and BART terminology. The main text follows Introduction, Methods, Simulation studies, Multiple-myeloma application, and Discussion. Detailed calculations, proofs, implementation qualifications and complete results are organized into Appendix A–F.

The revision assembled the main article, displays and web appendices into one document. The current version separates them into three documents, each with its editable `.tex` file alongside the PDF:

- `main`: the article, including the main equations and three conditional/calibration theorems.
- `main_tables_figures`: seven tables and Figures 1--2.
- `appendix`: computation, ESS derivation, conditional theory, survival calculations, and complete simulation/application results.

The approved explanatory material is retained in `content_migration/`. These notes are the basis of the redraft, rather than additional sections to insert verbatim. The handoff is a description of manuscript changes, not a study-running guide.

## 2. Model specification: clarify the data groups and conditioning

The underlying historical-plus-discrepancy construction is retained. The revision introduced the control model as `y_i = f(x_i) + D_i g(x_i) + error`. The new main text writes the three outcome models explicitly: historical controls have mean `f(x)`, trial controls have mean `f(x)+g(x)`, and trial-treated patients have a separately fitted mean `h(x)`. Individual outcomes, covariate vectors, sample sizes and residual variances use subscripts `1` (trial treatment), `2` (trial control) and `3` (historical control), with separate within-group indices. Trial-standardized estimands average over `\mathbf{X}_1` and `\mathbf{X}_2`, with one contribution per participant. This also makes clear that the control-only equation is not the complete model for all trial participants.

The revised explanation distinguishes the following:

- Both control sources inform `f` in the final two-arm fit; only trial-control partial residuals directly inform `g`.
- The preliminary historical fit used for ESS is different: its likelihood contains historical outcomes only.
- The latent `z_jl` is a terminal-node mixture indicator, not an observed source indicator. The shared `w` is a prior mixture probability, not a historical-patient likelihood weight.
- The `f` and `g` trees have separate structures. Sharing candidate variables or cutpoints does not make the discrepancy trees inherit the historical trees.

The single-arm explanation is brief in the main text. Its target control surface remains `f+g`; without trial controls, the discrepancy receives no outcome-likelihood information. This property was already present in the revision. What is now made explicit is the current implementation: the reported single-arm/application fits fix `w` and fix `tau_0^2 = min(s_0^2, tau_1^2/2)`, whereas the revision's description treated the discrepancy scale as a prior draw. Calibration still averages over an untruncated scale law. Those specifications must not be described as identical.

## 3. ESS: rebuild the argument and distinguish the averaging operations

The ESS section now begins with its reference experiment. Independent normal observations with known residual variance and a flat prior on their mean give a posterior mean variance of `sigma_2^2/m`. Matching that variance to uncertainty about the historical-data-informed target control mean defines an equivalent count `m`; averaging those conditional counts over the calibration scale law then defines ESS. This makes explicit that both sides concern uncertainty about a mean, not individual outcome variability versus mean uncertainty.

Definitions now use `\coloneqq`, whereas derived equalities retain `=`. Bold `\mathbf{X}_a` and `\mathbf{Y}_a` distinguish group data collections from scalar simulation predictors. Notation has been simplified: tree structures and terminal-node parameters are written as `\mathcal T_j^f, \theta^f_{jl}` and `\mathcal T_j^g, \theta^g_{jl}`, and the separate treatment surface is `h(x)`. `V_f` replaces `V_mu^f`, `V_g` denotes the specified prior variance of the standardized discrepancy mean, and `ESS_tau0` replaces `ESS_0` to identify the all-spike calibration reference. Conditioning is written using `\mathbf{Y}_3`, `\mathbf{X}_3`, and `\mathbf{X}_1,\mathbf{X}_2`. The derivation explains where the leaf proportions `a_jl` come from and retains covariance among predictions when calculating the variance of a subgroup or overall mean.

There is an important distinction from the revision's general definition. Its Appendix D first defined a working ESS by averaging reciprocal conditional variances over leaf indicators and tree structures. The new general definition first integrates the indicators and trees **inside the discrepancy variance**, and then averages the matched counts over the spike variance **outside the reciprocal**:

\[
V_g(w,\tau_0^2)
=\{w\tau_0^2+(1-w)\tau_1^2\}
 E_{\mathcal T^g}\!\left[\sum_{j,\ell}a_{j\ell}^2\right],
\qquad
\mathrm{ESS}(w,s_0^2)
=E_{\tau_0^2}\!\left[\frac{\sigma_2^2}{V_f+V_g(w,\tau_0^2)}\right].
\]

These orders of averaging are generally different. For the all-spike reference, the resulting formula retains the revision's implemented expected-tree-weight coefficient before inversion; that numerical convention is not newly introduced. The new walkthrough derives it through the variance definition rather than leaving the relationship implicit.

The ceiling `sigma_2^2/V_f` is retained and derived as the zero-spike-scale limit. The text also distinguishes all-spike calibration from allowing slab assignments in the final fit. It does not claim a global maximum-borrowing property over an untruncated scale law that can exceed the slab variance.

Several qualifications already appeared in the revision and are retained: this is not exact marginal ELIR; subgroup ESS values are not additive; and survival reference units are uncensored normal log-time observations. The redraft develops their interpretation more directly. It also separates the full-variance mathematical ESS from the saved block-average selector, finite grid and ceiling override. The current numerical selector has not been replaced as part of manuscript drafting.

## 4. Posterior computation and theoretical scope

The computation walkthrough has been shortened and reorganized around partial residuals, marginal leaf factors, the Metropolis–Hastings ratio, and conditional parameter updates. Both `M_f` and `M_g` now use the same precision-sum notation `A` and `B`. The density-ratio representation for `M_g` and the auxiliary residual-mean variance notation have been removed. This is an algebraic and presentational rewrite, not a change to the tree acceptance calculation.

The text states explicitly that discrepancy tree proposals integrate out both the terminal-node parameter and `z_jl`, conditional on the shared hyperparameters. After a tree is retained, including after rejection of a proposed move, its indicators and node parameters are sampled conditionally. Appendix A provides the normal integral, sampling steps and acceptance probability.

The theoretical presentation is more compact and narrower. Calibration is now Theorem 1; bounded additional shrinkage relative to the slab-only leaf estimate is Theorem 2; conditional learning of a leaf discrepancy is Theorem 3. The revision's extended fixed-partition ensemble results, detailed detection-threshold discussion and block-rule asymptotics have not all been carried into the new theorem statements. This is a reduction in scope, not a new proof of full-model robustness. The revision already qualified its theory; the redraft preserves the distinctions between leaf inference, latent-component selection, and inference with unknown tree structures.

## 5. Simulation reporting: use the current study results

The numerical results have been replaced by summaries from the four current simulation folders, rather than copied from the revision's tables. The revision's main two-arm study reported 500 replicates; the current completed results contain 100 replicates per setting. The current reporting set comprises 414 result files and 606 performance rows after including supplementary survival RMST estimands.

The main differences collaborators should track are:

| Item | Revision manuscript | Current manuscript |
|---|---|---|
| Scenario numbering | Sc4 regional; Sc5 global | Sc4 global; Sc5 regional |
| Two-arm reporting | Included additional nested-model and tuning comparisons | Primary target 100; targets 75, 50, f90, f50 and f25 in sensitivity results |
| Extra experiments | CAHB transfer, development/oracle analyses, and `H_g=10`, `w=0`, `w=1` comparisons | Not carried into the current empirical comparison |
| Single-arm reporting | Single-arm results and sensitivities in the revision study | Gaussian `n_1=200`; survival `n_1=30,200`; Sc1/Sc2; five lrcBART configurations plus four comparators |
| Minimum-scale label | Earlier full-borrowing terminology | `s0min`, explicitly a small-scale prior benchmark |
| Decision frequency | Described through interval exclusion of the null | Reports the actual saved, uncalibrated decision rules; Gaussian benefit threshold 0.5, survival ratio threshold 1, with method-specific rules documented |

The single-arm lrcBART configurations are fixed `w=1` with target 100, f90 and s0min, plus fixed `w=0.9` with target 100 and f90. The two-arm primary configuration retains `H_f=50`, `H_g=5` and an updated `w` under a Beta(1,1) prior. MAP and PSCL are reported for Gaussian two-arm simulations, not added to the survival comparison.

The new tables include Monte Carlo standard errors in the appendix and distinguish requested ESS, capped targets, evaluated curve ESS and ceilings. The available runs are alternative-hypothesis simulations; the draft does not claim to provide null Type I error or calibrated-power results. Design-specific generators, comparator implementations and survival numerical conventions are documented from the current code. Differences from the revision's numerical estimates must not be interpreted as a paired demonstration of improved performance or as a change caused solely by rewriting the manuscript.

## 6. Application: retain the cohorts, change the reporting and local ESS target

The application retains 30 EloKRd patients and 253 UCMM historical controls. Harmonization was already part of the revision. The current analysis reports **harmonized coding only**, removing the original-coding sensitivity comparison. PFS and OS five-year RMST results remain, with arm summaries, ratios, differences and intervals. The main lrcBART setting remains `H_f=10`, `H_g=5`, fixed `w=1`, and requested ESS 30. The appendix reports all 58 method/endpoint rows, including the `H_f=50` and `w=0.9` sensitivities. Application labels 100/75/50 correspond to counts 30/22.5/15; they are absolute counts in the simulations.

Data-construction reporting is more explicit: the UCMM selection flow, harmonized category values, differing age reference dates, endpoint definitions, and the high-risk-cytogenetics fallback are described. The postinduction nature of ASCT and descriptive interpretation of the comparison are retained qualifications, not newly discovered features. The draft also avoids equating the available 30-person extract with the full published trial enrollment without clinical confirmation.

The application results were refreshed from the rerun case-study outputs. Table 6 and Table S10 distinguish the hypothetical EloKRd control RMST, the historical-model RMST standardized to EloKRd covariates, and the observed unadjusted UCMM RMST. KM reports only the observed UCMM quantity; adjusted methods report the two model-based control quantities. Updated contrasts, sensitivity ranges and diagnostics replace the earlier numerical values. Simulation results and the preliminary subgroup ESS calculation were retained.

### Figure 2: preserve the revision's visual style, update its inputs

Figure 2 retains the revision's two panels against `X_5`: the probability of a small discrepancy and posterior mean discrepancy, with a shared scenario legend, distinct point shapes, explicit region labels and dashed truth references. It now uses the first 40 replicates of the current Gaussian primary analyses for three regional shifts (0.5, 1 and 2) in current Sc5. The global-shift curve is omitted. Profile-level summaries were recovered with the saved datasets, scales and seeds and checked against the saved aggregate results. The old figure's numerical curves were not copied.

The original MAPBART manuscript's Figure 3 DAGs are migrated as Figure 1, illustrating the three simulation source-selection mechanisms. The DAG precedes the regional plot in order of first citation. Its structure is retained, with notation defined locally in the caption and a reference in the simulation-design text.

### Revision Figure 2: replace it with Table 7

The revision's Figure 2 showed pointwise ESS at 30 trial profiles against separate predictors, with both PFS and OS. It used an inverse-variance plug-in with the fitted mean spike variance in the denominator. That display is removed.

Table 7 instead reports PFS ESS for the **mean control log-time within joint covariate subgroups**, using the settled expectation-of-conditional-ratios definition. Its final design is:

- Columns: age band, sex, race, Hispanic/Latino status, cytogenetic risk, ASCT, EloKRd count, UCMM count, and historical prior ESS.
- Age bands: under 50, 50–59, 60–69, and 70+. Exact ages remain in model predictions.
- Rows: the 22 combinations observed in EloKRd, ordered hierarchically by the displayed covariates. Combinations not represented in EloKRd are omitted, including historical-only combinations.
- Numerical ESS: the 22 combinations represented in EloKRd, reported to one decimal place. Every displayed row has an ESS. Displayed counts sum to 30 EloKRd and 98 UCMM patients; all 253 UCMM patients remain in the historical fit.

Counts describe patients within the combination. Each ESS uses information from the **full historical cohort**, not just the UCMM patients counted in that row. ESS values range from 5.4 to 23.4 and do not generally sum to overall ESS. OS remains in the outcome tables but is omitted from this local-information table.

Hierarchical comparator labels are `hier. LM` and `hier. AFT`; prior settings appear in table footnotes rather than method-name parentheses. Where both settings are reported, table footnotes identify the first row as 0.05 and the second as 0.5 within each reporting block.

Comparator display labels use -NP for trial-only controls, -CP for complete pooling, and -PP for source-indicator-adjusted pooling. In single-arm settings, -CP denotes use of all historical controls without a discrepancy term. These are display changes; fitted methods and numerical results are unchanged.

## 7. Methodological qualifications

The manuscript has been revised using the current project’s results. The existing analysis code, saved study results and cleaned datasets were not changed during redrafting.

One methodological issue remains: the implemented ESS calibration does not fully match the settled derivation. Differences concern the use of block-specific rather than full-variance variances, the spike-variance distribution, and the tree-generation rules. The present manuscript retains the implementation and describes these differences explicitly. Aligning it with the derivation would require a separate methodological revision and rerunning affected analyses.

The editorial review also makes the modified supplementary simulation RMST functionals explicit: the two-arm summaries substitute a residual standard deviation under the stated rule, and simulation ratios floor the control denominator. These are preserved results, not merely integration approximations. The primary median-survival analyses and application analytic RMST summaries are unaffected. Author details, funding, conflicts of interest and the code/data-access statement remain to be supplied.
