# Summary

The revision has repaired the central narrative problem identified in Round 1. The title, abstract, model, estimand, posterior computation, primary simulation, and discussion now describe the joint Gaussian treatment model. The single-arm myeloma analysis and the earlier Gaussian/survival simulations are explicitly identified as separate-treatment variants. The fresh study is also described with the correct five condition names, shared seeds for the first three conditions, source-specific outcome construction, comparator specifications, and retained results for all 200 datasets.

The paper is now reviewable as a focused methodological manuscript. Its current evidence supports a qualified demonstration of the proposed formulation, not a general performance claim. The remaining obstacles for a top statistics journal are substantive: 18 of 200 Joint LRC method--dataset groups meet the stated ATE diagnostic flag, the minimum ATE bulk ESS is about 35, and discrepancy diagnostics have bulk ESS near 4 with R-hat near 2. In addition, the three methods differ in trial-arm variance sharing, range inputs, calibration and mean structure, so the modest pooled RMSE differences cannot identify which component produced them. The manuscript acknowledges both limitations, but the evidence and claims should be made conditional on them at the level of the primary conclusion.

# What the paper does well

- The contribution is now intelligible: a separately regularized source-discrepancy forest is embedded in a likelihood that jointly uses treated trial subjects, concurrent controls and external controls. The introduction also correctly says that a source indicator in one BART forest can represent regional differences, so the distinction is not incorrectly reduced to “two forests versus one.”

- The identification boundary is unusually clear. The manuscript states that Delta = psi only under the common-effect restriction and gives the g(x) + b, psi - b likelihood invariance when concurrent controls are absent.

- The joint validation is reproducible on paper. Appendix G gives the 300-trial/300-external design, source probability, intercept, regions, heterogeneous and null truths, seeds, chain lengths, calibration stages, and the retained method--dataset metrics.

- The results are candid. The five-condition table retains the unfavorable heterogeneous and null results, reports coverage counts and ATE flags, and the text says that the 10% criterion was not met. The discussion also states that the comparisons do not isolate a benefit caused solely by joint treatment modeling.

- The application is now properly scoped. It reports the current 30 EloKRd and 200 UCMM cohort and calls the analysis a separate-treatment single-arm illustration; it does not present it as evidence for the joint treatment model or as a joint survival analysis.

# Resolution check

1. **Resolved — Round 1 Major 1, primary model absent.** The title and abstract now name the joint method, and the introduction states the joint contribution (`main_critical_final.tex`, lines 5--24). The primary model is defined as $Y_i=f(x_i)+S_i g(x_i)+A_i psi + epsilon_i$ with both randomized arms sharing the trial surface (lines 26--38). The primary simulation is included immediately after the fixed-region result (`source/joint_results_main.tex`, lines 1--26), while the separate-treatment survival and single-arm material is explicitly moved to a variant section (main, lines 218--231).

2. **Resolved — Round 1 Major 2, estimand and identification.** The main text defines the trial-standardized ATE, states Delta = psi under the common-effect restriction, identifies the heterogeneous condition as a deliberate violation, and states the no-concurrent-control invariance (main, lines 40--47). Appendix G repeats the group means and the realized ATE target and gives the individual heterogeneous effect and its trial-profile mean (lines 4--17 and 91--99). This is sufficient to prevent the heterogeneous condition or the single-arm analysis from being read as evidence for an unrestricted treatment-effect model.

3. **Partially resolved — Round 1 Major 3, comparator confounding.** Appendix G now provides the requested comparator table, including mean model, trial variance, range inputs and calibration (`source/appendix_joint_study.tex`, lines 70--85), and the main result says that the comparison concerns complete fitted formulations (`source/joint_results_main.tex`, lines 14--20). The discussion explicitly declines to attribute the result solely to joint treatment modeling (main, lines 279--284). The remaining problem is empirical: no matched-variance or matched-prior comparison isolates the discrepancy forest, shared trial variance, or the joint treatment update. The manuscript should therefore keep the pooled reductions as descriptive evidence for these complete formulations rather than as evidence for a mechanism-specific gain.

4. **Partially resolved — Round 1 Major 4, condition-specific evidence and convergence.** The five-condition table, coverage counts, bootstrap qualification, and ATE flags now appear in the primary results (`source/joint_results_main.tex`, lines 16--26; `source/generated/joint_main_table.tex`). Appendix G additionally reports the maximum R-hat and minimum bulk ESS for both ATE and discrepancy summaries (lines 105--140). This resolves the presentation failure. It does not resolve the underlying evidence limitation: the Joint LRC ATE flag affects 18/200 groups, its minimum bulk ESS is 35.2, and discrepancy diagnostics are much worse. A reader still needs the primary conclusion to be explicitly conditional on these unsettled chains, rather than treating the retained RMSE and coverage table as settled operating characteristics.

5. **Resolved — Round 1 Major 5, application and legacy evidence.** The application is now labelled as a separate-treatment illustration and says that the absence of concurrent controls precludes validation of the joint decomposition (main, lines 218--231). The text reports the current 30-plus-200 cohort (line 235), distinguishes the 100-replicate separate-treatment record from the primary study (line 219), and states that no joint-treatment survival analysis is reported (line 226). The earlier Round 1 wording that called this an older 30-plus-253 cohort was corrected: that cohort belongs only to earlier local revision artifacts.

# Major comments

## 1. The primary numerical evidence remains convergence-limited

**QUOTE OR LOCATION.** `source/appendix_joint_study.tex`, lines 107--140, says that ATE flags use R-hat > 1.05 or bulk ESS below 400, then reports 18/200 flagged Joint LRC groups, minimum bulk ESS 35.2, and discrepancy maximum R-hat near 2 with minimum bulk ESS near 4. `source/joint_results_main.tex`, lines 22--26, acknowledges that these diagnostics limit interpretation.

**ISSUE AND CONSEQUENCE.** The disclosure is good, but the primary table still presents RMSE, coverage and pooled paired comparisons as the validation evidence without a pre-specified analysis of how flagged method--dataset groups affect those summaries. A flag is not proof that a particular posterior mean is unusable, but 18 flagged Joint LRC groups and extremely weak discrepancy diagnostics are too frequent to treat the comparison as a settled operating-characteristic result. This is especially important because the apparent gain is modest and the 10% criterion is not met.

**CONCRETE FIX.** Before making a performance conclusion, either complete and independently verify the flagged ATE chains under a pre-specified rule, or add a clearly labelled flagged-chain sensitivity analysis with the exact denominator and rule retained. If no new chain analysis is authorized, make the primary table explicitly descriptive and state in the abstract, results and conclusion that the comparative evidence remains provisional because of the observed ATE diagnostics. Keep the discrepancy summaries explicitly diagnostic; do not use them as evidence of regional recovery.

## 2. The central comparative claim must remain about complete formulations

**QUOTE OR LOCATION.** `source/joint_results_main.tex`, lines 14--20, and `source/appendix_joint_study.tex`, lines 70--85, correctly state that Joint LRC, Joint source-BART and Separated source-BART differ in mean model, trial variance, range inputs and calibration. The main discussion (`main_critical_final.tex`, lines 280--288) also says that the comparison does not isolate a benefit caused solely by joint treatment modeling.

**ISSUE AND CONSEQUENCE.** This is now an explicit limitation, but it remains the main interpretive constraint on the proposed contribution. The pooled RMSE reductions of 4.64% and 6.9% cannot establish that the spike-and-slab discrepancy forest, the shared trial variance, or joint estimation of psi is responsible. Without a matched specification, the paper cannot support a mechanism-specific superiority claim even if all chains were stable.

**CONCRETE FIX.** Carry the “complete fitted formulations” qualification into the abstract and the first sentence introducing the pooled comparison. Describe the result as a comparator-specific pattern under the stated prior, calibration and variance choices. Either add a pre-specified matched-variance/prior sensitivity study before making a mechanism claim, or make the contribution explicitly a model proposal with bounded evidence rather than an empirically demonstrated improvement.

# Minor comments

1. **Abstract strength.** The abstract says that the results “identify circumstances in which joint modeling is useful” (main, line 11), although the paired comparison is modest, the 10% gate is unmet and the ATE diagnostics are unsettled. Change “identify” to “suggest” or “illustrate under the stated formulation.”

2. **Literature positioning.** The introduction’s discussion of Zhou and Ji and BCF is substantially more careful (main, line 20). Keep the claim limited to the explicitly regularized source-discrepancy prior within the joint trial likelihood and control-mean reference. Do not imply that source-indicator BART cannot represent regional discrepancies or that two forests alone establish novelty.

3. **Diagnostic denominator.** The table caption correctly says that flags count method--dataset groups rather than chains (`source/generated/joint_main_table.tex`). Repeat “method--dataset groups” in the main results sentence whenever the 18/200, 12/200 and 21/200 values are quoted, so readers do not read them as 18/600 chain failures.

4. **ATE versus discrepancy diagnostics.** The appendix distinguishes them, but the main result could state once that the ATE flags are the treatment-effect diagnostic and the near-2 R-hat/near-4 ESS values concern discrepancy contrasts. This prevents the two forms of instability from being conflated.

5. **Null interpretation.** The statement that 36--39 of 40 intervals cover in the five-condition table should remain paired with the existing sentence that 40 noisy null datasets do not establish nominal error control (`source/joint_results_main.tex`, line 18). Avoid describing these counts as calibration evidence.

6. **Heterogeneous target wording.** In `source/joint_results_main.tex`, line 12, call $1+0.5(X_5 - bar X_5,trial)$ an individual treatment effect and state that its trial-profile average is one. Appendix G already supplies this precise qualification.

7. **Calibration boundary.** Appendix G, lines 64--68, appropriately separates training-only calibration from fit-specific outcome ranges and residual scales. The main ESS discussion should retain the sentence that the requested count of 100 is a control-mean prior-setting reference, not treatment information (`main_critical_final.tex`, line 187), wherever the number 100 is mentioned.

8. **Application transition.** The separate-treatment survival/single-arm section is correctly scoped, but the transition after the joint Gaussian results should repeat that it is a different likelihood. This will prevent a reader scanning from the primary table directly into the application from treating the application as joint-model validation.

9. **Reproducibility access.** Appendix G, line 103, says that the original validation archive contains raw chains while the accompanying source package does not redistribute them. Identify the archive location and checksum in the final reproducibility materials, or state precisely what a reader can obtain from the supplied package.

10. **Scope of empirical generalization.** The conclusion appropriately calls for a randomized-trial application and larger null/heterogeneous studies (main, line 288). Keep this future-work sentence adjacent to the statement that the current Gaussian validation does not transfer to survival or single-arm analyses.

# Suggestions to enhance the paper

- Retain the current five-condition table as the central numerical display and make the paired pooled comparison visibly secondary to it.

- A compact graphical model linking $f$, $g$, $S$ and A psi would help readers distinguish source membership from treatment assignment without adding methodological claims.

- If the authors do not add a matched comparator study, title the empirical result as a validation of the complete implementations and reserve “benefit” for the compatible-condition descriptive pattern.

- Preserve the explicit limitation that the single-arm application demonstrates standardized external-control prediction, not identification of psi or the untreated trial-control discrepancy.

# Source-dependent flags

- The primary-source note verifies the limited Zhou--Ji and BCF positioning, but it is not an exhaustive novelty search. Any broader claim about priority or closest precedent still requires a complete literature review.

- The deidentified merged n230 application RData is locally available, but the raw 797-patient UCMM source remains private and 58 posterior result RData files are absent from this checkout. The application assessment therefore remains dependent on the available merged data, generated JSON/manifests and source scripts; no new fit was run.

- The primary validation is Gaussian with a common treatment coefficient. It supplies no direct evidence for a joint survival model or for treatment-effect identification in the single-arm application.

- All 200 datasets and all three method rows are retained, but the observed convergence flags and missing matched-mechanism sensitivity limit the strength of any claim about general superiority.

# Verdict

The revision resolves the Round 1 coherence failures and now presents a scientifically honest, bounded joint-model paper. For a top statistics journal, the remaining convergence limitation and the inability to attribute the small performance difference to a specific model component require another substantive revision or an explicit downgrade of the empirical claim. The manuscript should be acceptable for further review after those evidence and wording boundaries are handled, but the current record does not support acceptance as a settled performance contribution.

VERDICT: major revision
MAJOR COUNT: 2

What bothered me most: the paper is now candid about unstable chains, yet its main validation table still invites readers to treat those same flagged draws as settled comparative evidence.




