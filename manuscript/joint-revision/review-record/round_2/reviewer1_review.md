# Round 2 methodologist review

## Summary

The revision addresses the main Round 1 defect. The primary manuscript now specifies
\[
Y_i=f(x_i)+S_i g(x_i)+A_i\psi+\epsilon_i,
\]
with both randomized arms sharing a trial variance, a separate external variance, and a common treatment effect. The estimand, treatment full conditional, treatment-adjusted forest updates, fixed-region theorem, Schur-complement distinction, and separate-treatment boundary are now stated consistently. The fixed-region proof is valid under its explicit fixed-design, known-variance, fixed-prior assumptions. The globally restricted tree weight is also distinguished from the forced-terminal calibration recursion.

Two issues remain material for a top statistics journal. The central joint empirical comparison still pools all method--dataset groups while 18/200 Joint LRC groups have an ATE convergence flag, and the paper does not show the existing chain-index sensitivity of the pooled comparison. The primary trial residual-scale prior is also estimated from a within-trial regression that omits the treatment indicator, so its empirical-Bayes provenance should be justified or made a sensitivity. These are bounded and implementable; they do not require reopening the resolved joint theorem or relabeling the control ESS as treatment information.

## What the paper does well

- The primary data structure now identifies the roles of external controls, concurrent controls, and treated participants. The text correctly states that concurrent controls are needed to separate (g) from a common (\psi), and that the heterogeneous condition violates the common-effect restriction.
- The treatment update in Eq. (\(\ref{eq:psi_update}\)) correctly uses only treated observations directly because (A_i=0) for controls, while uncertainty in (f+g) is propagated through the joint posterior. The appendix correctly distinguishes this full-conditional variance from marginal treatment uncertainty.
- The fixed-region theorem and proof in `source/appendix_theory.tex`, Section C.joint, give a sufficient full-rank design condition and a valid finite-mixture concentration argument. The proof does not claim tree-boundary recovery, ensemble consistency, or treatment-risk improvement.
- The ESS section now explicitly makes the control-mean quantity a prior-scale reference, supplies the treatment/forest Schur complement, and separates marginal ELIR from the working variance calibration.
- The primary Gaussian study is now separated from the older survival and single-arm likelihoods. Its modest, comparator-specific findings and poor discrepancy diagnostics are reported without a general-superiority claim.

## Resolution check

1. **Round 1 Major 1 — Resolved.** The main likelihood and estimand are joint; the separate-treatment survival and single-arm material is labeled as a variant in `main_critical_final.tex`, Section `sec:variant`.

2. **Round 1 Major 2 — Partially resolved.** The joint conditionals, shared trial variance, all-group range rule, treatment prior variance, and source-specific variance updates are documented. The remaining issue is the treatment-omitting preliminary regression used to set the trial residual-scale prior; see Major comment 1 below.

3. **Round 1 Major 3 — Resolved.** The fixed-region joint theorem and proof include (f_j), (g_j), and (\psi), require all three group types in each region, and state the limits of the result. The control risk theorem is explicitly retained as a benchmark.

4. **Round 1 Major 4 — Partially resolved.** Web Appendix G is now the primary joint experiment, all conditions and flags are retained, and superiority is qualified. The remaining concern is whether the flagged-chain empirical summaries are stable enough for the displayed pooled comparison; see Major comment 2.

5. **Round 1 Major 5 — Resolved.** The manuscript states that the control ESS is not a joint treatment ESS, gives the Schur complement, and presents ELIR as a distinct diagnostic.

## Major comments

### 1. The empirical-Bayes trial variance prior is documented but its treatment-scale provenance remains unresolved

**Location and quote.** `source/appendix_computation.tex`, Section A.1, states that the joint wrapper estimates the source residual scale “by regressing centered outcomes on covariates within source” and that “this preliminary fit does not remove the treatment coefficient.” The primary model has one trial variance for treated and control observations (`main_critical_final.tex`, lines 38 and 86).

**Problem and consequence.** For the randomized source, a covariate-only preliminary regression leaves the treatment mean difference in its residuals. The resulting empirical-Bayes scale is therefore partly driven by the treatment effect and allocation, while the final likelihood estimates (\psi) jointly. This is a defensible conditioned-on prior choice, but it is not a neutral residual-variance estimate and can change shrinkage, interval widths, and the tree likelihood. The fixed-region theorem conditions on known variance and does not justify this data-dependent scale choice.

**Required fix.** State explicitly that the trial variance prior scale is an empirical-Bayes choice based on a misspecified treatment-omitting preliminary regression, and explain why it is preferred. At minimum, add a sensitivity using a preliminary regression that includes (A_i), or report a source-verification result showing that the chosen scale is practically unchanged. If no such sensitivity is available, retain the choice but move its effect into the limitations and avoid describing the fixed-region result as validating the fitted variance hierarchy. This is a prior-provenance correction, not a request for a new treatment model.

### 2. The primary ATE comparison remains computationally provisional and should show the existing chain-sensitivity evidence

**Location and quote.** `source/joint_results_main.tex` reports 18/200 Joint LRC method--dataset convergence flags, a minimum ATE bulk ESS of about 35, and says that “unsettled treatment-effect and discrepancy chains limit the conclusions.” It nevertheless reports pooled RMSE and paired bootstrap comparisons over all datasets. The Discussion repeats that the numerical record is unsettled (`main_critical_final.tex`, lines 279--288).

**Problem and consequence.** The point estimates and intervals from flagged groups remain part of the primary comparison. The reported paired interval quantifies replicate resampling uncertainty, not Monte Carlo uncertainty from these flagged chains. Without showing how the existing retained chains change the pooled direction, the reader cannot tell whether the 4.64% or 6.9% reductions are stable features of the saved runs. This is especially relevant because the discrepancy diagnostics are much worse than the ATE diagnostics.

**Required fix.** Use the existing chain-specific summaries to add a compact sensitivity display: report the pooled constant-condition reduction for each chain index, the condition-specific flag counts, and state that these are descriptive chain-sensitivity checks rather than convergence proof. Keep all datasets and the original pooled rows; do not select a converged subset or replace the primary estimates. If the authors do not add this display, the results section should call the pooled reductions provisional numerical patterns and avoid using them as evidence of a stable performance advantage.

## Minor comments

1. The theorem’s “all three group counts tend to infinity” condition is sufficient rather than minimal. Keep that wording, and make clear that the result is for a common (\psi), not for the heterogeneous stress condition.

2. The theorem proof uses the smallest group-by-region precision to lower-bound (X^T W X). State once that variances are known and fixed in this argument; estimated residual variances are only part of the fitted implementation.

3. Appendix C.4 still describes how a separate independent treatment estimator inherits a control-risk improvement under exact agreement. Label this explicitly as a separate-treatment contrast benchmark so it cannot be read as a risk result for the joint (\psi).

4. The global restricted tree weight is now correctly stated, including the (1-p_k(d)) factor for an exhausted terminal node. Add one sentence that its normalized marginal tree law is not the ordinary recursively normalized Galton--Watson prior, and that the calibration simulator intentionally targets a different forced-terminal recursion.

5. The abstract says “several treatment-effect chains” are unsettled, while the tables report method--dataset groups flagged by an aggregate over three chains. Use one denominator and label the unit consistently.

6. The comparator table appropriately notes that trial-variance sharing, ranges, and calibration differ. The Results should repeat that these differences prevent attributing the pooled RMSE pattern solely to joint treatment modeling.

7. The ELIR appendix correctly requires the joint mixture law. Make the distinction from the control ESS visible at the first point where ELIR is mentioned, rather than relying on Appendix B.9 for the qualification.

8. The separate-treatment survival appendix says to use the Gaussian updates in Appendix A. Point the reader to the separate-treatment subsection specifically, since the first part of Appendix A now describes the joint sampler.

## Suggestions to enhance the paper

Add a short table beside the primary results with three columns: estimand, model likelihood, and computational qualification. This would keep the joint Gaussian evidence, separate-treatment survival evidence, and control-mean ESS reference from being conflated. A one-paragraph description of the existing chain-index sensitivity is sufficient; no new fit or reclassification is needed. The conclusion should retain the current narrow claim: the joint formulation is a coherent primary model with a modest, comparator-specific Gaussian pattern, while treatment-effect computation, discrepancy recovery, and joint ESS remain open questions.

## Source-dependent flags

- The review assumes the source verification recorded in the Round 1 audit: (S_i=1) for both randomized arms, treatment-adjusted (f)- and (g)-residuals, the shared trial variance, and the exact treatment full conditional. These should remain checked against the executable source after any further edit.
- The tree prior is now an explicit globally restricted finite-support weight. Its MH ratios are coherent for that target, while the calibration simulator uses forced termination when no split remains. This is correctly disclosed as a calibration mismatch and should not be silently described as one normalized recursive prior.
- The joint validation table retains flagged groups and discrepancy diagnostics. Any chain-sensitivity display must use the saved chains and preserve the current rows; deleting or reclassifying flagged groups would change the estimand of the empirical summary.
- The theorem assumes fixed prior hyperparameters and known positive variances. It does not cover the empirical-Bayes residual scales, estimated (w) and (\tau_0^2), learned tree partitions, censoring, or heterogeneous treatment effects.

## Verdict

The revision resolves the primary model, identification, theorem-scope, and ESS-interpretation defects from Round 1. It is close to a defensible top-journal methods presentation, but the empirical-Bayes trial variance provenance and the stability of the central ATE comparison need explicit treatment. Both can be addressed with bounded manuscript and saved-output checks; no new broad simulation program is required for this review.

What bothered me most is that the revised scientific claims are now appropriately cautious, but the displayed pooled ATE improvement still rests on flagged chains without showing the already-available chain-index sensitivity.

VERDICT: MAJOR REVISION
MAJOR COUNT: 2
