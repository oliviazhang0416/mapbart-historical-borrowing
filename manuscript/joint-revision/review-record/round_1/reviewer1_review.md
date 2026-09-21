# Round 1 methodologist review

## Summary

The manuscript is careful about the scope of its current control-source theory, its flat-baseline risk benchmark, and the distinction between calibration ESS and the fitting prior. Those qualifications are strengths. In its present form, however, it is a separate-treatment paper: the likelihood in `main_critical_final.tex`, Eq. (1), fits (f_{m trt}) independently, while the borrowing model is learned from the two control sources. The requested primary analysis is instead the joint randomized-trial model

\[
y_i=f(x_i)+D_i g(x_i)+A_i\psi+\epsilon_i,
\qquad \epsilon_i\sim N(0,\sigma^2_{g_i}),
\]

with (D_i=1) for both randomized arms and (A_i) indicating treatment. That change is substantive. It changes the estimand, likelihood, variance model, full conditionals, prior-scale provenance, ESS interpretation, theoretical assumptions, and the evidence that can support the main conclusion. The manuscript cannot be accepted as a joint-treatment paper by relabeling the present treatment surface.

The documented Web Appendix H is a useful primary evidence source for the revised paper. It gives the joint coefficient update, uses a shared randomized-trial variance and an external-control variance, and retains all 200 datasets and convergence flags. Its pooled constant-effect improvement is modest (4.64% relative to separated source-BART and 6.9% relative to joint source-BART); the paired interval against separated source-BART includes zero, and the discrepancy chains have very poor diagnostics. These results support a qualified exploratory claim, not a general joint-model superiority claim. The fixed-region control theorems and the flat-baseline control risk theorem remain valid as benchmarks under their stated assumptions, but they do not establish a joint treatment-effect guarantee.

## What the paper does well

- The control-source construction is stated coherently, including the distinction between (f), (f+g), and the total discrepancy. The exact Gaussian leaf factors and the collapsed tree-move argument in `source/appendix_computation.tex`, Sections A.2--A.5, are useful building blocks for a joint derivation.
- The manuscript explicitly says that its regional results do not automatically transfer to learned forests, shared hyperparameters, data-dependent scales, or censored outcomes (`main_critical_final.tex`, lines 246--248). That restraint should be preserved when the treatment arm is added.
- The ESS section correctly treats its quantity as information about a specified current-control mean and distinguishes the untruncated calibration law from the truncated fitting law (`main_critical_final.tex`, lines 154--174; `source/appendix_ess.tex`, Sections B.5 and B.7).
- Web Appendix H supplies the essential joint Gaussian algebra: with (\beta_{\rm trt}\sim N(0,100)), (V_\beta=(100^{-1}+\sum_i A_i^2/\sigma_{g_i}^2)^{-1}) and mean (V_\beta\sum_i A_i r_i/\sigma_{g_i}^2), followed by treatment-adjusted tree residuals. Its explicit statement that the calibration is a prior-scale procedure rather than an ESS=100 interpretation for the joint model is scientifically appropriate.
- `03-theory/appendix_new.tex`, Sections D.2--D.3, makes the important distinction between conditional Gaussian precision calculations and marginal ELIR for the joint ensemble distribution. That material can safely support a clearly labeled information diagnostic, provided it is not presented as a validated joint-model ESS or an operating-characteristic result.

## Major comments

### 1. The primary likelihood and estimand are not the required joint-treatment model

**Location and quote.** `main_critical_final.tex`, Eq. (1), lines 31--42, states
“(Y_{\mathrm{trt},i}\sim N\{f_{\mathrm{trt}}(x_{\mathrm{trt},i}),\sigma_{\mathrm{trt}}^2\})” and then defines the contrast from the separately fitted (f_{\mathrm{trt}}). Lines 60--68 explicitly call the treatment fit “separate from the control model.” The same separation is used in `source/appendix_survival.tex`, Sections D.1--D.4, and throughout the application.

**Problem and consequence.** This likelihood does not let treated outcomes inform the common treatment coefficient or the joint decomposition of the randomized-trial mean into (f+g) and treatment effect. It therefore does not define or analyze the requested primary model. The present Gaussian and survival claims, ESS statements, and application estimands refer to a different analysis. In particular, replacing (f_{\rm trt}) by a symbol (\psi) without changing the residual updates would leave the scientific target and posterior unchanged.

**Required fix.** Rewrite the primary model and estimand around (y=f+Dg+A\psi+\epsilon). Define (D=1) for every randomized observation, distinguish the external-control mean (f), randomized-control mean (f+g), and treated mean (f+g+\psi), and state the overlap and randomized-assignment assumptions needed for a treatment interpretation. If the primary model uses a common (\psi), its standardized ATE is (\psi); if a heterogeneous (\psi(x)) is intended, give that function its own model and estimand. Move the current separate-treatment survival and single-arm analysis to an explicitly labeled variant. The abstract, contribution paragraph, figures/tables, discussion, and generated results must be rewritten so that old separate-treatment evidence is not used as evidence for the joint primary model.

### 2. The joint sampler, variance model, and prior-scale provenance must be derived and made internally consistent

**Location and quote.** `source/appendix_computation.tex`, Eqs. (A.2)--(A.5), updates (f) from controls only, (g) from randomized controls only, and then states at lines 132--133 that “the treatment sum of trees is fitted separately.” The main text at line 68 gives a separate (\sigma_{\rm trt}^2). By contrast, Web Appendix H.2 uses one variance for both randomized arms, a separate external variance, a (N(0,100)) common treatment coefficient, and subtracts (A_i\beta_{\rm trt}) before every tree update.

**Problem and consequence.** These are different posteriors. In the joint model, the (f)-tree likelihood uses all observations after removing the current (Dg) and (A\psi) contributions; the (g)-tree likelihood uses all randomized observations after removing (f) and (A\psi); and the treatment coefficient is coupled to both the randomized treated and randomized control residuals through the weighted normal update. A separate treatment variance, a treatment-only range scale, or control-only partial residuals silently changes the target. The choice between a shared randomized-arm variance and arm-specific variances also changes the likelihood and the treatment full conditional.

**Required fix.** Add the complete joint conditional algebra. For a common coefficient (\psi\sim N(0,v_\psi)), the conditional update must be
\[
 V_\psi=(v_\psi^{-1}+\sum_i A_i^2/\sigma_{g_i}^2)^{-1},\qquad
 m_\psi=V_\psi\sum_i A_i\{y_i-f(x_i)-D_i g(x_i)\}/\sigma_{g_i}^2.
\]
State the exact residuals and marginal leaf factors for both tree ensembles, including the treatment subtraction, and state the variance updates for the selected shared or arm-specific structure. Document centering and whether the (N(0,100)) number is a variance. Reconcile every terminal-node and residual-variance range rule with the actual joint fitting outcomes. Web Appendix H says its joint range rule uses all three groups, while the current appendix gives separate historical and treatment ranges; both cannot describe one primary implementation. This is a manuscript fix, but the executable source and generated outputs must then be re-audited before any results are retained.

### 3. The current theorems do not transfer to the joint primary model; a bounded fixed-region theorem is feasible, but no BART or risk guarantee follows

**Location and quote.** `main_critical_final.tex`, lines 187--248 and `source/appendix_theory.tex`, Sections C.1--C.4, assume only historical and randomized controls with means (f_j) and (f_j+g_j). The text correctly says that the full learned-ensemble result remains an objective. Neither theorem contains treated observations or a treatment coefficient.

**Problem and consequence.** The current learning theorem can support a control-benchmark statement, and the flat-baseline risk theorem is explicitly a one-region current-control benchmark. Neither proves consistency of (\psi), identification of (g) in the presence of (\psi), or an ATE risk improvement. Identification fails if randomized controls are absent: within a region, treated and external means identify (f) and (f+g+\psi), but (g) and (\psi) can be shifted in opposite directions without changing the likelihood. The theorem must therefore make concurrent randomized controls and region overlap explicit.

**Required fix.** Add either a new theorem or an explicit benchmark boundary. A directly supportable fixed-region result is: in each fixed region (j), observe independent means (f_j) from external controls, (f_j+g_j) from randomized controls, and (f_j+g_j+\psi) from treated randomized subjects; assume all three region counts grow, residual variances are known and positive, and conditional on each finite spike/slab state the prior precision is fixed and positive definite. For parameter vector (\vartheta=(f_1,\ldots,f_J,g_1,\ldots,g_J,\psi)), the design rows are (e_{f_j}), (e_{f_j}+e_{g_j}), and (e_{f_j}+e_{g_j}+e_\psi). They have full column rank when all three counts are positive. With (G_n=X^TWX), (Q_z) the fixed prior precision, (V_z=(Q_z+G_n)^{-1}), and (m_z=V_zX^TWy), (V_z\to0), the bias is (-V_zQ_z\vartheta^*\to0), and (V_zG_nV_z\preceq V_z\to0). A finite mixture over (z) then gives posterior concentration in probability and (B_{c,j}\to\mathbb{I}(|g_j^*|<c)) away from (|g_j^*|=c), while (\psi\) is consistent. This is a fixed-region Gaussian result only; it does not establish learned-tree posterior concentration, data-dependent-scale behavior, censored-outcome validity, or risk improvement. If the authors do not add and prove this bounded result, they should retain the current theorems as control benchmarks and remove any implication that they analyze the joint treatment effect.

### 4. The evidence for the revised primary model is currently absent from the manuscript and is not strong enough for an unqualified superiority claim

**Location and quote.** The main simulation design and results (`main_critical_final.tex`, lines 262--299; `source/appendix_simulation.tex`, Sections E.1--E.5) describe 100-replicate studies with a separately fitted treated surface. The primary joint evidence is instead in the separate Web Appendix H, which reports 40 datasets per condition, five conditions, and joint LRC, joint source-BART, and separated source-BART.

**Problem and consequence.** The old Gaussian/survival tables cannot validate the new primary model. Web Appendix H is suitable evidence, but it reports modest gains and serious discrepancy-chain instability: pooled constant-effect RMSE reduction is 4.64% versus separated source-BART with a paired interval ([-0.011089,0.002817]), and 6.9% versus joint source-BART; discrepancy ESS minima are about 4 and maximum (\widehat R) values are about 2. The ATE convergence flags affect 18/200 joint-LRC fits, 12/200 joint-source fits, and 21/200 separated fits. The heterogeneous condition also uses a covariate-varying treatment truth while the documented joint fit uses a common (\beta_{\rm trt}), so it is a deliberate misspecification/sensitivity condition rather than evidence for a correctly specified heterogeneous treatment model. The null results are too small to establish nominal error control.

**Required fix and remaining empirical work.** Make Web Appendix H the primary joint-treatment experiment, retain every dataset and convergence flag, and state the reductions as exploratory and comparator-specific. Do not reclassify failed chains or substitute the old 100-replicate results. Label the old separate-treatment Gaussian, survival, and single-arm studies as a separate variant, or remove their use in the primary conclusion. If the paper seeks a broad claim about heterogeneous joint treatment effects, additional correctly specified treatment-function experiments and targeted diagnosis of the flagged chains remain empirical work; they cannot be replaced by prose changes. No new performance claim should be made from the current evidence.

### 5. The existing ESS is a control-prior calibration, not a joint-treatment ESS

**Location and quote.** `main_critical_final.tex`, lines 96--174, defines ESS for the current-control mean (\mu=f+g) and calibrates (s_0^2) with (w=1), an untruncated scale law, (V_f) from an external fit, and the randomized-control variance. `source/appendix_ess.tex`, Sections B.5--B.7, explicitly records the differences between calibration and fitting. Web Appendix H.2 likewise states that its training-only calibration is “not an ESS=100 interpretation for the joint model.”

**Problem and consequence.** In the joint model, treated outcomes inform (\psi), and the treatment contrast is coupled to the control surfaces through the shared likelihood and variance. Calling the selected scale “ESS 100” for the joint ATE would attribute information to the prior that the current definition does not measure. The quantity also depends on whether randomized treatment and control share a variance, on the standardization population, and on data-dependent range and residual scales. The D.2--D.3 ELIR derivation is mathematically useful but does not by itself produce a joint ESS or a performance guarantee.

**Required fix.** State prominently that the existing control ESS is a prior-setting reference for (f+g), not a joint ATE ESS. If a joint information number is reported, define its target, reference likelihood, conditioning/mixing law, and treatment-coefficient contribution, then derive or label it as a diagnostic. Safe incorporation of D.2--D.3 is as an appendix distinction between conditional working precision and marginal ELIR; it should not be used to justify an ESS=100 interpretation or to tune a threshold after seeing performance. Preserve the fixed-parameter ESS theorem as a benchmark and disclose all data-dependent scale choices.

For a useful fixed-forest Gaussian check, let (Z_0) be the design matrix for all (f)- and (g)-leaf coefficients, let (A) be the treatment column, (W) the source-precision matrix, and (Lambda_0) the conditional leaf-prior covariance. The joint conditional precision is
\[
\begin{pmatrix}
\Lambda_0^{-1}+Z_0^TWZ_0 & Z_0^TWA\\
A^TWZ_0 & v_\psi^{-1}+A^TWA
\end{pmatrix}.
\]
Consequently, after integrating the forest coefficients, the conditional marginal variance of (\psi) is
\[
\left[v_\psi^{-1}+A^TWA-A^TWZ_0(\Lambda_0^{-1}+Z_0^TWZ_0)^{-1}Z_0^TWA\right]^{-1},
\]
whereas ((v_\psi^{-1}+A^TWA)^{-1}) is only the full-conditional variance given the current forest. The latter is smaller and must not be reported as posterior precision or a joint ESS. Mixing over tree states, indicators, scales, and hyperparameters adds another total-variance term. This identity is a safe appendix diagnostic, not a joint-ensemble ESS theorem.

## Minor comments

1. **Notation.** The manuscript uses source labels, (D), (A), and (g_i) in several incompatible roles across the current theory and the joint documentation. Define one source indicator (S_i) or (g_i), set (D_i=\mathbb I\{S_i=\mathrm{RCT}\}), and reserve (A_i) for treatment throughout. The control-only theorem may use a different symbol in its local notation.

2. **Treatment estimand.** The current contrast in Eq. (2) is a difference of two independently estimated surfaces. Under a common (\psi), write the standardized ATE and its target directly; under (\psi(x)), specify the treatment-function prior and the trial-profile average. Do not use the old (f_{\rm trt}-f-g) formula as if it were the joint estimand.

3. **Identification statement.** Add a short rank/overlap condition to the methods and state that (n_{1,j}>0) is needed to separate (g_j) from a common (\psi) in the fixed-region benchmark. External support alone does not identify the treatment effect or the discrepancy outside the observed randomized-control support.

4. **Variance notation.** The main model has (\sigma_{\rm trt}^2,\sigma_1^2,\sigma_2^2), whereas Web Appendix H shares the randomized-arm variance. State one primary choice and update it consistently in Gaussian and censored-log-time calculations.

5. **Prior scale units.** Write (\psi\sim N(0,100)) as a variance statement and distinguish it from a standard deviation of 100. Explain whether the outcome centering and range rule are applied before or after the joint treatment term is introduced.

6. **Simulation labels.** The present tables call all two-arm results primary even though they are generated and fitted under the separate-treatment model. Add method/version labels to prevent the old generated tables from being read as joint evidence.

7. **Heterogeneity condition.** In Web Appendix H, identify the heterogeneous-treatment setting as misspecified for the common-(\beta_{\rm trt}) fit, and report it as a stress test. It cannot establish performance for a heterogeneous joint model.

8. **Single-arm scope.** A single-arm treated sample still cannot identify (g) or jointly separate (g) from a treatment effect. Retain the prior-driven qualification, but make clear that the single-arm variant does not validate the concurrent-control joint primary model.

9. **Risk language.** Keep “strict MSE improvement” attached to the flat-baseline, one-region control benchmark. Do not carry it into the joint ATE discussion without a new risk calculation.

10. **Diagnostics.** Report treatment-coefficient diagnostics separately from tree/discrepancy diagnostics. A satisfactory scalar ATE (\widehat R) does not establish exploration of (g), as the Web Appendix H discrepancy diagnostics demonstrate.

## Suggestions to enhance the paper

Add a one-page model-comparison table distinguishing the joint primary model from the separate-treatment variant by likelihood, estimand, variance structure, priors, calibration target, and supporting evidence. Put the joint conditional updates and the bounded fixed-region theorem in a new computation/theory appendix. In the results, lead with Web Appendix H’s five conditions and its explicit uncertainty, then preserve the older Gaussian, survival, and single-arm results as a clearly labeled secondary variant. The discussion should say that concurrent randomized controls identify (g) only where their support overlaps the external data and that the present joint evidence is exploratory.

## Source-dependent flags

- This review verifies the mathematical specifications documented in the manuscript, Web Appendix H, `source/appendix_*`, and `03-theory/appendix_new.tex`; it does not certify an executable sampler. Before using any joint result, verify in the source that (D_i=1) is applied to both randomized arms, the same (\beta) full conditional is used in every tree and variance update, the selected randomized-variance structure is actually sampled, and the treatment-adjusted residual is passed to all (f)- and (g)-tree moves.
- The current generated tables and figures describe the separate-treatment study. They must not be relabeled as joint results without a hash-preserving provenance check against the Web Appendix H outputs and the joint configuration.
- Web Appendix H states that its range prior uses all three groups in a joint fit and that its calibration excludes treated outcomes. Those two choices are defensible, but they must be stated as conditioned-on prior-setting choices rather than silently merged with the current appendix’s separate treatment range and residual-scale rules.
- The fixed-region theorem above assumes known positive variances and fixed prior hyperparameters. It does not cover the fitted BART ensemble, estimated (w), estimated (\tau_0^2), data-dependent scales, survival censoring, or treatment-function heterogeneity. Any broader theorem would require additional proof.
- D.2--D.3’s ELIR calculations require the joint mixing law and retained cross-tree dependence. Independently recombining fitted draws or inserting the control ESS calibration into that formula would change the target distribution.
- In the fixed-forest identity above, the Schur-complement subtraction is nonnegative, so the marginal treatment variance is at least the full-conditional variance. Any implementation or table using the full-conditional quantity as if it integrated out (f,g) is a source-dependent correctness error.

## Verdict

The paper has a potentially publishable joint-treatment extension, but the current manuscript is not yet a valid presentation of that extension. The necessary manuscript corrections are concrete: replace the primary likelihood and estimand, derive the joint sampler, qualify the ESS, add or clearly separate a joint fixed-region theorem, and move the old evidence to a variant. The modest and diagnostically unstable Web Appendix H results also require restrained claims; a broad joint-method superiority statement would require further correctly specified empirical evidence and targeted computational diagnosis.

What bothered me most is that the manuscript’s strongest theory and most polished simulation narrative currently describe a separate treatment fit, while the requested primary scientific claim is about a jointly identified treatment coefficient.

VERDICT: MAJOR REVISION
MAJOR COUNT: 5
