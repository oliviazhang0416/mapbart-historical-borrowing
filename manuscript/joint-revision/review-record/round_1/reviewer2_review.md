# Summary

The manuscript has a potentially useful practical problem: an RCT control
surface and an external-control surface can agree in some covariate regions and
disagree in others, so a single borrowing weight is hard to interpret. The
LRC-BART control decomposition, local discrepancy summaries, exact leaf
marginals, and the distinction between working prior ESS and the fitting
distribution are potentially valuable contributions.

The current working manuscript is not yet reviewable as a paper centered on
the newly selected joint treatment model. Its main model and most of its
simulation and application narrative still use a separately fitted treatment
surface, while the joint model appears as an additional Gaussian experiment in
`05-writing/sections_1_2.tex`, `sections_3_5.tex`, and Web Appendix H. The
joint experiment is informative but exploratory: the pooled gain is modest,
the prespecified 10% criterion is not met, and treatment-effect convergence
checks fail for 18/200 joint-LRC fits. The heterogeneous and null conditions
also favor the separated baseline over joint LRC. The application remains a
single-arm, separate-treatment illustration and therefore cannot validate the
joint treatment model or identify an untreated control discrepancy.

The paper can become a coherent contribution, but it requires a substantial
reorganization around the joint estimand and an explicit separation between
the joint randomized-trial evidence and the legacy separate-treatment
simulation/application record.

# What the paper does well

- The local borrowing problem is concrete and important. The manuscript
  distinguishes prognostic adjustment from residual source discrepancy and
  gives practitioners a compatibility map rather than only a global borrowing
  scalar.
- The control model is described with useful operational detail: source-specific
  outcomes, shallow discrepancy trees, a spike-and-slab leaf prior, residual
  variance handling, minimum leaf sizes, and a reproducible calibration
  procedure.
- The manuscript is unusually candid about theoretical scope. It states that
  the fixed-partition results do not establish consistency for the full
  posterior over learned tree partitions and that the single-arm discrepancy
  is prior-driven.
- Web Appendix H records the fresh joint design with 200 datasets, five
  conditions, shared replicate seeds for the three constant-effect conditions,
  30,000 burn-in and 12,000 retained draws, all three comparators, and the
  saved-chain diagnostics. This is a strong provenance base for a qualified
  exploratory result.
- The validation record reports the unfavorable conditions rather than hiding
  them: joint LRC is worse under heterogeneous effects and the null, and its
  two-variable discrepancy diagnostics are unstable. The paired bootstrap
  interval against separated source-BART includes zero.

# Major comments

## 1. The selected primary model is absent from the working manuscript’s main story

**QUOTE OR LOCATION.** In `05-writing/journal-revision/main_critical_final.tex`,
the data model has a separate `f_{\mathrm{trt}}(x)` and the text says “The
treatment surface $f_{\mathrm{trt}}$ is fitted separately.” It repeats that
“The treatment fit is separate from the control model.” The joint model instead
appears in `05-writing/sections_1_2.tex`, subsection “Joint treatment modeling
for randomized trials,” and is evaluated as an “additional” experiment in
`05-writing/sections_3_5.tex`.

**ISSUE AND CONSEQUENCE.** This is a scope failure, not a local wording issue.
The title, abstract, estimand, prior description, posterior computation,
simulation comparators, application, and discussion describe the separate
treatment model, while the requested primary method is
$y_i=f(x_i)+D_i g(x_i)+A_i\psi+\varepsilon_i$. A reader cannot tell which
model the paper proposes, which posterior is the contribution, or whether the
reported operating characteristics support the claimed method. It also makes
the contribution appear to be an optional add-on rather than the selected
method.

**CONCRETE FIX.** Reorganize the manuscript around the joint model from the
start. Define $D_i=1$ for every RCT subject (treated or control), $A_i=1$ only
for treated RCT subjects, and $D_i=0,A_i=0$ for external controls. State that
the RCT-control mean is $f+g$, the treated mean is $f+g+\psi$ under the common
effect model, and the standardized ATE is $\psi$. Move the old separate
treatment formulation to a clearly labelled comparator or legacy variant.
Rewrite the abstract, introduction, contribution, posterior updates, and
discussion so the 200-dataset joint experiment is the primary evidence. State
explicitly that survival and the single-arm application remain separate-model
extensions unless a joint survival model is actually fitted.

## 2. The joint estimand, treatment-effect assumption, and identification boundary need one coherent statement

**QUOTE OR LOCATION.** `05-writing/sections_1_2.tex` defines the joint model
and says “the standardized ATE is $\beta_{\rm trt}$,” while Web Appendix H
introduces the same coefficient as $\beta_{\rm trt}\sim\Normal(0,100)$.
The fresh design nevertheless includes heterogeneous individual effects
$1+0.5(X_{5i}-\bar X_{5,\mathrm{trial}})$ and a null condition. The old
`main_critical_final.tex` instead defines a contrast involving separate
$f_{\mathrm{trt}}$ and $f+g$ surfaces.

**ISSUE AND CONSEQUENCE.** The common-effect model identifies a scalar
coefficient under its structural restriction, but the heterogeneous condition
does not satisfy that restriction. Calling the fitted coefficient an ATE in
all five conditions is acceptable only if the target and misspecification are
made explicit; otherwise the reader may interpret the heterogeneous result as
evidence for a general treatment-effect model. The single-arm text also needs
to distinguish two boundaries: without concurrent controls, the untreated
trial control and discrepancy are not identified, and the joint coefficient
$\psi$ is not estimated by the joint RCT model because the necessary control
arm is absent.

**CONCRETE FIX.** Add a short estimand box before the model: define the
standardized trial ATE, explain that the primary joint likelihood imposes a
constant treatment effect $\psi$, and say that the heterogeneous condition
tests robustness to violation of that restriction while scoring the posterior
mean against the realized trial-profile ATE. Keep the single-arm statement
precise: treated outcomes alone do not identify either the untreated trial
control surface or a trial/external discrepancy; do not imply that they
identify $\psi$. Report the common-effect prior variance in the same notation
and units as the outcome.

## 3. The comparator differences confound the joint architecture with residual-variance handling

**QUOTE OR LOCATION.** Web Appendix H.2 states that both joint models share one
residual variance across the two randomized trial arms, whereas separated
source-BART estimates the treated-arm variance independently and keeps
source-specific control variances. `sections_3_5.tex` correctly says the
comparisons concern the “full fitted formulations,” but the surrounding result
language still presents the numbers as a gain from joint modeling.

**ISSUE AND CONSEQUENCE.** Joint LRC, joint source-BART, and separated
source-BART differ in more than whether treatment is fitted jointly. They also
differ in trial-arm variance sharing, range-rule inputs, treatment/control
prior scales, and the LRC-only calibration. The pooled RMSE comparison cannot
attribute a difference to the discrepancy architecture or to joint treatment
borrowing alone. This matters especially because the apparent gain is only
4.64% versus the separated baseline and disappears or reverses in several
conditions.

**CONCRETE FIX.** Put a compact comparator specification table in the main
text or Web Appendix H with, for each method, mean structure, treatment
parameterization, trial-arm variance structure, external variance, tree/range
inputs, and calibration. Change all interpretation to “full formulation” or
“joint LRC under these variance and prior choices,” unless matched-variance
sensitivity evidence is added. Do not claim that the experiment isolates the
benefit of joint treatment modeling. Keep the source-specific variance
qualification visible in the main results paragraph.

## 4. The fresh joint evidence is exploratory and must be presented condition by condition, with convergence failures attached to the estimates

**QUOTE OR LOCATION.** `sections_3_5.tex` reports pooled constant-effect RMSE
reductions of 4.6% and 6.9% and says the result is exploratory. The frozen
validation summary reports joint-LRC RMSE 0.1750, 0.2068, 0.2151, 0.2092,
and 0.2298 across Compatible, One-variable, Two-variable, Heterogeneous,
and Null conditions. Treatment-effect convergence checks fail in 18/200
joint-LRC fits; the minimum ATE bulk ESS is about 35, and 51 method-dataset
groups are flagged across the three methods.

**ISSUE AND CONSEQUENCE.** A pooled favorable number can be read as the main
result even though the method is worse than separated source-BART under
heterogeneous effects (0.2092 versus 0.1917) and the null (0.2298 versus
0.2100), and worse under the two-variable condition (0.2151 versus 0.2113).
The treatment estimate is not uniformly better, the primary 10% gate is not
met, and unsettled chains make a settled performance claim premature. The
current manuscript’s older 100-replicate simulations cannot be used to fill
this evidence gap because they evaluate the separate-treatment model.

**CONCRETE FIX.** Make the five-condition table the main joint result, with
coverage counts and convergence flags adjacent to each method or clearly
cross-referenced. Keep the pooled constant-effect comparison as a prespecified
summary, but label it secondary to the per-condition display and retain its
bootstrap interval. State in the abstract and conclusion that the evidence is
exploratory, the 10% criterion was not met, and no general treatment-effect
gain was shown. Do not call the regional discrepancy recovered: the saved
discrepancy diagnostics are unstable and the two-variable contrast is far
below the generated contrast.

## 5. The application and legacy simulation record cannot serve as direct evidence for the joint method

**QUOTE OR LOCATION.** The current journal-revision base describes a single-arm
EloKRd application with a separately fitted treatment surface and the current
30-plus-200 n230 cohort. The older 30-plus-253 cohort appears only in the
earlier local revision application and HR artifacts. The upstream current
empirical record has current RMST outputs in the n230 generated files. The
joint validation conclusion explicitly says that the current application
remains a separate-treatment variant and that the fresh joint batch is a
different empirical record.

**ISSUE AND CONSEQUENCE.** Presenting the application immediately after the
joint randomized-trial evidence suggests that it validates the joint method,
although it has no randomized control arm and no fitted $A_i\psi$ joint model.
It also risks mixing n283 RMST/ESS values with the n230 application outputs.
The reader needs to know what practical workflow is actually available for a
single-arm trial and which claims are unsupported by that example.

**CONCRETE FIX.** Label the application “legacy separate-treatment, single-arm
illustration” or move it to a clearly separated extension section. If it is
retained in the primary paper, use the current n230 cohort and generated RMST
outputs consistently, and state that it demonstrates the control-borrowing
component only; it is not joint-treatment evidence and cannot identify
$\psi$ or the untreated trial-control discrepancy. Keep the older n283
analysis as an archived provenance record. Make the 500-replicate original
simulation and the 100-replicate upstream current simulation visibly distinct
from the 200-dataset joint validation.

# Minor comments

1. **Notation collision.** In the joint sections, use a source indicator such
   as $D_i$ consistently and reserve $g(x)$ for the discrepancy. Define $A_i$
   before the first equation and state that treated RCT subjects have
   $D_i=A_i=1$, RCT controls have $D_i=1,A_i=0$, and external controls have
   $D_i=A_i=0$.

2. **Variance labels.** The joint model’s common randomized-arm variance and
   separate external variance should be stated in the model equation and the
   comparator table, not only in Web Appendix H. Avoid calling all three
   methods “source-specific variance” when the joint fits share the two trial
   arms.

3. **Prior-scale language.** The joint LRC calibration is a practical
   training-only scale procedure. Keep the explicit statement that it is not
   an ESS=100 interpretation for the joint posterior; do not let the old
   two-arm ESS language imply otherwise.

4. **Convergence denominators.** Report the denominator for every failure rate
   (18/200 method-dataset groups, not 18/600 chains) and distinguish the
   ATE diagnostics from the much worse discrepancy diagnostics. A chain-level
   failure table would make the statement reproducible.

5. **Bootstrap interpretation.** The shared-seed block bootstrap interval is
   appropriate for the paired pooled comparison, but it is not a posterior
   interval and does not repair failed MCMC. State this next to the interval,
   not only in Web Appendix H.

6. **Design labels.** Use condition names in the joint study (Compatible,
   One-variable, Two-variable, Heterogeneous, Null) and keep Sc1--Sc5 labels
   for the separate original simulation. This prevents the 40-per-condition
   joint experiment from being mistaken for the older 100- or 500-replicate
   record.

7. **Interaction claim.** The checkerboard condition is a useful stress test,
   but neither marginal split signal nor a low discrepancy contrast establishes
   interaction recovery. Report the generated contrast, estimated contrast,
   and diagnostics together whenever the condition is discussed.

8. **Practical workflow.** Give one reproducible command path or archive
   manifest for the joint fit, including code/library commit, seeds, calibration
   inputs, chain lengths, and the saved-trace summary. The current report has
   these ingredients across several files but a practitioner should not have
   to reconstruct them from the validation directory.

# Suggestions to enhance the paper

- Add a one-panel graphical model showing the three roles of $f$, $g$, and
  $\psi$ and the distinction between RCT source membership and treatment.
- Frame the contribution as a joint treatment/control extension with a
  locally regularized control discrepancy, then state that the main evidence
  is a bounded randomized-trial validation rather than a superiority study.
- Retain the application because it is practically relevant, but use it to
  demonstrate profile-standardized external-control prediction and prior ESS,
  not to support the joint-treatment claim.
- Put the full per-condition joint table and convergence flags before the
  legacy separate-treatment operating-characteristic tables. This would make
  the scientific hierarchy visible without discarding the earlier record.

# Source-dependent flags

- Claims about novelty relative to MAP, commensurate, power-prior, BART source
  indicators, BCF, CAHB, and related local-borrowing methods require checking
  the cited papers directly. This review does not treat the manuscript’s
  “closest” labels as independently verified novelty findings.
- The deidentified merged n230 RData is locally available at
  `../yunxuan-repo/lrcbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n230.RData`.
  The raw 797-patient UCMM source remains private, and the 58 posterior result
  RData files are absent from this checkout. Cohort counts and RMST provenance
  are therefore accepted as documented by the available merged file, generated
  JSON/manifests, and source scripts; no new fit was run.
- The joint validation is a Gaussian experiment. It provides no evidence that
  the joint treatment model works for the survival application or for the
  single-arm setting.

# Verdict

The local method and the frozen validation record support a potentially useful
paper, but the working manuscript currently presents the wrong primary model
and overconnects a separate-treatment application and legacy simulations to a
new joint-treatment claim. The needed changes are substantial but concrete:
make the joint estimand and likelihood primary, report fair full-formulation
comparisons, attach convergence qualifications to the condition-specific
results, and reclassify the single-arm application as separate-model evidence.

VERDICT: major revision  
MAJOR COUNT: 5

What bothered me most: the manuscript asks the reader to evaluate a joint
treatment model while its main methods, application, and most operating
characteristics still describe a separately fitted treatment model.
