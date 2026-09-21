# Model walkthrough

Status: revised and settled with the user on 2026-09-16. This version supersedes the initial model walkthrough saved on the same date.

Purpose: authoritative content reference for the model section of the updated manuscript. This approved revision separates the model specification from the derivations in [posterior computation](posterior_computation_walkthrough.md), clarifies the sum-of-trees and variance priors, and states the survival, overlap, and single-arm assumptions. Implementation points at the end remain unresolved; saving this text does not change analysis code or reported results.

## 1. Define the two control populations and their outcome surfaces

Use group labels $\mathrm{trt}$, $1$ and $2$ for randomized treatment, randomized control and RWD control, respectively. For group $a\in\{\mathrm{trt},1,2\}$, let $n_a$ denote the sample size, $Y_{a,i}$ the outcome and $x_{a,i}$ the covariate vector for participant $i=1,\ldots,n_a$, and $\sigma_a^2$ the residual variance.
For $a\in\{\mathrm{trt},1,2\}$, write $\mathbf{X}_a\coloneqq (x_{a,1},\ldots,x_{a,n_a})$ and $\mathbf{Y}_a\coloneqq (Y_{a,1},\ldots,Y_{a,n_a})$. The trial size is $N\coloneqq n_{\mathrm{trt}}+n_{1}$. Trial standardization averages over groups $\mathrm{trt}$ and $1$, with one contribution per participant. In a single-arm trial, $n_{1}=0$, so only group $\mathrm{trt}$ contributes. For Gaussian outcomes, conditional independence is assumed under
\[
\begin{aligned}
Y_{\mathrm{trt},i}\mid x_{\mathrm{trt},i}&\sim N\{f_{\mathrm{trt}}(x_{\mathrm{trt},i}),\sigma_{\mathrm{trt}}^2\},\\
Y_{1,i}\mid x_{1,i}&\sim N\{f(x_{1,i})+g(x_{1,i}),\sigma_1^2\},\\
Y_{2,i}\mid x_{2,i}&\sim N\{f(x_{2,i}),\sigma_2^2\}.
\end{aligned}
\]

The observations are conditionally independent given the model parameters.

Here:

- \(f(x)\) is the RWD-control conditional mean.
- \(f(x)+g(x)\) is the RCT-control conditional mean.
- \(g(x)\) is the difference between the two control means at the same covariate profile.
- \(\sigma_{2}^2\) and \(\sigma_{1}^2\) describe individual-outcome variation around the respective means.

Group membership is observed and is shown directly in the subscripts; it is distinct from the latent leaf indicator.

**Both control sources contribute to estimating \(f\) in the final joint model.** RWD outcomes directly inform \(f\), while RCT-control outcomes inform \(f+g\). RWD information therefore also influences estimation of \(g\) through the jointly estimated \(f\).

This final fit differs from the preliminary ESS calculation, whose outcome likelihood for \(f\) uses only \(\mathbf{Y}_{2}\).

## 2. Represent the outcome surface and discrepancy using separate sums of trees

Suppressing the fixed outcome-centering constant in the notation, write

\[
\boxed{
f(x)\coloneqq \sum_{h^f=1}^{H_f}\sum_{\ell\in\mathcal T_{h^f}^f}\theta^f_{h^f\ell}\mathbf 1\{x\in\ell\},
\qquad
g(x)\coloneqq \sum_{h^g=1}^{H_g}\sum_{\ell\in\mathcal T_{h^g}^g}\theta^g_{h^g\ell}\mathbf 1\{x\in\ell\}.
}
\]

Here \(\mathcal T_{h^f}^f\) and \(\mathcal T_{h^g}^g\) denote prognostic and discrepancy tree structures, and \(\theta^f_{h^f\ell}\) and \(\theta^g_{h^g\ell}\) their terminal-node parameters. In each sum, \(\ell\) indexes a terminal-node region. Each node contributes its parameter when the profile lies in that region. The centering constant is restored when forming outcome predictions.

The two sums of trees have separate tree structures. Their priors favor smaller trees through depth-dependent splitting probabilities:

\[
P(\text{split at depth }d)
=
\alpha_f(1+d)^{-\beta_f}
\quad\text{or}\quad
\alpha_g(1+d)^{-\beta_g}.
\]

The current settings are

\[
(\alpha_f,\beta_f)=(0.95,2),
\qquad
(\alpha_g,\beta_g)=(0.5,3).
\]

Thus, the discrepancy sum-of-trees has a stronger prior preference for shallow trees, allowing the source difference to be simpler than the underlying outcome surface.

Candidate splitting variables and cutpoints can be shared while the realized trees remain distinct. The agreed ESS revision calls for a consistent candidate-cutpoint rule based on \(\mathbf{X}_{2}\) and \(\mathbf{X}_{\mathrm{trt}},\mathbf{X}_{1}\); implementation alignment is still pending.

## 3. Specify the leaf priors and explain local borrowing

The \(f\)-leaf parameters have normal priors:

\[
\theta^f_{h^f\ell}\mid \mathcal T_{h^f}^f
\sim N(0,\lambda_f^2),
\]

independently conditional on the tree structures and fixed prior scale.

For each \(g\) leaf,

\[
z_{h^g\ell}\mid w
\sim\operatorname{Bernoulli}(w),
\]

and

\[
\boxed{
\theta^g_{h^g\ell}\mid z_{h^g\ell},\tau_0^2,\tau_1^2
\sim
\begin{cases}
N(0,\tau_0^2),&z_{h^g\ell}=1,\\
N(0,\tau_1^2),&z_{h^g\ell}=0,
\end{cases}
\qquad
0<\tau_0^2<\tau_1^2.
}
\]

The indicators are independent conditional on \(w\), and the leaf parameters are independent conditional on their indicators, scales, and tree structures.

The narrow spike shrinks a leaf contribution strongly toward zero. The broader slab allows larger positive or negative source differences. Both components remain centered at zero.

Consequently, the prior encourages the RCT-control mean to remain near \(f(x)\), while allowing departures supported by the trial controls.

Several distinctions matter:

- The spike is continuous, so a spike leaf can have a nonzero contribution.
- \(z_{h^g\ell}\) is specific to a leaf; the group subscripts identify treatment and control sources.
- \(w\) is a shared prior allocation probability, not a patient-specific borrowing weight.
- At a profile \(x\), the discrepancy \(g(x)\) combines one leaf contribution from every discrepancy tree.

The posterior-computation section explains how the data update the trees, indicators, and leaf values.

## 4. Specify the remaining priors and the role of ESS calibration

For the default two-arm model,

\[
w\sim\operatorname{Beta}(1,1),
\]

and

\[
\boxed{
\tau_0^2
\sim
\text{scaled-Inv-}\chi^2(\nu_0,s_0^2)
\text{ truncated to }(0,\tau_1^2),
\qquad \nu_0=3.
}
\]

The quantities have distinct roles:

| Quantity | Role |
| --- | --- |
| \(\lambda_f^2\) | Prior variance of an \(f\)-leaf contribution |
| \(\tau_0^2\) | Spike variance for a \(g\) leaf |
| \(\tau_1^2\) | Fixed, broader slab variance |
| \(w\) | Common prior probability of a spike allocation |
| \(s_0^2\) | Scale parameter governing the prior for \(\tau_0^2\) |

The implementation sets \(\lambda_f^2\) and \(\tau_1^2\) using the outcome range and the respective numbers of trees.

For the two control-source residual variances, the implemented prior can be written

\[
\sigma_s^2
\sim
\operatorname{InvGamma}
\left(
\frac{\nu_\sigma}{2},
\frac{\nu_\sigma\lambda_{\sigma,s}}{2}
\right),
\qquad s=2,3,
\]

independently, using the shape–scale convention. The current implementation uses \(\nu_\sigma=3\), with positive variance-scale settings \(\lambda_{\sigma,s}\) established before sampling.

**ESS calibration selects \(s_0^2\).** In the default two-arm fit, the data subsequently update \(\tau_0^2\), \(w\), the leaf allocations, and the outcome surfaces. The ESS target therefore specifies a prior calibration rather than a fixed amount of posterior borrowing.

The separate ESS section retains the distinction between its untruncated calibration scale law and the truncated scale prior used for fitting.

## 5. Define the treatment contrast over the trial population

Fit the treated outcomes separately:

\[
Y_{\mathrm{trt},i}\mid x_{\mathrm{trt},i}
\sim N\!\left(f_{\mathrm{trt}}(x_{\mathrm{trt},i}),\sigma_{\mathrm{trt}}^2\right),
\]

using a treatment-arm BART model for $f_{\mathrm{trt}}$.

For Gaussian outcomes, the trial-standardized mean contrast is

\[
\boxed{
\Delta
\coloneqq 
\frac1N\sum_{a\in\{\mathrm{trt},1\}}\sum_{i=1}^{n_a}
\left[
f_{\mathrm{trt}}(x_{a,i})-f(x_{a,i})-g(x_{a,i})
\right].
}
\]

The average uses all trial covariate profiles, including treated subjects and controls when both are present. Evaluating this expression at each posterior draw gives the posterior distribution of the contrast.

The ability to estimate the source difference from outcomes depends on covariate support. Where both control sources are represented, their data inform the comparison between \(f(x)\) and \(f(x)+g(x)\). Outside their shared support, predictions depend more heavily on the sum-of-trees structure, priors, and extrapolation.

## 6. Extend the model to survival outcomes

For survival outcomes, the Gaussian model applies to latent log event time:

\[
\begin{aligned}
\log T_{2,i}&=f(x_{2,i})+\epsilon_{2,i},\\
\log T_{1,i}&=f(x_{1,i})+g(x_{1,i})+\epsilon_{1,i},\quad\text{randomized controls}.
\end{aligned}
\]

Accordingly, \(f\) and \(f+g\) describe conditional means of **log event time**, and \(g\) describes the difference between those log-time means.

Under conditional noninformative censoring, an uncensored subject contributes its event-time density, while a subject censored at \(C_i\) contributes the probability of surviving beyond \(C_i\). The posterior-computation section handles censored log times through data augmentation.

For an RCT-control profile,

\[
S_{1}(t\mid x)
=
1-\Phi\!\left(
\frac{\log t-f(x)-g(x)}{\sigma_{1}}
\right).
\]

The treatment survival function uses \(f_{\mathrm{trt}}(x)\) and \(\sigma_{\mathrm{trt}}\). Standardize each arm to the same trial profiles:

\[
\bar S_a(t)
\coloneqq 
\frac1N\sum_{s\in\{\mathrm{trt},1\}}\sum_{i=1}^{n_s}S_a(t\mid x_{s,i}),
\qquad a\in\{\mathrm{trt},1\}.
\]

The population median survival ratio compares the medians obtained from \(\bar S_{\mathrm{trt}}\) and \(\bar S_{1}\). An RMST ratio compares their integrals up to a prespecified horizon. Both are functionals of standardized survival curves, rather than simply transformations of the average log-time contrast.

## 7. State the additional assumptions for a single-arm trial

When \(n_{1}=0\), there are no observed RCT-control outcomes.

The RWD outcomes estimate \(f\), and the trial-control prediction remains

\[
f(x)+g(x).
\]

Conditional on the specified hyperparameters, \(g\) receives no outcome-likelihood update and remains governed by its prior. The separately fitted treated outcomes estimate \(f_{\mathrm{trt}}\); they do not identify the untreated trial discrepancy.

The current primary single-arm configuration fixes \(w=1\), with \(w<1\) used for sensitivity analysis. Even with \(w=1\), the discrepancy remains uncertain because its spike variance is positive.

The single-arm implementation also assumes

\[
\boxed{\sigma_{1}^2=\sigma_{2}^2}
\]

for the unobserved trial-control residual variance. This assumption matters for ESS reference units and counterfactual survival predictions, and cannot be checked using RCT-control outcomes in a single-arm trial.

The resulting treatment comparison depends on the specified discrepancy and variance assumptions, alongside the information supplied by the RWD.

## Implementation points retained with this review

- **Single-arm spike scale:** both single-arm scripts currently fix \(\tau_0^2\) at \(\min(s_0^2,\tau_1^2/2)\). They therefore do not propagate a spike-scale hyperprior in the same way as the two-arm model. This remains unresolved. [Single-arm setting](/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-gaussian-single-arm/lrcBART.R:335), [sampler initialization](/Users/oliviazhang/Desktop/lrcBART/clrcbart.cpp:121).
- **Tree support:** the shared candidate-cutpoint and tree-support specification agreed during the ESS discussion has not yet been aligned across all stages.
- **Additional survival RMST output:** the two-arm survival script currently uses an adaptive control-variance choice that can substitute treatment variance. Those outputs cannot automatically be described as RMST calculations from the control survival model above. This needs resolution when settling the reported analyses. [RMST setting](/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-survival/lrcBART.R:18), [variance selection](/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-sim-survival/rmst_helpers.R:132).
