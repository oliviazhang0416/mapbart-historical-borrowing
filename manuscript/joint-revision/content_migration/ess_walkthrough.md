# Prior effective sample size for the RWD-informed prior on the current-control mean

## Notation and objective

Let groups \(\mathrm{trt}\), \(1\), and \(2\) denote randomized treatment, randomized control, and RWD control, respectively. The current-control mean is defined on the \(n_{1}\) randomized-control covariate profiles \(\mathbf X_{1}=(x_{1,1},\ldots,x_{1,n_{1}})\). Historical outcomes and covariates \((\mathbf Y_{2},\mathbf X_{2})\) inform the prognostic surface \(f\), whereas the discrepancy surface \(g\) describes the current-minus-historical control difference.

The goal is to express the information in the prior informed by historical real-world data (the RWD-informed prior) for the current-control mean as an equivalent number of observations from a prespecified reference likelihood. The construction:

1. identifies the induced RWD-informed prior and its conditional marginal variance;
2. derives the variance contributions from \(f\) and \(g\);
3. matches inverse RWD-informed-prior variance to Fisher information from reference observations; and
4. averages the conditional equivalent counts over the calibration distribution of the spike variance.

## 1. RWD-informed prior and variance decomposition

Define the historical-control, discrepancy, and current-control means over \(\mathbf X_{1}\) by

\[
\mu_f\coloneqq \frac1{n_{1}}\sum_{i=1}^{n_{1}}f(x_{1,i}),
\qquad
\mu_g\coloneqq \frac1{n_{1}}\sum_{i=1}^{n_{1}}g(x_{1,i}),
\qquad
\mu\coloneqq\mu_f+\mu_g.
\]

A preliminary historical fit produces the posterior of \(f\) from \((\mathbf Y_{2},\mathbf X_{2})\), evaluated at \(\mathbf X_{1}\). Independently, the specified prior for \(g\) induces a distribution for \(\mu_g\). Together they induce the RWD-informed prior

\[
\pi
(\mu\mid\mathbf Y_{2},\mathbf X_{2},\mathbf X_{1},w,\tau_0^2).
\]

Conditional on \(w\), \(\tau_0^2\), the covariates, and all other fixed prior settings, the historical posterior for \(f\) and the prior for \(g\) are independent. Hence

\[
\begin{aligned}
V_{\pi}(w,\tau_0^2)
&\coloneqq
\operatorname{Var}_{\pi}
(\mu\mid\mathbf Y_{2},\mathbf X_{2},\mathbf X_{1},w,\tau_0^2)\\
&=V_f+V_g(w,\tau_0^2),
\end{aligned}
\]

where

\[
V_f\coloneqq
\operatorname{Var}(\mu_f\mid\mathbf Y_{2},\mathbf X_{2},\mathbf X_{1})
\]

and

\[
V_g(w,\tau_0^2)\coloneqq
\operatorname{Var}(\mu_g\mid\mathbf X_{1},w,\tau_0^2).
\]

\(V_f\) is the uncertainty remaining in the historical surface averaged over the randomized-control profiles after observing the historical data. \(V_g\) is the uncertainty introduced when moving from the historical-control mean to the current-control mean.

## 2. Calculating \(V_f\)

For posterior draw \(b=1,\ldots,B\) from the preliminary historical fit, calculate

\[
\mu_f^{(b)}\coloneqq
\frac1{n_{1}}\sum_{i=1}^{n_{1}}f^{(b)}(x_{1,i}).
\]

Then estimate \(V_f\) by

\[
\widehat V_f
=
\frac1{B-1}\sum_{b=1}^B
(\mu_f^{(b)}-\bar\mu_f)^2,
\qquad
\bar\mu_f\coloneqq B^{-1}\sum_{b=1}^B\mu_f^{(b)}.
\]

Equivalently,

\[
V_f
=\frac1{n_{1}^2}
\sum_{i=1}^{n_{1}}\sum_{i'=1}^{n_{1}}
\operatorname{Cov}
\{f(x_{1,i}),f(x_{1,i'})
\mid\mathbf Y_{2},\mathbf X_{2},\mathbf X_{1}\}.
\]

The variance of posterior averages retains tree uncertainty and covariance among predictions at different randomized-control profiles. Averaging pointwise variances would omit the covariance terms.

## 3. Calculating \(V_g(w,\tau_0^2)\)

Write the discrepancy surface as

\[
g(x)\coloneqq
\sum_{h^g=1}^{H_g}
\sum_{\ell\in\mathcal T_{h^g}^g}
\theta^g_{h^g\ell}\mathbf1\{x\in\ell\}.
\]

For terminal node \(\ell\) of discrepancy tree \(h^g\), define its randomized-control-profile weight

\[
a_{h^g\ell}\coloneqq
\frac1{n_{1}}\sum_{i=1}^{n_{1}}
\mathbf1\{x_{1,i}\in\ell\text{ of tree }h^g\}.
\]

Regrouping the profile averages gives

\[
\mu_g
=\sum_{h^g,\ell}a_{h^g\ell}\theta^g_{h^g\ell}.
\]

The leaf prior is

\[
z_{h^g\ell}\mid w\sim\operatorname{Bernoulli}(w),
\qquad
\theta^g_{h^g\ell}\mid z_{h^g\ell},\tau_0^2,\tau_1^2
\sim
\begin{cases}
N(0,\tau_0^2),&z_{h^g\ell}=1,\\
N(0,\tau_1^2),&z_{h^g\ell}=0.
\end{cases}
\]

Conditional on the trees and indicators, independence and the zero prior means give

\[
\operatorname{Var}(\mu_g\mid\mathcal T^g,\mathbf z,w,\tau_0^2)
=
\sum_{h^g,\ell}a_{h^g\ell}^2
\{z_{h^g\ell}\tau_0^2+(1-z_{h^g\ell})\tau_1^2\}.
\]

The conditional mean remains zero when indicators and trees are integrated out, so the between-mean terms in both applications of the law of total variance vanish. Therefore

\[
\boxed{
V_g(w,\tau_0^2)
=
\{w\tau_0^2+(1-w)\tau_1^2\}
E_{\mathcal T^g}
\left(\sum_{h^g,\ell}a_{h^g\ell}^2\right).
}
\]

The expectation is with respect to the discrepancy-tree prior used for calibration. It is evaluated by prior-tree simulation.

## 4. Information in the RWD-informed prior and reference likelihood

The proposed method measures information in the RWD-informed prior by inverse conditional marginal variance:

\[
\mathcal I_{\pi}(w,\tau_0^2)
\coloneqq
\frac1{V_{\pi}(w,\tau_0^2)}
=\frac1{V_f+V_g(w,\tau_0^2)}.
\]

This quantity measures the concentration of the RWD-informed prior through its marginal variance.

Use one normal current-control observation as the reference information unit:

\[
Y\mid\mu\sim N(\mu,\sigma_{1}^2),
\]

where \(\sigma_{1}^2\) is fixed for the theoretical calibration. Its expected Fisher information is

\[
\mathcal I_{\mathrm{ref}}(\mu)
\coloneqq i_F(\mu)
\coloneqq
E_{Y\mid\mu}
\left[-\frac{\partial^2}{\partial\mu^2}
\log p(Y\mid\mu)\right]
=\frac1{\sigma_{1}^2}.
\]

Fisher information is additive for independent observations, so \(m\) reference observations supply \(m/\sigma_{1}^2\).

## 5. Conditional matching and prior ESS

For fixed \(w\) and \(\tau_0^2\), define \(m(w,\tau_0^2)\) by matching the information in the RWD-informed prior to the information in \(m(w,\tau_0^2)\) reference observations:

\[
\mathcal I_{\pi}(w,\tau_0^2)
=m(w,\tau_0^2)\mathcal I_{\mathrm{ref}}(\mu).
\]

Thus

\[
\boxed{
m(w,\tau_0^2)
=\frac{\sigma_{1}^2}{V_f+V_g(w,\tau_0^2)}.
}
\]

The scale \(s_0^2\) indexes the calibration distribution of \(\tau_0^2\). The proposed prior ESS averages the conditional equivalent counts:

\[
\boxed{
\operatorname{ESS}(w,s_0^2)
\coloneqq
E_{\tau_0^2}\{m(w,\tau_0^2)\}
=E_{\tau_0^2}
\left[
\frac{\sigma_{1}^2}{V_f+V_g(w,\tau_0^2)}
\right].
}
\]

The result can be noninteger because it is an information equivalent. It is conditional on the selected randomized-control profiles, reference likelihood, historical dataset, covariate-dependent tree prior, and fixed calibration settings.

## 6. All-spike calibration and information ceiling

For prior-scale selection, set \(w=1\) and use the untruncated distribution

\[
\tau_0^2\sim\operatorname{scaled\text{-}Inv}\chi^2(\nu_0,s_0^2).
\]

Then

\[
V_g(1,\tau_0^2)
=\tau_0^2
E_{\mathcal T^g}
\left(\sum_{h^g,\ell}a_{h^g\ell}^2\right)
\]

and

\[
\operatorname{ESS}_{\tau_0}(s_0^2)
\coloneqq\operatorname{ESS}(1,s_0^2)
=E_{\tau_0^2}
\left[
\frac{\sigma_{1}^2}{
V_f+\tau_0^2E_{\mathcal T^g}
(\sum_{h^g,\ell}a_{h^g\ell}^2)}
\right].
\]

Under \(0<V_f<\infty\), \(0<\sigma_{1}^2<\infty\), \(\nu_0>0\), and a finite discrepancy tree almost surely, this curve is continuous, strictly decreasing, and strictly convex in \(s_0^2>0\). Its limits are

\[
\lim_{s_0^2\downarrow0}
\operatorname{ESS}_{\tau_0}(s_0^2)
=\frac{\sigma_{1}^2}{V_f},
\qquad
\lim_{s_0^2\to\infty}
\operatorname{ESS}_{\tau_0}(s_0^2)=0.
\]

Hence

\[
\operatorname{ESS}_{\mathrm{ceiling}}
\coloneqq\frac{\sigma_{1}^2}{V_f}.
\]

As \(s_0^2\downarrow0\), discrepancy variation vanishes but historical-posterior uncertainty remains. The ceiling therefore depends on the information supplied by the historical data through \(V_f\); it is not a historical head count.

## 7. Relationship to established prior-ESS definitions

All four definitions compare a target-prior information quantity with a calibrated or reference information quantity, but they use different information measures and calibrations.

| Method | Target information | Reference or calibration | ESS interpretation |
| --- | --- | --- | --- |
| Neuenschwander et al. (2010) | Inverse marginal variance \(1/V_{\mathrm{target}}\) of a partial-pooling MAP prior | Complete-pooling MAP prior with variance \(V_0\) and assigned ESS \(N_H\) | \(N_HV_0/V_{\mathrm{target}}\): information retained relative to complete pooling |
| Morita–Thall–Müller (2008) | Log-prior curvature at the target-prior mean | Expected log-posterior curvature from a same-mean vague prior plus \(m\) hypothetical observations | The \(m\) giving the closest curvature match |
| ELIR | Pointwise log-prior curvature | Expected Fisher information in one current-study observation | Prior expectation of the local information ratio |
| Proposed method | Inverse conditional marginal variance \(1/\{V_f+V_g(w,\tau_0^2)\}\) | Expected Fisher information \(1/\sigma_{1}^2\) in one normal current-control observation | Mean, over \(\tau_0^2\), of the conditional variance-matched counts |

Neuenschwander’s complete-pooling row supplies the ESS anchor for its partial-pooling prior. Morita matches target-prior curvature to the expected curvature of a posterior formed from a low-information prior and hypothetical data. ELIR averages local ratios of prior curvature to one-observation Fisher information over the target prior. The proposed method averages conditional inverse-variance matches over the spike-variance calibration distribution.

## 8. Implementation and interpretation

Numerical calibration substitutes a fixed estimate \(\widehat\sigma_{1}^2\) for the theoretical \(\sigma_{1}^2\). The theoretical curve estimates \(V_f\) from all preliminary posterior draws, simulates discrepancy trees from the calibration prior, and uses common chi-square draws across the scale grid. The reported selector additionally constructs block-specific estimates of \(V_f\) and applies its stated 0.95 acceptance rule. That selection procedure is an implementation rule rather than a change to the ESS definition.

The calibration uses an untruncated spike-variance distribution and its specified tree support. The two-arm fit uses a truncated spike-variance distribution, and single-arm fits may fix \(\tau_0^2\). The reported prior ESS therefore corresponds to the calibration specification.

An ESS of \(m\) means that the RWD-informed prior has, on average over the spike-variance calibration distribution, the same conditional precision about the specified current-control mean as \(m\) independent observations from the normal reference likelihood. Pointwise ESS values do not generally sum or average to the overall ESS because predictions at different profiles are correlated and the inverse-variance transformation is nonlinear.

For a single-arm analysis, the current-control standardization profiles are specified separately.

## References

1. Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ. *Summarizing historical information on controls in clinical trials.* Clinical Trials. 2010;7(1):5–18.
2. Morita S, Thall PF, Müller P. *Determining the effective sample size of a parametric prior.* Biometrics. 2008;64(2):595–602.
3. Neuenschwander B, Weber S, Schmidli H, O’Hagan A. *Predictively consistent prior effective sample sizes.* Biometrics. 2020;76(2):578–587.
