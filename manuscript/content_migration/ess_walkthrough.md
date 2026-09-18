# Prior effective sample size for the RWD-informed trial-control mean

## Notation and objective

Let groups \(1\), \(2\), and \(3\) denote randomized-trial treatment, randomized-trial control, and historical control. The prior-ESS target is defined on the \(n_2\) randomized-control covariate profiles \(\mathbf X_2=(x_{2,1},\ldots,x_{2,n_2})\). Historical outcomes and covariates \((\mathbf Y_3,\mathbf X_3)\) inform the prognostic surface \(f\), whereas the discrepancy surface \(g\) describes the current-minus-historical control difference.

The goal is to express the prior information about the standardized current-control mean as an equivalent number of observations from a prespecified reference likelihood. The construction:

1. identifies the induced target prior and its conditional marginal variance;
2. derives the variance contributions from \(f\) and \(g\);
3. matches inverse target variance to Fisher information from reference observations; and
4. averages the conditional equivalent counts over the calibration distribution of the spike variance.

## 1. Target prior and variance decomposition

Define the historical-control, discrepancy, and current-control means over \(\mathbf X_2\) by

\[
\mu_f\coloneqq \frac1{n_2}\sum_{i=1}^{n_2}f(x_{2,i}),
\qquad
\mu_g\coloneqq \frac1{n_2}\sum_{i=1}^{n_2}g(x_{2,i}),
\qquad
\mu\coloneqq\mu_f+\mu_g.
\]

A preliminary historical fit produces the posterior of \(f\) from \((\mathbf Y_3,\mathbf X_3)\), evaluated at \(\mathbf X_2\). Independently, the specified prior for \(g\) induces a distribution for \(\mu_g\). Together they induce the target prior

\[
\pi_{\mathrm{target}}
(\mu\mid\mathbf Y_3,\mathbf X_3,\mathbf X_2,w,\tau_0^2).
\]

This notation names the induced prior; it does not introduce a separate modeling distribution. Conditional on \(w\), \(\tau_0^2\), the covariates, and all other fixed prior settings, the historical posterior for \(f\) and the prior for \(g\) are independent. Hence

\[
\begin{aligned}
V_{\mathrm{target}}(w,\tau_0^2)
&\coloneqq
\operatorname{Var}
(\mu\mid\mathbf Y_3,\mathbf X_3,\mathbf X_2,w,\tau_0^2)\\
&=V_f+V_g(w,\tau_0^2),
\end{aligned}
\]

where

\[
V_f\coloneqq
\operatorname{Var}(\mu_f\mid\mathbf Y_3,\mathbf X_3,\mathbf X_2)
\]

and

\[
V_g(w,\tau_0^2)\coloneqq
\operatorname{Var}(\mu_g\mid\mathbf X_2,w,\tau_0^2).
\]

\(V_f\) is the uncertainty remaining in the standardized historical surface after observing the historical data. \(V_g\) is the uncertainty introduced when transporting that surface to the current-control population.

## 2. Calculating \(V_f\)

For posterior draw \(b=1,\ldots,B\) from the preliminary historical fit, calculate

\[
\mu_f^{(b)}\coloneqq
\frac1{n_2}\sum_{i=1}^{n_2}f^{(b)}(x_{2,i}).
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
=\frac1{n_2^2}
\sum_{i=1}^{n_2}\sum_{i'=1}^{n_2}
\operatorname{Cov}
\{f(x_{2,i}),f(x_{2,i'})
\mid\mathbf Y_3,\mathbf X_3,\mathbf X_2\}.
\]

The variance of posterior averages retains tree uncertainty and covariance among predictions at different target profiles. Averaging pointwise variances would omit the covariance terms.

## 3. Calculating \(V_g(w,\tau_0^2)\)

Write the discrepancy surface as

\[
g(x)\coloneqq
\sum_{h^g=1}^{H_g}
\sum_{\ell\in\mathcal T_{h^g}^g}
\theta^g_{h^g\ell}\mathbf1\{x\in\ell\}.
\]

For terminal node \(\ell\) of discrepancy tree \(h^g\), define its target-profile weight

\[
a_{h^g\ell}\coloneqq
\frac1{n_2}\sum_{i=1}^{n_2}
\mathbf1\{x_{2,i}\in\ell\text{ of tree }h^g\}.
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

## 4. Target and reference information

The proposed method measures target-prior information by inverse conditional marginal variance:

\[
\mathcal I_{\mathrm{target}}(w,\tau_0^2)
\coloneqq
\frac1{V_{\mathrm{target}}(w,\tau_0^2)}
=\frac1{V_f+V_g(w,\tau_0^2)}.
\]

This is a global concentration measure. If the target prior is normal, it also equals its constant negative log-prior curvature; exact normality is not required for the variance-based definition.

Use one normal current-control observation as the reference information unit:

\[
Y\mid\mu\sim N(\mu,\sigma_2^2),
\]

where \(\sigma_2^2\) is fixed for the theoretical calibration. Its expected Fisher information is

\[
\mathcal I_{\mathrm{ref}}(\mu)
\coloneqq i_F(\mu)
\coloneqq
E_{Y\mid\mu}
\left[-\frac{\partial^2}{\partial\mu^2}
\log p(Y\mid\mu)\right]
=\frac1{\sigma_2^2}.
\]

No reference prior is required. Independence makes the reference information additive, so \(m\) observations supply \(m/\sigma_2^2\).

## 5. Conditional matching and prior ESS

For fixed \(w\) and \(\tau_0^2\), define \(m(w,\tau_0^2)\) by writing the target information on the left and reference information on the right:

\[
\mathcal I_{\mathrm{target}}(w,\tau_0^2)
=m(w,\tau_0^2)\mathcal I_{\mathrm{ref}}(\mu).
\]

Thus

\[
\boxed{
m(w,\tau_0^2)
=\frac{\sigma_2^2}{V_f+V_g(w,\tau_0^2)}.
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
\frac{\sigma_2^2}{V_f+V_g(w,\tau_0^2)}
\right].
}
\]

The result can be noninteger because it is an information equivalent. It is conditional on the selected target profiles, reference likelihood, historical dataset, covariate-dependent tree prior, and fixed calibration settings.

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
\frac{\sigma_2^2}{
V_f+\tau_0^2E_{\mathcal T^g}
(\sum_{h^g,\ell}a_{h^g\ell}^2)}
\right].
\]

Under \(0<V_f<\infty\), \(0<\sigma_2^2<\infty\), \(\nu_0>0\), and a finite discrepancy tree almost surely, this curve is continuous, strictly decreasing, and strictly convex in \(s_0^2>0\). Its limits are

\[
\lim_{s_0^2\downarrow0}
\operatorname{ESS}_{\tau_0}(s_0^2)
=\frac{\sigma_2^2}{V_f},
\qquad
\lim_{s_0^2\to\infty}
\operatorname{ESS}_{\tau_0}(s_0^2)=0.
\]

Hence

\[
\operatorname{ESS}_{\mathrm{ceiling}}
\coloneqq\frac{\sigma_2^2}{V_f}.
\]

As \(s_0^2\downarrow0\), discrepancy variation vanishes but historical-posterior uncertainty remains. The ceiling therefore depends on the information supplied by the historical data through \(V_f\); it is not a historical head count.

## 7. Relationship to established prior-ESS definitions

All four definitions compare a target-prior information quantity with a calibrated or reference information quantity, but they use different information measures and calibrations.

| Method | Target information | Reference or calibration | ESS interpretation |
| --- | --- | --- | --- |
| Neuenschwander et al. (2010) | Inverse marginal variance \(1/V_{\mathrm{target}}\) of a partial-pooling MAP prior | Complete-pooling MAP prior with variance \(V_0\) and assigned ESS \(N_H\) | \(N_HV_0/V_{\mathrm{target}}\): information retained relative to complete pooling |
| Morita–Thall–Müller (2008) | Log-prior curvature at the target-prior mean | Expected log-posterior curvature from a same-mean vague prior plus \(m\) hypothetical observations | The \(m\) giving the closest curvature match |
| ELIR | Pointwise log-prior curvature | Expected Fisher information in one current-study observation | Prior expectation of the local information ratio |
| Proposed method | Inverse conditional marginal variance \(1/\{V_f+V_g(w,\tau_0^2)\}\) | Expected Fisher information \(1/\sigma_2^2\) in one normal current-control observation | Mean, over \(\tau_0^2\), of the conditional variance-matched counts |

Neuenschwander’s complete-pooling row is an ESS anchor rather than a hypothetical-likelihood reference. Morita’s reference object is a posterior formed from a low-information prior and hypothetical data. ELIR and the proposed method use one-observation Fisher information directly, but ELIR averages local curvature ratios over the target prior, whereas the proposed method averages conditional inverse-variance matches over the spike-variance calibration distribution. The proposed ESS therefore does not inherit ELIR’s predictive-consistency property.

## 8. Implementation and interpretation

Numerical calibration substitutes a fixed estimate \(\widehat\sigma_2^2\) for the theoretical \(\sigma_2^2\). The theoretical curve estimates \(V_f\) from all preliminary posterior draws, simulates discrepancy trees from the calibration prior, and uses common chi-square draws across the scale grid. The reported selector additionally constructs block-specific estimates of \(V_f\) and applies its stated 0.95 acceptance rule. That selection procedure is an implementation rule rather than a change to the ESS definition.

The calibration distribution for \(\tau_0^2\), the calibration tree support, and the final fitting prior need not coincide. In particular, calibration uses an untruncated spike-variance distribution, while the two-arm fit uses a truncated distribution and single-arm fits may fix \(\tau_0^2\). Reported ESS must therefore be labeled as the calibration definition rather than an exact post-fit information measure.

An ESS of \(m\) means that the full historical-data-informed target prior has, on average over the spike-variance calibration distribution, the same conditional precision about the specified current-control mean as \(m\) independent observations from the normal reference likelihood. It is not the number of historical patients retained, a compatibility diagnostic, or an additive patient count. Overall and subgroup ESS values are generally nonadditive because their standardized means are correlated and the inverse-variance transformation is nonlinear.

The group-2 target requires prespecified current-control profiles. When \(n_2=0\), the same definition can be used only after specifying a separate target profile set; it cannot be formed from an empty \(\mathbf X_2\).

## References

1. Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ. *Summarizing historical information on controls in clinical trials.* Clinical Trials. 2010;7(1):5–18.
2. Morita S, Thall PF, Müller P. *Determining the effective sample size of a parametric prior.* Biometrics. 2008;64(2):595–602.
3. Neuenschwander B, Weber S, Schmidli H, O’Hagan A. *Predictively consistent prior effective sample sizes.* Biometrics. 2020;76(2):578–587.
