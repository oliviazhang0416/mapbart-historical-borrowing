# ESS walkthrough

Notation: \(1\) denotes trial treatment, \(2\) trial control and \(3\) historical control. Standardization uses \(\mathbf{X}_1\) and \(\mathbf{X}_2\), with \(N\coloneqq n_1+n_2\). Each participant contributes once. In a single-arm trial, \(n_2=0\), so standardization uses only \(\mathbf{X}_1\).

Status: settled with the user on 2026-09-16, including the ESS ceiling definition and limiting derivation.

Purpose: authoritative content reference for the future ESS sections of the main manuscript and appendix. The walkthrough below preserves the agreed explanation. Its proposed covariate/cutpoint revision is not an assertion that the analysis code has already been changed.

---

The aim is to find a reference sample size \(m\) such that **the posterior variance of the unknown mean under a flat-prior reference analysis equals the RWD-informed prior variance of the trial control mean**. We then interpret that matched sample size as a conditional prior ESS.

The distribution obtained from historical outcomes is a posterior relative to those outcomes. Combined with the discrepancy prior, it supplies prior information for the trial control mean.

1. **Derive the reference posterior variance.**

   Consider \(m\) hypothetical independent observations,

   \[
   Y_j^{\mathrm{ref}}\mid\mu
   \overset{\mathrm{iid}}{\sim}
   N(\mu,\sigma_2^2),
   \qquad j=1,\ldots,m,
   \]

   with known \(\sigma_2^2\) and a flat baseline prior \(p(\mu)\propto1\).

   By Bayes' rule,

   \[
   p(\mu\mid Y_{1:m}^{\mathrm{ref}})
   \propto
   \exp\left\{
   -\frac{1}{2\sigma_2^2}
   \sum_{j=1}^{m}(Y_j^{\mathrm{ref}}-\mu)^2
   \right\}.
   \]

   Since

   \[
   \sum_j(Y_j^{\mathrm{ref}}-\mu)^2
   =
   \sum_j(Y_j^{\mathrm{ref}}-\bar Y^{\mathrm{ref}})^2
   +
   m(\mu-\bar Y^{\mathrm{ref}})^2,
   \]

   dropping the term that does not depend on \(\mu\) gives

   \[
   p(\mu\mid Y_{1:m}^{\mathrm{ref}})
   \propto
   \exp\left\{
   -\frac{(\mu-\bar Y^{\mathrm{ref}})^2}
   {2\sigma_2^2/m}
   \right\}.
   \]

   Therefore,

   \[
   \boxed{
   \mu\mid Y_{1:m}^{\mathrm{ref}}
   \sim N\!\left(\bar Y^{\mathrm{ref}},\frac{\sigma_2^2}{m}\right).
   }
   \]

   Thus, \(\sigma_2^2/m\) is the **posterior variance of the unknown mean**, rather than the variance of an individual outcome.

2. **Construct the RWD-informed distribution for the trial control mean.**

   Let \(\mathbf{X}_3\) and \(\mathbf{Y}_3\) denote the historical covariates and outcomes, and \(\mathbf{X}_1,\mathbf{X}_2\) the trial covariates.

   The trial control surface is \(f(x)+g(x)\). Over the \(N\) trial profiles,

   \[
   \mu
   \coloneqq 
   \frac1N\sum_{s=1}^2\sum_{i=1}^{n_s}\{f(x_{s,i})+g(x_{s,i})\}
   =
   \mu_f+\mu_g.
   \]

   In the proposed revision, use \(\mathbf{X}_3\) and \(\mathbf{X}_1,\mathbf{X}_2\) to construct candidate cutpoints consistently across the preliminary \(f\) fit, the prior-tree calculation for \(g\), and the final joint fit. This does not require \(f\) and \(g\) to share realized trees.

   The preliminary fit uses only \(\mathbf{Y}_3\) in its outcome likelihood:

   \[
   p(f\mid \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2)
   \propto
   p(\mathbf{Y}_3\mid f,\mathbf{X}_3)
   \,
   p(f\mid \mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2).
   \]

   Combine this distribution with an independent discrepancy prior for \(g\). The covariates, slab variance \(\tau_1^2\), and other specified reference settings are held fixed throughout.

3. **Solve the reference sample-size matching problem.**

   Initially hold \(w\) and \(\tau_0^2\) fixed. Find \(m\) satisfying

   \[
   \underbrace{\frac{\sigma_2^2}{m}}_{
   \text{reference posterior variance}
   }
   =
   \underbrace{
   \operatorname{Var}
   \left(
   \mu\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,
   w,\tau_0^2
   \right)
   }_{
   \text{RWD-informed prior variance for the trial mean}
   }.
   \]

   Solving gives

   \[
   \boxed{
   m(w,\tau_0^2)
   =
   \frac{\sigma_2^2}
   {\operatorname{Var}
   \left(
   \mu\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,
   w,\tau_0^2
   \right)}.
   }
   \]

4. **Interpret the matched sample size as conditional prior ESS.**

   Suppose the solution is \(m=20\). The RWD-informed distribution for \(\mu\) then has the same variance as the reference posterior obtained from 20 independent normal observations under a flat baseline prior.

   Under this variance-matching definition, we therefore assign the RWD-informed distribution a **prior ESS of 20 reference observations**.

   More generally, \(m(w,\tau_0^2)\) is its **conditional prior ESS**, conditional on \(w\) and \(\tau_0^2\). It may be noninteger because it represents an information equivalent, rather than an actual patient count.

5. **Define how conditional ESS values are summarized over the spike variance.**

   For a candidate \(s_0^2\), specify

   \[
   \tau_0^2
   \sim
   \text{scaled-Inv-}\chi^2(\nu_0,s_0^2).
   \]

   Different values of \(\tau_0^2\) produce different matched sample sizes \(m(w,\tau_0^2)\). The calibration chooses their mean:

   \[
   \boxed{
   \mathrm{ESS}(w,s_0^2)
   \coloneqq 
   E_{\tau_0^2}[m(w,\tau_0^2)].
   }
   \]

   Equivalently, the unexpanded definition is

   \[
   \mathrm{ESS}(w,s_0^2)
   =
   E_{\tau_0^2}
   \left[
   \frac{\sigma_2^2}
   {\operatorname{Var}
   \left(
   \mu\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,
   w,\tau_0^2
   \right)}
   \right].
   \]

   **Averaging these conditional matches is an additional calibration choice.** It does not follow automatically from variance matching.

6. **Calculate \(V_f\).**

   Define

   \[
   \boxed{
   V_f
   \coloneqq 
   \operatorname{Var}
   (\mu_f\mid \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2).
   }
   \]

   For each preliminary posterior draw, calculate

   \[
   \mu_f^{(b)}
   \coloneqq 
   \frac1N\sum_{s=1}^2\sum_{i=1}^{n_s} f^{(b)}(x_{s,i}).
   \]

   The variance across these averages estimates \(V_f\). It includes uncertainty in the \(f\) trees and leaf parameters, including dependence among predictions at different trial profiles.

7. **Calculate \(V_g\), including the averages over indicators and trees.**

   For fixed discrepancy trees,

   \[
   \mu_g
   =
   \sum_{j,\ell}a_{j\ell}\theta^g_{j\ell},
   \qquad
   a_{j\ell}
   \coloneqq 
   \frac{\#\{\text{trial profiles in leaf }(j,\ell)\}}{N}.
   \]

   Under the mixture prior,

   \[
   z_{j\ell}\mid w\sim\operatorname{Bernoulli}(w),
   \]

   \[
   \theta^g_{j\ell}\mid z_{j\ell},\tau_0^2
   \sim
   N\!\left(
   0,\,
   z_{j\ell}\tau_0^2+(1-z_{j\ell})\tau_1^2
   \right).
   \]

   Conditional independence gives

   \[
   \operatorname{Var}(\mu_g\mid\mathcal T^g,z,\tau_0^2)
   =
   \sum_{j,\ell}a_{j\ell}^2
   \{z_{j\ell}\tau_0^2+(1-z_{j\ell})\tau_1^2\}.
   \]

   First average over \(z\). Every conditional component has mean zero, so the between-component mean term in the law of total variance vanishes. Using \(E(z_{j\ell}\mid w)=w\),

   \[
   \operatorname{Var}(\mu_g\mid\mathcal T^g,w,\tau_0^2)
   =
   \{w\tau_0^2+(1-w)\tau_1^2\}
   \sum_{j,\ell}a_{j\ell}^2.
   \]

   Next average over the specified tree prior, independent of \(w\) and the scales. The conditional mean is again zero for every configuration. Consequently,

   \[
   \boxed{
   V_g(w,\tau_0^2)
   =
   \{w\tau_0^2+(1-w)\tau_1^2\}
   E_{\mathcal T^g}\left[\sum_{j,\ell}a_{j\ell}^2\right],
   }
   \]

   where

   \[
   V_g(w,\tau_0^2)
   \coloneqq 
   \operatorname{Var}
   (\mu_g\mid \mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,w,\tau_0^2).
   \]

   Thus, \(V_g\) is the prior variance of the average discrepancy after integrating out leaf parameters, indicators, and trees.

8. **Substitute the variance decomposition and select the all-spike reference.**

   In the calibration reference, the RWD-updated \(f\) distribution and discrepancy prior are independent conditional on the specified settings. Therefore,

   \[
   \operatorname{Var}
   \left(
   \mu\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,
   w,\tau_0^2
   \right)
   =
   V_f+V_g(w,\tau_0^2).
   \]

   Hence the conditional prior ESS is

   \[
   m(w,\tau_0^2)
   =
   \frac{\sigma_2^2}{V_f+V_g(w,\tau_0^2)}.
   \]

   Here, \(V_g\) adds uncertainty about transferring the historical surface to the trial population. It discounts the information attributed to the historical fit.

   For spike-scale calibration, choose the all-spike reference \(w=1\), retaining the specified tree prior:

   \[
   V_g(1,\tau_0^2)
   =
   \tau_0^2
   E_{\mathcal T^g}\left[\sum_{j,\ell}a_{j\ell}^2\right].
   \]

   Its mean conditional ESS is therefore

   \[
   \boxed{
   \mathrm{ESS}_{\tau_0}(s_0^2)
   \coloneqq 
   E_{\tau_0^2}[m(1,\tau_0^2)]
   =
   E_{\tau_0^2}
   \left[
   \frac{\sigma_2^2}
   {V_f+V_g(1,\tau_0^2)}
   \right].
   }
   \]

9. **Define and derive the ESS ceiling.**

   The ceiling is the largest prior ESS available under the all-spike calibration reference when the discrepancy contributes no uncertainty. Hold \(V_f>0\), \(\sigma_2^2\), \(\nu_0>0\), and the specified tree prior fixed.

   First, generate the scaled inverse-chi-square prior as

   \[
   \tau_0^2=\frac{\nu_0s_0^2}{U},
   \qquad U\sim\chi^2_{\nu_0}.
   \]

   For every fixed \(U>0\), \(\tau_0^2\) approaches zero as \(s_0^2\) approaches zero. This takes a limit through positive prior scales.

   Second, under the all-spike reference,

   \[
   V_g(1,\tau_0^2)
   =
   \tau_0^2
   E_{\mathcal T^g}\left[\sum_{j,\ell}a_{j\ell}^2\right].
   \]

   The tree-weight expectation is fixed and finite: the nonnegative leaf weights sum to one within each tree, so their squared sum across the sum-of-trees is at most \(H_g\). Consequently, \(V_g(1,\tau_0^2)\) approaches zero, and the conditional ESS approaches

   \[
   \frac{\sigma_2^2}{V_f+V_g(1,\tau_0^2)}
   \;\longrightarrow\;
   \frac{\sigma_2^2}{V_f}.
   \]

   Third, passing this limit through the expectation requires justification. Every conditional ESS satisfies

   \[
   0\leq
   \frac{\sigma_2^2}{V_f+V_g(1,\tau_0^2)}
   \leq
   \frac{\sigma_2^2}{V_f}.
   \]

   Using the common \(U\) representation above, the quantities being averaged converge almost surely and are bounded by the same finite constant. Bounded convergence therefore gives

   \[
   \boxed{
   \mathrm{ESS}_{\mathrm{ceiling}}
   \coloneqq 
   \lim_{s_0^2\downarrow0}
   \mathrm{ESS}_{\tau_0}(s_0^2)
   =
   \frac{\sigma_2^2}{V_f}.
   }
   \]

   This explains both the limit and the upper bound. Removing discrepancy uncertainty leaves the historical-fit uncertainty \(V_f\); it does not make the trial control mean known perfectly. Adding nonnegative discrepancy variance increases the denominator and reduces ESS.

   The all-spike condition matters. If \(w<1\) remains fixed, the slab still contributes uncertainty as the spike variance vanishes, so the limit is generally smaller than \(\sigma_2^2/V_f\).

   For example, \(\sigma_2^2=4\) and \(V_f=0.05\) give a ceiling of 80. Under no discrepancy, the historical-data-informed distribution for the trial control mean has the same variance as the flat-prior reference posterior based on 80 independent observations. The f90, f50, and f25 requested targets are then 72, 40, and 20; an absolute target of 100 is capped at 80.

   The ceiling depends on the historical fit, target covariate profiles, and reference variance. It need not equal either cohort's patient count. It bounds this prior calibration quantity, not information after incorporating trial outcomes. A finite positive scale grid generally only approaches the limiting ceiling. For survival, the reference observations concern uncensored log event times, not patients under the trial's censoring pattern.

The all-spike choice concerns calibration; the final joint fit can update both \(f\) and \(g\) using RCT outcomes. Also, averaging conditional ESS values differs from matching the fully marginalized prior variance. The latter would instead give

\[
\frac{\sigma_2^2}
{V_f+E_{\tau_0^2}[V_g(1,\tau_0^2)]}.
\]

In the current implementation, \(\sigma_2^2\) is estimated from RCT controls in two-arm analyses and then held fixed for calibration. The calibration scale law is untruncated, whereas fitting truncates \(\tau_0^2<\tau_1^2\).

**Open numerical implementation issue.** The saved ceiling uses the variance from all preliminary standardized-mean draws, but the saved central ESS curve averages conditional-ratio calculations using separate block estimates of that variance. The block variances are rescaled to have mean equal to the full-draw variance. Because inversion is nonlinear, the averaged curve's zero-scale limit need not equal the saved ceiling and can exceed it. See the Gaussian two-arm [variance calculation](../../lrcbart-sim-gaussian/ess_local/ess_cal.R:210), [ceiling calculation](../../lrcbart-sim-gaussian/ess_local/ess_cal.R:260), and [curve averaging](../../lrcbart-sim-gaussian/ess_local/ess_cal.R:320). Alignment of the central numerical curve with the same \(V_f\) remains open; saving this derivation does not change calibration code, targets, or results.
