---
title: "Comparison of Prior Effective Sample Size Definitions"
subtitle: "Neuenschwander variance calibration, Morita–Thall–Müller, ELIR, and the proposed variance-matching definition"
date: "2026-09-18"
lang: en
---

<style>
body {
  max-width: 1800px;
  margin: 0 auto;
  padding: 1.5rem;
  color: #17202a;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  line-height: 1.4;
}

table {
  width: 100%;
  border-collapse: collapse;
  margin: 1.1rem 0 1.6rem;
}

table:first-of-type {
  display: table;
  table-layout: fixed;
  overflow: visible;
  font-size: 0.8rem;
}

table:first-of-type th,
table:first-of-type td {
  min-width: 0;
  max-width: none;
  padding: 0.45rem;
  overflow-wrap: anywhere;
  word-break: normal;
  hyphens: auto;
}

table:first-of-type math {
  font-size: 0.9em;
}

th,
td {
  padding: 0.65rem;
  border: 1px solid #d5d8dc;
  vertical-align: top;
}

thead {
  background: #eaf2f8;
}

table:first-of-type thead tr:first-child th {
  background: #d6eaf8;
  text-align: center;
  vertical-align: middle;
  font-size: 0.82rem;
  letter-spacing: 0.02em;
}

table:first-of-type thead tr:nth-child(2) th {
  background: #eaf2f8;
}

table:first-of-type tbody tr:nth-child(1),
table:first-of-type tbody tr:nth-child(2) {
  background: #f4f9fc;
}

table:first-of-type tbody tr:nth-child(3) {
  background: #fbf7ed;
}

table:first-of-type tbody tr:nth-child(4) {
  background: #f1f8f3;
}

table:first-of-type tbody tr:nth-child(5) {
  background: #f7f2fa;
}


@media print {
  @page { size: A3 landscape; margin: 10mm; }
  body { max-width: none; padding: 0; }
  table:first-of-type { font-size: 7.4pt; }
}
</style>

This table puts four prior-ESS approaches into a common notation for a scalar parameter. Throughout, \(\theta_h\) denotes a historical-data parameter or estimand, with \(h\) indexing historical studies when applicable, and \(\theta_\star\) denotes the corresponding current-study parameter or estimand. The two Neuenschwander rows make the calibration step explicit: the complete-pooling MAP prior supplies the known ESS anchor, and the partial-pooling MAP prior is the prior whose ESS is inferred.

<table class="ess-comparison">
<colgroup>
<col style="width: 9%">
<col style="width: 11%">
<col style="width: 11%">
<col style="width: 10%">
<col style="width: 12%">
<col style="width: 12%">
<col style="width: 12%">
<col style="width: 11%">
<col style="width: 12%">
</colgroup>
<thead>
<tr>
<th rowspan="2">Method / case</th>
<th colspan="3">Target</th>
<th colspan="3">Reference</th>
<th colspan="2">ESS</th>
</tr>
<tr>
<th>Target prior</th>
<th>Information measure</th>
<th>Information quantity</th>
<th>Information source</th>
<th>Information measure</th>
<th>Information quantity</th>
<th>Matching / calibration</th>
<th>Result</th>
</tr>
</thead>
<tbody>
<tr>
<td>**Neuenschwander et al. (2010): complete pooling**</td>
<td>$\pi_0(\theta_\star)\coloneqq p(\theta_\star\mid D_H,\tau=0)$. <br>where $D_H\coloneqq(Y_1,\ldots,Y_H)$ and $\pi_0$ is the complete-pooling MAP prior.</td>
<td>Inverse marginal variance:<br>$V_0\coloneqq\operatorname{Var}_{\pi_0}(\theta_\star)$. <br>where $V_0$ is the complete-pooling prior variance; for a normal prior, inverse variance equals log-prior curvature.</td>
<td>$\mathcal I_0\coloneqq1/V_0$. <br>where $\mathcal I_0$ is the complete-pooling prior information.</td>
<td>–</td>
<td>–</td>
<td>–</td>
<td>Calibration anchor: $ESS_0\coloneqq N_H$. <br>where $N_H$ is the historical sample size assigned to the complete-pooling prior.</td>
<td>$ESS_0=N_H$. <br>where $ESS_0$ is the anchor ESS.</td>
</tr>
<tr>
<td>**Neuenschwander et al. (2010): partial pooling**</td>
<td>$\pi_{\mathrm{target}}(\theta_\star)\coloneqq p(\theta_\star\mid D_H)$. <br>where $\pi_{\mathrm{target}}$ is the partial-pooling MAP prior.</td>
<td>Inverse marginal variance:<br>$V_{\mathrm{target}}\coloneqq\operatorname{Var}_{\pi_{\mathrm{target}}}(\theta_\star)$. <br>where $V_{\mathrm{target}}$ is the partial-pooling prior variance and the same measure is used for both MAP priors.</td>
<td>$\mathcal I_{\mathrm{target}}\coloneqq1/V_{\mathrm{target}}$. <br>where $\mathcal I_{\mathrm{target}}$ is information in the partial-pooling MAP prior.</td>
<td>–</td>
<td>–</td>
<td>–</td>
<td>$\mathcal I_{\mathrm{target}}$<br>$\displaystyle =\frac{ESS_{\mathrm{target}}}{N_H}\mathcal I_0$. <br>where $ESS_{\mathrm{target}}/N_H$ is the information retained relative to complete pooling.</td>
<td>$\displaystyle ESS_{\mathrm{target}}=N_H\frac{V_0}{V_{\mathrm{target}}}$. <br>where $V_0/V_{\mathrm{target}}$ is the relative-precision factor.</td>
</tr>
<tr>
<td>**Morita–Thall–Müller (2008)**</td>
<td>$\pi_{\mathrm{target}}(\theta_\star)$. <br>where $\theta_\star$ is the current-study parameter and $\pi_{\mathrm{target}}$ is the supplied prior whose ESS is sought; the definition does not prescribe how this prior is constructed.</td>
<td>Log-prior curvature at $\bar\theta_\star$:<br>$\bar\theta_\star\coloneqq E_{\pi_{\mathrm{target}}}(\theta_\star)$<br>$\displaystyle i\{q(\theta_\star)\}\coloneqq-\frac{d^2}{d\theta_\star^2}\log q(\theta_\star)$. <br>where $q\coloneqq\pi_{\mathrm{target}}$.</td>
<td>$\mathcal I_{\mathrm{target}}\coloneqq i\{\pi_{\mathrm{target}}(\bar\theta_\star)\}$<br>$\displaystyle =-\left.\frac{d^2}{d\theta_\star^2}\log\pi_{\mathrm{target}}(\theta_\star)\right\rvert_{\theta_\star=\bar\theta_\star}$. <br>where $\mathcal I_{\mathrm{target}}$ is target-prior curvature evaluated at its mean.</td>
<td>$\pi_{\epsilon,m}(\theta_\star\mid Y_m)$<br>$\displaystyle \propto\pi_\epsilon(\theta_\star)\prod_{j=1}^m f(Y_j\mid\theta_\star)$. <br>where $\pi_\epsilon$ is a same-mean $\epsilon$-information prior and $Y_m$ is a hypothetical sample of size $m$.</td>
<td>Expected log-posterior curvature at $\bar\theta_\star$:<br>$\displaystyle i\{q(\theta_\star)\}\coloneqq-\frac{d^2}{d\theta_\star^2}\log q(\theta_\star)$. <br>where $q\coloneqq\pi_{\epsilon,m}$ and the curvature is averaged over $Y_m$.</td>
<td>$\mathcal I_{\mathrm{ref}}(m)\coloneqq i\{\pi_\epsilon(\bar\theta_\star)\}$<br>$\displaystyle {}+E_{Y_m}\{i_F(Y_m;\bar\theta_\star)\}$. <br>where $i_F(Y_m;\theta_\star)$ is observed likelihood information and the expectation uses the target prior-predictive distribution.</td>
<td>$\mathcal I_{\mathrm{target}}\approx\mathcal I_{\mathrm{ref}}(m)$. <br>where $m$ is varied to obtain the closest information match.</td>
<td>$\displaystyle ESS_{\mathrm{MTM}}\coloneqq\arg\min_{m\ge0}d(m)$. <br>where $d(m)\coloneqq\lvert\mathcal I_{\mathrm{target}}-\mathcal I_{\mathrm{ref}}(m)\rvert$; the original definition uses integer $m$ and interpolation permits a noninteger ESS.</td>
</tr>
<tr>
<td>**ELIR (Neuenschwander et al., 2020)**</td>
<td>$\pi_{\mathrm{target}}(\theta_\star)$. <br>where $\theta_\star$ is the current-study parameter and $\pi_{\mathrm{target}}$ is the prior whose ESS is sought; ELIR does not restrict how this prior is constructed.</td>
<td>Pointwise log-prior curvature:<br>$\displaystyle i\{\pi(\theta_\star)\}\coloneqq-\frac{d^2}{d\theta_\star^2}\log\pi(\theta_\star)$. <br>where the curvature is evaluated locally at each $\theta_\star$.</td>
<td>$\mathcal I_{\mathrm{target}}(\theta_\star)\coloneqq i\{\pi_{\mathrm{target}}(\theta_\star)\}$<br>$\displaystyle =-\frac{d^2}{d\theta_\star^2}\log\pi_{\mathrm{target}}(\theta_\star)$. <br>where $\mathcal I_{\mathrm{target}}(\theta_\star)$ is local target-prior information.</td>
<td>One information unit $Y_1\sim f(\cdot\mid\theta_\star)$; no reference prior. <br>where $Y_1$ is one unit under the current-study sampling model $f$.</td>
<td>Expected Fisher information from one observation:<br>$\displaystyle i_F(\theta_\star)\coloneqq E_{Y_1\mid\theta_\star}\!\left[-\frac{d^2}{d\theta_\star^2}\log f(Y_1\mid\theta_\star)\right]$. <br>where the expectation is under the current-study sampling model.</td>
<td>$\mathcal I_{\mathrm{ref}}(\theta_\star)\coloneqq i_F(\theta_\star)$<br>$\displaystyle =E_{Y_1\mid\theta_\star}\!\left[-\frac{d^2}{d\theta_\star^2}\log f(Y_1\mid\theta_\star)\right]$. <br>where $i_F(\theta_\star)$ is expected Fisher information in one unit.</td>
<td>$\mathcal I_{\mathrm{target}}(\theta_\star)$<br>$=m(\theta_\star)\mathcal I_{\mathrm{ref}}(\theta_\star)$. <br>where $m(\theta_\star)$ is the local equivalent sample size.</td>
<td>$\displaystyle ESS_{\mathrm{ELIR}}\coloneqq E_{\pi_{\mathrm{target}}}\{m(\theta_\star)\}$<br>$\displaystyle =E_{\pi_{\mathrm{target}}}\!\left\{\frac{\mathcal I_{\mathrm{target}}(\theta_\star)}{\mathcal I_{\mathrm{ref}}(\theta_\star)}\right\}$. <br>where the expectation averages local ESS over the target prior; the original paper denotes $m(\theta_\star)$ by $r(\theta_\star)$.</td>
</tr>
<tr>
<td>**Proposed conditional variance-matching ESS**</td>
<td>$\theta_\star\coloneqq\theta_h+\delta_\star$,<br>$\displaystyle \pi_{\mathrm{target}}(\theta_\star\mid w,\tau_0^2)$<br>$\displaystyle \coloneqq\int p_g(\theta_\star-\theta_h\mid w,\tau_0^2)$<br>$\displaystyle {}\times p(\theta_h\mid D_H)\,d\theta_h$. <br>where $\theta_h$ and $\delta_\star$ are the standardized historical-control and discrepancy estimands induced by the posterior of $f$ and prior for $g$, respectively.</td>
<td>Inverse conditional marginal variance:<br>$V_{\mathrm{target}}(w,\tau_0^2)\coloneqq\operatorname{Var}_{\pi_{\mathrm{target}}}(\theta_\star\mid w,\tau_0^2)$. <br>where $V_{\mathrm{target}}(w,\tau_0^2)$ is the conditional target-prior variance.</td>
<td>$\displaystyle \mathcal I_{\mathrm{target}}(w,\tau_0^2)\coloneqq1/V_{\mathrm{target}}(w,\tau_0^2)$. <br>where $\mathcal I_{\mathrm{target}}(w,\tau_0^2)$ is the target information conditional on the discrepancy-prior settings.</td>
<td>One information unit $Y_1\sim f(\cdot\mid\theta_\star)$; no reference prior. <br>where $Y_1$ is one unit under the current-study sampling model $f$.</td>
<td>Expected Fisher information from one observation:<br>$\displaystyle i_F(\theta_\star)\coloneqq E_{Y_1\mid\theta_\star}\!\left[-\frac{d^2}{d\theta_\star^2}\log f(Y_1\mid\theta_\star)\right]$. <br>where the expectation is under the current-study sampling model.</td>
<td>$\mathcal I_{\mathrm{ref}}(\theta_\star)\coloneqq i_F(\theta_\star)$<br>$\displaystyle =E_{Y_1\mid\theta_\star}\!\left[-\frac{d^2}{d\theta_\star^2}\log f(Y_1\mid\theta_\star)\right]$. <br>where $i_F(\theta_\star)$ is expected Fisher information in one reference observation.</td>
<td>$\mathcal I_{\mathrm{target}}(w,\tau_0^2)$<br>$=m(w,\tau_0^2;\theta_\star)\mathcal I_{\mathrm{ref}}(\theta_\star)$. <br>where $m(w,\tau_0^2;\theta_\star)$ is the conditional equivalent reference sample size.</td>
<td>$\displaystyle m(w,\tau_0^2;\theta_\star)=\frac{\mathcal I_{\mathrm{target}}(w,\tau_0^2)}{i_F(\theta_\star)}$<br>$\displaystyle ESS(w,s_0^2)\coloneqq E_{\tau_0^2}\{m(w,\tau_0^2;\theta_\star)\}$. <br>where the first equality solves the matching equation and the expectation uses the calibration distribution indexed by $s_0^2$; for the normal reference, $i_F$ is constant and the dependence on $\theta_\star$ drops out.</td>
</tr>
</tbody>
</table>

## Normal–normal comparison

The following derivations isolate the reference information used by each method under normal sampling models. Each subsection starts from that method's reference construction and ends with its reference information quantity.

### Neuenschwander et al. (2010): complete-pooling calibration

Neuenschwander's normal random-effects model starts from (H) historical trial estimates,

$$
Y_h\mid\theta_h\sim N(\theta_h,\sigma_h^2),
\qquad h=1,\ldots,H,
$$

with exchangeable historical and new-trial parameters

$$
\theta_1,\ldots,\theta_H,\theta_\star
\mid\mu,\tau^2
\overset{\mathrm{iid}}{\sim}
N(\mu,\tau^2).
$$

For fixed \(\tau\), first integrate out \(\theta_h\):

$$
Y_h\mid\mu,\tau
\sim N(\mu,\sigma_h^2+\tau^2).
$$

Under the paper's locally uniform prior for \(\mu\), define

$$
w_h(\tau)\coloneqq\frac1{\sigma_h^2+\tau^2},
\qquad
W_\tau\coloneqq\sum_{h=1}^H w_h(\tau),
\qquad
\widehat\mu_\tau
\coloneqq\frac{\sum_{h=1}^H w_h(\tau)Y_h}{W_\tau}.
$$

Normal conjugacy gives the posterior distribution of the common location:

$$
\mu\mid D_H,\tau
\sim N\!\left(\widehat\mu_\tau,\frac1{W_\tau}\right).
$$

Next use \(\theta_\star\mid\mu,\tau\sim N(\mu,\tau^2)\) and integrate over the posterior of \(\mu\). The meta-analytic predictive distribution is

$$
\theta_\star\mid D_H,\tau
\sim
N\!\left(\widehat\mu_\tau,V_\tau\right),
$$

where

$$
V_\tau
\coloneqq\frac1{W_\tau}+\tau^2
=\frac{1}{\displaystyle\sum_{h=1}^H
\frac1{\sigma_h^2+\tau^2}}
+\tau^2.
$$

If the variance of historical trial estimate \(h\) is \(\sigma_h^2=\sigma^2/n_h\), then

$$
w_h(\tau)
=\frac{1}{\sigma^2/n_h+\tau^2}
=\frac{n_h}{\sigma^2+n_h\tau^2},
$$

and, for possibly unequal trial sizes,

$$
V_\tau
=\frac{1}{\displaystyle\sum_{h=1}^H
\frac{n_h}{\sigma^2+n_h\tau^2}}
+\tau^2.
$$

If all \(H\) historical trials have size \(n_h=n\), define the total historical sample size as \(N_H\coloneqq Hn\). Then

$$
V_\tau
=\frac{\sigma^2+n\tau^2}{Hn}+\tau^2
=\frac{\sigma^2}{N_H}
+\tau^2\left(1+\frac1H\right).
$$

Complete pooling is the special case \(\tau=0\). Hence

$$
V_0
\coloneqq\operatorname{Var}(\theta_\star\mid D_H,\tau=0)
=\frac{1}{\displaystyle\sum_{h=1}^H\frac1{\sigma_h^2}},
\qquad
\mathcal I_0\coloneqq\frac1{V_0}
=\sum_{h=1}^H\frac1{\sigma_h^2}.
$$

When \(\sigma_h^2=\sigma^2/n_h\), define \(N_H\coloneqq\sum_{h=1}^Hn_h\). Then

$$
V_0=\frac{\sigma^2}{N_H},
\qquad
\mathcal I_0=\frac{N_H}{\sigma^2}.
$$

The complete-pooling prior is assigned its historical count \(N_H\). Therefore, its implied information per reference patient is

$$
\boxed{
\mathcal I_{\mathrm{ref},1}
\coloneqq\frac{\mathcal I_0}{N_H}
=\frac1{N_HV_0}
=\frac1{\sigma^2}.
}
$$

### Morita–Thall–Müller (2008): reference-posterior curvature

Let the normal reference likelihood and the same-mean \(\epsilon\)-information prior be

$$
Y_j\mid\theta_\star\overset{\mathrm{iid}}{\sim}N(\theta_\star,\sigma^2),
\qquad j=1,\ldots,m,
\qquad
\pi_\epsilon\coloneqq N(\mu_0,V_\epsilon).
$$

The resulting reference posterior is normal. Its variance is

$$
V_{\epsilon,m}
\coloneqq\left(\frac1{V_\epsilon}+\frac{m}{\sigma^2}\right)^{-1}.
$$

For a normal posterior, log-posterior curvature equals inverse variance. Therefore, Morita's reference information quantity is

$$
\boxed{
\mathcal I_{\mathrm{ref}}(m)
\coloneqq\frac1{V_{\epsilon,m}}
=\frac1{V_\epsilon}+\frac{m}{\sigma^2}.
}
$$

The posterior variance does not depend on the realized hypothetical observations, so no further averaging over their values is needed in this normal case.

### ELIR: one-observation Fisher information

For one observation from the normal reference likelihood,

$$
Y_1\mid\theta_\star\sim N(\theta_\star,\sigma^2),
$$

the expected Fisher information is

$$
\boxed{
\mathcal I_{\mathrm{ref}}(\theta_\star)
\coloneqq i_F(\theta_\star)
=E_{Y_1\mid\theta_\star}\!\left[
-\frac{d^2}{d\theta_\star^2}\log f(Y_1\mid\theta_\star)
\right]
=\frac1{\sigma^2}.
}
$$

It is constant in \(\theta_\star\). Consequently, \(m\) independent reference observations provide \(m/\sigma^2\) information. ELIR uses no reference prior.

### Proposed conditional variance-matching definition

The proposed method uses the same one-observation reference source as ELIR:

$$
Y_1\mid\theta_\star\sim N(\theta_\star,\sigma^2).
$$

Thus its reference information quantity is also

$$
\boxed{
\mathcal I_{\mathrm{ref}}(\theta_\star)
\coloneqq i_F(\theta_\star)
=\frac1{\sigma^2},
}
$$

and \(m\) independent reference observations provide

$$
m\,\mathcal I_{\mathrm{ref}}(\theta_\star)
=\frac{m}{\sigma^2}.
$$

No reference prior is required. In the manuscript specialization, the table's \(\theta_\star\) is \(\mu\), and the reference variance is \(\sigma_2^2\), so \(i_F(\theta_\star)=1/\sigma_2^2\).

## References

1. Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ. *Summarizing historical information on controls in clinical trials.* Clinical Trials. 2010;7(1):5–18. [doi:10.1177/1740774509356002](https://doi.org/10.1177/1740774509356002).
2. Morita S, Thall PF, Müller P. *Determining the effective sample size of a parametric prior.* Biometrics. 2008;64(2):595–602. [doi:10.1111/j.1541-0420.2007.00888.x](https://doi.org/10.1111/j.1541-0420.2007.00888.x).
3. Neuenschwander B, Weber S, Schmidli H, O’Hagan A. *Predictively consistent prior effective sample sizes.* Biometrics. 2020;76(2):578–587. [doi:10.1111/biom.13252](https://doi.org/10.1111/biom.13252).
