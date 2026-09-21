---
pagetitle: "Section 2.5: theoretical objectives and candidate results"
output:
  html_document:
    toc: true
    toc_depth: 3
    math_method: mathml
---

# Section 2.5: agreed theoretical objectives and candidate results

19 September 2026. Development note. The fixed-region statements and proofs below have been incorporated into Section 2.5 and Appendix C. They remain mathematical drafts requiring human review, as recorded in [HANDOFF.md](../HANDOFF.md). The broader full-model theorem remains a research objective.

## Agreed objectives

1. Establish learning of the scientifically meaningful total source discrepancy and of regions of practical agreement.
2. Characterize the estimation benefit of historical borrowing under agreement and the estimation error it can introduce under conflict, including implications for the standardized treatment contrast.

The existing fixed-leaf results are supporting results about the mixture mechanism. The proposed main results should concern total control means and total discrepancies, with historical uncertainty included.

## Candidate Theorem 2: learning local agreement and the current-control surface

Write the true historical mean as $f_*$, the true current-control mean as $m_{1,*}$, and $g_*=m_{1,*}-f_*$. Let $Q$ be a prespecified target covariate distribution on a region supported by both control populations.

The desired full-model conclusion is joint posterior consistency of $f$ and $m_1=f+g$ in $L^2(Q)$. Explicitly, for every $\epsilon>0$,

$$
\Pi\{\|f-f_*\|_{L^2(Q)}+\|f+g-m_{1,*}\|_{L^2(Q)}>\epsilon\mid\mathcal D\}
\longrightarrow 0
$$

in probability under the true sampling distribution. This implies consistency of the total $g$. If $Q\{|g_*|=c\}=0$, the intended local-summary consequence is

$$
\int\left|B_c(x)-\mathbf{1}\{|g_*(x)|<c\}\right|\,dQ(x)
\longrightarrow 0.
$$

In words, the posterior probability of practical agreement approaches one where the true source difference is smaller than $c$ and zero where it is larger, with the approximation error averaged over the target population. This states the result directly in terms of posterior evidence, without assigning labels or introducing a 50% decision cutoff. The displayed full-model target is population-weighted; pointwise or uniform conclusions require separate arguments.

With analogous $L^2(Q)$ consistency of the treatment surface, Cauchy--Schwarz yields consistency of the $Q$-standardized treatment contrast. The empirical-trial standardization used in the manuscript needs its corresponding design argument; it is not interchangeable with fixed $Q$ without justification.

### Assumptions and proof status

This is a full-model theorem target, not an established theorem for the implemented lrcBART prior. The proof must specify normal sampling, overlap (for example bounded density ratios of $Q$ relative to both control designs), truth regularity/approximability, tree-prior support and complexity, prior scales, and growth of both control samples. A first full-model argument may assume comparable sample growth; extreme historical-to-trial imbalance requires an explicit extension.

One must verify these properties for the actual depth prior, admissible splitting rules, hyperpriors, and data-dependent scale calibration. Merely assuming joint posterior contraction and proving the displayed consequence would establish a useful corollary, but would not establish the proposed main result. Existing BART results cannot be transferred automatically: Ročková and Saha study the tree prior and use a modification to obtain optimal concentration ([On Theory for BART](https://proceedings.mlr.press/v89/rockova19a.html)).

For clarity, the local-summary implication follows by excluding the set where $\bigl||g_*|-c\bigr|\leq\eta$. Outside that set, a posterior draw can put $|g|$ on the opposite side of $c$ from $|g_*|$ only if $|g-g_*|\geq\eta$. The integrated absolute error of $B_c$ relative to the true agreement indicator is at most

$$
Q\{\bigl||g_*|-c\bigr|\leq\eta\}
+\Pi(\|g-g_*\|_{L^2(Q)}>\epsilon\mid\mathcal D)
+\epsilon^2/\eta^2.
$$

First take the large-sample limit, then $\epsilon$ down to zero, then $\eta$ down to zero. This argument uses posterior probabilities and does not assume convergence of posterior second moments.

A tractable first version uses finitely many fixed regions, constant true means within each region, known residual variances, fixed positive discrepancy scales, and observations from both sources increasing in every target region. It jointly estimates both unknown source means in all regions. In that finite-dimensional model, normal-mixture calculations establish consistency and practical-agreement learning directly. This is more informative about borrowing than fixing $f$ and all other contributions, but remains a regional benchmark rather than full BART theory.

There is no outcome-based discrepancy-learning claim in a single-arm trial under the separate treatment model.

### First rigorous version: joint learning in fixed regions

This version establishes learning of total regional source means and discrepancies, including uncertainty in both sources. Its partition is fixed in advance. It does not assert the full-model target above for learned tree ensembles.

**Theorem (Learning local practical agreement in a joint regional model).** Let $R_1,\ldots,R_J$ be a fixed finite partition of the target covariate domain. Within each region $j$, suppose independent historical-control and trial-control outcomes follow

$$
Y_{2,ji}\sim N(f_j^*,\sigma_2^2),\qquad
Y_{1,ji}\sim N(f_j^*+g_j^*,\sigma_1^2),
$$

where the positive residual variances are known. Fit the corresponding model with independent priors across regions and independent $f_j$ and $g_j$ within each region,

$$
f_j\sim N(0,\lambda_j^2),\qquad
g_j\sim wN(0,\tau_0^2)+(1-w)N(0,\tau_1^2),
$$

with all prior parameters fixed, $0<\lambda_j^2<\infty$, $0<w<1$ and $0<\tau_0^2<\tau_1^2<\infty$. If $n_{1,j}$ and $n_{2,j}$ both tend to infinity for every $j$, then for every $\epsilon>0$,

$$
\Pi\left\{\max_{1\leq j\leq J}
\bigl(|f_j-f_j^*|+|g_j-g_j^*|\bigr)>\epsilon
\mid\mathcal D\right\}\longrightarrow 0
$$

almost surely under the true sampling law. Consequently, for a fixed $c>0$ such that $|g_j^*|\ne c$ for every target region,

$$
B_{c,j}\mathrel{:=}\Pi(|g_j|<c\mid\mathcal D)
\longrightarrow
\begin{cases}
1,& |g_j^*|<c,\\
0,& |g_j^*|>c,
\end{cases}
$$

almost surely, simultaneously over all $j$. The current-control means $f_j+g_j$ are also consistently learned. For $x\in R_j$, the total regional function is $g(x)=g_j$ and $B_c(x)=B_{c,j}$.

**Proof.** Fix a region and put $v_s=\sigma_s^2/n_{s,j}$. The vector of sample statistics

$$
\widehat\beta_j=(\bar Y_{2,j},\bar Y_{1,j}-\bar Y_{2,j})^\top
$$

has mean $\beta_j^*=(f_j^*,g_j^*)^\top$ and covariance

$$
\Sigma_j=\begin{pmatrix}v_2&-v_2\\-v_2&v_1+v_2\end{pmatrix}.
$$

By the strong law, $\widehat\beta_j$ tends to $\beta_j^*$ almost surely, and $\Sigma_j$ tends to the zero matrix. Conditional on discrepancy component $k$, the prior covariance is $P_{jk}=\operatorname{diag}(\lambda_j^2,\tau_k^2)$. Normal updating gives posterior covariance and mean

$$
C_{jk}=(\Sigma_j^{-1}+P_{jk}^{-1})^{-1},\qquad
m_{jk}=\widehat\beta_j-C_{jk}P_{jk}^{-1}\widehat\beta_j.
$$

Since $0\prec C_{jk}\preceq\Sigma_j$ in the positive-semidefinite order, $C_{jk}$ tends to zero. Because $P_{jk}$ is fixed and positive definite, $m_{jk}$ tends to $\beta_j^*$ almost surely. Each of the two normal posterior components therefore concentrates at the true pair of means, and so does their posterior mixture, whatever its mixing probabilities. A finite union bound proves the simultaneous result across regions. Positive separation of every $|g_j^*|$ from $c$ then yields the agreement-probability limits. The result for $f_j+g_j$ follows by the triangle inequality.

This proof does not hold historical predictions or other regional means fixed. It also does not require recovery of a spike/slab allocation. The mixture components overlap and either component can learn the true discrepancy. The role of this theorem is to justify the scientific posterior summary; the estimation benefit specific to borrowing is the subject of the next result.

## Candidate Theorem 3: regional borrowing benefit and conflict risk

### Exact benchmark model

Consider one fixed region, with independent sample means

$$
\bar Y_1\sim N(\mu_1,v_1),\qquad
\bar Y_2\sim N(\mu_2,v_2),\qquad
v_s=\sigma_s^2/n_s>0.
$$

Let $d=\mu_1-\mu_2$. Use a flat prior on $\mu_2$ and the independent discrepancy prior

$$
d\sim wN(0,t_0)+(1-w)N(0,t_1),\qquad
0<w<1,\quad 0<t_0<t_1<\infty.
$$

Here $t_k$ denotes a variance, not a standard deviation. All scales and $w$ are fixed. The flat historical-mean prior is a deliberate benchmark simplification; the fitted BART model has proper Gaussian leaf priors. This benchmark puts the mixture on the total regional discrepancy, not on an individual contribution within an arbitrary multi-tree decomposition.

Set $D=\bar Y_1-\bar Y_2$, $V=v_1+v_2$, and

$$
\begin{aligned}
\gamma(D)&=\frac{w\phi(D;0,V+t_0)}{w\phi(D;0,V+t_0)+(1-w)\phi(D;0,V+t_1)},\\[4pt]
a_k&=\frac{v_1}{V+t_k},\\[4pt]
a(D)&=\gamma(D)a_0+\{1-\gamma(D)\}a_1.
\end{aligned}
$$

The posterior mean of the current-control mean is exactly

$$
\widehat\mu_1=\bar Y_1-a(D)D
=\{1-a(D)\}\bar Y_1+a(D)\bar Y_2.
$$

Both sources' sampling uncertainties are included. The coefficient $a(D)$ describes the estimator in this benchmark; it is not $w$, a patient likelihood weight, $B_c$, or the manuscript's prior ESS.

### Conclusions and derivation

1. **Adaptive borrowing coefficient.** The coefficient $a(D)$ decreases strictly with $|D|>0$. It is bounded between $a_1$ and $a_0$ and tends to $a_1$ under extreme observed disagreement. This result describes the current-control estimate rather than a conditional discrepancy-leaf update.

2. **Strict improvement in squared-error risk under exact agreement.** When the true $d=0$, the posterior-mean estimator is unbiased and has mean squared error strictly less than $v_1$, the risk of the trial-only sample mean. This is a repeated-sampling MSE result, not a claim that every realized posterior interval is shorter.

   To prove it, put $a_*=v_1/V$ and $T=(v_2\bar Y_1+v_1\bar Y_2)/V$. Under normal sampling, $T$ and $D$ are independent, $\operatorname{Var}(T)=v_1v_2/V$, and $E(T)=\mu_1-a_*d$. Since $\widehat\mu_1=T+\{a_*-a(D)\}D$, at $d=0$ symmetry gives zero bias and

   $$
   R(0)=\frac{v_1v_2}{V}
   +E_0[\{a_*-a(D)\}^2D^2]
   <\frac{v_1v_2}{V}+a_*^2V=v_1.
   $$

   The strict inequality uses $0<a(D)<a_*$ and $P(D\ne0)=1$. Continuity of $R(d)$ also gives an open neighborhood of zero with lower MSE than trial-only, although its size depends on the specified parameters.

3. **Exact conflict-risk characterization.** For any fixed true discrepancy $d$,

   $$
   \operatorname{Bias}_d(\widehat\mu_1)=-E_d[a(D)D],
   $$

   $$
   R(d)=\frac{v_1v_2}{V}
   +E_{D\sim N(d,V)}[\{(a_*-a(D))D-a_*d\}^2].
   $$

   These one-dimensional expectations quantify the tradeoff using both sample sizes, both observation variances, discrepancy size, and prior parameters. Do not substitute the data-dependent $a(D)$ into a fixed-weight bias-variance formula. An equivalent identity, obtained by normal integration by parts, is

   $$
   R(d)=v_1+E_d[a(D)^2D^2]-2v_1E_d[a(D)+Da'(D)].
   $$

4. **Limit of protection at fixed sample size.** Since $a_1>0$, the shift from the trial-only estimate approaches $-a_1D$ as $|D|$ increases, and does not disappear. For fixed $v_1$, $v_2$, $t_0$, $t_1$ and $w$, as $|d|$ goes to infinity,

   $$
   \begin{aligned}
   \operatorname{Bias}_d(\widehat\mu_1)&=-a_1d+o(1),\\[4pt]
   R(d)&=a_1^2d^2+(1-a_1)^2v_1+a_1^2v_2+o(1).
   \end{aligned}
   $$

   Gaussian tail attenuation of $\gamma$ makes all the remaining polynomial-times-$\gamma$ expectations vanish. Thus no uniform finite bound on borrowing bias or risk over arbitrarily large discrepancies follows for this fixed-variance normal-mixture benchmark. Increasing trial size is a different limit: with fixed positive scales and fixed true $d$, $a(D)\leq v_1/t_0$, so its borrowing-induced bias vanishes as $v_1$ goes to zero (for bounded $v_2$).

These identities and the agreement-risk inequality have been derived here for the benchmark. They are not yet results for random tree ensembles, estimated shared hyperparameters, proper baseline shrinkage, or empirical scale selection. Deterministic normal quadrature checked the two risk identities in 15 parameter/discrepancy configurations; maximum relative numerical difference was 1.4e-14. This is an algebra check, not a substitute for the derivations above.

### Connection to the treatment contrast

For fixed regional target weights $p_j$, the control contribution to contrast-estimation error is minus the weighted sum of regional control errors. Its MSE contains squared weighted bias and all weighted covariance terms. Under independent regions and fixed independent regional priors, those covariances vanish. If every region has exact agreement, unbiased regional estimators and the strict regional MSE improvement yield an improvement for the standardized control mean, and also for a contrast using the same independent treatment estimator. In conflicting regions the bias terms must be retained. Shared $w$, shared scales, and overlapping trees can induce dependence and require covariance accounting rather than automatic regional additivity.

## Intended manuscript development

The first objective is a learning theorem for total functions, with learning of posterior practical-agreement probabilities and estimand consistency as consequences. The second is an estimation-risk theorem, initially with an exact regional benchmark and then an extension as far as the actual joint model permits. A universal guarantee of no bias under arbitrary conflict, or of improved estimation for every possible discrepancy, is not a proposed claim.

The regional risk calculation is a useful foundation, but its novelty alone is limited by its reduction to a normal commensurate model. The contribution should come from a justified connection to covariate-specific total discrepancies and the standardized estimand, and from a defensible extension to the model used in the paper.
