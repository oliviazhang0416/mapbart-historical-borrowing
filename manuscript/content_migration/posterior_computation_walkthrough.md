# Posterior computation walkthrough

Notation: \(1\) denotes trial treatment, \(2\) trial control and \(3\) historical control. Standardization uses \(\mathbf{X}_1\) and \(\mathbf{X}_2\), with \(N\coloneqq n_1+n_2\). Each participant contributes once. In a single-arm trial, \(n_2=0\), so standardization uses only \(\mathbf{X}_1\).

Status: settled with the user on 2026-09-16.

Purpose: authoritative content reference for the posterior-computation section of the updated manuscript. This preserves the approved walkthrough of marginal leaf likelihoods, their role in Metropolis–Hastings tree updates, and the subsequent conditional leaf draws. Both leaf factors use A and B; the derivation omits density-ratio notation and does not introduce v.

No analysis code or original manuscript sources were changed when saving this section. Existing implementation-alignment issues remain recorded in [implementation qualifications](implementation_notes.md).

## 1. Update one tree conditional on the rest of the model

Hold all other trees, their leaf parameters, and the variance and mixture parameters fixed.

For tree \(j\), form the partial residuals

\[
\begin{aligned}
r_{3,i}^f&\coloneqq Y_{3,i}-f^{(-j)}(x_{3,i}),\\
r_{2,i}^f&\coloneqq Y_{2,i}-f^{(-j)}(x_{2,i})-g(x_{2,i}),\\
r_{2,i}^g&\coloneqq Y_{2,i}-f(x_{2,i})-g^{(-j)}(x_{2,i}).
\end{aligned}
\]

The RCT residuals use only trial controls. Here, \(f^{(-j)}\) and \(g^{(-j)}\) exclude the tree being updated. These residuals remain fixed when comparing candidate structures for that tree.

In the generic leaf integral below, source labels are suppressed on the residuals only. Within each candidate leaf, integrate out its parameter \(\theta\). For a \(g\) leaf, also sum over its indicator \(z\). The resulting **marginal leaf likelihood** is

\[
\boxed{
p(\mathbf r_\ell\mid \mathcal T_j,\text{rest})
=
\left[\prod_{i\in\ell}\phi(r_i;0,\sigma_i^2)\right]M_\ell,
}
\]

where \(\sigma_i^2\) is the appropriate source variance and \(M_\ell\) is the leaf factor derived below. The final argument of \(\phi\) denotes variance.

## 2. Derive the f-leaf factor

For an \(f\) leaf,

\[
r_i\mid\theta\sim N(\theta,\sigma_i^2),
\qquad
\theta\sim N(0,\lambda_f^2).
\]

Let \(n_{2,\ell},n_{3,\ell}\) be the numbers of RCT and RWD controls in the leaf, and \(S_{2,\ell},S_{3,\ell}\) their respective partial-residual sums. Define

\[
A\coloneqq \frac{n_{2,\ell}}{\sigma_2^2}+\frac{n_{3,\ell}}{\sigma_3^2},
\qquad
B\coloneqq \frac{S_{2,\ell}}{\sigma_2^2}+\frac{S_{3,\ell}}{\sigma_3^2}.
\]

Expanding the normal likelihood and multiplying by the prior leaves the following integral after factoring out the observation-level term:

\[
M_f
=
\frac{1}{\sqrt{2\pi\lambda_f^2}}
\int
\exp\!\left[
-\frac12\left(A+\lambda_f^{-2}\right)\theta^2+B\theta
\right]d\theta.
\]

Completing the square gives a normal kernel with precision \(A+\lambda_f^{-2}\) and mean \(B/(A+\lambda_f^{-2})\). Integrating that kernel yields

\[
\boxed{
M_f
=
(1+\lambda_f^2A)^{-1/2}
\exp\!\left\{
\frac{\lambda_f^2B^2}{2(1+\lambda_f^2A)}
\right\}.
}
\]

Both control sources contribute, weighted by their residual precisions.

## 3. Derive the g-leaf factor using the same normal integral

For a \(g\) leaf, only RCT-control partial residuals contribute:

\[
r_i\mid\theta\sim N(\theta,\sigma_2^2).
\]

Accordingly, define

\[
A\coloneqq \frac{n_{2,\ell}}{\sigma_2^2},
\qquad
B\coloneqq \frac{S_{2,\ell}}{\sigma_2^2},
\]

where \(n_{2,\ell}\) and \(S_{2,\ell}\) are the count and sum of the RCT-control partial residuals in this leaf.

The leaf prior is

\[
\theta\mid z\sim
\begin{cases}
N(0,\tau_0^2),&z=1,\quad P(z=1)=w,\\
N(0,\tau_1^2),&z=0,\quad P(z=0)=1-w.
\end{cases}
\]

For each component, apply the normal integration from step 2, replacing \(\lambda_f^2\) by the corresponding component variance. Averaging the two results over \(z\) gives

\[
\boxed{
\begin{aligned}
M_g={}&
w(1+\tau_0^2A)^{-1/2}
\exp\!\left\{
\frac{\tau_0^2B^2}{2(1+\tau_0^2A)}
\right\}\\
&+(1-w)(1+\tau_1^2A)^{-1/2}
\exp\!\left\{
\frac{\tau_1^2B^2}{2(1+\tau_1^2A)}
\right\}.
\end{aligned}
}
\]

Each component integrates out \(\theta\); the weighted sum integrates out \(z\). The current values of \(w,\tau_0^2,\tau_1^2\) remain fixed during this calculation.

Thus, \(M_f\) and \(M_g\) use the same normal-integration formula, with different data contributions and leaf priors:

| Leaf | \(A\) | \(B\) | Leaf prior |
| --- | --- | --- | --- |
| \(f\) | \(n_{2,\ell}/\sigma_2^2+n_{3,\ell}/\sigma_3^2\) | \(S_{2,\ell}/\sigma_2^2+S_{3,\ell}/\sigma_3^2\) | One normal component |
| \(g\) | \(n_{2,\ell}/\sigma_2^2\) | \(S_{2,\ell}/\sigma_2^2\) | Spike-and-slab normal mixture |

The residual sums are computed from the corresponding tree's partial residuals.

For a \(g\) leaf without RCT controls, \(A=B=0\), so the formula automatically gives

\[
M_g=w+(1-w)=1.
\]

## 4. Use the leaf factors in the MH tree update

Conditional on the shared hyperparameters, the leaf priors factorize. The marginal likelihood for the whole tree is therefore

\[
p(\mathbf r\mid \mathcal T_j,\text{rest})
=
\left[\prod_i\phi(r_i;0,\sigma_i^2)\right]
\prod_{\ell\in \mathcal T_j}M_\ell.
\]

The observation-level product is identical under the current and proposed partitions: the residuals are fixed, and every included observation appears once. It therefore cancels from the likelihood ratio.

For a proposal \(\mathcal T_j\rightarrow \mathcal T_j^*\), the MH acceptance probability is

\[
\boxed{
\alpha
=
\min\!\left\{
1,\;
\frac{p(\mathcal T_j^*)}{p(\mathcal T_j)}
\frac{q(\mathcal T_j\mid \mathcal T_j^*)}{q(\mathcal T_j^*\mid \mathcal T_j)}
\frac{\prod_{\ell\in \mathcal T_j^*}M_\ell}
     {\prod_{\ell\in \mathcal T_j}M_\ell}
\right\}.
}
\]

Here, \(p(\mathcal T)\) is the tree prior and \(q\) is the proposal probability, with the fixed conditioning suppressed.

For a grow move splitting parent \(P\) into children \(L\) and \(R\), unchanged leaves cancel, leaving the likelihood factor

\[
\boxed{\frac{M_LM_R}{M_P}.}
\]

Use \(M_f\) for an \(f\)-tree update and \(M_g\) for a \(g\)-tree update.

## 5. Draw the leaf values after updating the structure

After accepting or rejecting the proposed tree:

- For an \(f\) tree, draw its leaf parameters from their normal conditional posteriors.
- For a \(g\) tree, draw each indicator from its posterior component probability, then draw its leaf parameter from the corresponding normal conditional posterior.

These draws are made for the retained tree even when the structural proposal is rejected. The sampler then proceeds to the next tree.

The leaf integrations condition on the current variances and mixture parameters; any updates to those parameters occur in separate steps. For survival outcomes, the same leaf calculations apply conditional on the current observed or imputed log event times.
