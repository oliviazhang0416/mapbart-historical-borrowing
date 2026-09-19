# Patient-profile prior ESS in the multiple-myeloma application

## Purpose and scope

Table 7 reports progression-free-survival (PFS) prior effective sample size for each of the 30 EloKRd patient covariate profiles. The information comes from the full 200-patient triplet UCMM cohort through the preliminary historical model. Overall survival remains in the outcome tables and is not included in this local-information table.

The six displayed covariate dimensions are continuous age, sex, race, Hispanic/Latino status, cytogenetic risk and ASCT. Race is represented in the model by Black and Other indicators with White as reference, so the six dimensions correspond to seven numeric predictors. ASCT is postinduction.

## 1. Profile-specific target

For EloKRd patient \(i\), let \(x_i\) be the complete model covariate vector and define the conditional mean control log-PFS

\[
\mu_i \coloneqq f(x_i)+g(x_i).
\]

The preliminary historical model is fitted to all UCMM outcomes and evaluated at \(x_i\). Its posterior uncertainty at that profile is

\[
V_f(x_i)
\coloneqq
\operatorname{Var}\{f(x_i)\mid
\mathbf Y_2,\mathbf X_2,x_i\}.
\]

The 4,000 retained preliminary posterior predictions give the sample estimate of \(V_f(x_i)\). The historical data are not filtered to patients whose covariates exactly match \(x_i\).

## 2. Discrepancy variance at one profile

At a single profile, exactly one terminal node contributes from each of the \(H_g\) discrepancy trees. The contributing node has weight one. Therefore

\[
V_g(w,\tau_0^2;x_i)
=H_g\{w\tau_0^2+(1-w)\tau_1^2\}.
\]

The calibration reference sets \(w=1\), giving

\[
V_g(1,\tau_0^2;x_i)=H_g\tau_0^2.
\]

This single-profile identity does not require prior-tree simulation because the sum of squared profile weights is exactly one in every tree.

## 3. Matching reference and prior ESS

The reference experiment for the same patient profile is

\[
Y_{ij}^{\mathrm{ref}}\mid x_i
\sim N(\mu_i,\sigma_1^2).
\]

The target prior and the reference experiment therefore concern the same scalar parameter, \(\mu_i\). One reference observation supplies Fisher information \(1/\sigma_1^2\). At fixed \(\tau_0^2\), matching that reference information to the inverse RWD-informed prior variance gives

\[
m_i(\tau_0^2)
=
\frac{\sigma_1^2}
{V_f(x_i)+H_g\tau_0^2}.
\]

The patient-profile ESS averages these matched counts over the calibration distribution:

\[
\boxed{
\operatorname{ESS}_{\tau_0}(s_0^2;x_i)
=E_{\tau_0^2}\left[
\frac{\sigma_1^2}
{V_f(x_i)+H_g\tau_0^2}
\right],
\qquad
\tau_0^2\sim
\operatorname{scaled\text{-}Inv}\chi^2(\nu_0,s_0^2).
}
\]

The application uses \(H_f=10\), \(H_g=5\), \(\nu_0=3\), \(s_0^2=0.00398107\), and the saved PFS residual variance \(\sigma_1^2=1.91808\). The single-arm reference uses the UCMM residual variance under \(\sigma_1^2=\sigma_2^2\). Every profile uses the same \(\sigma_1^2\), \(H_g\), and scale law. Differences across rows arise from \(V_f(x_i)\).

## 4. Interpretation

If a row has prior ESS 20, the UCMM-informed prior supplies, on average over the spike-variance calibration law, the same conditional precision about mean control log-PFS at that covariate profile as approximately 20 independent normal-reference observations with that same profile.

The unit is one uncensored normal-reference observation on the log-PFS scale. It is not an equivalent number of patients under the trial's censoring distribution. The value describes prior precision and does not establish compatibility between UCMM outcomes and unobserved EloKRd control outcomes. It also does not quantify a realized posterior gain from borrowing.

The full UCMM fit informs every row. The displayed UCMM count summarizes how many historical patients match the row's age group and categorical covariates, but the ESS is not computed only from those patients. A high ESS therefore need not track the displayed count.

## 5. Relationship to the overall ESS

The standardized overall current-control target is

\[
\bar\mu=\frac1{30}\sum_{i=1}^{30}\mu_i.
\]

Its variance includes covariance among predictions:

\[
\operatorname{Var}(\bar\mu)
=\frac1{30^2}\left[
\sum_i\operatorname{Var}(\mu_i)
+\sum_{i\ne i'}\operatorname{Cov}(\mu_i,\mu_{i'})
\right].
\]

The same fitted functions inform all profiles, so these covariance terms generally do not vanish. The variance-to-ESS transformation is also nonlinear. The 30 pointwise ESS values therefore do not sum or generally average to the overall ESS. The overall ESS must be calculated from the variance of \(\bar\mu\).

## 6. Table 7 specification

Table 7 contains 30 rows ordered by age group and exact age. Age group and UCMM cohort count are displayed on every row. Its columns are:

1. age group;
2. age;
3. sex;
4. race;
5. Hispanic/Latino status;
6. cytogenetic risk;
7. ASCT;
8. UCMM cohort count for the displayed covariate combination;
9. prior ESS.

The displayed categories are under 50, 50--59, 60--69 and 70 years or older. A row's UCMM count is the number of triplet-treated historical controls matching its age group, sex, race, Hispanic/Latino status, cytogenetic risk and ASCT. Counts repeat when EloKRd rows share a combination and are not additive down the table. Age is displayed to one decimal year, while calculations use exact age. Prior ESS is displayed to one decimal place. The observed pointwise values range from 3.8 to 24.5. The table note states that the combination counts are descriptive, the full historical cohort informs each profile, and the common calibration settings, reference unit, nonaggregation property, compatibility limitation and postinduction status of ASCT continue to apply.

## 7. Reproducible calculation

The saved calibration object contains the 30 values of \(V_f(x_i)\) in the same order as the EloKRd rows in the merged analysis extract. For each profile, draw

\[
U_b\sim\chi^2_3,
\qquad
\tau_{0b}^2=\frac{3s_0^2}{U_b},
\]

and compute

\[
\widehat{\operatorname{ESS}}_i
=\frac1B\sum_{b=1}^B
\frac{\widehat\sigma_1^2}
{\widehat V_f(x_i)+H_g\tau_{0b}^2}.
\]

The manuscript build uses \(B=100{,}000\) common scale draws for all 30 profiles. Reusing the draws makes cross-profile comparisons depend only on the saved profile-specific historical uncertainty.
