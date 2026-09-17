# PFS ESS in the real-data application

Notation: \(1\) denotes trial treatment, \(2\) trial control and \(3\) historical control. Standardization uses \(\mathbf{X}_1\) and \(\mathbf{X}_2\), with \(N\coloneqq n_1+n_2\). Each participant contributes once. In a single-arm trial, \(n_2=0\), so standardization uses only \(\mathbf{X}_1\).

Status: settled with the user, including interpretation and nonadditivity. The latest display decision replaces Figure 2 entirely with Table 7: the 22 joint covariate combinations observed in EloKRd, separate EloKRd and UCMM counts, and numerical PFS ESS for the 22 combinations represented in EloKRd. Every displayed combination has a defined empirical trial-subgroup target and numerical ESS. This supersedes the earlier dot-plot, pointwise-row and separate-panel designs.

Purpose: authoritative content reference for the ESS subsection of the future real-data application. This extends the [settled ESS definition](ess_walkthrough.md) to individual trial covariate profiles and subgroup means. The pointwise derivation remains valid background; the agreed main table targets subgroup means, as specified in section 6. The remaining application content is not yet settled. Saving this text does not recalculate study results or regenerate the manuscript table.

Approved style: a table of joint covariates, EloKRd and UCMM counts, and numerical ESS, as specified in section 6. The user requested removal of the synthetic preview files; this written specification remains the design reference. The approved bands are under 50, 50–59, 60–69, and 70+ (left-closed intervals at 50, 60 and 70); exact ages remain in predictions. The manuscript build includes Table 7.

---

1. **Start from the saved ESS definition and identify the target.**

   The saved overall target is the trial-standardized mean control outcome,

   \[
   \mu\coloneqq \frac1N\sum_{i=1}^{N}\{f(x_{1,i})+g(x_{1,i})\},
   \]

   with the general ESS definition

   \[
   \mathrm{ESS}(w,s_0^2)
   =E_{\tau_0^2}\left[
   \frac{\sigma_2^2}
   {\operatorname{Var}(\mu\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,w,\tau_0^2)}
   \right].
   \]

   For PFS, the surfaces describe mean log event time. The historical outcome records \(\mathbf{Y}_3\) include the observed PFS times and censoring indicators, handled through the survival likelihood.

   A pointwise display changes the scalar target to

   \[
   \mu(x_{1,i})\coloneqq f(x_{1,i})+g(x_{1,i}),
   \]

   the mean log-PFS time under control at trial patient \(i\)'s complete covariate profile. It retains the same preliminary historical fit, covariate conditioning, fixed reference variance, and globally calibrated PFS prior settings. Trial profiles are used because the question concerns information about control outcomes for the trial population.

2. **Apply the same variance-matching definition at a profile.**

   The pointwise extension of the saved definition is

   \[
   \boxed{
   \mathrm{ESS}(w,s_0^2;x_{1,i})
   =E_{\tau_0^2}\left[
   \frac{\sigma_2^2}
   {\operatorname{Var}(\mu(x_{1,i})\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,w,\tau_0^2)}
   \right].
   }
   \]

   As in the saved definition, this reference combines the historical-outcome posterior for \(f\) with an independent prior for \(g\). It is not the final joint posterior after trial-control outcomes have been incorporated.

   Define

   \[
   V_f(x_{1,i})\coloneqq \operatorname{Var}\{f(x_{1,i})\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2\},
   \]

   \[
   V_g(w,\tau_0^2;x_{1,i})
   \coloneqq \operatorname{Var}\{g(x_{1,i})\mid
   \mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,w,\tau_0^2\}.
   \]

   The variance defining \(V_g\) is taken under the discrepancy prior. Their sum is the variance in the denominator. At a single profile, exactly one leaf contributes from each of the \(H_g\) discrepancy trees. The contributing leaf has weight one, and all other leaves have weight zero. Thus the sum of squared weights is \(H_g\) for every tree configuration. Integrating the mixture indicators and trees gives

   \[
   V_g(w,\tau_0^2;x_{1,i})
   =H_g\{w\tau_0^2+(1-w)\tau_1^2\}.
   \]

   For the saved all-spike calibration reference, set \(w=1\), retaining the specified tree prior. Consequently,

   \[
   V_g(1,\tau_0^2;x_{1,i})=H_g\tau_0^2,
   \]

   and

   \[
   \boxed{
   \mathrm{ESS}_{\tau_0}(s_0^2;x_{1,i})
   =E_{\tau_0^2}\left[
   \frac{\sigma_2^2}{V_f(x_{1,i})+V_g(1,\tau_0^2;x_{1,i})}
   \right]
   =E_{\tau_0^2}\left[
   \frac{\sigma_2^2}{V_f(x_{1,i})+H_g\tau_0^2}
   \right].
   }
   \]

   The expectation uses the same untruncated scaled-inverse-chi-squared calibration law with parameters \(\nu_0,s_0^2\) as the saved overall ESS. Use the same globally calibrated PFS \(s_0^2\) for every profile; do not recalibrate it separately for each patient or subgroup. Estimate \(V_f(x_{1,i})\) from the variance of the preliminary posterior predictions at that profile, and average the displayed ratio over draws from the calibration scale law.

3. **Interpret one point.**

   Suppose

   \[
   \mathrm{ESS}_{\tau_0}(s_0^2;x_{1,i})=20.
   \]

   The RWD-informed prior information about the mean log-PFS time under control at this complete covariate profile is worth approximately 20 reference control observations, under the all-spike calibration reference.

   More precisely, at each fixed \(\tau_0^2\), find the number of independent normal-reference observations whose flat-prior posterior variance for the mean equals the RWD-informed prior variance of \(f(x_{1,i})+g(x_{1,i})\). Average those matched counts over the calibration distribution of \(\tau_0^2\). A value of 20 is that average; it does not assert that the fully scale-marginalized prior has variance \(\sigma_2^2/20\).

   The information comes from the full RWD dataset through the fitted model. It is not restricted to historical patients with exactly the same covariates or membership in the displayed subgroup. The units are uncensored log-PFS reference observations, not actual trial patients with the trial's censoring pattern. In this single-arm application, the reference residual variance uses the assumption \(\sigma_2^2=\sigma_3^2\).

   Higher values indicate more precise RWD-informed prior estimates of the control mean at those profiles. With the same \(H_g\), reference variance, and scale law across profiles, this variation comes through \(V_f(x_{1,i})\). It does not establish outcome compatibility between the RWD and unobserved trial controls, or quantify a realized posterior gain from borrowing.

4. **Explain the pointwise display considered earlier.**

   The earlier display calculated one pointwise ESS for each of the 30 trial profiles and displayed the same 30 values against age, sex, race, ethnicity, high-risk cytogenetics, and ASCT status in six panels. This is retained as a possible diagnostic, not the selected main Figure 2. The PFS-only scope remains agreed; omit OS points and the PFS-versus-OS legend.

   Each dot uses the patient's full observed predictor vector. For example, a dot under ASCT: Yes uses that patient's actual age, sex, race, ethnicity, and cytogenetic risk together with ASCT status. Other predictors are not replaced by averages or held at reference values.

   A panel groups or positions these profile-specific values by its horizontal-axis predictor. It does not display a separate ESS for a subgroup mean. Differences between groups may reflect their other covariates as well, and cannot be attributed solely to the displayed predictor. The vertical value for a given patient is identical across the six panels; any horizontal jitter only separates overlapping points.

5. **Distinguish pointwise, subgroup, and overall ESS.**

   For a subgroup \(G\) with \(N_G\) trial patients, its mean control log-PFS target is

   \[
   \mu_G\coloneqq \frac1{N_G}\sum_{i\in G}\{f(x_{1,i})+g(x_{1,i})\}.
   \]

   Its ESS must be calculated by substituting \(\mu_G\) into the saved variance-matching definition, retaining the same reference settings and scale law. Adding the patient-profile ESS values within \(G\) does not generally give the ESS of \(\mu_G\).

   Likewise, subgroup ESS values are not generally additive to the overall ESS. For the two disjoint ASCT groups, let \(p\) be the proportion of trial patients with ASCT. Then

   \[
   \mu=p\mu_{\mathrm{ASCT}}+(1-p)\mu_{\mathrm{no\ ASCT}}.
   \]

   At each fixed spike variance in the same calibration reference, suppressing the common conditioning only to shorten this identity,

   \[
   \begin{aligned}
   \operatorname{Var}(\mu)
   ={}&p^2\operatorname{Var}(\mu_{\mathrm{ASCT}})
   +(1-p)^2\operatorname{Var}(\mu_{\mathrm{no\ ASCT}})\\
   &+2p(1-p)\operatorname{Cov}
   (\mu_{\mathrm{ASCT}},\mu_{\mathrm{no\ ASCT}}).
   \end{aligned}
   \]

   The overall ESS is obtained by matching this overall variance and then averaging the matched counts over the scale law. Adding subgroup ESS values does not account for the weights, shared-model covariance, or the inverse-variance transformation. Thus there is no general identity equating their sum to overall ESS; equality can occur in special cases.

   Patient-profile ESS values likewise neither sum nor generally average to overall ESS. Their sum is not constrained to the overall calibration target. Calculating ESS for a subgroup or the whole trial requires the variance of that target mean, including dependence among predictions at different profiles.

6. **Report joint-subgroup ESS and observed cohort counts in Table 7.**

   The user selected exactly option 1's visual style, with option 4's subgroup-mean ESS target. Each row represents one distinct observed combination of **age band, sex, race, Hispanic/Latino status, cytogenetic risk, and ASCT status**. Its main message is: **the historical cohort supplies this much prior information for estimating the mean control outcome in each joint subgroup**.

   The RWD supplies the information. The trial covariates define the target subgroup. Include every distinct combination represented in EloKRd once. Count UCMM patients in those same combinations. Omit historical-only combinations and combinations absent from both cohorts. Patients matching on all five categorical predictors and the age band belong to the same row. Use their actual complete trial profiles when averaging the control surface, including each patient's exact age. Age bands define display groups only; they do not replace continuous age in the model or predictions, and band midpoints are not substituted for actual ages. Do not restrict the historical fit to matching historical subgroup members.

   Apply the saved all-spike calibration definition to the subgroup target from section 5:

   \[
   \boxed{
   \mathrm{ESS}_{\tau_0}(s_0^2;G)
   =E_{\tau_0^2}\left[
   \frac{\sigma_2^2}
   {\operatorname{Var}(\mu_G\mid
   \mathbf{Y}_3,\mathbf{X}_3,\mathbf{X}_1,\mathbf{X}_2,w=1,\tau_0^2)}
   \right].
   }
   \]

   Retain the same globally calibrated PFS scale law, historical fit, and fixed reference variance for all subgroups. Compute the variance of the subgroup mean, including dependence among predictions. For f, this requires taking the variance across posterior draws of the within-subgroup average. For g, it requires the subgroup leaf weights in the saved tree-prior calculation; the single-profile simplification H_g * tau_0^2 is not generally the variance of a subgroup mean.

   An ESS of 20 for a row means that the full RWD-informed prior supplies information equivalent to approximately 20 reference control observations for estimating mean control log-PFS among trial patients with that complete combination of categories, averaged over the calibration scale law. The normal-reference units, single-arm variance assumption, distinction from fully marginalized variance, and nonadditivity qualifications above continue to apply. Differences between subgroup scores are not isolated effects of any one predictor. A combination represented by one patient is still a valid row, for which the subgroup-mean ESS reduces to pointwise ESS; no minimum subgroup size has been imposed by this design choice.

   Agreed display:

   - Replace Figure 2 entirely with Table 7, PFS only; no dot plot remains.
   - Columns: Age band, Sex, Race, Hispanic/Latino, Cytogenetic risk, ASCT, EloKRd n, UCMM n, RWD prior ESS.
   - Include only the 22 combinations observed in EloKRd; omit all combinations with zero EloKRd patients. Order hierarchically by the displayed covariates: ascending age bands; Female, Male; White, Black, Other; No, Yes; Standard, High; No, Yes.
   - Age bands are under 50, 50–59, 60–69 and 70+. Exact ages remain in predictions.
   - Counts refer to observed patients within each displayed combination and sum to 30 EloKRd and 98 UCMM patients. The remaining 155 UCMM patients are omitted from the display but retained in the full historical model fit.
   - Report numerical ESS to one decimal place for the 22 combinations represented in EloKRd. Every displayed combination has trial profiles and a numerical ESS. Do not invent representative ages.
   - Each ESS uses the full historical fit and the unchanged subgroup-mean definition. It is not computed only from the historical patients counted in that row, and subgroup ESS values are not additive.
   - Retain the equivalent uncensored normal-reference units and the postinduction-ASCT qualification in the table note.

   The pointwise derivation remains explanatory background. The earlier separate-panel choice is superseded, and the synthetic preview files have been removed at the user's request. No supplementary figure has been committed by this choice.

## Implementation boundary for the revised figure

The existing revision figure computes the pointwise quantity

\[
\frac{\widehat\sigma_2^2}
{\widehat V_f(x_{1,i})+H_g\,\overline{\tau_0^2}},
\]

where the scale samples supplied to the map come from the fitted control model. This puts the mean scale inside the denominator and does not implement the settled expectation of conditional ratios under the calibration scale law. The relevant source is `revision/01-code/lrcbart/R/lrcbart.R`, function `lrc_ess_map`, called by `revision/04-application/code/run_lrcbart.R`; `make_figures.R` currently displays both PFS and OS.

The implemented main table preserves the 22 PFS subgroup-mean ESS values under the settled scale-averaging definition and adds observed counts from the harmonized analysis extract. The observed-combination table does not extend the ESS calculation to hypothetical profiles. Existing candidate-cutpoint/tree-support and single-arm fixed-scale differences remain implementation issues. No canonical study fit, analysis script or cleaned dataset was changed. See source/generated/subgroup_provenance.json for ESS provenance and source/generated/subgroup_table_validation.json for count and grid validation.
