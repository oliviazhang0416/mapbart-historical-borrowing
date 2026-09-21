# Current handoff — joint-primary manuscript, 20 September 2026

Read README.md and review-record/FINAL-ASSESSMENT.md first. The current main is main_critical_final.tex/pdf, with main.tex an entry-point wrapper. The immutable upstream input is retained in the parent critical-revisor round_0 and upstream Git commit51b269deff1bd3dccf79b655092feaa242e8078c. The older upstream handoff is in content_migration/HANDOFF-UPSTREAM-20260920.md.

The primary fitted model is Y=f(x)+Sg(x)+Apsi+error: S=1 for both trial arms; A=1 only for treatment. The implemented effect is constant, psi~N(0,100), and the trial arms share residual variance. Do not present a heterogeneous treatment forest or joint survival analysis as fitted. The source discrepancy and effect are not separately identified in a single-arm trial.

Two critical review rounds, the joint primary rewrite, a house-voice pass, clean/ tracked builds, numerical checks and visual QA are complete. Main13pages, tables4pages, appendix76pages. Author field blank. All600 joint method-dataset rows and diagnostic flags retained, with chain-index sensitivity. No new MCMC fits were run. The53-file original frozen manifest verifies unchanged. The n230 application remains a separate-treatment illustration; absent raw posterior files prevent independent reproduction of that analysis.

Keep the surviving empirical and theory limitations explicit. Next work should address flagged fits, matched variance/prior choices, the empirical-Bayes treatment-omitting variance-scale rule, calibration-versus-fitting priors, null/heterogeneity assessment and a direct joint-model application. Current control-mean ESS is not treatment ESS. A fixed-region Gaussian theorem does not establish learned-BART consistency.

The current package is copied into revision/manuscript-and-slides/paper-current, RISW-2026/sources/paper/journal-revision and yunxuan-repo/manuscript/joint-revision. Root main_draft.pdf/main_tracked.pdf in manuscript-and-slides and RISW-2026 point to this revision. Original root copies were preserved under reports/pre-joint-primary-20260920. The repository branch is reconcile/yunxuan-20260920; this branch carries the coauthor-review package and targets upstream `revision/lrc-bart-20260917`. Yunxuan's original manuscript directory files remain unchanged.

The slide decks were preserved and were not rewritten during this manuscript pass. Their next scientific revision must follow the chosen joint-primary scope. Do not claim that a source-code build, checksum or completed review resolves the remaining scientific concerns.
