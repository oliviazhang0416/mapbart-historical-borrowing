# Joint-treatment manuscript assessment

The main paper has been substantially rewritten around the joint Gaussian likelihood. Yunxuan's compact 13-page main-paper structure, separate displays and detailed supplement are retained. The current output contains 13 main pages, 4 pages of main tables and a 76-page supplement. Author fields are blank.

The model, estimand, treatment-adjusted forest updates, common treatment prior, trial variance and identification boundary are consistent with the frozen joint implementation. A fixed-region joint-learning theorem and proof have been added; this is not a learned-BART consistency theorem. The control-mean ESS reference is distinguished from treatment uncertainty, with a Schur-complement calculation and marginal-ELIR discussion.

The five-condition 200-dataset Gaussian experiment is primary. All 600 method-dataset rows and diagnostic flags remain. Direct recomputation verified all15 RMSE and coverage rows plus flag counts and Monte Carlo errors. Saved-chain RMSE sensitivity retains all120 Compatible/One-variable/Two-variable datasets for each chain index. It supports the observed pooled direction but does not establish convergence. The original53 frozen code/output hashes are unchanged. No new MCMC fit was run.

The earlier100-replicate simulations and n230 single-arm survival illustration remain clearly separate-treatment variants. Stale sensitivity ranges and ESS summaries were reconciled against n230 producing summaries. The main application numbers have not been independently refitted because58 underlying posterior output files are absent. Old n283 HR outputs were not substituted.

Two rounds of internal LLM-assisted technical review are complete, with dispositions in each round ledger. Both second-round reviewers retained substantive concerns about convergence and unmatched comparator choices. Their requested bounded presentation remedies and saved-chain sensitivity were incorporated; the evidence limitations remain open. These internal editing audits are not coauthor review or journal peer review and do not establish journal readiness.

Before top-statistics-journal submission:

1. Fix the estimand and model configuration, resolve flagged ATE/discrepancy fits under a uniform rule and rerun the locked comparison if needed.
2. Separate architectural effects from variance-sharing and prior-range choices through matched sensitivities, including trial residual-scale estimation with treatment included.
3. Align calibration and fitting tree/scale priors or explicitly validate the current scale-selection reference; do not interpret a control target of100 as treatment ESS.
4. Expand null and heterogeneous-effect assessment with Monte Carlo uncertainty and a prespecified decision rule if decision claims are desired.
5. Supply an application that directly assesses the joint model and complete a focused prior-art comparison. The current fixed-region theory and single-arm illustration cannot stand in for these.

Build and integrity: all3 clean PDFs compile with recorder and zero undefined-reference/citation, overfull, duplicate-label or duplicate-destination warnings. Main, table and appendix redlines also compile with no such warnings; table/appendix PDFs use additions-only color fallback, with complete unified patches. Pure prose pass removed12 redundant bold lead-ins and passed protected mathematical, numerical, citation and structural checks. Visual QA notes and full-page renders are retained in the revision workspace.
