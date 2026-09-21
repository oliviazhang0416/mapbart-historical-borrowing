# Joint-treatment manuscript revision — 20 September 2026

The current main paper is `main_critical_final.tex` / `main_critical_final.pdf` (13 pages). `main.tex` is an entry-point wrapper; `main.pdf`, `main_draft.tex` and `main_draft.pdf` are convenience copies. Author names and affiliations remain blank.

The companion `main_tables_figures.pdf` contains four main tables. `appendix.pdf` contains the computational, information, theoretical and empirical details, including clearly separate treatment-model variants.

## Build

From this directory, run `python3 source/scripts/build_review_package.py`. From elsewhere, use the script’s full path. It builds all three documents with recorder output, resolves external references through repeated passes, retains logs under `build/`, and reports undefined references/citations and layout warnings. TeX Live and Python 3 are required. The older inherited `build_pdfs.py` is retained for lineage; use `build_review_package.py` for this revision.

## Changes from Yunxuan's revision

The primary likelihood is now `Y = f(x) + S g(x) + A psi + error`, where both trial arms have S=1 and only treated trial participants have A=1. The paper defines the common treatment effect, treatment-adjusted tree updates, shared trial residual variance, and the lack of joint treatment/source identification without concurrent controls. It includes a fixed-region joint-learning theorem and proof, and distinguishes the control-mean prior information reference from uncertainty about psi.

The 200-dataset joint Gaussian study is the primary numerical evidence. Its condition-specific results, comparator differences, convergence flags, paired comparisons and saved-chain sensitivity are reported. The earlier 100-replicate experiments and the n230 single-arm survival application remain separate-treatment analyses. Their results are retained with their own likelihoods and provenance. The existing survival illustration is not a joint-model application.

## Review and tracked changes

`review-record/` contains two rounds of internal LLM-assisted technical review and their resolution ledgers. These are editing audits, not coauthor review or journal peer review. The immutable comparison baseline is Yunxuan's upstream commit 51b269deff1bd3dccf79b655092feaa242e8078c, preserved in Git history and in the unchanged parent `manuscript/` files. No original numerical fit was overwritten.

`main_tracked.pdf`, `main_tables_figures_tracked.pdf` and `appendix_tracked.pdf` compare with that baseline. The main uses colored additions/deletions. Tabular deletion markup broke alignment in the table and appendix redlines, so those PDFs use additions-only color markup. Full deletions and equation changes remain in all three `*_tracked.patch` files. Redline PDFs are reading aids; the clean PDFs and exact patches are authoritative. `tracked-source/` retains the reference/assets needed to compile the redlines.

## Reproducibility and submission status

`reproducibility/joint-validation/` contains source-only R/C++ code, the frozen job ledger, summaries for all 600 method-dataset rows, chain-index sensitivity, checksums and a no-fit numerical check. Run `python3 reproducibility/joint-validation/verify_joint_validation.py`. The 1,800 original raw chain files remain in the parent revision workspace at `09-method-exploration/validation/results/raw/`; they are not redistributed here. Original application posterior result files are unavailable in the upstream checkout; application numbers were reconciled against its saved summaries and sources, without a new fit.

This is a manuscript for coauthor review, not a submission-ready claim. Before a top-statistics-journal submission, resolve flagged treatment and discrepancy fits, assess matched variance/prior specifications (including a treatment-adjusted variance prior scale), align or validate the calibration reference against the fitting prior, strengthen null and heterogeneity assessment, and add an application that directly evaluates the joint randomized-trial model. No general superiority, joint-survival validation or treatment-ESS interpretation is claimed.

## Coauthor review

Please start with [REVIEW-FOR-YUNXUAN.md](REVIEW-FOR-YUNXUAN.md). It links the six clean/tracked PDFs and lists the scientific decisions we would like to settle before further experiments. The frozen simulation source is unchanged; the build script additionally locates TeX programs on PATH for portability.
