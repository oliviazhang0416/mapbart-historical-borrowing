# LRC-BART revision handoff — 2026-09-25

Branch `revision/lrc-bart-20260925`, based on upstream `revision/lrc-bart-20260917`
(`51b269deff1bd3dccf79b655092feaa242e8078c`). Prepared by Yunxuan.

## 1. What this branch does

It migrates the content of **PR#2** into the code.

PR#2 ("Make joint treatment modeling primary in the LRC-BART manuscript",
`koaeraser:reconcile/yunxuan-20260920` → `revision/lrc-bart-20260917`) is a
manuscript-only change: it makes joint modeling of treated trial participants,
concurrent controls and external controls the primary formulation of the paper,
but leaves the implementation on the branch untouched.

This branch carries that formulation into the code. The shared C++ model and the
Gaussian and survival two-arm joint fitting workflows are implemented. Gaussian
and survival single-arm generation, fitting, reporting and plotting are migrated
under PR#2's separate-treatment formulation, as are the case-study scripts.

Diff against the base: **66 files changed, +3274 / −1504.** No files added, none
removed. `manuscript/` is deliberately untouched — manuscript migration remains
deferred, so the paper on this branch is still the base text, not PR#2's.

Scope note: within this branch, joint fits (`joint_model = TRUE`) were run to
completion for two-arm Gaussian *and* two-arm survival — 56 complete
100-replicate result sets each, one chain, 1,000 burn-in, 1,000 retained draws.
Single-arm and the case study remain non-joint (`joint_model = FALSE`).

## 2. Joint LRC-BART performs about the same as the non-joint version

Joint LRC-BART on this branch is **nearly identical in performance** to the
previous version on `revision/lrc-bart-20260917`, which does *not* add joint
modeling of the treatment effect.

That is the central finding, and it cuts against PR#2's premise. If modelling
the treatment effect jointly leaves the operating characteristics essentially
where they were, the joint layer is not earning its place, and **the case for
adding it has not been made.** The burden here is on the joint formulation: it
adds model complexity, a second set of variance components and a slower fit, so
parity with the simpler model is an argument against it, not a neutral result.

Strength of evidence, stated honestly: this comparison is qualitative, not a
recomputed head-to-head. The baseline `LRC-BART_` results available on
`revision/lrc-bart-20260917` cover two-arm Gaussian at sc1c / sc1c-i / sc1c-ii
(six ESS targets each) and two-arm survival at sc1c only. A matched, fully
recomputed comparison over those overlapping scenarios is the obvious next step
and has not been run.

## 3. The claimed advantage over the BART counterparts is unclear

PR#2's own validation summary does not establish it. From
`joint-revision/reproducibility/joint-validation/VERIFY-OUTPUT.txt` (RMSE,
n = 40 per scenario):

| scenario | JOINT-LRC-01 | JOINT-SOURCE-01 | SEPARATED-SOURCE-01 | Δ LRC vs separated | MCSE |
|---|---|---|---|---|---|
| compatible | 0.175 | 0.213 | 0.206 | **−0.031** | 0.017 |
| heterogeneous | 0.209 | 0.211 | 0.192 | **+0.018** | 0.021 |
| null | 0.230 | 0.215 | 0.210 | **+0.020** | 0.029 |
| one | 0.207 | 0.216 | 0.211 | −0.004 | 0.018 |
| two | 0.215 | 0.215 | 0.211 | +0.004 | 0.018 |

Reading this table:

- Only the **compatible** scenario favours joint LRC-BART, and its margin
  (−0.031) is about 1.5 MCSE.
- In **heterogeneous** and **null** the joint LRC fit is *worse* than the
  separated baseline, by roughly one MCSE in each case.
- In **one** and **two** the gaps are ≤ 0.004 — far inside Monte Carlo noise.
- Coverage does not separate the methods: 188 / 200, 191 / 200 and 188 / 200
  against a nominal 190.

So the advantage is confined to the one scenario where the external controls are
compatible, and is reversed where they are not. A single favourable cell across
five scenarios, at this MCSE and n = 40, does not support a general claim of
superiority over the BART counterparts.

Two further caveats on the same evidence:

- ATE draws flagged for `R-hat > 1.05` or bulk ESS < 400 run 18 / 200, 12 / 200
  and 21 / 200 for the three methods. The rates are comparable, so flags are not
  a knock against LRC specifically — but they are not negligible for any method.
- The `chain-sensitivity/` companion, which is where the headline "LRC percentage
  reduction" figures come from, states outright that it "does not filter
  convergence flags, refit models or establish convergence." Reductions derived
  that way should not be read as established.

## 4. PR#2's evidence covers two-arm Gaussian only

PR#2's joint validation is, by its own README, a **"200-dataset joint Gaussian
validation"**: 1,800 jobs (200 datasets × 3 methods × 3 chains) at 30,000 burn-in
and 12,000 retained draws.

There is no joint validation for survival, for either single-arm setting, or for
the application. PR#2 nonetheless makes joint modeling the primary formulation of
the whole paper, whose tables draw on all four simulation projects and the
case study. **The conclusion is considerably broader than the evidence
supporting it.**

## What would settle the three points above

1. A matched joint-vs-non-joint recomputation on the overlapping scenarios
   (Gaussian sc1c / sc1c-i / sc1c-ii, survival sc1c), same data, same ESS
   targets, with convergence flags filtered rather than ignored.
2. Joint validation extended past two-arm Gaussian, or the manuscript's claims
   narrowed to the setting actually validated.
3. Re-derived comparator reductions with convergence filtering, and MCSE
   reported alongside every margin.

## Not in this branch

`res/`, `inserts/` and the simulation `data/` directories are excluded — they are
regenerable. Rebuild with each project's `data_gen_p10.R`, then `run_all.R`.
Patient-level inputs are excluded. The existing `.gitignore` already enforces
all of this; the published `merged_elokrd_ucmm_n230.RData` cohort and the
vendored `psrwe/data/` files are retained from the base branch.
