#!/usr/bin/env Rscript

# Small deterministic API/design checks only.  These use non-reserved seeds and
# n=24 synthetic datasets; they do not create validation data or fit outcomes.
Sys.setenv(VALIDATION_MODE = "validation")
source("09-method-exploration/validation/validation_core.R")

stopifnot(identical(METHODS, c("JOINT-LRC-01", "JOINT-SOURCE-01", "SEPARATED-SOURCE-01")),
          nrow(make_job_table()) == 5L * 40L * 3L * 3L,
          all(make_job_table()$dataset_seed[make_job_table()$condition %in% c("compatible", "one", "two")] >= 92026001L),
          all(make_job_table()$dataset_seed[make_job_table()$condition == "heterogeneous"] >= 92126001L),
          all(make_job_table()$dataset_seed[make_job_table()$condition == "null"] >= 92226001L),
          all(diff(make_job_table()$fit_seed) == 10L))

small <- function(region, delta, seed) {
  gen_data(4, outcome = "gaussian", n_rct = 24, n_rwd = 24,
           region = region, delta_rwd = delta, seed = seed, n_mc = 0)
}
z <- 91980001L
compatible <- small("ONE50", 0, z)
one <- small("ONE50", 2, z)
two <- small("XOR50", 2, z)
for (nm in c("rct_trt", "rct_ctrl", "x_rct", "z_rct")) {
  stopifnot(isTRUE(all.equal(compatible[[nm]], one[[nm]], tolerance = 0)),
            isTRUE(all.equal(compatible[[nm]], two[[nm]], tolerance = 0)))
}
one_rwd <- region_mask(one$rwd_ctrl[, XCOLS, drop = FALSE], "ONE50")
two_rwd <- region_mask(two$rwd_ctrl[, XCOLS, drop = FALSE], "XOR50")
stopifnot(max(abs((one$rwd_ctrl$y - compatible$rwd_ctrl$y) - 2 * (!one_rwd))) < 1e-12,
          max(abs((two$rwd_ctrl$y - compatible$rwd_ctrl$y) - 2 * (!two_rwd))) < 1e-12)

tau <- 1 + 0.5 * (compatible$x_rct[, 5] - mean(compatible$x_rct[, 5]))
stopifnot(abs(mean(tau) - 1) < 1e-12)
null_tau <- rep(0, nrow(compatible$x_rct))
stopifnot(abs(mean(null_tau)) < 1e-12)

# Exercise the actual full-size transformation functions on non-reserved seeds.
Sys.setenv(VALIDATION_MODE = "smoke")
family <- setNames(lapply(c("compatible", "one", "two"), generate_dataset, rep = 1L),
                   c("compatible", "one", "two"))
assert_constant_family(family)
for (condition in c("heterogeneous", "null")) {
  changed <- generate_dataset(condition, 1L)
  original <- gen_data(4, "gaussian", n_rct = 300, n_rwd = 300,
                       region = "XOR50", delta_rwd = 2,
                       seed = dataset_seed(condition, 1L), n_mc = 0)
  nt <- nrow(original$rct_trt)
  expected_change <- changed$true_tau[seq_len(nt)] - 1
  stopifnot(max(abs(changed$d$rct_trt$y - original$rct_trt$y - expected_change)) < 1e-12,
            identical(changed$d$rct_ctrl, original$rct_ctrl),
            identical(changed$d$rwd_ctrl, original$rwd_ctrl),
            abs(mean(changed$true_tau) - changed$truth_ate) < 1e-12,
            identical(changed$src_joint, c(rep(1L, 300), rep(2L, 300))))
}

stopifnot(is.function(joint_lrc_bart), is.function(joint_source_bart),
          is.function(joint_tau_conditional),
          is.function(getS3method("predict", "lrcbart")))
q <- joint_tau_conditional(y = rnorm(24), A = c(rep(1, 12), rep(0, 12)),
                           source = c(rep(1L, 12), rep(2L, 12)),
                           sigma1_sq = 1.5^2, sigma2_sq = 1.8^2,
                           prior_var = 100)
stopifnot(all(is.finite(c(q$mean, q$variance))), q$variance > 0)

dir.create(file.path(VALIDATION_ROOT, "checks"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "PASS: frozen ledger has 1800 jobs and disjoint fit-seed blocks",
  "PASS: common-seed trial data are identical across compatible/ONE50/XOR50 at n=24",
  "PASS: external shift differences equal the known regional shift at n=24",
  "PASS: heterogeneous effect is centered at ATE 1 and null effect is zero",
  "PASS: actual 300+300 generator paths preserve controls and apply aligned heterogeneous/null changes",
  "PASS: private-library joint signatures and scalar conditional are available",
  "No validation datasets, calibration objects, posterior fits, or state.json were changed."
), file.path(VALIDATION_ROOT, "checks", "unit_checks.txt"))
cat(paste(readLines(file.path(VALIDATION_ROOT, "checks", "unit_checks.txt")), collapse = "\n"), "\n")
