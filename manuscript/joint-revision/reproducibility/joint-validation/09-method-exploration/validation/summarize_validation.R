#!/usr/bin/env Rscript

source("09-method-exploration/validation/validation_core.R")
suppressPackageStartupMessages(library(posterior))

mode <- tolower(Sys.getenv("VALIDATION_MODE", "validation"))
if (!mode %in% c("validation", "smoke")) stop("VALIDATION_MODE must be validation or smoke")
mode_dir <- file.path(VALIDATION_ROOT, if (mode == "smoke") "smoke" else "results")
jobs_path <- file.path(mode_dir, "jobs_selected.csv")
if (!file.exists(jobs_path)) stop("missing selected jobs: ", jobs_path)
jobs <- read.csv(jobs_path, stringsAsFactors = FALSE)
if (!all(c(1, 2, 3) %in% jobs$chain)) stop("canonical summaries require all three chains")

raw_path <- function(j) file.path(mode_dir, "raw",
  sprintf("%s_r%02d_c%d_%s.rds", j$condition, j$rep, j$chain, j$method))
load_raw <- function(condition, rep, chain, method) {
  p <- raw_path(data.frame(condition = condition, rep = rep, chain = chain, method = method))
  if (!file.exists(p)) stop("missing raw result: ", p)
  z <- readRDS(p)
  expected_job <- jobs[jobs$condition == condition & jobs$rep == rep &
                       jobs$chain == chain & jobs$method == method, , drop = FALSE]
  cached <- load_dataset(mode_dir, condition, rep)
  stopifnot(nrow(expected_job) == 1L,
            identical(z$method, method), identical(z$condition, condition),
            z$rep == rep, z$chain == chain,
            z$dataset_seed == expected_job$dataset_seed_effective,
            z$seed == expected_job$fit_seed_effective,
            identical(z$data_sha256, cached$data_sha256),
            identical(z$code_hash, CODE_HASH),
            identical(z$library_hash, LIBRARY_HASH),
            identical(z$config_hash, method_config_hash(method, z$n_burn, z$n_draw)),
            length(z$ate) == z$n_draw, length(z$g_contrast) == z$n_draw,
            z$interval_level == .95, z$variance_model == "source-specific")
  if (mode == "validation") stopifnot(z$n_burn == 30000L, z$n_draw == 12000L)
  z
}

groups <- unique(jobs[c("condition", "rep", "method")])
metric_rows <- list(); canonical_rows <- list(); conv_rows <- list()
for (i in seq_len(nrow(groups))) {
  g <- groups[i, , drop = FALSE]
  ch <- sort(jobs$chain[jobs$condition == g$condition & jobs$rep == g$rep & jobs$method == g$method])
  if (!identical(ch, 1:3)) stop("not exactly three chains for ", paste(g, collapse = "/"))
  pp <- lapply(ch, function(k) load_raw(g$condition, g$rep, k, g$method))
  draw_counts <- vapply(pp, function(z) length(z$ate), integer(1))
  stopifnot(length(unique(draw_counts)) == 1L)
  n_draw <- draw_counts[1]
  for (field in c("data_sha256", "truth", "truth_control", "true_g",
                  "dataset_seed", "batch_started_utc", "preselection_utc",
                  "primary_candidate_id", "primary_baseline_id", "config_hash"))
    stopifnot(identical(pp[[1]][[field]], pp[[2]][[field]]),
              identical(pp[[1]][[field]], pp[[3]][[field]]))
  arr <- array(NA_real_, dim = c(n_draw, 3, 2),
               dimnames = list(iteration = NULL, chain = paste0("c", 1:3),
                               variable = c("ate", "g_contrast")))
  for (k in 1:3) {
    arr[, k, 1] <- pp[[k]]$ate[seq_len(n_draw)]
    arr[, k, 2] <- pp[[k]]$g_contrast[seq_len(n_draw)]
  }
  cv <- as.data.frame(summarise_draws(as_draws_array(arr)))
  cv$method <- g$method; cv$scenario <- g$condition; cv$rep <- g$rep
  conv_rows[[length(conv_rows) + 1L]] <- cv

  ate <- as.numeric(arr[, , "ate"])
  lower <- unname(quantile(ate, .025)); upper <- unname(quantile(ate, .975))
  truth <- pp[[1]]$truth
  ctrl_surface <- Reduce("+", lapply(pp, `[[`, "control_surface")) / 3
  truth_ctrl <- pp[[1]]$truth_control
  gmean <- Reduce("+", lapply(pp, `[[`, "g_mean")) / 3
  gtrue <- pp[[1]]$true_g
  arow <- data.frame(method = g$method, scenario = g$condition, rep = g$rep,
                     ate_mean = mean(ate), truth = truth, lower = lower, upper = upper,
                     max_rhat = cv$rhat[cv$variable == "ate"],
                     min_bulk_ess = cv$ess_bulk[cv$variable == "ate"],
                     dataset_seed = pp[[1]]$dataset_seed,
                     fit_seed_c1 = pp[[1]]$seed, fit_seed_c2 = pp[[2]]$seed,
                     fit_seed_c3 = pp[[3]]$seed,
                     batch_started_utc = pp[[1]]$batch_started_utc,
                     data_sha256 = pp[[1]]$data_sha256,
                     data_hash = pp[[1]]$data_sha256,
                     config_hash = pp[[1]]$config_hash,
                     code_hash = pp[[1]]$code_hash,
                     library_hash = pp[[1]]$library_hash,
                     variance_model = pp[[1]]$variance_model,
                     interval_level = pp[[1]]$interval_level,
                     preselection_utc = pp[[1]]$preselection_utc,
                     primary_candidate_id = pp[[1]]$primary_candidate_id,
                     primary_baseline_id = pp[[1]]$primary_baseline_id,
                     source1_controls_n = pp[[1]]$source_counts[1],
                     source2_external_n = pp[[1]]$source_counts[2],
                     source1_joint_n = pp[[1]]$source_counts_joint[1],
                     source2_joint_n = pp[[1]]$source_counts_joint[2],
                     stringsAsFactors = FALSE)
  canonical_rows[[length(canonical_rows) + 1L]] <- arow
  metric_rows[[length(metric_rows) + 1L]] <- data.frame(
    method = g$method, scenario = g$condition, rep = g$rep, n_draw = n_draw,
    ate_mean = mean(ate), truth = truth, lower = lower, upper = upper,
    ate_bias = mean(ate) - truth, ate_rmse_draws = sqrt(mean((ate - truth)^2)),
    covered = lower <= truth && upper >= truth, width = upper - lower,
    trial_control_rmse = sqrt(mean((ctrl_surface - truth_ctrl)^2)),
    g_rmse = sqrt(mean((gmean - gtrue)^2)),
    g_contrast_mean = mean(vapply(pp, function(z) mean(z$g_contrast), numeric(1))),
    ate_rhat = cv$rhat[cv$variable == "ate"],
    ate_bulk_ess = cv$ess_bulk[cv$variable == "ate"],
    g_rhat = cv$rhat[cv$variable == "g_contrast"],
    g_bulk_ess = cv$ess_bulk[cv$variable == "g_contrast"],
    source1_controls_n = pp[[1]]$source_counts[1], source2_external_n = pp[[1]]$source_counts[2],
    source1_joint_n = pp[[1]]$source_counts_joint[1], source2_joint_n = pp[[1]]$source_counts_joint[2],
    data_sha256 = pp[[1]]$data_sha256, data_hash = pp[[1]]$data_sha256,
    config_hash = pp[[1]]$config_hash, code_hash = pp[[1]]$code_hash,
    library_hash = pp[[1]]$library_hash, stringsAsFactors = FALSE)
}

canonical <- do.call(rbind, canonical_rows)
metrics <- do.call(rbind, metric_rows)
conv <- do.call(rbind, conv_rows)
write.csv(canonical, file.path(mode_dir, "canonical.csv"), row.names = FALSE)
write.csv(metrics, file.path(mode_dir, "per_dataset.csv"), row.names = FALSE)
write.csv(conv, file.path(mode_dir, "convergence.csv"), row.names = FALSE)

summary_rows <- lapply(split(metrics, list(metrics$method, metrics$scenario), drop = TRUE), function(x) {
  data.frame(method = x$method[1], scenario = x$scenario[1], n = nrow(x),
             ate_bias = mean(x$ate_bias), ate_rmse = sqrt(mean(x$ate_bias^2)),
             coverage = mean(x$covered), coverage_count = sum(x$covered),
             width = mean(x$width), trial_control_rmse = mean(x$trial_control_rmse),
             g_rmse = mean(x$g_rmse), max_ate_rhat = max(x$ate_rhat),
             min_ate_bulk_ess = min(x$ate_bulk_ess), max_g_rhat = max(x$g_rhat),
             min_g_bulk_ess = min(x$g_bulk_ess), stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, summary_rows), file.path(mode_dir, "summary.csv"), row.names = FALSE)

conv_summary <- do.call(rbind, lapply(split(conv, list(conv$method, conv$variable), drop = TRUE), function(x) {
  data.frame(method = x$method[1], variable = x$variable[1], n = nrow(x),
             rhat_median = median(x$rhat, na.rm = TRUE), rhat_max = max(x$rhat, na.rm = TRUE),
             bulk_ess_min = min(x$ess_bulk, na.rm = TRUE),
             bulk_ess_median = median(x$ess_bulk, na.rm = TRUE), stringsAsFactors = FALSE)
}))
write.csv(conv_summary, file.path(mode_dir, "convergence_summary.csv"), row.names = FALSE)

writeLines(c(
  paste0("Validation summary mode=", mode),
  paste0("Canonical rows: ", nrow(canonical), "; expected methods: ", paste(METHODS, collapse = ", ")),
  "canonical.csv is directly consumable by 09-method-exploration/evaluate.R.",
  "ATE convergence is reported in canonical.csv; discrepancy convergence is retained separately in convergence.csv.",
  "No promising designation is assigned by this script; the primary agent applies the preregistered gates."
), file.path(mode_dir, "REPORT.md"))
print(do.call(rbind, summary_rows))
