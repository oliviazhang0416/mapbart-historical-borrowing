# LRC-BART ADDITION START
# Collect available harmonized case-study results for PFS and OS.

projectDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/lrcbart-case-study-mm"
resultDir <- file.path(projectDir, "res")
merged_file <- Sys.getenv(
  "MERGED_FILE",
  unset = file.path(projectDir, "data_cleaned", "merged_elokrd_ucmm_n283.RData")
)
data_tag <- regmatches(
  basename(merged_file), regexpr("n[0-9]+", basename(merged_file))
)
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(merged_file))

outcomes <- unique(toupper(strsplit(Sys.getenv("OUTCOMES", "PFS,OS"), "[,[:space:]]+")[[1]]))
outcomes <- outcomes[nzchar(outcomes)]
selected_configs <- unique(strsplit(Sys.getenv("LRC_CONFIGS", "default,Hf50,w0.9,Hf50_w0.9"),
                                   "[,[:space:]]+")[[1]])
selected_configs <- selected_configs[nzchar(selected_configs)]
failed_methods <- strsplit(Sys.getenv("FAILED_METHODS", ""), ",", fixed = TRUE)[[1]]
stopifnot(length(outcomes) > 0L, all(outcomes %in% c("PFS", "OS")),
          length(selected_configs) > 0L,
          all(selected_configs %in% c("default", "Hf50", "w0.9", "Hf50_w0.9")))

rows <- list()
arm_columns <- function(estimate, prefix) {
  values <- if (is.null(estimate)) rep(NA_real_, 3L) else
    as.numeric(estimate[c("est", "lo", "hi")])
  stats::setNames(as.list(values), paste0(prefix, c("_est", "_lower", "_upper")))
}

add_result <- function(method, outcome, path, coding = "harmonized",
                       config = "default", target = NA_character_,
                       prior = NA_real_) {
  script <- c(KM = "km.R", `AFT-CP` = "AFTv2.R", `Standard BART` = "BARTv2.R",
              HierAFT = "hierAFT.R", `LRC-BART` = "lrcBART.R")[[method]]
  if (paste(outcome, script, sep = ":") %in% failed_methods) return(invisible(NULL))
  if (!file.exists(path)) return(invisible(NULL))
  result <- readRDS(path)
  calibration <- result$calibration
  diagnostics <- result$diagnostics
  sampler <- result$sampler
  ucmm_est <- result$rmst_ucmm_est
  # Older KM/AFT/BART files already contain this quantity as rmst_ctrl_est.
  # Older hierarchical/LRC files need a new fit to save the external component.
  if (is.null(ucmm_est) && method %in% c("KM", "AFT-CP", "Standard BART"))
    ucmm_est <- result$rmst_ctrl_est
  ucmm_population <- result$rmst_ucmm_population
  if (is.null(ucmm_population))
    ucmm_population <- if (method == "KM") "UCMM" else "EloKRd"
  rows[[length(rows) + 1L]] <<- data.frame(
    outcome = outcome, method = method, coding = coding,
    config = config, target = target, prior = prior,
    target_N = if (!is.null(result$settings$target_N))
      as.numeric(result$settings$target_N) else NA_real_,
    tau_rmst = if (!is.null(result$tau_rmst))
      as.numeric(result$tau_rmst) else NA_real_,
    interval_type = if (method == "KM") "confidence" else "credible",
    contrast_control = if (method == "KM") "UCMM (unadjusted)" else
      "Hypothetical control (EloKRd)",
    arm_columns(result$rmst_trt_est, "rmst_trt"),
    arm_columns(if (method == "KM") NULL else result$rmst_ctrl_est,
                "rmst_hyp_ctrl"),
    arm_columns(ucmm_est, "rmst_ucmm"),
    rmst_ucmm_population = ucmm_population,
    estimate = as.numeric(result$delta_hat),
    lower = as.numeric(result$ci_95[1]),
    upper = as.numeric(result$ci_95[2]),
    difference = if (!is.null(result$delta_diff))
      as.numeric(result$delta_diff) else NA_real_,
    difference_lower = if (!is.null(result$ci_diff_95))
      as.numeric(result$ci_diff_95[1]) else NA_real_,
    difference_upper = if (!is.null(result$ci_diff_95))
      as.numeric(result$ci_diff_95[2]) else NA_real_,
    s0_sq = if (!is.null(calibration$s0_sq))
      as.numeric(calibration$s0_sq) else NA_real_,
    ess_tau0 = if (!is.null(calibration$ess_tau0))
      as.numeric(calibration$ess_tau0) else NA_real_,
    ess_realized = if (!is.null(sampler$ess_realized))
      as.numeric(sampler$ess_realized) else NA_real_,
    ess_ceiling = if (!is.null(calibration$ceiling))
      as.numeric(calibration$ceiling) else NA_real_,
    rhat = if (!is.null(diagnostics["rhat"]))
      as.numeric(diagnostics["rhat"]) else NA_real_,
    ess_bulk = if (!is.null(diagnostics["ess_bulk"]))
      as.numeric(diagnostics["ess_bulk"]) else NA_real_,
    ess_tail = if (!is.null(diagnostics["ess_tail"]))
      as.numeric(diagnostics["ess_tail"]) else NA_real_,
    file = basename(path), stringsAsFactors = FALSE
  )
  invisible(NULL)
}

for (outcome in outcomes) {
  add_result(
    "KM", outcome,
    file.path(resultDir, paste0("KM_results_", outcome, "_", data_tag,
                               ".RData")),
    coding = "none"
  )
  add_result(
    "AFT-CP", outcome,
    file.path(resultDir, paste0("AFTv2_results_", outcome, "_", data_tag,
                               "_w1.RData"))
  )
  add_result(
    "Standard BART", outcome,
    file.path(resultDir, paste0("BARTv2_results_", outcome, "_", data_tag,
                               ".RData"))
  )
  for (prior_value in c(0.05, 0.5))
    add_result(
      "HierAFT", outcome,
      file.path(resultDir, paste0("HierAFT_results_", outcome, "_", data_tag,
                                 "_prior", prior_value, ".RData")),
      config = paste0("prior", prior_value), prior = prior_value
    )

  lrc_configs <- data.frame(
    config = c("default", "Hf50", "w0.9", "Hf50_w0.9"),
    suffix = c("", "_Hf50", "_w0.9", "_Hf50_w0.9"),
    stringsAsFactors = FALSE
  )
  lrc_configs <- lrc_configs[lrc_configs$config %in% selected_configs, , drop = FALSE]
  for (config_index in seq_len(nrow(lrc_configs)))
    for (target_name in c("100", "75", "50", "f90", "f50", "f25"))
      add_result(
        "LRC-BART", outcome,
        file.path(
          resultDir,
          paste0("LRC-BART_results_", outcome, "_", data_tag,
                 lrc_configs$suffix[config_index], "_N", target_name,
                 ".RData")
        ),
        config = lrc_configs$config[config_index], target = target_name
      )
}

if (!length(rows)) stop("No completed case-study result files found")
results_table <- do.call(rbind, rows)
rownames(results_table) <- NULL
write.csv(results_table, file.path(resultDir, "results_table.csv"),
          row.names = FALSE)
saveRDS(results_table, file.path(resultDir, "results_table.RData"))

display_table <- results_table[, c(
  "outcome", "method", "config", "target", "contrast_control",
  "estimate", "lower", "upper",
  "s0_sq", "ess_tau0", "ess_realized"
)]
numeric_columns <- vapply(display_table, is.numeric, logical(1))
display_table[numeric_columns] <- lapply(
  display_table[numeric_columns], function(x) round(x, 3)
)
print(display_table, row.names = FALSE)

format_rmst <- function(prefix) {
  est <- results_table[[paste0(prefix, "_est")]]
  lo <- results_table[[paste0(prefix, "_lower")]]
  hi <- results_table[[paste0(prefix, "_upper")]]
  ifelse(is.na(est), NA_character_, sprintf("%.3f [%.3f, %.3f]", est, lo, hi))
}
arm_table <- results_table[, c(
  "outcome", "method", "config", "target", "tau_rmst", "interval_type",
  "rmst_ucmm_population"
)]
arm_table$treatment <- format_rmst("rmst_trt")
arm_table$hypothetical_control <- format_rmst("rmst_hyp_ctrl")
arm_table$ucmm_control <- format_rmst("rmst_ucmm")
cat("\nRMST in years: estimate [95% interval]\n",
    "Model UCMM controls are standardized to EloKRd covariates.\n",
    "KM reports unadjusted UCMM RMST and has no hypothetical control.\n",
    "AFT-CP and Standard BART have identical UCMM and hypothetical controls.\n",
    sep = "")
print(arm_table, row.names = FALSE)
cat("Saved", file.path(resultDir, "results_table.csv"), "\n")
# LRC-BART ADDITION END
