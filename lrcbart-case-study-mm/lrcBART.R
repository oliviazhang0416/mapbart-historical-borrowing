rm(list = ls())


# Resolve the repository root: walk up from this script's own location first,
# then from the working directory. Works under Rscript, source(), and R CMD.
.lrcRoot <- local({
  .up <- function(d) {
    d <- tryCatch(normalizePath(d, winslash = "/", mustWork = TRUE),
                  error = function(e) NA_character_)
    if (is.na(d)) return(NA_character_)
    while (!file.exists(file.path(d, ".lrcbart-root")) && dirname(d) != d) d <- dirname(d)
    if (file.exists(file.path(d, ".lrcbart-root"))) d else NA_character_
  }
  cand <- character(0)
  for (i in seq_len(sys.nframe())) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(of) && nzchar(of)) cand <- c(cand, dirname(of))
  }
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) cand <- c(cand, dirname(sub("^--file=", "", m[1])))
  cand <- c(cand, getwd())
  hit <- NA_character_
  for (p in cand) if (is.na(hit)) hit <- .up(p)
  if (is.na(hit))
    stop("lrcbart repository root (.lrcbart-root marker) not found from: ",
         paste(unique(cand), collapse = ", "))
  hit
})
# LRC-BART MODIFICATION START
# This file began as the case study's original MAP-BART method file. It retains the
# censored-AFT treatment fit, five-year RMST estimand, result structure, and
# user-facing script role while replacing the MAP-BART control fit with the
# approved single-arm LRC-BART f(X)+g(X) fit.
mainDir <- .lrcRoot
projDir <- file.path(
  mainDir, "lrcbart-case-study-mm"
)

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppEigen)
  library(survival)
})

OUTCOME <- Sys.getenv("OUTCOME", unset = "PFS")
stopifnot(OUTCOME %in% c("PFS", "OS"))
merged_file <- Sys.getenv(
  "MERGED_FILE",
  unset = file.path(projDir, "data_cleaned", "merged_elokrd_ucmm_n283.RData")
)
if (!file.exists(merged_file) &&
    file.exists(file.path(projDir, "data_cleaned", basename(merged_file))))
  merged_file <- file.path(projDir, "data_cleaned", basename(merged_file))
if (!file.exists(merged_file))
  stop("Merged file not found: ", merged_file,
       " (run private_data/data_merge.R first)")
load(merged_file)
merged_file <- normalizePath(merged_file)
data_tag <- regmatches(
  basename(merged_file), regexpr("n[0-9]+", basename(merged_file))
)
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(merged_file))
harmonized_covariates <- c(
  "age", "male", "race_Black", "race_Other", "hispanic",
  "high_risk_cyto", "asct"
)
if (!all(harmonized_covariates %in% names(merged)))
  stop("Merged data lacks harmonized columns; run private_data/data_merge.R")
X_all <- as.matrix(merged[, harmonized_covariates, drop = FALSE])
storage.mode(X_all) <- "double"
Y <- if (OUTCOME == "PFS") merged$pfs_months / 12 else merged$os_months / 12
event <- as.integer(if (OUTCOME == "PFS") merged$pfs_status else merged$os_status)
trt <- as.integer(merged$trt)
trt_idx <- which(trt == 1L)
ctrl_idx <- which(trt == 0L)

cat(sprintf("=== lrcBART.R: OUTCOME = %s; coding = harmonized ===\n", OUTCOME))
cat(sprintf("Merged data: %d subjects (%d EloKRd, %d UCMM)\n",
            length(trt), length(trt_idx), length(ctrl_idx)))
cat("Columns:", paste(colnames(X_all), collapse = ", "), "\n")

# The revision application fits four chains with 2,000 warm-up and 2,000
# retained draws. Stage 1 retains 4,000 draws after the same warm-up.
n_chains <- 4L
n_burn <- 2000L
n_draw <- 2000L
stage1_draw <- 4000L
H_g <- 5L
tau_rmst <- 5

# One method file owns the reported configuration and three approved
# sensitivities. The reported H_f=10, w=1 setting is called `default`.
configuration_table <- data.frame(
  config = c("default", "Hf50", "w0.9", "Hf50_w0.9"),
  H_f = c(10L, 50L, 10L, 50L),
  w = c(1, 1, 0.9, 0.9),
  stringsAsFactors = FALSE
)
requested_configs <- strsplit(
  Sys.getenv("LRC_CONFIGS", unset = paste(configuration_table$config,
                                           collapse = ",")),
  "[, ]+"
)[[1]]
requested_configs <- requested_configs[nzchar(requested_configs)]
bad_configs <- setdiff(requested_configs, configuration_table$config)
if (length(bad_configs))
  stop("Unknown LRC_CONFIGS value(s): ", paste(bad_configs, collapse = ", "))
configuration_table <- configuration_table[
  match(requested_configs, configuration_table$config), , drop = FALSE
]

target_names <- c("100", "75", "50", "f90", "f50", "f25")
target_labels <- stats::setNames(
  c("100% of EloKRd n (N=30)", "75% of EloKRd n (N=22.5)",
    "50% of EloKRd n (N=15)", "90% ceiling", "50% ceiling",
    "25% ceiling"),
  target_names
)

# Load the active standalone sources directly. No C++ copy is kept here.
lrc_cpp_cache <- file.path(tempdir(), "sourceCpp_case_study_lrcBART")
dir.create(lrc_cpp_cache, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file.path(mainDir, "lrcBART", "clrcbart.cpp"),
                cacheDir = lrc_cpp_cache)
Rcpp::sourceCpp(file.path(mainDir, "lrcBART", "cess.cpp"),
                cacheDir = lrc_cpp_cache)

treatment_cpp_cache <- file.path(tempdir(), "sourceCpp_case_study_treatment")
dir.create(treatment_cpp_cache, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file.path(mainDir, "aBART", "cabart.cpp"),
                cacheDir = treatment_cpp_cache)
source(file.path(mainDir, "bartModelMatrix.R"), local = TRUE)

resultDir <- file.path(projDir, "res")
cacheDir <- file.path(resultDir, "cache")
essDir <- file.path(projDir, "res", "ess")
dir.create(cacheDir, recursive = TRUE, showWarnings = FALSE)
dir.create(essDir, recursive = TRUE, showWarnings = FALSE)

# #------------------------------------------------
# #------------ ESS Curve Construction ------------
# #------------------------------------------------
# Construct or load one ESS curve for each requested H_f. ess_cal.R contains
# the explicit s_0^2 range and returns the six selected target rows.
calibrations <- list()
for (H_f_value in unique(configuration_table$H_f)) {
  ess_data_file <- merged_file
  ess_outcome <- OUTCOME
  ess_data_tag <- data_tag
  ess_X_all <- X_all
  ess_time <- Y
  ess_status <- event
  ess_trt <- trt
  ess_H_g <- H_g
  ess_H_f <- H_f_value
  ess_n_burn <- n_burn
  ess_n_draw <- stage1_draw
  ess_seed <- 1L
  ess_checkpoint_dir <- essDir
  source(file.path(projDir, "ess_local", "ess_cal.R"), local = TRUE)
  calibrations[[as.character(H_f_value)]] <- ess_calibration
}

# #------------------------------------------------
# #---------------- RCT Treatment -----------------
# #------------------------------------------------
# Fit the inherited 50-tree EloKRd treatment AFT-BART once per chain and cache
# it for all LRC-BART configurations and targets.
treatment_settings <- list(
  H_f = 50L, n_chains = n_chains, n_burn = n_burn,
  n_draw = n_draw, coding = "harmonized", outcome = OUTCOME
)
treatment_cache_file <- file.path(
  cacheDir, paste0("lrcbart_treatment_", tolower(OUTCOME), "_", data_tag,
                   ".RData")
)
treatment_cache_valid <- FALSE
if (file.exists(treatment_cache_file) &&
    file.info(treatment_cache_file)$mtime >= file.info(merged_file)$mtime) {
  treatment_cached <- readRDS(treatment_cache_file)
  treatment_cache_valid <- identical(treatment_cached$settings,
                                     treatment_settings) &&
    identical(treatment_cached$data_file, normalizePath(merged_file))
}

fit_treatment_chain <- function(chain_id) {
  x_train_frame <- data.frame(X = X_all[trt_idx, , drop = FALSE])
  x_test_frame <- data.frame(X = X_all[trt_idx, , drop = FALSE])
  y_train <- log(Y[trt_idx])
  event_train <- as.integer(event[trt_idx])
  numcut <- 100L
  temp <- bartModelMatrix(
    x_train_frame, numcut, usequants = FALSE,
    xinfo = matrix(0, 0, 0), rm.const = TRUE
  )
  x_train <- t(temp$X)
  x_test <- bartModelMatrix(x_test_frame)
  x_test <- t(x_test[, temp$rm.const, drop = FALSE])
  p <- nrow(x_train)
  n <- ncol(x_train)
  np <- ncol(x_test)
  offset <- mean(y_train)
  y_centered <- y_train - offset
  aft <- survival::survreg(
    survival::Surv(exp(y_centered), event_train) ~ .,
    data = data.frame(t(x_train)), dist = "lognormal"
  )
  sigest <- aft$scale
  sigdf <- 3
  lambda <- sigest^2 * stats::qchisq(0.10, sigdf) / sigdf
  tau <- diff(range(y_centered)) / (4 * sqrt(50L))
  grp <- temp$grp
  if (is.null(grp)) grp <- seq_len(p)
  seed <- 1000L + chain_id
  fit <- cabart(
    1L, n, p, np, x_train, y_centered, event_train, x_test,
    50L, temp$numcut, n_draw, n_burn, 1L,
    2, 0.95, offset, tau, sigdf, lambda, sigest,
    rep(1, n), FALSE, 0, 1, grp, 0.5, 1, p, FALSE,
    10000L, temp$xinfo, seed
  )
  if (n_burn > 0L)
    fit$sigma <- fit$sigma[-seq_len(n_burn)]
  list(mu = fit$yhat.test, sigma = fit$sigma)
}

if (treatment_cache_valid) {
  treatment_chains <- treatment_cached$fit
} else {
  treatment_chains <- lapply(seq_len(n_chains), fit_treatment_chain)
  saveRDS(
    list(fit = treatment_chains, settings = treatment_settings,
         data_file = normalizePath(merged_file)),
    treatment_cache_file
  )
}

compute_pop_rmst <- function(mu_draws, sigma_draws, tau) {
  lt <- log(tau)
  output <- numeric(nrow(mu_draws))
  for (draw_id in seq_len(nrow(mu_draws))) {
    mu <- mu_draws[draw_id, ]
    sigma <- sigma_draws[draw_id]
    output[draw_id] <- mean(
      exp(mu + sigma^2 / 2) *
        pnorm((lt - mu - sigma^2) / sigma) +
        tau * pnorm((mu - lt) / sigma)
    )
  }
  output
}

arm_summary <- function(values) {
  interval <- quantile(values, c(0.025, 0.975), na.rm = TRUE)
  c(est = median(values, na.rm = TRUE),
    lo = unname(interval[1]), hi = unname(interval[2]))
}

treatment_rmst <- lapply(
  treatment_chains,
  function(fit) compute_pop_rmst(fit$mu, fit$sigma, tau_rmst)
)

# #------------------------------------------------
# #------------- ESS Target Fits ------------------
# #------------------------------------------------
for (configuration_index in seq_len(nrow(configuration_table))) {
  configuration <- configuration_table[configuration_index, ]
  calibration <- calibrations[[as.character(configuration$H_f)]]

  for (target_name in target_names) {
    target_row <- calibration$targets[
      calibration$targets$target_name == target_name, , drop = FALSE
    ]
    if (nrow(target_row) != 1L)
      stop("Missing ESS target ", target_name)

    cat(sprintf(
      "\n%s %s: H_f=%d, H_g=%d, w=%g, %s, s0^2=%g\n",
      OUTCOME, configuration$config, configuration$H_f, H_g,
      configuration$w, target_labels[target_name], target_row$s0_sq
    ))

    control_rmst <- vector("list", n_chains)
    ucmm_rmst <- vector("list", n_chains)
    ratio_chains <- matrix(NA_real_, n_draw, n_chains)
    difference_chains <- matrix(NA_real_, n_draw, n_chains)
    sigma_control_chains <- matrix(NA_real_, n_draw, n_chains)
    g_draws_all <- NULL
    profile_ess_sum <- numeric(length(trt_idx))
    realized_ess <- numeric(n_chains)
    tau0_mean <- numeric(n_chains)
    w_mean <- numeric(n_chains)
    spike_fraction <- numeric(n_chains)
    acceptance <- vector("list", n_chains)

    x_train_control <- X_all[ctrl_idx, , drop = FALSE]
    x_test_control <- X_all[trt_idx, , drop = FALSE]
    y_control <- log(Y[ctrl_idx])
    status_control <- as.integer(event[ctrl_idx])
    y_center <- mean(y_control)
    y_centered <- y_control - y_center
    outcome_range <- diff(range(y_centered))

    sigma_fit <- try(survival::survreg(
      survival::Surv(exp(y_control), status_control) ~ .,
      data = data.frame(x_train_control), dist = "lognormal"
    ), silent = TRUE)
    sigma_hat <- if (!inherits(sigma_fit, "try-error") &&
        is.finite(sigma_fit$scale) && sigma_fit$scale > 0)
      sigma_fit$scale else stats::sd(y_centered)
    if (!is.finite(sigma_hat) || sigma_hat <= 0) sigma_hat <- 1
    if (!is.finite(outcome_range) || outcome_range <= 0)
      outcome_range <- 2 * sigma_hat

    cutpoints <- vector("list", ncol(x_train_control))
    for (predictor_id in seq_len(ncol(x_train_control))) {
      unique_x <- sort(unique(x_train_control[, predictor_id]))
      if (length(unique_x) < 2L)
        stop("LRC-BART requires two observed values for predictor ",
             colnames(x_train_control)[predictor_id])
      if (length(unique_x) <= 100L) {
        cutpoints[[predictor_id]] <-
          (unique_x[-1] + unique_x[-length(unique_x)]) / 2
      } else {
        x_range <- range(unique_x)
        cutpoints[[predictor_id]] <-
          x_range[1] + diff(x_range) * seq_len(100L) / 101
      }
    }

    for (chain_id in seq_len(n_chains)) {
      hyper <- list(
        H_f = as.integer(configuration$H_f), H_g = as.integer(H_g),
        alpha_f = 0.95, beta_f = 2, alpha_g = 0.5, beta_g = 3,
        lambda_f_sq = (outcome_range /
                         (4 * sqrt(configuration$H_f)))^2,
        nu_sigma = 3,
        lambda_sigma1 = sigma_hat^2 * stats::qchisq(0.10, 3) / 3,
        lambda_sigma2 = sigma_hat^2 * stats::qchisq(0.10, 3) / 3,
        nu0 = 3, s0_sq = as.numeric(target_row$s0_sq),
        tau1_sq = outcome_range^2 / (4 * 2^2 * H_g),
        a_w = 1, b_w = 1, update_tau0 = FALSE,
        augmentation = TRUE, p_grow = 0.3, p_prune = 0.3,
        n_min = 5L, n_min_g = 0L,
        sigma1_sq_init = sigma_hat^2,
        sigma2_sq_init = sigma_hat^2,
        w_fixed = configuration$w
      )
      control <- list(
        n_burn = n_burn, n_draw = n_draw, thin = 1L,
        verbose = FALSE, keep_trees = TRUE, keep_train = FALSE,
        keep_test = TRUE,
        seed = as.integer(100L * chain_id + 10L * configuration_index)
      )
      fit_control <- clrcbart(
        t(x_train_control), y_centered,
        rep.int(2L, length(ctrl_idx)), status_control,
        y_control - y_center, t(x_test_control),
        cutpoints, hyper, control
      )
      control_prediction <- fit_control$f_test + fit_control$g_test + y_center
      sigma_control <- sqrt(fit_control$sigma1_sq)
      control_rmst[[chain_id]] <- compute_pop_rmst(
        control_prediction, sigma_control, tau_rmst
      )
      # Both controls use EloKRd covariates. The UCMM component excludes g
      # and uses the external-control residual variance.
      ucmm_rmst[[chain_id]] <- compute_pop_rmst(
        fit_control$f_test + y_center,
        sqrt(fit_control$sigma2_sq), tau_rmst
      )
      ratio_chains[, chain_id] <-
        treatment_rmst[[chain_id]] / control_rmst[[chain_id]]
      difference_chains[, chain_id] <-
        treatment_rmst[[chain_id]] - control_rmst[[chain_id]]
      sigma_control_chains[, chain_id] <- sigma_control
      g_draws_all <- rbind(g_draws_all, fit_control$g_test)

      Vg <- cess_forest_gvar(
        fit_control$g_draws, fit_control$cutpoints,
        t(x_test_control), fit_control$tau0_sq, fit_control$tau1_sq
      )
      realized_ess[chain_id] <- mean(
        fit_control$sigma1_sq / (calibration$V_mu_f + Vg)
      )
      for (profile_id in seq_along(trt_idx)) {
        Vg_profile <- cess_forest_gvar(
          fit_control$g_draws, fit_control$cutpoints,
          t(x_test_control[profile_id, , drop = FALSE]),
          fit_control$tau0_sq, fit_control$tau1_sq
        )
        profile_ess_sum[profile_id] <- profile_ess_sum[profile_id] + mean(
          fit_control$sigma1_sq /
            (calibration$V_profile_f[profile_id] + Vg_profile)
        )
      }
      tau0_mean[chain_id] <- mean(fit_control$tau0_sq)
      w_mean[chain_id] <- mean(fit_control$w)
      spike_fraction[chain_id] <- mean(
        fit_control$K0 / pmax(fit_control$L_g, 1)
      )
      acceptance[[chain_id]] <- fit_control$accept
    }

    post_ratio <- as.vector(ratio_chains)
    post_difference <- as.vector(difference_chains)
    rmst_trt <- unlist(treatment_rmst, use.names = FALSE)
    rmst_ctrl <- unlist(control_rmst, use.names = FALSE)
    rmst_ucmm <- unlist(ucmm_rmst, use.names = FALSE)
    sigma_trt <- unlist(lapply(treatment_chains, `[[`, "sigma"),
                        use.names = FALSE)
    sigma_ctrl <- as.vector(sigma_control_chains)
    g_summary <- data.frame(
      profile = seq_along(trt_idx),
      mean = colMeans(g_draws_all),
      sd = apply(g_draws_all, 2, stats::sd),
      lo = apply(g_draws_all, 2, stats::quantile, 0.025),
      hi = apply(g_draws_all, 2, stats::quantile, 0.975)
    )
    ess_map <- data.frame(
      profile = seq_along(trt_idx),
      ess = profile_ess_sum / n_chains
    )

    diagnostic <- c(rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_)
    if (requireNamespace("posterior", quietly = TRUE)) {
      diagnostic <- c(
        rhat = posterior::rhat(log(ratio_chains)),
        ess_bulk = posterior::ess_bulk(log(ratio_chains)),
        ess_tail = posterior::ess_tail(log(ratio_chains))
      )
    }

    results <- list(
      post_ratio = post_ratio,
      delta_hat = median(post_ratio, na.rm = TRUE),
      ci_95 = quantile(post_ratio, c(0.025, 0.975), na.rm = TRUE),
      delta_diff = median(post_difference, na.rm = TRUE),
      ci_diff_95 = quantile(post_difference, c(0.025, 0.975), na.rm = TRUE),
      tau_rmst = tau_rmst,
      rmst_ctrl = rmst_ctrl, rmst_trt = rmst_trt,
      rmst_ucmm = rmst_ucmm,
      rmst_trt_est = arm_summary(rmst_trt),
      rmst_ctrl_est = arm_summary(rmst_ctrl),
      rmst_ucmm_est = arm_summary(rmst_ucmm),
      rmst_ucmm_population = "EloKRd",
      sigma_trt_est = arm_summary(sigma_trt),
      sigma_ctrl_est = arm_summary(sigma_ctrl),
      ratio_chains = ratio_chains,
      difference_chains = difference_chains,
      g_summary = g_summary, ess_map = ess_map,
      calibration = list(
        target_name = target_name,
        target = target_row$target,
        requested = target_row$requested,
        capped = target_row$capped,
        s0_sq = target_row$s0_sq,
        ess_tau0 = target_row$ess_tau0,
        ceiling = calibration$ceiling,
        V_mu_f = calibration$V_mu_f,
        sigma1_sq = calibration$sigma1_sq,
        c_g = calibration$c_g
      ),
      diagnostics = diagnostic,
      sampler = list(
        ess_realized = mean(realized_ess),
        tau0_sq_post = mean(tau0_mean), w = mean(w_mean),
        spike_fraction = mean(spike_fraction), acceptance = acceptance
      ),
      settings = list(
        outcome = OUTCOME, data_tag = data_tag, coding = "harmonized",
        config = configuration$config, H_f = configuration$H_f,
        H_g = H_g, w = configuration$w,
        n_chains = n_chains, n_burn = n_burn, n_draw = n_draw,
        target = target_name,
        target_N = target_row$target, tau_rmst = tau_rmst
      )
    )

    config_suffix <- switch(
      configuration$config,
      default = "",
      Hf50 = "_Hf50",
      `w0.9` = "_w0.9",
      Hf50_w0.9 = "_Hf50_w0.9"
    )
    result_file <- file.path(
      resultDir,
      paste0("LRC-BART_results_", OUTCOME, "_", data_tag,
             config_suffix, "_N", target_name, ".RData")
    )
    saveRDS(results, result_file)
    cat(sprintf(
      "  UCMM control standardized to EloKRd: RMST %.3f [%.3f, %.3f] years\n",
      results$rmst_ucmm_est["est"], results$rmst_ucmm_est["lo"],
      results$rmst_ucmm_est["hi"]
    ))
    cat(sprintf(
      "Saved %s: RMST ratio %.3f [%.3f, %.3f], realized ESS %.1f\n",
      basename(result_file), results$delta_hat, results$ci_95[1],
      results$ci_95[2], results$sampler$ess_realized
    ))
  }
}

cat("\n=== DONE ===\n")
# LRC-BART MODIFICATION END
