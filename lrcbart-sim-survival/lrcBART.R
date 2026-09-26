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
# This file began as the survival two-arm mBART.R. It retains that script's
# replicate, median-ratio/RMST metric, threshold, and result-file roles.
# The joint revision replaces independent treatment/control fits with one
# treatment-aware censored AFT model and retains control-only ESS calibration.
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-survival")
n_replicates <- 100L
# JOINT LRC-BART MODIFICATION START
dir.create(file.path(projectDir, "res"), recursive = TRUE, showWarnings = FALSE)
# JOINT LRC-BART MODIFICATION END
data_folder <- "data"

library(Rcpp)
library(RcppEigen)
library(survival)
source(file.path(projectDir, "rmst_helpers.R"))
rmst_tau <- 3
# JOINT LRC-BART MODIFICATION START
rmst_control_sigma <- "own"  # both trial arms use the same sigma1 draw
# JOINT LRC-BART MODIFICATION END

# JOINT LRC-BART MODIFICATION START
sc <- 3
variant <- ""
# JOINT LRC-BART MODIFICATION END
hypo <- "alternative"
cor <- if (sc == 3) 0.5 else 1
delta_rwd <- 1
region <- "X5"

# The default value runs the six target fits. run_all.R overrides it with w0,
# w1, or Hg10 only for a sensitivity invocation.
lrc_config <- "default"
stopifnot(lrc_config %in% c("default", "w0", "w1", "Hg10"))
# LRC-BART MODIFICATION START
# The primary configuration uses the concise inherited result naming. Only
# optional sensitivity configurations receive an identifying suffix.
lrc_result_config <- if (lrc_config == "default") "" else
  paste0("_", lrc_config)
# LRC-BART MODIFICATION END

# JOINT LRC-BART MODIFICATION START
# Sampling schedule matches lrcbart-historical-borrowing; AFT fitting remains joint.
n_chains <- 1L
joint_n_burn <- 1000L
joint_n_draw <- 1000L
joint_thin <- 1L
g_sweeps <- 1L
calibration_n_burn <- 1000L
calibration_n_draw <- 1000L
stopifnot(n_chains >= 1L, n_chains == as.integer(n_chains),
          joint_n_burn >= 0L, joint_n_burn == as.integer(joint_n_burn),
          joint_n_draw >= 4L, joint_n_draw == as.integer(joint_n_draw),
          joint_thin >= 1L, joint_thin == as.integer(joint_thin),
          g_sweeps >= 1L, g_sweeps == as.integer(g_sweeps))
if (!requireNamespace("digest", quietly = TRUE))
  stop("Joint result provenance requires the digest package")
# JOINT LRC-BART MODIFICATION END
alpha_beta <- data.frame(alpha = 0.95, beta = 2)
ntree <- 50L
threshold <- 0.95

# JOINT LRC-BART MODIFICATION START
stopifnot(length(sc) == 1L, sc %in% 1:3,
          length(variant) == 1L, variant %in% c("", "a", "b", "c", "c-i", "c-ii"),
          sc == 1 || variant == "", hypo %in% c("null", "alternative"),
          length(delta_rwd) == 1L, is.finite(delta_rwd))
if (variant == "c-ii" && hypo != "null")
  stop("Sc1c-ii requires hypo = 'null'")
if (variant == "c-i" && hypo != "alternative")
  stop("Sc1c-i requires hypo = 'alternative'; use Sc1c-ii for the null study")
if (variant == "") {
  delta_rwd <- 0
  region <- "none"
} else if (variant == "a") {
  region <- "none"
} else if (variant == "b") {
  stopifnot(region %in% c("X5", "X7"))
} else {
  stopifnot(region == "X5X7", delta_rwd == 2)
}
scenario_id <- paste0("sc", sc, variant)
if (sc == 3) scenario_id <- paste0(scenario_id, "_cor", cor)
if (variant == "b") scenario_id <- paste0(scenario_id, "_", region)
if (nzchar(variant)) scenario_id <- paste0(scenario_id, "_d", delta_rwd)
# JOINT LRC-BART MODIFICATION END

sample_files <- list.files(file.path(projectDir, data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (!length(sample_files))
  stop("No data files found in ", file.path(projectDir, data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
data_prefix <- paste0("data_p", p_obs, "_", scenario_id, "_", hypo, "_")
replicate_files <- Sys.glob(file.path(projectDir, data_folder,
                                      paste0(data_prefix, "*.RData")))
if (length(replicate_files) < n_replicates)
  stop("Expected ", n_replicates, " simulated data files for ", scenario_id,
       "; found ", length(replicate_files))

# Load the active standalone source directly. No model or ESS C++ copy lives
# inside this research subproject.
lrc_cpp_cache <- file.path(tempdir(), "sourceCpp_lrcBART")
dir.create(lrc_cpp_cache, showWarnings = FALSE, recursive = TRUE)
Rcpp::sourceCpp(file.path(.lrcRoot, "lrcBART", "clrcbart.cpp"),
                cacheDir = lrc_cpp_cache)
Rcpp::sourceCpp(file.path(.lrcRoot, "lrcBART", "cess.cpp"),
                cacheDir = lrc_cpp_cache)

# JOINT LRC-BART MODIFICATION START
# The joint call replaces the independent treated-arm fit and its cache.
model_sources <- c(list.files(file.path(.lrcRoot, "lrcBART"),
                              pattern = "\\.(cpp|h)$", full.names = TRUE, recursive = TRUE),
                   file.path(projectDir, "lrcBART.R"),
                   file.path(projectDir, "rmst_helpers.R"),
                   file.path(projectDir, "ess_local", "ess_cal.R"))
source_hash <- digest::digest(unname(tools::md5sum(sort(model_sources))),
                              algo = "sha256")
# Reporting-only storage: chain draws and intermediate summaries remain in memory.
# No chain checkpoint or diagnostic files are written.
# JOINT LRC-BART MODIFICATION END

H_g <- if (lrc_config == "Hg10") 10L else 5L
w_fixed <- if (lrc_config == "w0") 0 else
  if (lrc_config == "w1") 1 else NA_real_
w_prior <- if (lrc_config == "Hg10") c(4, 1) else c(1, 1)
target_names <- if (lrc_config == "default")
  c("100", "75", "50", "f90", "f50", "f25") else "100"

set.seed(6)
seed <- sample.int(10000000L, n_replicates, replace = FALSE)
results <- list()
decisions <- list()
result_columns <- c(
  "bias", "sd", "rmse", "w1distance", "w2distance", "ci",
  "coverage", "tp_calibrated", "tp", "fp", "pehe_subj",
  "bias_subj", "bias.trt.sigma", "sd.trt.sigma",
  "w2distance.trt.sigma", "bias.ctrl.sigma", "sd.ctrl.sigma",
  "w2distance.ctrl.sigma", "bias.trt.median.pop",
  "w2distance.trt.median.pop", "bias.ctrl.median.pop",
  "w2distance.ctrl.median.pop", "s0_sq", "ess_prior",
  "ess_realized", "ess_ceiling", "ess_target",
  "ess_target_requested", "ess_capped", "w", "spike_frac",
  "map_mean", "map_region", "map_outside", "g_region",
  "g_outside", "ess_region", "ess_outside",
  "estimate", "lower", "upper", "median_interval_excludes_one", "p_above_eff_star",
  "psi_hat", "psi_truth", "psi_bias", "psi_lower", "psi_upper",
  "psi_interval_excludes_zero", "rmst_lower", "rmst_upper", "rmst_interval_excludes_one",
  "rmst_tau", "rmst_true", "rmst_hat", "bias_rmst", "sd_rmst",
  "ci_rmst", "coverage_rmst", "decision_rmst", "tp_rmst",
  "rmse_rmst", "w1distance_rmst", "w2distance_rmst",
  "tp_calibrated_rmst", "bias.trt.rmst.pop",
  "w2distance.trt.rmst.pop", "bias.ctrl.rmst.pop",
  "w2distance.ctrl.rmst.pop", "bias_subj_rmst", "pehe_subj_rmst"
)
for (target_name in target_names) {
  results[[target_name]] <- data.frame(
    alpha = 0.95, beta = 2, H = ntree, c = cor,
    iteration = seq_len(n_replicates), stringsAsFactors = FALSE
  )
  for (result_name in result_columns)
    results[[target_name]][[result_name]] <- NA_real_
  results[[target_name]]$lrc_config <- lrc_config
  results[[target_name]]$target <- target_name
  # JOINT LRC-BART ADDITION START
  results[[target_name]]$scenario_id <- scenario_id
  results[[target_name]]$hypothesis <- hypo
  results[[target_name]]$joint_model <- TRUE
  results[[target_name]]$effect_estimand <- "trial-standardized population median ratio"
  results[[target_name]]$n_chains <- n_chains
  results[[target_name]]$n_burn <- joint_n_burn
  results[[target_name]]$n_draw <- joint_n_draw
  results[[target_name]]$thin <- joint_thin
  results[[target_name]]$g_sweeps <- g_sweeps
  results[[target_name]]$tau_prior_var <- 100
  results[[target_name]]$data_hash <- NA_character_
  results[[target_name]]$source_hash <- source_hash
  results[[target_name]]$calibration_hash <- NA_character_
  # A previous aggregate is not a completed summary of this invocation.
  old_summary <- file.path(projectDir, "res", paste0("LRC-BART-joint_p", p_obs,
    "_", scenario_id, lrc_result_config, "_N", target_name, "_", hypo, ".RData"))
  if (file.exists(old_summary)) unlink(old_summary)
  # JOINT LRC-BART ADDITION END
  decisions[[target_name]] <- numeric(n_replicates)
}

for (replicate_id in seq_len(n_replicates)) {
  data_file <- file.path(
    projectDir, data_folder,
    paste0("data_p", p_obs, "_", scenario_id, "_", hypo, "_",
           replicate_id, ".RData")
  )
  dat <- readRDS(data_file)
  this_seed <- seed[replicate_id]
  # JOINT LRC-BART ADDITION START
  stopifnot(identical(dat$scenario_id, scenario_id),
            identical(dat$hypothesis, hypo))
  data_hash <- digest::digest(file = data_file, algo = "sha256")
  # JOINT LRC-BART ADDITION END

  # #------------------------------------------------
  # #------------ ESS Curve Construction ------------
  # #------------------------------------------------
  # Build the replicate-specific ESS curve and target table once.
  ess_data_file <- data_file
  ess_H_g <- H_g
  ess_H_f <- ntree
  # JOINT LRC-BART MODIFICATION START
  ess_n_burn <- calibration_n_burn
  ess_n_draw <- calibration_n_draw
  # JOINT LRC-BART MODIFICATION END
  ess_seed <- this_seed
  ess_checkpoint_dir <- file.path(projectDir, "res", "ess")
  source(file.path(projectDir, "ess_local", "ess_cal.R"), local = TRUE)
  calibration <- ess_calibration

  # JOINT LRC-BART MODIFICATION START
  # Calibrate once per dataset; every target and chain shares this checkpoint.
  calibration_hash <- digest::digest(calibration, algo = "sha256")
  # JOINT LRC-BART MODIFICATION END

  # #------------------------------------------------
  # #--------------- ESS Calibration ----------------
  # #------------------------------------------------
  for (target_index in seq_along(target_names)) {
    target_name <- target_names[target_index]
    target_row <- calibration$targets[
      calibration$targets$target_name == target_name, , drop = FALSE
    ]
    if (nrow(target_row) != 1L)
      stop("Missing ESS calibration target ", target_name)

    # #------------------------------------------------
    # #-------------- Joint LRC-BART Fit --------------
    # #------------------------------------------------
    # JOINT LRC-BART MODIFICATION START
    # Joint likelihood: all trial rows and external controls, in original order.
    # Model: f(X) + I(trial)*g(X) + A*psi.
    # Model applies to latent log event times, with censoring augmentation.
    x_names <- paste0("X", seq_len(p_obs))
    d <- dat$X[, "D"]
    z <- dat$X[, "Z"]
    stopifnot(all(d %in% c(0, 1)), all(z %in% c(0, 1)),
              all(z[d == 0] == 0), any(d == 1 & z == 0), any(d == 1 & z == 1))
    idx_ctrl <- rep(TRUE, length(z))  # inherited name now selects all training rows
    treatment_joint <- as.numeric(z)
    x_train_control <- as.matrix(
      dat$X[idx_ctrl, x_names, drop = FALSE]
    )
    x_test_control <- as.matrix(
      dat$X[d == 1, x_names, drop = FALSE]
    )
    storage.mode(x_train_control) <- "double"
    storage.mode(x_test_control) <- "double"
    # LRC-BART MODIFICATION START
    # Censoring is represented by status=0 and a lower bound on latent log
    # survival time. clrcbart performs the truncated-normal augmentation.
    stopifnot(length(dat$y) == length(d), all(is.finite(dat$y)), all(dat$y > 0),
              length(dat$delta) == length(d), all(dat$delta %in% c(0, 1)))
    # Only observed times and statuses enter the likelihood, never y_true.
    y_control <- log(as.numeric(dat$y[idx_ctrl]))
    status_control <- as.integer(dat$delta[idx_ctrl])
    source_control <- ifelse(d[idx_ctrl] == 1, 1L, 2L)
    y_center <- mean(y_control)
    yc <- y_control - y_center
    censor_control <- y_control - y_center
    outcome_range <- diff(range(yc))
    lambda_f_sq <- (outcome_range / (4 * sqrt(ntree)))^2

    y_sigma1 <- yc[source_control == 1L]
    x_sigma1 <- x_train_control[source_control == 1L, , drop = FALSE]
    event_sigma1 <- status_control[source_control == 1L]
    sigma_hat1 <- if (length(y_sigma1) < 2L) 1 else
      if (any(event_sigma1 == 0L)) sd(y_sigma1[event_sigma1 == 1L]) else sd(y_sigma1)
    if (length(y_sigma1) > ncol(x_sigma1) + 5L) {
      if (any(event_sigma1 == 0L)) {
        sigma_fit1 <- try(survreg(
          Surv(exp(y_sigma1), event_sigma1) ~ .,
          data = data.frame(x_sigma1), dist = "lognormal"), silent = TRUE)
        if (!inherits(sigma_fit1, "try-error") &&
            is.finite(sigma_fit1$scale) && sigma_fit1$scale > 0)
          sigma_hat1 <- sigma_fit1$scale
      } else {
        sigma_fit1 <- try(lm(y_sigma1 ~ ., data = data.frame(y_sigma1, x_sigma1)), silent = TRUE)
        if (!inherits(sigma_fit1, "try-error") &&
            is.finite(summary(sigma_fit1)$sigma) && summary(sigma_fit1)$sigma > 0)
          sigma_hat1 <- summary(sigma_fit1)$sigma
      }
    }

    y_sigma2 <- yc[source_control == 2L]
    x_sigma2 <- x_train_control[source_control == 2L, , drop = FALSE]
    event_sigma2 <- status_control[source_control == 2L]
    sigma_hat2 <- if (length(y_sigma2) < 2L) 1 else
      if (any(event_sigma2 == 0L)) sd(y_sigma2[event_sigma2 == 1L]) else sd(y_sigma2)
    if (length(y_sigma2) > ncol(x_sigma2) + 5L) {
      if (any(event_sigma2 == 0L)) {
        sigma_fit2 <- try(survreg(
          Surv(exp(y_sigma2), event_sigma2) ~ .,
          data = data.frame(x_sigma2), dist = "lognormal"), silent = TRUE)
        if (!inherits(sigma_fit2, "try-error") &&
            is.finite(sigma_fit2$scale) && sigma_fit2$scale > 0)
          sigma_hat2 <- sigma_fit2$scale
      } else {
        sigma_fit2 <- try(lm(y_sigma2 ~ ., data = data.frame(y_sigma2, x_sigma2)), silent = TRUE)
        if (!inherits(sigma_fit2, "try-error") &&
            is.finite(summary(sigma_fit2)$sigma) && summary(sigma_fit2)$sigma > 0)
          sigma_hat2 <- summary(sigma_fit2)$sigma
      }
    }

    if (!is.finite(outcome_range) || outcome_range <= 0 ||
        !is.finite(sigma_hat1) || sigma_hat1 <= 0 ||
        !is.finite(sigma_hat2) || sigma_hat2 <= 0)
      stop("Joint AFT fitting requires positive range and estimable source residual scales")
    qchi <- qchisq(0.10, df = 3) / 3
    tau1_sq <- outcome_range^2 / (4 * 2^2 * H_g)
    cutpoints <- vector("list", ncol(x_train_control))
    for (j in seq_len(ncol(x_train_control))) {
      unique_x <- sort(unique(x_train_control[, j]))
      if (length(unique_x) < 2L)
        stop("LRC-BART requires at least two observed values for predictor ", j)
      if (length(unique_x) <= 100L) {
        cutpoints[[j]] <-
          (unique_x[-1] + unique_x[-length(unique_x)]) / 2
      } else {
        x_range <- range(unique_x)
        cutpoints[[j]] <-
          x_range[1] + diff(x_range) * seq_len(100L) / 101
      }
    }
    hyper <- list(
      H_f = as.integer(ntree), H_g = as.integer(H_g),
      alpha_f = 0.95, beta_f = 2, alpha_g = 0.5, beta_g = 3,
      lambda_f_sq = lambda_f_sq, nu_sigma = 3,
      lambda_sigma1 = sigma_hat1^2 * qchi,
      lambda_sigma2 = sigma_hat2^2 * qchi,
      nu0 = 3, s0_sq = as.numeric(target_row$s0_sq),
      tau1_sq = tau1_sq, tau_prior_var = 100,
      a_w = w_prior[1], b_w = w_prior[2], update_tau0 = TRUE,
      augmentation = TRUE, p_grow = 0.3, p_prune = 0.3,
      n_min = 5L, n_min_g = 5L,
      sigma1_sq_init = sigma_hat1^2,
      sigma2_sq_init = sigma_hat2^2,
      w_fixed = w_fixed
    )
    control <- list(
      n_burn = as.integer(joint_n_burn), n_draw = as.integer(joint_n_draw),
      thin = as.integer(joint_thin), g_sweeps = as.integer(g_sweeps),
      verbose = FALSE, keep_trees = TRUE, keep_train = FALSE, keep_test = TRUE
    )
    chain_results <- vector("list", n_chains)
    for (chain_id in seq_len(n_chains)) {
      # Stable target/configuration seeds, independent of loop ordering.
      seed_key <- list(this_seed, scenario_id, hypo, replicate_id, target_name, lrc_config)
      seed_base <- strtoi(substr(digest::digest(seed_key, algo = "sha256"), 1, 7), 16L)
      control$seed <- as.integer((seed_base + chain_id) %% 2147483646 + 1)
      fit_chain <- clrcbart(
        t(x_train_control), yc, as.integer(source_control),
        status_control, censor_control, t(x_test_control),
        cutpoints, hyper, control, treatment = treatment_joint
      )
      stopifnot(isTRUE(fit_chain$joint_model), length(fit_chain$tau) == joint_n_draw)
      fit_chain$control_test <- fit_chain$f_test + fit_chain$g_test + y_center
      # Consume forests for inherited borrowing summaries before discarding them.
      fit_chain$Vg <- cess_forest_gvar(fit_chain$g_draws, cutpoints,
        t(x_test_control), fit_chain$tau0_sq, fit_chain$tau1_sq)
      fit_chain$Vg_region <- fit_chain$Vg_outside <- rep(NA_real_, joint_n_draw)
      R <- dat$region_rct
      if (!is.null(R) && any(R) && any(!R)) {
        if (target_name == "100") {
          fit_chain$Vg_region <- cess_forest_gvar(fit_chain$g_draws, cutpoints,
            t(x_test_control[R, , drop = FALSE]), fit_chain$tau0_sq, fit_chain$tau1_sq)
          fit_chain$Vg_outside <- cess_forest_gvar(fit_chain$g_draws, cutpoints,
            t(x_test_control[!R, , drop = FALSE]), fit_chain$tau0_sq, fit_chain$tau1_sq)
        }
      }
      fit_chain <- fit_chain[c("tau", "control_test", "g_test", "sigma1_sq",
        "sigma2_sq", "tau0_sq", "tau1_sq", "w", "K0", "L_g", "Vg",
        "Vg_region", "Vg_outside")]
      chain_results[[chain_id]] <- fit_chain
      cat("Completed", scenario_id, "target", target_name,
          "replicate", replicate_id, "chain", chain_id, "of", n_chains, "\n")
    }
    # Pool all complete chains in memory for the final reporting calculations.
    stopifnot(length(chain_results) == n_chains,
      all(vapply(chain_results, function(x) length(x$tau) == joint_n_draw, logical(1))))
    ctrl <- list()
    for (field in c("control_test", "g_test"))
      ctrl[[field]] <- do.call(rbind, lapply(chain_results, `[[`, field))
    for (field in c("tau", "sigma1_sq", "sigma2_sq", "w", "K0", "L_g",
                    "Vg", "Vg_region", "Vg_outside"))
      ctrl[[field]] <- unlist(lapply(chain_results, `[[`, field), use.names = FALSE)
    rm(chain_results, fit_chain)
    # JOINT LRC-BART MODIFICATION END

    # JOINT LRC-BART MODIFICATION START
    # #------------------------------------------------
    # #--------------- RCT Control --------------------
    # #------------------------------------------------
    # Control log-time means: f(X) + g(X), with centering restored.
    post_ctrl <- ctrl$control_test

    # #------------------------------------------------
    # #-------------- RCT Treatment -------------------
    # #------------------------------------------------
    # Treatment log-time means: f(X) + g(X) + psi.
    post_trt <- sweep(post_ctrl, 1L, ctrl$tau, "+")
    # Both trial arms share the same posterior residual variance.

    # #------------------------------------------------
    # #---------- Calculate Population Median Ratio ----------
    # #------------------------------------------------
    nd <- length(ctrl$tau)
    sigma_trt <- sigma_ctrl <- sqrt(ctrl$sigma1_sq)
    # Solve both mixture medians on log time without the legacy 0.05/30 clips.
    median_survival_trt <- pop_median_draws(post_trt, sigma_trt, log_scale = TRUE)
    median_survival_ctrl <- pop_median_draws(post_ctrl, sigma_ctrl, log_scale = TRUE)
    post_samples <- median_survival_trt / median_survival_ctrl
    if (any(!is.finite(post_samples)) ||
        max(abs(log(post_samples) - ctrl$tau)) > 1e-7)
      stop("Joint population median ratio is inconsistent with exp(psi)")
    eff <- exp(dat$treat_eff_true)
    # JOINT LRC-BART MODIFICATION END
    eff_star <- exp(dat$treat_eff_star)
    delta_hat <- median(post_samples, na.rm = TRUE)
    decision <- mean(post_samples > eff_star, na.rm = TRUE)
    threshold_file <- file.path(
      projectDir, "res",
      paste0("LRC-BART-joint_p", p_obs, "_", scenario_id,
             lrc_result_config, "_N", target_name, "_threshold.RData")
    )
    threshold_this <- if (hypo == "alternative" && file.exists(threshold_file))
      readRDS(threshold_file) else threshold

    row <- results[[target_name]][replicate_id, , drop = FALSE]
    # JOINT LRC-BART MODIFICATION START
    row$estimate <- delta_hat
    row$psi_hat <- mean(ctrl$tau)
    row$psi_truth <- dat$mean_log_time_effect
    row$psi_bias <- row$psi_hat - row$psi_truth
    row$psi_lower <- unname(quantile(ctrl$tau, 0.025))
    row$psi_upper <- unname(quantile(ctrl$tau, 0.975))
    row$psi_interval_excludes_zero <- as.numeric(row$psi_lower > 0 || row$psi_upper < 0)
    row$lower <- exp(row$psi_lower)
    row$upper <- exp(row$psi_upper)
    row$median_interval_excludes_one <- as.numeric(row$lower > 1 || row$upper < 1)
    row$p_above_eff_star <- decision
    row$data_hash <- data_hash
    row$calibration_hash <- calibration_hash
    row$bias <- delta_hat - eff
    # JOINT LRC-BART MODIFICATION END
    row$sd <- sd(post_samples, na.rm = TRUE)
    row$rmse <- (delta_hat - eff)^2
    row$w1distance <- mean(abs(post_samples - eff), na.rm = TRUE)
    row$w2distance <- sqrt(mean((post_samples - eff)^2, na.rm = TRUE))
    row$ci <- row$upper - row$lower
    row$coverage <- as.numeric(
      row$lower <= eff && row$upper >= eff
    )
    row$tp_calibrated <- as.numeric(decision > threshold_this)
    row$tp <- as.numeric(decision > 0.95)
    row$fp <- if (hypo == "null") row$tp else NA_real_

    #------------------------------------------
    #------ Calculate subject-wise ------------
    #------------------------------------------
    true_i <- exp(dat$eff_i)
    delta_i <- apply(exp(post_trt - post_ctrl), 2, median, na.rm = TRUE)
    row$bias_subj <- mean(delta_i - true_i)
    row$pehe_subj <- sqrt(mean((delta_i - true_i)^2))

    #-------------------------------------
    #------ Variance estimation ----------
    #-------------------------------------
    row$bias.trt.sigma <- mean(sigma_trt) - dat$sigma_rct
    row$sd.trt.sigma <- sd(sigma_trt)
    row$w2distance.trt.sigma <-
      sqrt(mean((sigma_trt - dat$sigma_rct)^2))
    row$bias.ctrl.sigma <- mean(sigma_ctrl) - dat$sigma_rct
    row$sd.ctrl.sigma <- sd(sigma_ctrl)
    row$w2distance.ctrl.sigma <-
      sqrt(mean((sigma_ctrl - dat$sigma_rct)^2))

    #--------------------------------------------------
    #------ Population median survival time -----------
    #--------------------------------------------------
    row$bias.trt.median.pop <-
      median(median_survival_trt, na.rm = TRUE) - dat$true_median_trt_pop
    row$w2distance.trt.median.pop <- sqrt(mean(
      (median_survival_trt - dat$true_median_trt_pop)^2, na.rm = TRUE
    ))
    row$bias.ctrl.median.pop <-
      median(median_survival_ctrl, na.rm = TRUE) - dat$true_median_ctrl_pop
    row$w2distance.ctrl.median.pop <- sqrt(mean(
      (median_survival_ctrl - dat$true_median_ctrl_pop)^2, na.rm = TRUE
    ))

    #-------------------------------------
    #---------- Calculate RMST -----------
    #-------------------------------------
    # JOINT LRC-BART MODIFICATION START
    # Match the generator's 200-point RMST grid. The machine floor only prevents
    # division by zero; do not impose the legacy 0.05 floor on joint predictions.
    rmst_metrics <- compute_rmst_metrics(
      post_trt, sigma_trt, post_ctrl, sigma_ctrl, dat,
      tau = rmst_tau, threshold = threshold_this,
      control_sigma = rmst_control_sigma, K = 200L,
      floor = .Machine$double.xmin, joint = TRUE
    )
    # JOINT LRC-BART MODIFICATION END
    for (rmst_name in names(rmst_metrics))
      row[[rmst_name]] <- rmst_metrics[[rmst_name]]
    # LRC-BART MODIFICATION END

    # JOINT LRC-BART MODIFICATION START
    Vg <- ctrl$Vg
    # JOINT LRC-BART MODIFICATION END
    row$s0_sq <- target_row$s0_sq
    row$ess_prior <- target_row$ess_tau0
    row$ess_realized <- mean(
      ctrl$sigma1_sq / (calibration$V_mu_f + Vg)
    )
    row$ess_ceiling <- calibration$ceiling
    row$ess_target <- target_row$target
    row$ess_target_requested <- target_row$requested
    row$ess_capped <- as.numeric(target_row$capped)
    row$w <- mean(ctrl$w)
    row$spike_frac <- mean(ctrl$K0 / pmax(ctrl$L_g, 1))
    borrow_draw <- abs(ctrl$g_test) < 0.5
    rct_control_profile <-
      dat$X[dat$X[, "D"] == 1, "Z"] == 0
    row$map_mean <- mean(
      borrow_draw[, rct_control_profile, drop = FALSE]
    )

    rct_region <- dat$region_rct
    if (!is.null(rct_region) && any(rct_region) && any(!rct_region)) {
      rct_control_region <- rct_region[rct_control_profile]
      if (any(rct_control_region) && any(!rct_control_region)) {
        row$map_region <- mean(
          borrow_draw[, rct_control_profile, drop = FALSE][,
            rct_control_region, drop = FALSE]
        )
        row$map_outside <- mean(
          borrow_draw[, rct_control_profile, drop = FALSE][,
            !rct_control_region, drop = FALSE]
        )
        row$g_region <- mean(
          ctrl$g_test[, rct_control_profile, drop = FALSE][,
            rct_control_region, drop = FALSE]
        )
        row$g_outside <- mean(
          ctrl$g_test[, rct_control_profile, drop = FALSE][,
            !rct_control_region, drop = FALSE]
        )
      }
      if (target_name == "100") {
        # JOINT LRC-BART MODIFICATION START
        Vg_region <- ctrl$Vg_region
        Vg_outside <- ctrl$Vg_outside
        # JOINT LRC-BART MODIFICATION END
        row$ess_region <- mean(
          ctrl$sigma1_sq /
            (calibration$V_mu_f_region[["region"]] + Vg_region)
        )
        row$ess_outside <- mean(
          ctrl$sigma1_sq /
            (calibration$V_mu_f_region[["outside"]] + Vg_outside)
        )
      }
    }

    results[[target_name]][replicate_id, ] <- row
    decisions[[target_name]][replicate_id] <- decision
    # JOINT LRC-BART ADDITION START
    # Save completed reporting rows, never posterior draws or chain summaries.
    report_file <- file.path(projectDir, "res", paste0("LRC-BART-joint_p", p_obs,
      "_", scenario_id, lrc_result_config, "_N", target_name, "_", hypo, ".RData"))
    report_rows <- results[[target_name]][seq_len(replicate_id), , drop = FALSE]
    attr(report_rows, "planned_replicates") <- n_replicates
    attr(report_rows, "complete") <- replicate_id == n_replicates
    report_tmp <- tempfile("report-", tmpdir = file.path(projectDir, "res"))
    saveRDS(report_rows, report_tmp)
    if (!file.rename(report_tmp, report_file)) stop("Cannot publish reporting results")
    rm(ctrl, post_ctrl, post_trt, post_samples, borrow_draw)
    # JOINT LRC-BART ADDITION END
    cat("Done lrcBART",
        if (nzchar(lrc_config)) paste0(" ", lrc_config) else "",
        " target ", target_name, " replicate ", replicate_id, " of ", n_replicates,
        "\n", sep = "")
  }
}

for (target_name in target_names) {
  threshold_file <- file.path(
    projectDir, "res",
    paste0("LRC-BART-joint_p", p_obs, "_", scenario_id,
           lrc_result_config, "_N", target_name, "_threshold.RData")
  )
  if (hypo == "null") {
    candidates <- seq(0.1, 0.999, by = 0.001)
    type1 <- numeric(length(candidates))
    for (candidate_index in seq_along(candidates))
      type1[candidate_index] <- mean(
        decisions[[target_name]] > candidates[candidate_index], na.rm = TRUE
      )
    calibrated_threshold <- candidates[which.min(abs(type1 - 0.05))]
    saveRDS(calibrated_threshold, threshold_file)
  }
  result_file <- file.path(
    projectDir, "res",
    paste0("LRC-BART-joint_p", p_obs, "_", scenario_id,
           lrc_result_config, "_N", target_name, "_", hypo, ".RData")
  )
  attr(results[[target_name]], "planned_replicates") <- n_replicates
  attr(results[[target_name]], "complete") <- TRUE
  saveRDS(results[[target_name]], result_file)
}
# LRC-BART MODIFICATION END
