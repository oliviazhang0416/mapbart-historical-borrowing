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
# censored AFT treatment fit, replicate, median-ratio/RMST metric, threshold,
# and result-file roles while replacing only the MAP-BART control fit with the
# approved LRC-BART f+g fit and replicate-specific ESS calibration.
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-survival")
n_replicates <- 100L
data_folder <- "data"

library(Rcpp)
library(RcppEigen)
library(survival)
source(file.path(projectDir, "rmst_helpers.R"))
rmst_tau <- 3
rmst_control_sigma <- "adaptive"

sc <- 3
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

# The revision's production workflow uses 1,000 post-warm-up draws after
# 1,000 burn-in iterations for LRC-BART.
ndpost <- 1000L
nskip <- 1000L
keepevery <- 1L
alpha_beta <- data.frame(alpha = 0.95, beta = 2)
ntree <- 50L
threshold <- 0.95

scenario_id <- paste0("sc", sc)
if (sc == 3) scenario_id <- paste0(scenario_id, "_cor", cor)
if (sc == 4) scenario_id <- paste0(scenario_id, "_d", delta_rwd)
if (sc == 5)
  scenario_id <- paste0(scenario_id, "_", region, "_d", delta_rwd)

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

treatment_cpp_cache <- file.path(tempdir(), "sourceCpp_lrcBART_treatment")
dir.create(treatment_cpp_cache, showWarnings = FALSE, recursive = TRUE)
Rcpp::sourceCpp(file.path(.lrcRoot, "aBART", "cabart.cpp"),
                cacheDir = treatment_cpp_cache)
source(file.path(.lrcRoot, "bartModelMatrix.R"), local = TRUE)

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

  # #------------------------------------------------
  # #------------ ESS Curve Construction ------------
  # #------------------------------------------------
  # Build the replicate-specific ESS curve and target table once.
  ess_data_file <- data_file
  ess_H_g <- H_g
  ess_H_f <- ntree
  ess_n_burn <- nskip
  ess_n_draw <- ndpost
  ess_seed <- this_seed
  ess_checkpoint_dir <- file.path(projectDir, "res", "ess")
  source(file.path(projectDir, "ess_local", "ess_cal.R"), local = TRUE)
  calibration <- ess_calibration

  # #------------------------------------------------
  # #---------------- RCT Treatment -----------------
  # #------------------------------------------------
  # Fit the inherited treated-arm BART once for this replicate, or load the
  # saved fit for a sensitivity invocation.
  # LRC-BART MODIFICATION START
  # Keep internal treatment-fit checkpoints separate from scientific results.
  treatment_cache_dir <- file.path(projectDir, "res", "cache")
  if (!dir.exists(treatment_cache_dir))
    dir.create(treatment_cache_dir, recursive = TRUE)
  treatment_cache_file <- file.path(
    treatment_cache_dir,
    paste0("lrcbart_trt_p", p_obs, "_", scenario_id, "_",
           hypo, "_", replicate_id, ".RData")
  )
  # LRC-BART MODIFICATION END
  treatment_settings <- list(
    ndpost = ndpost, nskip = nskip, keepevery = keepevery,
    ntree = ntree, alpha = 0.95, beta = 2,
    seed = as.integer(this_seed)
  )
  treatment_cache_valid <- FALSE
  if (file.exists(treatment_cache_file)) {
    treatment_cached <- readRDS(treatment_cache_file)
    treatment_cache_valid <-
      identical(treatment_cached$settings, treatment_settings) &&
      identical(normalizePath(treatment_cached$data_file),
                normalizePath(data_file)) &&
      file.info(treatment_cache_file)$mtime >= file.info(data_file)$mtime
  }

  if (treatment_cache_valid) {
    treatment <- treatment_cached$fit
  } else {
    if (lrc_config != "default")
      stop("Run the default lrcBART.R configuration before ", lrc_config,
           " so the saved treatment draws are available")

    d <- dat$X[, "D"]
    z <- dat$X[, "Z"]
    x_names <- paste0("X", seq_len(p_obs))
    # LRC-BART MODIFICATION START
    # Retain the original survival project's treated-arm AFT-BART fit.
    y.train <- log(dat$y[d == 1 & z == 1])
    event <- as.integer(dat$delta[d == 1 & z == 1])
    x.train <- data.frame(
      X = as.matrix(dat$X[d == 1 & z == 1, x_names, drop = FALSE])
    )
    x.test <- data.frame(
      X = as.matrix(dat$X[d == 1, x_names, drop = FALSE])
    )
    ntype <- 1
    sparse <- FALSE
    theta <- 0
    omega <- 1
    a <- 0.5
    b <- 1
    augment <- FALSE
    rho <- NULL
    xinfo <- matrix(0.0, 0, 0)
    usequants <- FALSE
    rm.const <- TRUE
    sigest <- NA
    sigdf <- 3
    sigquant <- 0.90
    k <- 2.0
    sigmaf <- NA
    lambda <- NA
    offset <- 0
    weights <- rep(1, length(y.train))
    numcut <- 100L
    printevery <- 10000L
    n <- length(y.train)

    temp <- bartModelMatrix(x.train, numcut, usequants = usequants,
                            xinfo = xinfo, rm.const = rm.const)
    x.train <- t(temp$X)
    numcut <- temp$numcut
    xinfo <- temp$xinfo
    x.test <- bartModelMatrix(x.test)
    x.test <- t(x.test[, temp$rm.const, drop = FALSE])
    rm.const <- temp$rm.const
    grp <- temp$grp
    p <- nrow(x.train)
    np <- ncol(x.test)
    if (!length(rho)) rho <- p
    if (!length(rm.const)) rm.const <- seq_len(p)
    if (!length(grp)) grp <- seq_len(p)
    y.train <- y.train - offset

    if (is.na(lambda)) {
      aft_fit <- survreg(
        Surv(exp(y.train), event) ~ .,
        data = data.frame(t(x.train)), dist = "lognormal"
      )
      sigest <- aft_fit$scale
      qchi <- qchisq(1.0 - sigquant, sigdf)
      lambda <- sigest^2 * qchi / sigdf
    } else {
      sigest <- sqrt(lambda)
    }
    tau <- if (is.na(sigmaf))
      (max(y.train) - min(y.train)) / (2 * k * sqrt(ntree)) else
        sigmaf / sqrt(ntree)

    set.seed(this_seed)
    treatment <- cabart(
      ntype, n, p, np, x.train, y.train, event, x.test,
      ntree, numcut, ndpost * keepevery, nskip, keepevery,
      2, 0.95, offset, tau, sigdf, lambda, sigest, weights,
      sparse, theta, omega, grp, a, b, rho, augment,
      printevery, xinfo, as.integer(this_seed)
    )
    if (nskip > 0)
      treatment$sigma <- treatment$sigma[-seq_len(nskip)]
    # LRC-BART MODIFICATION END
    saveRDS(
      list(fit = treatment, settings = treatment_settings,
           data_file = normalizePath(data_file)),
      treatment_cache_file
    )
  }

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
    # #----------- RCT + Historical Control -----------
    # #------------------------------------------------
    # LRC-BART control fit, written inline in the inherited replicate/target
    # loop so the data preparation and C++ call remain visible together.
    x_names <- paste0("X", seq_len(p_obs))
    d <- dat$X[, "D"]
    z <- dat$X[, "Z"]
    idx_ctrl <- z == 0
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
    sigma_hat1 <- 1
    if (length(y_sigma1) > ncol(x_sigma1) + 5L) {
      sigma_fit1 <- try(survreg(
        Surv(exp(y_sigma1 + y_center),
             status_control[source_control == 1L]) ~ .,
        data = data.frame(x_sigma1), dist = "lognormal"
      ), silent = TRUE)
      if (!inherits(sigma_fit1, "try-error") &&
          is.finite(sigma_fit1$scale) && sigma_fit1$scale > 0)
        sigma_hat1 <- sigma_fit1$scale
    } else if (is.finite(sd(y_sigma1)) && sd(y_sigma1) > 0) {
      sigma_hat1 <- sd(y_sigma1)
    }

    y_sigma2 <- yc[source_control == 2L]
    x_sigma2 <- x_train_control[source_control == 2L, , drop = FALSE]
    sigma_hat2 <- 1
    if (length(y_sigma2) > ncol(x_sigma2) + 5L) {
      sigma_fit2 <- try(survreg(
        Surv(exp(y_sigma2 + y_center),
             status_control[source_control == 2L]) ~ .,
        data = data.frame(x_sigma2), dist = "lognormal"
      ), silent = TRUE)
      if (!inherits(sigma_fit2, "try-error") &&
          is.finite(sigma_fit2$scale) && sigma_fit2$scale > 0)
        sigma_hat2 <- sigma_fit2$scale
    } else if (is.finite(sd(y_sigma2)) && sd(y_sigma2) > 0) {
      sigma_hat2 <- sd(y_sigma2)
    }

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
      tau1_sq = tau1_sq,
      a_w = w_prior[1], b_w = w_prior[2], update_tau0 = TRUE,
      augmentation = TRUE, p_grow = 0.3, p_prune = 0.3,
      n_min = 5L, n_min_g = 5L,
      sigma1_sq_init = sigma_hat1^2,
      sigma2_sq_init = sigma_hat2^2,
      w_fixed = w_fixed
    )
    control <- list(
      n_burn = as.integer(nskip), n_draw = as.integer(ndpost), thin = 1L,
      verbose = FALSE, keep_trees = TRUE, keep_train = FALSE,
      keep_test = TRUE, seed = as.integer(this_seed + 100L + target_index)
    )
    ctrl <- clrcbart(
      t(x_train_control), yc, as.integer(source_control),
      status_control, censor_control, t(x_test_control),
      cutpoints, hyper, control
    )
    ctrl$f_test <- ctrl$f_test + y_center
    ctrl$control_test <- ctrl$f_test + ctrl$g_test
    ctrl$y_center <- y_center
    ctrl$H_g <- H_g
    ctrl$cutpoints <- cutpoints
    ctrl$calibration <- calibration
    # LRC-BART MODIFICATION END

    #-------------------------------------
    #---------- Calculate ATE ------------
    #-------------------------------------
    post_trt <- treatment$yhat.test
    post_ctrl <- ctrl$control_test
    nd <- min(nrow(post_trt), nrow(post_ctrl))
    post_trt <- post_trt[seq_len(nd), , drop = FALSE]
    post_ctrl <- post_ctrl[seq_len(nd), , drop = FALSE]
    # LRC-BART MODIFICATION START
    # Preserve the original survival project's population median-ratio
    # estimand, now using the LRC-BART control log-time surface.
    sigma_trt <- treatment$sigma[seq_len(nd)]
    sigma_ctrl <- sqrt(ctrl$sigma1_sq[seq_len(nd)])
    median_survival_trt <- pop_median_draws(post_trt, sigma_trt)
    median_survival_ctrl <- pop_median_draws(post_ctrl, sigma_ctrl)
    post_samples <- median_survival_trt / median_survival_ctrl
    eff <- exp(dat$treat_eff_true)
    eff_star <- exp(dat$treat_eff_star)
    delta_hat <- median(post_samples, na.rm = TRUE)
    decision <- mean(post_samples > eff_star, na.rm = TRUE)
    threshold_file <- file.path(
      projectDir, "res",
      paste0("LRC-BART_p", p_obs, "_", scenario_id,
             lrc_result_config, "_N", target_name, "_threshold.RData")
    )
    threshold_this <- if (hypo == "alternative" && file.exists(threshold_file))
      readRDS(threshold_file) else threshold

    row <- results[[target_name]][replicate_id, , drop = FALSE]
    row$bias <- delta_hat - eff
    row$sd <- sd(post_samples, na.rm = TRUE)
    row$rmse <- (delta_hat - eff)^2
    row$w1distance <- mean(abs(post_samples - eff), na.rm = TRUE)
    row$w2distance <- sqrt(mean((post_samples - eff)^2, na.rm = TRUE))
    row$ci <- diff(quantile(post_samples, c(0.025, 0.975), na.rm = TRUE))
    row$coverage <- as.numeric(
      quantile(post_samples, 0.025, na.rm = TRUE) <= eff &&
        quantile(post_samples, 0.975, na.rm = TRUE) >= eff
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
    rmst_metrics <- compute_rmst_metrics(
      post_trt, sigma_trt, post_ctrl, sigma_ctrl, dat,
      tau = rmst_tau, threshold = threshold_this,
      control_sigma = rmst_control_sigma
    )
    for (rmst_name in names(rmst_metrics))
      row[[rmst_name]] <- rmst_metrics[[rmst_name]]
    # LRC-BART MODIFICATION END

    Vg <- cess_forest_gvar(
      ctrl$g_draws, ctrl$cutpoints,
      t(dat$X[dat$X[, "D"] == 1,
              paste0("X", seq_len(p_obs)), drop = FALSE]),
      ctrl$tau0_sq, ctrl$tau1_sq
    )
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
        xrct <- dat$X[
          dat$X[, "D"] == 1,
          paste0("X", seq_len(p_obs)), drop = FALSE
        ]
        Vg_region <- cess_forest_gvar(
          ctrl$g_draws, ctrl$cutpoints,
          t(xrct[rct_region, , drop = FALSE]),
          ctrl$tau0_sq, ctrl$tau1_sq
        )
        Vg_outside <- cess_forest_gvar(
          ctrl$g_draws, ctrl$cutpoints,
          t(xrct[!rct_region, , drop = FALSE]),
          ctrl$tau0_sq, ctrl$tau1_sq
        )
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
    cat("Done lrcBART",
        if (nzchar(lrc_config)) paste0(" ", lrc_config) else "",
        " target ", target_name, " replicate ", replicate_id, " of ", n_replicates,
        "\n", sep = "")
  }
}

for (target_name in target_names) {
  threshold_file <- file.path(
    projectDir, "res",
    paste0("LRC-BART_p", p_obs, "_", scenario_id,
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
    paste0("LRC-BART_p", p_obs, "_", scenario_id,
           lrc_result_config, "_N", target_name, "_", hypo, ".RData")
  )
  saveRDS(results[[target_name]], result_file)
}
# LRC-BART MODIFICATION END
