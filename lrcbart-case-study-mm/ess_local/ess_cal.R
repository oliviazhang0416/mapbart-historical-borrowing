# LRC-BART ADDITION START
# Censored-survival single-arm ESS curve for one case-study endpoint and H_f.
# This file is sourced by lrcBART.R and returns `ess_calibration`.

ess_required <- c(
  "ess_data_file", "ess_outcome", "ess_data_tag", "ess_X_all", "ess_time",
  "ess_status", "ess_trt", "ess_H_g", "ess_H_f", "ess_n_burn",
  "ess_n_draw", "ess_seed", "ess_checkpoint_dir"
)
for (ess_name in ess_required)
  if (!exists(ess_name, inherits = TRUE))
    stop("ess_cal.R requires `", ess_name, "` from lrcBART.R")
if (!exists("clrcbart", mode = "function", inherits = TRUE) ||
    !exists("cess_cg", mode = "function", inherits = TRUE) ||
    !exists("cess_ess_tau0", mode = "function", inherits = TRUE))
  stop("lrcBART.R must load clrcbart.cpp and cess.cpp before ess_cal.R")

ess_alpha_g <- 0.5
ess_beta_g <- 3
ess_nu0 <- 3
ess_q <- 0.95

# #------------------------------------------------
# #------------ ESS Curve Construction ------------
# #------------------------------------------------
# This is the supplied s_0^2 range used to construct the calibration curve.
ess_grid <- 10^seq(-6, 1, length.out = 141)
ess_n_blocks <- 20L
ess_n_tau <- 4000L
ess_n_cg <- 2000L

dir.create(ess_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
# H_g=5 is the application default and is intentionally omitted from names.
ess_checkpoint <- file.path(
  ess_checkpoint_dir,
  paste0("ess_", tolower(ess_outcome), "_", ess_data_tag,
         "_Hf", ess_H_f, ".RData")
)
ess_data_mtime <- file.info(ess_data_file)$mtime
ess_settings <- list(
  outcome = ess_outcome, H_g = as.integer(ess_H_g),
  H_f = as.integer(ess_H_f), alpha_g = ess_alpha_g,
  beta_g = ess_beta_g, nu0 = ess_nu0, q = ess_q,
  grid = as.numeric(ess_grid), n_blocks = as.integer(ess_n_blocks),
  n_tau = as.integer(ess_n_tau), n_cg = as.integer(ess_n_cg),
  n_burn = as.integer(ess_n_burn), n_draw = as.integer(ess_n_draw),
  seed = as.integer(ess_seed), augmentation = TRUE,
  # JOINT LRC-BART ADDITION START
  data_hash = unname(tools::md5sum(ess_data_file)),
  code_hash = tools::md5sum(c(
    file.path(mainDir, "lrcBART", "clrcbart.cpp"),
    file.path(mainDir, "lrcBART", "cess.cpp"),
    list.files(file.path(mainDir, "lrcBART", "include.lrc"),
               full.names = TRUE, recursive = TRUE),
    file.path(projDir, "ess_local", "ess_cal.R")
  )),
  # JOINT LRC-BART ADDITION END
  coding = "harmonized"
)

ess_target_names <- c("100", "75", "50", "f90", "f50", "f25")
ess_use_checkpoint <- FALSE
if (file.exists(ess_checkpoint) &&
    file.info(ess_checkpoint)$mtime >= ess_data_mtime) {
  ess_cached <- tryCatch(readRDS(ess_checkpoint), error = function(e) NULL)
  ess_use_checkpoint <- is.list(ess_cached) && identical(ess_cached$settings, ess_settings) &&
    identical(normalizePath(ess_cached$data_file),
              normalizePath(ess_data_file)) &&
    identical(as.character(ess_cached$targets$target_name), ess_target_names)
}

if (ess_use_checkpoint) {
  ess_calibration <- ess_cached
} else {
  ess_X_all <- as.matrix(ess_X_all)
  storage.mode(ess_X_all) <- "double"
  ess_idx_rwd <- ess_trt == 0L
  ess_idx_rct <- ess_trt == 1L
  if (!any(ess_idx_rwd) || !any(ess_idx_rct))
    stop("Case-study calibration requires UCMM controls and EloKRd profiles")

  # Stage 1: plain censored AFT-BART on UCMM, evaluated at EloKRd profiles.
  ess_stage1_y <- log(as.numeric(ess_time[ess_idx_rwd]))
  ess_stage1_status <- as.integer(ess_status[ess_idx_rwd])
  ess_stage1_x <- ess_X_all[ess_idx_rwd, , drop = FALSE]
  ess_stage1_xtest <- ess_X_all[ess_idx_rct, , drop = FALSE]
  ess_stage1_center <- mean(ess_stage1_y)
  ess_stage1_yc <- ess_stage1_y - ess_stage1_center
  ess_stage1_censor <- ess_stage1_y - ess_stage1_center
  ess_stage1_range <- diff(range(ess_stage1_yc))

  ess_stage1_lm <- try(survival::survreg(
    survival::Surv(exp(ess_stage1_y), ess_stage1_status) ~ .,
    data = data.frame(ess_stage1_x), dist = "lognormal"
  ), silent = TRUE)
  ess_stage1_sigma_hat <- if (!inherits(ess_stage1_lm, "try-error") &&
      is.finite(ess_stage1_lm$scale) && ess_stage1_lm$scale > 0)
    ess_stage1_lm$scale else stats::sd(ess_stage1_yc)
  if (!is.finite(ess_stage1_sigma_hat) || ess_stage1_sigma_hat <= 0)
    ess_stage1_sigma_hat <- 1
  if (!is.finite(ess_stage1_range) || ess_stage1_range <= 0)
    ess_stage1_range <- 2 * ess_stage1_sigma_hat

  ess_stage1_cutpoints <- vector("list", ncol(ess_stage1_x))
  for (ess_j in seq_len(ncol(ess_stage1_x))) {
    ess_u <- sort(unique(ess_stage1_x[, ess_j]))
    if (length(ess_u) < 2L)
      stop("LRC-BART requires two observed values for predictor ",
           colnames(ess_stage1_x)[ess_j])
    if (length(ess_u) <= 100L) {
      ess_stage1_cutpoints[[ess_j]] <-
        (ess_u[-1] + ess_u[-length(ess_u)]) / 2
    } else {
      ess_r <- range(ess_u)
      ess_stage1_cutpoints[[ess_j]] <-
        ess_r[1] + diff(ess_r) * seq_len(100L) / 101
    }
  }
  ess_qchi <- stats::qchisq(0.10, 3) / 3
  ess_stage1_hyper <- list(
    H_f = as.integer(ess_H_f), H_g = 0L,
    alpha_f = 0.95, beta_f = 2, alpha_g = ess_alpha_g,
    beta_g = ess_beta_g,
    lambda_f_sq = (ess_stage1_range / (4 * sqrt(ess_H_f)))^2,
    nu_sigma = 3,
    lambda_sigma1 = ess_stage1_sigma_hat^2 * ess_qchi,
    lambda_sigma2 = ess_stage1_sigma_hat^2 * ess_qchi,
    nu0 = ess_nu0, s0_sq = 1, tau1_sq = 1,
    a_w = 1, b_w = 1, update_tau0 = FALSE,
    augmentation = TRUE, p_grow = 0.3, p_prune = 0.3,
    n_min = 5L, n_min_g = 0L,
    sigma1_sq_init = ess_stage1_sigma_hat^2,
    sigma2_sq_init = ess_stage1_sigma_hat^2,
    w_fixed = NA_real_
  )
  ess_stage1_control <- list(
    n_burn = as.integer(ess_n_burn), n_draw = as.integer(ess_n_draw),
    thin = 1L, verbose = FALSE, keep_trees = FALSE,
    keep_train = FALSE, keep_test = TRUE, seed = as.integer(ess_seed)
  )
  ess_stage1_fit <- clrcbart(
    t(ess_stage1_x), ess_stage1_yc,
    rep.int(2L, length(ess_stage1_yc)), ess_stage1_status,
    ess_stage1_censor, t(ess_stage1_xtest),
    ess_stage1_cutpoints, ess_stage1_hyper, ess_stage1_control
  )
  ess_stage1_f_test <- ess_stage1_fit$f_test + ess_stage1_center

  ess_mu_draws <- rowMeans(ess_stage1_f_test)
  ess_V_whole <- stats::var(ess_mu_draws)
  if (!is.finite(ess_V_whole) || ess_V_whole <= 0)
    stop("Stage-1 mean draw has non-positive variance")
  ess_n_blocks <- min(as.integer(ess_n_blocks),
                      floor(length(ess_mu_draws) / 2))
  if (ess_n_blocks < 2L) stop("At least four Stage-1 draws are required")
  ess_block <- cut(seq_along(ess_mu_draws), breaks = ess_n_blocks,
                   labels = FALSE)
  ess_V_raw <- as.numeric(tapply(ess_mu_draws, ess_block, stats::var))
  if (any(!is.finite(ess_V_raw)) || mean(ess_V_raw) <= 0)
    ess_V_raw <- rep(ess_V_whole, ess_n_blocks)
  ess_V_blocks <- ess_V_raw * ess_V_whole / mean(ess_V_raw)

  # No RCT controls are observed, so the RWD residual variance supplies the
  # single-arm sigma_1^2 calibration quantity.
  ess_sigma1_sq <- mean(ess_stage1_fit$sigma2_sq)
  ess_sigma1_blocks <- rep(ess_sigma1_sq, length(ess_V_blocks))

  ess_cg_grid <- rbind(ess_stage1_x, ess_stage1_xtest)
  ess_cg_cutpoints <- vector("list", ncol(ess_cg_grid))
  for (ess_j in seq_len(ncol(ess_cg_grid))) {
    ess_u <- sort(unique(ess_cg_grid[, ess_j]))
    if (length(ess_u) <= 100L) {
      ess_cg_cutpoints[[ess_j]] <-
        (ess_u[-1] + ess_u[-length(ess_u)]) / 2
    } else {
      ess_r <- range(ess_u)
      ess_cg_cutpoints[[ess_j]] <-
        ess_r[1] + diff(ess_r) * seq_len(100L) / 101
    }
  }
  ess_c_g <- cess_cg(
    t(ess_stage1_xtest), ess_cg_cutpoints,
    ess_alpha_g, ess_beta_g, as.integer(ess_n_cg), 12L,
    as.integer(ess_seed + 21L)
  )
  ess_ceiling <- ess_sigma1_sq / ess_V_whole

  set.seed(ess_seed + 31L)
  ess_chi <- stats::rchisq(ess_n_tau, ess_nu0)
  ess_tau0_blocks <- matrix(
    NA_real_, nrow = length(ess_V_blocks), ncol = length(ess_grid)
  )
  for (ess_k in seq_along(ess_grid)) {
    ess_tau0_sq <- ess_nu0 * ess_grid[ess_k] / ess_chi
    ess_tau0_blocks[, ess_k] <- cess_ess_tau0(
      ess_V_blocks, ess_sigma1_blocks, ess_c_g,
      as.integer(ess_H_g), ess_tau0_sq
    )
  }

  ess_requested <- c(
    30, 22.5, 15,
    0.90 * ess_ceiling, 0.50 * ess_ceiling, 0.25 * ess_ceiling
  )
  ess_fraction <- c(1, 0.75, 0.50, 0.90, 0.50, 0.25)
  ess_targets <- data.frame(
    target_name = ess_target_names,
    target = pmin(ess_requested, ess_ceiling),
    requested = ess_requested, fraction = ess_fraction,
    s0_sq = NA_real_, ess_tau0 = NA_real_,
    capped = ess_requested >= ess_ceiling,
    boundary = NA_character_, stringsAsFactors = FALSE
  )
  for (ess_k in seq_along(ess_target_names)) {
    ess_G <- colMeans(ess_tau0_blocks <= ess_targets$target[ess_k])
    if (ess_targets$capped[ess_k]) {
      ess_pick <- 1L
      ess_targets$boundary[ess_k] <- "lower"
    } else {
      ess_ok <- which(ess_G >= ess_q)
      if (!length(ess_ok)) {
        ess_pick <- length(ess_grid)
        ess_targets$boundary[ess_k] <- "upper"
      } else {
        ess_pick <- ess_ok[1]
        if (ess_pick == 1L) ess_targets$boundary[ess_k] <- "lower"
      }
    }
    ess_targets$s0_sq[ess_k] <- ess_grid[ess_pick]
    ess_targets$ess_tau0[ess_k] <- mean(ess_tau0_blocks[, ess_pick])
  }

  ess_calibration <- list(
    data_file = normalizePath(ess_data_file), data_mtime = ess_data_mtime,
    outcome = ess_outcome, data_tag = ess_data_tag,
    settings = ess_settings, H_f = as.integer(ess_H_f),
    H_g = as.integer(ess_H_g), alpha_g = ess_alpha_g,
    beta_g = ess_beta_g, nu0 = ess_nu0, q = ess_q,
    grid = ess_grid, targets = ess_targets,
    s0_sq = stats::setNames(ess_targets$s0_sq, ess_targets$target_name),
    ceiling = ess_ceiling, V_mu_f = ess_V_whole,
    V_profile_f = apply(ess_stage1_f_test, 2, stats::var),
    ess_tau0_blocks = ess_tau0_blocks,
    ess_tau0_grid = colMeans(ess_tau0_blocks),
    sigma1_sq = ess_sigma1_sq,
    c_g = ess_c_g,
    checkpoint = ess_checkpoint
  )
  saveRDS(ess_calibration, ess_checkpoint)
}
# LRC-BART ADDITION END
