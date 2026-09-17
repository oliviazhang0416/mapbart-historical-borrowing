# LRC-BART MODIFICATION START
# Replicate-specific censored-survival single-arm LRC-BART ESS calibration.
#
# This file is sourced from inside lrcBART.R's replicate loop. It follows the
# original project's procedural script style: the calculation is performed
# directly and returned in `ess_calibration`, without private helper functions
# or a second replicate loop.

ess_required <- c("ess_data_file", "ess_H_g", "ess_H_f", "ess_n_burn",
                  "ess_n_draw", "ess_seed", "ess_checkpoint_dir")
for (ess_name in ess_required) {
  if (!exists(ess_name, inherits = TRUE))
    stop("ess_cal.R requires `", ess_name, "` from lrcBART.R")
}
if (!exists("clrcbart", mode = "function", inherits = TRUE) ||
    !exists("cess_cg", mode = "function", inherits = TRUE) ||
    !exists("cess_ess_tau0", mode = "function", inherits = TRUE))
  stop("lrcBART.R must load clrcbart.cpp and cess.cpp before ess_cal.R")
if (!file.exists(ess_data_file)) stop("Data file not found: ", ess_data_file)

ess_alpha_g <- 0.5
ess_beta_g <- 3
ess_nu0 <- 3
ess_q <- 0.95
ess_grid <- 10^seq(-6, 1, length.out = 141)
ess_n_blocks <- 20L
ess_n_tau <- 4000L
ess_n_cg <- 2000L

dir.create(ess_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
ess_stem <- sub("\\.RData$", "", basename(ess_data_file), ignore.case = TRUE)
ess_checkpoint <- file.path(
  ess_checkpoint_dir, paste0("ess_", ess_stem, "_Hg", ess_H_g, ".RData")
)
ess_data_mtime <- file.info(ess_data_file)$mtime
ess_settings <- list(
  H_g = as.integer(ess_H_g), H_f = as.integer(ess_H_f),
  alpha_g = ess_alpha_g, beta_g = ess_beta_g, nu0 = ess_nu0, q = ess_q,
  grid = as.numeric(ess_grid), n_blocks = as.integer(ess_n_blocks),
  n_tau = as.integer(ess_n_tau), n_cg = as.integer(ess_n_cg),
  n_burn = as.integer(ess_n_burn), n_draw = as.integer(ess_n_draw),
  seed = as.integer(ess_seed), augmentation = TRUE
)

ess_use_checkpoint <- FALSE
if (file.exists(ess_checkpoint) &&
    file.info(ess_checkpoint)$mtime >= ess_data_mtime) {
  ess_cached <- readRDS(ess_checkpoint)
  ess_use_checkpoint <- identical(ess_cached$settings, ess_settings) &&
    identical(normalizePath(ess_cached$data_file),
              normalizePath(ess_data_file)) &&
    all(c("100", "f90", "s0min") %in%
          names(ess_cached$s0_sq))
}

if (ess_use_checkpoint) {
  ess_calibration <- ess_cached
} else {
  ess_dat <- readRDS(ess_data_file)
  ess_x_names <- grep("^X[0-9]+$", colnames(ess_dat$X), value = TRUE)
  if (!length(ess_x_names))
    stop("No X1, X2, ... predictors found in ", ess_data_file)
  ess_x_all <- as.matrix(ess_dat$X[, ess_x_names, drop = FALSE])
  storage.mode(ess_x_all) <- "double"
  ess_d <- ess_dat$X[, "D"]
  ess_z <- ess_dat$X[, "Z"]
  ess_idx_rwd <- ess_d == 0 & ess_z == 0
  ess_idx_rct <- ess_d == 1
  if (!any(ess_idx_rwd) || !any(ess_idx_rct))
    stop("Survival single-arm calibration requires RWD controls and RCT profiles")
  if (any(ess_d == 1 & ess_z == 0))
    stop("Survival single-arm calibration received RCT controls")

  # Stage 1: plain censored AFT-BART on the RWD controls, evaluated at all
  # RCT profiles. Status 0 invokes latent log-time augmentation.
  ess_stage1_y <- log(as.numeric(ess_dat$y[ess_idx_rwd]))
  ess_stage1_status <- as.integer(ess_dat$delta[ess_idx_rwd])
  ess_stage1_x <- ess_x_all[ess_idx_rwd, , drop = FALSE]
  ess_stage1_xtest <- ess_x_all[ess_idx_rct, , drop = FALSE]
  ess_stage1_center <- mean(ess_stage1_y)
  ess_stage1_yc <- ess_stage1_y - ess_stage1_center
  ess_stage1_censor <- ess_stage1_y - ess_stage1_center
  ess_stage1_range <- diff(range(ess_stage1_yc))
  ess_stage1_sigma_hat <- 1
  if (length(ess_stage1_yc) > ncol(ess_stage1_x) + 5L) {
    ess_stage1_lm <- try(survival::survreg(
      survival::Surv(exp(ess_stage1_y), ess_stage1_status) ~ .,
      data = data.frame(ess_stage1_x), dist = "lognormal"
    ), silent = TRUE)
    if (!inherits(ess_stage1_lm, "try-error") &&
        is.finite(ess_stage1_lm$scale) && ess_stage1_lm$scale > 0)
      ess_stage1_sigma_hat <- ess_stage1_lm$scale
  } else if (is.finite(stats::sd(ess_stage1_yc)) &&
             stats::sd(ess_stage1_yc) > 0) {
    ess_stage1_sigma_hat <- stats::sd(ess_stage1_yc)
  }
  if (!is.finite(ess_stage1_range) || ess_stage1_range <= 0)
    ess_stage1_range <- 2 * ess_stage1_sigma_hat
  ess_stage1_cutpoints <- vector("list", ncol(ess_stage1_x))
  for (ess_j in seq_len(ncol(ess_stage1_x))) {
    ess_u <- sort(unique(ess_stage1_x[, ess_j]))
    if (length(ess_u) < 2L)
      stop("LRC-BART requires at least two observed values for predictor ", ess_j)
    if (length(ess_u) <= 100L) {
      ess_stage1_cutpoints[[ess_j]] <-
        (ess_u[-1] + ess_u[-length(ess_u)]) / 2
    } else {
      ess_r <- range(ess_u)
      ess_stage1_cutpoints[[ess_j]] <-
        ess_r[1] + diff(ess_r) * seq_len(100L) / 101
    }
  }
  ess_stage1_hyper <- list(
    H_f = as.integer(ess_H_f), H_g = 0L,
    alpha_f = 0.95, beta_f = 2, alpha_g = 0.5, beta_g = 3,
    lambda_f_sq = (ess_stage1_range / (4 * sqrt(ess_H_f)))^2,
    nu_sigma = 3,
    lambda_sigma1 = ess_stage1_sigma_hat^2 * stats::qchisq(0.10, 3) / 3,
    lambda_sigma2 = ess_stage1_sigma_hat^2 * stats::qchisq(0.10, 3) / 3,
    nu0 = 3, s0_sq = 1, tau1_sq = 1,
    a_w = 1, b_w = 1, update_tau0 = FALSE,
    augmentation = TRUE, p_grow = 0.3, p_prune = 0.3,
    n_min = 5L, n_min_g = 0L,
    sigma1_sq_init = ess_stage1_sigma_hat^2,
    sigma2_sq_init = ess_stage1_sigma_hat^2,
    w_fixed = NA_real_
  )
  ess_stage1_control <- list(
    n_burn = as.integer(ess_n_burn), n_draw = as.integer(ess_n_draw),
    thin = 1L, verbose = FALSE, keep_trees = FALSE, keep_train = FALSE,
    keep_test = TRUE, seed = as.integer(ess_seed)
  )
  ess_stage1_fit <- clrcbart(
    t(ess_stage1_x), ess_stage1_yc,
    rep.int(2L, length(ess_stage1_yc)), ess_stage1_status,
    ess_stage1_censor,
    t(ess_stage1_xtest), ess_stage1_cutpoints,
    ess_stage1_hyper, ess_stage1_control
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
  # No RCT controls exist. In single-arm mode the standalone sampler uses the
  # RWD residual variance for both sigma_1^2 and sigma_2^2.
  ess_sigma1_sq <- mean(ess_stage1_fit$sigma2_sq)
  ess_sigma1_blocks <- rep(ess_sigma1_sq, length(ess_V_blocks))

  ess_region_rct <- ess_dat$region_rct
  ess_V_region <- c(region = NA_real_, outside = NA_real_)
  if (!is.null(ess_region_rct) &&
      length(ess_region_rct) == ncol(ess_stage1_f_test) &&
      any(ess_region_rct) && any(!ess_region_rct)) {
    ess_V_region["region"] <- stats::var(rowMeans(
      ess_stage1_f_test[, ess_region_rct, drop = FALSE]
    ))
    ess_V_region["outside"] <- stats::var(rowMeans(
      ess_stage1_f_test[, !ess_region_rct, drop = FALSE]
    ))
  }

  ess_cg_grid <- rbind(ess_x_all[ess_idx_rwd, , drop = FALSE],
                       ess_x_all[ess_idx_rct, , drop = FALSE])
  ess_cg_cutpoints <- vector("list", ncol(ess_cg_grid))
  for (ess_j in seq_len(ncol(ess_cg_grid))) {
    ess_u <- sort(unique(ess_cg_grid[, ess_j]))
    if (length(ess_u) < 2L)
      stop("LRC-BART requires at least two observed values for predictor ", ess_j)
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
    t(ess_x_all[ess_idx_rct, , drop = FALSE]), ess_cg_cutpoints,
    ess_alpha_g, ess_beta_g, as.integer(ess_n_cg), 12L,
    as.integer(ess_seed + 21L)
  )
  ess_ceiling <- ess_sigma1_sq / ess_V_whole

  set.seed(ess_seed + 31L)
  ess_chi <- stats::rchisq(ess_n_tau, ess_nu0)
  ess_tau0_blocks <- matrix(NA_real_, nrow = length(ess_V_blocks),
                            ncol = length(ess_grid))
  for (ess_k in seq_along(ess_grid)) {
    ess_tau0_sq <- ess_nu0 * ess_grid[ess_k] / ess_chi
    ess_tau0_blocks[, ess_k] <- cess_ess_tau0(
      ess_V_blocks, ess_sigma1_blocks, ess_c_g,
      as.integer(ess_H_g), ess_tau0_sq
    )
  }

  ess_target_names <- c("100", "f90", "s0min")
  # s0min is the fixed near-complete-pooling benchmark: it uses the smallest
  # supplied s0^2 value (1e-6) instead of applying the target probability rule.
  ess_requested <- c(100, 0.90 * ess_ceiling, NA_real_)
  ess_fraction <- c(NA_real_, 0.90, NA_real_)
  ess_targets <- data.frame(
    target_name = ess_target_names,
    target = NA_real_, requested = ess_requested,
    fraction = ess_fraction, s0_sq = NA_real_, ess_tau0 = NA_real_,
    capped = FALSE, boundary = NA_character_, stringsAsFactors = FALSE
  )
  for (ess_k in seq_along(ess_target_names)) {
    if (ess_target_names[ess_k] == "s0min") {
      ess_pick <- 1L
      ess_targets$target[ess_k] <- NA_real_
      ess_targets$capped[ess_k] <- FALSE
      ess_targets$boundary[ess_k] <- "lower"
    } else {
      ess_targets$target[ess_k] <- min(ess_requested[ess_k], ess_ceiling)
      ess_targets$capped[ess_k] <- ess_requested[ess_k] >= ess_ceiling
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
    }
    ess_targets$s0_sq[ess_k] <- ess_grid[ess_pick]
    ess_targets$ess_tau0[ess_k] <- mean(ess_tau0_blocks[, ess_pick])
  }
  ess_s0_sq <- stats::setNames(ess_targets$s0_sq,
                               ess_targets$target_name)

  ess_calibration <- list(
    data_file = normalizePath(ess_data_file), data_mtime = ess_data_mtime,
    settings = ess_settings,
    H_f = as.integer(ess_H_f), H_g = as.integer(ess_H_g),
    alpha_g = ess_alpha_g, beta_g = ess_beta_g,
    nu0 = ess_nu0, q = ess_q,
    grid = ess_grid, targets = ess_targets, s0_sq = ess_s0_sq,
    ceiling = ess_ceiling, V_mu_f = ess_V_whole,
    V_mu_f_region = ess_V_region,
    V_blocks = ess_V_blocks, V_blocks_raw = ess_V_raw,
    ess_tau0_blocks = ess_tau0_blocks,
    ess_tau0_grid = colMeans(ess_tau0_blocks),
    sigma1_sq = ess_sigma1_sq,
    sigma1_sq_blocks = ess_sigma1_blocks, c_g = ess_c_g,
    stage1_mu_draws = ess_mu_draws,
    stage1_sigma2_sq = ess_stage1_fit$sigma2_sq,
    checkpoint = ess_checkpoint
  )
  saveRDS(ess_calibration, ess_checkpoint)
}
# LRC-BART MODIFICATION END
