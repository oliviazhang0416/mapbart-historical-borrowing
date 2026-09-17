# LRC-BART MODIFICATION START
# Replicate-specific MAP ESS calibration.
#
# MAP.R sources this file inside its inherited target/replicate loop. The
# calculation is procedural and returns `map_calibration`; it does not define
# a second workflow function or its own replicate loop.

map_required <- c("map_data_file", "map_target_Ns", "map_seed",
                  "map_checkpoint_dir")
for (map_name in map_required) {
  if (!exists(map_name, inherits = TRUE))
    stop("ess_cal_map.R requires `", map_name, "` from MAP.R")
}
if (!requireNamespace("RBesT", quietly = TRUE))
  stop("Package 'RBesT' is required for MAP ESS calibration")
if (!file.exists(map_data_file)) stop("Data file not found: ", map_data_file)

map_target_Ns <- as.integer(map_target_Ns)
map_S_grid <- seq(0.05, 0.45, by = 0.05)
map_tilde_nu <- 3
map_beta_prior <- 10
map_sigma_ref <- NULL
if (!length(map_target_Ns) || any(map_target_Ns <= 0L))
  stop("map_target_Ns must contain positive integers")

dir.create(map_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
map_checkpoint_file <- file.path(
  map_checkpoint_dir,
  paste0("map_ess_", tools::file_path_sans_ext(basename(map_data_file)),
         "_checkpoint.RData")
)
map_settings <- list(
  target_Ns = map_target_Ns, S_grid = map_S_grid,
  tilde_nu = map_tilde_nu, beta_prior = map_beta_prior,
  sigma_ref = map_sigma_ref, seed = as.integer(map_seed)
)

map_use_checkpoint <- FALSE
if (file.exists(map_checkpoint_file)) {
  map_cached <- readRDS(map_checkpoint_file)
  map_use_checkpoint <- identical(map_cached$settings, map_settings) &&
    identical(normalizePath(map_cached$data_file),
              normalizePath(map_data_file)) &&
    file.info(map_checkpoint_file)$mtime >= file.info(map_data_file)$mtime
}

if (map_use_checkpoint) {
  map_calibration <- map_cached
} else {
  map_dat <- readRDS(map_data_file)
  map_idx <- map_dat$X[, "D"] == 0 & map_dat$X[, "Z"] == 0
  map_y <- as.numeric(map_dat$y[map_idx])
  if (length(map_y) < 2L || any(!is.finite(map_y)))
    stop("MAP calibration needs at least two finite RWD-control outcomes")
  map_n_RWD <- length(map_y)
  map_ybar <- mean(map_y)
  map_sd_y <- stats::sd(map_y)
  map_sigma_ref_used <- if (is.null(map_sigma_ref))
    map_sd_y else as.numeric(map_sigma_ref)
  map_historical_data <- data.frame(
    y = map_ybar, n = map_n_RWD,
    y.se = map_sd_y / sqrt(map_n_RWD), study = factor(1L)
  )

  # These are RBesT posterior-sampler settings for one dataset. They are
  # independent of the simulation replicate count.
  map_old_options <- options(
    RBesT.MC.chains = 1,
    RBesT.MC.warmup = 500,
    RBesT.MC.iter = 2000,
    RBesT.MC.thin = 1,
    RBesT.MC.control = list(
      adapt_delta = 0.95, max_treedepth = 10, stepsize = 0.1
    ),
    mc.cores = 1L
  )
  set.seed(as.integer(map_seed))
  map_ESS <- numeric(length(map_S_grid))
  for (map_k in seq_along(map_S_grid)) {
    map_prior <- RBesT::gMAP(
      formula = cbind(y, y.se) ~ 1 | study,
      data = map_historical_data,
      tau.dist = "InvGamma",
      tau.prior = c(map_tilde_nu / 2,
                    map_tilde_nu * map_S_grid[map_k] / 2),
      beta.prior = map_beta_prior,
      family = gaussian
    )
    map_mixture <- RBesT::automixfit(map_prior, Nc = 1:3, thresh = 0)
    map_ESS[map_k] <- as.numeric(
      RBesT::ess(map_mixture, sigma = map_sigma_ref_used)
    )
  }
  options(map_old_options)

  map_targets <- data.frame(
    target_N = map_target_Ns,
    s0_sq = NA_real_, ess_prior = NA_real_
  )
  for (map_k in seq_along(map_target_Ns)) {
    map_eligible <- which(map_ESS <= map_target_Ns[map_k])
    if (!length(map_eligible))
      stop("No MAP s^2 grid value attains ESS <= ", map_target_Ns[map_k],
           " for ", basename(map_data_file), "; extend the grid upward")
    map_pick <- map_eligible[1]
    map_targets$s0_sq[map_k] <- map_S_grid[map_pick]
    map_targets$ess_prior[map_k] <- map_ESS[map_pick]
  }

  map_calibration <- list(
    data_file = normalizePath(map_data_file), settings = map_settings,
    n_RWD = map_n_RWD, ybar = map_ybar, sd_y = map_sd_y,
    sigma_ref = map_sigma_ref_used,
    S_grid = map_S_grid, ESS = map_ESS, targets = map_targets,
    checkpoint_file = map_checkpoint_file
  )
  saveRDS(map_calibration, map_checkpoint_file)
}
# LRC-BART MODIFICATION END
