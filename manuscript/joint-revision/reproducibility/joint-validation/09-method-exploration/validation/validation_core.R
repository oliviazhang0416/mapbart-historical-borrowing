# Shared design, data construction, calibration, fitting, and cache helpers for
# the frozen validation queue.  This file has no top-level numerical run.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")

REVISION_ROOT <- normalizePath(".", mustWork = TRUE)
VALIDATION_ROOT <- file.path(REVISION_ROOT, "09-method-exploration", "validation")
JOINT_LIBRARY <- file.path(REVISION_ROOT, "09-method-exploration", "joint-model", "private-library")
.libPaths(c(normalizePath(JOINT_LIBRARY, mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(lrcbart))

source(file.path(REVISION_ROOT, "01-code", "R", "dgp.R"))
source(file.path(REVISION_ROOT, "08-interaction-pilot", "20260918", "code", "scenarios.R"))

XCOLS <- paste0("X", 1:10)
METHODS <- c("JOINT-LRC-01", "JOINT-SOURCE-01", "SEPARATED-SOURCE-01")
CONDITIONS <- c("compatible", "one", "two", "heterogeneous", "null")
N_BURN_DEFAULT <- 30000L
N_DRAW_DEFAULT <- 12000L
H_F_DEFAULT <- 50L
H_G_LRC <- 5L

hash_files <- function(paths) {
  paths <- paths[file.exists(paths)]
  paths <- paths[!file.info(paths)$isdir]
  if (!length(paths)) return(NA_character_)
  h <- tools::md5sum(paths)
  paste(names(h), as.character(h), sep = "=", collapse = ";")
}

hash_object <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

CODE_HASH <- hash_files(c(
  file.path(VALIDATION_ROOT, "validation_core.R"),
  file.path(VALIDATION_ROOT, "run_validation.R"),
  file.path(VALIDATION_ROOT, "summarize_validation.R"),
  file.path(VALIDATION_ROOT, "checks.R"),
  file.path(REVISION_ROOT, "01-code", "R", "dgp.R"),
  file.path(REVISION_ROOT, "08-interaction-pilot", "20260918", "code", "scenarios.R"),
  list.files(file.path(REVISION_ROOT, "09-method-exploration", "joint-model", "lrcbart-joint"),
             recursive = TRUE, full.names = TRUE)))
LIBRARY_HASH <- hash_files(list.files(file.path(JOINT_LIBRARY, "lrcbart"),
                                      recursive = TRUE, full.names = TRUE))
PRIOR_CONFIG <- list(H_f = H_F_DEFAULT, H_g_lrc = H_G_LRC, alpha_g = 0.5,
                     beta_g = 3, tau_prior_var = 100,
                     separated_source_variances = TRUE,
                     calibration = "original training-only routine; N_target=100 is retained for LRC calibration only")

method_config_hash <- function(method, n_burn, n_draw) {
  digest::digest(c(PRIOR_CONFIG, list(method = method, n_burn = n_burn,
                                      n_draw = n_draw,
                                      variance_model = "source-specific")),
                 algo = "sha256", serialize = TRUE)
}

condition_spec <- function(condition) {
  if (condition == "compatible") return(list(region = "ONE50", delta_rwd = 0, seed_base = 92026000L))
  if (condition == "one") return(list(region = "ONE50", delta_rwd = 2, seed_base = 92026000L))
  if (condition == "two") return(list(region = "XOR50", delta_rwd = 2, seed_base = 92026000L))
  if (condition == "heterogeneous") return(list(region = "XOR50", delta_rwd = 2, seed_base = 92126000L))
  if (condition == "null") return(list(region = "XOR50", delta_rwd = 2, seed_base = 92226000L))
  stop("unknown condition: ", condition)
}

region_mask <- function(x, region) {
  x <- as.matrix(x)
  switch(region,
         ONE50 = x[, 5] > 2,
         XOR50 = (x[, 5] > 2) == (x[, 7] > 2),
         stop("unknown validation region: ", region))
}

dataset_seed <- function(condition, rep, mode = tolower(Sys.getenv("VALIDATION_MODE", "validation"))) {
  if (mode == "smoke") {
    smoke_base <- switch(condition,
                         compatible = 91990000L, one = 91990000L, two = 91990000L,
                         heterogeneous = 91991000L, null = 91992000L,
                         stop("unknown smoke condition: ", condition))
    return(smoke_base + as.integer(rep))
  }
  condition_spec(condition)$seed_base + as.integer(rep)
}

make_job_table <- function() {
  rows <- list(); k <- 0L
  for (condition in CONDITIONS) for (rep in 1:40) for (chain in 1:3) {
    for (method in METHODS) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        job_id = k, condition = condition, rep = rep, chain = chain,
        method = method, dataset_seed = dataset_seed(condition, rep, mode = "validation"),
        fit_seed = 94000000L + 10L * k, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

assert_same_manifest <- function(path, expected) {
  if (!file.exists(path)) {
    write.csv(expected, path, row.names = FALSE)
    return(invisible(TRUE))
  }
  got <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(identical(names(got), names(expected)), nrow(got) == nrow(expected))
  for (nm in names(expected)) stopifnot(identical(as.character(got[[nm]]), as.character(expected[[nm]])))
  invisible(TRUE)
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, tmp)
  stopifnot(file.rename(tmp, path))
  invisible(path)
}

generate_dataset <- function(condition, rep) {
  sp <- condition_spec(condition)
  seed <- dataset_seed(condition, rep)
  d <- gen_data(4, outcome = "gaussian", n_rct = 300, n_rwd = 300,
                region = sp$region, delta_rwd = sp$delta_rwd,
                seed = seed, n_mc = 0)
  stopifnot(nrow(d$x_rct) == 300L, nrow(d$rwd_ctrl) == 300L,
            ncol(d$x_rct) == 10L, d$constants$sigma_rct == 1.5,
            d$constants$sigma_rwd == 1.8)
  region_rct <- as.logical(d$region_rct)
  true_g <- if (condition == "compatible") rep(0, 300) else -2 * (!region_rct)
  true_tau <- rep(1, 300)
  if (condition == "heterogeneous") {
    true_tau <- 1 + 0.5 * (d$x_rct[, 5] - mean(d$x_rct[, 5]))
    ntrt <- nrow(d$rct_trt)
    d$rct_trt$y <- d$rct_trt$y + true_tau[seq_len(ntrt)] - 1
    d$mu1_rct <- d$mu0_rct + true_tau
    d$truth <- list(ate = mean(true_tau), trt = mean(d$mu1_rct), ctrl = mean(d$mu0_rct))
  } else if (condition == "null") {
    d$rct_trt$y <- d$rct_trt$y - 1
    true_tau <- rep(0, 300)
    d$mu1_rct <- d$mu0_rct
    d$truth <- list(ate = 0, trt = mean(d$mu1_rct), ctrl = mean(d$mu0_rct))
  }
  stopifnot(abs(mean(true_tau) - unname(d$truth$ate)) < 1e-12)
  xc <- as.matrix(d$rct_ctrl[, XCOLS, drop = FALSE])
  xw <- as.matrix(d$rwd_ctrl[, XCOLS, drop = FALSE])
  x <- rbind(xc, xw)
  y <- c(d$rct_ctrl$y, d$rwd_ctrl$y)
  src <- rep(c(1L, 2L), c(nrow(xc), nrow(xw)))
  rct <- rbind(d$rct_trt, d$rct_ctrl)
  x_joint <- rbind(as.matrix(rct[, XCOLS, drop = FALSE]), xw)
  y_joint <- c(as.numeric(rct$y), as.numeric(d$rwd_ctrl$y))
  src_joint <- c(rep(1L, 300L), rep(2L, 300L))
  A_joint <- c(as.numeric(d$z_rct), rep(0, 300L))
  stopifnot(max(abs(x_joint[seq_len(300), , drop = FALSE] - d$x_rct)) < 1e-12,
            identical(as.integer(A_joint[seq_len(300)]), as.integer(d$z_rct)))
  list(condition = condition, rep = as.integer(rep), seed = seed, d = d,
       x = x, y = y, src = src, x_joint = x_joint, y_joint = y_joint,
       src_joint = src_joint, A_joint = A_joint,
       region_rct = region_rct, true_g = true_g,
       true_tau = true_tau, truth_ate = unname(d$truth$ate),
       source_counts = as.integer(table(factor(src, levels = 1:2))),
       source_counts_joint = as.integer(table(factor(src_joint, levels = 1:2))),
       region = sp$region, delta_rwd = sp$delta_rwd)
}

assert_constant_family <- function(family) {
  stopifnot(all(c("compatible", "one", "two") %in% names(family)))
  a <- family$compatible$d; b <- family$one$d; c <- family$two$d
  for (nm in c("rct_trt", "rct_ctrl", "x_rct", "z_rct")) {
    stopifnot(isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 0)),
              isTRUE(all.equal(a[[nm]], c[[nm]], tolerance = 0)))
  }
  one_rwd <- region_mask(b$rwd_ctrl[, XCOLS, drop = FALSE], "ONE50")
  two_rwd <- region_mask(c$rwd_ctrl[, XCOLS, drop = FALSE], "XOR50")
  stopifnot(max(abs((b$rwd_ctrl$y - a$rwd_ctrl$y) - 2 * (!one_rwd))) < 1e-12,
            max(abs((c$rwd_ctrl$y - a$rwd_ctrl$y) - 2 * (!two_rwd))) < 1e-12,
            all(family$compatible$true_g == 0),
            all(family$one$true_g == -2 * (!family$one$region_rct)),
            all(family$two$true_g == -2 * (!family$two$region_rct)))
  invisible(TRUE)
}

calibration_for <- function(dat, seed, n_burn = 1000L, n_draw = 1000L) {
  xc <- as.matrix(dat$d$rct_ctrl[, XCOLS, drop = FALSE])
  xw <- as.matrix(dat$d$rwd_ctrl[, XCOLS, drop = FALSE])
  xt <- dat$d$x_rct
  np <- bart_fit(dat$d$rct_ctrl$y, xc, H_f = H_F_DEFAULT,
                 n_burn = n_burn, n_draw = n_draw, seed = seed + 100001L)
  stage <- bart_fit(dat$d$rwd_ctrl$y, xw, H_f = H_F_DEFAULT,
                    n_burn = n_burn, n_draw = n_draw, seed = seed + 200001L)
  cal <- lrc_ess_calibrate(dat$d$rwd_ctrl$y, xw, xt, N_target = 100,
                            fit_stage1 = stage, sigma1_sq = mean(np$sigma1_sq),
                            H_f = H_F_DEFAULT, H_g = H_G_LRC,
                            alpha_g = 0.5, beta_g = 3,
                            n_burn = n_burn, n_draw = n_draw,
                            seed = seed + 300001L, warn = FALSE)
  cal$fit_stage1 <- NULL
  list(cal = cal, seeds = c(rct_stage = seed + 100001L,
                            rwd_stage = seed + 200001L,
                            calibration = seed + 300001L),
       n_burn = n_burn, n_draw = n_draw)
}

dataset_path <- function(mode_dir, condition, rep) {
  file.path(mode_dir, "data", sprintf("%s_r%02d.rds", condition, rep))
}

ensure_dataset <- function(mode_dir, condition, rep, calib_burn = 1000L,
                            calib_draw = 1000L) {
  path <- dataset_path(mode_dir, condition, rep)
  key <- paste(condition, rep, dataset_seed(condition, rep), calib_burn, calib_draw, sep = "|")
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$cache_key, key)) stop("dataset cache mismatch: ", path)
    if (!identical(old$code_hash, CODE_HASH) || !identical(old$library_hash, LIBRARY_HASH))
      stop("dataset code/library hash mismatch: ", path)
    return(old)
  }
  dat <- generate_dataset(condition, rep)
  cal <- calibration_for(dat, dat$seed, calib_burn, calib_draw)
  data_sha256 <- hash_object(dat)
  out <- list(cache_key = key, data_sha256 = data_sha256, data_hash = data_sha256,
              code_hash = CODE_HASH,
              library_hash = LIBRARY_HASH, data = dat, calibration = cal,
              created = format(Sys.time(), tz = "UTC"))
  atomic_save_rds(out, path)
  out
}

ensure_data_for_jobs <- function(mode_dir, jobs) {
  pairs <- unique(jobs[c("condition", "rep")])
  for (i in seq_len(nrow(pairs))) {
    condition <- pairs$condition[i]; rep <- pairs$rep[i]
    if (condition %in% c("compatible", "one", "two")) {
      fam <- lapply(c("compatible", "one", "two"), function(cc)
        ensure_dataset(mode_dir, cc, rep))
      names(fam) <- c("compatible", "one", "two")
      assert_constant_family(lapply(fam, `[[`, "data"))
    } else {
      ensure_dataset(mode_dir, condition, rep)
    }
    message("Data/calibration ready: ", condition, " rep=", rep)
  }
  invisible(TRUE)
}

load_dataset <- function(mode_dir, condition, rep) {
  p <- dataset_path(mode_dir, condition, rep)
  if (!file.exists(p)) stop("dataset cache missing; call ensure_data_for_jobs first: ", p)
  readRDS(p)
}

compact_fit <- function(fit, keep_tau = FALSE) {
  tau <- if (keep_tau) fit$tau else NULL
  fit$draws <- list(f = NULL, g = NULL, tau = tau)
  fit$f_train <- NULL; fit$g_train <- NULL
  fit$calib <- NULL; fit$x_train <- NULL
  fit$y <- NULL; fit$status <- NULL
  fit
}

fit_joint_method <- function(dat, cal, method, seed, n_burn, n_draw) {
  xfit <- dat$x_joint
  source_mode <- method == "JOINT-SOURCE-01"
  if (source_mode) xfit <- cbind(xfit, source_rct = as.numeric(dat$src_joint == 1L))
  fit <- if (method == "JOINT-LRC-01") {
    joint_lrc_bart(dat$y_joint, xfit, dat$src_joint, treatment = dat$A_joint,
                   H_f = H_F_DEFAULT, H_g = H_G_LRC, alpha_g = 0.5, beta_g = 3,
                   s0_sq = cal$cal$s0_sq, g_sweeps = 5,
                   calib = cal$cal, n_burn = n_burn, n_draw = n_draw,
                   n_min = 5, seed = seed, keep_trees = TRUE, keep_train = TRUE)
  } else {
    joint_source_bart(dat$y_joint, dat$x_joint, dat$src_joint,
                      treatment = dat$A_joint,
                      H_f = H_F_DEFAULT, s0_sq = cal$cal$s0_sq,
                      calib = cal$cal, n_burn = n_burn, n_draw = n_draw,
                      n_min = 5, seed = seed, keep_trees = TRUE, keep_train = TRUE)
  }
  x1 <- if (source_mode) cbind(dat$d$x_rct, source_rct = 1) else dat$d$x_rct
  x0 <- if (source_mode) cbind(dat$d$x_rct, source_rct = 0) else dat$d$x_rct
  pred1 <- predict(fit, x1, arm = "rct_control")
  pred0 <- predict(fit, x0, arm = "rwd_control")
  trt <- rowMeans(predict(fit, x1, treatment = rep(1, 300)))
  compact <- compact_fit(fit, keep_tau = TRUE)
  list(fit = compact, pred1 = pred1, pred0 = pred0, trt_mean = trt,
       sigma1_sq = fit$sigma1_sq, sigma2_sq = fit$sigma2_sq)
}

fit_separated_source <- function(dat, cal, seed, n_burn, n_draw) {
  xt <- dat$d$x_rct
  trt_fit <- bart_fit(dat$d$rct_trt$y,
                      as.matrix(dat$d$rct_trt[, XCOLS, drop = FALSE]),
                      H_f = H_F_DEFAULT, n_burn = n_burn, n_draw = n_draw,
                      seed = seed + 1L, keep_trees = TRUE, keep_train = TRUE)
  xctrl <- cbind(dat$x, source_rct = as.numeric(dat$src == 1L))
  ctrl_fit <- lrc_bart(dat$y, xctrl, dat$src, H_f = H_F_DEFAULT, H_g = 0,
                       s0_sq = cal$cal$s0_sq, calib = cal$cal,
                       n_burn = n_burn, n_draw = n_draw, seed = seed + 2L,
                       n_min = 5, keep_trees = TRUE, keep_train = TRUE)
  pred1 <- predict(ctrl_fit, cbind(xt, source_rct = 1), arm = "rct_control")
  pred0 <- predict(ctrl_fit, cbind(xt, source_rct = 0), arm = "rwd_control")
  trt <- rowMeans(predict(trt_fit, xt))
  list(fit = list(trt = compact_fit(trt_fit), ctrl = compact_fit(ctrl_fit)),
       pred1 = pred1, pred0 = pred0, trt_mean = trt,
       sigma1_sq = ctrl_fit$sigma1_sq, sigma2_sq = ctrl_fit$sigma2_sq)
}

fit_validation_job <- function(dat, cal, method, seed, n_burn, n_draw) {
  ans <- if (method %in% c("JOINT-LRC-01", "JOINT-SOURCE-01"))
    fit_joint_method(dat, cal, method, seed, n_burn, n_draw) else
    fit_separated_source(dat, cal, seed, n_burn, n_draw)
  pred1 <- ans$pred1; pred0 <- ans$pred0
  stopifnot(nrow(pred1) == n_draw, ncol(pred1) == 300L,
            nrow(pred0) == n_draw, ncol(pred0) == 300L)
  ctrl <- rowMeans(pred1); ate <- ans$trt_mean - ctrl
  if (method %in% c("JOINT-LRC-01", "JOINT-SOURCE-01"))
    stopifnot(max(abs(ate - ans$fit$tau)) < 1e-10)
  gdraw <- pred1 - pred0
  R <- dat$region_rct
  gcontrast <- rowMeans(gdraw[, R, drop = FALSE]) - rowMeans(gdraw[, !R, drop = FALSE])
  trace <- cbind(trt_mean = ans$trt_mean, control_mean = ctrl,
                 ate = ate, g_contrast = gcontrast)
  colnames(trace) <- c("trt_mean", "control_mean", "ate", "g_contrast")
  ans$pred1 <- NULL; ans$pred0 <- NULL
  list(method = method, condition = dat$condition, rep = dat$rep,
       seed = seed, n_burn = n_burn, n_draw = n_draw, trace = trace,
       ate = ate, control_mean = ctrl, trt_mean = ans$trt_mean,
       control_surface = colMeans(pred1), external_surface = colMeans(pred0),
       g_mean = colMeans(gdraw), g_contrast = gcontrast,
       sigma1_sq = ans$sigma1_sq, sigma2_sq = ans$sigma2_sq,
       truth = dat$truth_ate, truth_control = as.numeric(dat$d$mu0_rct),
       true_g = dat$true_g, region_rct = R, source_counts = dat$source_counts,
       source_counts_joint = dat$source_counts_joint,
       source = dat$src_joint, source_used = dat$src_joint,
       fit_source_arrays = if (method == "SEPARATED-SOURCE-01")
         list(control = dat$src, treatment = rep(1L, nrow(dat$d$rct_trt))) else
         list(joint = dat$src_joint),
       controls_source = dat$src, treatment = dat$A_joint,
       fit = ans$fit, dataset_seed = dat$seed,
       data_sha256 = dat$data_sha256, data_hash = dat$data_sha256,
       code_hash = CODE_HASH, library_hash = LIBRARY_HASH,
       prior_config = PRIOR_CONFIG,
       source_variance_config = "separate sigma1_sq and sigma2_sq; source codes preserved",
       calibration_seeds = cal$seeds, data_cache_key = paste(dat$condition, dat$rep, dat$seed, sep = "|"))
}
