rm(list=ls())

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
# devtools::install_github("olssol/psrwe")
library(miceadds)
library(dplyr)
# LRC-BART MODIFICATION START
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-gaussian")
n_replicates <- 100L
# JOINT LRC-BART MODIFICATION START
dir.create(file.path(projectDir, "res"), recursive = TRUE, showWarnings = FALSE)
# JOINT LRC-BART MODIFICATION END
data_folder <- "data"
# psrwe (Wang et al., GPL >= 3) is vendored at the repository root so this
# comparator runs from a clean clone. PSRWE_DIR overrides it if needed.
psrweDir <- Sys.getenv("PSRWE_DIR", unset = file.path(.lrcRoot, "psrwe"))
if (!dir.exists(psrweDir))
  stop("psrwe not found at ", psrweDir,
       ". Expected the vendored copy at <repo>/psrwe, or set PSRWE_DIR.")
# LRC-BART MODIFICATION END

old_wd <- getwd()
# LRC-BART MODIFICATION START
setwd(psrweDir)
# LRC-BART MODIFICATION END
source.all("R/")
setwd(old_wd)

# JOINT LRC-BART MODIFICATION START
sc <- 3
variant <- ""
# JOINT LRC-BART MODIFICATION END
# LRC-BART MODIFICATION START
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(projectDir, data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(projectDir, data_folder))
# LRC-BART MODIFICATION END
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# LRC-BART ADDITION START
delta_rwd <- 0
region <- "none"
# The scenario filename tag is assembled inline after `cor` is set.
# LRC-BART ADDITION END

threshold <- 0.95

# LRC-BART MODIFICATION START
# Strength of unmeasured confounding
# JOINT LRC-BART MODIFICATION START
if (sc %in% c(1, 2)){
# JOINT LRC-BART MODIFICATION END
  cor <- 1
}
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[3]
}

# Validate the configured replicate count against this scenario
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
data_prefix <- paste0("data_p", p_obs, "_", scenario_id, "_", hypo, "_")
replicate_files <- Sys.glob(file.path(projectDir, data_folder,
                                      paste0(data_prefix, "*.RData")))
if (length(replicate_files) < n_replicates)
  stop("Expected ", n_replicates, " simulated data files for ", scenario_id,
       "; found ", length(replicate_files))
# LRC-BART MODIFICATION END

# LRC-BART MODIFICATION START
# Auto-load calibrated threshold when running under alternative
if (hypo == "alternative") {
  threshold_file <- file.path(projectDir, "res",
    paste0("PSCL_p", p_obs, "_", scenario_id, "_threshold.RData"))

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}
# LRC-BART MODIFICATION END

set.seed(6)
seed <- sample(1:10000, n_replicates*4, replace = F)
ee <- 1

res <- data.frame(c = rep(cor, each = n_replicates))
res$bias <- NA
res$sd <- NA
res$rmse <- NA
res$w1distance <- NA
res$w2distance <- NA
res$ci <- NA
res$coverage <- NA
res$tp_calibrated <- NA
res$tp <- NA
res$fp <- NA
res$pehe_subj <- NA  # Not applicable for PSCL
res$bias_subj <- NA  # Not applicable for PSCL
res$bias.trt.sigma <- NA  # Not applicable for PSCL
res$sd.trt.sigma <- NA
res$w2distance.trt.sigma <- NA
res$bias.ctrl.sigma <- NA
res$sd.ctrl.sigma <- NA
res$w2distance.ctrl.sigma <- NA
res$bias.trt.median.pop <- NA  # Not applicable for PSCL
res$w2distance.trt.median.pop <- NA
res$bias.ctrl.median.pop <- NA  # Not applicable for PSCL
res$w2distance.ctrl.median.pop <- NA
res$iteration <- rep(seq_len(n_replicates), times = nrow(res) / n_replicates)

# For threshold calibration under H0
if (hypo == "null") {
  all_decisions <- numeric(n_replicates * length(cor))
  decision_idx <- 1
}

cat("Loaded threshold from H0 run:", threshold)

for (c_idx in seq_along(cor)){
  c <- cor[c_idx]
for (replicate_id in seq_len(n_replicates)){

  set.seed(seed[ee])
  ee <- ee + 1

  # Calculate row index directly
  row_idx <- (c_idx - 1) * n_replicates + replicate_id

  # Use tryCatch to compute results, return list on success or NULL on error
  iter_result <- tryCatch({

    # LRC-BART MODIFICATION START
    # Each generated file contains the current replicate's RCT and RWD rows.
    filename_rct <- file.path(projectDir, data_folder,
                              paste0(data_prefix, replicate_id, ".RData"))
    filename_rwd <- filename_rct
    # LRC-BART MODIFICATION END
    
    # Read RCT data from current replicate
    data_rct_full <- readRDS(filename_rct)
    rct_idx <- data_rct_full$X[,"D"] == 1

    # Read RWD data from current replicate
    data_rwd_full <- readRDS(filename_rwd)
    rwd_idx <- data_rwd_full$X[,"D"] == 0

    # Combine RCT and RWD from current replicate
    data <- cbind(
      rbind(data_rct_full$X[rct_idx, ],
            data_rwd_full$X[rwd_idx, ]),
      data.frame("y" = c(data_rct_full$y[rct_idx],
                         data_rwd_full$y[rwd_idx]))
    )

    # Also need filename variable for reading true values later
    filename <- filename_rct

    ps <- psrwe_est(data, v_covs = paste("X", 1:10, sep = ""),
                    v_grp = "D", cur_grp_level = "1",
                    v_arm = "Z", ctl_arm_level = "0")

    borrow <- psrwe_borrow(ps, total_borrow = 100,
                           method = "distance")
    fit <- psrwe_compl(borrow, v_outcome = "y",
                       outcome_type = "continuous")
    fit_ci <- psrwe_outana(fit)

    #-------------------------------------
    #---------- Calculate ATE ------------
    #-------------------------------------
    # JOINT LRC-BART MODIFICATION START
    eff <- readRDS(filename)$treat_eff_true
    # JOINT LRC-BART MODIFICATION END
    eff_star <- readRDS(filename)$treat_eff_star

    # Posterior estimate (point estimate for PSCL)
    delta_hat <- fit$Effect$Overall_Estimate$Mean
    print(delta_hat)
    
    # Within replicate metrics
    bias <- delta_hat - eff
    w2distance <- NA

    # Across replicate metrics
    rmse <- (delta_hat - eff)^2

    # Other
    SD <- fit$Effect$Overall_Estimate$StdErr
    w1distance <- NA
    ci_lower <- fit_ci$CI$Effect$Overall_Estimate$Lower
    ci_upper <- fit_ci$CI$Effect$Overall_Estimate$Upper
    ci <- ci_upper - ci_lower
    coverage <- ci_lower <= eff & ci_upper >= eff

    if (hypo == "null") {
      fp <- as.numeric(ci_lower > eff_star)
      tp <- NA
      tp_calibrated <- NA
    } else {
      tp <- as.numeric(ci_lower > eff_star)
      fp <- NA
      tp_calibrated <- NA
    }

    # Return results as a list
    list(bias = bias, sd = SD, rmse = rmse, w1distance = w1distance,
         w2distance = w2distance, ci = ci, coverage = coverage,
         tp_calibrated = tp_calibrated, tp = tp, fp = fp)

  }, error = function(e) {
    cat(paste0("\nError in iteration ", replicate_id, " for cor = ", c, ":\n"))
    cat(paste0("  ", e$message, "\n"))
    cat("  Skipping this iteration and setting results to NA\n\n")
    NULL  # Return NULL on error
  })

  # Assign results outside tryCatch (if successful)
  if (!is.null(iter_result)) {
    res[row_idx, "bias"] <- iter_result$bias
    res[row_idx, "sd"] <- iter_result$sd
    res[row_idx, "rmse"] <- iter_result$rmse
    res[row_idx, "w1distance"] <- iter_result$w1distance
    res[row_idx, "w2distance"] <- iter_result$w2distance
    res[row_idx, "ci"] <- iter_result$ci
    res[row_idx, "coverage"] <- iter_result$coverage
    res[row_idx, "tp_calibrated"] <- iter_result$tp_calibrated
    res[row_idx, "tp"] <- iter_result$tp
    res[row_idx, "fp"] <- iter_result$fp
  }
  # If iter_result is NULL (error), res stays NA (default)

  print(paste0("Done for replicate ", replicate_id, " cor = ", c))
}
  print(paste0("Done for cor = ", c))
}

# Note: PSCL method doesn't support threshold calibration as it doesn't provide posterior samples
# No threshold calibration section needed for PSCL

# Note: RMSE is already set to NA for each replicate above since PSCL
# doesn't provide posterior samples needed to calculate it properly

# LRC-BART MODIFICATION START
saveRDS(res, file = file.path(projectDir, "res",
  paste0("PSCL_p", p_obs, "_", scenario_id, "_", hypo, ".RData")))
# LRC-BART MODIFICATION END
