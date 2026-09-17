rm(list=ls())
# LRC-BART MODIFICATION START
mainDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing"
projectDir <- file.path(mainDir, "lrcbart-sim-gaussian")
n_replicates <- 100L
# LRC-BART MODIFICATION END
library(RBesT)

# LRC-BART MODIFICATION START
data_folder <- "data"
# LRC-BART MODIFICATION END

sc <- 2
# LRC-BART MODIFICATION START
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(projectDir, data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(projectDir, data_folder))
# LRC-BART MODIFICATION END
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# LRC-BART MODIFICATION START
# The inherited 100/75/50 MAP targets remain here. Their s^2 values are now
# calibrated and checkpointed for each generated replicate.
target_Ns <- c(100L, 75L, 50L)
# LRC-BART MODIFICATION END

# LRC-BART ADDITION START
delta_rwd <- 0
region <- "none"
# The scenario filename tag is assembled inline after `cor` is set.
# LRC-BART ADDITION END

threshold <- 0.95

# LRC-BART MODIFICATION START
# Strength of unmeasured confounding
if (sc %in% c(1, 2, 4, 5)){
  cor <- 1
}
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[1]
}

# Validate the configured replicate count against this scenario
scenario_id <- paste0("sc", sc)
if (sc == 3) scenario_id <- paste0(scenario_id, "_cor", cor)
if (sc == 4) scenario_id <- paste0(scenario_id, "_d", delta_rwd)
if (sc == 5)
  scenario_id <- paste0(scenario_id, "_", region, "_d", delta_rwd)
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
    paste0("MAP_p", p_obs, "_", scenario_id, "_threshold.RData"))

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}
# LRC-BART MODIFICATION END

# LRC-BART MODIFICATION START
for (.idx in seq_along(target_Ns)) {
  target_N <- target_Ns[.idx]

  cat("\n========== Running with N_target =", target_N, "==========\n")
# LRC-BART MODIFICATION END

  set.seed(6)
  seed <- sample(1:10000, n_replicates*4, replace = F)
  ee <- 1

  res <- data.frame(
    dist = "HalfNormal",
    c = cor
  )
  res <- res[rep(1:nrow(res), each = n_replicates), ]
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
  res$pehe_subj <- NA  # Not applicable for MAP
  res$bias_subj <- NA  # Not applicable for MAP
  res$bias.trt.sigma <- NA  # MAP method doesn't estimate these directly
  res$sd.trt.sigma <- NA
  res$w2distance.trt.sigma <- NA
  res$bias.ctrl.sigma <- NA
  res$sd.ctrl.sigma <- NA
  res$w2distance.ctrl.sigma <- NA
  res$bias.trt.median.pop <- NA  # Not applicable for MAP
  res$w2distance.trt.median.pop <- NA
  res$bias.ctrl.median.pop <- NA  # Not applicable for MAP
  res$w2distance.ctrl.median.pop <- NA
  res$iteration <- rep(seq_len(n_replicates), times = nrow(res) / n_replicates)

  # For threshold calibration under H0
  if (hypo == "null") {
    all_decisions <- numeric(n_replicates * length(cor))
    decision_idx <- 1
  }

  cat("Loaded threshold from H0 run:", threshold)

  for (c in cor){
    for (replicate_id in seq_len(n_replicates)){
      
      # LRC-BART MODIFICATION START
      this_seed <- seed[ee]
      set.seed(this_seed)
      # LRC-BART MODIFICATION END
      ee <- ee + 1
      
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
      data <- data.frame(
        rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
              data_rwd_full$X[rwd_idx, paste0("X",1:10)])
      )
      data$D <- c(data_rct_full$X[rct_idx, "D"],
                  data_rwd_full$X[rwd_idx, "D"])
      data$Z <- c(data_rct_full$X[rct_idx, "Z"],
                  data_rwd_full$X[rwd_idx, "Z"])
      data$Y <- c(data_rct_full$y[rct_idx],
                  data_rwd_full$y[rwd_idx])
      
      # Also need filename variable for reading true values later
      filename <- filename_rct
      
      print(dim(data))

      # LRC-BART ADDITION START
      # Calibrate once for this replicate; the checkpoint is reused by the
      # other two target passes through this inherited outer target loop.
      map_data_file <- filename_rwd
      map_target_Ns <- target_Ns
      map_seed <- this_seed
      map_checkpoint_dir <- file.path(projectDir, "res", "ess")
      source(file.path(projectDir, "ess_local", "ess_cal_map.R"),
             local = TRUE)
      target_row <- map_calibration$targets[
        map_calibration$targets$target_N == target_N, , drop = FALSE]
      if (nrow(target_row) != 1L)
        stop("Missing MAP ESS calibration target ", target_N)
      prior <- target_row$s0_sq
      # LRC-BART ADDITION END
      
      # Control
      historical_mean <- mean(data[data$D == 0 & data$Z == 0, "Y"])  
      historical_sd   <- sd(data[data$D == 0 & data$Z == 0, "Y"])  
      historical_n    <- nrow(data[data$D == 0 & data$Z == 0, ])
      
      historical_data <- data.frame(
        y = historical_mean, 
        n = historical_n,
        y.se = historical_sd / sqrt(historical_n), 
        study = factor(1:length(historical_mean)) 
      )
      
      map_prior <- gMAP(
        formula = cbind(y, y.se) ~ 1 | study,  
        data = historical_data,
        beta.prior = 10,  
        # tau.dist = "HalfNormal",
        # tau.prior = c(0, 1),
        tau.dist = "InvGamma",
        tau.prior = c(3/2, 3*prior/2),
        family = gaussian
      )
      
      # samples <- unlist(map_prior$fit@sim$samples[[1]]["theta_pred"])
      map_mixture <- automixfit(map_prior)
      print(ess(map_mixture, sigma = 1.5))
      
      new_mean <- mean(data[data$D == 1 & data$Z == 0, "Y"])  
      new_sd   <- sd(data[data$D == 1 & data$Z == 0, "Y"])  
      new_n    <- nrow(data[data$D == 1 & data$Z == 0, ])
      
      posterior_ctrl <- postmix(map_mixture, 
                                m = new_mean, 
                                se = new_sd / sqrt(new_n))
      
      # Treatment
      treatment_mean <- mean(data[data$D == 1 & data$Z == 1, "Y"])  
      treatment_sd   <- sd(data[data$D == 1 & data$Z == 1, "Y"])  
      treatment_n    <- nrow(data[data$D == 1 & data$Z == 1, ])
      
      # Define a weakly informative normal prior for treatment
      weak_prior <- mixnorm(c(1,              #w
                              0,              #mean
                              10              #se
      ), param = 'ms')
      posterior_trt <- postmix(weak_prior, 
                               m = treatment_mean, 
                               se = treatment_sd / sqrt(treatment_n))
      
      #-------------------------------------
      #---------- Calculate ATE ------------
      #-------------------------------------
      samples_trt <- rmix(posterior_trt, 1500)
      samples_ctrl <- rmix(posterior_ctrl, 1500)
      
      post_samples <- samples_trt - samples_ctrl
      eff <- readRDS(filename)$treat_eff
      eff_star <- readRDS(filename)$treat_eff_star
      
      # Posterior estimate
      delta_hat <- mean(post_samples, na.rm = TRUE)
      
      # Within replicate metrics
      bias <- delta_hat - eff
      w2distance <- sqrt(mean((post_samples - eff)^2, na.rm = TRUE))
      
      # Across replicate metrics
      rmse <- (delta_hat - eff)^2
      
      # Other
      SD <- sqrt(var(post_samples, na.rm = TRUE))
      w1distance <- mean(abs(post_samples - eff), na.rm = TRUE)
      coverage <- quantile(post_samples, probs = c(0.025), na.rm = TRUE) <= eff &
        quantile(post_samples, probs = c(0.975), na.rm = TRUE) >= eff
      ci <- diff(quantile(post_samples, c(0.025, 0.975), na.rm = TRUE))
      
      # True Positive/Power
      decision <- mean(post_samples > eff_star, na.rm = TRUE)
      
      tp_calibrated <- as.numeric(decision > threshold)
      tp <- as.numeric(decision > 0.95)
      
      if (hypo == "null") {
        all_decisions[decision_idx] <- decision
        decision_idx <- decision_idx + 1
      }
      
      # False Positive
      if (hypo == "null") {
        fp <- tp
      } else {
        fp <- NA
      }
      
      res[which(res$c==c & res$iteration==replicate_id), "bias"] <- bias
      res[which(res$c==c & res$iteration==replicate_id), "sd"] <- SD
      res[which(res$c==c & res$iteration==replicate_id), "rmse"] <- rmse
      res[which(res$c==c & res$iteration==replicate_id), "w1distance"] <- w1distance
      res[which(res$c==c & res$iteration==replicate_id), "w2distance"] <- w2distance
      res[which(res$c==c & res$iteration==replicate_id), "ci"] <- ci
      res[which(res$c==c & res$iteration==replicate_id), "coverage"] <- coverage
      res[which(res$c==c & res$iteration==replicate_id), "tp_calibrated"] <- tp_calibrated
      res[which(res$c==c & res$iteration==replicate_id), "tp"] <- tp
      res[which(res$c==c & res$iteration==replicate_id), "fp"] <- fp
      
      # Subject-wise calculations - not applicable for MAP (operates at group level)
      res[which(res$c==c & res$iteration==replicate_id), "pehe_subj"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "bias_subj"] <- NA
      
      # Variance estimation for treatment and control groups
      # MAP method doesn't model residual variance, so these remain NA
      res[which(res$c==c & res$iteration==replicate_id), "bias.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "sd.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "w2distance.trt.sigma"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "bias.ctrl.sigma"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "sd.ctrl.sigma"] <- NA
      res[which(res$c==c & res$iteration==replicate_id), "w2distance.ctrl.sigma"] <- NA
      
      #--------------------------------------------------
      #------ Population mean estimation ----------------
      #--------------------------------------------------
      true_mean_trt_pop <- readRDS(filename)$true_mean_trt_pop
      true_mean_ctrl_pop <- readRDS(filename)$true_mean_ctrl_pop
      
      # Calculate population means
      mean_trt_pop_hat <- mean(samples_trt)
      mean_ctrl_pop_hat <- mean(samples_ctrl)
      
      # Calculate metrics
      bias_trt_median_pop <- mean_trt_pop_hat - true_mean_trt_pop
      w2distance_trt_median_pop <- sqrt(mean((samples_trt - true_mean_trt_pop)^2))
      bias_ctrl_median_pop <- mean_ctrl_pop_hat - true_mean_ctrl_pop
      w2distance_ctrl_median_pop <- sqrt(mean((samples_ctrl - true_mean_ctrl_pop)^2))
      
      res[which(res$c==c & res$iteration==replicate_id), "bias.trt.median.pop"] <- bias_trt_median_pop
      res[which(res$c==c & res$iteration==replicate_id), "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
      res[which(res$c==c & res$iteration==replicate_id), "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
      res[which(res$c==c & res$iteration==replicate_id), "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop
      
      print(paste0("Done for replicate ", replicate_id, " cor = ", c))
    }
    print(paste0("Done for cor = ", c))
  }

  # Threshold calibration under null hypothesis
  if (hypo == "null") {
    cat("\n=== Threshold Calibration (H0) ===\n")
    cat("Original threshold: ", threshold, "\n")

    # Grid search for exact type I error = 0.05 (more precise)
    candidate_thresholds <- seq(0.1, 0.999, by = 0.001)
    type1_errors <- sapply(candidate_thresholds, function(thresh) {
      mean(all_decisions > thresh, na.rm = TRUE)
    })

    best_idx <- which.min(abs(type1_errors - 0.05))
    calibrated_threshold_grid <- candidate_thresholds[best_idx]
    actual_type1_error <- type1_errors[best_idx]

    cat("Selected threshold: ", calibrated_threshold_grid, "\n")
    cat("Actual Type I error: ", actual_type1_error, "\n")

    # LRC-BART MODIFICATION START
    # Save calibrated threshold for use in alternative hypothesis run.
    threshold_file <- file.path(projectDir, "res",
      paste0("MAP_p", p_obs, "_", scenario_id, "_threshold.RData"))
    saveRDS(calibrated_threshold_grid, file = threshold_file)
    # LRC-BART MODIFICATION END
  }

  # LRC-BART MODIFICATION START
  # Save results with the inherited target label in the filename.
  saveRDS(res, file = file.path(projectDir, "res",
    paste0("MAP_p", p_obs, "_", scenario_id, "_N", target_N,
           "_", hypo, ".RData")))
  # LRC-BART MODIFICATION END

} # End of target_N loop
