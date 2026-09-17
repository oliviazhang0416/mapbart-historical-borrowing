rm(list = ls())
# LRC-BART MODIFICATION START
mainDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing"
projectDir <- file.path(mainDir, "lrcbart-sim-gaussian")
n_replicates <- 100L
# LRC-BART MODIFICATION END
.scpp_cache <- file.path(tempdir(), "scpp_BARTv3")
dir.create(.scpp_cache, showWarnings = FALSE, recursive = TRUE)
library(Rcpp)
library(RcppEigen)
library(gtools)

# LRC-BART MODIFICATION START
data_folder <- "data"
# LRC-BART MODIFICATION END

sc <- 3
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
if (sc %in% c(1, 2, 4, 5)){
  cor <- 1
}
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[3]
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
    paste0("BARTv3_p", p_obs, "_", scenario_id, "_threshold.RData"))

  if (file.exists(threshold_file)) {
    threshold <- readRDS(threshold_file)
  }
}
# LRC-BART MODIFICATION END

ndpost=1500L
nskip=100L
keepevery=1L

# alpha_beta <- data.frame(alpha = c(0.9, 0.95, 0.98),  # -> bigger trees
#                          beta = c(2.5, 2, 1.5))       # -> bigger trees
# ntree <- c(1, 5, 10, 50, 100, 200)
alpha_beta <- data.frame(alpha = c(0.95),  # -> bigger trees
                         beta = c(2))       # -> bigger trees
ntree <- c(50)

# LRC-BART MODIFICATION START
sourceCpp("/Users/oliviazhang/Desktop/wBART/cwbart.cpp", cacheDir = .scpp_cache)
# LRC-BART MODIFICATION END

set.seed(6)
seed <- sample(1:10000, n_replicates*4, replace = F)
ee <- 1

res <- expand.grid(
  id = 1:nrow(alpha_beta),
  H = ntree,
  c = cor
)

# Merge with alpha-beta combinations
res <- merge(res,
             alpha_beta,
             by.x = "id",
             by.y = "row.names",
             all.x = TRUE)
res <- res[ , c("alpha", "beta", "H", "c")]
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
res$pehe_subj <- NA
res$bias_subj <- NA
res$bias.trt.sigma <- NA
res$sd.trt.sigma <- NA
res$w2distance.trt.sigma <- NA
res$bias.ctrl.sigma <- NA
res$sd.ctrl.sigma <- NA
res$w2distance.ctrl.sigma <- NA
res$bias.trt.median.pop <- NA
res$w2distance.trt.median.pop <- NA
res$bias.ctrl.median.pop <- NA
res$w2distance.ctrl.median.pop <- NA
res$iteration <- rep(seq_len(n_replicates), times = nrow(res) / n_replicates)

# For threshold calibration under H0
if (hypo == "null") {
  all_decisions <- numeric(nrow(alpha_beta) * length(cor) * length(ntree) * n_replicates)
  decision_idx <- 1
}

cat("Loaded threshold from H0 run:", threshold)


for (which in 1:nrow(alpha_beta)){
  alpha <- alpha_beta[which, "alpha"]
  beta <- alpha_beta[which, "beta"]
for (c in cor){
for (H in ntree){
  for (replicate_id in seq_len(n_replicates)){
    this_seed <- seed[ee]   # one seed per iteration; both arms share it
    set.seed(this_seed)
    ee <- ee + 1


    data = list()

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
    data$X <- rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
                    data_rwd_full$X[rwd_idx, paste0("X",1:10)])
    data$D <- c(data_rct_full$X[rct_idx, "D"],
                data_rwd_full$X[rwd_idx, "D"])
    data$Z <- c(data_rct_full$X[rct_idx, "Z"],
                data_rwd_full$X[rwd_idx, "Z"])
    data$Y <- c(data_rct_full$y[rct_idx],
                data_rwd_full$y[rwd_idx])

    # Also need filename variable for reading true values later
    filename <- filename_rct

    #------------------------------------------------
    #----------- RCT + Historical Control -----------
    #------------------------------------------------
    # v3: include data source indicator D in the control arm model
    x.train = data.frame(X = as.matrix(data$X[data$Z == 0, ]),
                         D = data$D[data$Z == 0])
    y.train = data$Y[data$Z == 0]
    x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]),
                        D = 1)

    sparse=FALSE
    theta=0
    omega=1
    a=0.5
    b=1
    augment=FALSE
    rho=NULL
    xinfo=matrix(0.0,0,0)
    usequants=FALSE
    cont=FALSE
    rm.const=TRUE
    sigest=NA
    sigdf=3
    sigquant=.90
    k=2.0
    power=beta
    base=alpha
    sigmaf=NA
    lambda=NA
    fmean=0
    w=rep(1,length(y.train))

    numcut=100L
    nkeeptrain=ndpost
    nkeeptest=ndpost
    nkeeptestmean=ndpost
    nkeeptreedraws=ndpost
    printevery=10000L
    transposed=FALSE

    n = length(y.train)

    if(!transposed) {
      # LRC-BART MODIFICATION START
      source("/Users/oliviazhang/Desktop/bartModelMatrix.R")
      # LRC-BART MODIFICATION END
      temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                             cont=cont, xinfo=xinfo, rm.const=rm.const)
      x.train = t(temp$X)
      numcut = temp$numcut
      xinfo = temp$xinfo
      if(length(x.test)>0) {
        x.test = bartModelMatrix(x.test)
        x.test = t(x.test[ , temp$rm.const])
      }
      rm.const <- temp$rm.const
      grp <- temp$grp
      rm(temp)
    }

    p = nrow(x.train)
    np = ncol(x.test)
    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    y.train = y.train-fmean
    if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
      nkeeptrain=ndpost
    }
    if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
      nkeeptest=ndpost
    }
    if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
      nkeeptestmean=ndpost
    }
    if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
      nkeeptreedraws=ndpost
    }

    nu=sigdf
    if(is.na(lambda)) {
      if(is.na(sigest)) {
        if(p < n) {
          df = data.frame(t(x.train),y.train)
          lmf = lm(y.train~.,df)
          sigest = summary(lmf)$sigma
        } else {
          sigest = sd(y.train)
        }
      }
      qchi = qchisq(1.0-sigquant,nu)
      lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
    } else {
      sigest=sqrt(c)
    }

    if(is.na(sigmaf)) {
      tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
    } else {
      tau = sigmaf/sqrt(ntree)
    }

    set.seed(this_seed)
    res_ctrl = cwbart(n,  #number of observations in training data
                      p,  #dimension of x
                      np, #number of observations in test data
                      x.train,   #pxn training data x
                      y.train,   #pxn training data x
                      x.test,   #p*np test data x
                      ntree,
                      numcut,
                      ndpost*keepevery,
                      nskip,
                      power,
                      base,
                      tau,
                      nu,
                      lambda,
                      sigest,
                      w,
                      sparse,
                      theta,
                      omega,
                      grp,
                      a,
                      b,
                      rho,
                      augment,
                      nkeeptrain,
                      nkeeptest,
                      nkeeptestmean,
                      nkeeptreedraws,
                      printevery,
                      xinfo
    )

    if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
    if(keepevery>1) res_ctrl$sigma = c(res_ctrl$sigma[nskip.],
                                       res_ctrl$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
    res_ctrl$sigma = res_ctrl$sigma[-(nskip.)]

    # #------------------------------------------------
    # #---------------- RCT Treatment -----------------
    # #------------------------------------------------
    x.train = data.frame(X = as.matrix(data$X[data$Z == 1, ]))
    y.train = data$Y[data$Z == 1]
    x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))

    sparse=FALSE
    theta=0
    omega=1
    a=0.5
    b=1
    augment=FALSE
    rho=NULL
    xinfo=matrix(0.0,0,0)
    usequants=FALSE
    cont=FALSE
    rm.const=TRUE
    sigest=NA
    sigdf=3
    sigquant=.90
    # k=2.0
    power=beta
    base=alpha
    k=3.0
    sigmaf=NA
    lambda=NA
    fmean=0
    w=rep(1,length(y.train))
    numcut=100L

    nkeeptrain=ndpost
    nkeeptest=ndpost
    nkeeptestmean=ndpost
    nkeeptreedraws=ndpost
    printevery=10000L
    transposed=FALSE

    n = length(y.train)

    if(!transposed) {
      # LRC-BART MODIFICATION START
      source("/Users/oliviazhang/Desktop/bartModelMatrix.R")
      # LRC-BART MODIFICATION END
      temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                             cont=cont, xinfo=xinfo, rm.const=rm.const)
      x.train = t(temp$X)
      numcut = temp$numcut
      xinfo = temp$xinfo
      if(length(x.test)>0) {
        x.test = bartModelMatrix(x.test)
        x.test = t(x.test[ , temp$rm.const])
      }
      rm.const <- temp$rm.const
      grp <- temp$grp
      rm(temp)
    }

    p = nrow(x.train)
    np = ncol(x.test)
    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    y.train = y.train-fmean
    if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
      nkeeptrain=ndpost
    }
    if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
      nkeeptest=ndpost
    }
    if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
      nkeeptestmean=ndpost
    }
    if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
      nkeeptreedraws=ndpost
    }

    nu=sigdf
    if(is.na(lambda)) {
      if(is.na(sigest)) {
        if(p < n) {
          df = data.frame(t(x.train),y.train)
          lmf = lm(y.train~.,df)
          sigest = summary(lmf)$sigma
        } else {
          sigest = sd(y.train)
        }
      }
      qchi = qchisq(1.0-sigquant,nu)
      lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
    } else {
      sigest=sqrt(lambda)
    }

    if(is.na(sigmaf)) {
      tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
    } else {
      tau = sigmaf/sqrt(ntree)
    }


    set.seed(this_seed)
    res_trt = cwbart(n,  #number of observations in training data
                     p,  #dimension of x
                     np, #number of observations in test data
                     x.train,   #pxn training data x
                     y.train,   #pxn training data x
                     x.test,   #p*np test data x
                     ntree,
                     numcut,
                     ndpost*keepevery,
                     nskip,
                     power,
                     base,
                     tau,
                     nu,
                     lambda,
                     sigest,
                     w,
                     sparse,
                     theta,
                     omega,
                     grp,
                     a,
                     b,
                     rho,
                     augment,
                     nkeeptrain,
                     nkeeptest,
                     nkeeptestmean,
                     nkeeptreedraws,
                     printevery,
                     xinfo
    )

    if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
    if(keepevery>1) res_trt$sigma = c(res_trt$sigma[nskip.],
                                      res_trt$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
    res_trt$sigma = res_trt$sigma[-(nskip.)]

    #-------------------------------------
    #---------- Calculate ATE ------------
    #-------------------------------------
    post_samples <- rowMeans(res_trt$yhat.test - res_ctrl$yhat.test)
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

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "bias"] <- bias

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "sd"] <- SD

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "rmse"] <- rmse

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w1distance"] <- w1distance

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w2distance"] <- w2distance

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "ci"] <- ci

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "coverage"] <- coverage

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "tp_calibrated"] <- tp_calibrated

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "tp"] <- tp

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "fp"] <- fp

    #------------------------------------------
    #------ Calculate subject-wise ------------
    #------------------------------------------
    post_samples_i <- res_trt$yhat.test - res_ctrl$yhat.test   # n_draws x n_subjects
    eff_i <- readRDS(filename)$eff_i

    # Posterior estimates per subject
    delta_hat_i <- apply(post_samples_i, 2, mean)    # n_subjects

    # Within replicate metrics
    bias_subj <- mean(delta_hat_i - eff_i)
    pehe_subj <- sqrt(mean((delta_hat_i - eff_i)^2))

    res[which(res$c==c & res$iteration==replicate_id), "bias_subj"] <- bias_subj
    res[which(res$c==c & res$iteration==replicate_id), "pehe_subj"] <- pehe_subj

    #-------------------------------------
    #------ Variance estimation ----------
    #-------------------------------------
    true_sigma_trt <- true_sigma_ctrl <- readRDS(filename)$sigma_rct

    # Treatment group
    sigma_samples_trt <- res_trt$sigma
    sigma_est_trt <- mean(sigma_samples_trt)
    sd_trt_sigma <- sd(sigma_samples_trt)
    w2distance_trt_sigma <- sqrt(mean((sigma_samples_trt - true_sigma_trt)^2, na.rm = TRUE))

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "bias.trt.sigma"] <- sigma_est_trt - true_sigma_trt
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "sd.trt.sigma"] <- sd_trt_sigma
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w2distance.trt.sigma"] <- w2distance_trt_sigma

    # Control group
    sigma_samples_ctrl <- res_ctrl$sigma
    sigma_est_ctrl <- mean(sigma_samples_ctrl)
    sd_ctrl_sigma <- sd(sigma_samples_ctrl)
    w2distance_ctrl_sigma <- sqrt(mean((sigma_samples_ctrl - true_sigma_trt)^2, na.rm = TRUE))

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "bias.ctrl.sigma"] <- sigma_est_ctrl - true_sigma_ctrl
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "sd.ctrl.sigma"] <- sd_ctrl_sigma
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w2distance.ctrl.sigma"] <- w2distance_ctrl_sigma

    #--------------------------------------------------
    #------ Population mean estimation ----------------
    #--------------------------------------------------
    true_mean_trt_pop <- readRDS(filename)$true_mean_trt_pop
    true_mean_ctrl_pop <- readRDS(filename)$true_mean_ctrl_pop

    # Calculate population means
    mean_trt_pop_hat <- mean(rowMeans(res_trt$yhat.test))
    mean_ctrl_pop_hat <- mean(rowMeans(res_ctrl$yhat.test))

    # Calculate metrics
    mu_pred_trt_samples <- rowMeans(res_trt$yhat.test)
    mu_pred_ctrl_samples <- rowMeans(res_ctrl$yhat.test)
    bias_trt_median_pop <- mean_trt_pop_hat - true_mean_trt_pop
    w2distance_trt_median_pop <- sqrt(mean((mu_pred_trt_samples - true_mean_trt_pop)^2))
    bias_ctrl_median_pop <- mean_ctrl_pop_hat - true_mean_ctrl_pop
    w2distance_ctrl_median_pop <- sqrt(mean((mu_pred_ctrl_samples - true_mean_ctrl_pop)^2))

    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "bias.trt.median.pop"] <- bias_trt_median_pop
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w2distance.trt.median.pop"] <- w2distance_trt_median_pop
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "bias.ctrl.median.pop"] <- bias_ctrl_median_pop
    res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
        "w2distance.ctrl.median.pop"] <- w2distance_ctrl_median_pop

    print(paste0("Done for replicate ", replicate_id, " cor = ", c))
  }
  print(paste0("Done for cor = ", c,
               " H = ", H,
               " alpha = ", alpha,
               " beta = ", beta))
}
}
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
    paste0("BARTv3_p", p_obs, "_", scenario_id, "_threshold.RData"))
  saveRDS(calibrated_threshold_grid, file = threshold_file)
  # LRC-BART MODIFICATION END
}

# LRC-BART MODIFICATION START
saveRDS(res, file = file.path(projectDir, "res",
  paste0("BARTv3_p", p_obs, "_", scenario_id, "_", hypo, ".RData")))
# LRC-BART MODIFICATION END
