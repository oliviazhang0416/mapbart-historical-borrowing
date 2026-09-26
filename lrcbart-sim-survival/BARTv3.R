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
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-survival")
n_replicates <- 100L
# JOINT LRC-BART MODIFICATION START
dir.create(file.path(projectDir, "res"), recursive = TRUE, showWarnings = FALSE)
# JOINT LRC-BART MODIFICATION END
# LRC-BART MODIFICATION END
.scpp_cache <- file.path(tempdir(), "scpp_BARTv3")
dir.create(.scpp_cache, showWarnings = FALSE, recursive = TRUE)
library(Rcpp)
library(RcppEigen)
library(gtools)
library(survival)
library(dplyr)

# LRC-BART MODIFICATION START
source(file.path(projectDir, "rmst_helpers.R"))
# LRC-BART MODIFICATION END
rmst_tau <- 3   # RMST restriction horizon (admin-censoring horizon); additional output
rmst_control_sigma <- "own"  # use the control model's own residual SD for control-arm RMST

# LRC-BART MODIFICATION START
data_folder <- "data"
# LRC-BART MODIFICATION END

# JOINT LRC-BART MODIFICATION START
sc <- 1
variant <- ""
# JOINT LRC-BART MODIFICATION END
# Infer p_obs from any data file in the folder (count columns named X1, X2, ...)
sample_files <- list.files(file.path(projectDir, data_folder),
                           pattern = "^data_.*\\.RData$", full.names = TRUE)
if (length(sample_files) == 0) stop("No data files found in ", file.path(projectDir, data_folder))
p_obs <- sum(grepl("^X\\d+$", colnames(readRDS(sample_files[1])$X)))
hypo <- "alternative"  # "null" or "alternative"

# LRC-BART ADDITION START
delta_rwd <- 0
region <- "none"
# LRC-BART ADDITION END

threshold <- 0.95

# Strength of unmeasured confounding
# LRC-BART MODIFICATION START
# JOINT LRC-BART MODIFICATION START
if (sc %in% c(1, 2)){
# JOINT LRC-BART MODIFICATION END
  cor <- 1
}
# LRC-BART MODIFICATION END
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)[2]
}

# Validate the configured replicate count against this scenario
# LRC-BART MODIFICATION START
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
# Auto-load calibrated threshold when running under alternative.
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
sourceCpp(file.path(.lrcRoot, "aBART", "cabart.cpp"), cacheDir = .scpp_cache)
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
res$rmst_tau <- NA
res$rmst_true <- NA
res$rmst_hat <- NA
res$bias_rmst <- NA
res$sd_rmst <- NA
res$ci_rmst <- NA
res$coverage_rmst <- NA
res$decision_rmst <- NA
res$tp_rmst <- NA
res$rmse_rmst <- NA
res$w1distance_rmst <- NA
res$w2distance_rmst <- NA
res$tp_calibrated_rmst <- NA
res$bias.trt.rmst.pop <- NA
res$w2distance.trt.rmst.pop <- NA
res$bias.ctrl.rmst.pop <- NA
res$w2distance.ctrl.rmst.pop <- NA
res$bias_subj_rmst <- NA
res$pehe_subj_rmst <- NA
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
        
        # Read RWD data from iteration 1
        data_rwd_full <- readRDS(filename_rwd)
        rwd_idx <- data_rwd_full$X[,"D"] == 0
        
        # Combine RCT (from current replicate_id) and RWD (from replicate_id 1)
        data$X <- rbind(data_rct_full$X[rct_idx, paste0("X",1:10)],
                        data_rwd_full$X[rwd_idx, paste0("X",1:10)])
        data$D <- c(data_rct_full$X[rct_idx, "D"],
                    data_rwd_full$X[rwd_idx, "D"])
        data$Z <- c(data_rct_full$X[rct_idx, "Z"],
                    data_rwd_full$X[rwd_idx, "Z"])
        data$Y <- c(data_rct_full$y[rct_idx],
                    data_rwd_full$y[rwd_idx])
        data$event <- c(data_rct_full$delta[rct_idx],
                        data_rwd_full$delta[rwd_idx])
        
        #------------------------------------------------
        #----------- RCT + Historical Control -----------
        #------------------------------------------------
        # v3: include data source indicator D in the control arm model
        x.train = data.frame(X = as.matrix(data$X[data$Z == 0, ]),
                             D = data$D[data$Z == 0])
        y.train = log(data$Y[data$Z == 0])
        event = as.integer(data$event[data$Z == 0])

        x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]),
                            D = 1)
        
        ntype=1
        sparse=FALSE
        theta=0
        omega=1
        a=0.5
        b=1
        augment=FALSE
        rho=NULL
        xinfo=matrix(0.0,0,0)
        usequants=FALSE
        rm.const=TRUE
        sigdf=3
        sigquant=.90
        k=2.0
        power=beta
        base=alpha
        sigmaf=NA
        lambda=NA
        offset=0
        w=rep(1,length(y.train))
        numcut=100L
        printevery=10000L
        transposed=FALSE
        
        n = length(y.train)
        
        if(!transposed) {
          # LRC-BART MODIFICATION START
          source(file.path(.lrcRoot, "bartModelMatrix.R"))
          # LRC-BART MODIFICATION END
          temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                                 xinfo=xinfo, rm.const=rm.const)
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
        
        y.train = y.train-offset
        if(is.na(lambda)) {
          # sigest = summary(lm(y.train~.,
          #                     data.frame(t(x.train),y.train)))$sigma
          aft_fit = survreg(Surv(exp(y.train), event) ~ .,
                            data = data.frame(t(x.train)),
                            dist = "lognormal")
          sigest = aft_fit$scale
          
          qchi = qchisq(1.0-sigquant,sigdf)
          lambda = (sigest*sigest*qchi)/sigdf #lambda parameter for sigma prior
        } else {
          sigest=sqrt(lambda)
        }
        
        if(is.na(sigmaf)) {
          tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
          tau = sigmaf/sqrt(ntree)
        }
        
        set.seed(this_seed)
        res_ctrl = cabart(ntype,
                          n,
                          p,
                          np,
                          x.train,
                          y.train,
                          event,
                          x.test,
                          ntree,
                          numcut,
                          ndpost*keepevery,
                          nskip,
                          keepevery,
                          power,
                          base,
                          offset,
                          tau,
                          sigdf,
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
                          printevery,
                          xinfo,
                          seed[ee]  # seed for random number generation
        )
        
        if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
        if(keepevery>1) res_ctrl$sigma = c(res_ctrl$sigma[nskip.],
                                           res_ctrl$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
        res_ctrl$sigma = res_ctrl$sigma[-(nskip.)]
        
        # Calculate median survival times for control
        mu_draws <- res_ctrl$yhat.test
        sig_draw <- res_ctrl$sigma
        w <- rep(1/ncol(mu_draws), ncol(mu_draws))
        
        S_mix_fun <- function(t, mu_vec, sig, w) {
          pnorm(log(t), mean = mu_vec, sd = sig, lower.tail = FALSE) |> as.numeric() |>
            {\(v) sum(w * v)}()
        }
        
        pop_median_root_ctrl <- rep(NA_real_, nrow(mu_draws))
        upper_default <- pmin(exp(apply(mu_draws, 1, max) + 6 * sig_draw), 1e10)
        lower <- min(data$Y[data$event == 1]) / 10
        
        for (m in 1:nrow(mu_draws)) {
          f <- function(t) S_mix_fun(t, mu_draws[m, ], sig_draw[m], w) - 0.5
          if (f(upper_default[m]) > 0) {
            warning(paste0("Control group: f(upper_default[", m, "]) > 0, setting median to NA"))
            pop_median_root_ctrl[m] <- NA
            next
          }
          lo <- max(1e-10, lower)
          fm_lo <- f(lo); fm_hi <- f(upper_default[m])
          if (fm_lo < 0) {
            warning(paste0("Control group: f(lo) < 0 for m=", m, ", setting median to lower bound"))
            pop_median_root_ctrl[m] <- lo
          } else {
            pop_median_root_ctrl[m] <- uniroot(f, interval = c(lo, upper_default[m]))$root
          }
        }
        
        median_survival_ctrl <- pop_median_root_ctrl
        
        #------------------------------------------------
        #---------------- RCT Treatment -----------------
        #------------------------------------------------
        x.train = data.frame(X = as.matrix(data$X[data$Z == 1, ]))
        y.train = log(data$Y[data$Z == 1])
        event = as.integer(data$event[data$Z == 1])
        x.test = data.frame(X = as.matrix(data$X[data$D == 1, ]))
        
        ntype=1
        sparse=FALSE
        theta=0
        omega=1
        a=0.5
        b=1
        augment=FALSE
        rho=NULL
        xinfo=matrix(0.0,0,0)
        usequants=FALSE
        rm.const=TRUE
        sigest=NA
        sigdf=3
        sigquant=.90
        k=2.0
        power=beta
        base=alpha
        sigmaf=NA
        lambda=NA
        offset=0
        w=rep(1,length(y.train))
        numcut=100L
        printevery=10000L
        transposed=FALSE
        
        n = length(y.train)
        
        if(!transposed) {
          # LRC-BART MODIFICATION START
          source(file.path(.lrcRoot, "bartModelMatrix.R"))
          # LRC-BART MODIFICATION END
          temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                                 xinfo=xinfo, rm.const=rm.const)
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
        
        y.train = y.train-offset
        
        if(is.na(lambda)) {
          # sigest = summary(lm(y.train~.,
          #                     data.frame(t(x.train),y.train)))$sigma
          aft_fit = survreg(Surv(exp(y.train), event) ~ .,
                            data = data.frame(t(x.train)),
                            dist = "lognormal")
          sigest = aft_fit$scale
          
          qchi = qchisq(1.0-sigquant,sigdf)
          lambda = (sigest*sigest*qchi)/sigdf #lambda parameter for sigma prior
        } else {
          sigest=sqrt(lambda)
        }
        
        if(is.na(sigmaf)) {
          tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
          tau = sigmaf/sqrt(ntree)
        }
        
        
        set.seed(this_seed)
        res_trt = cabart(ntype,
                         n,
                         p,
                         np,
                         x.train,
                         y.train,
                         event,
                         x.test,
                         ntree,
                         numcut,
                         ndpost*keepevery,
                         nskip,
                         keepevery,
                         power,
                         base,
                         offset,
                         tau,
                         sigdf,
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
                         printevery,
                         xinfo,
                         seed[ee]  # seed for random number generation
        )
        
        if(nskip>0){nskip.=1:nskip }  else{nskip. = 0 }
        if(keepevery>1) res_trt$sigma = c(res_trt$sigma[nskip.],
                                          res_trt$sigma[nskip+seq(1, ndpost*keepevery, keepevery)])
        res_trt$sigma = res_trt$sigma[-(nskip.)]
        
        # Calculate median survival times for treatment
        mu_draws <- res_trt$yhat.test
        sig_draw <- res_trt$sigma
        w <- rep(1/ncol(mu_draws), ncol(mu_draws))
        
        pop_median_root_trt <- rep(NA_real_, nrow(mu_draws))
        upper_default <- pmin(exp(apply(mu_draws, 1, max) + 6 * sig_draw), 1e10)
        
        for (m in 1:nrow(mu_draws)) {
          f <- function(t) S_mix_fun(t, mu_draws[m, ], sig_draw[m], w) - 0.5
          if (f(upper_default[m]) > 0) {
            warning(paste0("Treatment group: f(upper_default[", m, "]) > 0, setting median to NA"))
            pop_median_root_trt[m] <- NA
            next
          }
          lo <- max(1e-10, lower)
          fm_lo <- f(lo); fm_hi <- f(upper_default[m])
          if (fm_lo < 0) {
            warning(paste0("Treatment group: f(lo) < 0 for m=", m, ", setting median to lower bound"))
            pop_median_root_trt[m] <- lo
          } else {
            pop_median_root_trt[m] <- uniroot(f, interval = c(lo, upper_default[m]))$root
          }
        }
        
        median_survival_trt <- pop_median_root_trt
        
        #-------------------------------------
        #---------- Calculate ATE ------------
        #-------------------------------------
        post_samples <- median_survival_trt / median_survival_ctrl
        data_tmp <- readRDS(filename_rct)
        eff <- exp(data_tmp$treat_eff_true)
        eff_star <- exp(data_tmp$treat_eff_star)

        # ---- RMST(tau)-ratio (additional output; median-ratio left unchanged) ----
        .rm <- compute_rmst_metrics(res_trt$yhat.test, res_trt$sigma,
                                    res_ctrl$yhat.test, tail(res_ctrl$sigma, length(res_trt$sigma)),
                                    data_tmp, tau = rmst_tau, threshold = threshold,
                                    control_sigma = rmst_control_sigma)
        .sel <- which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id)
        for (.nm in names(.rm)) res[.sel, .nm] <- .rm[[.nm]]

        # Posterior estimate
        delta_hat <- median(post_samples, na.rm = TRUE)
        
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
            "coverage"] <- coverage
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
            "ci"] <- ci
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
            "tp_calibrated"] <- tp_calibrated
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
            "tp"] <- tp
        
        res[which(res$alpha==alpha & res$beta==beta & res$c==c & res$H==H & res$iteration==replicate_id),
            "fp"] <- fp
        
        #------------------------------------------
        #------ Calculate subject-wise ------------
        #------------------------------------------
        post_samples_i <- exp(res_trt$yhat.test) / exp(res_ctrl$yhat.test)   # n_draws x n_subjects
        eff_i <- exp(readRDS(filename_rct)$eff_i)

        # Posterior estimates per subject
        delta_hat_i <- apply(post_samples_i, 2, median)    # n_subjects
        
        # post_samples_1i <- exp(res_trt$yhat.test)      # n_draws x n_subjects
        # post_samples_0i <- exp(res_ctrl$yhat.test)     # n_draws x n_subjects 
        # eff_i <- exp(readRDS(filename_rct)$eff_i)
        # 
        # # Posterior estimates per subject
        # delta_hat_1i <- apply(post_samples_1i, 2, median)    # n_subjects
        # delta_hat_0i <- apply(post_samples_0i, 2, median)    # n_subjects
        # delta_hat_i <- delta_hat_1i / delta_hat_0i           # n_subjects
        
        # Within replicate metrics
        bias_subj <- mean(delta_hat_i - eff_i)
        pehe_subj <- sqrt(mean((delta_hat_i - eff_i)^2))
        
        res[which(res$c==c & res$iteration==replicate_id), "bias_subj"] <- bias_subj
        res[which(res$c==c & res$iteration==replicate_id), "pehe_subj"] <- pehe_subj

        #-------------------------------------
        #------ Variance estimation ----------
        #-------------------------------------
        true_sigma_trt <- true_sigma_ctrl <- readRDS(filename_rct)$sigma_rct
        
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
        #------ Population median survival time -----------
        #--------------------------------------------------
        true_median_trt_pop <- readRDS(filename_rct)$true_median_trt_pop
        true_median_ctrl_pop <- readRDS(filename_rct)$true_median_ctrl_pop

        # Use already calculated population medians
        median_trt_pop_hat <- median(median_survival_trt, na.rm = TRUE)
        median_ctrl_pop_hat <- median(median_survival_ctrl, na.rm = TRUE)

        # Calculate metrics
        bias_trt_median_pop <- median_trt_pop_hat - true_median_trt_pop
        w2distance_trt_median_pop <- sqrt(mean((median_survival_trt - true_median_trt_pop)^2, na.rm = TRUE))
        bias_ctrl_median_pop <- median_ctrl_pop_hat - true_median_ctrl_pop
        w2distance_ctrl_median_pop <- sqrt(mean((median_survival_ctrl - true_median_ctrl_pop)^2, na.rm = TRUE))

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
  # Save the null-run threshold under the same scenario identifier.
  threshold_file <- file.path(projectDir, "res",
    paste0("BARTv3_p", p_obs, "_", scenario_id, "_threshold.RData"))
  # LRC-BART MODIFICATION END

  saveRDS(calibrated_threshold_grid, file = threshold_file)
}


# LRC-BART MODIFICATION START
saveRDS(res, file = file.path(projectDir, "res",
  paste0("BARTv3_p", p_obs, "_", scenario_id, "_", hypo, ".RData")))
# LRC-BART MODIFICATION END
