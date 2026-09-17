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
library(MASS)
library(Rlab)

# LRC-BART MODIFICATION START
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-gaussian")
n_replicates <- 100L
# LRC-BART MODIFICATION END

sc <- 3
hypo <- "alternative"

# LRC-BART ADDITION START
# Sc4 is the global RWD shift; Sc5 is regional partial compatibility.
delta_rwd <- 1
region <- "X5"
if (sc <= 3) {
  delta_rwd <- 0
  region <- "none"
}
if (sc == 4) region <- "none"
stopifnot(sc %in% 1:5, region %in% c("none", "X5", "X7"),
          length(delta_rwd) == 1, is.finite(delta_rwd))

# Region indicators and scenario filename tags are assembled inline where
# they are used below.
# LRC-BART ADDITION END

n <- 300

r <- 2
n1 <- (r/(r+1)) * n
n2 <- (1/(r+1)) * n
n3 <- n

# LRC-BART MODIFICATION START
# Retain the original two-pool design while making both exact samples available
# with overwhelming probability under the source-assignment mechanisms.
n_rct_pool <- 4 * (n1 + n2)
n_rwd_pool <- 4 * n3
# LRC-BART MODIFICATION END

p_obs <- 10
conf_cols <- c(5, 6)
extreme_cols <- c(1)

# LRC-BART MODIFICATION START
if (sc %in% c(1, 2, 4, 5)){
  cor <- 1
}
# LRC-BART MODIFICATION END
if (sc == 3){
  cor <- c(-0.5, 0, 0.5)
}

set.seed(6)
seed <- sample(1:10000, n_replicates*length(cor), replace = F)
ee <- 1

seed_rwd <- 2026

# RWD-generation mode: TRUE = one frozen pool (files tagged _fz); FALSE = fresh RWD each replicate
# LRC-BART MODIFICATION START
rwd_frozen <- FALSE
# LRC-BART MODIFICATION END
set.seed(seed_rwd)
seed_rwd_replicate <- sample(1:10000, n_replicates*length(cor), replace = F)   # per-replicate RWD seeds

sd_x <- 1
sd_U <- 1
m_x <- c(1, -1, 0.5, -1, 2, -2, 2, -3, 1, 0)

if (sc == 2) {
  beta_D   <- -1.20
  sel_temp <- 0.5
}
if (sc == 3) {
  beta_D   <- -1.20
  sel_temp <- 0.5
}
if (hypo == "null") {
  eff <- 0
} else {
  eff <- 1.0
}
eff_star <- 0.5

# LRC-BART MODIFICATION START
if (sc %in% c(1, 4, 5)) {
  b0_rct <- 3.0
  b0_rwd <- 3.0
  
  sd_Y_rct <- 1.5
  sd_Y_rwd <- 1.8
  
} else if (sc == 2) {
  b0_rct <- 2.5
  b0_rwd <- 2.5
  
  sd_Y_rct <- 1.5
  sd_Y_rwd <- 1.8
  
} else if (sc == 3) {
  b0_rct <- 2.5
  b0_rwd <- 2.5
  
  sd_Y_rct <- 1.5
  sd_Y_rwd <- 1.8
}
# LRC-BART MODIFICATION END

beta_rct <- c(-0.50, -0.75, -0.50, -0.50,
              -1.35, -0.80, -0.01, -0.01, -0.01, -0.01)
beta_rwd <- beta_rct

b0_rct <- b0_rct - 2.58   # was -0.54; extra -2.04 offsets the added 2.0 bump (holds Y location; R2 unaffected)
b0_rwd <- b0_rwd - 2.58

gamma <- 0
eta <- 0
#===========================================================================
#============================ Generate RWD_pool ============================
#===========================================================================
gen_rwd_pool <- function(seed_base) {
  set.seed(seed_base)
  
  X_rwd_pool <- matrix(rep(NA, n_rwd_pool * p_obs), ncol = p_obs)
  if (sc == 3) {
    U_rwd_pool <- matrix(rep(NA, n_rwd_pool * length(conf_cols)), ncol = length(conf_cols))
  }
  
  for (i in 1:n_rwd_pool) {
    #========================== Generate X ==========================
    X_rwd_pool[i, 1:p_obs] <- rnorm(p_obs, mean = m_x[1:p_obs], sd = sd_x)
    
    for (col in extreme_cols) {
      if (runif(1) < 0.15) {
        X_rwd_pool[i, col] <- rnorm(1, mean = sample(c(-3, 3), 1), sd = 0.5)
      }
    }
    
    #========================== Generate U ==========================
    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        U_rwd_pool[i, j] <- rnorm(1, mean = m_x[conf_cols[j]], sd = sd_x)
      }
    }
  }
  
  #========================== Generate Y ==========================
  set.seed(seed_base + 1000)
  
  if (sc == 3) {
    W_rwd_pool <- U_rwd_pool
  } else {
    W_rwd_pool <- X_rwd_pool[, conf_cols]
  }
  
  y_rwd_pool <- rep(NA, n_rwd_pool)
  lp_rwd_pool <- rep(NA, n_rwd_pool)
  
  for (i in 1:n_rwd_pool) {
    x_vec <- X_rwd_pool[i, ]
    x_vec[conf_cols] <- W_rwd_pool[i, ]
    
    lp <- b0_rwd +
      beta_rwd[1] * x_vec[1]^2 +
      beta_rwd[2] * exp(-(x_vec[2] + 1)^2 / 2) +
      beta_rwd[3] * abs(x_vec[3] - 0.5) +
      beta_rwd[4] * abs(x_vec[4] + 1) +
      beta_rwd[5] * tanh(x_vec[5] - 2) +
      beta_rwd[6] * tanh(x_vec[6] + 2) +
      2.0 * ( (x_vec[5] - 2)^2 / (1 + ((x_vec[5] - 2)/1.5)^2) +
                (x_vec[6] + 2)^2 / (1 + ((x_vec[6] + 2)/1.5)^2) ) +
      beta_rwd[7] * x_vec[7] +
      beta_rwd[8] * x_vec[8] +
      beta_rwd[9] * x_vec[9] +
      beta_rwd[10] * x_vec[10]
    
    lp_rwd_pool[i] <- lp
    y_rwd_pool[i] <- lp + rnorm(n = 1, sd = sd_Y_rwd)
  }
  
  X_rwd_pool_by_c <- vector("list", length(cor))
  for (ci in seq_along(cor)) {
    Xc <- X_rwd_pool
    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        U_std <- (U_rwd_pool[,j] - mean(U_rwd_pool[,j])) / sd(U_rwd_pool[,j])
        set.seed(seed_base + j + 100)
        Xc[, conf_cols[j]] <- cor[ci] * U_std + sqrt(1 - cor[ci]^2) * rnorm(n_rwd_pool, mean = 0, sd = sd_U)
      }
    }
    X_rwd_pool_by_c[[ci]] <- Xc
  }
  
  for (.nm in c("X_rwd_pool","U_rwd_pool","C_rwd_pool","entry_rwd_pool","admin_rwd_pool",
                "W_rwd_pool","y_rwd_pool","lp_rwd_pool","yobs_rwd_pool","censor_rwd_pool",
                "label_rwd_pool","X_rwd_pool_by_c"))
    if (exists(.nm, inherits = FALSE)) assign(.nm, get(.nm), envir = .GlobalEnv)
}

# Frozen mode: build the pool once, up front (per-replicate mode rebuilds it inside the loop).
if (rwd_frozen) gen_rwd_pool(seed_rwd)

for (ci in seq_along(cor)) {
  c <- cor[ci]
  if (rwd_frozen) X_rwd_pool <- X_rwd_pool_by_c[[ci]]

  #===========================================================================
  #============================ Generate RCT_pool ============================
  #===========================================================================
  for (replicate_id in seq_len(n_replicates)){

    # Per-replicate mode: fresh RWD pool for each replicate (before set.seed(seed[ee]) so RCT is unchanged).
    if (!rwd_frozen) { gen_rwd_pool(seed_rwd_replicate[ee]); X_rwd_pool <- X_rwd_pool_by_c[[ci]] }

    set.seed(seed[ee])
    ee <- ee + 1

    #========================== Generate X ==========================
    X_rct_pool <- matrix(rep(NA, n_rct_pool * p_obs), ncol = p_obs)
    if (sc == 3) {
      U_rct_pool <- matrix(rep(NA, n_rct_pool * length(conf_cols)), ncol = length(conf_cols))
    }

    for (i in 1:n_rct_pool) {
      X_rct_pool[i, 1:p_obs] <- rnorm(p_obs, mean = m_x[1:p_obs], sd = sd_x)

      for (col in extreme_cols) {
        if (runif(1) < 0.15) {
          X_rct_pool[i, col] <- rnorm(1, mean = sample(c(-3, 3), 1), sd = 0.5)
        }
      }

      if (sc == 3) {
        for (j in 1:length(conf_cols)) {
          U_rct_pool[i, j] <- rnorm(1, mean = m_x[conf_cols[j]], sd = sd_x)
        }
      }
    }

    if (sc == 3) {
      for (j in 1:length(conf_cols)) {
        U_std <- (U_rct_pool[,j] - mean(U_rct_pool[,j])) / sd(U_rct_pool[,j])
        X_rct_pool[, conf_cols[j]] <- c * U_std + sqrt(1 - c^2) * rnorm(n_rct_pool, mean = 0, sd = sd_U)
        print(cor(X_rct_pool[, conf_cols[j]], U_rct_pool[,j]))
      }
    }

    # Combine RCT and RWD
    X <- rbind(X_rct_pool, X_rwd_pool)
    if (sc == 3) {
      U <- rbind(U_rct_pool, U_rwd_pool)
    }

    #========================== Generate D ==========================
    n_total_pool <- n_rct_pool + n_rwd_pool
    # LRC-BART MODIFICATION START
    if (sc %in% c(1, 4, 5)) {
      D <- rbern(n_total_pool, prob = (n1+n2)/(n1+n2+n3))
    }
    # LRC-BART MODIFICATION END
    if (sc == 2) {
      if (length(beta_D) > 1) {
        lp <- X[, conf_cols] %*% beta_D
      } else {
        lp <- beta_D * rowSums(X[, conf_cols])
      }
      thres <- quantile(lp, (n1+n2)/(n1+n2+n3))
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }
    if (sc == 3) {
      if (length(beta_D) > 1) {
        lp <- U[, 1:length(conf_cols)] %*% beta_D
      } else {
        lp <- beta_D * rowSums(U[, 1:length(conf_cols)])
      }
      thres <- quantile(lp, (n1+n2)/(n1+n2+n3))
      D <- rbern(n_total_pool, plogis((lp - thres) / sel_temp))
    }

    #===========================================================================
    #============================== Generate RCT ===============================
    #===========================================================================

    idx_rct_selected <- which(D == 1 & (1:n_total_pool) <= n_rct_pool)

    # LRC-BART MODIFICATION START
    if (length(idx_rct_selected) < (n1+n2))
      stop("RCT pool did not contain 300 source-eligible subjects")
    idx_rct <- sample(idx_rct_selected, n1+n2, replace = FALSE)
    # LRC-BART MODIFICATION END

    n_rct <- length(idx_rct)
    X_rct_sel <- X[idx_rct, , drop = FALSE]
    if (sc == 3) {
      W_rct_sel <- matrix(U[idx_rct, ], ncol = length(conf_cols))
    } else {
      W_rct_sel <- X_rct_sel[, conf_cols]
    }

    #========================== Generate Z ==========================
    Z_rct <- rbern(n_rct, prob = r/(r+1))

    #========================== Generate Y ==========================
    y_rct <- rep(NA, n_rct)
    lp_rct <- rep(NA, n_rct)

    true_mean_trt <- rep(NA, n_rct)
    true_mean_ctrl <- rep(NA, n_rct)
    eff_i_rct <- rep(NA, n_rct)

    for (i in 1:n_rct) {
      x_vec <- X_rct_sel[i, ]
      x_vec[conf_cols] <- W_rct_sel[i, ]
      modifier <- gamma * sum(x_vec[conf_cols])

      lp <- b0_rct +
        beta_rct[1] * x_vec[1]^2 +
        beta_rct[2] * exp(-(x_vec[2] + 1)^2 / 2) +
        beta_rct[3] * abs(x_vec[3] - 0.5) +
        beta_rct[4] * abs(x_vec[4] + 1) +
        beta_rct[5] * tanh(x_vec[5] - 2) +
        beta_rct[6] * tanh(x_vec[6] + 2) +
        2.0 * ( (x_vec[5] - 2)^2 / (1 + ((x_vec[5] - 2)/1.5)^2) +
                  (x_vec[6] + 2)^2 / (1 + ((x_vec[6] + 2)/1.5)^2) ) +
        beta_rct[7] * x_vec[7] +
        beta_rct[8] * x_vec[8] +
        beta_rct[9] * x_vec[9] +
        beta_rct[10] * x_vec[10]

      eff_i <- eff

      if(Z_rct[i] == 1){
        lp_rct[i] <- lp + modifier + eff_i
        y_rct[i] <- lp + modifier + eff_i + rnorm(n = 1, sd = sd_Y_rct)
      } else {
        lp_rct[i] <- lp + modifier
        y_rct[i] <- lp + modifier + rnorm(n = 1, sd = sd_Y_rct)
      }

      eff_i_rct[i] <- eff_i
      true_mean_trt[i] <- lp + modifier + eff_i
      true_mean_ctrl[i] <- lp + modifier
    }

    mean_trt_pop <- mean(true_mean_trt)
    mean_ctrl_pop <- mean(true_mean_ctrl)

    #===========================================================================
    #============================== Generate RWD ===============================
    #===========================================================================
    idx_rwd_selected <- which(D == 0 & (1:n_total_pool) > n_rct_pool)
    pool_indices_available <- idx_rwd_selected - n_rct_pool

    # LRC-BART MODIFICATION START
    if (length(pool_indices_available) < n3)
      stop("RWD pool did not contain 300 source-eligible subjects")
    pool_indices <- sample(pool_indices_available, n3, replace = FALSE)
    # LRC-BART MODIFICATION END

    n_rwd <- length(pool_indices)
    X_rwd <- X_rwd_pool[pool_indices, , drop = FALSE]
    if (sc == 3) {
      U_rwd <- U_rwd_pool[pool_indices, , drop = FALSE]
    }

    y_rwd <- y_rwd_pool[pool_indices]
    lp_rwd <- lp_rwd_pool[pool_indices]
    Z_rwd <- rep(0, n_rwd)

    # LRC-BART ADDITION START
    # The shift is deterministic and is applied after all outcome noise has
    # been drawn, so matched scenarios consume the same random-number stream.
    if (region == "X5") {
      region_rwd <- X_rwd[, 5] > 2
      region_rct <- X_rct_sel[, 5] > 2
    } else if (region == "X7") {
      region_rwd <- X_rwd[, 7] <= 2
      region_rct <- X_rct_sel[, 7] <= 2
    } else {
      region_rwd <- rep(FALSE, n_rwd)
      region_rct <- rep(FALSE, n_rct)
    }
    shift_rwd <- if (sc == 4) rep(delta_rwd, n_rwd) else
      if (sc == 5) delta_rwd * (!region_rwd) else rep(0, n_rwd)
    shift_at_rct <- if (sc == 4) rep(delta_rwd, n_rct) else
      if (sc == 5) delta_rwd * (!region_rct) else rep(0, n_rct)
    y_rwd <- y_rwd + shift_rwd
    lp_rwd <- lp_rwd + shift_rwd
    true_mean_rwd_at_rct <- true_mean_ctrl + shift_at_rct
    # LRC-BART MODIFICATION START
    this_scenario_id <- paste0("sc", sc)
    # LRC-BART MODIFICATION END
    if (sc == 3)
      this_scenario_id <- paste0(this_scenario_id, "_cor", c)
    if (sc == 4)
      this_scenario_id <- paste0(this_scenario_id, "_d", delta_rwd)
    if (sc == 5)
      this_scenario_id <- paste0(
        this_scenario_id, "_", region, "_d", delta_rwd
      )
    # LRC-BART ADDITION END

    n_rct_trt <- sum(Z_rct == 1)
    n_rct_ctrl <- sum(Z_rct == 0)
    cat("n_rct_trt =", n_rct_trt, ", n_rct_ctrl =", n_rct_ctrl, ", n_rwd =", n_rwd,"\n")

    X_rct_df <- data.frame(X_rct_sel[, 1:p_obs])
    colnames(X_rct_df) <- paste0("X", 1:p_obs)
    X_rwd_df <- data.frame(X_rwd[, 1:p_obs])
    colnames(X_rwd_df) <- paste0("X", 1:p_obs)
    X_rct_df$D <- 1
    X_rwd_df$D <- 0
    X_rct_df$Z <- Z_rct
    X_rwd_df$Z <- Z_rwd

    if (sc == 3) {
      U_rct_df <- data.frame(W_rct_sel)
      colnames(U_rct_df) <- paste0("U", 1:length(conf_cols))
      U_rwd_df <- data.frame(U_rwd)
      colnames(U_rwd_df) <- paste0("U", 1:length(conf_cols))
      U_combined <- rbind(U_rct_df, U_rwd_df)
    } else {
      U_combined <- NULL
    }

    cat("eff =", eff,
        "mean(eff_i_rct) =", mean(eff_i_rct),
        "eff_true =", mean_trt_pop - mean_ctrl_pop,"\n")

    # Output
    out <- list(X = rbind(X_rct_df, X_rwd_df),
                U = U_combined,
                y = c(y_rct, y_rwd),
                treat_eff = eff,
                treat_eff_true = mean_trt_pop - mean_ctrl_pop,
                treat_eff_star = eff_star,
                gamma_eff = gamma,  # modifier coefficient
                eta_eff = if (exists("eta")) eta else 0,  # HTE coefficient for treatment effect
                sigma_rct = sd_Y_rct,
                sigma_rwd = sd_Y_rwd,
                true_mean_trt = true_mean_trt,
                true_mean_ctrl = true_mean_ctrl,
                true_mean_trt_pop = mean_trt_pop,
                true_mean_ctrl_pop = mean_ctrl_pop,
                eff_i = eff_i_rct,
                lp = c(lp_rct, lp_rwd),
                n_rct = n_rct,
                n_rwd = n_rwd,
                # LRC-BART ADDITION START
                scenario_id = this_scenario_id,
                delta_rwd = delta_rwd,
                region = region,
                region_rct = region_rct,
                true_mean_rwd_at_rct = true_mean_rwd_at_rct
                # LRC-BART ADDITION END
                )

    # Save with scenario identifiers
    # LRC-BART MODIFICATION START
    saveRDS(out, file = file.path(
      projectDir, "data",
      paste0("data_p", p_obs, "_", this_scenario_id, "_", hypo,
             if (rwd_frozen) "_fz" else "", "_", replicate_id, ".RData")
    ))
    # LRC-BART MODIFICATION END
  }
}

# LRC-BART MODIFICATION START
# The inherited trailing "Testing code" block is intentionally omitted. The
# production generator ends after writing the requested replicate datasets.
# LRC-BART MODIFICATION END
