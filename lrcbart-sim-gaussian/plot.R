rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
# LRC-BART MODIFICATION START
# The inherited Sc1-Sc3 result plot is extended below with the revision's Sc4
# global-shift and Sc5 regional-partial-compatibility configurations. All
# result paths point to this migrated project; the former MAP-BART series is
# read from the default LRC-BART target files.
mainDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/"
# LRC-BART MODIFICATION END

p_obs <- 10
# Pipeline selection; standalone runs use the same alternative-hypothesis default.
plot_context <- NULL
# ---- Local plot selection and result readers ----
# Plot selection is independent of which comparison methods were fitted today.
`%||%` <- function(x, default) if (is.null(x)) default else x
plot_context <- plot_context %||% list()
plot_hypothesis <- plot_context$hypothesis %||% "alternative"
stopifnot(length(plot_hypothesis) == 1L, plot_hypothesis %in% c("null", "alternative"))
plot_saved <- character()

plot_scenario_id <- function(cfg) {
  tag <- paste0("sc", cfg$sc)
  if (cfg$sc == 3) tag <- paste0(tag, "_cor", cfg$cor)
  if (cfg$sc == 4) tag <- paste0(tag, "_d", cfg$delta_rwd)
  if (cfg$sc == 5) tag <- paste0(tag, "_", cfg$region, "_d", cfg$delta_rwd)
  tag
}
selected_sc <- function(default) {
  if (is.null(plot_context$scenarios)) return(default)
  unique(vapply(plot_context$scenarios, function(x) as.integer(x$sc), integer(1)))
}
selected_configurations <- function(sc, default) {
  if (is.null(plot_context$scenarios) || sc <= 3) return(default)
  chosen <- Filter(function(x) x$sc == sc, plot_context$scenarios)
  vapply(chosen, function(x) if (sc == 4) paste0("d", x$delta_rwd) else
    paste0(x$region, "_d", x$delta_rwd), character(1))
}
selected_correlations <- function(sc, default) {
  if (is.null(plot_context$scenarios) || sc != 3) return(default)
  unique(vapply(Filter(function(x) x$sc == sc, plot_context$scenarios), `[[`, numeric(1), "cor"))
}
plot_method_setting <- function(method, name, default)
  plot_context$methods[[method]][[name]] %||% default

# Preserve the familiar ordering while retaining custom scenario settings.
plot_factor <- function(x, levels) {
  actual <- unique(as.character(x[!is.na(x)]))
  factor(x, levels = unique(c(levels[levels %in% actual], actual)))
}

plot_path_allowed <- function(path, check_failures = FALSE) {
  name <- basename(path)
  if (!grepl(paste0("_", plot_hypothesis, "(_|\\.)"), name)) return(FALSE)
  if (!is.null(plot_context$scenarios)) {
    tags <- vapply(plot_context$scenarios, plot_scenario_id, character(1))
    if (!any(vapply(tags, function(tag) grepl(paste0("_", tag, "_"), name, fixed = TRUE), logical(1))))
      return(FALSE)
  }
  if (!is.null(plot_context$n_T) && !is.na(plot_context$n_T) &&
      !grepl(paste0("_n", plot_context$n_T, "_"), name, fixed = TRUE)) return(FALSE)
  if (check_failures) for (failure in plot_context$failed_results) {
    method <- if (failure$method == "lrcBART") "LRC-BART" else failure$method
    if (!startsWith(name, paste0(method, "_")) ||
        !grepl(paste0("_", failure$scenario, "_"), name, fixed = TRUE) ||
        !grepl(paste0("_", failure$hypothesis, "_"), paste0(sub("\\.RData$", "", name), "_"), fixed = TRUE)) next
    if (!is.na(failure$n_T) && !grepl(paste0("_n", failure$n_T, "_"), name, fixed = TRUE)) next
    if (method != "LRC-BART") return(FALSE)
    tail <- strsplit(name, paste0("_", failure$scenario, "_"), fixed = TRUE)[[1]][2]
    if (failure$config == "default" && startsWith(tail, "N") && !startsWith(tail, "Ns0min_")) return(FALSE)
    if (failure$config == "s0min" && startsWith(tail, "Ns0min_")) return(FALSE)
    if (startsWith(tail, paste0(failure$config, "_N"))) return(FALSE)
  }
  TRUE
}
read_plot_data <- function(path) {
  if (!plot_path_allowed(path)) stop("Dataset excluded by the current plot selection")
  readRDS(path)
}
read_plot_result <- function(path) {
  if (!plot_path_allowed(path, TRUE)) stop("Result excluded by the current plot selection or a failed fit")
  result <- readRDS(path)
  if (is.data.frame(result) && "iteration" %in% names(result) && !is.null(plot_context$n_replicates))
    result <- result[!is.na(result$iteration) & result$iteration <= plot_context$n_replicates, , drop = FALSE]
  if (is.data.frame(result) && !nrow(result)) stop("No requested replicates in result")
  if (plot_hypothesis == "null" && is.data.frame(result)) {
    if ("fp" %in% names(result)) result$tp <- result$fp
    if ("fp_rmst" %in% names(result)) result$tp_rmst <- result$fp_rmst
  }
  result
}
save_pipeline_plot <- function(filename, ...) {
  if (plot_hypothesis == "null") filename <- sub("(\\.[^.]+)$", "_null\\1", filename)
  if (!is.null(plot_context$output_dir)) filename <- file.path(plot_context$output_dir, basename(filename))
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, ...)
  plot_saved <<- c(plot_saved, filename)
  invisible(filename)
}
pipeline_table <- function(d, ...) {
  if (plot_hypothesis == "null") {
    names(d)[names(d) == "Power"] <- "Type I error"
    names(d)[names(d) == "Power_calib"] <- "Calibrated Type I error"
    # Optional null-file columns can otherwise duplicate the main null summary.
    d <- d[, !duplicated(names(d)), drop = FALSE]
  }
  gridExtra::tableGrob(d, ...)
}
finish_pipeline_plot <- function() {
  if (!length(plot_saved)) stop("No figures were generated for the requested plot selection")
  invisible(plot_saved)
}
lrc_plot_configurations <- function(single_arm, sc) {
  defaults <- if (single_arm) c("default", "w0.9", "s0min") else c("default", "w1", "w0")
  configs <- plot_context$lrc_configs %||% defaults
  # Match the runner's sensitivity scope, including in standalone plot runs.
  sensitivity_scenarios <- plot_context$lrc_sensitivity_scenarios %||% 1:3
  if (!single_arm && !sc %in% sensitivity_scenarios)
    configs <- configs[configs == "default"]
  do.call(rbind, lapply(configs, function(config) {
    if (single_arm) {
      stopifnot(config %in% c("default", "w0.9", "s0min"))
      target <- if (config == "s0min") "s0min" else c("100", "f90")
      suffix <- if (config == "w0.9") "_w0.9" else ""
      label <- if (config == "s0min") "LRC-BART(s₀²=min)" else
        paste0("LRC-BART(N=", target, ",w=", if (config == "w0.9") "0.9" else "1", ")")
    } else {
      stopifnot(config %in% c("default", "w0", "w1", "Hg10"))
      target <- if (config == "default") c("100", "75", "50", "f90", "f50", "f25") else "100"
      suffix <- if (config == "default") "" else paste0("_", config)
      label <- paste0("LRC-BART(N=", target,
        if (config == "default") "" else paste0(",", if (config == "Hg10") "Hg=10" else sub("w", "w=", config)), ")")
    }
    data.frame(config = config, suffix = suffix, target = target, label = label, stringsAsFactors = FALSE)
  }))
}
weight_labels <- function(method, weights) {
  if (identical(as.numeric(weights), 1)) method else paste0(method, "(w=", weights, ")")
}
# ---- End local plot helpers ----


# LRC-BART MODIFICATION START
# hierLM still keys on a single calibrated prior value.
lmv4_prior_vals <- plot_method_setting("hierLM.R", "prior_vals", 0.25)

# MAP and lrcBART now key on N_target rather than the raw s^2 value: each
# analysis file loops over multiple N_target multipliers from its
# calibration .RData and writes one result file per target (with the
# integer target N in the filename).  List the N values to read here.
map_target_Ns <- plot_method_setting("MAP.R", "target_Ns", c(100, 75, 50))
# LRC-BART targets and fixed-weight variants are selected inside each scenario.
# LRC-BART MODIFICATION END

# Helper function to add separator rows between groups
add_configuration_separators <- function(df, configuration_col = "Configuration_Label") {
  if (nrow(df) == 0 || !configuration_col %in% colnames(df)) return(df)

  # Get unique values in order of appearance
  unique_vals <- unique(as.character(df[[configuration_col]]))

  result <- data.frame()
  for (i in seq_along(unique_vals)) {
    group_data <- df[as.character(df[[configuration_col]]) == unique_vals[i], ]
    result <- rbind(result, group_data)

    # Add separator row after each group (except the last)
    if (i < length(unique_vals)) {
      sep_row <- group_data[1, ]
      sep_row[] <- ""
      result <- rbind(result, sep_row)
    }
  }

  return(result)
}

for (sc in selected_sc(1:5)) {
lrcbart_configurations <- lrc_plot_configurations(FALSE, sc)
tryCatch({

if (sc %in% c(1, 2, 4, 5)){

  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # LRC-BART ADDITION START
  # Sc1/Sc2 use the Base label. Sc4 rows distinguish the global
  # RWD shift, while Sc5 rows distinguish both compatible region and shift.
  if (sc %in% c(1, 2)) configurations <- c("Base")
  if (sc == 4) configurations <- c("d1", "d2")
  if (sc == 5)
    configurations <- c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2")
  # LRC-BART ADDITION END

  for (configuration in selected_configurations(sc, configurations)) {

    # Initialize empty lists for this configuration
    configuration_res_list_ATE <- list()
    configuration_res_list_sigma <- list()

    # LRC-BART ADDITION START
    # Construct the exact scenario tag written by the migrated method files.
    if (sc %in% c(1, 2)) {
      scenario_suffix <- paste0("_sc", sc)
      configuration_label <- "Default"
    }
    if (sc == 4) {
      scenario_suffix <- paste0("_sc4_", configuration)
      configuration_label <- paste0("delta RWD = ", sub("^d", "", configuration))
    }
    if (sc == 5) {
      scenario_suffix <- paste0("_sc5_", configuration)
      configuration_parts <- strsplit(configuration, "_d", fixed = TRUE)[[1]]
      configuration_label <- paste0("region = ", configuration_parts[1],
                            ", delta RWD = ", configuration_parts[2])
    }
    file_suffix <- paste0(scenario_suffix, paste0("_", plot_hypothesis, ".RData"))
    # LRC-BART ADDITION END

    # Try to read LMv1 results
    tryCatch({
      LMv1_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-gaussian/res/LMv1_p",p_obs,file_suffix))
      LMv1_res_ATE <- LMv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv1_res_sigma <- LMv1_res[, c("c","iteration",
                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                "pehe_subj","bias_subj")]
      LMv1_res_ATE$Method <- "LMv1"
      LMv1_res_sigma$Method <- "LMv1"
      LMv1_res_ATE$Configuration <- configuration_label
      LMv1_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv1_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv1_res_sigma
    }, error = function(e) {
    })

    # Try to read LMv2 results
    tryCatch({
      LMv2_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-gaussian/res/LMv2_p",p_obs,file_suffix))
      LMv2_res_ATE <- LMv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv2_res_sigma <- LMv2_res[, c("c","iteration",
                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                  "pehe_subj","bias_subj")]
      LMv2_res_ATE$Method <- "LMv2"
      LMv2_res_sigma$Method <- "LMv2"
      LMv2_res_ATE$Configuration <- configuration_label
      LMv2_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv2_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv2_res_sigma
    }, error = function(e) {
    })

    # Try to read LMv3 results
    tryCatch({
      LMv3_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-gaussian/res/LMv3_p",p_obs,file_suffix))
      LMv3_res_ATE <- LMv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv3_res_sigma <- LMv3_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      LMv3_res_ATE$Method <- "LMv3"
      LMv3_res_sigma$Method <- "LMv3"
      LMv3_res_ATE$Configuration <- configuration_label
      LMv3_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv3_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv3_res_sigma
    }, error = function(e) {
    })

    # Try to read hierLM results for each prior value
    for (prior_val in lmv4_prior_vals) {
      tryCatch({
        # Construct file path with prior.
        lmv4_file <- paste0(mainDir,"/lrcbart-sim-gaussian/res/hierLM_p",p_obs,
                            scenario_suffix,"_prior",prior_val,
                            paste0("_", plot_hypothesis, ".RData"))

        if (file.exists(lmv4_file)) {
          hierLM_res <- read_plot_result(lmv4_file)
          hierLM_res_ATE <- hierLM_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          hierLM_res_sigma <- hierLM_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          hierLM_res_ATE$Method <- paste0("hierLM(", prior_val, ")")
          hierLM_res_sigma$Method <- paste0("hierLM(", prior_val, ")")
          hierLM_res_ATE$Configuration <- configuration_label
          hierLM_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- hierLM_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- hierLM_res_sigma
        }
      }, error = function(e) {
      })
    }

    # Try to read MAP results for each N_target value (ATE only, no sigma data)
    for (target_N in map_target_Ns) {
      tryCatch({
        # Construct file path with target N.
        map_file <- paste0(mainDir,"/lrcbart-sim-gaussian/res/MAP_p",p_obs,
                           scenario_suffix,"_N",target_N,
                           paste0("_", plot_hypothesis, ".RData"))

        if (file.exists(map_file)) {
          MAP_res <- read_plot_result(map_file)
          MAP_res_ATE <- MAP_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          MAP_res_sigma <- MAP_res[, c("c","iteration",
                                      "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                      "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                      "bias.trt.median.pop","w2distance.trt.median.pop",
                                      "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                      "pehe_subj","bias_subj")]
          MAP_res_ATE$Method <- paste0("MAP(N=", target_N, ")")
          MAP_res_sigma$Method <- paste0("MAP(N=", target_N, ")")
          MAP_res_ATE$Configuration <- configuration_label
          MAP_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- MAP_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- MAP_res_sigma
        }
      }, error = function(e) {
      })
    }

    # Try to read PSCL results (ATE only, no sigma data, may have NA values)
    # Keep all rows including NAs - plotting functions will handle with na.rm = TRUE
    tryCatch({
      PSCL_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-gaussian/res/PSCL_p",p_obs,file_suffix))
      PSCL_res_ATE <- PSCL_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      PSCL_res_ATE$Method <- "PSCL"
      PSCL_res_ATE$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- PSCL_res_ATE
      # PSCL doesn't have sigma data, so no sigma entry
    }, error = function(e) {
    })

    # Try to read BARTv1 / BARTv2 results
    for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
    tryCatch({
      BART_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-gaussian/res/",bver,"_p",p_obs,file_suffix))
      BART_res <- BART_res[BART_res$alpha == 0.95 & BART_res$beta == 2, ]
      BART_res_ATE <- BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      BART_res_sigma <- BART_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      BART_res_ATE$Method <- bver
      BART_res_sigma$Method <- bver
      BART_res_ATE$Configuration <- configuration_label
      BART_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- BART_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- BART_res_sigma
    }, error = function(e) {
    })
    }

    # Try to read lrcBART results for each N_target value
    for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

      tryCatch({
        # Construct file path with target N.
        lrcbart_file <- paste0(mainDir,"/lrcbart-sim-gaussian/res/LRC-BART_p",p_obs,
                               scenario_suffix,lrc_suffix,"_N",target_N,
                               paste0("_", plot_hypothesis, ".RData"))

        if (file.exists(lrcbart_file)) {
          lrcBART_res <- read_plot_result(lrcbart_file)
          lrcBART_res <- lrcBART_res[lrcBART_res$alpha == 0.95 & lrcBART_res$beta == 2, ]
          lrcBART_res_ATE <- lrcBART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          lrcBART_res_sigma <- lrcBART_res[, c("c","iteration",
                                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                "pehe_subj","bias_subj")]
          lrcBART_res_ATE$Method <- lrc_label
          lrcBART_res_sigma$Method <- lrc_label
          lrcBART_res_ATE$Configuration <- configuration_label
          lrcBART_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- lrcBART_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- lrcBART_res_sigma
        }
      }, error = function(e) {
        # Silently skip if file doesn't exist
      })
    }

    # Combine results for this configuration (only if we have data)
    if (length(configuration_res_list_ATE) > 0) {
      configuration_res_ATE <- do.call(rbind, configuration_res_list_ATE)
      all_res_ATE <- rbind(all_res_ATE, configuration_res_ATE)
    }

    if (length(configuration_res_list_sigma) > 0) {
      configuration_res_sigma <- do.call(rbind, configuration_res_list_sigma)
      all_res_sigma <- rbind(all_res_sigma, configuration_res_sigma)
    }
  }

  # Assign to res and res_sigma for plotting
  res <- all_res_ATE
  res_sigma <- all_res_sigma

  # Load Type I error from null hypothesis files (if current hypo is not null)
  all_null_FP <- data.frame()
  for (configuration in selected_configurations(sc, configurations)) {
    # LRC-BART ADDITION START
    # Match the alternative-file configuration tag when optional null files
    # are present.
    if (sc %in% c(1, 2)) {
      scenario_suffix <- paste0("_sc", sc)
      configuration_label <- "Default"
    }
    if (sc == 4) {
      scenario_suffix <- paste0("_sc4_", configuration)
      configuration_label <- paste0("delta RWD = ", sub("^d", "", configuration))
    }
    if (sc == 5) {
      scenario_suffix <- paste0("_sc5_", configuration)
      configuration_parts <- strsplit(configuration, "_d", fixed = TRUE)[[1]]
      configuration_label <- paste0("region = ", configuration_parts[1],
                            ", delta RWD = ", configuration_parts[2])
    }
    null_file_suffix <- paste0(scenario_suffix, "_null.RData")
    # LRC-BART ADDITION END

    # Try to load null results for each method
    methods_list <- list(
      list(name = "LMv1", file_prefix = "LMv1"),
      list(name = "LMv2", file_prefix = "LMv2"),
      list(name = "LMv3", file_prefix = "LMv3"),
      list(name = "hierLM", file_prefix = "hierLM"),
      list(name = "MAP", file_prefix = "MAP"),
      list(name = "PSCL", file_prefix = "PSCL"),
      list(name = "BARTv1", file_prefix = "BARTv1"),
      list(name = "BARTv2", file_prefix = "BARTv2"),
      list(name = "BARTv3", file_prefix = "BARTv3"),
      list(name = "lrcBART", file_prefix = "LRC-BART")
    )

    for (method_info in methods_list) {
      tryCatch({
        null_file <- paste0(mainDir,"/lrcbart-sim-gaussian/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
        if (file.exists(null_file)) {
          null_res <- read_plot_result(null_file)

          # For BART and lrcBART, filter by alpha and beta
          if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3", "lrcBART")) {
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          }

          # Calculate mean FP (Type I error) - use na.rm = TRUE for PSCL which may have NAs
          if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              Configuration = configuration_label,
              Method = method_info$name,
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {
        # Silently skip if null file doesn't exist
      })
    }
  }
}
if (sc == 3){
  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # Loop through the base configuration and correlations
  configurations <- c("Base")
  cor_vals <- c(-0.5, 0, 0.5)

  # Track which files were found
  files_checked <- list()
  files_found <- list()

  for (configuration in selected_configurations(sc, configurations)) {
    for (cor_val in selected_correlations(sc, cor_vals)) {

      # Initialize empty lists for this configuration and correlation
      configuration_res_list_ATE <- list()
      configuration_res_list_sigma <- list()

      # Construct file suffix based on configuration
      # For sc == 3: include "cor" prefix before correlation value
      if (configuration == "Base") {
        file_suffix <- paste0("_sc", sc, "_cor", cor_val, paste0("_", plot_hypothesis, ".RData"))
        configuration_label <- configuration
      } else {
        file_suffix <- paste0("_sc", sc, configuration, "_cor", cor_val, paste0("_", plot_hypothesis, ".RData"))
        configuration_label <- configuration
      }

      # Try to read LMv1 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-gaussian/res/LMv1_p",p_obs,file_suffix)
        files_checked[[length(files_checked) + 1]] <- file_path
        if (file.exists(file_path)) {
          LMv1_res <- read_plot_result(file_path)
          files_found[[length(files_found) + 1]] <- file_path
          LMv1_res_ATE <- LMv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv1_res_sigma <- LMv1_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          LMv1_res_ATE$Method <- "LMv1"
          LMv1_res_sigma$Method <- "LMv1"
          LMv1_res_ATE$Configuration <- configuration_label
          LMv1_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv1_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv1_res_sigma
        }
      }, error = function(e) {
      })

      # Try to read LMv2 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-gaussian/res/LMv2_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          LMv2_res <- read_plot_result(file_path)
          LMv2_res_ATE <- LMv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv2_res_sigma <- LMv2_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          LMv2_res_ATE$Method <- "LMv2"
          LMv2_res_sigma$Method <- "LMv2"
          LMv2_res_ATE$Configuration <- configuration_label
          LMv2_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv2_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv2_res_sigma
        }
      }, error = function(e) {
      })

      # Try to read LMv3 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-gaussian/res/LMv3_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          LMv3_res <- read_plot_result(file_path)
          LMv3_res_ATE <- LMv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LMv3_res_sigma <- LMv3_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          LMv3_res_ATE$Method <- "LMv3"
          LMv3_res_sigma$Method <- "LMv3"
          LMv3_res_ATE$Configuration <- configuration_label
          LMv3_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv3_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv3_res_sigma
        }
      }, error = function(e) {
      })

      # Try to read hierLM results for each prior value
      for (prior_val in lmv4_prior_vals) {
        tryCatch({
          # Construct file path with prior
          if (configuration == "Base") {
            lmv4_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/hierLM_p",p_obs,"_sc",sc,"_cor",cor_val,"_prior",prior_val,paste0("_", plot_hypothesis, ".RData"))
          } else {
            lmv4_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/hierLM_p",p_obs,"_sc",sc,configuration,"_cor",cor_val,"_prior",prior_val,paste0("_", plot_hypothesis, ".RData"))
          }

          if (file.exists(lmv4_file)) {
            hierLM_res <- read_plot_result(lmv4_file)
            hierLM_res_ATE <- hierLM_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            hierLM_res_sigma <- hierLM_res[, c("c","iteration",
                                          "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                          "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                          "bias.trt.median.pop","w2distance.trt.median.pop",
                                          "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                          "pehe_subj","bias_subj")]
            hierLM_res_ATE$Method <- paste0("hierLM(", prior_val, ")")
            hierLM_res_sigma$Method <- paste0("hierLM(", prior_val, ")")
            hierLM_res_ATE$Configuration <- configuration_label
            hierLM_res_sigma$Configuration <- configuration_label
            configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- hierLM_res_ATE
            configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- hierLM_res_sigma
          }
        }, error = function(e) {
        })
      }

      # Try to read MAP results for each N_target value (ATE only, no sigma data)
      for (target_N in map_target_Ns) {
        tryCatch({
          # Construct file path with target N
          if (configuration == "Base") {
            map_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,"_cor",cor_val,"_N",target_N,paste0("_", plot_hypothesis, ".RData"))
          } else {
            map_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/MAP_p",p_obs,"_sc",sc,configuration,"_cor",cor_val,"_N",target_N,paste0("_", plot_hypothesis, ".RData"))
          }

          if (file.exists(map_file)) {
            MAP_res <- read_plot_result(map_file)
            MAP_res_ATE <- MAP_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            MAP_res_sigma <- MAP_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
            MAP_res_ATE$Method <- paste0("MAP(N=", target_N, ")")
            MAP_res_sigma$Method <- paste0("MAP(N=", target_N, ")")
            MAP_res_ATE$Configuration <- configuration_label
            MAP_res_sigma$Configuration <- configuration_label
            configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- MAP_res_ATE
            configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- MAP_res_sigma
          }
        }, error = function(e) {
        })
      }

      # Try to read PSCL results (ATE only, no sigma data, may have NA values)
      # Keep all rows including NAs - plotting functions will handle with na.rm = TRUE
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-gaussian/res/PSCL_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          PSCL_res <- read_plot_result(file_path)
          PSCL_res_ATE <- PSCL_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          PSCL_res_ATE$Method <- "PSCL"
          PSCL_res_ATE$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- PSCL_res_ATE
          # PSCL doesn't have sigma data, so no sigma entry
        }
      }, error = function(e) {
      })

      # Try to read BARTv1 / BARTv2 results
      for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-gaussian/res/",bver,"_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          BART_res <- read_plot_result(file_path)
          BART_res <- BART_res[BART_res$alpha == 0.95 & BART_res$beta == 2, ]
          BART_res_ATE <- BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          BART_res_sigma <- BART_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          BART_res_ATE$Method <- bver
          BART_res_sigma$Method <- bver
          BART_res_ATE$Configuration <- configuration_label
          BART_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- BART_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- BART_res_sigma
        }
      }, error = function(e) {
      })
      }

      # Try to read lrcBART results for each N_target value
      for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

        tryCatch({
          # Construct file path with target N
          if (configuration == "Base") {
            lrcbart_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/LRC-BART_p",p_obs,"_sc",sc,"_cor",cor_val,lrc_suffix,"_N",target_N,paste0("_", plot_hypothesis, ".RData"))
          } else {
            lrcbart_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/LRC-BART_p",p_obs,"_sc",sc,configuration,"_cor",cor_val,lrc_suffix,"_N",target_N,paste0("_", plot_hypothesis, ".RData"))
          }

          if (file.exists(lrcbart_file)) {
            lrcBART_res <- read_plot_result(lrcbart_file)
            lrcBART_res <- lrcBART_res[lrcBART_res$alpha == 0.95 & lrcBART_res$beta == 2, ]
            lrcBART_res_ATE <- lrcBART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            lrcBART_res_sigma <- lrcBART_res[, c("c","iteration",
                                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                  "pehe_subj","bias_subj")]
            lrcBART_res_ATE$Method <- lrc_label
            lrcBART_res_sigma$Method <- lrc_label
            lrcBART_res_ATE$Configuration <- configuration_label
            lrcBART_res_sigma$Configuration <- configuration_label
            configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- lrcBART_res_ATE
            configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- lrcBART_res_sigma
          }
        }, error = function(e) {
          # Silently skip if file doesn't exist
        })
      }

      # Combine results for this configuration and correlation (only if we have data)
      if (length(configuration_res_list_ATE) > 0) {
        configuration_res_ATE <- do.call(rbind, configuration_res_list_ATE)
        all_res_ATE <- rbind(all_res_ATE, configuration_res_ATE)
      }

      if (length(configuration_res_list_sigma) > 0) {
        configuration_res_sigma <- do.call(rbind, configuration_res_list_sigma)
        all_res_sigma <- rbind(all_res_sigma, configuration_res_sigma)
      }
    }
  }

  # Assign to res and res_sigma for plotting
  if (nrow(all_res_ATE) > 0) {
    res <- all_res_ATE
  } else {
    stop("No results files found for scenario 3")
  }

  if (nrow(all_res_sigma) > 0) {
    res_sigma <- all_res_sigma
  } else {
    res_sigma <- data.frame()
  }

  # Load Type I error from null hypothesis files
  all_null_FP <- data.frame()
  methods_list <- list(
    list(name = "LMv1", file_prefix = "LMv1"),
    list(name = "LMv2", file_prefix = "LMv2"),
    list(name = "LMv3", file_prefix = "LMv3"),
    list(name = "hierLM", file_prefix = "hierLM"),
    list(name = "MAP", file_prefix = "MAP"),
    list(name = "PSCL", file_prefix = "PSCL"),
    list(name = "BARTv1", file_prefix = "BARTv1"),
    list(name = "BARTv2", file_prefix = "BARTv2"),
    list(name = "BARTv3", file_prefix = "BARTv3"),
    list(name = "lrcBART", file_prefix = "LRC-BART")
  )

  for (configuration in selected_configurations(sc, configurations)) {
    for (cor_val in selected_correlations(sc, cor_vals)) {
      # Construct null file suffix based on configuration
      if (configuration == "Base") {
        null_file_suffix <- paste0("_sc", sc, "_cor", cor_val, "_null.RData")
        configuration_label <- configuration
      } else {
        null_file_suffix <- paste0("_sc", sc, configuration, "_cor", cor_val, "_null.RData")
        configuration_label <- configuration
      }

      for (method_info in methods_list) {
        tryCatch({
          null_file <- paste0(mainDir,"lrcbart-sim-gaussian/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
          if (file.exists(null_file)) {
            null_res <- read_plot_result(null_file)

            # For BART and lrcBART, filter by alpha and beta
            if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3", "lrcBART")) {
              null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            }

            # Calculate mean FP (Type I error) - use na.rm = TRUE for PSCL which may have NAs
            if ("fp" %in% colnames(null_res) && nrow(null_res) > 0) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                Configuration = configuration_label,
                c = cor_val,
                Method = method_info$name,
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {
          # Silently skip if null file doesn't exist
        })
      }
    }
  }
}

# Check if we have any data to plot
if (nrow(res) == 0) {
  stop("No results data found. Please check that result files exist and paths are correct.")
}

# Create plotting-row labels for all scenarios.
if (sc %in% 1:5) {
  configuration_levels <- "Default"
  if (sc == 4)
    configuration_levels <- paste0("delta RWD = ", c(1, 2))
  if (sc == 5)
    configuration_levels <- c(
      paste0("region = X5, delta RWD = ", c(0.5, 1, 2)),
      paste0("region = X7, delta RWD = ", c(1, 2))
    )
  if (nrow(res) > 0 && "Configuration" %in% colnames(res)) {
    res$Configuration_Label <- if (sc <= 3)
      plot_factor(res$Configuration, levels = configuration_levels) else
      plot_factor(res$Configuration, levels = configuration_levels)
  }
  if (nrow(res_sigma) > 0 && "Configuration" %in% colnames(res_sigma)) {
    res_sigma$Configuration_Label <- if (sc <= 3)
      plot_factor(res_sigma$Configuration, levels = configuration_levels) else
      plot_factor(res_sigma$Configuration, levels = configuration_levels)
  }
}

# Create method levels dynamically.  hierLM still labels by prior value;
# MAP and LRC-BART label by N_target.
lmv4_methods  <- paste0("hierLM(",      lmv4_prior_vals, ")")
map_methods   <- paste0("MAP(N=",     map_target_Ns,   ")")
lrcbart_methods <- lrcbart_configurations$label
all_method_levels <- c("LMv1", "LMv2", "LMv3", "BARTv1", "BARTv2", "BARTv3", "PSCL", map_methods, lmv4_methods, lrcbart_methods)

if (nrow(res) > 0) {
  res$Method <- factor(res$Method, levels = all_method_levels)
}
if (nrow(res_sigma) > 0) {
  res_sigma$Method <- factor(res_sigma$Method, levels = all_method_levels)
}

# Define specific colors for each method (default ggplot2 colors)
n_methods <- length(all_method_levels)
default_colors <- scales::hue_pal()(n_methods)
method_colors <- setNames(default_colors, all_method_levels)

if (sc %in% c(1, 2, 4, 5)){

  summary_table <- res %>%
    filter(!is.na(bias)) %>%
    group_by(Configuration_Label, Configuration, Method) %>%
    summarise(
      Bias = sprintf("%.2f", mean(bias, na.rm = TRUE)),
      SD = sprintf("%.2f", mean(sd, na.rm = TRUE)),
      RMSE = sprintf("%.2f", sqrt(mean(rmse, na.rm = TRUE))),
      W1Distance = sprintf("%.2f", mean(w1distance, na.rm = TRUE)),
      W2Distance = sprintf("%.2f", mean(w2distance, na.rm = TRUE)),
      CI_length = sprintf("%.2f", mean(ci, na.rm = TRUE)),
      CI_coverage = sprintf("%.2f", mean(coverage, na.rm = TRUE)),
      Power = sprintf("%.2f", mean(tp, na.rm = TRUE)),
      Power_calib = sprintf("%.2f", mean(tp_calibrated, na.rm = TRUE)),
      N = n(),
      N_bias_ge1 = sum(abs(bias) >= 1, na.rm = TRUE),
      .groups = "drop"
    )

  # Add Type I error from null files if available
  if (nrow(all_null_FP) > 0) {
    summary_table <- summary_table %>%
      left_join(all_null_FP, by = c("Configuration", "Method")) %>%
      rename(`Type I error` = Type_I_error)
  }

  # Remove Configuration column (keep only Configuration_Label)
  summary_table <- summary_table %>% dplyr::select(-Configuration)

  # Add separator rows between configurations and create table
  summary_table_sep <- add_configuration_separators(summary_table, "Configuration_Label")
  table_grob <- pipeline_table(summary_table_sep, rows = NULL)

  res <- res %>%
    filter(!is.na(bias)) %>%
    pivot_longer(cols = c(bias, sd, rmse, w1distance, w2distance),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))

  facet_plot <- ggplot(res, aes(x = Method, y = value, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Treatment Effect Estimation Performance Across Configurations",
      x = NULL,
      y = NULL
    ) +
    facet_grid(Configuration_Label ~ metric, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray50"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  # Variance estimation plots (only if we have sigma data)
  if (nrow(res_sigma) > 0) {
    res_sigma_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
      pivot_longer(cols = c(bias.trt.sigma, sd.trt.sigma, w2distance.trt.sigma,
                            bias.ctrl.sigma, sd.ctrl.sigma, w2distance.ctrl.sigma),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.sigma", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "sd" ~ "SD",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        # Order group factor with Treatment first (on top)
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_sigma <- res_sigma_long %>%
      group_by(Configuration_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      # Reorder columns: Bias, SD, W2Distance
      dplyr::select(Configuration_Label, Method, group, Bias, SD, W2Distance)

    # Split into treatment and control tables with clear labels
    summary_table_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_trt_sep <- add_configuration_separators(summary_table_trt, "Configuration_Label")
    summary_table_ctrl_sep <- add_configuration_separators(summary_table_ctrl, "Configuration_Label")
    table_grob_trt <- pipeline_table(summary_table_trt_sep, rows = NULL)
    table_grob_ctrl <- pipeline_table(summary_table_ctrl_sep, rows = NULL)

    # Combine treatment and control tables horizontally
    table_grob_sigma <- wrap_plots(table_grob_trt, table_grob_ctrl, ncol = 2)

    facet_plot_sigma <- ggplot(res_sigma_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Variance Estimation Performance Across Configurations",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + Configuration_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )

    # Check if subject-wise columns exist
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
      # Subject-wise bias and RMSE plots
      res_subj_long <- res_sigma %>%
        filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
        pivot_longer(cols = c(bias_subj, pehe_subj),
                     names_to = "metric",
                     values_to = "value") %>%
        mutate(
          metric_type = case_when(
            metric == "bias_subj" ~ "Bias",
            metric == "pehe_subj" ~ "PEHE",
            TRUE ~ metric
          )
        )

      if (nrow(res_subj_long) > 0) {
        summary_table_subj <- res_subj_long %>%
          group_by(Configuration_Label, Method, metric_type) %>%
          summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
          pivot_wider(names_from = metric_type, values_from = Mean)

        # Add separator rows between configurations and create table
        summary_table_subj_sep <- add_configuration_separators(summary_table_subj, "Configuration_Label")
        table_grob_subj <- pipeline_table(summary_table_subj_sep, rows = NULL)

        facet_plot_subj <- ggplot(res_subj_long, aes(x = Method, y = value, fill = Method)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
          geom_boxplot(alpha = 0.7) +
          stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
          scale_fill_manual(values = method_colors) +
          labs(
            title = "Subject-wise Performance Across Configurations",
            x = NULL,
            y = NULL
          ) +
          facet_grid(Configuration_Label ~ metric_type, scales = "free_y") +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray90", color = "gray50"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
          )
      } else {
        facet_plot_subj <- NULL
        table_grob_subj <- NULL
      }
    } else {
      facet_plot_subj <- NULL
      table_grob_subj <- NULL
    }

    # Population mean plots
    res_pop_mean_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, map_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
      pivot_longer(cols = c(bias.trt.median.pop, w2distance.trt.median.pop,
                            bias.ctrl.median.pop, w2distance.ctrl.median.pop),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.median\\.pop", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_pop_mean <- res_pop_mean_long %>%
      group_by(Configuration_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean)

    # Split into treatment and control tables
    summary_table_pop_mean_trt <- summary_table_pop_mean %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_mean_ctrl <- summary_table_pop_mean %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_pop_mean_trt_sep <- add_configuration_separators(summary_table_pop_mean_trt, "Configuration_Label")
    summary_table_pop_mean_ctrl_sep <- add_configuration_separators(summary_table_pop_mean_ctrl, "Configuration_Label")
    table_grob_pop_mean_trt <- pipeline_table(summary_table_pop_mean_trt_sep, rows = NULL)
    table_grob_pop_mean_ctrl <- pipeline_table(summary_table_pop_mean_ctrl_sep, rows = NULL)
    table_grob_pop_mean <- wrap_plots(table_grob_pop_mean_trt, table_grob_pop_mean_ctrl, ncol = 2)

    facet_plot_pop_mean <- ggplot(res_pop_mean_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Mean Performance",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + Configuration_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )

    # LRC-BART MODIFICATION START
    # Size every plot and table panel from its actual facet/row count. Fixed
    # table weights caused long Sc4/Sc5 tables to spill into adjacent figures.
    n_configuration_rows <- length(unique(res$Configuration_Label))
    ate_plot_height <- max(6, 3 * n_configuration_rows)
    ate_table_height <- max(3, 0.30 * (nrow(summary_table_sep) + 1))
    pop_plot_height <- max(7, 5 * n_configuration_rows)
    pop_table_height <- max(
      3,
      0.30 * (max(nrow(summary_table_pop_mean_trt_sep),
                  nrow(summary_table_pop_mean_ctrl_sep)) + 1)
    )
    sigma_plot_height <- max(7, 5 * n_configuration_rows)
    sigma_table_height <- max(
      3,
      0.30 * (max(nrow(summary_table_trt_sep),
                  nrow(summary_table_ctrl_sep)) + 1)
    )

    # Combine ATE, population mean, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      subj_plot_height <- max(6, 3 * n_configuration_rows)
      subj_table_height <- max(
        3, 0.30 * (nrow(summary_table_subj_sep) + 1)
      )
      final_heights <- c(
        ate_plot_height, ate_table_height,
        pop_plot_height, pop_table_height,
        subj_plot_height, subj_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    } else {
      final_heights <- c(
        ate_plot_height, ate_table_height,
        pop_plot_height, pop_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    }
    # LRC-BART MODIFICATION END
  } else {
    # Only ATE plot if no sigma data
    n_configuration_rows <- length(unique(res$Configuration_Label))
    final_heights <- c(
      max(6, 3 * n_configuration_rows),
      max(3, 0.30 * (nrow(summary_table_sep) + 1))
    )
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = final_heights)
  }

  # Dynamically calculate dimensions based on the number of configurations
  # and method series that are actually reported.
  n_reporting_methods <- length(unique(as.character(
    res$Method[!is.na(res$Method)]
  )))
  n_reporting_methods <- max(1, n_reporting_methods)

  # Calculate total dimensions
  # LRC-BART ADDITION START
  # More reported methods require wider boxplot panels. Table heights above
  # are based on their actual method/configuration row counts.
  total_height <- sum(final_heights) + 2
  total_width <- max(14, 1.4 * n_reporting_methods)
  # LRC-BART ADDITION END

  save_pipeline_plot(paste0(file.path(mainDir),"/lrcbart-sim-gaussian/inserts/p",p_obs,"_sc",sc,"_all_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}
if (sc == 3){

  # Check if res has data and required columns
  if (nrow(res) == 0) {
    stop("No data remaining for sc == 3 after filtering. All iterations may have |bias| >= 1. Check the filtering threshold or data quality.")
  }
  if (!"Configuration_Label" %in% colnames(res)) {
    stop("Configuration_Label column not found in res. Check data loading.")
  }

  summary_table <- res %>%
    filter(!is.na(bias)) %>%
    group_by(Configuration_Label, Configuration, c, Method) %>%
    summarise(
      Bias = sprintf("%.2f", mean(bias, na.rm = TRUE)),
      SD = sprintf("%.2f", mean(sd, na.rm = TRUE)),
      RMSE = sprintf("%.2f", sqrt(mean(rmse, na.rm = TRUE))),
      W1Distance = sprintf("%.2f", mean(w1distance, na.rm = TRUE)),
      W2Distance = sprintf("%.2f", mean(w2distance, na.rm = TRUE)),
      CI_length = sprintf("%.2f", mean(ci, na.rm = TRUE)),
      CI_coverage = sprintf("%.2f", mean(coverage, na.rm = TRUE)),
      Power = sprintf("%.2f", mean(tp, na.rm = TRUE)),
      Power_calib = sprintf("%.2f", mean(tp_calibrated, na.rm = TRUE)),
      N = n(),
      N_bias_ge1 = sum(abs(bias) >= 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Order by configuration first, then correlation
    arrange(Configuration_Label, c, Method)

  # Add Type I error from null files if available
  if (nrow(all_null_FP) > 0) {
    summary_table <- summary_table %>%
      left_join(all_null_FP, by = c("Configuration", "c", "Method")) %>%
      rename(`Type I error` = Type_I_error)
  }

  # Create combined Configuration + Correlation label for separator function
  summary_table$Config_Cor_Label <- paste0(summary_table$Configuration_Label, ", rho=", summary_table$c)

  # Remove Configuration column (keep only Configuration_Label) and rename c to Correlation
  summary_table <- summary_table %>% dplyr::select(-Configuration) %>% rename(Correlation = c)

  # Add separator rows between each configuration-correlation combination and create table
  summary_table_sep <- add_configuration_separators(summary_table, "Config_Cor_Label")
  # Remove the Config_Cor_Label column from display (it's redundant with Configuration_Label + Correlation)
  summary_table_sep <- summary_table_sep %>% dplyr::select(-Config_Cor_Label)
  table_grob <- pipeline_table(summary_table_sep, rows = NULL)

  # Create combined Configuration + Correlation label for faceting
  # Order: group by configuration first, then by correlation within each configuration
  if (nrow(res) > 0) {
    res$Config_Cor_Label <- plot_factor(paste0(res$Configuration_Label, ", rho=", res$c),
                                   levels = as.vector(t(outer(c("Default"),
                                                              c(-0.5, 0, 0.5),
                                                              function(s, r) paste0(s, ", rho=", r)))))
  }
  if (nrow(res_sigma) > 0) {
    res_sigma$Config_Cor_Label <- plot_factor(paste0(res_sigma$Configuration_Label, ", rho=", res_sigma$c),
                                         levels = as.vector(t(outer(c("Default"),
                                                                    c(-0.5, 0, 0.5),
                                                                    function(s, r) paste0(s, ", rho=", r)))))
  }

  res <- res %>%
    filter(!is.na(bias)) %>%
    pivot_longer(cols = c(bias, sd, rmse, w1distance, w2distance),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))

  # Determine number of unique combinations that exist
  n_configuration_cor <- length(unique(res$Config_Cor_Label))

  facet_plot <- ggplot(res, aes(x = Method, y = value, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_boxplot(alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "Treatment Effect Estimation Performance Across Configurations and Correlations",
      x = NULL,
      y = NULL
    ) +
    facet_grid(Config_Cor_Label ~ metric, scales = "free_y") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13, face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = "gray50"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  # Variance estimation plots (only if we have sigma data)
  if (nrow(res_sigma) > 0) {
    res_sigma_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
      pivot_longer(cols = c(bias.trt.sigma, sd.trt.sigma, w2distance.trt.sigma,
                            bias.ctrl.sigma, sd.ctrl.sigma, w2distance.ctrl.sigma),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.sigma", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "sd" ~ "SD",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_sigma <- res_sigma_long %>%
      group_by(Config_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      dplyr::select(Config_Cor_Label, c, Method, group, Bias, SD, W2Distance) %>%
      rename(Correlation = c) %>%
      # Order by configuration first, then correlation
      arrange(Config_Cor_Label, Method)

    summary_table_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_trt_sep <- add_configuration_separators(summary_table_trt, "Config_Cor_Label")
    summary_table_ctrl_sep <- add_configuration_separators(summary_table_ctrl, "Config_Cor_Label")
    table_grob_trt <- pipeline_table(summary_table_trt_sep, rows = NULL)
    table_grob_ctrl <- pipeline_table(summary_table_ctrl_sep, rows = NULL)
    table_grob_sigma <- wrap_plots(table_grob_trt, table_grob_ctrl, ncol = 2)

    facet_plot_sigma <- ggplot(res_sigma_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Variance Estimation Performance Across Configurations and Correlations",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + Config_Cor_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )

    # Subject-wise plots
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
      res_subj_long <- res_sigma %>%
        filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
        pivot_longer(cols = c(bias_subj, pehe_subj),
                     names_to = "metric",
                     values_to = "value") %>%
        mutate(
          metric_type = case_when(
            metric == "bias_subj" ~ "Bias",
            metric == "pehe_subj" ~ "PEHE",
            TRUE ~ metric
          )
        )

      if (nrow(res_subj_long) > 0) {
        summary_table_subj <- res_subj_long %>%
          group_by(Config_Cor_Label, c, Method, metric_type) %>%
          summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
          pivot_wider(names_from = metric_type, values_from = Mean) %>%
          rename(Correlation = c) %>%
          # Order by configuration first, then correlation
          arrange(Config_Cor_Label, Method)

        # Add separator rows between configurations and create table
        summary_table_subj_sep <- add_configuration_separators(summary_table_subj, "Config_Cor_Label")
        table_grob_subj <- pipeline_table(summary_table_subj_sep, rows = NULL)

        facet_plot_subj <- ggplot(res_subj_long, aes(x = Method, y = value, fill = Method)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
          geom_boxplot(alpha = 0.7) +
          stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
          scale_fill_manual(values = method_colors) +
          labs(
            title = "Subject-wise Performance Across Configurations and Correlations",
            x = NULL,
            y = NULL
          ) +
          facet_grid(Config_Cor_Label ~ metric_type, scales = "free_y") +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "none",
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray90", color = "gray50"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
          )
      } else {
        facet_plot_subj <- NULL
        table_grob_subj <- NULL
      }
    } else {
      facet_plot_subj <- NULL
      table_grob_subj <- NULL
    }

    # Population mean plots
    res_pop_mean_long <- res_sigma %>%
      filter(Method %in% c("LMv1", "LMv2", "LMv3", lmv4_methods, map_methods, "BARTv1", "BARTv2", lrcbart_methods)) %>%
      pivot_longer(cols = c(bias.trt.median.pop, w2distance.trt.median.pop,
                            bias.ctrl.median.pop, w2distance.ctrl.median.pop),
                   names_to = "metric",
                   values_to = "value") %>%
      mutate(
        group = ifelse(grepl("trt", metric), "Treatment", "Control"),
        metric_type = gsub("\\.(trt|ctrl)\\.median\\.pop", "", metric),
        metric_type = case_when(
          metric_type == "bias" ~ "Bias",
          metric_type == "w2distance" ~ "W2Distance",
          TRUE ~ metric_type
        ),
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_pop_mean <- res_pop_mean_long %>%
      group_by(Config_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      rename(Correlation = c) %>%
      # Order by configuration first, then correlation
      arrange(Config_Cor_Label, Method)

    summary_table_pop_mean_trt <- summary_table_pop_mean %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_mean_ctrl <- summary_table_pop_mean %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_pop_mean_trt_sep <- add_configuration_separators(summary_table_pop_mean_trt, "Config_Cor_Label")
    summary_table_pop_mean_ctrl_sep <- add_configuration_separators(summary_table_pop_mean_ctrl, "Config_Cor_Label")
    table_grob_pop_mean_trt <- pipeline_table(summary_table_pop_mean_trt_sep, rows = NULL)
    table_grob_pop_mean_ctrl <- pipeline_table(summary_table_pop_mean_ctrl_sep, rows = NULL)
    table_grob_pop_mean <- wrap_plots(table_grob_pop_mean_trt, table_grob_pop_mean_ctrl, ncol = 2)

    facet_plot_pop_mean <- ggplot(res_pop_mean_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Mean Performance Across Configurations and Correlations",
        x = NULL,
        y = NULL
      ) +
      facet_grid(group + Config_Cor_Label ~ metric_type, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 13, face = "bold"),
        legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "gray90", color = "gray50"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )

    # LRC-BART MODIFICATION START
    # Size every plot and table panel from its actual facet/row count. Fixed
    # table weights caused the multi-correlation tables to overlap figures.
    n_configuration_rows <- length(unique(res$Config_Cor_Label))
    ate_plot_height <- max(6, 3 * n_configuration_rows)
    ate_table_height <- max(3, 0.30 * (nrow(summary_table_sep) + 1))
    pop_plot_height <- max(7, 5 * n_configuration_rows)
    pop_table_height <- max(
      3,
      0.30 * (max(nrow(summary_table_pop_mean_trt_sep),
                  nrow(summary_table_pop_mean_ctrl_sep)) + 1)
    )
    sigma_plot_height <- max(7, 5 * n_configuration_rows)
    sigma_table_height <- max(
      3,
      0.30 * (max(nrow(summary_table_trt_sep),
                  nrow(summary_table_ctrl_sep)) + 1)
    )

    # Combine ATE, population mean, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      subj_plot_height <- max(6, 3 * n_configuration_rows)
      subj_table_height <- max(
        3, 0.30 * (nrow(summary_table_subj_sep) + 1)
      )
      final_heights <- c(
        ate_plot_height, ate_table_height,
        pop_plot_height, pop_table_height,
        subj_plot_height, subj_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    } else {
      final_heights <- c(
        ate_plot_height, ate_table_height,
        pop_plot_height, pop_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_mean / table_grob_pop_mean / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    }
    # LRC-BART MODIFICATION END
  } else {
    # Only ATE plot if no sigma data
    n_configuration_rows <- length(unique(res$Config_Cor_Label))
    final_heights <- c(
      max(6, 3 * n_configuration_rows),
      max(3, 0.30 * (nrow(summary_table_sep) + 1))
    )
    final_plot <- facet_plot / table_grob +
      plot_layout(heights = final_heights)
  }

  # Dynamically calculate dimensions based on the number of
  # configuration/correlation rows and method series actually reported.
  n_reporting_methods <- length(unique(as.character(
    res$Method[!is.na(res$Method)]
  )))
  n_reporting_methods <- max(1, n_reporting_methods)

  # Calculate total dimensions. Table panels already scale from the number of
  # reporting rows, so use their combined requested height directly.
  total_height <- sum(final_heights) + 2
  total_width <- max(14, 1.4 * n_reporting_methods)

  save_pipeline_plot(paste0(file.path(mainDir),"lrcbart-sim-gaussian/inserts/p",p_obs,"_sc",sc,"_all_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}

}, error = function(e) {
  message(sprintf("Skipped sc = %d: %s", sc, conditionMessage(e)))
})
}

finish_pipeline_plot()
