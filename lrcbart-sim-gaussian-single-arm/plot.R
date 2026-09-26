# JOINT LRC-BART MODIFICATION START
# Matching original single-arm plot, adapted to the agreed Sc1 descendants.
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
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
# LRC-BART MODIFICATION START
# The inherited Gaussian result plot is retained for the two approved
# single-arm scenarios and their available comparator and LRC-BART results.
mainDir <- .lrcRoot
# LRC-BART MODIFICATION END

p_obs <- 10
# Pipeline selection; standalone runs use the same alternative-hypothesis default.
plot_context <- list(hypothesis = if ("--null" %in% commandArgs(TRUE)) "null" else "alternative")
# ---- Local plot selection and result readers ----
# Plot selection is independent of which comparison methods were fitted today.
`%||%` <- function(x, default) if (is.null(x)) default else x
plot_context <- plot_context %||% list()
plot_hypothesis <- plot_context$hypothesis %||% "alternative"
stopifnot(length(plot_hypothesis) == 1L, plot_hypothesis %in% c("null", "alternative"))
plot_saved <- character()

plot_scenario_id <- function(cfg) {
  variant <- cfg$variant %||% ""
  tag <- paste0("sc", cfg$sc, variant)
  if (cfg$sc == 3) tag <- paste0(tag, "_cor", cfg$cor)
  if (variant == "b") tag <- paste0(tag, "_", cfg$region)
  if (nzchar(variant)) tag <- paste0(tag, "_d", cfg$delta_rwd)
  tag
}
selected_sc <- function(default) {
  if (is.null(plot_context$scenarios)) return(default)
  chosen <- vapply(plot_context$scenarios, function(x)
    paste0(x$sc, x$variant %||% ""), character(1))
  default[default %in% chosen]
}
selected_configurations <- function(sc, default) {
  if (is.null(plot_context$scenarios) || sc %in% c("1", "2", "3")) return(default)
  chosen <- Filter(function(x) paste0(x$sc, x$variant %||% "") == sc,
                   plot_context$scenarios)
  unique(vapply(chosen, function(x) if (sc == "1b")
    paste0(x$region, "_d", x$delta_rwd) else paste0("d", x$delta_rwd), character(1)))
}
# JOINT LRC-BART MODIFICATION END
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
  if (identical(attr(result, "complete"), FALSE))
    stop("Incomplete reporting file: ", basename(path))
  message("Plot input: ", basename(path), " (", nrow(result), " replicates)")
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
lrc_plot_configurations <- function(single_arm) {
  defaults <- if (single_arm) c("default", "w0.9", "s0min") else "default"
  configs <- plot_context$lrc_configs %||% defaults
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
# Retained single-arm hierarchical-model prior scale.
lmv4_prior_vals <- plot_method_setting("hierLM.R", "prior_vals", 0.05)

# MAP and lrcBART now key on N_target rather than the raw s^2 value: each
# analysis file loops over multiple N_target multipliers from its
# calibration .RData and writes one result file per target (with the
# integer target N in the filename).  List the N values to read here.
map_target_Ns <- plot_method_setting("MAP.R", "target_Ns", numeric(0))
# Standalone default; run_all.R overrides this from the ess calibration file.
# Listed in decreasing N (Inf first) for LRC-BART curve/legend order.
lmv2_rwd_w_vals <- plot_method_setting("LMv2.R", "rwd_w_vals", 1)
lrcbart_configurations <- lrc_plot_configurations(TRUE)
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

for (sc in selected_sc(c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2"))) {
tryCatch({

if (sc %in% c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2")){

  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  configurations <- if (sc == "1a") c("d1", "d2") else
    if (sc == "1b") c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2") else
    if (sc %in% c("1c", "1c-i", "1c-ii")) "d2" else "Base"

  for (configuration in selected_configurations(sc, configurations)) {

    # Initialize empty lists for this configuration
    configuration_res_list_ATE <- list()
    configuration_res_list_sigma <- list()

    # LRC-BART MODIFICATION START
    scenario_suffix <- paste0("_sc", sc,
      if (sc %in% c("1", "2")) "" else paste0("_", configuration))
    configuration_label <- if (sc %in% c("1", "2")) "Default" else configuration
    file_suffix <- paste0(scenario_suffix, paste0("_", plot_hypothesis, ".RData"))
    # LRC-BART MODIFICATION END

    # Try to read LMv2 results
    for (.rwd_w in lmv2_rwd_w_vals) {
    tryCatch({
      LMv2_file <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/res/LMv2_p",p_obs,
                          scenario_suffix,paste0("_", plot_hypothesis, "_w", .rwd_w, ".RData"))
      if (!file.exists(LMv2_file)) stop("LMv2 result is not available")
      LMv2_res <- read_plot_result(LMv2_file)
      LMv2_res_ATE <- LMv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      LMv2_res_sigma <- LMv2_res[, c("c","iteration",
                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                  "pehe_subj","bias_subj")]
      LMv2_res_ATE$Method <- weight_labels("LMv2", lmv2_rwd_w_vals)[match(.rwd_w, lmv2_rwd_w_vals)]
      LMv2_res_sigma$Method <- weight_labels("LMv2", lmv2_rwd_w_vals)[match(.rwd_w, lmv2_rwd_w_vals)]
      LMv2_res_ATE$Configuration <- configuration_label
      LMv2_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LMv2_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LMv2_res_sigma
    }, error = function(e) {
    })

    }

    # Try to read hierLM results for each prior value
    for (prior_val in lmv4_prior_vals) {
      tryCatch({
        # Construct file path with prior.
        lmv4_file <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/res/hierLM_p",p_obs,
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

    # Try to read BARTv2 results
    for (bver in "BARTv2") {
    tryCatch({
      BART_file <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/res/",bver,
                          "_p",p_obs,file_suffix)
      if (!file.exists(BART_file)) stop("BARTv2 result is not available")
      BART_res <- read_plot_result(BART_file)
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

    # Try to read every approved single-arm LRC-BART configuration.
    for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      tryCatch({
        target_N <- lrcbart_configurations$target[lrc_index]
        lrcbart_file <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/res/LRC-BART_p",p_obs,
                               scenario_suffix,
                               lrcbart_configurations$suffix[lrc_index],
                               "_N",target_N,
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
          lrcBART_res_ATE$Method <- lrcbart_configurations$label[lrc_index]
          lrcBART_res_sigma$Method <- lrcbart_configurations$label[lrc_index]
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
    scenario_suffix <- paste0("_sc", sc,
      if (sc %in% c("1", "2")) "" else paste0("_", configuration))
    configuration_label <- if (sc %in% c("1", "2")) "Default" else configuration
    null_file_suffix <- paste0(scenario_suffix, "_null.RData")

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
        null_file <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
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
# Check if we have any data to plot
if (nrow(res) == 0) {
  stop("No results data found. Please check that result files exist and paths are correct.")
}

# Create plotting-row labels for all scenarios.
if (sc %in% c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2")) {
  configuration_levels <- "Default"
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
all_method_levels <- c(weight_labels("LMv2", lmv2_rwd_w_vals), "BARTv2", lmv4_methods, lrcbart_methods)

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

if (sc %in% c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2")){

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

  save_pipeline_plot(paste0(file.path(mainDir),"/lrcbart-sim-gaussian-single-arm/inserts/p",p_obs,"_sc",sc,"_all_results.jpg"),
         width = total_width,
         height = total_height,
         final_plot,
         limitsize = FALSE)
}
}, error = function(e) {
  message(sprintf("Skipped sc = %s: %s", sc, conditionMessage(e)))
})
}

finish_pipeline_plot()

# JOINT LRC-BART MODIFICATION END
