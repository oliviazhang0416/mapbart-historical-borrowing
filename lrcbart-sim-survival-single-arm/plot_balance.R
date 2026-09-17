rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(survival)
library(survminer)
library(cowplot)
library(patchwork)
library(ggh4x)
library(ggExtra)
# LRC-BART MODIFICATION START
# Retain the inherited survival balance plots for the migrated Sc1/Sc2 data.
mainDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing"
projectDir <- file.path(mainDir, "lrcbart-sim-survival-single-arm")
n_replicates <- 100L

data_folder <- "data"
n_T <- 30L
stopifnot(n_T %in% c(30L, 200L))
size_suffix <- paste0("_n", n_T)
p_obs <- 10L
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

# LRC-BART MODIFICATION END

# LRC-BART MODIFICATION START
# Load the final requested replicate for each approved scenario. Each saved
# object already contains the treated trial arm and external controls.
all_data <- data.frame()
for (sc in selected_sc(1:2)) {
  for (hypo in plot_hypothesis) {
    filename <- file.path(
      projectDir, data_folder,
      paste0("data_p", p_obs, size_suffix, "_sc", sc, "_", hypo, "_",
             n_replicates, ".RData")
    )
    data_full <- tryCatch(read_plot_data(filename), error = function(e) NULL)
    if (is.null(data_full)) {
      message("Data file not found or unreadable: ", filename)
      next
    }

    data <- cbind(
      data_full$X,
      data.frame(y = data_full$y, delta = data_full$delta)
    )
    data <- data[
      (data$D == 1 & data$Z == 1) | (data$D == 0 & data$Z == 0),
      , drop = FALSE
    ]
    data$U1 <- NA_real_
    data$U2 <- NA_real_
    data$Scenario <- sc
    data$Configuration <- "Default"
    data$Correlation <- NA_real_
    data$Hypothesis <- hypo
    all_data <- rbind(all_data, data)
  }
}
# LRC-BART MODIFICATION END

if (nrow(all_data) == 0) {
  stop("No data files were found. Please check that the data files exist for the specified scenarios.")
}

# Check which hypotheses were loaded
hypotheses_loaded <- unique(all_data$Hypothesis)
cat("Successfully loaded data for hypothesis:", paste(hypotheses_loaded, collapse = ", "), "\n")

# Define new color palette for groups
group_colors <- c("External Control" = "lightblue", "Treatment" = "red")

# Recreate Group factor (single-arm: no RCT control)
all_data$Group <- factor(ifelse(all_data$Z == 1, "Treatment", "External Control"),
                     levels = c("External Control", "Treatment"))

# Create factor for Hypothesis
all_data$Hypothesis <- factor(all_data$Hypothesis, levels = c("null", "alternative"))

# Create combined Group-Hypothesis variable for stratification
all_data$Group_Hypo <- interaction(all_data$Group, all_data$Hypothesis, sep = " - ")

# LRC-BART MODIFICATION START
# The migrated design has one configuration and no configuration variable.
all_data$Scenario_Label <- factor(
  paste0("Scenario ", all_data$Scenario),
  levels = c("Scenario 1", "Scenario 2")
)

# Reorder data so External Control is plotted last (on top).
all_data_reordered <- all_data %>% arrange(desc(Group == "External Control"))
all_data_reordered$Scenario_F <- factor(
  paste0("Scenario ", all_data_reordered$Scenario),
  levels = c("Scenario 1", "Scenario 2")
)
all_data_reordered$Configuration_F <- plot_factor(
  all_data_reordered$Configuration, levels = "Default"
)
# LRC-BART MODIFICATION END

# Create Correlation label for column faceting
all_data_reordered$Cor_Label <- ifelse(is.na(all_data_reordered$Correlation),
                                        "",
                                        paste0("ρ=", all_data_reordered$Correlation))

# Kaplan-Meier survival curves for outcome
# Create separate lists for null and alternative plots
km_plots_null <- list()
km_plots_alt <- list()

for (scenario in levels(all_data$Scenario_Label)) {
  df <- all_data[all_data$Scenario_Label == scenario, ]

  # Get available hypotheses for this scenario
  hypos_in_scenario <- unique(df$Hypothesis)
  hypos_in_scenario <- hypos_in_scenario[!is.na(hypos_in_scenario)]

  # Create null plot if available
  if ("null" %in% hypos_in_scenario) {
    df_null <- df[df$Hypothesis == "null", ]

    fit_null <- survfit(Surv(y, delta) ~ Group, data = df_null)

    p_null <- ggsurvplot(fit_null, data = df_null,
                         palette = c("lightblue", "red"),
                         legend = "none",
                         xlab = "Time",
                         ylab = "Survival Probability",
                         title = paste0(scenario, "\n(H0: Null)"),
                         ggtheme = theme_classic(base_size = 9),
                         conf.int = FALSE,
                         risk.table = FALSE)$plot +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    km_plots_null[[length(km_plots_null) + 1]] <- p_null
  }

  # Create alternative plot if available
  if ("alternative" %in% hypos_in_scenario) {
    df_alt <- df[df$Hypothesis == "alternative", ]

    fit_alt <- survfit(Surv(y, delta) ~ Group, data = df_alt)

    p_alt <- ggsurvplot(fit_alt, data = df_alt,
                        palette = c("lightblue", "red"),
                        legend = "none",
                        xlab = "Time",
                        ylab = "Survival Probability",
                        title = paste0(scenario, "\n(H1: Alternative)"),
                        ggtheme = theme_classic(base_size = 9),
                        conf.int = FALSE,
                        risk.table = FALSE)$plot +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    km_plots_alt[[length(km_plots_alt) + 1]] <- p_alt
  }
}

# Dynamically determine layout for KM plots
n_scenarios <- length(levels(all_data$Scenario_Label))

# Determine number of rows needed
if (length(km_plots_null) > 0 && length(km_plots_alt) > 0) {
  km_nrow <- 2
} else if (length(km_plots_null) > 0 || length(km_plots_alt) > 0) {
  km_nrow <- 1
} else {
  stop("No KM plots were created")
}

# Calculate km_ncol based on actual number of plots (not just n_scenarios)
# Use the maximum of null and alt plot counts to ensure we have enough columns
n_km_plots <- max(length(km_plots_null), length(km_plots_alt))
km_ncol <- n_km_plots  # Use actual number of plots for columns

# Use patchwork to combine KM plots instead of ggarrange
# This is compatible with the patchwork layout used later
if (length(km_plots_null) > 0 && length(km_plots_alt) > 0) {
  # Both null and alternative: create two rows
  # Each row gets the full number of columns
  p_km <- wrap_plots(km_plots_null, ncol = length(km_plots_null)) /
          wrap_plots(km_plots_alt, ncol = length(km_plots_alt))
} else if (length(km_plots_null) > 0) {
  # Only null plots
  p_km <- wrap_plots(km_plots_null, ncol = length(km_plots_null))
} else if (length(km_plots_alt) > 0) {
  # Only alternative plots
  p_km <- wrap_plots(km_plots_alt, ncol = length(km_plots_alt))
} else {
  stop("No KM plots were created")
}

# Density plots of observed survival times (censored vs uncensored)
# Commented out: these plots are not included in the final figure
# all_data_reordered$Event_Status <- factor(ifelse(all_data_reordered$delta == 1, "Event", "Censored"),
#                                            levels = c("Event", "Censored"))
#
# # Create combined facet variable for Scenario and Event Status
# all_data_reordered$Facet_Label <- paste0(all_data_reordered$Scenario_Label, "\n", all_data_reordered$Event_Status)
#
# p_density <- ggplot(all_data_reordered, aes(x = y, fill = Group)) +
#   geom_density(alpha = 0.6) +
#   facet_wrap(~Facet_Label, nrow = n_scenarios, ncol = 2, scales = "free_y") +
#   labs(title = "Distribution of Observed Times by Event Status", x = "Time", y = "Density") +
#   scale_fill_manual(values = group_colors, name = "Group") +
#   theme_classic(base_size = 9) +
#   theme(legend.position = "none",
#         strip.text = element_text(size = 7, face = "bold"),
#         strip.background = element_rect(fill = "grey90", color = "black"),
#         panel.spacing = unit(0.3, "lines"))

# Create a dummy plot with all three groups for legend
dummy_data <- data.frame(
  x = rep(1, 2),
  Group = factor(c("External Control", "Treatment"),
                 levels = c("External Control", "Treatment"))
)

legend_plot <- ggplot(dummy_data, aes(x = x, fill = Group)) +
  geom_bar() +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Extract legend
legend <- get_legend(legend_plot)

# Reshape data to long format for combined covariate plot
# For sc == 1 and 2, only X5, X6 are available (confounders)
# For sc == 3, X5, U1, X6, U2 are available

library(tidyr)

# Create long format data for covariates
# Handle sc == 1, sc == 2, and sc == 3 separately
sc1_data <- all_data_reordered[all_data_reordered$Scenario == 1, ]
sc2_data <- all_data_reordered[all_data_reordered$Scenario == 2, ]
sc3_data <- all_data_reordered[all_data_reordered$Scenario == 3, ]

# For sc == 1: reshape X5, X6 and include outcome y
if (nrow(sc1_data) > 0) {
  sc1_long <- sc1_data %>%
    dplyr::select(Scenario_F, Configuration_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc1_long$Covariate <- factor(sc1_long$Covariate, levels = c("X5", "X6"))
} else {
  sc1_long <- NULL
}

# For sc == 2: reshape X5, X6 and include outcome y
if (nrow(sc2_data) > 0) {
  sc2_long <- sc2_data %>%
    dplyr::select(Scenario_F, Configuration_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc2_long$Covariate <- factor(sc2_long$Covariate, levels = c("X5", "X6"))
} else {
  sc2_long <- NULL
}

# For sc == 3: reshape X5, U1, X6, U2 in the specified order and include outcome y
if (nrow(sc3_data) > 0 && "U1" %in% colnames(sc3_data)) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_F, Configuration_F, Cor_Label, Group, X5, U1, X6, U2, y) %>%
    pivot_longer(cols = c(X5, U1, X6, U2),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(sc3_long$Covariate, levels = c("X5", "U1", "X6", "U2"))
  has_U_plot <- TRUE
} else if (nrow(sc3_data) > 0) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_F, Configuration_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(sc3_long$Covariate, levels = c("X5", "X6"))
  has_U_plot <- FALSE
} else {
  sc3_long <- NULL
  has_U_plot <- FALSE
}

# Load ggh4x for nested faceting with visual separators
library(ggh4x)

# Create covariate vs outcome plot for sc == 1 with marginal distributions
if (!is.null(sc1_long) && nrow(sc1_long) > 0) {
  # Filter out rows with NA in Value or y
  sc1_long_clean <- sc1_long[!is.na(sc1_long$Value) & !is.na(sc1_long$y) & sc1_long$y > 0, ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc1_main <- ggplot(sc1_long_clean, aes(x = Value, y = log(y), color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "log(Y)") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and configuration)
  p_cov_sc1_top <- ggplot(sc1_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 1", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour (3 covariate rows)
  p_cov_sc1 <- wrap_elements(p_cov_sc1_top / p_cov_sc1_main + plot_layout(heights = c(1, 1)))
  has_sc1_plot <- TRUE
} else {
  has_sc1_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 2 with marginal distributions
if (!is.null(sc2_long) && nrow(sc2_long) > 0) {
  # Filter out rows with NA in Value or y
  sc2_long_clean <- sc2_long[!is.na(sc2_long$Value) & !is.na(sc2_long$y) & sc2_long$y > 0, ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc2_main <- ggplot(sc2_long_clean, aes(x = Value, y = log(y), color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "log(Y)") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and configuration)
  p_cov_sc2_top <- ggplot(sc2_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 2", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour (3 covariate rows)
  p_cov_sc2 <- wrap_elements(p_cov_sc2_top / p_cov_sc2_main + plot_layout(heights = c(1, 1)))
  has_sc2_plot <- TRUE
} else {
  has_sc2_plot <- FALSE
}

# Combine sc1 and sc2 plots vertically (since each now has histogram + outcome density)
if (has_sc1_plot && has_sc2_plot) {
  p_cov_sc12 <- p_cov_sc1 / p_cov_sc2
  has_sc12_plot <- TRUE
} else if (has_sc1_plot) {
  p_cov_sc12 <- p_cov_sc1
  has_sc12_plot <- TRUE
} else if (has_sc2_plot) {
  p_cov_sc12 <- p_cov_sc2
  has_sc12_plot <- TRUE
} else {
  has_sc12_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 3 with marginal distributions
if (!is.null(sc3_long) && nrow(sc3_long) > 0) {
  # Filter out rows with NA in Value or y
  sc3_long_clean <- sc3_long[!is.na(sc3_long$Value) & !is.na(sc3_long$y) & sc3_long$y > 0, ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc3_main <- ggplot(sc3_long_clean, aes(x = Value, y = log(y), color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "log(Y)") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Top marginal: Covariate distribution histogram (one per covariate and configuration)
  p_cov_sc3_top <- ggplot(sc3_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenario 3", x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combine with patchwork: histogram grid above main contour (6 covariate rows)
  p_cov_sc3 <- wrap_elements(p_cov_sc3_top / p_cov_sc3_main + plot_layout(heights = c(2, 3)))
  has_sc3_plot <- TRUE
} else {
  has_sc3_plot <- FALSE
}

# Dynamically set heights based on number of scenarios
km_height <- max(5.0, km_nrow * 1.0)  # Scale KM plot height with number of rows
# density_height <- max(10.0, n_scenarios * 2.0)  # Scale density plot height with scenarios (commented out)

# Calculate covariate plot heights based on number of rows (Configuration + Covariate)
n_sc1_configuration <- length(unique(sc1_data$Configuration_F[!is.na(sc1_data$Configuration_F)]))
n_sc2_configuration <- length(unique(sc2_data$Configuration_F[!is.na(sc2_data$Configuration_F)]))
n_sc3_configuration <- length(unique(sc3_data$Configuration_F[!is.na(sc3_data$Configuration_F)]))
# sc1 and sc2 are side by side, so use max of their heights
cov_height_sc12 <- max(n_sc1_configuration, n_sc2_configuration, n_sc3_configuration) * 10
cov_height_sc3 <- max(n_sc1_configuration, n_sc2_configuration, n_sc3_configuration) * 10

# Create a plot from the legend
legend_plot_final <- ggplot(dummy_data, aes(x = x, fill = Group)) +
  geom_bar() +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_void() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Arrange the plots using patchwork
# p_km is already a patchwork object, no conversion needed
# Note: p_density has been removed from the layout (density plots commented out)
if (has_sc12_plot && has_sc3_plot) {
  plot <- legend_plot_final / p_km / p_cov_sc12 / p_cov_sc3 +
    plot_layout(heights = c(0.5, km_height, cov_height_sc12, cov_height_sc3))
  total_height <- 0.1 + km_height + cov_height_sc12 + cov_height_sc3 + 1
} else if (has_sc12_plot) {
  plot <- legend_plot_final / p_km / p_cov_sc12 +
    plot_layout(heights = c(0.5, km_height, cov_height_sc12))
  total_height <- 0.1 + km_height + cov_height_sc12 + 1
} else if (has_sc3_plot) {
  plot <- legend_plot_final / p_km / p_cov_sc3 +
    plot_layout(heights = c(0.5, km_height, cov_height_sc3))
  total_height <- 0.1 + km_height + cov_height_sc3 + 1
} else {
  plot <- legend_plot_final / p_km +
    plot_layout(heights = c(0.5, km_height))
  total_height <- 0.1 + km_height + 1
}

# Calculate total plot dimensions dynamically
total_width <- km_ncol * 3 + 2  # Scale width with number of columns, +2 for margins

save_pipeline_plot(file.path(projectDir, "inserts", paste0("balance_all_scenarios", size_suffix, ".jpg")),
       width = total_width,
       height = total_height,
       plot,
       limitsize = FALSE)

finish_pipeline_plot()
