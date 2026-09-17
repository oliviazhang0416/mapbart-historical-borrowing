rm(list=ls())
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggh4x)
library(ggExtra)
library(tidyr)
# LRC-BART MODIFICATION START
# The inherited balance plot is retained for the two approved single-arm
# scenarios. Its only displayed groups are RCT treatment and RWD control.
mainDir <- "/Users/oliviazhang/Desktop/lrcbart-historical-borrowing/"
data_folder <- "data"
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

iter_rct <- 50
iter_rwd <- 50

# Load data for all scenarios and configurations
# Try to load both null and alternative data
all_data <- data.frame()

for (sc in selected_sc(1:2)) {
  if (sc %in% c(1, 2, 4, 5)){
    # LRC-BART ADDITION START
    # Non-correlation scenarios use their scenario-specific shift/region
    # configurations as plotting rows.
    if (sc %in% c(1, 2)) configurations <- c("Base")
    if (sc == 4) configurations <- c("d1", "d2")
    if (sc == 5)
      configurations <- c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2")
    # LRC-BART ADDITION END
    cor <- 1
  } else {
    # For scenario 3, configurations with correlation
    configurations <- c("Base")
    cor <- c(-0.5, 0, 0.5)
  }

  if (sc %in% c(1, 2, 4, 5)) {
    for (configuration in selected_configurations(sc, configurations)) {
      for (c in selected_correlations(sc, cor)) {
        # Try both null and alternative
        for (hypo in plot_hypothesis) {
          # LRC-BART ADDITION START
          # Construct the exact Sc1/Sc2/Sc4/Sc5 tag written by data_gen.
          if (sc %in% c(1, 2)) {
            scenario_tag <- paste0("sc", sc)
            configuration_label <- "Base"
          }
          if (sc == 4) {
            scenario_tag <- paste0("sc4_", configuration)
            configuration_label <- paste0("delta RWD = ", sub("^d", "", configuration))
          }
          if (sc == 5) {
            scenario_tag <- paste0("sc5_", configuration)
            configuration_parts <- strsplit(configuration, "_d", fixed = TRUE)[[1]]
            configuration_label <- paste0("region = ", configuration_parts[1],
                                  ", delta RWD = ", configuration_parts[2])
          }
          filename_rct <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,
                                 "/data_p",p_obs,"_",scenario_tag,"_",hypo,
                                 "_",iter_rct,".RData")
          filename_rwd <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,
                                 "/data_p",p_obs,"_",scenario_tag,"_",hypo,
                                 "_",iter_rwd,".RData")
          # LRC-BART ADDITION END

          # Read RCT data and extract RCT portion (D==1)
          data_rct_full <- if (file.exists(filename_rct)) tryCatch({
            read_plot_data(filename_rct)
          }, error = function(e) {
            message(paste0("Error reading RCT file: ", filename_rct))
            return(NULL)
          }) else NULL

          # Read RWD data and extract RWD portion (D==0)
          data_rwd_full <- if (file.exists(filename_rwd)) tryCatch({
            read_plot_data(filename_rwd)
          }, error = function(e) {
            message(paste0("Error reading RWD file: ", filename_rwd))
            return(NULL)
          }) else NULL

          # Combine RCT and RWD portions
          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            # Extract RCT portion (D==1)
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            # Add NA columns for U (sc == 1, 2 don't have U)
            data_rct$U1 <- NA
            data_rct$U2 <- NA
            data_rct <- data_rct[data_rct$D == 1, ]

            # Extract RWD portion (D==0)
            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            # Add NA columns for U (sc == 1, 2 don't have U)
            data_rwd$U1 <- NA
            data_rwd$U2 <- NA
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            # Combine RCT and RWD
            data <- rbind(data_rct, data_rwd)

            # Add scenario and configuration labels
            data$Scenario <- sc
            data$Configuration <- configuration_label
            data$Correlation <- NA
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  } else {
    # Scenario 3: loop through configurations and correlations
    for (configuration in selected_configurations(sc, configurations)) {
      for (c in selected_correlations(sc, cor)) {
        # Try both null and alternative
        for (hypo in plot_hypothesis) {
          # Construct filename: data_p10_sc3F_cor0.5_alternative_1.RData
          if (configuration == "Base") {
            filename_rct <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
            configuration_label <- configuration
          } else {
            filename_rct <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,configuration,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/lrcbart-sim-gaussian-single-arm/",data_folder,"/data_p",p_obs,"_sc",sc,configuration,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
            configuration_label <- configuration
          }

          # Read RCT data and extract RCT portion (D==1)
          data_rct_full <- if (file.exists(filename_rct)) tryCatch({
            read_plot_data(filename_rct)
          }, error = function(e) {
            message(paste0("Error reading RCT file: ", filename_rct))
            return(NULL)
          }) else NULL

          # Read RWD data and extract RWD portion (D==0)
          data_rwd_full <- if (file.exists(filename_rwd)) tryCatch({
            read_plot_data(filename_rwd)
          }, error = function(e) {
            message(paste0("Error reading RWD file: ", filename_rwd))
            return(NULL)
          }) else NULL

          # Combine RCT and RWD portions
          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            # Extract RCT portion (D==1)
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y))
            # Add U columns if available (sc == 3), otherwise add NA columns
            if (!is.null(data_rct_full$U)) {
              data_rct <- cbind(data_rct, data_rct_full$U)
            } else {
              data_rct$U1 <- NA
              data_rct$U2 <- NA
            }
            data_rct <- data_rct[data_rct$D == 1, ]

            # Extract RWD portion (D==0)
            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y))
            # Add U columns if available (sc == 3), otherwise add NA columns
            if (!is.null(data_rwd_full$U)) {
              data_rwd <- cbind(data_rwd, data_rwd_full$U)
            } else {
              data_rwd$U1 <- NA
              data_rwd$U2 <- NA
            }
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            # Combine RCT and RWD
            data <- rbind(data_rct, data_rwd)

            # Add scenario label, configuration, and correlation
            data$Scenario <- sc
            data$Configuration <- configuration_label
            data$Correlation <- c
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  }
}

# Check if any data was loaded
if (nrow(all_data) == 0) {
  stop("No data files were found. Please check that the data files exist for the specified scenarios.")
}

# Check which hypotheses were loaded
hypotheses_loaded <- unique(all_data$Hypothesis)
cat("Successfully loaded data for hypothesis:", paste(hypotheses_loaded, collapse = ", "), "\n")

# Define new color palette for groups
group_colors <- c("External Control" = "lightblue", "Treatment" = "red")

# Recreate Group factor
all_data$Group <- factor(ifelse(all_data$Z == 1, "Treatment",
                                "External Control"),
                         levels = c("External Control", "Treatment"))

# Create factor for Hypothesis
all_data$Hypothesis <- factor(all_data$Hypothesis, levels = c("null", "alternative"))

# Create combined Group-Hypothesis variable for stratification
all_data$Group_Hypo <- interaction(all_data$Group, all_data$Hypothesis, sep = " - ")

# Create combined scenario label for faceting
all_data$Scenario_Label <- NA

# Get unique combinations that actually exist in the data
existing_noncor <- unique(all_data[all_data$Scenario %in% c(1, 2, 4, 5),
                                   c("Scenario", "Configuration")])
existing_sc3 <- unique(all_data[all_data$Scenario == 3, c("Scenario", "Configuration", "Correlation")])

# LRC-BART ADDITION START
# For non-correlation scenarios, use Base for Sc1/Sc2 and show the actual
# shift/region configuration for Sc4/Sc5.
if (nrow(existing_noncor) > 0) {
  for (i in 1:nrow(existing_noncor)) {
    sc <- existing_noncor$Scenario[i]
    configuration <- existing_noncor$Configuration[i]
    idx <- which(all_data$Scenario == sc & all_data$Configuration == configuration)
    all_data$Scenario_Label[idx] <- if (sc %in% c(1, 2))
      paste0("Scenario ", sc, ".", configuration) else
      paste0("Scenario ", sc, "\n", configuration)
  }
}
# LRC-BART ADDITION END

# For scenario 3: create labels like "Scenario 3.Base\nρ=0.5"
if (nrow(existing_sc3) > 0) {
  for (i in 1:nrow(existing_sc3)) {
    configuration <- existing_sc3$Configuration[i]
    c <- existing_sc3$Correlation[i]
    idx <- which(all_data$Scenario == 3 & all_data$Configuration == configuration & all_data$Correlation == c)
    all_data$Scenario_Label[idx] <- paste0("Scenario 3.", configuration, "\nρ=", c)
  }
}

# Define factor levels in order (only for scenarios that exist)
# For scenario 3, create all combinations of configurations and correlations
sc3_labels <- expand.grid(configuration = c("Base"),
                          cor = c(-0.5, 0, 0.5))
sc3_labels <- paste0("Scenario 3.", sc3_labels$configuration, "\nρ=", sc3_labels$cor)

all_scenario_labels <- c(paste0("Scenario 1.", c("Base")),
                         paste0("Scenario 2.", c("Base")),
                         sc3_labels,
                         paste0("Scenario 4\ndelta RWD = ", c(1, 2)),
                         paste0("Scenario 5\nregion = X5, delta RWD = ",
                                c(0.5, 1, 2)),
                         paste0("Scenario 5\nregion = X7, delta RWD = ",
                                c(1, 2)))
# Keep only labels that actually exist in the data
existing_labels <- all_scenario_labels[all_scenario_labels %in% unique(all_data$Scenario_Label)]
all_data$Scenario_Label <- plot_factor(all_data$Scenario_Label, levels = existing_labels)

# Reorder data so External Control is plotted last (on top)
all_data_reordered <- all_data %>% arrange(desc(Group == "External Control"))

# Create factor for Scenario to ensure proper ordering
all_data_reordered$Scenario_F <- factor(paste0("Scenario ", all_data_reordered$Scenario),
                                         levels = paste0("Scenario ", 1:5))

# Create factor for Configuration to ensure proper ordering
all_data_reordered$Configuration_F <- plot_factor(all_data_reordered$Configuration,
                                            levels = c(
                                              "Base",
                                              paste0("delta RWD = ", c(1, 2)),
                                              paste0("region = X5, delta RWD = ",
                                                     c(0.5, 1, 2)),
                                              paste0("region = X7, delta RWD = ",
                                                     c(1, 2))
                                            ))

# Create Correlation label for column faceting
all_data_reordered$Cor_Label <- ifelse(is.na(all_data_reordered$Correlation),
                                        "",
                                        paste0("ρ=", all_data_reordered$Correlation))

# Outcome distribution plots - separate by hypothesis
# Create separate lists for null and alternative plots
outcome_plots_null <- list()
outcome_plots_alt <- list()

n_scenarios <- length(levels(all_data$Scenario_Label))

for (scenario in levels(all_data$Scenario_Label)) {
  df <- all_data[all_data$Scenario_Label == scenario, ]

  # Get available hypotheses for this scenario
  hypos_in_scenario <- unique(df$Hypothesis)
  hypos_in_scenario <- hypos_in_scenario[!is.na(hypos_in_scenario)]

  # Create null plot if available
  if ("null" %in% hypos_in_scenario) {
    df_null <- df[df$Hypothesis == "null", ]

    p_null <- ggplot(df_null, aes(x = y, color = Group, fill = Group)) +
      geom_density(alpha = 0.3, linewidth = 1.2) +
      labs(title = paste0(scenario, "\n(H0: Null)"), x = "Y", y = "Density") +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    outcome_plots_null[[length(outcome_plots_null) + 1]] <- p_null
  }

  # Create alternative plot if available
  if ("alternative" %in% hypos_in_scenario) {
    df_alt <- df[df$Hypothesis == "alternative", ]

    p_alt <- ggplot(df_alt, aes(x = y, color = Group, fill = Group)) +
      geom_density(alpha = 0.3, linewidth = 1.2) +
      labs(title = paste0(scenario, "\n(H1: Alternative)"), x = "Y", y = "Density") +
      scale_fill_manual(values = group_colors, name = "Group") +
      scale_color_manual(values = group_colors, name = "Group") +
      theme_classic(base_size = 9) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9, face = "bold", hjust = 0.5))

    outcome_plots_alt[[length(outcome_plots_alt) + 1]] <- p_alt
  }
}

# Determine number of rows needed for outcome plots
if (length(outcome_plots_null) > 0 && length(outcome_plots_alt) > 0) {
  outcome_nrow <- 2
} else if (length(outcome_plots_null) > 0 || length(outcome_plots_alt) > 0) {
  outcome_nrow <- 1
} else {
  stop("No outcome plots were created")
}

# Calculate outcome_ncol based on actual number of plots
n_outcome_plots <- max(length(outcome_plots_null), length(outcome_plots_alt))
outcome_ncol <- n_outcome_plots

# Use patchwork to combine outcome plots
if (length(outcome_plots_null) > 0 && length(outcome_plots_alt) > 0) {
  p_outcome <- wrap_plots(outcome_plots_null, ncol = length(outcome_plots_null)) /
          wrap_plots(outcome_plots_alt, ncol = length(outcome_plots_alt))
} else if (length(outcome_plots_null) > 0) {
  p_outcome <- wrap_plots(outcome_plots_null, ncol = length(outcome_plots_null))
} else if (length(outcome_plots_alt) > 0) {
  p_outcome <- wrap_plots(outcome_plots_alt, ncol = length(outcome_plots_alt))
} else {
  stop("No outcome plots were created")
}

# Create a dummy plot with both single-arm groups for the shared legend.
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

# Create long format data for covariates
# Handle sc == 1, sc == 2, and sc == 3 separately
sc1_data <- all_data_reordered[all_data_reordered$Scenario == 1, ]
sc2_data <- all_data_reordered[all_data_reordered$Scenario == 2, ]
sc3_data <- all_data_reordered[all_data_reordered$Scenario == 3, ]
# LRC-BART ADDITION START
sc45_data <- all_data_reordered[all_data_reordered$Scenario %in% c(4, 5), ]
# LRC-BART ADDITION END

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

# LRC-BART ADDITION START
# Sc4 and Sc5 retain the inherited measured-covariate display while separating
# the global-shift and regional-partial-compatibility configurations by row.
if (nrow(sc45_data) > 0) {
  sc45_long <- sc45_data %>%
    dplyr::select(Scenario_F, Configuration_F, Cor_Label, Group, X5, X6, y) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc45_long$Covariate <- factor(sc45_long$Covariate,
                                levels = c("X5", "X6"))
} else {
  sc45_long <- NULL
}
# LRC-BART ADDITION END

# Create covariate vs outcome plot for sc == 1 with marginal distributions
if (!is.null(sc1_long) && nrow(sc1_long) > 0) {
  # Filter out rows with NA in Value or y
  sc1_long_clean <- sc1_long[!is.na(sc1_long$Value) & !is.na(sc1_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc1_main <- ggplot(sc1_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
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

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc1 <- wrap_elements(p_cov_sc1_top / p_cov_sc1_main + plot_layout(heights = c(1, 1)))
  has_sc1_plot <- TRUE
} else {
  has_sc1_plot <- FALSE
}

# Create covariate vs outcome plot for sc == 2 with marginal distributions
if (!is.null(sc2_long) && nrow(sc2_long) > 0) {
  # Filter out rows with NA in Value or y
  sc2_long_clean <- sc2_long[!is.na(sc2_long$Value) & !is.na(sc2_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc2_main <- ggplot(sc2_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
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

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc2 <- wrap_elements(p_cov_sc2_top / p_cov_sc2_main + plot_layout(heights = c(1, 1)))
  has_sc2_plot <- TRUE
} else {
  has_sc2_plot <- FALSE
}

# Combine sc1 and sc2 plots vertically
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
  sc3_long_clean <- sc3_long[!is.na(sc3_long$Value) & !is.na(sc3_long$y), ]

  # Main contour plot: Covariate (x-axis) vs Outcome (y-axis)
  p_cov_sc3_main <- ggplot(sc3_long_clean, aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Configuration_F + Covariate ~ Cor_Label, scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
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

  # Combine with patchwork: histogram grid above main contour
  p_cov_sc3 <- wrap_elements(p_cov_sc3_top / p_cov_sc3_main + plot_layout(heights = c(2, 3)))
  has_sc3_plot <- TRUE
} else {
  has_sc3_plot <- FALSE
}

# LRC-BART ADDITION START
# Create the matching covariate/outcome display for Sc4 and Sc5.
if (!is.null(sc45_long) && nrow(sc45_long) > 0) {
  sc45_long_clean <- sc45_long[
    !is.na(sc45_long$Value) & !is.na(sc45_long$y), ]

  p_cov_sc45_main <- ggplot(sc45_long_clean,
                            aes(x = Value, y = y, color = Group)) +
    geom_point(alpha = 0.1, size = 0.5, na.rm = TRUE) +
    geom_density_2d(alpha = 0.8, linewidth = 0.5, na.rm = TRUE) +
    facet_nested(Scenario_F + Configuration_F + Covariate ~ Cor_Label,
                 scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    labs(x = "Covariate Value", y = "Y") +
    scale_color_manual(values = group_colors, name = "Group") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA,
                                      linewidth = 0.5))

  p_cov_sc45_top <- ggplot(sc45_long_clean, aes(x = Value, fill = Group)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.3,
                   alpha = 0.6, position = "identity", na.rm = TRUE) +
    facet_nested(Scenario_F + Configuration_F + Covariate ~ Cor_Label,
                 scales = "free",
                 nest_line = element_line(linewidth = 1, color = "black")) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "Covariate Distributions: Scenarios 4 and 5",
         x = "Value", y = "Density") +
    theme_classic(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          panel.spacing.x = unit(0.3, "lines"),
          panel.spacing.y = unit(0.3, "lines"),
          panel.border = element_rect(color = "grey50", fill = NA,
                                      linewidth = 0.5))

  p_cov_sc45 <- wrap_elements(
    p_cov_sc45_top / p_cov_sc45_main + plot_layout(heights = c(1, 1)))
  has_sc45_plot <- TRUE
} else {
  has_sc45_plot <- FALSE
}
# LRC-BART ADDITION END

# Dynamically set heights based on number of scenarios
outcome_height <- max(5.0, outcome_nrow * 1.0)

# Calculate covariate plot heights based on number of rows (Configuration + Covariate)
n_sc1_configuration <- length(unique(sc1_data$Configuration_F[!is.na(sc1_data$Configuration_F)]))
n_sc2_configuration <- length(unique(sc2_data$Configuration_F[!is.na(sc2_data$Configuration_F)]))
n_sc3_configuration <- length(unique(sc3_data$Configuration_F[!is.na(sc3_data$Configuration_F)]))
n_sc45_configuration <- length(unique(sc45_data$Configuration_F[!is.na(sc45_data$Configuration_F)]))
# sc1 and sc2 are side by side, so use max of their heights
cov_height_sc12 <- max(n_sc1_configuration, n_sc2_configuration, n_sc3_configuration) * 10
cov_height_sc3 <- max(n_sc1_configuration, n_sc2_configuration, n_sc3_configuration) * 10
cov_height_sc45 <- max(1, n_sc45_configuration) * 10

# Create a plot from the legend
legend_plot_final <- ggplot(dummy_data, aes(x = x, fill = Group)) +
  geom_bar() +
  scale_fill_manual(values = group_colors, name = "Group") +
  theme_void() +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9))

# Arrange the plots using patchwork (no KM plots for linear models).
plot_parts <- list(legend_plot_final, p_outcome)
plot_heights <- c(0.5, outcome_height)
if (has_sc12_plot) {
  plot_parts[[length(plot_parts) + 1]] <- p_cov_sc12
  plot_heights <- c(plot_heights, cov_height_sc12)
}
if (has_sc3_plot) {
  plot_parts[[length(plot_parts) + 1]] <- p_cov_sc3
  plot_heights <- c(plot_heights, cov_height_sc3)
}
if (has_sc45_plot) {
  plot_parts[[length(plot_parts) + 1]] <- p_cov_sc45
  plot_heights <- c(plot_heights, cov_height_sc45)
}
plot <- wrap_plots(plot_parts, ncol = 1, heights = plot_heights)
total_height <- sum(plot_heights) + 1

# Calculate total plot dimensions dynamically
total_width <- outcome_ncol * 3 + 2  # Scale width with number of columns, +2 for margins

# Include data_folder in output filename to distinguish plots from different data sources
output_suffix <- ifelse(data_folder == "data", "", paste0("_", data_folder))
save_pipeline_plot(paste0(file.path(mainDir),"/lrcbart-sim-gaussian-single-arm/inserts/balance_all_scenarios", output_suffix, ".jpg"),
       width = total_width,
       height = total_height,
       plot,
       limitsize = FALSE)

finish_pipeline_plot()
