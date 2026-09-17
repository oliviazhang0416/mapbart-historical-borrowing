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
library(survival)
library(survminer)
library(cowplot)
library(patchwork)
library(tidyr)

# LRC-BART MODIFICATION START
# The inherited survival balance plot is extended with the revision's Sc4
# global-shift and Sc5 regional-partial-compatibility configurations.
mainDir <- .lrcRoot
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

# Load data for all scenarios (without configuration distinction)
all_data <- data.frame()

for (sc in selected_sc(1:5)) {
  if (sc %in% c(1, 2, 4, 5)){
    # LRC-BART ADDITION START
    if (sc %in% c(1, 2)) configurations <- c("Base")
    if (sc == 4) configurations <- c("d1", "d2")
    if (sc == 5)
      configurations <- c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2")
    # LRC-BART ADDITION END
    cor <- 1
  } else {
    configurations <- c("Base")
    cor <- c(-0.5, 0, 0.5)
  }

  if (sc %in% c(1, 2, 4, 5)) {
    for (configuration in selected_configurations(sc, configurations)) {
      for (c in selected_correlations(sc, cor)) {
        for (hypo in plot_hypothesis) {
          # LRC-BART ADDITION START
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
          filename_rct <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,
                                 "/data_p",p_obs,"_",scenario_tag,"_",hypo,
                                 "_",iter_rct,".RData")
          filename_rwd <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,
                                 "/data_p",p_obs,"_",scenario_tag,"_",hypo,
                                 "_",iter_rwd,".RData")
          # LRC-BART ADDITION END

          data_rct_full <- if (file.exists(filename_rct)) tryCatch({
            read_plot_data(filename_rct)
          }, error = function(e) {
            return(NULL)
          }) else NULL

          data_rwd_full <- if (file.exists(filename_rwd)) tryCatch({
            read_plot_data(filename_rwd)
          }, error = function(e) {
            return(NULL)
          }) else NULL

          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y,
                                        "delta" = data_rct_full$delta))
            data_rct$U1 <- NA
            data_rct$U2 <- NA
            data_rct <- data_rct[data_rct$D == 1, ]

            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y,
                                        "delta" = data_rwd_full$delta))
            data_rwd$U1 <- NA
            data_rwd$U2 <- NA
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            data <- rbind(data_rct, data_rwd)
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
    for (configuration in selected_configurations(sc, configurations)) {
      for (c in selected_correlations(sc, cor)) {
        for (hypo in plot_hypothesis) {
          if (configuration == "Base") {
            filename_rct <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,"/data_p",p_obs,"_sc",sc,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
          } else {
            filename_rct <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,"/data_p",p_obs,"_sc",sc,configuration,"_cor",c,"_",hypo,"_",iter_rct,".RData")
            filename_rwd <- paste0(mainDir,"/lrcbart-sim-survival/",data_folder,"/data_p",p_obs,"_sc",sc,configuration,"_cor",c,"_",hypo,"_",iter_rwd,".RData")
          }

          data_rct_full <- if (file.exists(filename_rct)) tryCatch({
            read_plot_data(filename_rct)
          }, error = function(e) {
            return(NULL)
          }) else NULL

          data_rwd_full <- if (file.exists(filename_rwd)) tryCatch({
            read_plot_data(filename_rwd)
          }, error = function(e) {
            return(NULL)
          }) else NULL

          if (!is.null(data_rct_full) && !is.null(data_rwd_full)) {
            data_rct <- cbind(data_rct_full$X,
                              data.frame("y" = data_rct_full$y,
                                        "delta" = data_rct_full$delta))
            if (!is.null(data_rct_full$U)) {
              data_rct <- cbind(data_rct, data_rct_full$U)
            } else {
              data_rct$U1 <- NA
              data_rct$U2 <- NA
            }
            data_rct <- data_rct[data_rct$D == 1, ]

            data_rwd <- cbind(data_rwd_full$X,
                              data.frame("y" = data_rwd_full$y,
                                        "delta" = data_rwd_full$delta))
            if (!is.null(data_rwd_full$U)) {
              data_rwd <- cbind(data_rwd, data_rwd_full$U)
            } else {
              data_rwd$U1 <- NA
              data_rwd$U2 <- NA
            }
            data_rwd <- data_rwd[data_rwd$D == 0, ]

            data <- rbind(data_rct, data_rwd)
            data$Scenario <- sc
            data$Configuration <- "Base"
            data$Correlation <- c
            data$Hypothesis <- hypo

            all_data <- rbind(all_data, data)
          }
        }
      }
    }
  }
}

if (nrow(all_data) == 0) {
  stop("No data files were found.")
}

# Define color palette
group_colors <- c("Control" = "darkblue", "External Control" = "lightblue", "Treatment" = "red")

# Create Group factor
all_data$Group <- factor(ifelse(all_data$Z == 1, "Treatment",
                            ifelse(all_data$Z == 0 & all_data$D == 1, "Control", "External Control")),
                     levels = c("Control", "External Control", "Treatment"))

all_data$Hypothesis <- factor(all_data$Hypothesis, levels = c("null", "alternative"))

# LRC-BART ADDITION START
# Create one panel label per scenario configuration.
all_data$Scenario_Label <- ifelse(
  all_data$Scenario == 3,
  paste0("Scenario 3\n(ρ=", all_data$Correlation, ")"),
  ifelse(all_data$Scenario %in% c(4, 5),
         paste0("Scenario ", all_data$Scenario, "\n", all_data$Configuration),
         paste0("Scenario ", all_data$Scenario))
)
all_data$Scenario_Label <- plot_factor(all_data$Scenario_Label,
                                   levels = c("Scenario 1", "Scenario 2",
                                              "Scenario 3\n(ρ=-0.5)", "Scenario 3\n(ρ=0)", "Scenario 3\n(ρ=0.5)",
                                              paste0("Scenario 4\ndelta RWD = ", c(1, 2)),
                                              paste0("Scenario 5\nregion = X5, delta RWD = ", c(0.5, 1, 2)),
                                              paste0("Scenario 5\nregion = X7, delta RWD = ", c(1, 2))))
# LRC-BART ADDITION END

# Reorder data so External Control is plotted last (on top)
all_data_reordered <- all_data %>% arrange(desc(Group == "External Control"))

# ============================================================================
# PLOT 1: Kaplan-Meier Survival Curves
# ============================================================================

# Create KM plots for each scenario configuration.
km_plots <- list()

for (scenario_label in levels(all_data$Scenario_Label)) {
  df <- all_data[all_data$Scenario_Label == scenario_label, ]
  if (nrow(df) == 0) next

  fit <- survfit(Surv(y, delta) ~ Group, data = df)

  p <- ggsurvplot(fit, data = df,
                  palette = c("darkblue", "lightblue", "red"),
                  legend = "none",
                  xlab = "Time",
                  ylab = "Survival Probability",
                  title = scenario_label,
                  ggtheme = theme_classic(base_size = 11),
                  conf.int = FALSE,
                  risk.table = FALSE)$plot +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

  km_plots[[length(km_plots) + 1]] <- p
}

# Combine KM plots (no legend)
p_km_final <- wrap_plots(km_plots, ncol = 4)

# Save KM survival plot
output_suffix <- ifelse(data_folder == "data", "", paste0("_", data_folder))
save_pipeline_plot(paste0(file.path(mainDir),"/lrcbart-sim-survival/inserts/balanced_1_surv.jpg"),
       width = 16,
       height = 3 * ceiling(length(km_plots) / 4),
       p_km_final,
       limitsize = FALSE)

# ============================================================================
# PLOT 2: Covariate Balance Plot
# ============================================================================

# Reshape data for covariate histograms
# For sc == 1 and 2: X5, X6
# For sc == 3: X5, U1, X6, U2

sc1_data <- all_data_reordered[all_data_reordered$Scenario == 1, ]
sc2_data <- all_data_reordered[all_data_reordered$Scenario == 2, ]
sc3_data <- all_data_reordered[all_data_reordered$Scenario == 3, ]
# LRC-BART ADDITION START
sc45_data <- all_data_reordered[all_data_reordered$Scenario %in% c(4, 5), ]
# LRC-BART ADDITION END

# Define covariate label mapping
cov_labels <- c("X5" = "X₅ (Measured)", "X6" = "X₆ (Measured)",
                "U1" = "U₁ (Unmeasured)", "U2" = "U₂ (Unmeasured)")

# Reshape to long format
if (nrow(sc1_data) > 0) {
  sc1_long <- sc1_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc1_long$Covariate <- factor(cov_labels[sc1_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc1_long <- NULL
}

if (nrow(sc2_data) > 0) {
  sc2_long <- sc2_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc2_long$Covariate <- factor(cov_labels[sc2_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc2_long <- NULL
}

if (nrow(sc3_data) > 0 && "U1" %in% colnames(sc3_data) && !all(is.na(sc3_data$U1))) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_Label, Group, X5, U1, X6, U2) %>%
    pivot_longer(cols = c(X5, U1, X6, U2),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(cov_labels[sc3_long$Covariate],
                                levels = c("X₅ (Measured)", "U₁ (Unmeasured)",
                                           "X₆ (Measured)", "U₂ (Unmeasured)"))
} else if (nrow(sc3_data) > 0) {
  sc3_long <- sc3_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc3_long$Covariate <- factor(cov_labels[sc3_long$Covariate],
                                levels = c("X₅ (Measured)", "X₆ (Measured)"))
} else {
  sc3_long <- NULL
}

# LRC-BART ADDITION START
if (nrow(sc45_data) > 0) {
  sc45_long <- sc45_data %>%
    dplyr::select(Scenario_Label, Group, X5, X6) %>%
    pivot_longer(cols = c(X5, X6),
                 names_to = "Covariate",
                 values_to = "Value")
  sc45_long$Covariate <- factor(cov_labels[sc45_long$Covariate],
                                levels = c("X₅ (Measured)",
                                           "X₆ (Measured)"))
} else {
  sc45_long <- NULL
}
# LRC-BART ADDITION END

# Combine all scenarios
all_long <- rbind(sc1_long, sc2_long, sc3_long, sc45_long)
all_long <- all_long[!is.na(all_long$Value), ]

# Create covariate balance histogram plot
p_balance <- ggplot(all_long, aes(x = Value, fill = Group)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, alpha = 0.6, position = "identity") +
  facet_grid(Covariate ~ Scenario_Label, scales = "free") +
  scale_fill_manual(values = group_colors, name = "Group") +
  labs(x = "Covariate Value", y = "Density") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

# Determine height based on number of covariates (no legend)
n_covariates <- length(unique(all_long$Covariate))
plot_height <- n_covariates * 2

# Save covariate balance plot for Sc1-Sc5 configurations.
save_pipeline_plot(paste0(file.path(mainDir),"/lrcbart-sim-survival/inserts/balanced_2_surv.jpg"),
       width = 24,
       height = plot_height,
       p_balance,
       limitsize = FALSE)

finish_pipeline_plot()
