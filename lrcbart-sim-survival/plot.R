rm(list=ls())

# Resolve the repository root by walking up from the working directory.
.lrcRoot <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (!file.exists(file.path(d, ".lrcbart-root")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, ".lrcbart-root")))
    stop("lrcbart repository root not found from ", getwd())
  d
})
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)
# LRC-BART MODIFICATION START
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-survival")
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

# hierAFT (parametric) keeps a fixed prior; LRC-BART is calibrated per target
# and result files use _N<target>_ suffixes.
aftv4_prior_vals <- plot_method_setting("hierAFT.R", "prior_vals", 0.25)
# LRC-BART targets and fixed-weight variants are selected inside each scenario.
# Backward-compatibility alias used by hierAFT loops below.
prior_vals       <- aftv4_prior_vals

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

# =============================================================================
# build_rmst_column(sc)  --  RMST(tau)-ratio column (right side of the unified
# figure). RMST performance facet + summary table, arm-specific population RMST
# facet + tables, and (if present) subject-wise RMST facet + table. Reads the
# RMST columns fresh from res/*.RData. Returns a patchwork stack, or NULL if
# absent. Uses `configuration` from the enclosing scope (matches the median column).
# =============================================================================
build_rmst_column <- function(sc) {
  lrcbart_configurations <- lrc_plot_configurations(FALSE, sc)
  RMST_COLS <- c("c","iteration","rmst_tau","rmst_true","rmst_hat",
                 "bias_rmst","sd_rmst","rmse_rmst","w1distance_rmst","w2distance_rmst",
                 "ci_rmst","coverage_rmst","tp_rmst","tp_calibrated_rmst",
                 "bias.trt.rmst.pop","w2distance.trt.rmst.pop",
                 "bias.ctrl.rmst.pop","w2distance.ctrl.rmst.pop",
                 "bias_subj_rmst","pehe_subj_rmst")
  resdir <- paste0(file.path(projectDir, "res"), "/")
  read_rmst <- function(file, method) {
    if (!file.exists(file)) return(NULL)
    d <- tryCatch(read_plot_result(file), error = function(e) NULL)
    if (is.null(d) || !"bias_rmst" %in% names(d)) return(NULL)
    out <- d[, intersect(RMST_COLS, names(d)), drop = FALSE]; out$Method <- method; out
  }
  # LRC-BART MODIFICATION START
  build_one <- function(scenario_tag, scenario_label, cor_tag = NULL) {
    suf <- function(extra = "")
      paste0("_", scenario_tag, extra, paste0("_", plot_hypothesis, ".RData"))
    fr <- list(read_rmst(paste0(resdir,"AFTv1_p",p_obs,suf()),"AFTv1"),
               read_rmst(paste0(resdir,"AFTv2_p",p_obs,suf()),"AFTv2"),
               read_rmst(paste0(resdir,"AFTv3_p",p_obs,suf()),"AFTv3"),
               read_rmst(paste0(resdir,"BARTv1_p",p_obs,suf()),"BARTv1"),
               read_rmst(paste0(resdir,"BARTv2_p",p_obs,suf()),"BARTv2"),
               read_rmst(paste0(resdir,"BARTv3_p",p_obs,suf()),"BARTv3"))
    for (pr in aftv4_prior_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"hierAFT_p",p_obs,suf(paste0("_prior",pr))), paste0("hierAFT(",pr,")"))))
    for (li in seq_len(nrow(lrcbart_configurations)))
      fr <- c(fr, list(read_rmst(paste0(resdir, "LRC-BART_p", p_obs,
        suf(paste0(lrcbart_configurations$suffix[li], "_N", lrcbart_configurations$target[li]))),
        lrcbart_configurations$label[li])))
    fr <- fr[!sapply(fr, is.null)]; if (!length(fr)) return(NULL)
    out <- bind_rows(fr)
    out$Configuration <- scenario_label
    out$cor <- if (is.null(cor_tag)) NA_character_ else as.character(cor_tag)
    out
  }
  if (sc == 1) dat <- build_one("sc1", "Base")
  if (sc == 2) dat <- build_one("sc2", "Base")
  if (sc == 3) dat <- bind_rows(lapply(c(-0.5, 0, 0.5), function(cv)
    build_one(paste0("sc3_cor", cv), "Base", cv)))
  if (sc == 4) dat <- bind_rows(lapply(c(1, 2), function(delta)
    build_one(paste0("sc4_d", delta), paste0("delta RWD = ", delta))))
  if (sc == 5) {
    tags <- c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2")
    dat <- bind_rows(lapply(tags, function(tag) {
      parts <- strsplit(tag, "_d", fixed = TRUE)[[1]]
      build_one(paste0("sc5_", tag),
                paste0("region = ", parts[1], ", delta RWD = ", parts[2]))
    }))
  }
  if (is.null(dat) || !nrow(dat)) return(NULL)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  grp <- if (sc == 3) c("Method", "Configuration", "cor") else
    c("Method", "Configuration")
  # Finalize a summary table the SAME way as the median-ratio tables: insert
  # separator rows, and for sc3 keep a plain "Default" label plus a separate
  # Correlation column, separating rows by the Configuration+rho combination
  # (the combined "..., rho=" label is used only to place separators, then dropped).
  fin_tab <- function(tb, label_col = "Configuration_Label") {
    if (sc == 3) {
      tb <- tb %>% arrange(.data[[label_col]], Correlation, Method) %>%
        mutate(.sep = paste0(.data[[label_col]], ", rho=", Correlation))
      add_configuration_separators(tb, ".sep") %>% dplyr::select(-.sep)
    } else add_configuration_separators(tb, label_col)
  }
  for (.cc in c("rmse_rmst","w1distance_rmst","w2distance_rmst","tp_calibrated_rmst"))
    if (!.cc %in% names(dat)) dat[[.cc]] <- NA_real_
  .fmt  <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", m) }
  .fmtR <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", sqrt(m)) }
  rmst_main <- dat %>% group_by(across(all_of(grp))) %>%
    summarise(Bias=.fmt(bias_rmst), SD=.fmt(sd_rmst), RMSE=.fmtR(rmse_rmst),
              W1Distance=.fmt(w1distance_rmst), W2Distance=.fmt(w2distance_rmst),
              CI_length=.fmt(ci_rmst), CI_coverage=.fmt(coverage_rmst),
              Power=.fmt(tp_rmst), Power_calib=.fmt(tp_calibrated_rmst),
              N=n(), N_bias_ge1=sum(abs(bias_rmst)>=1,na.rm=TRUE), .groups="drop") %>%
    mutate(Configuration_Label=factor(paste0("Config ", Configuration)), .before=1)
  if (sc==3) rmst_main <- rmst_main %>% rename(Correlation=cor) %>% relocate(Correlation, .after=Configuration_Label)
  table_grob <- pipeline_table(fin_tab(rmst_main), rows = NULL)
  long <- dat %>% pivot_longer(c(bias_rmst,sd_rmst,ci_rmst,coverage_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
  long$metric <- factor(long$metric, levels=c("bias_rmst","sd_rmst","ci_rmst","coverage_rmst"), labels=c("Bias","SD","CI length","Coverage"))
  thm <- theme_minimal() + theme(plot.title=element_text(size=14,face="bold",hjust=0.5), axis.text.x=element_text(size=9,angle=45,hjust=1), legend.position="none", strip.text=element_text(size=11,face="bold"), strip.background=element_rect(fill="gray90",color="gray50"), panel.border=element_rect(color="black",fill=NA,linewidth=0.5))
  p <- ggplot(long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("RMST(τ)-ratio performance — Scenario ",sc," (τ = ",sprintf("%.1f",mean(dat$rmst_tau,na.rm=TRUE)),")"), x=NULL, y=NULL) + thm
  p <- if (sc==3) p + facet_grid(cor~metric,scales="free_y") else
    p + facet_grid(Configuration~metric,scales="free_y")
  arm_long <- dat %>% pivot_longer(c(bias.trt.rmst.pop,w2distance.trt.rmst.pop,bias.ctrl.rmst.pop,w2distance.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control")), Metric=ifelse(grepl("^bias",metric),"Bias","W2Distance"))
  arm_summary <- arm_long %>% group_by(across(all_of(c(grp, "Arm", "Metric")))) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=Metric, values_from=m)
  mk_arm <- function(g) {
    t <- arm_summary %>% filter(Arm==g) %>% dplyr::select(-Arm) %>%
      mutate(Configuration_Label=factor(paste0("Config ", Configuration)), Group=g)
    if (sc==3) t %>% rename(Correlation=cor) %>% dplyr::select(Group, Configuration_Label, Correlation, Method, Bias, W2Distance)
    else       t %>% dplyr::select(Group, Configuration_Label, Method, Bias, W2Distance)
  }
  arm_tab_trt <- mk_arm("Treatment"); arm_tab_ctrl <- mk_arm("Control")
  table_grob_arm <- wrap_plots(pipeline_table(fin_tab(arm_tab_trt),rows=NULL), pipeline_table(fin_tab(arm_tab_ctrl),rows=NULL), ncol=2)
  # Extra fixed-weight methods need table space, including Sc3 separator rows.
  n_configuration_rows <- nrow(unique(dat[, c("Configuration", "cor")]))
  table_height <- function(d) max(3, 0.35 * (nrow(fin_tab(d)) + 1))
  rmst_heights <- c(max(6, 3 * n_configuration_rows), table_height(rmst_main),
                    max(7, 5 * n_configuration_rows),
                    max(table_height(arm_tab_trt), table_height(arm_tab_ctrl)))
  p_arm <- ggplot(arm_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("Arm-specific population RMST(τ) — Scenario ",sc,"  (estimate vs DGP truth, per arm)"), x=NULL,y=NULL) + thm
  p_arm <- if (sc==3) p_arm + facet_grid(Arm+cor~Metric,scales="free_y") else
    p_arm + facet_grid(Arm+Configuration~Metric,scales="free_y")
  has_subj <- all(c("bias_subj_rmst","pehe_subj_rmst") %in% names(dat))
  if (has_subj) {
    subj_long <- dat %>% pivot_longer(c(bias_subj_rmst,pehe_subj_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
    subj_long$metric <- factor(subj_long$metric, levels=c("bias_subj_rmst","pehe_subj_rmst"), labels=c("Bias","PEHE"))
    subj_tab <- subj_long %>% group_by(across(all_of(c(grp, "metric")))) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=metric, values_from=m) %>%
      mutate(Configuration_Label=factor(paste0("Config ", Configuration)))
    subj_tab <- if (sc==3) subj_tab %>% rename(Correlation=cor) %>% dplyr::select(Configuration_Label, Correlation, Method, Bias, PEHE)
                else        subj_tab %>% dplyr::select(Configuration_Label, Method, Bias, PEHE)
    table_grob_subj <- pipeline_table(fin_tab(subj_tab), rows = NULL)
    p_subj <- ggplot(subj_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
      geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
      labs(title=paste0("Subject-wise RMST(τ)-ratio — Scenario ",sc,"  (per-subject estimate vs truth)"), x=NULL,y=NULL) + thm
    p_subj <- if (sc==3) p_subj + facet_grid(cor~metric,scales="free_y") else
      p_subj + facet_grid(Configuration~metric,scales="free_y")
    rmst_heights <- c(rmst_heights, max(6, 3 * n_configuration_rows),
                      table_height(subj_tab), 7, 3)
    column <- p / table_grob / p_arm / table_grob_arm / p_subj / table_grob_subj /
             patchwork::plot_spacer() / patchwork::plot_spacer() +
             plot_layout(heights = rmst_heights)
    attr(column, "required_height") <- sum(rmst_heights)
    return(column)
  }
  # LRC-BART MODIFICATION END
  rmst_heights <- c(rmst_heights, 7, 3)
  column <- p / table_grob / p_arm / table_grob_arm /
    patchwork::plot_spacer() / patchwork::plot_spacer() +
    plot_layout(heights = rmst_heights)
  attr(column, "required_height") <- sum(rmst_heights)
  column
}

for (sc in selected_sc(1:5)) {
lrcbart_configurations <- lrc_plot_configurations(FALSE, sc)
tryCatch({

if (sc %in% c(1, 2, 4, 5)){

  # Initialize empty data frames
  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()

  # LRC-BART ADDITION START
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

    # Try to read AFTv1 results
    tryCatch({
      AFTv1_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-survival/res/AFTv1_p",p_obs,file_suffix))
      AFTv1_res_ATE <- AFTv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv1_res_sigma <- AFTv1_res[, c("c","iteration",
                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                "pehe_subj","bias_subj")]
      AFTv1_res_ATE$Method <- "AFTv1"
      AFTv1_res_sigma$Method <- "AFTv1"
      AFTv1_res_ATE$Configuration <- configuration_label
      AFTv1_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv1_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv1_res_sigma
    }, error = function(e) {})

    # Try to read AFTv2 results
    tryCatch({
      AFTv2_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-survival/res/AFTv2_p",p_obs,file_suffix))
      AFTv2_res_ATE <- AFTv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv2_res_sigma <- AFTv2_res[, c("c","iteration",
                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                  "pehe_subj","bias_subj")]
      AFTv2_res_ATE$Method <- "AFTv2"
      AFTv2_res_sigma$Method <- "AFTv2"
      AFTv2_res_ATE$Configuration <- configuration_label
      AFTv2_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv2_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv2_res_sigma
    }, error = function(e) {})

    # Try to read AFTv3 results
    tryCatch({
      AFTv3_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-survival/res/AFTv3_p",p_obs,file_suffix))
      AFTv3_res_ATE <- AFTv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
      AFTv3_res_sigma <- AFTv3_res[, c("c","iteration",
                                    "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                    "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                    "bias.trt.median.pop","w2distance.trt.median.pop",
                                    "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                    "pehe_subj","bias_subj")]
      AFTv3_res_ATE$Method <- "AFTv3"
      AFTv3_res_sigma$Method <- "AFTv3"
      AFTv3_res_ATE$Configuration <- configuration_label
      AFTv3_res_sigma$Configuration <- configuration_label
      configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv3_res_ATE
      configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv3_res_sigma
    }, error = function(e) {})

    # Try to read hierAFT results for each prior value
    for (prior_val in prior_vals) {
      tryCatch({
        aftv4_file_suffix <- paste0(scenario_suffix, "_prior", prior_val, paste0("_", plot_hypothesis, ".RData"))
        file_path <- paste0(mainDir,"/lrcbart-sim-survival/res/hierAFT_p",p_obs,aftv4_file_suffix)
        if (file.exists(file_path)) {
          hierAFT_res <- read_plot_result(file_path)
          hierAFT_res_ATE <- hierAFT_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          hierAFT_res_sigma <- hierAFT_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          hierAFT_res_ATE$Method <- paste0("hierAFT(", prior_val, ")")
          hierAFT_res_sigma$Method <- paste0("hierAFT(", prior_val, ")")
          hierAFT_res_ATE$Configuration <- configuration_label
          hierAFT_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- hierAFT_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- hierAFT_res_sigma
        }
      }, error = function(e) {})
    }

    # Try to read BARTv1 / BARTv2 results
    for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
    tryCatch({
      BART_res <- read_plot_result(paste0(mainDir,"/lrcbart-sim-survival/res/",bver,"_p",p_obs,file_suffix))
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
    }, error = function(e) {})
    }

    # Try to read LRC-BART results for each N_target value
    for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

      tryCatch({
        lrcbart_file_suffix <- paste0(scenario_suffix, lrc_suffix,"_N", target_N, paste0("_", plot_hypothesis, ".RData"))
        file_path <- paste0(mainDir,"/lrcbart-sim-survival/res/LRC-BART_p",p_obs,lrcbart_file_suffix)
        if (file.exists(file_path)) {
          LRC_BART_res <- read_plot_result(file_path)
          LRC_BART_res <- LRC_BART_res[LRC_BART_res$alpha == 0.95 & LRC_BART_res$beta == 2, ]
          LRC_BART_res_ATE <- LRC_BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          LRC_BART_res_sigma <- LRC_BART_res[, c("c","iteration",
                                                "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                "bias.trt.median.pop","w2distance.trt.median.pop",
                                                "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                "pehe_subj","bias_subj")]
          LRC_BART_res_ATE$Method <- lrc_label
          LRC_BART_res_sigma$Method <- lrc_label
          LRC_BART_res_ATE$Configuration <- configuration_label
          LRC_BART_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LRC_BART_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LRC_BART_res_sigma
        }
      }, error = function(e) {})
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
      list(name = "AFTv1", file_prefix = "AFTv1"),
      list(name = "AFTv2", file_prefix = "AFTv2"),
      list(name = "AFTv3", file_prefix = "AFTv3"),
      list(name = "BARTv1", file_prefix = "BARTv1"),
      list(name = "BARTv2", file_prefix = "BARTv2"),
      list(name = "BARTv3", file_prefix = "BARTv3")
    )

    for (method_info in methods_list) {
      tryCatch({
        null_file <- paste0(mainDir,"/lrcbart-sim-survival/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
        if (file.exists(null_file)) {
          null_res <- read_plot_result(null_file)

          # For BART, filter by alpha and beta
          if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3")) {
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          }

          # Calculate mean FP (Type I error)
          if ("fp" %in% colnames(null_res)) {
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

    # Try to load hierAFT null results for each prior value
    for (prior_val in prior_vals) {
      tryCatch({
        aftv4_null_suffix <- paste0(scenario_suffix, "_prior", prior_val, "_null.RData")
        null_file <- paste0(mainDir,"/lrcbart-sim-survival/res/hierAFT_p", p_obs, aftv4_null_suffix)
        if (file.exists(null_file)) {
          null_res <- read_plot_result(null_file)
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              Configuration = configuration_label,
              Method = paste0("hierAFT(", prior_val, ")"),
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {})
    }

    # Try to load LRC-BART null results for each N_target value
    for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

      tryCatch({
        lrcbart_null_suffix <- paste0(scenario_suffix, lrc_suffix,"_N", target_N, "_null.RData")
        null_file <- paste0(mainDir,"/lrcbart-sim-survival/res/LRC-BART_p", p_obs, lrcbart_null_suffix)
        if (file.exists(null_file)) {
          null_res <- read_plot_result(null_file)
          null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
          if ("fp" %in% colnames(null_res)) {
            mean_FP <- mean(null_res$fp, na.rm = TRUE)
            all_null_FP <- rbind(all_null_FP, data.frame(
              Configuration = configuration_label,
              Method = lrc_label,
              Type_I_error = sprintf("%.2f", mean_FP)
            ))
          }
        }
      }, error = function(e) {})
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

      # Try to read AFTv1 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-survival/res/AFTv1_p",p_obs,file_suffix)
        files_checked[[length(files_checked) + 1]] <- file_path
        if (file.exists(file_path)) {
          AFTv1_res <- read_plot_result(file_path)
          files_found[[length(files_found) + 1]] <- file_path
          AFTv1_res_ATE <- AFTv1_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv1_res_sigma <- AFTv1_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv1_res_ATE$Method <- "AFTv1"
          AFTv1_res_sigma$Method <- "AFTv1"
          AFTv1_res_ATE$Configuration <- configuration_label
          AFTv1_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv1_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv1_res_sigma
        }
      }, error = function(e) {})

      # Try to read AFTv2 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-survival/res/AFTv2_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          AFTv2_res <- read_plot_result(file_path)
          AFTv2_res_ATE <- AFTv2_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv2_res_sigma <- AFTv2_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv2_res_ATE$Method <- "AFTv2"
          AFTv2_res_sigma$Method <- "AFTv2"
          AFTv2_res_ATE$Configuration <- configuration_label
          AFTv2_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv2_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv2_res_sigma
        }
      }, error = function(e) {})

      # Try to read AFTv3 results
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-survival/res/AFTv3_p",p_obs,file_suffix)
        if (file.exists(file_path)) {
          AFTv3_res <- read_plot_result(file_path)
          AFTv3_res_ATE <- AFTv3_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
          AFTv3_res_sigma <- AFTv3_res[, c("c","iteration",
                                        "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                        "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                        "bias.trt.median.pop","w2distance.trt.median.pop",
                                        "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                        "pehe_subj","bias_subj")]
          AFTv3_res_ATE$Method <- "AFTv3"
          AFTv3_res_sigma$Method <- "AFTv3"
          AFTv3_res_ATE$Configuration <- configuration_label
          AFTv3_res_sigma$Configuration <- configuration_label
          configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- AFTv3_res_ATE
          configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- AFTv3_res_sigma
        }
      }, error = function(e) {})

      # Try to read hierAFT results for each prior value
      for (prior_val in prior_vals) {
        tryCatch({
          aftv4_file_suffix <- paste0("_sc", sc,
            if (configuration == "Base") "" else configuration,
            "_cor", cor_val, "_prior", prior_val, paste0("_", plot_hypothesis, ".RData"))
          file_path <- paste0(mainDir,"lrcbart-sim-survival/res/hierAFT_p",p_obs,aftv4_file_suffix)
          if (file.exists(file_path)) {
            hierAFT_res <- read_plot_result(file_path)
            hierAFT_res_ATE <- hierAFT_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            hierAFT_res_sigma <- hierAFT_res[, c("c","iteration",
                                          "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                          "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                          "bias.trt.median.pop","w2distance.trt.median.pop",
                                          "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                          "pehe_subj","bias_subj")]
            hierAFT_res_ATE$Method <- paste0("hierAFT(", prior_val, ")")
            hierAFT_res_sigma$Method <- paste0("hierAFT(", prior_val, ")")
            hierAFT_res_ATE$Configuration <- configuration_label
            hierAFT_res_sigma$Configuration <- configuration_label
            configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- hierAFT_res_ATE
            configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- hierAFT_res_sigma
          }
        }, error = function(e) {})
      }

      # Try to read BARTv1 / BARTv2 results
      for (bver in c("BARTv1", "BARTv2", "BARTv3")) {
      tryCatch({
        file_path <- paste0(mainDir,"lrcbart-sim-survival/res/",bver,"_p",p_obs,file_suffix)
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
      }, error = function(e) {})
      }

      # Try to read LRC-BART results for each N_target value
      for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

        tryCatch({
          lrcbart_file_suffix <- paste0("_sc", sc,
            if (configuration == "Base") "" else configuration,
            "_cor", cor_val, lrc_suffix,"_N", target_N, paste0("_", plot_hypothesis, ".RData"))
          file_path <- paste0(mainDir,"lrcbart-sim-survival/res/LRC-BART_p",p_obs,lrcbart_file_suffix)
          if (file.exists(file_path)) {
            LRC_BART_res <- read_plot_result(file_path)
            LRC_BART_res <- LRC_BART_res[LRC_BART_res$alpha == 0.95 & LRC_BART_res$beta == 2, ]
            LRC_BART_res_ATE <- LRC_BART_res[, c("c","iteration","bias","sd","rmse","w1distance","w2distance","ci","coverage","tp","fp","tp_calibrated")]
            LRC_BART_res_sigma <- LRC_BART_res[, c("c","iteration",
                                                  "bias.trt.sigma","sd.trt.sigma","w2distance.trt.sigma",
                                                  "bias.ctrl.sigma","sd.ctrl.sigma","w2distance.ctrl.sigma",
                                                  "bias.trt.median.pop","w2distance.trt.median.pop",
                                                  "bias.ctrl.median.pop","w2distance.ctrl.median.pop",
                                                  "pehe_subj","bias_subj")]
            LRC_BART_res_ATE$Method <- lrc_label
            LRC_BART_res_sigma$Method <- lrc_label
            LRC_BART_res_ATE$Configuration <- configuration_label
            LRC_BART_res_sigma$Configuration <- configuration_label
            configuration_res_list_ATE[[length(configuration_res_list_ATE) + 1]] <- LRC_BART_res_ATE
            configuration_res_list_sigma[[length(configuration_res_list_sigma) + 1]] <- LRC_BART_res_sigma
          }
        }, error = function(e) {})
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
    list(name = "AFTv1", file_prefix = "AFTv1"),
    list(name = "AFTv2", file_prefix = "AFTv2"),
    list(name = "AFTv3", file_prefix = "AFTv3"),
    list(name = "BARTv1", file_prefix = "BARTv1"),
    list(name = "BARTv2", file_prefix = "BARTv2"),
    list(name = "BARTv3", file_prefix = "BARTv3")
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
          null_file <- paste0(mainDir,"lrcbart-sim-survival/res/", method_info$file_prefix, "_p", p_obs, null_file_suffix)
          if (file.exists(null_file)) {
            null_res <- read_plot_result(null_file)

            # For BART, filter by alpha and beta
            if (method_info$name %in% c("BARTv1", "BARTv2", "BARTv3")) {
              null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            }

            # Calculate mean FP (Type I error)
            if ("fp" %in% colnames(null_res)) {
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

      # Try to load hierAFT null results for each prior value
      for (prior_val in prior_vals) {
        tryCatch({
          aftv4_null_suffix <- paste0("_sc", sc,
            if (configuration == "Base") "" else configuration,
            "_cor", cor_val, "_prior", prior_val, "_null.RData")
          null_file <- paste0(mainDir,"lrcbart-sim-survival/res/hierAFT_p", p_obs, aftv4_null_suffix)
          if (file.exists(null_file)) {
            null_res <- read_plot_result(null_file)
            if ("fp" %in% colnames(null_res)) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                Configuration = configuration_label,
                c = cor_val,
                Method = paste0("hierAFT(", prior_val, ")"),
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {})
      }

      # Try to load LRC-BART null results for each N_target value
      for (lrc_index in seq_len(nrow(lrcbart_configurations))) {
      target_N <- lrcbart_configurations$target[lrc_index]
      lrc_suffix <- lrcbart_configurations$suffix[lrc_index]
      lrc_label <- lrcbart_configurations$label[lrc_index]

        tryCatch({
          lrcbart_null_suffix <- paste0("_sc", sc,
            if (configuration == "Base") "" else configuration,
            "_cor", cor_val, lrc_suffix,"_N", target_N, "_null.RData")
          null_file <- paste0(mainDir,"lrcbart-sim-survival/res/LRC-BART_p", p_obs, lrcbart_null_suffix)
          if (file.exists(null_file)) {
            null_res <- read_plot_result(null_file)
            null_res <- null_res[null_res$alpha == 0.95 & null_res$beta == 2, ]
            if ("fp" %in% colnames(null_res)) {
              mean_FP <- mean(null_res$fp, na.rm = TRUE)
              all_null_FP <- rbind(all_null_FP, data.frame(
                Configuration = configuration_label,
                c = cor_val,
                Method = lrc_label,
                Type_I_error = sprintf("%.2f", mean_FP)
              ))
            }
          }
        }, error = function(e) {})
      }
    }
  }
}

# Check if we have any data to plot
if (nrow(res) == 0) {
  stop("No results data found. Please check that result files exist and paths are correct.")
}

# LRC-BART MODIFICATION START
# Create labels for all migrated scenarios.
if (sc %in% 1:5) {
  if (nrow(res) > 0 && "Configuration" %in% colnames(res)) {
    res$Configuration_Label <- plot_factor(
      res$Configuration,
      levels = unique(res$Configuration)
    )
  }
  if (nrow(res_sigma) > 0 && "Configuration" %in% colnames(res_sigma)) {
    res_sigma$Configuration_Label <- plot_factor(
      res_sigma$Configuration,
      levels = unique(res_sigma$Configuration)
    )
  }
}
# LRC-BART MODIFICATION END

# Create dynamic method levels including hierAFT and LRC-BART with prior values
# hierAFT still labels by prior value; LRC-BART labels by N_target.
aftv4_methods <- paste0("hierAFT(",     aftv4_prior_vals, ")")
lrcbart_methods <- lrcbart_configurations$label
all_method_levels <- c("AFTv1", "AFTv2", "AFTv3", "BARTv1", "BARTv2", "BARTv3", aftv4_methods, lrcbart_methods)

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
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

  # Variance estimation plots (only if we have sigma data)
  if (nrow(res_sigma) > 0) {
    res_sigma_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
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
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # Check if subject-wise columns exist
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
      # Subject-wise bias and RMSE plots
      res_subj_long <- res_sigma %>%
        filter(Method %in% all_method_levels) %>%
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
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
          )
      } else {
        facet_plot_subj <- NULL
        table_grob_subj <- NULL
      }
    } else {
      facet_plot_subj <- NULL
      table_grob_subj <- NULL
    }

    # Population median plots
    res_pop_median_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
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

    summary_table_pop_median <- res_pop_median_long %>%
      group_by(Configuration_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean)

    # Split into treatment and control tables
    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_pop_median_trt_sep <- add_configuration_separators(summary_table_pop_median_trt, "Configuration_Label")
    summary_table_pop_median_ctrl_sep <- add_configuration_separators(summary_table_pop_median_ctrl, "Configuration_Label")
    table_grob_pop_median_trt <- pipeline_table(summary_table_pop_median_trt_sep, rows = NULL)
    table_grob_pop_median_ctrl <- pipeline_table(summary_table_pop_median_ctrl_sep, rows = NULL)
    table_grob_pop_median <- wrap_plots(table_grob_pop_median_trt, table_grob_pop_median_ctrl, ncol = 2)

    facet_plot_pop_median <- ggplot(res_pop_median_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Median Performance",
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
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # LRC-BART MODIFICATION START
    # Size plot and table panels from their actual configuration/row counts.
    n_configuration_rows <- length(unique(res$Configuration_Label))
    ate_plot_height <- max(6, 3 * n_configuration_rows)
    ate_table_height <- max(3, 0.30 * (nrow(summary_table_sep) + 1))
    pop_plot_height <- max(7, 5 * n_configuration_rows)
    pop_table_height <- max(
      3, 0.30 * (max(nrow(summary_table_pop_median_trt_sep),
                      nrow(summary_table_pop_median_ctrl_sep)) + 1)
    )
    sigma_plot_height <- max(7, 5 * n_configuration_rows)
    sigma_table_height <- max(
      3, 0.30 * (max(nrow(summary_table_trt_sep),
                      nrow(summary_table_ctrl_sep)) + 1)
    )

    # Combine ATE, population median, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      subj_plot_height <- max(6, 3 * n_configuration_rows)
      subj_table_height <- max(3, 0.30 * (nrow(summary_table_subj_sep) + 1))
      final_heights <- c(
        ate_plot_height, ate_table_height, pop_plot_height, pop_table_height,
        subj_plot_height, subj_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    } else {
      final_heights <- c(
        ate_plot_height, ate_table_height, pop_plot_height, pop_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_sigma / table_grob_sigma +
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

  # ---- Append the RMST(tau)-ratio column (right) beside the median-ratio column (left) ----
  rmst_col <- build_rmst_column(sc)
  if (!is.null(rmst_col)) {
    final_plot <- final_plot | rmst_col
  }

  # LRC-BART MODIFICATION START
  # Width follows the number of methods actually present. Table panels above
  # follow their actual row counts, preventing table/figure overlap.
  n_reporting_methods <- length(unique(as.character(
    res$Method[!is.na(res$Method)]
  )))
  n_reporting_methods <- max(1, n_reporting_methods)
  total_height <- max(sum(final_heights), attr(rmst_col, "required_height") %||% 0) + 2
  total_width <- max(14, 1.4 * n_reporting_methods)
  if (!is.null(rmst_col)) total_width <- 2 * total_width
  # LRC-BART MODIFICATION END

  save_pipeline_plot(file.path(projectDir, "inserts",
                   paste0("p", p_obs, "_sc", sc, "_all_results.jpg")),
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
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )

  # Variance estimation plots (only if we have sigma data)
  if (nrow(res_sigma) > 0) {
    res_sigma_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
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
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # Subject-wise plots
    if ("pehe_subj" %in% colnames(res_sigma) && "bias_subj" %in% colnames(res_sigma)) {
      res_subj_long <- res_sigma %>%
        filter(Method %in% all_method_levels) %>%
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
            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
          )
      } else {
        facet_plot_subj <- NULL
        table_grob_subj <- NULL
      }
    } else {
      facet_plot_subj <- NULL
      table_grob_subj <- NULL
    }

    # Population median plots
    res_pop_median_long <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
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

    summary_table_pop_median <- res_pop_median_long %>%
      group_by(Config_Cor_Label, c, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      rename(Correlation = c) %>%
      # Order by configuration first, then correlation
      arrange(Config_Cor_Label, Method)

    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between configurations and create tables
    summary_table_pop_median_trt_sep <- add_configuration_separators(summary_table_pop_median_trt, "Config_Cor_Label")
    summary_table_pop_median_ctrl_sep <- add_configuration_separators(summary_table_pop_median_ctrl, "Config_Cor_Label")
    table_grob_pop_median_trt <- pipeline_table(summary_table_pop_median_trt_sep, rows = NULL)
    table_grob_pop_median_ctrl <- pipeline_table(summary_table_pop_median_ctrl_sep, rows = NULL)
    table_grob_pop_median <- wrap_plots(table_grob_pop_median_trt, table_grob_pop_median_ctrl, ncol = 2)

    facet_plot_pop_median <- ggplot(res_pop_median_long, aes(x = Method, y = value, fill = Method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
      geom_boxplot(alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black", color = "black") +
      scale_fill_manual(values = method_colors) +
      labs(
        title = "Population Median Performance Across Configurations and Correlations",
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
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)
      )

    # LRC-BART MODIFICATION START
    # Size plot and table panels from their actual correlation/row counts.
    n_configuration_rows <- length(unique(res$Config_Cor_Label))
    ate_plot_height <- max(6, 3 * n_configuration_rows)
    ate_table_height <- max(3, 0.30 * (nrow(summary_table_sep) + 1))
    pop_plot_height <- max(7, 5 * n_configuration_rows)
    pop_table_height <- max(
      3, 0.30 * (max(nrow(summary_table_pop_median_trt_sep),
                      nrow(summary_table_pop_median_ctrl_sep)) + 1)
    )
    sigma_plot_height <- max(7, 5 * n_configuration_rows)
    sigma_table_height <- max(
      3, 0.30 * (max(nrow(summary_table_trt_sep),
                      nrow(summary_table_ctrl_sep)) + 1)
    )

    # Combine ATE, population median, subject-wise, and variance plots with tables
    if (!is.null(facet_plot_subj)) {
      subj_plot_height <- max(6, 3 * n_configuration_rows)
      subj_table_height <- max(3, 0.30 * (nrow(summary_table_subj_sep) + 1))
      final_heights <- c(
        ate_plot_height, ate_table_height, pop_plot_height, pop_table_height,
        subj_plot_height, subj_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_subj / table_grob_subj / facet_plot_sigma / table_grob_sigma +
        plot_layout(heights = final_heights)
    } else {
      final_heights <- c(
        ate_plot_height, ate_table_height, pop_plot_height, pop_table_height,
        sigma_plot_height, sigma_table_height
      )
      final_plot <- facet_plot / table_grob / facet_plot_pop_median / table_grob_pop_median / facet_plot_sigma / table_grob_sigma +
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

  # ---- Append the RMST(tau)-ratio column (right) beside the median-ratio column (left) ----
  rmst_col <- build_rmst_column(sc)
  if (!is.null(rmst_col)) {
    final_plot <- final_plot | rmst_col
  }

  n_reporting_methods <- length(unique(as.character(
    res$Method[!is.na(res$Method)]
  )))
  n_reporting_methods <- max(1, n_reporting_methods)
  total_height <- max(sum(final_heights), attr(rmst_col, "required_height") %||% 0) + 2
  total_width <- max(14, 1.4 * n_reporting_methods)
  if (!is.null(rmst_col)) total_width <- 2 * total_width

  save_pipeline_plot(file.path(projectDir, "inserts",
                   paste0("p", p_obs, "_sc", sc, "_all_results.jpg")),
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
