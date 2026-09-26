# JOINT LRC-BART MODIFICATION START
# Matching original survival single-arm plot, adapted to the agreed Sc1 descendants.
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
# Retain the inherited survival result layout while reading the migrated
# single-arm files for the requested trial size.
mainDir <- .lrcRoot
projectDir <- file.path(mainDir, "lrcbart-sim-survival-single-arm")

p_obs <- 10L
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

n_T <- 30L
stopifnot(n_T %in% c(30L, 200L))
size_suffix <- paste0("_n", n_T)
resultDir <- file.path(projectDir, "res", paste0("n", n_T))
insertDir <- file.path(projectDir, "inserts")
dir.create(insertDir, recursive = TRUE, showWarnings = FALSE)

aftv4_prior_vals <- plot_method_setting("hierAFT.R", "prior_vals", 0.05)
aftv2_rwd_w_vals <- plot_method_setting("AFTv2.R", "rwd_w_vals", 1)
n_rwd_nominal <- 300
lrcbart_configurations <- lrc_plot_configurations(TRUE)
prior_vals <- aftv4_prior_vals
# LRC-BART MODIFICATION END

# Helper function to add separator rows between configurations
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
# figure). Contains the SAME elements as the former standalone plot_rmst.R:
# RMST performance facet + summary table, arm-specific population RMST facet +
# tables, and (if present) subject-wise RMST facet + table. Reads the RMST
# columns fresh from res/*.RData. Returns a patchwork stack, or NULL if absent.
# =============================================================================
build_rmst_column <- function(sc) {
  RMST_COLS <- c("c","iteration","rmst_tau","rmst_true","rmst_hat",
                 "bias_rmst","sd_rmst","rmse_rmst","w1distance_rmst","w2distance_rmst",
                 "ci_rmst","coverage_rmst","tp_rmst","tp_calibrated_rmst",
                 "bias.trt.rmst.pop","sd.trt.rmst.pop","w2distance.trt.rmst.pop",
                 "bias.ctrl.rmst.pop","sd.ctrl.rmst.pop","w2distance.ctrl.rmst.pop",
                 "bias_subj_rmst","pehe_subj_rmst")
  resdir <- paste0(resultDir, "/")
  read_rmst <- function(file, method) {
    if (!file.exists(file)) return(NULL)
    d <- tryCatch(read_plot_result(file), error = function(e) NULL)
    if (is.null(d) || !"bias_rmst" %in% names(d)) return(NULL)
    out <- d[, intersect(RMST_COLS, names(d)), drop = FALSE]; out$Method <- method; out
  }
  build_one <- function(scenario_tag, configuration_label) {
    suf <- function(extra = "")
      paste0("_", scenario_tag, extra, paste0("_", plot_hypothesis, ".RData"))
    fr <- list()
    for (.aw in aftv2_rwd_w_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"AFTv2_p",p_obs,size_suffix,
                                        sub("\\.RData$", paste0("_w",.aw,".RData"), suf())),
                                 weight_labels("AFTv2", aftv2_rwd_w_vals)[match(.aw, aftv2_rwd_w_vals)])))
    fr <- c(fr, list(read_rmst(paste0(resdir,"BARTv2_p",p_obs,size_suffix,suf()),"BARTv2")))
    for (pr in aftv4_prior_vals)
      fr <- c(fr, list(read_rmst(paste0(resdir,"hierAFT_p",p_obs,size_suffix,suf(paste0("_prior",pr))), paste0("HierAFT(",pr,")"))))
    for (lrc_index in seq_len(nrow(lrcbart_configurations)))
      fr <- c(fr, list(read_rmst(
        paste0(resdir, "LRC-BART_p", p_obs, size_suffix, "_", scenario_tag,
               lrcbart_configurations$suffix[lrc_index], "_N",
               lrcbart_configurations$target[lrc_index],
               paste0("_", plot_hypothesis, ".RData")),
        lrcbart_configurations$label[lrc_index]
      )))
    fr <- fr[!sapply(fr, is.null)]; if (!length(fr)) return(NULL)
    out <- bind_rows(fr); out$Configuration <- configuration_label; out
  }
  configurations <- if (sc == "1a") c("d1", "d2") else
    if (sc == "1b") c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2") else
    if (sc %in% c("1c", "1c-i", "1c-ii")) "d2" else "Base"
  dat <- bind_rows(lapply(selected_configurations(sc, configurations), function(configuration) {
    build_one(paste0("sc", sc, if (sc %in% c("1", "2")) "" else paste0("_", configuration)),
      if (configuration == "Base") "Default" else configuration)
  }))
  if (is.null(dat) || !nrow(dat)) return(NULL)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  grp <- c("Configuration", "Method")
  # Main table: SAME columns as the median main summary table. RMSE = sqrt(mean(.))
  # like the median table; metrics fall back to "-" if absent in old result files.
  for (.cc in c("rmse_rmst","w1distance_rmst","w2distance_rmst","tp_calibrated_rmst",
                "sd.trt.rmst.pop","sd.ctrl.rmst.pop"))
    if (!.cc %in% names(dat)) dat[[.cc]] <- NA_real_
  .fmt  <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", m) }
  .fmtR <- function(x) { m <- mean(x, na.rm = TRUE); if (!is.finite(m)) "-" else sprintf("%.2f", sqrt(m)) }
  rmst_main <- dat %>% group_by(Configuration, Method) %>%
    summarise(Bias=.fmt(bias_rmst), SD=.fmt(sd_rmst), RMSE=.fmtR(rmse_rmst),
              W1Distance=.fmt(w1distance_rmst), W2Distance=.fmt(w2distance_rmst),
              CI_length=.fmt(ci_rmst), CI_coverage=.fmt(coverage_rmst),
              Power=.fmt(tp_rmst), Power_calib=.fmt(tp_calibrated_rmst),
              N=n(), N_bias_ge1=sum(abs(bias_rmst)>=1,na.rm=TRUE), .groups="drop") %>%
    mutate(Configuration_Label=factor(Configuration), .before=1)
  table_grob <- pipeline_table(
    add_configuration_separators(rmst_main, "Configuration_Label"),
    rows = NULL
  )
  long <- dat %>% pivot_longer(c(bias_rmst,sd_rmst,ci_rmst,coverage_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
  long$metric <- factor(long$metric, levels=c("bias_rmst","sd_rmst","ci_rmst","coverage_rmst"), labels=c("Bias","SD","CI length","Coverage"))
  thm <- theme_minimal() + theme(plot.title=element_text(size=14,face="bold",hjust=0.5), axis.text.x=element_text(size=9,angle=45,hjust=1), legend.position="none", strip.text=element_text(size=11,face="bold"), strip.background=element_rect(fill="gray90",color="gray50"), panel.border=element_rect(color="black",fill=NA,linewidth=0.5))
  p <- ggplot(long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("RMST(τ)-ratio performance — Scenario ",sc," (τ = ",sprintf("%.1f",mean(dat$rmst_tau,na.rm=TRUE)),")"), x=NULL, y=NULL) + thm
  p <- p + facet_grid(Configuration~metric,scales="free_y")
  arm_long <- dat %>% pivot_longer(c(bias.trt.rmst.pop,w2distance.trt.rmst.pop,bias.ctrl.rmst.pop,w2distance.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control")), Metric=ifelse(grepl("^bias",metric),"Bias","W2Distance"))
  # Arm tables use the same columns as the median population-median tables.
  arm_summary <- arm_long %>% group_by(Configuration, Method, Arm, Metric) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=Metric, values_from=m)
  # SD column: mean posterior SD of the per-arm estimate (matches the ATE table's
  # SD = mean(sd_rmst); uses per-replicate posterior SD cols sd.{trt,ctrl}.rmst.pop).
  arm_sd <- dat %>% pivot_longer(c(sd.trt.rmst.pop,sd.ctrl.rmst.pop), names_to="metric", values_to="value") %>% filter(!is.na(value)) %>%
    mutate(Arm=factor(ifelse(grepl("\\.trt\\.",metric),"Treatment","Control"), levels=c("Treatment","Control"))) %>%
    group_by(Configuration, Method, Arm) %>% summarise(SD=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop")
  arm_summary <- dplyr::left_join(arm_summary, arm_sd, by=c("Configuration","Method","Arm"))
  arm_tab_trt <- arm_summary %>% filter(Arm=="Treatment") %>% dplyr::select(-Arm) %>% mutate(Configuration_Label=factor(Configuration), Group="Treatment") %>% dplyr::select(Group, Configuration_Label, Method, Bias, SD, W2Distance)
  arm_tab_ctrl <- arm_summary %>% filter(Arm=="Control") %>% dplyr::select(-Arm) %>% mutate(Configuration_Label=factor(Configuration), Group="Control") %>% dplyr::select(Group, Configuration_Label, Method, Bias, SD, W2Distance)
  table_grob_arm <- wrap_plots(
    pipeline_table(add_configuration_separators(arm_tab_trt,"Configuration_Label"),rows=NULL),
    pipeline_table(add_configuration_separators(arm_tab_ctrl,"Configuration_Label"),rows=NULL),
    ncol=2
  )
  p_arm <- ggplot(arm_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
    geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
    labs(title=paste0("Arm-specific population RMST(τ) — Scenario ",sc,"  (estimate vs DGP truth, per arm)"), x=NULL,y=NULL) + thm
  p_arm <- p_arm + facet_grid(Configuration+Arm~Metric,scales="free_y")
  # Table-slot height grows with the number of methods so all rows show in full.
  .tbl_h <- max(2, length(unique(as.character(dat$Method))) * 0.16)
  has_subj <- all(c("bias_subj_rmst","pehe_subj_rmst") %in% names(dat))
  if (has_subj) {
    subj_long <- dat %>% pivot_longer(c(bias_subj_rmst,pehe_subj_rmst), names_to="metric", values_to="value") %>% filter(!is.na(value))
    subj_long$metric <- factor(subj_long$metric, levels=c("bias_subj_rmst","pehe_subj_rmst"), labels=c("Bias","PEHE"))
    # Subject-wise table uses the same columns as the median subject-wise table.
    subj_tab <- subj_long %>% group_by(Configuration, Method, metric) %>% summarise(m=sprintf("%.2f",mean(value,na.rm=TRUE)),.groups="drop") %>% pivot_wider(names_from=metric, values_from=m) %>%
      mutate(Configuration_Label=factor(Configuration)) %>% dplyr::select(Configuration_Label, Method, Bias, PEHE)
    table_grob_subj <- pipeline_table(add_configuration_separators(subj_tab,"Configuration_Label"), rows = NULL)
    p_subj <- ggplot(subj_long, aes(Method,value,fill=Method)) + geom_hline(yintercept=0,linetype="dashed",color="black",alpha=0.6) +
      geom_boxplot(alpha=0.7,outlier.size=0.6) + stat_summary(fun=mean,geom="point",shape=23,size=2.5,fill="black") +
      labs(title=paste0("Subject-wise RMST(τ)-ratio — Scenario ",sc,"  (per-subject estimate vs truth)"), x=NULL,y=NULL) + thm
    p_subj <- p_subj + facet_grid(Configuration~metric,scales="free_y")
    # Blank bottom section (matches the median column's variance/sigma section) so
    # the two columns have identical row layout and align horizontally.
    return(p / table_grob / p_arm / table_grob_arm / p_subj / table_grob_subj /
             patchwork::plot_spacer() / patchwork::plot_spacer() +
             plot_layout(heights = c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h, 4, .tbl_h)))
  }
  p / table_grob / p_arm / table_grob_arm /
    patchwork::plot_spacer() / patchwork::plot_spacer() +
    plot_layout(heights = c(4, .tbl_h, 4, .tbl_h, 4, .tbl_h))
}

# LRC-BART MODIFICATION START
# Read the supported survival single-arm comparators and LRC-BART
# configurations while retaining the inherited reporting sections below.
for (sc in selected_sc(c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2"))) {
tryCatch({

if (sc %in% c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2")) {

  all_res_ATE <- data.frame()
  all_res_sigma <- data.frame()
  configurations <- if (sc == "1a") c("d1", "d2") else
    if (sc == "1b") c("X5_d0.5", "X5_d1", "X5_d2", "X7_d1", "X7_d2") else
    if (sc %in% c("1c", "1c-i", "1c-ii")) "d2" else "Base"
  for (configuration in selected_configurations(sc, configurations)) {
  scenario_tag <- paste0("sc", sc, if (sc %in% c("1", "2")) "" else paste0("_", configuration))
  configuration_label <- if (configuration == "Base") "Default" else configuration

  method_files <- data.frame(
    method = c(
      weight_labels("AFTv2", aftv2_rwd_w_vals), "BARTv2", paste0("HierAFT(", prior_vals, ")"),
      lrcbart_configurations$label
    ),
    path = c(
      file.path(resultDir, paste0("AFTv2_p", p_obs, size_suffix, "_", scenario_tag,
                                  paste0("_", plot_hypothesis, "_w", aftv2_rwd_w_vals, ".RData"))),
      file.path(resultDir, paste0("BARTv2_p", p_obs, size_suffix, "_", scenario_tag,
                                  paste0("_", plot_hypothesis, ".RData"))),
      file.path(resultDir, paste0("hierAFT_p", p_obs, size_suffix, "_", scenario_tag,
                                  "_prior", prior_vals,
                                  paste0("_", plot_hypothesis, ".RData"))),
      file.path(resultDir, paste0(
        "LRC-BART_p", p_obs, size_suffix, "_", scenario_tag,
        lrcbart_configurations$suffix, "_N",
        lrcbart_configurations$target, paste0("_", plot_hypothesis, ".RData")
      ))
    ),
    stringsAsFactors = FALSE
  )

  ate_columns <- c(
    "c", "iteration", "bias", "sd", "rmse", "w1distance",
    "w2distance", "ci", "coverage", "tp", "fp", "tp_calibrated"
  )
  sigma_columns <- c(
    "c", "iteration", "bias.trt.sigma", "sd.trt.sigma",
    "w2distance.trt.sigma", "bias.ctrl.sigma", "sd.ctrl.sigma",
    "w2distance.ctrl.sigma", "bias.trt.median.pop",
    "sd.trt.median.pop", "w2distance.trt.median.pop",
    "bias.ctrl.median.pop", "sd.ctrl.median.pop",
    "w2distance.ctrl.median.pop", "pehe_subj", "bias_subj"
  )

  for (method_index in seq_len(nrow(method_files))) {
    result_file <- method_files$path[method_index]
    if (!file.exists(result_file)) next
    method_result <- tryCatch(read_plot_result(result_file), error = function(e) NULL)
    if (is.null(method_result)) next
    if (all(c("alpha", "beta") %in% names(method_result)))
      method_result <- method_result[
        method_result$alpha == 0.95 & method_result$beta == 2,
        , drop = FALSE
      ]

    if (all(ate_columns %in% names(method_result))) {
      method_ate <- method_result[, ate_columns, drop = FALSE]
      method_ate$Method <- method_files$method[method_index]
      method_ate$Configuration <- configuration_label
      all_res_ATE <- rbind(all_res_ATE, method_ate)
    }
    if (all(sigma_columns %in% names(method_result))) {
      method_sigma <- method_result[, sigma_columns, drop = FALSE]
      method_sigma$Method <- method_files$method[method_index]
      method_sigma$Configuration <- configuration_label
      all_res_sigma <- rbind(all_res_sigma, method_sigma)
    }
  }

  }
  res <- all_res_ATE
  res_sigma <- all_res_sigma
  all_null_FP <- data.frame()
}
# LRC-BART MODIFICATION END

# Check if we have any data to plot
if (nrow(res) == 0) {
  stop("No results data found. Please check that result files exist and paths are correct.")
}

# Create the single retained configuration label.
if (sc %in% c("1", "1a", "1b", "1c", "1c-i", "1c-ii", "2")) {
  if (nrow(res) > 0 && "Configuration" %in% colnames(res))
    res$Configuration_Label <- plot_factor(res$Configuration, levels = "Default")
  if (nrow(res_sigma) > 0 && "Configuration" %in% colnames(res_sigma))
    res_sigma$Configuration_Label <- plot_factor(res_sigma$Configuration,
                                             levels = "Default")
}

# Create dynamic method levels for the supported methods.
aftv2_methods <- weight_labels("AFTv2", aftv2_rwd_w_vals)
aftv4_methods <- paste0("HierAFT(", aftv4_prior_vals, ")")
lrcbart_methods <- lrcbart_configurations$label
all_method_levels <- c(aftv2_methods, "BARTv2", aftv4_methods,
                       lrcbart_methods)

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

  # Add separator rows between subscenarios and create table
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
        group = factor(group, levels = c("Treatment", "Control"))
      )

    summary_table_sigma <- res_sigma_long %>%
      group_by(Configuration_Label, Method, group, metric_type) %>%
      summarise(Mean = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop") %>%
      pivot_wider(names_from = metric_type, values_from = Mean) %>%
      # Reorder columns: Bias, SD, W2Distance
      dplyr::select(Configuration_Label, Method, group, Bias, SD, W2Distance)

    # Split into treatment and control tables
    summary_table_sigma_trt <- summary_table_sigma %>%
      filter(group == "Treatment") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Treatment", .before = 1)

    summary_table_sigma_ctrl <- summary_table_sigma %>%
      filter(group == "Control") %>%
      dplyr::select(-group) %>%
      mutate(Group = "Control", .before = 1)

    # Add separator rows between subscenarios and create tables
    summary_table_sigma_trt_sep <- add_configuration_separators(summary_table_sigma_trt, "Configuration_Label")
    summary_table_sigma_ctrl_sep <- add_configuration_separators(summary_table_sigma_ctrl, "Configuration_Label")
    table_grob_sigma_trt <- pipeline_table(summary_table_sigma_trt_sep, rows = NULL)
    table_grob_sigma_ctrl <- pipeline_table(summary_table_sigma_ctrl_sep, rows = NULL)
    table_grob_sigma <- wrap_plots(table_grob_sigma_trt, table_grob_sigma_ctrl, ncol = 2)

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

        # Add separator rows between subscenarios and create table
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

    # SD column: mean posterior SD of the per-arm estimate (matches the ATE table's
    # SD = mean(sd); uses the per-replicate posterior SD cols sd.{trt,ctrl}.median.pop).
    sd_pop_median <- res_sigma %>%
      filter(Method %in% all_method_levels) %>%
      pivot_longer(cols = c(sd.trt.median.pop, sd.ctrl.median.pop),
                   names_to = "metric", values_to = "value") %>%
      mutate(group = factor(ifelse(grepl("trt", metric), "Treatment", "Control"),
                            levels = c("Treatment", "Control"))) %>%
      group_by(Configuration_Label, Method, group) %>%
      summarise(SD = sprintf("%.2f", mean(value, na.rm = TRUE)), .groups = "drop")
    summary_table_pop_median <- dplyr::left_join(summary_table_pop_median, sd_pop_median,
      by = c("Configuration_Label", "Method", "group"))

    # Split into treatment and control tables (Bias, SD, W2Distance)
    summary_table_pop_median_trt <- summary_table_pop_median %>%
      filter(group == "Treatment") %>%
      mutate(Group = "Treatment") %>%
      dplyr::select(Group, Configuration_Label, Method, Bias, SD, W2Distance)

    summary_table_pop_median_ctrl <- summary_table_pop_median %>%
      filter(group == "Control") %>%
      mutate(Group = "Control") %>%
      dplyr::select(Group, Configuration_Label, Method, Bias, SD, W2Distance)

    # Add separator rows between subscenarios and create tables
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
    # Match the other migrated folders by sizing each table from its actual
    # row count and each plot from its number of configuration rows.
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
      3, 0.30 * (max(nrow(summary_table_sigma_trt_sep),
                      nrow(summary_table_sigma_ctrl_sep)) + 1)
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
  n_reporting_methods <- length(unique(as.character(
    res$Method[!is.na(res$Method)]
  )))
  n_reporting_methods <- max(1, n_reporting_methods)
  total_height <- sum(final_heights) + 2
  total_width <- max(14, 1.4 * n_reporting_methods)
  if (!is.null(rmst_col)) total_width <- 2 * total_width
  # LRC-BART MODIFICATION END

  save_pipeline_plot(file.path(insertDir, paste0("p", p_obs, size_suffix, "_sc", sc, "_all_results.jpg")),
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
