# Master script: Gaussian single-arm historical borrowing pipeline.
#
#   1. Generate data with data_gen_p10.R
#   2. Plot generated-data balance with plot_balance.R and plot_balance_simple.R
#   3. Run LM-CP, BART-CP, HierLM, and LRC-BART
#   4. Plot replicate-specific ESS calibration
#   5. Plot simulation results
#
# Each analysis script owns its replicate loop. This runner changes header
# assignments in memory and evaluates each script in an isolated environment;
# it never rewrites the source files.

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
script_dir <- file.path(mainDir, "lrcbart-sim-gaussian-single-arm")
# LRC-BART MODIFICATION END

# ---- Pipeline flags ----
# All pipeline stages are on by default.
run_data_gen     <- TRUE
run_plot_balance <- TRUE
run_analysis     <- TRUE
run_ess_plot     <- TRUE
run_plot         <- TRUE

# LRC-BART ADDITION START
# The revision reports the fixed-w=0.9 and near-complete-pooling fits with the
# primary fixed-w=1 fits, so these supported single-arm sensitivities run by
# default after the default LRC-BART call.
run_lrc_sensitivities <- TRUE
sensitivity_configs <- c("w0.9", "s0min")
# LRC-BART ADDITION END

parallel_methods <- TRUE
use_parallel <- isTRUE(parallel_methods) && .Platform$OS.type == "unix"
n_cores <- NULL  # NULL => all parallel methods; otherwise cap each batch

# ---- Simulation design ----
p_obs <- 10L
# One folder-local setting controls data generation, analysis, and output names.
n_replicates <- 100L
hypo_values <- "alternative"

# LRC-BART MODIFICATION START
# Without RCT controls, Sc3's control discrepancy is not identifiable. The
# approved single-arm analysis therefore runs only Sc1 and Sc2.
data_configs <- list(list(sc = 1), list(sc = 2))
analysis_scenarios <- list(
  list(sc = 1, cor = 1),
  list(sc = 2, cor = 1)
)
# LRC-BART MODIFICATION END

# ---- Analysis files ----
files <- c("LMv2.R", "BARTv2.R", "hierLM.R", "lrcBART.R")

# Stan and standalone LRC-BART compilation run in the parent R process.
serial_files <- c("LMv2.R", "hierLM.R", "lrcBART.R")

seed_rct <- 123L
seed_rwd <- 456L
seed_method <- 789L

common_overrides <- list(
  seed_method = seed_method,
  n_replicates = n_replicates
)
method_overrides <- list(
  "LMv2.R" = list(rwd_w_vals = 1),
  "hierLM.R" = list(prior_vals = 0.05)
)
data_gen_overrides <- list(
  seed_rct = seed_rct,
  seed_rwd = seed_rwd,
  rwd_frozen = FALSE,
  p_obs = p_obs,
  n_replicates = n_replicates
)
# The runner owns execution using its configured method list and stage flags.
# ---- Execution helpers for this subproject ----
`%||%` <- function(x, default) if (is.null(x)) default else x

format_r_value <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

apply_overrides <- function(text, overrides) {
  seed_value <- overrides$seed_rct %||% overrides$seed_method
  rewrite <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr[[1]], as.name("set.seed")) && length(expr) >= 2L &&
        is.numeric(expr[[2]]) && !is.null(seed_value)) expr[[2]] <- seed_value
    if (identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]])) {
      name <- as.character(expr[[2]])
      if (name %in% names(overrides)) {
        expr[3] <- list(parse(text = format_r_value(overrides[[name]]))[[1]])
        return(expr)
      }
    }
    for (i in seq_along(expr)[-1L]) expr[i] <- list(rewrite(expr[[i]]))
    expr
  }
  paste(vapply(as.list(parse(text = text)), function(expr)
    paste(deparse(rewrite(expr), width.cutoff = 500L), collapse = "\n"), character(1)),
    collapse = "\n")
}

pipeline_result <- function(file, status = "OK", elapsed_sec = 0)
  data.frame(file = file, status = status, elapsed_sec = elapsed_sec,
             stringsAsFactors = FALSE)

run_one <- function(file, overrides = list()) {
  explicit_config <- overrides$lrc_config
  overrides <- modifyList(overrides, method_overrides[[file]] %||% list())
  # Explicit dispatch of a sensitivity must not collapse into a default override.
  if (file == "lrcBART.R" && !is.null(explicit_config)) overrides$lrc_config <- explicit_config
  start <- Sys.time()
  cat("Running:", file, "\n")
  status <- tryCatch({
    text <- paste(readLines(file.path(script_dir, file), warn = FALSE), collapse = "\n")
    text <- apply_overrides(text, overrides)
    eval(parse(text = text), new.env(parent = globalenv()))
    "OK"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  cat(sprintf(">>> %s: %s (%.1fs)\n", file, status, elapsed))
  pipeline_result(file, status, elapsed)
}

run_methods_parallel <- function(files, overrides, n_cores = NULL) {
  if (!length(files)) return(NULL)
  n_cores <- n_cores %||% length(files)
  if (length(n_cores) != 1L || !is.numeric(n_cores) || !is.finite(n_cores) ||
      n_cores < 1 || n_cores != floor(n_cores))
    stop("n_cores must be NULL or a positive integer")
  if (length(files) > n_cores) {
    batches <- split(files, ceiling(seq_along(files) / n_cores))
    return(do.call(rbind, lapply(batches, run_methods_parallel,
                                 overrides = overrides, n_cores = n_cores)))
  }
  status_dir <- tempfile("pipeline-workers-")
  dir.create(status_dir)
  on.exit(unlink(status_dir, recursive = TRUE), add = TRUE)
  worker <- function(file) {
    result <- tryCatch(run_one(file, overrides), error = function(e)
      pipeline_result(file, paste("FAILED:", conditionMessage(e)), NA_real_))
    marker <- file.path(status_dir, paste0(file, ".done"))
    writeLines(c(gsub("[\r\n]+", " ", result$status), format(result$elapsed_sec)),
               paste0(marker, ".tmp"))
    if (!file.rename(paste0(marker, ".tmp"), marker)) stop("Cannot publish worker status")
    invisible(NULL)
  }
  jobs <- setNames(lapply(files, function(file)
    parallel::mcparallel(worker(file), name = file)), files)
  pending <- files
  while (length(pending)) {
    Sys.sleep(0.2)
    for (file in pending) {
      marker <- file.path(status_dir, paste0(file, ".done"))
      if (file.exists(marker)) {
        pending <- setdiff(pending, file)
      } else {
        state <- suppressWarnings(system2("ps", c("-o", "state=", "-p",
                                      as.character(jobs[[file]]$pid)),
                                           stdout = TRUE, stderr = FALSE))
        state <- trimws(paste(state, collapse = ""))
        if (!nzchar(state) || startsWith(state, "Z")) pending <- setdiff(pending, file)
      }
    }
  }
  results <- lapply(files, function(file) {
    marker <- file.path(status_dir, paste0(file, ".done"))
    if (!file.exists(marker)) return(pipeline_result(file, "FAILED: worker died", NA_real_))
    lines <- tryCatch(readLines(marker, warn = FALSE), error = function(e) character())
    if (length(lines) != 2L) return(pipeline_result(file, "FAILED: malformed worker status", NA_real_))
    pipeline_result(file, lines[1], suppressWarnings(as.numeric(lines[2])))
  })
  suppressWarnings(try(parallel::mccollect(jobs, wait = FALSE), silent = TRUE))
  do.call(rbind, results)
}

finish_pipeline <- function(log) {
  results <- if (length(log)) do.call(rbind, log) else data.frame()
  cat("\nPipeline summary\n")
  print(results, row.names = FALSE)
  if (nrow(results) && any(grepl("^(FAILED|BLOCKED):", results$status)))
    stop("Pipeline finished with failed or blocked work; see the summary above.", call. = FALSE)
  invisible(results)
}

scenario_id <- function(cfg) {
  tag <- paste0("sc", cfg$sc)
  if (cfg$sc == 3) tag <- paste0(tag, "_cor", cfg$cor)
  if (cfg$sc == 4) tag <- paste0(tag, "_d", cfg$delta_rwd)
  if (cfg$sc == 5) tag <- paste0(tag, "_", cfg$region, "_d", cfg$delta_rwd)
  tag
}

data_key <- function(cfg, hypothesis, size = NA_integer_)
  paste(size, hypothesis, cfg$sc, cfg$delta_rwd %||% 0, cfg$region %||% "none", sep = "|")

# ---- Pipeline stages for this subproject ----
# Called with the enclosing runner's configuration and stage flags.
run_simulation_pipeline <- function() {
  log <- list()
  failed_data <- character()
  failed_results <- list()
  record <- function(stage, result, detail = "") {
    if (is.null(result)) return(invisible(NULL))
    result$stage <- stage
    result$detail <- detail
    log[[length(log) + 1L]] <<- result
    invisible(result)
  }
  blocked <- function(stage, file, detail)
    record(stage, pipeline_result(file, "BLOCKED: data generation failed"), detail)
  runner_env <- environment(run_simulation_pipeline)
  sizes <- get0("n_T_values", runner_env, inherits = FALSE, ifnotfound = NA_integer_)
  lrc_configs <- c("default", if (run_lrc_sensitivities) sensitivity_configs)
  # Method-specific numerical settings remain defined once in each runner.
  plot_methods <- method_overrides
  plot_context_for <- function(hypothesis, size, scenarios) list(
    hypothesis = hypothesis, scenarios = scenarios, n_replicates = n_replicates,
    n_T = size, methods = plot_methods, lrc_configs = lrc_configs,
    failed_results = failed_results
  )
  allowed_scenarios <- function(hypothesis, size) Filter(function(cfg)
    !data_key(cfg, hypothesis, size) %in% failed_data, analysis_scenarios)
  plot_calls <- function(scripts, stage) {
    for (size in sizes) for (hypothesis in hypo_values) {
      scope <- allowed_scenarios(hypothesis, size)
      detail <- paste("n_T=", size, hypothesis)
      for (script in scripts) {
        if (!length(scope)) {
          blocked(stage, script, detail)
          next
        }
        overrides <- list(p_obs = p_obs, n_replicates = n_replicates,
                          iter_rct = n_replicates, iter_rwd = n_replicates,
                          plot_context = plot_context_for(hypothesis, size, scope))
        if (!is.na(size)) overrides$n_T <- size
        record(stage, run_one(script, overrides), detail)
      }
    }
  }
  if (run_data_gen) {
    for (size in sizes) for (hypothesis in hypo_values) for (cfg in data_configs) {
      overrides <- modifyList(common_overrides, data_gen_overrides)
      overrides <- modifyList(overrides, c(list(hypo = hypothesis), cfg))
      if (!is.na(size)) overrides$n_T <- size
      result <- run_one(paste0("data_gen_p", p_obs, ".R"), overrides)
      key <- data_key(cfg, hypothesis, size)
      record("data", result, key)
      if (any(result$status != "OK")) failed_data <- c(failed_data, key)
    }
  }
  if (run_plot_balance) plot_calls(c("plot_balance.R", "plot_balance_simple.R"), "balance")

  if (run_analysis) {
    for (size in sizes) for (hypothesis in hypo_values) for (cfg in analysis_scenarios) {
      key <- data_key(cfg, hypothesis, size)
      detail <- paste(size, hypothesis, scenario_id(cfg))
      if (key %in% failed_data) {
        for (file in files) blocked("analysis", file, detail)
        next
      }
      overrides <- modifyList(common_overrides, c(list(hypo = hypothesis), cfg))
      if (!is.na(size)) overrides$n_T <- size
      parallel_files <- if (use_parallel) setdiff(files, serial_files) else character()
      sequential_files <- if (use_parallel) intersect(files, serial_files) else files
      result <- rbind(
        if (length(parallel_files)) run_methods_parallel(parallel_files, overrides, n_cores),
        if (length(sequential_files)) do.call(rbind, lapply(sequential_files, run_one, overrides = overrides))
      )
      record("analysis", result, detail)
      mark_failed <- function(rows, config = "default") {
        for (file in rows$file[rows$status != "OK"])
          failed_results[[length(failed_results) + 1L]] <<- list(
            method = sub("\\.R$", "", file), scenario = scenario_id(cfg),
            hypothesis = hypothesis, n_T = size, config = config)
      }
      mark_failed(result)
      if (run_lrc_sensitivities && "lrcBART.R" %in% files) {
        primary_ok <- any(result$file == "lrcBART.R" & result$status == "OK")
        for (config in sensitivity_configs) {
          if (primary_ok) {
            sensitivity <- run_one("lrcBART.R", modifyList(overrides, list(lrc_config = config)))
          } else {
            sensitivity <- pipeline_result("lrcBART.R", "BLOCKED: primary LRC-BART fit failed")
          }
          record("sensitivity", sensitivity, paste(detail, config))
          mark_failed(sensitivity, config)
        }
      }
    }
  }
  if (run_ess_plot) {
    for (size in sizes) for (hypothesis in hypo_values) for (cfg in analysis_scenarios) {
      detail <- paste(size, hypothesis, scenario_id(cfg))
      if (data_key(cfg, hypothesis, size) %in% failed_data) {
        blocked("ESS", "ess_local/ess_plot.R", detail)
        next
      }
      # A failed primary fit may have left an older checkpoint at this path.
      failed_primary <- any(vapply(failed_results, function(x)
        x$method == "lrcBART" && x$config == "default" && x$scenario == scenario_id(cfg) &&
        x$hypothesis == hypothesis && identical(x$n_T, size), logical(1)))
      if (failed_primary) {
        record("ESS", pipeline_result("ess_local/ess_plot.R", "BLOCKED: primary LRC-BART fit failed"), detail)
        next
      }
      res_dir <- file.path(script_dir, "res", "ess")
      size_tag <- if (is.na(size)) "" else paste0("_n", size)
      checkpoint <- file.path(res_dir, paste0("ess_data_p", p_obs, size_tag, "_",
        scenario_id(cfg), "_", hypothesis, "_", n_replicates, "_Hg5.RData"))
      if (!file.exists(checkpoint)) {
        record("ESS", pipeline_result("ess_local/ess_plot.R", "SKIPPED: completed checkpoint missing"), detail)
        next
      }
      overrides <- list(projectDir = script_dir, resDir = res_dir,
                        scenario = scenario_id(cfg), cal_path = checkpoint)
      if (!is.na(size)) overrides$n_T <- size
      record("ESS", run_one("ess_local/ess_plot.R", overrides), detail)
    }
  }
  if (run_plot) plot_calls("plot.R", "results")
  finish_pipeline(log)
}
# ---- End pipeline stages ----
stopifnot(length(n_replicates) == 1L, is.finite(n_replicates),
          n_replicates >= 1L, n_replicates == floor(n_replicates),
          length(hypo_values) > 0L, all(hypo_values %in% c("null", "alternative")))
pipeline_results <- run_simulation_pipeline()
