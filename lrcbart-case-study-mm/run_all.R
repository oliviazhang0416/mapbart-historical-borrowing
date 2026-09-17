# Master script: EloKRd plus UCMM myeloma application.
#
#   1. Clean and merge the private source data
#   2. Run both PFS and OS analyses
#   3. Construct the ESS plots after all requested LRC-BART fits finish
#   4. Build the result table

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
projectDir <- file.path(mainDir, "lrcbart-case-study-mm")

# All pipeline stages are on by default.
run_data_prep <- TRUE
run_analysis <- TRUE
run_ess_plot <- TRUE
run_summary <- TRUE

outcomes <- toupper(strsplit(
  Sys.getenv("OUTCOMES", unset = "PFS,OS"), "[, ]+"
)[[1]])
outcomes <- outcomes[nzchar(outcomes)]
stopifnot(length(outcomes) >= 1L, all(outcomes %in% c("PFS", "OS")))

merged_file <- Sys.getenv(
  "MERGED_FILE",
  unset = file.path(projectDir, "data_cleaned",
                     "merged_elokrd_ucmm_n283.RData")
)
if (!grepl("^/", merged_file))
  merged_file <- file.path(projectDir, "data_cleaned", basename(merged_file))
data_tag <- regmatches(
  basename(merged_file), regexpr("n[0-9]+", basename(merged_file))
)
if (!length(data_tag)) data_tag <- sub("\\.RData$", "", basename(merged_file))

method_files <- c("km.R", "AFTv2.R", "BARTv2.R", "hierAFT.R", "lrcBART.R")
script_override <- Sys.getenv("SCRIPTS", unset = "")
if (nzchar(script_override)) {
  method_files <- strsplit(script_override, "[, ]+")[[1]]
  method_files <- method_files[nzchar(method_files)]
}
if (!length(method_files)) stop("SCRIPTS must select at least one analysis method")
unknown_files <- setdiff(method_files,
                         c("km.R", "AFTv2.R", "BARTv2.R", "hierAFT.R",
                           "lrcBART.R"))
if (length(unknown_files))
  stop("Unknown SCRIPTS value(s): ", paste(unknown_files, collapse = ", "))

lrc_configs <- Sys.getenv(
  "LRC_CONFIGS", unset = "default,Hf50,w0.9,Hf50_w0.9"
)

# ---- Execution helpers for this subproject ----
pipeline_result <- function(file, status = "OK", elapsed_sec = 0)
  data.frame(file = file, status = status, elapsed_sec = elapsed_sec,
             stringsAsFactors = FALSE)

finish_pipeline <- function(log) {
  results <- if (length(log)) do.call(rbind, log) else data.frame()
  cat("\nPipeline summary\n")
  print(results, row.names = FALSE)
  if (nrow(results) && any(grepl("^(FAILED|BLOCKED):", results$status)))
    stop("Pipeline finished with failed or blocked work; see the summary above.", call. = FALSE)
  invisible(results)
}

quote_environment <- function(env) {
  if (!length(env)) return(character())
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*=", env))) stop("Invalid environment assignment")
  keys <- sub("=.*$", "", env)
  values <- substring(env, nchar(keys) + 2L)
  paste0(keys, "=", vapply(values, shQuote, character(1)))
}

# ---- Pipeline stages for this subproject ----
run_script <- function(path, env = character(), args = character()) {
  start <- Sys.time()
  cat("Running:", basename(path), "\n")
  status <- tryCatch({
    exit_code <- system2(file.path(R.home("bin"), "Rscript"),
                         c(shQuote(path), vapply(args, shQuote, character(1))),
                         env = quote_environment(env))
    if (identical(exit_code, 0L)) "OK" else paste("FAILED: exit code", exit_code)
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  cat(sprintf(">>> %s: %s (%.1fs)\n", basename(path), status, elapsed))
  pipeline_result(basename(path), status, elapsed)
}

run_case_study_pipeline <- function() {
  log <- list()
  failed_methods <- character()
  record <- function(stage, result, detail = "") {
    result$stage <- stage
    result$detail <- detail
    log[[length(log) + 1L]] <<- result
    result
  }
  blocked <- function(stage, file, detail = "data preparation failed")
    record(stage, pipeline_result(file, paste("BLOCKED:", detail)))
  ready <- TRUE
  if (run_data_prep) {
    cleaning <- lapply(c("data_cleaning_ucmm.R", "data_cleaning_elokrd.R"), function(file)
      record("data", run_script(file.path(projectDir, "private_data", file))))
    ready <- all(vapply(cleaning, function(x) x$status == "OK", logical(1)))
    if (ready) {
      merged <- record("data", run_script(file.path(projectDir, "private_data", "data_merge.R")))
      ready <- merged$status == "OK"
    } else blocked("data", "data_merge.R")
  }
  if (!file.exists(merged_file)) {
    record("data", pipeline_result("merged data", "FAILED: requested merged file is missing"))
    ready <- FALSE
  }
  configs <- unique(strsplit(lrc_configs, "[,[:space:]]+")[[1]])
  configs <- configs[nzchar(configs)]
  if (!length(configs) || any(!configs %in% c("default", "Hf50", "w0.9", "Hf50_w0.9")))
    stop("Invalid LRC_CONFIGS selection")
  if (run_analysis) for (outcome in outcomes) for (method in method_files) {
    if (!ready) {
      blocked("analysis", method)
      next
    }
    result <- record("analysis", run_script(file.path(projectDir, method), env = c(
      paste0("OUTCOME=", outcome), paste0("MERGED_FILE=", merged_file),
      "AFTV2_RWD_W=1", "HIERAFT_CONFIGS=0.05,0.5",
      paste0("LRC_CONFIGS=", paste(configs, collapse = ","))
    )), outcome)
    if (result$status != "OK") failed_methods <- c(failed_methods, paste(outcome, method, sep = ":"))
  }
  if (run_ess_plot) for (outcome in outcomes) {
    if (!ready || paste(outcome, "lrcBART.R", sep = ":") %in% failed_methods) {
      blocked("ESS", "ess_local/ess_plot.R", paste(outcome, "data preparation or LRC-BART fit failed"))
      next
    }
    hf <- unique(ifelse(grepl("Hf50", configs), 50L, 10L))
    paths <- file.path(projectDir, "res", "ess", paste0("ess_", tolower(outcome),
      "_", data_tag, "_Hf", hf, ".RData"))
    if (!any(file.exists(paths))) {
      record("ESS", pipeline_result("ess_plot.R", "SKIPPED: completed checkpoint missing"), outcome)
      next
    }
    record("ESS", run_script(file.path(projectDir, "ess_local", "ess_plot.R"),
      args = c("--outcome", outcome, "--data-tag", data_tag, "--hf", paste(hf, collapse = ","))), outcome)
  }
  if (run_summary) {
    if (!ready) {
      blocked("summary", "summarize.R")
    } else {
      record("summary", run_script(file.path(projectDir, "summarize.R"), env = c(
        paste0("MERGED_FILE=", merged_file), paste0("OUTCOMES=", paste(outcomes, collapse = ",")),
        paste0("LRC_CONFIGS=", paste(configs, collapse = ",")),
        paste0("FAILED_METHODS=", paste(failed_methods, collapse = ","))
      )))
    }
  }
  finish_pipeline(log)
}
# ---- End pipeline stages ----
pipeline_results <- run_case_study_pipeline()
