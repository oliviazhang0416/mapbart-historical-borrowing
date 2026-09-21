#!/usr/bin/env Rscript

# Resumable validation runner.  It writes only below validation/results (or
# validation/smoke when VALIDATION_MODE=smoke).  The full job ledger is frozen
# before any filters are applied.

source("09-method-exploration/validation/validation_core.R")

mode <- tolower(Sys.getenv("VALIDATION_MODE", "validation"))
if (!mode %in% c("validation", "smoke")) stop("VALIDATION_MODE must be validation or smoke")
n_burn <- as.integer(Sys.getenv("VALIDATION_N_BURN", as.character(N_BURN_DEFAULT)))
n_draw <- as.integer(Sys.getenv("VALIDATION_N_DRAW", as.character(N_DRAW_DEFAULT)))
workers <- min(6L, as.integer(Sys.getenv("VALIDATION_WORKERS", "3")))
stopifnot(n_burn > 0L, n_draw > 0L, workers >= 1L, workers <= 6L)
if (mode == "validation" && (n_burn != 30000L || n_draw != 12000L))
  stop("Frozen validation requires 30000 burn-in and 12000 retained draws; short runs require smoke mode")
mode_dir <- file.path(VALIDATION_ROOT, if (mode == "smoke") "smoke" else "results")
dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

batch_started_path <- file.path(mode_dir, "batch_started_utc.txt")
if (file.exists(batch_started_path)) {
  batch_started_utc <- trimws(readLines(batch_started_path, warn = FALSE)[1])
} else {
  batch_started_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  tmp_batch <- paste0(batch_started_path, ".tmp-", Sys.getpid())
  writeLines(batch_started_utc, tmp_batch)
  stopifnot(file.rename(tmp_batch, batch_started_path))
}
stopifnot(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", batch_started_utc))
preselection_utc <- Sys.getenv("VALIDATION_PRESELECTION_UTC", "")
if (mode == "validation" && !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$", preselection_utc))
  stop("validation requires VALIDATION_PRESELECTION_UTC in whole-second UTC format")
if (mode == "smoke") preselection_utc <- batch_started_utc
primary_candidate_id <- Sys.getenv("VALIDATION_PRIMARY_CANDIDATE_ID", "JOINT-LRC-01")
primary_baseline_id <- Sys.getenv("VALIDATION_PRIMARY_BASELINE_ID", "SEPARATED-SOURCE-01")

if (mode == "validation") {
  frozen_state <- jsonlite::fromJSON("09-method-exploration/state.json")
  stopifnot(identical(preselection_utc, frozen_state$candidate_frozen_at_utc),
            identical(primary_candidate_id, frozen_state$validation_primary_candidate),
            identical(primary_baseline_id, frozen_state$validation_baseline_id),
            preselection_utc < batch_started_utc)
}

jobs <- make_job_table()
assert_same_manifest(file.path(VALIDATION_ROOT, "jobs.csv"), jobs)
writeLines(c(
  "Frozen validation queue; this file is generated before filtering.",
  paste0("mode=", mode, "; output=", normalizePath(mode_dir, mustWork = FALSE)),
  paste0("n_burn=", n_burn, "; n_draw=", n_draw, "; workers=", workers),
  paste0("batch_started_utc=", batch_started_utc),
  paste0("preselection_utc=", preselection_utc),
  paste0("primary_candidate_id=", primary_candidate_id, "; primary_baseline_id=", primary_baseline_id),
  "methods: JOINT-LRC-01, JOINT-SOURCE-01, SEPARATED-SOURCE-01",
  "conditions: compatible, one (ONE50 d2), two (XOR50 d2), heterogeneous (XOR50 d2), null (XOR50 d2)",
  "constant seeds 92026001:92026040; heterogeneous 92126001:92126040; null 92226001:92226040",
  "calibration: original training-only routine, H_f=50, H_g=5, N_target=100, 1000 burn/1000 draws",
  "joint treatment prior: tau ~ N(0,100); no ESS=100 interpretation is claimed for joint models",
  "This file records the numerical run; scientific review follows completion."
), file.path(mode_dir, "CONFIG.txt"))

csv_filter <- function(name, allowed) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) return(allowed)
  val <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  if (any(!val %in% allowed)) stop("invalid values in ", name)
  val
}
conditions <- csv_filter("VALIDATION_CONDITIONS", CONDITIONS)
methods <- csv_filter("VALIDATION_METHODS", METHODS)
chains <- as.integer(csv_filter("VALIDATION_CHAINS", as.character(1:3)))
reps <- as.integer(csv_filter("VALIDATION_REPS", as.character(1:40)))
stopifnot(all(chains %in% 1:3), all(reps %in% 1:40))

selected <- jobs[jobs$condition %in% conditions & jobs$method %in% methods &
                 jobs$chain %in% chains & jobs$rep %in% reps, , drop = FALSE]
limit <- as.integer(Sys.getenv("VALIDATION_LIMIT", "0"))
if (limit > 0L) selected <- selected[seq_len(min(limit, nrow(selected))), , drop = FALSE]
if (nrow(selected) == 0L) stop("filter selected no jobs")
selected$dataset_seed_effective <- vapply(seq_len(nrow(selected)), function(i)
  dataset_seed(selected$condition[i], selected$rep[i]), integer(1))
selected$fit_seed_effective <- if (mode == "smoke")
  98000000L + 10L * selected$job_id else selected$fit_seed
write.csv(selected, file.path(mode_dir, "jobs_selected.csv"), row.names = FALSE)

# Data and calibration are prepared serially so that identical constant-source
# datasets cannot race into separate cache files.  Fits are the only parallel
# stage.  For a selected constant condition, all three paired source
# conditions for that replicate are built and checked.
ensure_data_for_jobs(mode_dir, selected)

raw_path <- function(j) file.path(mode_dir, "raw",
  sprintf("%s_r%02d_c%d_%s.rds", j$condition, j$rep, j$chain, j$method))

fit_job <- function(i) {
  j <- selected[i, , drop = FALSE]
  path <- raw_path(j)
  if (file.exists(path)) {
    old <- readRDS(path)
    cached <- load_dataset(mode_dir, j$condition, j$rep)
    fields <- c("method", "condition", "rep", "chain", "seed", "n_burn", "n_draw",
                "batch_started_utc", "preselection_utc", "config_hash",
                "variance_model", "interval_level", "primary_candidate_id",
                "primary_baseline_id")
    got <- old[fields]
    want <- list(method = j$method, condition = j$condition,
                 rep = as.integer(j$rep), chain = as.integer(j$chain),
                 seed = as.integer(j$fit_seed_effective), n_burn = n_burn, n_draw = n_draw,
                 batch_started_utc = batch_started_utc,
                 preselection_utc = preselection_utc,
                 config_hash = method_config_hash(j$method, n_burn, n_draw),
                 variance_model = "source-specific", interval_level = 0.95,
                 primary_candidate_id = primary_candidate_id,
                 primary_baseline_id = primary_baseline_id)
    if (!identical(as.character(got$method), want$method) ||
        !identical(as.character(got$condition), want$condition) ||
        !identical(as.integer(got$rep), want$rep) ||
        !identical(as.integer(got$chain), want$chain) ||
        !identical(as.integer(got$seed), want$seed) ||
        !identical(as.integer(got$n_burn), want$n_burn) ||
        !identical(as.integer(got$n_draw), want$n_draw) ||
        !identical(as.character(got$batch_started_utc), want$batch_started_utc) ||
        !identical(as.character(got$preselection_utc), want$preselection_utc) ||
        !identical(as.character(got$config_hash), want$config_hash) ||
        !identical(as.character(got$variance_model), want$variance_model) ||
        abs(as.numeric(got$interval_level) - want$interval_level) > 1e-12 ||
        !identical(as.character(got$primary_candidate_id), want$primary_candidate_id) ||
        !identical(as.character(got$primary_baseline_id), want$primary_baseline_id))
      stop("fit cache configuration mismatch: ", path)
    if (!identical(old$data_sha256, cached$data_sha256) ||
        !identical(old$code_hash, CODE_HASH) ||
        !identical(old$library_hash, LIBRARY_HASH))
      stop("fit cache data/code/library hash mismatch: ", path)
    return(TRUE)
  }
  cached <- load_dataset(mode_dir, j$condition, j$rep)
  payload <- fit_validation_job(cached$data, cached$calibration, j$method,
                                as.integer(j$fit_seed_effective), n_burn, n_draw)
  payload$data_hash <- cached$data_hash
  payload$data_sha256 <- cached$data_sha256
  payload$code_hash <- cached$code_hash
  payload$library_hash <- cached$library_hash
  payload$config_hash <- method_config_hash(j$method, n_burn, n_draw)
  payload$prior_config <- c(PRIOR_CONFIG, list(method = j$method,
                                               H_g = if (j$method == "JOINT-LRC-01") H_G_LRC else 0L))
  payload$source_variance_config <- "separate sigma1_sq and sigma2_sq; source codes preserved"
  payload$variance_model <- "source-specific"
  payload$interval_level <- 0.95
  payload$preselection_utc <- preselection_utc
  payload$primary_candidate_id <- primary_candidate_id
  payload$primary_baseline_id <- primary_baseline_id
  payload$condition <- j$condition; payload$rep <- as.integer(j$rep)
  payload$chain <- as.integer(j$chain); payload$seed <- as.integer(j$fit_seed_effective)
  payload$n_burn <- n_burn; payload$n_draw <- n_draw
  payload$batch_started_utc <- batch_started_utc
  payload$job_id <- as.integer(j$job_id)
  atomic_save_rds(payload, path)
  message("Saved job ", j$job_id, "/1800: ", j$condition, " rep=", j$rep,
          " chain=", j$chain, " ", j$method)
  TRUE
}

message("Validation runner prepared ", nrow(selected), " jobs in mode=", mode,
        "; workers=", workers)
ans <- parallel::mclapply(seq_len(nrow(selected)), fit_job, mc.cores = workers,
                           mc.preschedule = FALSE)
stopifnot(all(vapply(ans, isTRUE, logical(1))))
message("Validation jobs complete")
