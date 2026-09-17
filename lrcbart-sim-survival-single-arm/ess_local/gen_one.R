#!/usr/bin/env Rscript
# LRC-BART MODIFICATION START
# Generate only replicate 1 for a chosen survival single-arm scenario and
# treated-trial size. This optional utility mirrors run_all.R's source-text
# overrides and is useful for a bounded ESS check.

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

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) default else args[index + 1L]
}

sc_arg <- as.integer(get_arg("--sc", "1"))
n_T_arg <- as.integer(get_arg("--n_T", "30"))
seed_rct <- as.integer(get_arg("--seed_rct", "123"))
seed_rwd <- as.integer(get_arg("--seed_rwd", "456"))
stopifnot(sc_arg %in% c(1L, 2L), n_T_arg %in% c(30L, 200L))

data_gen <- paste0(
  .lrcRoot,
  "lrcbart-sim-survival-single-arm/data_gen_p10.R"
)
source_text <- paste(readLines(data_gen, warn = FALSE), collapse = "\n")
source_text <- sub("(?m)^sc <- [^\\n]*$", paste0("sc <- ", sc_arg),
                   source_text, perl = TRUE)
source_text <- sub("(?m)^n_T <- [^\\n]*$", paste0("n_T <- ", n_T_arg, "L"),
                   source_text, perl = TRUE)
source_text <- sub("(?m)^seed_rwd <- [^\\n]*$",
                   paste0("seed_rwd <- ", seed_rwd, "L"),
                   source_text, perl = TRUE)
source_text <- gsub("set\\.seed\\(\\s*6L?\\s*\\)",
                    paste0("set.seed(", seed_rct, "L)"),
                    source_text, perl = TRUE)
source_text <- sub("(?m)^n_replicates <- [^\\n]*$", "n_replicates <- 1L",
                   source_text, perl = TRUE)
eval(parse(text = source_text), envir = new.env(parent = globalenv()))
# LRC-BART MODIFICATION END
