# LRC-BART MODIFICATION START
# Minimal mean-curve companion to ess_plot.R. No model is refitted here.

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Package 'ggplot2' is required for ESS plotting")

scriptDir <- tryCatch({
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
  if (length(file_argument))
    dirname(normalizePath(sub("^--file=", "", file_argument[1])))
  else dirname(normalizePath(sys.frames()[[1]]$ofile))
}, error = function(error) getwd())

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) default else args[index + 1L]
}

projectDir <- dirname(scriptDir)
resDir <- get_arg("--resdir", file.path(projectDir, "res", "ess"))
scenario <- get_arg("--scenario", NA_character_)
cal_path <- get_arg("--res", NA_character_)
if (is.na(cal_path)) {
  candidates <- list.files(
    resDir, pattern = "^ess_data_.*_Hg[0-9]+\\.RData$", full.names = TRUE
  )
  if (!is.na(scenario) && nzchar(scenario))
    candidates <- candidates[grepl(paste0("_", scenario, "_"),
                                    basename(candidates), fixed = TRUE)]
  if (!length(candidates))
    stop("No LRC-BART ESS checkpoint found in ", resDir)
  cal_path <- candidates[which.max(file.info(candidates)$mtime)]
}

cal <- readRDS(cal_path)

mean_data <- data.frame(s0_sq = cal$grid, ESS_tau0 = cal$ess_tau0_grid)
selected <- cal$targets[, c("target_name", "s0_sq", "ess_tau0")]
plot <- ggplot2::ggplot(mean_data,
                        ggplot2::aes(x = s0_sq, y = ESS_tau0)) +
  ggplot2::geom_line(linewidth = 0.9, colour = "#0072B2") +
  ggplot2::geom_point(
    data = selected,
    ggplot2::aes(x = s0_sq, y = ess_tau0, shape = target_name),
    inherit.aes = FALSE, size = 2.2
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = expression(s[0]^2), y = expression(ESS[tau[0]]),
                shape = "Target") +
  ggplot2::theme_classic(base_size = 11)

default_out_path <- file.path(projectDir, "inserts", sub("\\.RData$", "_simple.pdf", basename(cal_path)))
if (identical(as.integer(cal$H_g), 5L))
  default_out_path <- sub("_Hg5_simple\\.pdf$", "_simple.pdf",
                          default_out_path)
out_path <- get_arg("--out", default_out_path)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(out_path, plot, width = 6.5, height = 4.2)
cat("Wrote", out_path, "\n")
# LRC-BART MODIFICATION END
