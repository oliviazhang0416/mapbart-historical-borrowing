# LRC-BART MODIFICATION START
# Plot one completed replicate-specific Gaussian single-arm ESS calibration.
# This script reads a checkpoint and performs no model fitting.

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
if (is.null(cal$ess_tau0_blocks) || is.null(cal$grid))
  stop("Checkpoint lacks ESS_tau0 draws or the s0^2 grid: ", cal_path)

curve_data <- data.frame(
  s0_sq = rep(cal$grid, each = nrow(cal$ess_tau0_blocks)),
  block = factor(rep(seq_len(nrow(cal$ess_tau0_blocks)),
                     times = length(cal$grid))),
  ESS_tau0 = as.vector(cal$ess_tau0_blocks)
)
mean_data <- data.frame(s0_sq = cal$grid, ESS_tau0 = cal$ess_tau0_grid)
selected <- cal$targets[, c("target_name", "target", "s0_sq", "ess_tau0")]
selected_targets <- selected[is.finite(selected$target), , drop = FALSE]

plot <- ggplot2::ggplot(
  curve_data,
  ggplot2::aes(x = s0_sq, y = ESS_tau0, group = block)
) +
  ggplot2::geom_line(alpha = 0.18, colour = "#7F8C8D") +
  ggplot2::geom_line(
    data = mean_data,
    ggplot2::aes(x = s0_sq, y = ESS_tau0, group = 1),
    linewidth = 1, colour = "#0072B2"
  ) +
  ggplot2::geom_hline(
    data = unique(selected_targets[c("target_name", "target")]),
    ggplot2::aes(yintercept = target, colour = target_name),
    linetype = "dashed", linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = selected,
    ggplot2::aes(x = s0_sq, y = ess_tau0, colour = target_name),
    inherit.aes = FALSE, size = 2.2
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = paste0("Single-arm LRC-BART ESS calibration: ",
                   tools::file_path_sans_ext(basename(cal$data_file))),
    subtitle = paste0("H_g = ", cal$H_g, "; ceiling = ",
                      format(round(cal$ceiling, 1), nsmall = 1)),
    x = expression(s[0]^2), y = expression(ESS[tau[0]]),
    colour = "Target"
  ) +
  ggplot2::theme_bw(base_size = 11)

default_out_path <- file.path(projectDir, "inserts", sub("\\.RData$", ".png", basename(cal_path)))
if (identical(as.integer(cal$H_g), 5L))
  default_out_path <- sub("_Hg5\\.png$", ".png", default_out_path)
out_path <- get_arg("--out", default_out_path)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(out_path, plot, width = 8, height = 5.5, dpi = 300)
cat("Wrote", out_path, "\n")
# LRC-BART MODIFICATION END
