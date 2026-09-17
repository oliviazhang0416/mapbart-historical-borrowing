# LRC-BART ADDITION START
# Overlay the completed H_f=10 and H_f=50 ESS curves for one endpoint.
# No model fitting is performed here.

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
get_arg <- function(flag, default) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) default else args[index + 1L]
}

projectDir <- file.path(.lrcRoot, "lrcbart-case-study-mm")
outcome <- toupper(get_arg("--outcome", "PFS"))
data_tag <- get_arg("--data-tag", "n283")
hf_values <- as.integer(strsplit(get_arg("--hf", "10,50"), ",", fixed = TRUE)[[1]])
stopifnot(length(hf_values) > 0L, all(hf_values %in% c(10L, 50L)))
resDir <- get_arg("--resdir", file.path(projectDir, "res", "ess"))
out_path <- get_arg(
  "--out",
  file.path(projectDir, "inserts", paste0("ess_curve_", tolower(outcome), "_", data_tag,
                           ".png"))
)

paths <- file.path(
  resDir,
  paste0("ess_", tolower(outcome), "_", data_tag, "_Hf", unique(hf_values),
         ".RData")
)
available <- file.exists(paths)
if (!any(available))
  stop("No completed ESS checkpoint found for ", outcome, " ", data_tag)
paths <- paths[available]
calibrations <- lapply(paths, readRDS)

colors <- c(`10` = "#0072B2", `50` = "#E69F00")
all_x <- unlist(lapply(calibrations, `[[`, "grid"))
all_y <- unlist(lapply(calibrations, function(x) x$ess_tau0_blocks))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
png(out_path, width = 8, height = 5.5, units = "in", res = 300)
par(mar = c(4.5, 4.6, 2.6, 1))
plot(
  NA, xlim = range(all_x), ylim = range(all_y, finite = TRUE), log = "x",
  xlab = expression(s[0]^2), ylab = expression(ESS[tau[0]]),
  main = paste0("LRC-BART ESS curve: ", outcome)
)
for (calibration in calibrations) {
  curve_color <- colors[as.character(calibration$H_f)]
  for (block_id in seq_len(nrow(calibration$ess_tau0_blocks)))
    lines(calibration$grid, calibration$ess_tau0_blocks[block_id, ],
          col = grDevices::adjustcolor(curve_color, alpha.f = 0.10), lwd = 0.7)
}
for (calibration in calibrations) {
  curve_color <- colors[as.character(calibration$H_f)]
  lines(calibration$grid, calibration$ess_tau0_grid,
        col = curve_color, lwd = 2.3, lty = 1)
  points(calibration$targets$s0_sq, calibration$targets$ess_tau0,
         pch = 19, cex = 0.85, col = curve_color)
}
legend(
  "topright",
  legend = c(
    vapply(calibrations, function(x) paste0("H_f = ", x$H_f, " mean"),
           character(1)),
    "Monte Carlo blocks", "Selected targets"
  ),
  col = c(
    vapply(calibrations, function(x) colors[as.character(x$H_f)],
           character(1)),
    "gray65", "black"
  ),
  lty = c(rep(1, length(calibrations)), 1, NA),
  lwd = c(rep(2.3, length(calibrations)), 0.8, NA),
  pch = c(rep(NA, length(calibrations) + 1L), 19),
  bty = "n"
)
mtext(
  paste0("H_g = 5; ceilings: ",
         paste(vapply(calibrations, function(x)
           paste0("H_f=", x$H_f, " ", round(x$ceiling, 1)), character(1)),
           collapse = "; ")),
  side = 3, line = 0.4, cex = 0.78
)
dev.off()
cat("Wrote", out_path, "\n")
# LRC-BART ADDITION END
