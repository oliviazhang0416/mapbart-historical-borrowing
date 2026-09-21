# Patient-profile PFS prior ESS for the 30 EloKRd covariate profiles.
library(jsonlite)

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", arg))
out <- normalizePath(file.path(dirname(script), ".."))
root <- normalizePath(file.path(out, "../.."))
data_file <- file.path(
  root, "lrcbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n230.RData"
)
calibration_file <- file.path(
  root, "lrcbart-case-study-mm/res/ess/ess_pfs_n230_Hf10.RData"
)

e <- new.env()
load(data_file, e)
dat <- e$merged
cal <- readRDS(calibration_file)
trial <- dat[dat$trt == 1, , drop = FALSE]
historical <- dat[dat$trt == 0, , drop = FALSE]
model_columns <- c(
  "age", "male", "race_Black", "race_Other", "hispanic",
  "high_risk_cyto", "asct"
)
age_breaks <- c(-Inf, 50, 60, 70, Inf)
age_labels <- c("<50", "50-59", "60-69", "70+")
trial_age_group <- cut(
  trial$age, age_breaks, right = FALSE, labels = age_labels
)
historical_age_group <- cut(
  historical$age, age_breaks, right = FALSE, labels = age_labels
)
make_display_profile <- function(x, age_group) {
  data.frame(
    Age_group = as.character(age_group),
    Sex = ifelse(x$male == 1, "Male", "Female"),
    Race = ifelse(
      x$race_Black == 1, "Black",
      ifelse(x$race_Other == 1, "Other", "White")
    ),
    Hispanic_Latino = ifelse(x$hispanic == 1, "Yes", "No"),
    Cytogenetic_risk = ifelse(x$high_risk_cyto == 1, "High", "Standard"),
    ASCT = ifelse(x$asct == 1, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}
trial_display <- make_display_profile(trial, trial_age_group)
historical_display <- make_display_profile(historical, historical_age_group)
combination_key <- function(x) apply(x, 1, paste, collapse = "|")
trial_combination_key <- combination_key(trial_display)
historical_combination_key <- combination_key(historical_display)
historical_combination_n <- table(historical_combination_key)
ucmm_combination_n <- as.integer(
  historical_combination_n[trial_combination_key]
)
ucmm_combination_n[is.na(ucmm_combination_n)] <- 0L
distinct_trial_combination <- !duplicated(trial_combination_key)
ucmm_displayed_combination_n <- sum(
  ucmm_combination_n[distinct_trial_combination]
)
stopifnot(
  sum(historical_combination_n) == nrow(historical),
  length(ucmm_combination_n) == nrow(trial),
  all(ucmm_combination_n >= 0L)
)

stopifnot(
  nrow(trial) == 30L,
  sum(dat$trt == 0) == 200L,
  !anyNA(trial[, model_columns]),
  nrow(unique(trial[, model_columns])) == 30L,
  length(cal$V_profile_f) == nrow(trial),
  identical(cal$outcome, "PFS"),
  cal$H_f == 10L,
  cal$H_g == 5L
)

target_row <- match("100", as.character(cal$targets$target_name))
stopifnot(length(target_row) == 1L, !is.na(target_row))
s0_sq <- as.numeric(cal$targets$s0_sq[target_row])
nu0 <- as.numeric(cal$nu0)
n_scale <- 100000L
scale_seed <- 460116L
set.seed(scale_seed)
tau0_sq <- nu0 * s0_sq / stats::rchisq(n_scale, nu0)

ess <- vapply(
  cal$V_profile_f,
  function(V_f) mean(cal$sigma1_sq / (V_f + cal$H_g * tau0_sq)),
  numeric(1)
)

profile_ess <- data.frame(
  Profile = seq_len(nrow(trial)),
  Age_group = trial_display$Age_group,
  Age_years = trial$age,
  Sex = trial_display$Sex,
  Race = trial_display$Race,
  Hispanic_Latino = trial_display$Hispanic_Latino,
  Cytogenetic_risk = trial_display$Cytogenetic_risk,
  ASCT = trial_display$ASCT,
  UCMM_covariate_combination_n = ucmm_combination_n,
  V_profile_f = as.numeric(cal$V_profile_f),
  Prior_ESS = ess,
  stringsAsFactors = FALSE
)
stopifnot(nrow(profile_ess) == 30L, all(is.finite(ess)), all(ess > 0))

write_json(
  profile_ess,
  file.path(out, "generated/profile_ess.json"),
  dataframe = "rows", pretty = TRUE, digits = 12
)

profile_table <- profile_ess[, c(
  "Age_group", "Age_years", "Sex", "Race", "Hispanic_Latino",
  "Cytogenetic_risk", "ASCT", "UCMM_covariate_combination_n", "Prior_ESS"
)]
profile_table <- profile_table[
  order(match(profile_table$Age_group, age_labels), profile_table$Age_years),
  , drop = FALSE
]
rownames(profile_table) <- NULL
write_json(
  profile_table,
  file.path(out, "generated/profile_table.json"),
  dataframe = "rows", pretty = TRUE, digits = 12
)

write_json(
  list(
    data_file = normalizePath(data_file),
    data_md5 = unname(tools::md5sum(data_file)),
    calibration_file = normalizePath(calibration_file),
    calibration_md5 = unname(tools::md5sum(calibration_file)),
    outcome = "PFS",
    historical_cohort = "UCMM verified triplet regimens",
    historical_n = sum(dat$trt == 0),
    trial_profiles = nrow(trial),
    unique_trial_profiles = nrow(unique(trial[, model_columns])),
    H_f = cal$H_f,
    H_g = cal$H_g,
    nu0 = nu0,
    s0_sq = s0_sq,
    sigma1_sq = cal$sigma1_sq,
    scale_draws = n_scale,
    scale_seed = scale_seed,
    UCMM_count_matching_variables = names(trial_display),
    displayed_covariate_combinations = sum(distinct_trial_combination),
    UCMM_patients_in_displayed_combinations = ucmm_displayed_combination_n,
    UCMM_patients_outside_displayed_combinations = nrow(historical) - ucmm_displayed_combination_n,
    ESS_range = range(ess)
  ),
  file.path(out, "generated/profile_provenance.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 12
)

write_json(
  list(
    rows = nrow(profile_table),
    unique_profiles = nrow(unique(trial[, model_columns])),
    finite_ESS_rows = sum(is.finite(profile_table$Prior_ESS)),
    PFS_only = TRUE,
    exact_age_used = TRUE,
    age_display_decimals = 1,
    historical_n = sum(dat$trt == 0),
    EloKRd_age_group_n = as.list(table(factor(
      trial_age_group, levels = age_labels
    ))),
    UCMM_covariate_combination_n_min = min(ucmm_combination_n),
    UCMM_covariate_combination_n_max = max(ucmm_combination_n),
    UCMM_covariate_combination_zero_rows = sum(ucmm_combination_n == 0L),
    UCMM_covariate_combination_row_sum = sum(ucmm_combination_n),
    displayed_covariate_combinations = sum(distinct_trial_combination),
    UCMM_patients_in_displayed_combinations = ucmm_displayed_combination_n,
    UCMM_patients_outside_displayed_combinations = nrow(historical) - ucmm_displayed_combination_n,
    ESS_min = min(ess),
    ESS_max = max(ess),
    V_profile_f_min = min(cal$V_profile_f),
    V_profile_f_max = max(cal$V_profile_f)
  ),
  file.path(out, "generated/profile_table_validation.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 12
)

header <- paste(c(
  "Age group", "Age", "Sex", "Race",
  "\\shortstack{Hispanic/\\\\Latino}",
  "\\shortstack{Cytogenetic\\\\risk}", "ASCT",
  "\\shortstack{UCMM\\\\cohort $n$}", "Prior ESS"
), collapse = " & ")
note <- paste0(
  "Each row is one EloKRd patient covariate profile. Age groups are $<50$, 50--59, 60--69 and $\\geq70$ years. ",
  "Age is rounded to one decimal year for display; calculations use exact age. UCMM cohort $n$ is the number of patients among all 200 triplet-treated UCMM controls matching that row's age group, sex, race, Hispanic/Latino status, cytogenetic risk and ASCT. Counts repeat when EloKRd rows share the same displayed combination and should not be summed down the table. ",
  "Prior ESS is informed by the full 200-patient triplet UCMM cohort and is expressed in uncensored normal-reference units for the conditional mean log-PFS at the same profile. ",
  "The same residual variance and globally selected PFS scale are used for all profiles. Values do not measure outcome compatibility and do not sum or average to overall ESS. ASCT is postinduction."
)

t <- c(
  "\\begingroup",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\renewcommand{\\arraystretch}{1.0}",
  "\\fontsize{9.5}{11}\\selectfont",
  "\\begin{longtable}{lrlllllrr}",
  "\\caption{UCMM-informed prior ESS for the 30 EloKRd patient covariate profiles, with UCMM covariate-combination counts (PFS).}\\label{tab:table7}\\\\",
  "\\toprule",
  paste0(header, " \\\\ \\midrule"),
  "\\endfirsthead",
  "\\multicolumn{9}{l}{Table \\thetable\\ (continued)}\\\\",
  "\\toprule",
  paste0(header, " \\\\ \\midrule"),
  "\\endhead",
  "\\midrule\\multicolumn{9}{r}{Continued on next page}\\\\\\endfoot",
  "\\bottomrule\\endlastfoot"
)
for (i in seq_len(nrow(profile_table))) {
  new_group <- i == 1L || profile_table$Age_group[i] != profile_table$Age_group[i - 1L]
  if (i > 1L && new_group) t <- c(t, "\\addlinespace[3pt]")
  age_group_cell <- switch(
    profile_table$Age_group[i],
    "<50" = "$<50$", "50-59" = "50--59",
    "60-69" = "60--69", "70+" = "$\\geq70$"
  )
  row <- c(
    age_group_cell,
    sprintf("%.1f", profile_table$Age_years[i]),
    profile_table$Sex[i],
    profile_table$Race[i],
    profile_table$Hispanic_Latino[i],
    profile_table$Cytogenetic_risk[i],
    profile_table$ASCT[i],
    profile_table$UCMM_covariate_combination_n[i],
    sprintf("%.1f", profile_table$Prior_ESS[i])
  )
  t <- c(t, paste0(paste(row, collapse = " & "), " \\\\"))
}
t <- c(
  t,
  "\\end{longtable}",
  paste0(
    "\\noindent\\begin{minipage}{\\linewidth}\\footnotesize ",
    note, "\\end{minipage}"
  ),
  "\\endgroup"
)
writeLines(t, file.path(out, "generated/table7.tex"))

cat(
  "Table 7:", nrow(profile_table), "EloKRd profiles; PFS prior ESS range",
  sprintf("%.6f to %.6f.\n", min(ess), max(ess))
)
