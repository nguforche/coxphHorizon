############################################################
# Simulation Demo —
# - Simulates a training cohort, fits Cox, exports artifacts
# - Simulates Mayo-like external cohort
# - Reports metrics across t = {1,3,6,9,12} months and their averages
#
############################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
})

source("cox_export.R")
source("cox_validate.R")

m2d <- function(m) round(m * 30.4375)

set.seed(2025)

# ---- Simulate training (Optum-like) ----
n_tr <- 6000
train <- tibble(
  age    = rnorm(n_tr, 60, 12),
  female = rbinom(n_tr, 1, 0.5),
  htn    = rbinom(n_tr, 1, 0.45),
  a1c    = rnorm(n_tr, 7.8, 1.4)
)
beta_true <- c(age = 0.02, female = -0.12, htn = 0.35, a1c = 0.22)
lp_tr     <- as.numeric(as.matrix(train[, names(beta_true)]) %*% beta_true)
lambda0   <- 6e-4
t_event   <- rexp(n_tr, rate = lambda0 * exp(lp_tr))
t_cens    <- runif(n_tr, min = m2d(6), max = m2d(30))
time_tr   <- pmin(t_event, t_cens)
status_tr <- as.integer(t_event <= t_cens)
train$time <- time_tr; train$status <- status_tr

# Fit Cox
fit <- coxph(Surv(time, status) ~ age + female + htn + a1c, data = train,
             ties = "breslow", x = TRUE, y = TRUE)

# Export Option B artifacts
horizons_month <- c(1,3,6,9,12)
horizons_days <- m2d(horizons_month)
export_dir <- "optum_exports_sim"



# This will add dcrc_cox_v1_thresholds.csv with columns: time, threshold, method, n_known
bundle <- cox_export_with_thresholds(fit, train_df = train,
                                    horizons_days = horizons_days,
                                    out_dir = export_dir, 
                                    prefix = "dcrc_cox_v1",
                                    method = "youden")

# ---- Simulate external validation data ----
n_ex <- 3500
ext <- tibble(
  age    = rnorm(n_ex, 62, 13),
  female = rbinom(n_ex, 1, 0.45),
  htn    = rbinom(n_ex, 1, 0.5),
  a1c    = rnorm(n_ex, 7.6, 1.5)
)
beta_true_ext <- c(age = 0.018, female = -0.10, htn = 0.32, a1c = 0.24)
lp_ex   <- as.numeric(as.matrix(ext[, names(beta_true_ext)]) %*% beta_true_ext)
lambda0_ext <- 6.5e-4
t_event_ex  <- rexp(n_ex, rate = lambda0_ext * exp(lp_ex))
t_cens_ex   <- runif(n_ex, min = m2d(6), max = m2d(36))
ext$time    <- pmin(t_event_ex, t_cens_ex)
ext$status  <- as.integer(t_event_ex <= t_cens_ex)



res <- validate(
  mayo_df = ext,
  coef_csv = file.path(export_dir, "dcrc_cox_v1_coefficients.csv"), 
  s0_csv   = file.path(export_dir, "dcrc_cox_v1_S0_table.csv"), 
  horizons_days = horizons_days,
  engine = c("uno", "timeROC")[2],     # "timeROC" is more robust (applies IPCW AUC)
  thresholds_csv = file.path(export_dir, "dcrc_cox_v1_thresholds.csv"),   # <-- use training thresholds
  compute_auc = TRUE
)


cat("\nMetrics by horizon (days):\n")
print(res$metrics_table)

cat("\nAverage metrics across t = {", paste(horizons_m, collapse = ","), "} months:\n", sep = "")
avg <- res$metrics_table %>%
  summarise(across(c(AUC,Sensitivity,Specificity,PPV,NPV,F1), \(x) mean(x, na.rm = TRUE)))
print(avg)






