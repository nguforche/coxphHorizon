############################################################
# Cox Export Toolkit 
# Exports ONLY:
#   - Cox coefficients (beta_hat)
#   - Baseline survival S0(t*) at selected horizons (days)
#   - Minimal metadata (variable names, encoding notes)
#
############################################################

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tibble)
  library(jsonlite) # for metadata JSON (optional)
})

# ---- Helpers ----

# Right-continuous step function from survfit baseline (S0)
make_S0_stepfun <- function(time_vec, surv_vec) {
  stats::stepfun(x = time_vec, y = c(1, surv_vec), right = TRUE)
}

# Months → days 
m2d <- function(m) round(m * 30.4375)

# ---- Core Export Function ----

#' Export Cox parameters 
#' @param fit coxph object (trained in a secure environment) 
#' @param horizons_days numeric vector of horizons (in days) to export S0(t*)
#' @param var_order optional explicit covariate order; defaults to names(coef(fit))
#' @param out_dir directory to write export artifacts (CSV/RDS/JSON)
#' @param prefix filename prefix, e.g. "dcrc_cox_v1"
#' @return list(beta_hat, S0_table, meta); also writes files to out_dir
cox_export <- function(fit,
                               horizons_days,
                               var_order = NULL,
                               out_dir = ".",
                               prefix = "cox_export") {
  stopifnot(inherits(fit, "coxph"))
  beta_hat <- stats::coef(fit)
  if (is.null(beta_hat)) stop("Cox model has no coefficients; check convergence.")
  if (is.null(var_order)) var_order <- names(beta_hat)

  # Baseline survival from survfit (lp=0)
  sf0 <- survival::survfit(fit)
  S0_fun <- make_S0_stepfun(sf0$time, sf0$surv)
  S0_vals <- S0_fun(horizons_days)

  S0_table <- tibble(time = as.numeric(horizons_days),
                     S0   = as.numeric(S0_vals)) %>%
    distinct(time, .keep_all = TRUE) %>%
    arrange(time)

  meta <- list(
    n_train          = fit$n,
    events_train     = fit$nevent,
    ties_method      = fit$method,
    variables        = unname(var_order),
    encoding_notes   = "Binary covariates must be encoded 0/1; continuous as-is. If any scaling/one-hot used, export means/SDs/levels here.",
    option           = "B (coefficients + S0 at horizons)",
    created_utc      = format(Sys.time(), tz = "UTC")
  )

  # ---- Write artifacts (no PHI) ----
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  coef_df <- tibble(variable = names(beta_hat), coef = as.numeric(beta_hat))
  readr::write_csv(coef_df, file.path(out_dir, sprintf("%s_coefficients.csv", prefix)))
  readr::write_csv(S0_table, file.path(out_dir, sprintf("%s_S0_table.csv", prefix)))
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
             file.path(out_dir, sprintf("%s_meta.json", prefix)))

  # Also provide a single RDS bundle for convenience (optional)
  saveRDS(list(beta_hat = beta_hat, S0_table = S0_table, meta = meta),
          file = file.path(out_dir, sprintf("%s_bundle.rds", prefix)))

  invisible(list(beta_hat = beta_hat, S0_table = S0_table, meta = meta))
}

# ----------compute optimal thresholds on TRAINING set  ----------

#' Compute optimal threshold per horizon on the training set 
#' @param fit   coxph trained model
#' @param train_df data frame with covariates + time/status 
#' @param covariates character vector (names in train_df) used in the Cox model
#' @param S0_table tibble(time,S0) from cox_export_optionB()
#' @param horizons_days numeric vector of horizons (days)
#' @param method "youden" or "closest.topleft"
#' @return tibble(time, threshold, method, n_known)
compute_thresholds <- function(fit, train_df, covariates, S0_table, horizons_days,
                                       method = c("youden","closest.topleft")) {
  method <- match.arg(method)
  # 1) compute lp and risk at each horizon
  beta_hat <- stats::coef(fit)
  stopifnot(all(covariates %in% names(train_df)))
  X <- as.matrix(train_df[, covariates, drop = FALSE])
  lp <- as.numeric(X %*% beta_hat)

  # helper to get survival => risk at t*
  predict_surv <- function(lp, t_star, S0_table) {
    S0_table <- S0_table[order(S0_table$time), ]
    idx <- max(which(S0_table$time <= t_star))
    if (is.infinite(idx)) idx <- 1L
    S0_t <- S0_table$S0[idx]
    S0_t ^ exp(lp)
  }

  rows <- lapply(horizons_days, function(t_star) {
    risk <- 1 - predict_surv(lp, t_star, S0_table)

    # build KNOWN set (label at horizon) from training data if time/status exist
    known <- rep(TRUE, length(risk))
    y_true <- NULL
    if (all(c("time","status") %in% names(train_df))) {
      time   <- as.numeric(train_df$time)
      status <- as.integer(train_df$status)
      known  <- (status == 1 & time <= t_star) | (time > t_star)
      y_true <- as.integer(status == 1 & time <= t_star)
    } else {
      # If time/status not present in the exported training frame,
      # we still compute thresholds via risk quantiles (fallback).
      # But best practice: pass time/status so thresholds align with your original training procedure.
      y_true <- rep(NA_integer_, length(risk))
    }

    rk <- risk[known]; yk <- y_true[known]
    # guard: need variation
    if (length(unique(na.omit(rk))) <= 1L || sum(known) == 0L ||
        (!all(is.na(yk)) && (sum(yk==1, na.rm=TRUE) == 0L || sum(yk==0, na.rm=TRUE) == 0L))) {
      return(tibble(time = t_star, threshold = NA_real_, method = method, n_known = sum(known)))
    }

    # threshold grid
    thr_unique <- sort(unique(rk))
    mids <- if (length(thr_unique) > 1L) (head(thr_unique,-1) + tail(thr_unique,-1))/2 else thr_unique
    thrs <- unique(c(min(thr_unique)-1e-12, mids, max(thr_unique)+1e-12))

    # If yk is NA (no time/status), default to the 80th percentile 
    if (all(is.na(yk))) {
      thr_star <- as.numeric(stats::quantile(rk, probs = 0.80, na.rm = TRUE))
      return(tibble(time = t_star, threshold = thr_star, method = paste0(method,"_no_labels"), n_known = sum(known)))
    }

    ss <- vapply(thrs, function(t0) {
      yp <- as.integer(rk >= t0)
      sens <- ifelse(sum(yk==1) > 0, sum(yp==1 & yk==1) / sum(yk==1), NA_real_)
      spec <- ifelse(sum(yk==0) > 0, sum(yp==0 & yk==0) / sum(yk==0), NA_real_)
      c(sens, spec)
    }, numeric(2))
    SENS <- ss[1,]; SPEC <- ss[2,]
    j <- if (method == "youden") which.max(SENS + SPEC - 1) else which.min(((1 - SPEC) - 0)^2 + (1 - SENS)^2)
    if (!is.finite(j) || length(j) < 1L) {
      thr_star <- NA_real_
    } else {
      thr_star <- as.numeric(thrs[j])
    }
    tibble(time = t_star, threshold = thr_star, method = method, n_known = sum(known))
  })

  dplyr::bind_rows(rows) |>
    dplyr::arrange(time)
}

# ---------- Export thresholds CSV alongside coeffs + S0_table ----------

#' Export Option B artifacts + thresholds
#' @param fit coxph model
#' @param train_df training data frame including time/status (preferred)
#' @param horizons_days numeric horizons
#' @param out_dir output dir
#' @param prefix file prefix
#' @param method threshold selection method
#' @return invisibly list(beta_hat, S0_table, thresholds)
cox_export_with_thresholds <- function(fit, train_df, horizons_days,
                                               out_dir = ".", prefix = "cox_export",
                                               method = c("youden","closest.topleft")) {
  method <- match.arg(method)

  # 1) export coeffs + S0_table
  exported <- cox_export(fit, horizons_days = horizons_days, out_dir = out_dir, prefix = prefix)
  beta_hat  <- exported$beta_hat
  S0_table  <- exported$S0_table
  covariates <- names(beta_hat)

  # 2) compute thresholds on training
  thresholds <- compute_thresholds(fit, train_df, covariates, S0_table, horizons_days, method = method)

  # 3) write thresholds
  readr::write_csv(thresholds, file.path(out_dir, sprintf("%s_thresholds.csv", prefix)))

  invisible(list(beta_hat = beta_hat, S0_table = S0_table, thresholds = thresholds, meta = exported$meta))
}




