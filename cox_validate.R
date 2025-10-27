############################################################
# Validation Toolkit —
# Uses ONLY:
#   - Exported coefficients (beta_hat)
#   - Exported baseline survival S0(t*) at selected horizons
# Computes:
#   - Predicted survival & risk at t*
#   - Time-dependent AUC (IPCW), Sens/Spec (Youden/closest), PPV, NPV, F1
#
# Dependencies: survival, timeROC, dplyr, tibble, purrr, readr, jsonlite
############################################################

suppressPackageStartupMessages({
  library(survival)
  library(timeROC)
  library(survAUC)  
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(jsonlite)
  library(tidyr)
  library(scales)
})

# ---- Utilities ----

m2d <- function(m) round(m * 30.4375)

# Force a single numeric scalar (NA if empty)
scalar_num <- function(x) {
  if (length(x) < 1L || all(is.na(x))) return(NA_real_)
  as.numeric(x[[1]])
}

# ---------- metrics at horizon using a PROVIDED threshold ----------
metrics_at_horizon_with_threshold <- function(time, status, risk, t_star, threshold, compute_auc = TRUE,
                                              engine = c("uno","timeROC"), iid = FALSE) {
  engine <- match.arg(engine)
  time   <- as.numeric(time)
  status <- as.integer(status)
  risk   <- as.numeric(risk)
  status[!(status %in% c(0L,1L))] <- 0L
  ok <- is.finite(time) & is.finite(status) & is.finite(risk)
  time <- time[ok]; status <- status[ok]; risk <- risk[ok]

  # AUC 
  AUC <- NA_real_
  if (isTRUE(compute_auc) && length(unique(na.omit(risk))) > 1L &&
      t_star > min(time) && t_star < max(time) &&
      sum(status == 1 & time <= t_star) > 0L && sum(time > t_star) > 0L) {

    if (engine == "uno") {
      auc_obj <- tryCatch(
        survAUC::AUC.uno(Surv(time, status), Surv(time, status), lpnew = risk, times = t_star),
        error = function(e) NULL
      )
      AUC <- if (is.null(auc_obj)) NA_real_ else as.numeric(auc_obj$iauc)
    } else {
      tr <- tryCatch(
        timeROC::timeROC(T = time, delta = status, marker = risk,
                         cause = 1, times = t_star, iid = iid, weighting = "marginal"),
        error = function(e) NULL
      )
      if (!is.null(tr)) {
        idx <- if (!is.null(tr$times)) which.min(abs(tr$times - t_star)) else NCOL(tr$FP)
        AUC <- suppressWarnings(as.numeric(tr$AUC[idx]))
      }
    }
  }

  # If no threshold provided (NA), return AUC-only row with NA for others
  if (!is.finite(threshold)) {
    return(list(t_days = t_star, AUC = AUC, Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }

  # Apply provided threshold
  y_pred <- as.integer(risk >= as.numeric(threshold))
  y_true <- as.integer(status == 1 & time <= t_star)
  known  <- (status == 1 & time <= t_star) | (time > t_star)

  if (!any(known)) {
    return(list(t_days = t_star, AUC = AUC, Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = as.numeric(threshold)))
  }

  TPn <- sum(y_pred[known] == 1 & y_true[known] == 1)
  FPn <- sum(y_pred[known] == 1 & y_true[known] == 0)
  TNn <- sum(y_pred[known] == 0 & y_true[known] == 0)
  FNn <- sum(y_pred[known] == 0 & y_true[known] == 1)

  sens <- ifelse((TPn + FNn) > 0, TPn / (TPn + FNn), NA_real_)
  spec <- ifelse((TNn + FPn) > 0, TNn / (TNn + FPn), NA_real_)
  ppv  <- ifelse((TPn + FPn) > 0, TPn / (TPn + FPn), NA_real_)
  npv  <- ifelse((TNn + FNn) > 0, TNn / (TNn + FNn), NA_real_)
  f1   <- ifelse((2*TPn + FPn + FNn) > 0, 2*TPn / (2*TPn + FPn + FNn), NA_real_)

  list(t_days = t_star, AUC = AUC, Sensitivity = sens, Specificity = spec,
       PPV = ppv, NPV = npv, F1 = f1, Threshold = as.numeric(threshold))
}

# ---------- join thresholds to horizons and compute metrics ----------
# thresholds_tbl must have columns: time, threshold
metrics_with_thresholds_over_horizons <- function(time, status, risk_by_t, thresholds_tbl,
                                                  engine = c("uno","timeROC"), iid = FALSE,
                                                  compute_auc = TRUE) {
  engine <- match.arg(engine)
  horizons_days <- as.numeric(names(risk_by_t))

  # Map threshold per horizon (exact match; if needed, use nearest)
  thr_map <- setNames(rep(NA_real_, length(horizons_days)), as.character(horizons_days))
  if (!is.null(thresholds_tbl) && all(c("time","threshold") %in% names(thresholds_tbl))) {
    # exact match first
    idx_exact <- match(horizons_days, thresholds_tbl$time)
    has_exact <- !is.na(idx_exact)
    thr_map[as.character(horizons_days[has_exact])] <- thresholds_tbl$threshold[idx_exact[has_exact]]
    # optional: nearest match for those missing
    if (any(!has_exact)) {
      for (t_miss in horizons_days[!has_exact]) {
        j <- which.min(abs(thresholds_tbl$time - t_miss))
        if (length(j) == 1L && is.finite(j)) {
          thr_map[as.character(t_miss)] <- thresholds_tbl$threshold[j]
        }
      }
    }
  }

  rows <- lapply(seq_along(horizons_days), function(i) {
    t_star <- horizons_days[i]
    risk   <- risk_by_t[[i]]
    thr_i  <- thr_map[[as.character(t_star)]]
    z <- metrics_at_horizon_with_threshold(time, status, risk, t_star,
                                           threshold = thr_i,
                                           compute_auc = compute_auc,
                                           engine = engine, iid = iid)
    tibble::tibble(
      t_days     = as.numeric(z$t_days),
      AUC        = as.numeric(z$AUC),
      Sensitivity= as.numeric(z$Sensitivity),
      Specificity= as.numeric(z$Specificity),
      PPV        = as.numeric(z$PPV),
      NPV        = as.numeric(z$NPV),
      F1         = as.numeric(z$F1),
      Threshold  = as.numeric(z$Threshold)
    )
  })

  dplyr::bind_rows(rows) |> dplyr::arrange(t_days)
}



#' Predict survival at t* (days) via Option B
#' S_i(t*) = S0(t*) ^ exp(lp_i)
#' @param lp numeric vector of linear predictors (X %*% beta_hat)
#' @param t_star single horizon (days)
#' @param S0_table tibble(time,S0) exported from Optum
#' @param interpolate if TRUE, interpolate log(S0) for non-listed t*
#' @return numeric vector of S_i(t*)
predict_surv <- function(lp, t_star, S0_table, interpolate = FALSE) {
  stopifnot(all(c("time","S0") %in% names(S0_table)))
  S0_table <- S0_table %>% arrange(time)
  if (isTRUE(interpolate)) {
    f <- stats::approxfun(x = S0_table$time,
                          y = log(pmax(S0_table$S0, .Machine$double.eps)),
                          rule = 2)
    S0_t <- exp(f(t_star))
  } else {
    idx <- max(which(S0_table$time <= t_star))
    if (is.infinite(idx)) idx <- 1L
    S0_t <- S0_table$S0[idx]
  }
  S0_t ^ exp(lp)
}

#' Time-dependent metrics at t*

# Compute metrics at a single horizon using Uno's AUC (fast)
metrics_at_horizon_uno <- function(time, status, risk, t_star,
                                   cutoff_strategy = c("youden","closest.topleft")) {
  cutoff_strategy <- match.arg(cutoff_strategy)

  # Degenerate marker guard
  if (length(unique(na.omit(risk))) <= 1L) {
    return(list(t_days = t_star, AUC = NA_real_,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }

  # AUC (Uno) — fast
  auc_obj <- tryCatch(
    survAUC::AUC.uno(Surv.rsp = Surv(time, status),
                     Surv.rsp.new = Surv(time, status),
                     lpnew = risk, times = t_star),
    error = function(e) NULL
  )
  AUC <- if (is.null(auc_obj)) NA_real_ else scalar_num(auc_obj$iauc)

  # Known label set at t*: event by t* vs observed event-free past t*
  known <- (status == 1 & time <= t_star) | (time > t_star)
  if (!any(known, na.rm = TRUE)) {
    return(list(t_days = t_star, AUC = AUC,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }
  y_true <- as.integer(status == 1 & time <= t_star)
  rk <- risk[known]; yk <- y_true[known]

  # If all risks identical in known set → no ROC
  if (length(unique(na.omit(rk))) <= 1L) {
    return(list(t_days = t_star, AUC = AUC,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }

  # Threshold grid: unique risks + small padding
  thr <- sort(unique(rk))
  # midpoints between unique thresholds to avoid ties at boundaries
  mid <- if (length(thr) > 1L) (head(thr, -1) + tail(thr, -1))/2 else thr
  thrs <- unique(c(min(thr) - 1e-12, mid, max(thr) + 1e-12))

  # Sens/Spec grid
  sens_spec <- vapply(thrs, function(t0) {
    yp <- as.integer(rk >= t0)
    sens <- ifelse(sum(yk == 1) > 0, sum(yp == 1 & yk == 1) / sum(yk == 1), NA_real_)
    spec <- ifelse(sum(yk == 0) > 0, sum(yp == 0 & yk == 0) / sum(yk == 0), NA_real_)
    c(sens, spec)
  }, numeric(2))
  SENS <- sens_spec[1, ]; SPEC <- sens_spec[2, ]
    
  j <- if (cutoff_strategy == "youden") {
    # maximize (sens + spec - 1)
    which.max(SENS + SPEC - 1)
  } else {
    # closest to (0,1)
    which.min(((1 - SPEC) - 0)^2 + (1 - SENS)^2)
  }
  if (!is.finite(j) || length(j) < 1L) {
    return(list(t_days = t_star, AUC = AUC,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }

  thr_star <- scalar_num(thrs[j])
  y_pred   <- as.integer(risk >= thr_star)

  # Confusion on KNOWN set
  TPn <- sum(y_pred[known] == 1 & y_true[known] == 1)
  FPn <- sum(y_pred[known] == 1 & y_true[known] == 0)
  TNn <- sum(y_pred[known] == 0 & y_true[known] == 0)
  FNn <- sum(y_pred[known] == 0 & y_true[known] == 1)

  sens <- ifelse((TPn + FNn) > 0, TPn / (TPn + FNn), NA_real_)
  spec <- ifelse((TNn + FPn) > 0, TNn / (TNn + FPn), NA_real_)
  ppv  <- ifelse((TPn + FPn) > 0, TPn / (TPn + FPn), NA_real_)
  npv  <- ifelse((TNn + FNn) > 0, TNn / (TNn + FNn), NA_real_)
  f1   <- ifelse((2*TPn + FPn + FNn) > 0, 2*TPn / (2*TPn + FPn + FNn), NA_real_)

  list(t_days = t_star, AUC = AUC,
       Sensitivity = sens, Specificity = spec,
       PPV = ppv, NPV = npv, F1 = f1, Threshold = thr_star)
}

#' Uses IPCW AUC (timeROC) + threshold selection (Youden or closest-to-(0,1))
metrics_at_horizon_timeROC <- function(time, status, risk, t_star,
                                       cutoff_strategy = c("youden","closest.topleft"),
                                       iid = FALSE,
                                       thresholds) {
  cutoff_strategy <- match.arg(cutoff_strategy)


  if (length(unique(na.omit(risk))) <= 1L) {
    return(list(t_days = t_star, AUC = NA_real_,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }

  tr <- tryCatch(
    timeROC::timeROC(T = time, delta = status, marker = risk,
                     cause = 1, times = t_star, iid = iid, weighting = "marginal"),
    error = function(e) NULL
  )
  
  if (is.null(tr)) {

    return(list(t_days = t_star, AUC = NA_real_,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }


  AUC <- scalar_num(tr$AUC[2])
  FP  <- tryCatch(tr$FP[,2], error=function(e) numeric(0))
  TP  <- tryCatch(tr$TP[,2], error=function(e) numeric(0))
   
  if (length(FP) == 0L || length(TP) == 0L) {
    return(list(t_days = t_star, AUC = AUC,
                Sensitivity = NA_real_, Specificity = NA_real_,
                PPV = NA_real_, NPV = NA_real_, F1 = NA_real_, Threshold = NA_real_))
  }


  j <- if (cutoff_strategy == "youden") which.max(TP - FP)
       else which.min((FP - 0)^2 + (1 - TP)^2)
  sens <- scalar_num(TP[j])
  spec <- scalar_num(1 - FP[j])
  thr_star <- scalar_num(thr[j])

print(sens)

  known <- (status == 1 & time <= t_star) | (time > t_star)
  y_true <- as.integer(status == 1 & time <= t_star)
  y_pred <- as.integer(risk >= thr_star)

  y_true_k <- y_true[known]; y_pred_k <- y_pred[known]
  TPn <- sum(y_pred_k == 1 & y_true_k == 1)
  FPn <- sum(y_pred_k == 1 & y_true_k == 0)
  TNn <- sum(y_pred_k == 0 & y_true_k == 0)
  FNn <- sum(y_pred_k == 0 & y_true_k == 1)

  ppv <- ifelse((TPn + FPn) > 0, TPn / (TPn + FPn), NA_real_)
  npv <- ifelse((TNn + FNn) > 0, TNn / (TNn + FNn), NA_real_)
  f1  <- ifelse((2*TPn + FPn + FNn) > 0, 2*TPn / (2*TPn + FPn + FNn), NA_real_)

  list(t_days = t_star, AUC = AUC,
       Sensitivity = sens, Specificity = spec, PPV = ppv, NPV = npv,
       F1 = f1, Threshold = thr_star)
}

# Vectorized over horizons; builds rows and binds (no mutate/rowwise)
metrics_over_horizons <- function(time, status, risk_by_t,
                                  cutoff_strategy = "youden",
                                  iid = FALSE, 
                                  engine = c("uno","timeROC")
                                 ) {
  engine <- match.arg(engine)

  horizons_days <- as.numeric(names(risk_by_t))
  res_list <- lapply(seq_along(horizons_days), function(i) {
    t_star <- horizons_days[i]
    risk   <- risk_by_t[[i]]
    
    #metrics_at_horizon_uno(time, status, risk, t_star, cutoff_strategy)
      
   if (engine == "uno") {
      metrics_at_horizon_uno(time, status, risk, t_star, cutoff_strategy)
    } else {
      metrics_at_horizon_timeROC(time, status, risk, t_star, cutoff_strategy, iid = iid)
    }

      
  })
    
  # Ensure every field is scalar
  out <- tibble::tibble(
    t_days = vapply(res_list, function(z) scalar_num(z$t_days), numeric(1)),
    AUC    = vapply(res_list, function(z) scalar_num(z$AUC), numeric(1)),
    Sensitivity = vapply(res_list, function(z) scalar_num(z$Sensitivity), numeric(1)),
    Specificity = vapply(res_list, function(z) scalar_num(z$Specificity), numeric(1)),
    PPV    = vapply(res_list, function(z) scalar_num(z$PPV), numeric(1)),
    NPV    = vapply(res_list, function(z) scalar_num(z$NPV), numeric(1)),
    F1     = vapply(res_list, function(z) scalar_num(z$F1), numeric(1)),
    Threshold = vapply(res_list, function(z) scalar_num(z$Threshold), numeric(1))
  )
  out[order(out$t_days), , drop = FALSE]
}


metrics_over_horizons1 <- function(time, status, risk_by_t,
                                  cutoff_strategy = "youden",
                                  engine = c("uno"),  # keep "uno" default
                                  iid = FALSE) {      # kept for API symmetry
  engine <- match.arg(engine)
  horizons_days <- as.numeric(names(risk_by_t))

  rows <- lapply(seq_along(horizons_days), function(i) {
    t_star <- horizons_days[i]
    risk   <- risk_by_t[[i]]

    if (engine == "uno") {
      z <- metrics_at_horizon_uno(time, status, risk, t_star, cutoff_strategy)
    } else {
      stop("Only engine='uno' is enabled in this build.")
    }

    # Build a one-row tibble with scalar fields
    tibble(
      t_days     = scalar_num(z$t_days),
      AUC        = scalar_num(z$AUC),
      Sensitivity= scalar_num(z$Sensitivity),
      Specificity= scalar_num(z$Specificity),
      PPV        = scalar_num(z$PPV),
      NPV        = scalar_num(z$NPV),
      F1         = scalar_num(z$F1),
      Threshold  = scalar_num(z$Threshold)
    )
  })

  bind_rows(rows) %>% arrange(t_days)
}

                  
# ---- Main validation function ----

#' Validate on test/external data
#' @param df data.frame with time (days), status {0,1}, and covariates
#' @param time_col, status_col column names for time/status
#' @param coef_csv path to exported coefficients CSV (variable, coef)
#' @param s0_csv path to exported S0 table CSV (time, S0)
#' @param horizons_days numeric vector of horizons (days) to evaluate
#' @param cutoff_strategy "youden" or "closest.topleft"
#' @param engine "uno" or "timeRoc"                  
#' @return list(metrics_table, risks_by_horizon, lp)
validate <- function(df,
                             time_col = "time",
                             status_col = "status",
                             coef_csv,
                             s0_csv,
                             horizons_days,
                             cutoff_strategy = "youden",   # if threshold is not available
                             engine = c("uno","timeROC"),  # only for AUC 
                             iid = FALSE,
                             thresholds_csv = NULL,        
                             compute_auc = TRUE) {

  engine <- match.arg(engine)
  coef_df  <- readr::read_csv(coef_csv, show_col_types = FALSE)
  S0_table <- readr::read_csv(s0_csv,   show_col_types = FALSE)
  thresh_df <- NULL
  if (!is.null(thresholds_csv) && file.exists(thresholds_csv)) {
    thresh_df <- readr::read_csv(thresholds_csv, show_col_types = FALSE)
  }

  stopifnot(all(coef_df$variable %in% names(df)))
  X <- as.matrix(df[, coef_df$variable, drop = FALSE])
  lp <- as.numeric(X %*% coef_df$coef)

  predict_surv <- function(lp, t_star, S0_table) {
    S0_table <- S0_table[order(S0_table$time), ]
    idx <- max(which(S0_table$time <= t_star))
    if (is.infinite(idx)) idx <- 1L
    S0_t <- S0_table$S0[idx]
    S0_t ^ exp(lp)
  }

  risk_by_t <- purrr::map(horizons_days, function(tt) {
    1 - predict_surv(lp, t_star = tt, S0_table = S0_table)
  })
  names(risk_by_t) <- as.character(horizons_days)

  if (!is.null(thresh_df)) {
    mt <- metrics_with_thresholds_over_horizons(
      time = df[[time_col]],
      status = df[[status_col]],
      risk_by_t = risk_by_t,
      thresholds_tbl = thresh_df,
      engine = engine,
      iid = iid,
      compute_auc = compute_auc
    )
  } else {
    # fallback to auto-thresholding (existing path)
    mt <- metrics_over_horizons(
      time = df[[time_col]],
      status = df[[status_col]],
      risk_by_t = risk_by_t,
      cutoff_strategy = cutoff_strategy,
      engine = engine,
      iid = iid
    )
  }

  list(metrics_table = mt, risks_by_horizon = risk_by_t, lp = lp)
}



