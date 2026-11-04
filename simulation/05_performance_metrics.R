# 05_performance_metrics.R ----------------------------------------------------
# Protocol sections 18-23: predictions in the fixed test population and
# prediction / clinical-utility metrics.

# --- Net benefit (protocol section 23) ----------------------------------------
net_benefit <- function(y, p, threshold) {
  n <- length(y)
  positive <- p >= threshold
  tp <- sum(positive & y == 1)
  fp <- sum(positive & y == 0)
  tp / n - fp / n * threshold / (1 - threshold)
}

# --- Calibration intercept and slope (protocol section 22) --------------------
# Mathematically equivalent, fast implementations of the two calibration fits.
#
# Calibration intercept: MLE of glm(Y ~ 1, offset = lp) is the value a such
# that sum(plogis(a + lp)) = sum(Y). Instead of glm's IRLS (slow on 1e5 rows)
# we solve this scalar equation with uniroot to 1e-10 precision; the solution
# is identical to the glm estimate.
calibration_intercept <- function(y, lp, tol = 1e-10) {
  target <- sum(y)
  f <- function(a) sum(plogis(a + lp)) - target
  uniroot(f, interval = c(-20, 20), tol = tol)$root
}

# Calibration slope: univariate logistic regression of Y on lp. Implemented
# with IRLS directly (weights w = p(1-p), working response z = eta + (y-p)/w)
# instead of glm(); converges to the identical MLE.
calibration_slope <- function(y, lp, tol = 1e-10, maxit = 100L) {
  b0 <- 0.0
  b1 <- 0.0
  for (i in seq_len(maxit)) {
    eta <- b0 + b1 * lp
    p <- plogis(eta)
    w <- p * (1 - p)
    w <- pmax(w, 1e-12)                    # guard against exact 0/1 probs
    z <- eta + (y - p) / w
    # Normal equations for Y ~ lp with weights w
    s_w <- sum(w); s_wx <- sum(w * lp); s_wxx <- sum(w * lp^2)
    s_wz <- sum(w * z); s_wxz <- sum(w * lp * z)
    det <- s_w * s_wxx - s_wx^2
    b0_new <- (s_wxx * s_wz - s_wx * s_wxz) / det
    b1_new <- (s_w * s_wxz - s_wx * s_wz) / det
    if (max(abs(b0_new - b0), abs(b1_new - b1)) < tol) {
      b0 <- b0_new; b1 <- b1_new
      break
    }
    b0 <- b0_new; b1 <- b1_new
  }
  c(b0 = b0, b1 = b1)
}

# --- All test-population metrics for one repetition ---------------------------
# p_ratio / p_component are 100k-length response predictions; only summary
# metrics are stored (protocol section 18).
compute_performance_metrics <- function(test_pop, p_ratio, p_component,
                                        quick = TRUE) {
  y <- test_pop$Y
  n_test <- length(y)

  # AUC (protocol section 19)
  auc_ratio <- as.numeric(pROC::auc(response = y, predictor = p_ratio, quiet = TRUE))
  auc_component <- as.numeric(pROC::auc(response = y, predictor = p_component, quiet = TRUE))
  delta_auc <- auc_component - auc_ratio

  # Brier score (protocol section 20)
  brier_ratio <- mean((y - p_ratio)^2)
  brier_component <- mean((y - p_component)^2)
  delta_brier <- brier_ratio - brier_component     # positive favors component

  # Log loss with clipping (protocol section 21)
  eps <- 1e-15
  p_ratio_clip <- pmin(pmax(p_ratio, eps), 1 - eps)
  p_component_clip <- pmin(pmax(p_component, eps), 1 - eps)
  logloss_ratio <- -mean(y * log(p_ratio_clip) + (1 - y) * log(1 - p_ratio_clip))
  logloss_component <- -mean(y * log(p_component_clip) + (1 - y) * log(1 - p_component_clip))
  delta_logloss <- logloss_ratio - logloss_component

  # Calibration intercept and slope (protocol section 22)
  lp_ratio <- qlogis(p_ratio_clip)
  lp_component <- qlogis(p_component_clip)
  if (quick) {
    cal_int_ratio <- calibration_intercept(y, lp_ratio)
    cal_int_component <- calibration_intercept(y, lp_component)
    s_ratio <- calibration_slope(y, lp_ratio)
    s_component <- calibration_slope(y, lp_component)
    cal_slope_ratio <- s_ratio[["b1"]]
    cal_slope_component <- s_component[["b1"]]
  } else {
    cal_int_ratio <- as.numeric(coef(glm(y ~ 1, offset = lp_ratio,
                                         family = binomial()))[1])
    cal_int_component <- as.numeric(coef(glm(y ~ 1, offset = lp_component,
                                             family = binomial()))[1])
    cal_slope_ratio <- as.numeric(coef(glm(y ~ lp_ratio, family = binomial()))["lp_ratio"])
    cal_slope_component <- as.numeric(coef(glm(y ~ lp_component, family = binomial()))["lp_component"])
  }

  # Net benefit at prespecified thresholds (protocol section 23)
  pi_test <- mean(y)
  out <- list()
  for (pt in c(0.10, 0.20, 0.30)) {
    tag <- sprintf("%03d", round(pt * 100))
    NB_ratio <- net_benefit(y, p_ratio, pt)
    NB_component <- net_benefit(y, p_component, pt)
    NB_all <- pi_test - (1 - pi_test) * pt / (1 - pt)
    NB_none <- 0
    out[[paste0("NB_ratio_", tag)]] <- NB_ratio
    out[[paste0("NB_component_", tag)]] <- NB_component
    out[[paste0("NB_all_", tag)]] <- NB_all
    out[[paste0("delta_NB_", tag)]] <- NB_component - NB_ratio
    out[[paste0("component_clinically_useful_", tag)]] <-
      as.integer(NB_component > max(NB_all, NB_none))
    out[[paste0("ratio_clinically_useful_", tag)]] <-
      as.integer(NB_ratio > max(NB_all, NB_none))
  }

  c(list(
    auc_ratio = auc_ratio, auc_component = auc_component, delta_auc = delta_auc,
    brier_ratio = brier_ratio, brier_component = brier_component, delta_brier = delta_brier,
    logloss_ratio = logloss_ratio, logloss_component = logloss_component,
    delta_logloss = delta_logloss,
    cal_intercept_ratio = cal_int_ratio, cal_intercept_component = cal_int_component,
    cal_slope_ratio = cal_slope_ratio, cal_slope_component = cal_slope_component
  ), out)
}

# --- Predictions on the fixed test population ---------------------------------
predict_test <- function(fit, test_pop) {
  p <- predict(fit, newdata = test_pop, type = "response")
  stopifnot(all(p > 0 & p < 1))
  p
}
