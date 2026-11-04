# 02_calibrate_parameters.R ---------------------------------------------------
# Protocol sections 10 and 11: calibrate baseline prognostic strength (gamma)
# and the per-scenario intercept (alpha) on large reference populations.

# --- Calibrate gamma so the covariate-only AUC hits a target in the
#     reference scenario (delta = 0, rho = 0.2, theta = 0.5).
calibrate_gamma <- function(target_auc, theta = 0.5, delta = 0, rho = 0.2,
                            n_cal = 5e5, master_seed = 20260731L,
                            gamma_grid = seq(0.1, 2.0, by = 0.02)) {
  set.seed(master_seed + 1L)          # deterministic calibration population
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
  x  <- MASS::mvrnorm(n_cal, c(0, 0), Sigma)
  XA <- x[, 1]; XB <- x[, 2]
  R  <- XA - XB
  Z1 <- rnorm(n_cal); Z2 <- rnorm(n_cal)

  auc_for_gamma <- function(g) {
    # Complete data-generating model (alpha irrelevant for AUC);
    # covariate-only score is gamma*Z1 + gamma*Z2.
    eta <- theta * R + delta * (XA + XB) + g * Z1 + g * Z2
    p   <- plogis(eta)
    Y   <- rbinom(n_cal, 1, p)
    fast_auc(g * (Z1 + Z2), Y)
  }

  # Coarse pass, then a fine pass around the best candidate.
  coarse_auc <- vapply(gamma_grid, auc_for_gamma, numeric(1))
  best_idx   <- which.min(abs(coarse_auc - target_auc))[1]
  best_g     <- gamma_grid[best_idx]
  fine_grid  <- seq(max(0.01, best_g - 0.02), best_g + 0.02, by = 0.002)
  fine_auc   <- vapply(fine_grid, auc_for_gamma, numeric(1))
  best_idx   <- which.min(abs(fine_auc - target_auc))[1]

  list(gamma = fine_grid[best_idx], auc = fine_auc[best_idx], target = target_auc)
}

# --- Calibrate the intercept alpha so marginal prevalence = 0.25 -------------
calibrate_intercept <- function(theta, delta, gamma_1, gamma_2, rho,
                                target_prev = 0.25, n_cal = 5e5, seed = 20260731L) {
  set.seed(seed)                     # deterministic per-scenario population
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
  x  <- MASS::mvrnorm(n_cal, c(0, 0), Sigma)
  R  <- x[, 1] - x[, 2]
  S  <- x[, 1] + x[, 2]
  Z1 <- rnorm(n_cal); Z2 <- rnorm(n_cal)

  # Linear predictor excluding the intercept (protocol section 11).
  eta_no_intercept <- theta * R + delta * S + gamma_1 * Z1 + gamma_2 * Z2

  target_function <- function(alpha, eta, target_prev) {
    mean(plogis(alpha + eta)) - target_prev
  }
  alpha <- uniroot(
    target_function,
    interval = c(-15, 15),
    eta = eta_no_intercept,
    target_prev = target_prev
  )$root

  prev_check <- mean(plogis(alpha + eta_no_intercept))
  list(alpha = alpha, prevalence_on_cal_population = prev_check,
       diff_from_target = prev_check - target_prev)
}

# --- Calibrate everything, with caching --------------------------------------
calibrate_parameters <- function(grid, master_seed = 20260731L, n_cal = 5e5,
                                 force = FALSE) {
  dir.create(sim_paths$summaries, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(sim_paths$summaries, "scenario_grid_calibrated.csv")

  # Cache key must cover every mechanism parameter that affects calibration:
  # theta, delta (hence beta_A/beta_B and D_true), rho, baseline setting
  # (which fixes gamma), mechanism, and the prevalence target. alpha does not
  # depend on n_train, so the training size is intentionally excluded.
  calib_key_cols <- c("scenario_group", "mechanism", "baseline_strength",
                      "theta", "delta", "beta_A", "beta_B", "D_true",
                      "rho", "prevalence_target")
  grid_key <- apply(grid[, calib_key_cols, drop = FALSE], 1,
                    paste, collapse = "|")

  if (!force && file.exists(cache_file)) {
    cached <- read.csv(cache_file, stringsAsFactors = FALSE)
    if (nrow(cached) == nrow(grid) &&
        all(cached$scenario_id == grid$scenario_id) &&
        "cache_key" %in% names(cached) &&
        all(cached$cache_key == grid_key) &&
        all(c("gamma_1", "gamma_2", "alpha") %in% names(cached))) {
      cat("Loaded cached calibration from", cache_file, "\n")
      grid$gamma_1 <- cached$gamma_1
      grid$gamma_2 <- cached$gamma_2
      grid$alpha   <- cached$alpha
      # Refresh derived columns (e.g. V_LP, whose definition changed) in the
      # cache file so the persisted grid is always consistent with the code.
      grid$cache_key <- grid_key
      write.csv(grid, cache_file, row.names = FALSE)
      grid$cache_key <- NULL
      return(grid)
    }
    cat("Cache mismatch (mechanism parameters changed) or incomplete; recalibrating.\n")
  }

  cat("Calibrating moderate baseline strength (covariate-only AUC ~ 0.72)...\n")
  gamma_moderate <- calibrate_gamma(0.72, master_seed = master_seed, n_cal = n_cal)
  cat(sprintf("  gamma_moderate = %.4f (achieved AUC %.4f)\n",
              gamma_moderate$gamma, gamma_moderate$auc))
  cat("Calibrating strong baseline strength (covariate-only AUC ~ 0.85)...\n")
  gamma_strong <- calibrate_gamma(0.85, master_seed = master_seed, n_cal = n_cal)
  cat(sprintf("  gamma_strong   = %.4f (achieved AUC %.4f)\n",
              gamma_strong$gamma, gamma_strong$auc))

  gamma_values <- data.frame(
    baseline_strength = c("moderate", "strong"),
    target_auc        = c(0.72, 0.85),
    gamma             = c(gamma_moderate$gamma, gamma_strong$gamma),
    achieved_auc      = c(gamma_moderate$auc, gamma_strong$auc)
  )
  write.csv(gamma_values, file.path(sim_paths$summaries, "calibrated_gamma_values.csv"),
            row.names = FALSE)

  # Same gamma for all core + mechanistic scenarios; strong gamma for STRONG.
  grid$gamma_1 <- ifelse(grid$baseline_strength == "strong",
                         gamma_strong$gamma, gamma_moderate$gamma)
  grid$gamma_2 <- grid$gamma_1

  cat("Calibrating per-scenario intercepts (prevalence ~ 0.25)...\n")
  cal_list <- lapply(seq_len(nrow(grid)), function(i) {
    s <- grid[i, ]
    cal <- calibrate_intercept(
      s$theta, s$delta, s$gamma_1, s$gamma_2, s$rho,
      target_prev = s$prevalence_target, n_cal = n_cal,
      seed = master_seed + s$scenario_id * 12345L
    )
    if (abs(cal$diff_from_target) > 1e-4) {
      warning(sprintf("Intercept %s deviates from target by %.6f",
                      s$scenario_label, cal$diff_from_target))
    }
    cat(sprintf("  [%02d] %-24s alpha = % .4f  (prev %.6f)\n",
                s$scenario_id, s$scenario_label, cal$alpha,
                cal$prevalence_on_cal_population))
    data.frame(scenario_id = s$scenario_id, scenario_label = s$scenario_label,
               alpha = cal$alpha,
               prevalence_on_cal_population = cal$prevalence_on_cal_population)
  })
  intercept_df <- do.call(rbind, cal_list)
  write.csv(intercept_df, file.path(sim_paths$summaries, "calibrated_intercepts.csv"),
            row.names = FALSE)

  grid$alpha <- intercept_df$alpha[match(grid$scenario_id, intercept_df$scenario_id)]
  grid$cache_key <- grid_key
  write.csv(grid, cache_file, row.names = FALSE)
  grid$cache_key <- NULL
  cat("Calibration complete; saved scenario grid with parameters.\n")
  grid
}
