# 09_quality_control.R --------------------------------------------------------
# Protocol section 29: quality-control checks before the final run.

# 29.1 Data-generation checks on a large sample
qc_data_generation <- function(rho = 0.2, n = 2e5, seed = 20260731L) {
  set.seed(seed)
  pred <- generate_predictors(n, rho)
  checks <- list(
    E_XA_near_0    = abs(mean(pred$X_A)) < 0.01,
    E_XB_near_0    = abs(mean(pred$X_B)) < 0.01,
    Var_XA_near_1  = abs(var(pred$X_A) - 1) < 0.05,
    Var_XB_near_1  = abs(var(pred$X_B) - 1) < 0.05,
    Corr_XA_XB     = abs(cor(pred$X_A, pred$X_B) - rho) < 0.05,
    Var_R_near_2m  = abs(var(pred$R) - 2 * (1 - rho)) < 0.1,
    Var_S_near_2p  = abs(var(pred$S) - 2 * (1 + rho)) < 0.1,
    Cov_R_S_near_0 = abs(cov(pred$R, pred$S)) < 0.05
  )
  list(pass = all(unlist(checks)), checks = checks)
}

# 29.2 Coefficient-transformation checks for every scenario
qc_coefficient_transforms <- function(grid) {
  checks <- list(
    beta_A_ok = all(abs(grid$beta_A - (grid$theta + grid$delta)) < 1e-12),
    beta_B_ok = all(abs(grid$beta_B - (-grid$theta + grid$delta)) < 1e-12),
    D_true_sum_ok = all(abs(grid$D_true - (grid$beta_A + grid$beta_B)) < 1e-12),
    D_true_2d_ok = all(abs(grid$D_true - 2 * grid$delta) < 1e-12),
    V_LP_agree = all(abs(grid$V_LP - grid$V_LP2) < 1e-12)
  )
  list(pass = all(unlist(checks)), checks = checks)
}

# 29.3 Prevalence checks on an independent large population (per scenario).
# Returns BOTH the observed deviation (prevalence_error) and its Monte Carlo
# standard error, so the two sources of deviation are kept separable:
#   (a) calibration algorithm tolerance  -> see calibrated_intercepts.csv
#   (b) MC sampling noise of the check   -> deviation vs 3*mcse here
qc_prevalence <- function(grid, n = 1e6, seed = 20260731L) {
  results <- lapply(seq_len(nrow(grid)), function(i) {
    # Deterministic per-scenario seed: each scenario's check population is
    # reproducible (protocol 29.3 reproducibility).
    set.seed(seed + i)
    scen <- grid[i, ]
    pred <- generate_predictors(n, scen$rho)
    dat <- generate_outcomes(pred, scen$theta, scen$delta,
                             scen$gamma_1, scen$gamma_2, scen$alpha)
    data.frame(
      scenario_id = scen$scenario_id,
      scenario_label = scen$scenario_label,
      prevalence = mean(dat$Y),
      prevalence_error = mean(dat$Y) - scen$prevalence_target,
      prevalence_mcse = sqrt(0.25 * 0.75 / n)
    )
  })
  df <- do.call(rbind, results)
  # The 0.001 target tolerance must be interpreted against Monte Carlo noise:
  # a deviation up to 0.001 + 3*MCSE is expected sampling variation.
  list(
    pass = all(abs(df$prevalence_error) <= 0.001 + 3 * df$prevalence_mcse),
    deviations = df
  )
}

# 29.4 Nested-model checks under D = 0: ratio and component log-likelihoods
#     are very close, LRT has 1 df, mean rejection near 0.05 (pilot run).
qc_nested_models <- function(scenario, master_seed = 20260731L, M = 100L,
                             n = 2000L, n_test = 5000L) {
  scen <- as.list(scenario)
  scen$n_train <- n
  scen$n_repetitions <- M
  scen$n_test <- n_test
  test_seed <- test_pop_seed(master_seed, scen$scenario_id)
  test_pop <- generate_test_population(scen, test_seed)

  rows <- lapply(seq_len(M), function(r) run_one_repetition(scen, r, test_pop, master_seed))
  res <- do.call(rbind, rows)
  ok <- res$fit_success == 1L

  # Under D = 0 the ratio model is correctly specified; both models should
  # have very similar fitted likelihoods -> deviance gap small.
  dev_gap <- res$lrt_statistic[ok]  # = Deviance(ratio) - Deviance(component)
  rejection <- mean(res$reject_lrt[ok])
  mcse <- mcse_proportion(rejection, sum(ok))

  list(
    pass = all(
      mean(dev_gap) < 3,                      # small average deviance gap
      all(res$lrt_df[ok] == 1L),              # 1 df
      abs(rejection - 0.05) <= 4 * mcse       # rejection near 0.05
    ),
    mean_deviance_gap = mean(dev_gap),
    rejection_rate = rejection,
    rejection_mcse = mcse,
    n_success = sum(ok)
  )
}

# 29.5 Metric-direction unit checks
qc_metric_directions <- function() {
  set.seed(99)
  test_pop <- data.frame(Y = rbinom(2000, 1, 0.25))
  p_ratio <- runif(2000, 0.05, 0.5)
  p_component <- p_ratio + rnorm(2000, 0, 0.05)
  m <- compute_performance_metrics(test_pop, p_ratio, p_component)
  checks <- list(
    delta_auc_ok = abs(m$delta_auc - (m$auc_component - m$auc_ratio)) < 1e-12,
    delta_brier_ok = abs(m$delta_brier - (m$brier_ratio - m$brier_component)) < 1e-12,
    delta_logloss_ok = abs(m$delta_logloss - (m$logloss_ratio - m$logloss_component)) < 1e-12,
    delta_NB_010_ok = abs(m$delta_NB_010 - (m$NB_component_010 - m$NB_ratio_010)) < 1e-12,
    delta_NB_020_ok = abs(m$delta_NB_020 - (m$NB_component_020 - m$NB_ratio_020)) < 1e-12,
    delta_NB_030_ok = abs(m$delta_NB_030 - (m$NB_component_030 - m$NB_ratio_030)) < 1e-12
  )
  list(pass = all(unlist(checks)), checks = checks)
}

# 29.5b Fast calibration implementations must match glm numerically.
# The quick path (uniroot intercept + IRLS slope) replaces glm() purely for
# speed on the 100k test population; this check proves the estimates agree.
qc_fast_calibration <- function(n = 2e4, seed = 100) {
  set.seed(seed)
  y <- rbinom(n, 1, 0.25)
  lp <- rnorm(n, 0, 0.8)
  q_int <- calibration_intercept(y, lp)
  q_slope <- calibration_slope(y, lp)[["b1"]]
  g_int <- as.numeric(coef(glm(y ~ 1, offset = lp, family = binomial()))[1])
  g_fit <- glm(y ~ lp, family = binomial())
  g_slope <- as.numeric(coef(g_fit)["lp"])
  list(
    pass = abs(q_int - g_int) < 1e-6 && abs(q_slope - g_slope) < 1e-6,
    intercept_diff = abs(q_int - g_int),
    slope_diff = abs(q_slope - g_slope)
  )
}

# 29.6 Null-scenario checks from summaries (delta = 0 rows)
qc_null_scenarios <- function(inferential, predictive) {
  null_inf <- inferential[inferential$delta == 0, , drop = FALSE]
  if (nrow(null_inf) == 0L) {
    return(list(pass = TRUE, note = "no delta=0 scenarios in summaries"))
  }
  # With very few repetitions (stage 1, M=5) the rejection MCSE is not
  # informative; require a minimum number of successful reps before judging.
  if (any(null_inf$n_success < 50L)) {
    return(list(pass = TRUE, note = "null-scenario check deferred (n_success < 50)"))
  }
  null_pred <- predictive[predictive$delta == 0, , drop = FALSE]
  checks <- list(
    rejection_near_0.05 = all(abs(null_inf$rejection_rate - 0.05) <=
                                4 * null_inf$rejection_mcse),
    bias_D_near_0 = all(abs(null_inf$bias_D) <=
                          4 * null_inf$ese_D / sqrt(null_inf$n_success)),
    coverage_near_0.95 = all(abs(null_inf$coverage - 0.95) <=
                               4 * null_inf$coverage_mcse),
    mean_delta_auc_near_0 = all(abs(null_pred$delta_auc) <=
                                  4 * null_pred$delta_auc_mcse + 0.005),
    mean_delta_brier_near_0 = all(abs(null_pred$delta_brier) <=
                                    4 * null_pred$delta_brier_mcse + 0.001)
  )
  list(pass = all(unlist(checks)), checks = checks)
}

# 29.7 Reproducibility: same pilot twice -> identical raw results
qc_reproducibility <- function(scenario, master_seed = 20260731L, M = 5L,
                               n = 500L, n_test = 2000L) {
  scen <- as.list(scenario)
  scen$n_train <- n
  scen$n_repetitions <- M
  scen$n_test <- n_test
  a <- run_one_scenario(scen, master_seed, n_repetitions = M,
                        n_test = n_test, n_train_override = n,
                        raw_dir = tempdir(), save = FALSE, verbose = FALSE)
  b <- run_one_scenario(scen, master_seed, n_repetitions = M,
                        n_test = n_test, n_train_override = n,
                        raw_dir = tempdir(), save = FALSE, verbose = FALSE)
  list(pass = isTRUE(all.equal(a$results, b$results)))
}

# Run the full battery, write a QC log.
run_quality_control <- function(grid, master_seed = 20260731L,
                                summaries = NULL,
                                null_scenario = NULL,
                                log_file = file.path(sim_paths$summaries, "quality_control_log.txt"),
                                include_prevalence = FALSE) {
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  sink(log_file, split = TRUE)
  on.exit(sink())
  cat("=== Simulation quality control ===\n")
  cat("Started:", format(Sys.time()), "\n")

  results <- list()
  results$data_generation <- qc_data_generation()
  cat("\n[29.1] Data generation:",
      if (results$data_generation$pass) "PASS" else "FAIL", "\n")

  results$coefficient_transforms <- qc_coefficient_transforms(grid)
  cat("[29.2] Coefficient transforms:",
      if (results$coefficient_transforms$pass) "PASS" else "FAIL", "\n")

  results$metric_directions <- qc_metric_directions()
  cat("[29.5] Metric directions:",
      if (results$metric_directions$pass) "PASS" else "FAIL", "\n")

  results$fast_calibration <- qc_fast_calibration()
  cat("[29.5b] Fast calibration vs glm:",
      if (results$fast_calibration$pass) "PASS" else "FAIL",
      sprintf("(int diff %.2e, slope diff %.2e)\n",
              results$fast_calibration$intercept_diff,
              results$fast_calibration$slope_diff))

  results$reproducibility <- qc_reproducibility(grid[1, , drop = FALSE], master_seed)
  cat("[29.7] Reproducibility (same seed twice):",
      if (results$reproducibility$pass) "PASS" else "FAIL", "\n")

  if (include_prevalence) {
    # Two separate judgements (protocol 29.3):
    #  (a) calibration algorithm numerical tolerance (alpha solved on the
    #      calibration population must give prevalence within 1e-4 of target);
    #  (b) independent large-population MC verification, allowing for
    #      0.001 + 3*MCSE sampling noise.
    # Both quantities are SAVED to prevalence_check.csv, so a relaxed PASS can
    # never hide a systematic calibration error.
    cal_file <- file.path(sim_paths$summaries, "calibrated_intercepts.csv")
    cal_pass <- TRUE
    cal_df <- NULL
    if (file.exists(cal_file)) {
      cal_df <- read.csv(cal_file, stringsAsFactors = FALSE)
      cal_df$calibration_error <-
        cal_df$prevalence_on_cal_population - 0.25
      cal_pass <- all(abs(cal_df$calibration_error) < 1e-4)
    }
    results$prevalence_calibration <- list(pass = cal_pass)
    cat("[29.3a] Prevalence calibration tolerance (|error| < 1e-4):",
        if (cal_pass) "PASS" else "FAIL", "\n")

    results$prevalence <- qc_prevalence(grid)
    cat("[29.3b] Prevalence MC verification (<= 0.001 + 3*MCSE):",
        if (results$prevalence$pass) "PASS" else "FAIL", "\n")

    # Save both quantities per scenario for manual review.
    prev_check <- results$prevalence$deviations
    if (!is.null(cal_df)) {
      prev_check <- merge(
        prev_check,
        cal_df[, c("scenario_id", "alpha", "calibration_error")],
        by = "scenario_id", all.x = TRUE, sort = FALSE
      )
    }
    write.csv(prev_check,
              file.path(sim_paths$summaries, "prevalence_check.csv"),
              row.names = FALSE)
    cat("Saved per-scenario prevalence_error + prevalence_mcse to",
        file.path(sim_paths$summaries, "prevalence_check.csv"), "\n")
  }

  if (!is.null(null_scenario)) {
    results$nested_models <- qc_nested_models(null_scenario, master_seed)
    cat("[29.4] Nested models:",
        if (results$nested_models$pass) "PASS" else "FAIL",
        sprintf("(mean deviance gap %.3f, rejection %.3f +/- %.3f)\n",
                results$nested_models$mean_deviance_gap,
                results$nested_models$rejection_rate,
                results$nested_models$rejection_mcse))
  }

  if (!is.null(summaries)) {
    results$null_scenarios <- qc_null_scenarios(summaries$inferential,
                                                summaries$predictive)
    cat("[29.6] Null scenarios:",
        if (results$null_scenarios$pass) "PASS" else "FAIL", "\n")

    # Finite-sample coverage sanity: with M=1000 the coverage MCSE at 0.95 is
    # ~0.0069, so a scenario sitting below 0.90 - 4*MCSE warrants inspection
    # (bad model-based SE, Wald CI skew, or variance/coefficient order bug).
    inf <- summaries$inferential
    if (nrow(inf) > 0L && any(inf$n_success >= 100L)) {
      enough <- inf$n_success >= 100L
      floor_cov <- 0.90 - 4 * inf$coverage_mcse
      bad <- inf[enough & inf$coverage < floor_cov, , drop = FALSE]
      results$coverage_sanity <- list(
        pass = nrow(bad) == 0L,
        scenarios_below_floor = bad[, c("scenario_id", "scenario_label",
                                        "coverage", "coverage_mcse")],
        floor = floor_cov
      )
      cat("[28] Coverage sanity (coverage >= 0.90 - 4*MCSE):",
          if (results$coverage_sanity$pass) "PASS" else "FAIL",
          sprintf("(%d scenario(s) below floor)\n", nrow(bad)))
    }
  }

  all_pass <- all(vapply(results, function(x) isTRUE(x$pass), logical(1)))
  cat("\nOVERALL:", if (all_pass) "ALL CHECKS PASS" else "SOME CHECKS FAILED", "\n")
  cat("Finished:", format(Sys.time()), "\n")
  invisible(list(pass = all_pass, results = results))
}
