# 11_diagnose_scenario.R ------------------------------------------------------
# Targeted diagnostic stage for the undercoverage observed in
# CORE_D0_N1000_RHO07 (null, n=1000, rho=0.7): rejection 0.079, coverage 0.921,
# SE ratio 0.935. Runs two scenarios x M repetitions and collects, per
# repetition:
#   - standardized statistic  Z_D = (D_hat - D_true) / se_D_hat
#   - Wald, likelihood-ratio (LRT) and score (Rao) tests for H0: D = 0
#   - theta/delta reparameterization (theta = (bA-bB)/2, delta = (bA+bB)/2)
#     with their model-based SEs
#   - vcov components Var(bA), Var(bB), Cov(bA, bB)
#   - condition number of the information matrix
#
# Key identity used for interval-construction diagnosis: for a null scenario,
# the profile-likelihood CI for D covers D_true = 0 if and only if the LRT
# does not reject, so profile coverage == 1 - LRT rejection rate.

run_diagnostic_repetition <- function(scenario, rep_id, test_pop, master_seed) {
  seed_id <- train_seed(master_seed, scenario$scenario_id, rep_id)
  set.seed(seed_id)
  row <- data.frame(
    repetition = rep_id, seed = seed_id,
    event_rate_train = NA_real_,
    b_A_hat = NA_real_, b_B_hat = NA_real_,
    theta_hat = NA_real_, delta_hat = NA_real_,
    D_hat = NA_real_, se_D_hat = NA_real_,
    var_bA = NA_real_, var_bB = NA_real_, cov_bAbB = NA_real_,
    se_theta = NA_real_, se_delta = NA_real_,
    z_D = NA_real_, wald_z = NA_real_,
    reject_wald = NA_integer_,
    lrt_stat = NA_real_, lrt_p = NA_real_, reject_lrt = NA_integer_,
    score_stat = NA_real_, score_p = NA_real_, reject_score = NA_integer_,
    kappa_info = NA_real_,
    converged = NA_integer_, fit_success = NA_integer_, failure_reason = "",
    stringsAsFactors = FALSE
  )
  tryCatch({
    train <- generate_train_data(scenario, seed_id)
    stopifnot(length(unique(train$Y)) == 2L)

    fit_ratio <- fit_ratio_model(train)
    fit_component <- fit_component_model(train)
    for (f in list(fit_ratio, fit_component)) {
      problem <- check_fit(f)
      if (!is.null(problem)) stop(problem)
    }

    # --- Three tests of H0: D = 0 ---
    lrt <- compute_lrt(fit_ratio, fit_component)
    rao <- anova(fit_ratio, fit_component, test = "Rao")
    score_stat <- rao$Rao[2]
    score_p <- rao[2, "Pr(>Chi)"]

    # --- Wald: D_hat / se ---
    D <- estimate_D(fit_component, scenario$D_true)
    wald_z <- D$D_hat / D$se_D_hat
    z_D <- (D$D_hat - scenario$D_true) / D$se_D_hat

    # --- theta/delta decomposition ---
    bA <- D$b_A_hat; bB <- D$b_B_hat
    theta_hat <- (bA - bB) / 2
    delta_hat <- (bA + bB) / 2
    V <- vcov(fit_component)
    var_bA <- V["X_A", "X_A"]; var_bB <- V["X_B", "X_B"]
    cov_bAbB <- V["X_A", "X_B"]
    se_theta <- 0.5 * sqrt(var_bA + var_bB - 2 * cov_bAbB)
    se_delta <- 0.5 * sqrt(var_bA + var_bB + 2 * cov_bAbB)

    row$event_rate_train <- mean(train$Y)
    row$b_A_hat <- bA; row$b_B_hat <- bB
    row$theta_hat <- theta_hat; row$delta_hat <- delta_hat
    row$D_hat <- D$D_hat; row$se_D_hat <- D$se_D_hat
    row$var_bA <- var_bA; row$var_bB <- var_bB; row$cov_bAbB <- cov_bAbB
    row$se_theta <- se_theta; row$se_delta <- se_delta
    row$z_D <- z_D; row$wald_z <- wald_z
    row$reject_wald <- as.integer(abs(wald_z) > qnorm(0.975))
    row$lrt_stat <- lrt$lrt_statistic; row$lrt_p <- lrt$lrt_p
    row$reject_lrt <- as.integer(lrt$lrt_p < 0.05)
    row$score_stat <- score_stat; row$score_p <- score_p
    row$reject_score <- as.integer(score_p < 0.05)
    row$kappa_info <- kappa(V)
    row$converged <- as.integer(fit_component$converged)
    row$fit_success <- 1L
  }, error = function(e) {
    row$fit_success <- 0L
    row$failure_reason <- conditionMessage(e)
  })
  row
}

# Run M repetitions of one scenario for diagnosis.
run_diagnostic_scenario <- function(scenario, master_seed, M = 5000L,
                                    n_test = 100000L, verbose = TRUE) {
  # CRITICAL: the main simulation runs under RNGkind("L'Ecuyer-CMRG") and the
  # per-repetition set.seed() calls produce different streams under the default
  # Mersenne-Twister. Repetition-level results are only comparable to Stage 4
  # raw data under the same RNG kind, so force it here as well.
  RNGkind("L'Ecuyer-CMRG")
  scen <- as.list(scenario)
  scen$n_test <- n_test
  test_pop <- generate_test_population(scen, test_pop_seed(master_seed, scen$scenario_id))
  rows <- vector("list", M)
  for (r in seq_len(M)) rows[[r]] <- run_diagnostic_repetition(scen, r, test_pop, master_seed)
  res <- do.call(rbind, rows)
  if (verbose) {
    cat(sprintf("[diag %s] %d reps, %d successful, %d failed\n",
                scen$scenario_label, M, sum(res$fit_success == 1L),
                sum(res$fit_success == 0L)))
  }
  res
}

# --- Scenario-level diagnostic summary ---------------------------------------
summarize_diagnostic <- function(res, scenario) {
  ok <- res$fit_success == 1L
  n_ok <- sum(ok)
  z <- res$z_D[ok]
  zc <- z - mean(z)
  skew <- mean(zc^3) / sd(z)^3
  ekurt <- mean(zc^4) / sd(z)^4 - 3

  mcse_prop <- function(p, m) sqrt(p * (1 - p) / m)
  wald_rej <- mean(res$reject_wald[ok])
  lrt_rej <- mean(res$reject_lrt[ok])
  score_rej <- mean(res$reject_score[ok])

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    n = scenario$n_train, rho = scenario$rho,
    n_success = n_ok, failure_rate = 1 - n_ok / nrow(res),
    # Z_D distribution
    z_mean = mean(z), z_sd = sd(z),
    z_skewness = skew, z_excess_kurtosis = ekurt,
    z_q2.5 = unname(quantile(z, 0.025)), z_q5 = unname(quantile(z, 0.05)),
    z_q95 = unname(quantile(z, 0.95)), z_q97.5 = unname(quantile(z, 0.975)),
    P_z_lt_neg1.96 = mean(z < -1.96), P_z_gt_1.96 = mean(z > 1.96),
    # three tests
    wald_rejection = wald_rej, wald_rejection_mcse = mcse_prop(wald_rej, n_ok),
    lrt_rejection = lrt_rej, lrt_rejection_mcse = mcse_prop(lrt_rej, n_ok),
    score_rejection = score_rej, score_rejection_mcse = mcse_prop(score_rej, n_ok),
    # D estimation
    bias_D = mean(res$D_hat[ok]) - scenario$D_true,
    ese_D = sd(res$D_hat[ok]), mean_se_D = mean(res$se_D_hat[ok]),
    se_ratio = mean(res$se_D_hat[ok]) / sd(res$D_hat[ok]),
    # vcov decomposition
    mean_var_bA = mean(res$var_bA[ok]), mean_var_bB = mean(res$var_bB[ok]),
    mean_cov_bAbB = mean(res$cov_bAbB[ok]),
    emp_var_D = var(res$D_hat[ok]),
    # theta/delta
    ese_theta = sd(res$theta_hat[ok]), mean_se_theta = mean(res$se_theta[ok]),
    ese_delta = sd(res$delta_hat[ok]), mean_se_delta = mean(res$se_delta[ok]),
    cor_theta_delta = cor(res$theta_hat[ok], res$delta_hat[ok]),
    kappa_mean = mean(res$kappa_info[ok]), kappa_p95 = unname(quantile(res$kappa_info[ok], 0.95)),
    stringsAsFactors = FALSE
  )
}

# --- Full diagnostic: abnormal scenario + matched control ---------------------
run_undercoverage_diagnosis <- function(grid, master_seed = 20260731L, M = 5000L,
                                        out_dir = sim_paths$summaries) {
  target <- grid[grid$scenario_label == "CORE_D0_N1000_RHO07", , drop = FALSE]
  control <- grid[grid$scenario_label == "CORE_D0_N1000_RHO02", , drop = FALSE]
  if (nrow(target) == 0L || nrow(control) == 0L) {
    stop("Diagnosis requires CORE_D0_N1000_RHO07 and CORE_D0_N1000_RHO02 in the grid")
  }
  cat("=== Undercoverage diagnostic: abnormal scenario vs matched control ===\n")
  cat("Target : CORE_D0_N1000_RHO07 (n=1000, rho=0.7, null)\n")
  cat("Control: CORE_D0_N1000_RHO02 (n=1000, rho=0.2, null)\n")

  res_t <- run_diagnostic_scenario(target[1, , drop = FALSE], master_seed, M = M)
  res_c <- run_diagnostic_scenario(control[1, , drop = FALSE], master_seed, M = M)

  summ_t <- summarize_diagnostic(res_t, target[1, , drop = FALSE])
  summ_c <- summarize_diagnostic(res_c, control[1, , drop = FALSE])
  summ <- rbind(summ_c, summ_t)
  rownames(summ) <- NULL

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(out_dir, "undercoverage_diagnosis.csv")
  write.csv(summ, outfile, row.names = FALSE)
  cat("Saved diagnosis to", outfile, "\n")

  print(summ, digits = 4)

  cat("\n--- Interpretation ---\n")
  for (i in seq_len(nrow(summ))) {
    s <- summ[i, ]
    cat(sprintf("\n[%s]  n_success=%d failure_rate=%.4f\n", s$scenario_label,
                s$n_success, s$failure_rate))
    cat(sprintf("  Z_D: mean=%.3f sd=%.3f skew=%.3f exkurt=%.3f  ",
                s$z_mean, s$z_sd, s$z_skewness, s$z_excess_kurtosis))
    cat(sprintf("P(Z<-1.96)=%.4f P(Z>1.96)=%.4f\n", s$P_z_lt_neg1.96, s$P_z_gt_1.96))
    cat(sprintf("  Rejection: Wald=%.4f LRT=%.4f Score=%.4f (profile coverage = 1 - LRT = %.4f)\n",
                s$wald_rejection, s$lrt_rejection, s$score_rejection,
                1 - s$lrt_rejection))
    cat(sprintf("  D: bias=%.4f ESE=%.4f meanSE=%.4f SEratio=%.3f\n",
                s$bias_D, s$ese_D, s$mean_se_D, s$se_ratio))
    cat(sprintf("  theta: ESE=%.4f meanSE=%.4f | delta: ESE=%.4f meanSE=%.4f | Cor(th,de)=%.3f\n",
                s$ese_theta, s$mean_se_theta, s$ese_delta, s$mean_se_delta,
                s$cor_theta_delta))
    cat(sprintf("  Var(bA)=%.5f Var(bB)=%.5f Cov=%.5f kappa mean=%.1f p95=%.1f\n",
                s$mean_var_bA, s$mean_var_bB, s$mean_cov_bAbB,
                s$kappa_mean, s$kappa_p95))
  }

  invisible(list(summary = summ, target_raw = res_t, control_raw = res_c))
}
