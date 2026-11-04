# 08_summarize_results.R ------------------------------------------------------
# Protocol sections 27 and 28: scenario-level summaries of inferential and
# predictive performance, with Monte Carlo standard errors.

# Data dictionary for the 64-column repetition-level output (protocol §25 and
# §34): every column name, its meaning, and its direction convention.
column_dictionary <- function() {
  data.frame(
    column = c(
      # scenario metadata
      "scenario_id", "scenario_label", "scenario_group", "baseline_strength",
      "theta", "delta", "D_true", "V_LP", "rho", "n_train",
      "repetition", "seed", "event_rate_train", "event_rate_test",
      # LRT
      "lrt_statistic", "lrt_df", "lrt_p", "reject_lrt",
      # D estimate
      "b_A_hat", "b_B_hat", "D_hat", "se_D_hat", "D_lower", "D_upper", "cover_D",
      # predictions
      "auc_covariate", "auc_ratio", "auc_component", "delta_auc",
      "brier_ratio", "brier_component", "delta_brier",
      "logloss_ratio", "logloss_component", "delta_logloss",
      "cal_intercept_ratio", "cal_intercept_component",
      "cal_slope_ratio", "cal_slope_component",
      # net benefit
      "NB_ratio_010", "NB_component_010", "NB_all_010", "delta_NB_010",
      "NB_ratio_020", "NB_component_020", "NB_all_020", "delta_NB_020",
      "NB_ratio_030", "NB_component_030", "NB_all_030", "delta_NB_030",
      # clinically useful flags
      "component_clinically_useful_010", "ratio_clinically_useful_010",
      "component_clinically_useful_020", "ratio_clinically_useful_020",
      "component_clinically_useful_030", "ratio_clinically_useful_030",
      # fit status
      "ratio_converged", "component_converged", "covariate_converged",
      "fit_success", "failure_reason", "error_message", "warning_message"
    ),
    meaning = c(
      "scenario index (1..24)", "scenario label, e.g. CORE_D020_N1000_RHO02",
      "CORE / MECH / STRONG", "moderate / strong baseline strength",
      "theta: ratio-direction coefficient", "delta: sum-direction (violation) coefficient",
      "D_true = beta_A + beta_B = 2*delta (formal constraint-violation parameter)",
      "standardized violation magnitude = |delta|*sqrt(2*(1+rho)) = |D|*sqrt((1+rho)/2)",
      "correlation of X_A, X_B", "training sample size",
      "repetition index", "repetition seed = master + id*1e5 + rep",
      "observed event rate in training sample", "observed event rate in test population",
      "LRT statistic (deviance difference)", "LRT df (must be 1)", "LRT p-value",
      "1 if lrt_p < 0.05 (reject fixed-ratio constraint)",
      "estimated beta_A (component model)", "estimated beta_B (component model)",
      "D_hat = b_A_hat + b_B_hat", "model-based SE of D_hat",
      "Wald CI lower bound for D", "Wald CI upper bound for D",
      "1 if CI covers D_true",
      "AUC of covariate-only model (Z1+Z2)", "AUC of fixed-ratio model",
      "AUC of unconstrained-component model",
      "auc_component - auc_ratio (positive = component better)",
      "Brier of ratio model", "Brier of component model",
      "brier_ratio - brier_component (positive = component better)",
      "log loss of ratio model (clipped 1e-15)", "log loss of component model",
      "logloss_ratio - logloss_component (positive = component better)",
      "calibration intercept of ratio model (ideal 0)",
      "calibration intercept of component model (ideal 0)",
      "calibration slope of ratio model (ideal 1)",
      "calibration slope of component model (ideal 1)",
      "net benefit of ratio model at threshold 0.10",
      "net benefit of component model at threshold 0.10",
      "net benefit of treat-all at threshold 0.10",
      "NB_component - NB_ratio at threshold 0.10 (positive = component better)",
      "same at threshold 0.20",
      "same at threshold 0.20", "same at threshold 0.20", "same at threshold 0.20",
      "same at threshold 0.30",
      "same at threshold 0.30", "same at threshold 0.30", "same at threshold 0.30",
      "1 if component NB > max(NB_all, NB_none) at 0.10",
      "1 if ratio NB > max(NB_all, NB_none) at 0.10",
      "same at 0.20", "same at 0.20",
      "same at 0.30", "same at 0.30",
      "glm convergence flag, ratio model", "glm convergence flag, component model",
      "glm convergence flag, covariate model",
      "1 = repetition succeeded, 0 = failed (never auto-replaced)",
      "failure reason (empty if success)", "error message (empty if success)",
      "captured glm warnings (empty if none)"
    ),
    stringsAsFactors = FALSE
  )
}

write_column_dictionary <- function(out_dir = sim_paths$summaries) {
  dict <- column_dictionary()
  write.csv(dict, file.path(out_dir, "column_dictionary.csv"), row.names = FALSE)
  invisible(dict)
}

mcse_proportion <- function(p, m) sqrt(p * (1 - p) / m)

summarize_one_scenario <- function(out) {
  scen <- out$scenario
  res <- out$results
  M <- nrow(res)
  ok <- res$fit_success == 1L
  n_ok <- sum(ok)

  # Inferential performance (section 27)
  rejection_rate <- mean(res$reject_lrt[ok])
  bias_D <- mean(res$D_hat[ok]) - scen$D_true
  ese_D <- sd(res$D_hat[ok])
  mean_model_se <- mean(res$se_D_hat[ok])
  se_ratio <- mean_model_se / ese_D
  coverage <- mean(res$cover_D[ok])

  # Predictive performance: mean, MCSE, SD, median, 2.5/97.5 pctiles (section 28)
  predictive_cols <- c(
    "auc_covariate", "auc_ratio", "auc_component", "delta_auc",
    "brier_ratio", "brier_component", "delta_brier",
    "logloss_ratio", "logloss_component", "delta_logloss",
    "cal_intercept_ratio", "cal_intercept_component",
    "cal_slope_ratio", "cal_slope_component",
    "NB_ratio_010", "NB_component_010", "NB_all_010", "delta_NB_010",
    "NB_ratio_020", "NB_component_020", "NB_all_020", "delta_NB_020",
    "NB_ratio_030", "NB_component_030", "NB_all_030", "delta_NB_030"
  )

  summ <- data.frame(
    scenario_id = scen$scenario_id,
    scenario_label = scen$scenario_label,
    scenario_group = scen$scenario_group,
    mechanism = scen$mechanism,
    baseline_strength = scen$baseline_strength,
    theta = scen$theta, delta = scen$delta, D_true = scen$D_true,
    # V_LP is recomputed from the raw data-generating parameters (protocol
    # section 12: |delta|*sqrt(2*(1+rho))), so a formula fix in the scenario
    # grid propagates to summaries without re-running the simulation.
    V_LP = abs(scen$delta) * sqrt(2 * (1 + scen$rho)),
    rho = scen$rho, n_train = scen$n_train,
    n_repetitions = M, n_success = n_ok,
    failure_rate = 1 - n_ok / M,
    rejection_rate = rejection_rate,
    rejection_mcse = mcse_proportion(rejection_rate, n_ok),
    bias_D = bias_D, ese_D = ese_D,
    mean_model_se = mean_model_se, se_ratio = se_ratio,
    coverage = coverage,
    coverage_mcse = mcse_proportion(coverage, n_ok),
    stringsAsFactors = FALSE
  )

  for (col in predictive_cols) {
    vals <- res[[col]][ok]
    summ[[col]] <- mean(vals)
    summ[[paste0(col, "_mcse")]] <- sd(vals) / sqrt(n_ok)
    summ[[paste0(col, "_sd")]] <- sd(vals)
    summ[[paste0(col, "_median")]] <- median(vals)
    summ[[paste0(col, "_p2.5")]] <- unname(quantile(vals, 0.025))
    summ[[paste0(col, "_p97.5")]] <- unname(quantile(vals, 0.975))
  }
  summ
}

# --- Centralized run log (protocol section 34) ----------------------------------
# Plain-text log with start/end time, machine specs, R and package versions,
# master seed, worker count, successful/failed repetitions and a reference to
# the output manifest (hashes).
write_run_log <- function(master_seed, stage, workers, M, n_test,
                          summaries = NULL, raw_dir = sim_paths$raw,
                          out_dir = sim_paths$summaries, start_time = NULL) {
  end_time <- Sys.time()
  lines <- c(
    "=== Simulation run log ===",
    sprintf("Start: %s", if (is.null(start_time)) "unknown" else format(start_time)),
    sprintf("End:   %s", format(end_time)),
    sprintf("Machine: %s", paste(Sys.info()[c("sysname", "release", "machine")],
                                 collapse = " ")),
    sprintf("R version: %s", R.version.string),
    sprintf("Master seed: %d", master_seed),
    sprintf("Stage: %d", stage),
    sprintf("Workers: %d", workers),
    sprintf("Repetitions per scenario (M): %d", M),
    sprintf("Test population size (n_test): %d", n_test),
    "Packages:"
  )
  si <- sessionInfo()
  if (!is.null(si$otherPkgs)) {
    for (p in names(si$otherPkgs)) {
      lines <- c(lines, sprintf("  %s %s", p, si$otherPkgs[[p]]$Version))
    }
  }
  if (!is.null(summaries)) {
    inf <- summaries$inferential
    lines <- c(lines,
      sprintf("Scenarios run: %d", nrow(inf)),
      sprintf("Successful repetitions: %d", sum(inf$n_success)),
      sprintf("Failed repetitions: %d", sum(inf$n_repetitions) - sum(inf$n_success))
    )
  }
  files <- list.files(raw_dir, pattern = "\\.rds$")
  lines <- c(lines,
    sprintf("Raw result files: %d", length(files)),
    "Output file hashes: see output_manifest.csv (regenerate with md5sum over simulation/)")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(out_dir, "run_log.txt"))
  invisible(lines)
}

# Aggregate all raw scenario files into inferential and predictive summaries.
aggregate_results <- function(raw_dir = sim_paths$raw,
                              summary_dir = sim_paths$summaries) {
  files <- list.files(raw_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) stop("No raw result files found in ", raw_dir)
  cat("Aggregating", length(files), "scenario result files\n")
  summaries <- lapply(files, function(f) summarize_one_scenario(readRDS(f)))
  summ <- do.call(rbind, summaries)
  summ <- summ[order(summ$scenario_id), , drop = FALSE]

  inferential_cols <- c(
    "scenario_id", "scenario_label", "scenario_group", "mechanism",
    "baseline_strength",
    "theta", "delta", "D_true", "V_LP", "rho", "n_train",
    "n_repetitions", "n_success", "failure_rate",
    "rejection_rate", "rejection_mcse",
    "bias_D", "ese_D", "mean_model_se", "se_ratio",
    "coverage", "coverage_mcse"
  )
  predictive_cols <- c(
    "scenario_id", "scenario_label", "scenario_group", "mechanism",
    "baseline_strength",
    "delta", "D_true", "V_LP", "rho", "n_train",
    "auc_covariate", "auc_ratio", "auc_component", "delta_auc",
    "delta_auc_mcse", "delta_auc_sd", "delta_auc_median",
    "delta_auc_p2.5", "delta_auc_p97.5",
    "delta_brier", "delta_brier_mcse", "delta_brier_sd",
    "delta_brier_median", "delta_brier_p2.5", "delta_brier_p97.5",
    "delta_logloss", "delta_logloss_mcse", "delta_logloss_sd",
    "delta_logloss_median", "delta_logloss_p2.5", "delta_logloss_p97.5",
    "cal_intercept_ratio", "cal_intercept_component",
    "cal_slope_ratio", "cal_slope_component",
    "NB_ratio_010", "NB_component_010", "NB_all_010", "delta_NB_010",
    "delta_NB_010_mcse", "delta_NB_010_sd", "delta_NB_010_median",
    "delta_NB_010_p2.5", "delta_NB_010_p97.5",
    "NB_ratio_020", "NB_component_020", "NB_all_020", "delta_NB_020",
    "delta_NB_020_mcse", "delta_NB_020_sd", "delta_NB_020_median",
    "delta_NB_020_p2.5", "delta_NB_020_p97.5",
    "NB_ratio_030", "NB_component_030", "NB_all_030", "delta_NB_030",
    "delta_NB_030_mcse", "delta_NB_030_sd", "delta_NB_030_median",
    "delta_NB_030_p2.5", "delta_NB_030_p97.5"
  )

  inferential <- summ[, inferential_cols, drop = FALSE]
  predictive <- summ[, predictive_cols, drop = FALSE]

  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(inferential, file.path(summary_dir, "summary_inferential.csv"),
            row.names = FALSE)
  write.csv(predictive, file.path(summary_dir, "summary_predictive.csv"),
            row.names = FALSE)
  cat("Saved summaries to", summary_dir, "\n")
  list(inferential = inferential, predictive = predictive)
}
