# 06_run_one_repetition.R -----------------------------------------------------
# Protocol sections 7, 14-25: run one repetition of one scenario and return a
# single row with every per-repetition output column listed in section 25.

# Deterministic seeds (protocol section 7):
#   train: seed_id = master_seed + scenario_id * 100000L + repetition_id
#   test population: a different deterministic seed space (see 07_run_scenarios.R)
train_seed <- function(master_seed, scenario_id, repetition_id) {
  master_seed + scenario_id * 100000L + repetition_id
}

# Empty row template: every output column of section 25 plus convergence flags.
blank_result_row <- function(scenario, repetition_id, seed_id) {
  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_label = scenario$scenario_label,
    scenario_group = scenario$scenario_group,
    baseline_strength = scenario$baseline_strength,
    theta = scenario$theta, delta = scenario$delta,
    D_true = scenario$D_true, V_LP = scenario$V_LP,
    rho = scenario$rho, n_train = scenario$n_train,
    repetition = repetition_id, seed = seed_id,
    event_rate_train = NA_real_, event_rate_test = NA_real_,
    lrt_statistic = NA_real_, lrt_df = NA_integer_, lrt_p = NA_real_,
    reject_lrt = NA_integer_,
    b_A_hat = NA_real_, b_B_hat = NA_real_, D_hat = NA_real_,
    se_D_hat = NA_real_, D_lower = NA_real_, D_upper = NA_real_,
    cover_D = NA_integer_,
    auc_covariate = NA_real_, auc_ratio = NA_real_, auc_component = NA_real_,
    delta_auc = NA_real_,
    brier_ratio = NA_real_, brier_component = NA_real_, delta_brier = NA_real_,
    logloss_ratio = NA_real_, logloss_component = NA_real_, delta_logloss = NA_real_,
    cal_intercept_ratio = NA_real_, cal_intercept_component = NA_real_,
    cal_slope_ratio = NA_real_, cal_slope_component = NA_real_,
    NB_ratio_010 = NA_real_, NB_component_010 = NA_real_, NB_all_010 = NA_real_,
    delta_NB_010 = NA_real_,
    NB_ratio_020 = NA_real_, NB_component_020 = NA_real_, NB_all_020 = NA_real_,
    delta_NB_020 = NA_real_,
    NB_ratio_030 = NA_real_, NB_component_030 = NA_real_, NB_all_030 = NA_real_,
    delta_NB_030 = NA_real_,
    component_clinically_useful_010 = NA_integer_, ratio_clinically_useful_010 = NA_integer_,
    component_clinically_useful_020 = NA_integer_, ratio_clinically_useful_020 = NA_integer_,
    component_clinically_useful_030 = NA_integer_, ratio_clinically_useful_030 = NA_integer_,
    ratio_converged = NA_integer_, component_converged = NA_integer_,
    covariate_converged = NA_integer_,
    fit_success = NA_integer_, failure_reason = "", error_message = "",
    warning_message = "",
    stringsAsFactors = FALSE
  )
}

# Run ONE repetition. Returns a one-row data.frame. Failures are RECORDED, not
# silently discarded (protocol section 24); they are never auto-replaced.
run_one_repetition <- function(scenario, repetition_id, test_pop, master_seed) {
  seed_id <- train_seed(master_seed, scenario$scenario_id, repetition_id)
  warnings_captured <- character()

  row <- tryCatch(
    withCallingHandlers(
      {
        set.seed(seed_id)                       # before generating the training sample
        train <- generate_train_data(scenario, seed_id)
        stopifnot(length(unique(train$Y)) == 2L) # both outcome classes present

        event_rate_train <- mean(train$Y)

        # --- Fit the three models ---
        fit_ratio <- fit_ratio_model(train)
        fit_component <- fit_component_model(train)
        fit_covariate <- fit_covariate_model(train)

        for (f in list(fit_ratio, fit_component, fit_covariate)) {
          problem <- check_fit(f)
          if (!is.null(problem)) stop(problem)
        }

        # --- Likelihood-ratio test (primary) ---
        lrt <- compute_lrt(fit_ratio, fit_component)

        # --- Estimate D and its uncertainty ---
        D <- estimate_D(fit_component, scenario$D_true)

        # --- Predictions in the fixed test population ---
        p_ratio <- predict_test(fit_ratio, test_pop)
        p_component <- predict_test(fit_component, test_pop)
        p_covariate <- predict_test(fit_covariate, test_pop)
        auc_covariate <- as.numeric(pROC::auc(
          response = test_pop$Y, predictor = p_covariate, quiet = TRUE
        ))

        # --- Performance metrics ---
        m <- compute_performance_metrics(test_pop, p_ratio, p_component)

        # --- Assemble the output row ---
        out <- blank_result_row(scenario, repetition_id, seed_id)
        out$event_rate_train <- event_rate_train
        out$event_rate_test <- mean(test_pop$Y)
        out$lrt_statistic <- lrt$lrt_statistic
        out$lrt_df <- lrt$lrt_df
        out$lrt_p <- lrt$lrt_p
        out$reject_lrt <- as.integer(lrt$lrt_p < 0.05)
        out$b_A_hat <- D$b_A_hat; out$b_B_hat <- D$b_B_hat
        out$D_hat <- D$D_hat; out$se_D_hat <- D$se_D_hat
        out$D_lower <- D$D_lower; out$D_upper <- D$D_upper
        out$cover_D <- D$cover_D
        out$auc_covariate <- auc_covariate
        for (nm in names(m)) out[[nm]] <- m[[nm]]
        out$ratio_converged <- as.integer(fit_ratio$converged)
        out$component_converged <- as.integer(fit_component$converged)
        out$covariate_converged <- as.integer(fit_covariate$converged)
        out$fit_success <- 1L
        out$failure_reason <- ""
        out
      },
      warning = function(w) {
        warnings_captured <<- c(warnings_captured, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      out <- blank_result_row(scenario, repetition_id, seed_id)
      out$fit_success <- 0L
      out$failure_reason <- conditionMessage(e)
      out$error_message <- conditionMessage(e)
      out
    }
  )
  if (length(warnings_captured) > 0) {
    row$warning_message <- paste(unique(warnings_captured), collapse = " | ")
  }
  row
}
