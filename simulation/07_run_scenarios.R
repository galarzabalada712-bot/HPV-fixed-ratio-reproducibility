# 07_run_scenarios.R ----------------------------------------------------------
# Protocol sections 13, 24 and 26: run all repetitions of one scenario
# (sequential inside the scenario), save raw results per scenario, and
# optionally parallelize ACROSS scenarios with future/furrr.

# Deterministic seed for the fixed test population, in a different seed space
# from training seeds (protocol section 7).
test_pop_seed <- function(master_seed, scenario_id) {
  master_seed + scenario_id * 10000000L + 999L
}

# Run all repetitions of ONE scenario against ONE fixed test population.
# Raw results are saved immediately to output/raw/<scenario_label>.rds
# (protocol section 26), before moving on to the next scenario.
run_one_scenario <- function(scenario, master_seed,
                             n_repetitions = NULL, n_test = NULL,
                             n_train_override = NULL,
                             raw_dir = sim_paths$raw, save = TRUE,
                             verbose = TRUE) {
  scen <- as.list(scenario)
  if (!is.null(n_repetitions)) scen$n_repetitions <- n_repetitions
  if (!is.null(n_test)) scen$n_test <- n_test
  if (!is.null(n_train_override)) scen$n_train <- n_train_override

  # One deterministic test population shared by both models across repetitions.
  test_seed <- test_pop_seed(master_seed, scen$scenario_id)
  test_pop <- generate_test_population(scen, test_seed)

  rows <- vector("list", scen$n_repetitions)
  n_fail <- 0L
  for (r in seq_len(scen$n_repetitions)) {
    row <- run_one_repetition(scen, r, test_pop, master_seed)
    n_fail <- n_fail + as.integer(row$fit_success == 0L)
    rows[[r]] <- row
  }
  results <- do.call(rbind, rows)

  out <- list(
    scenario = as.data.frame(scen, stringsAsFactors = FALSE),
    test_seed = test_seed,
    n_repetitions = scen$n_repetitions,
    results = results
  )

  if (save) {
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
    file <- file.path(raw_dir, paste0(scen$scenario_label, ".rds"))
    saveRDS(out, file)
  }

  if (verbose) {
    cat(sprintf("[%s] %d reps, %d successful, %d failed\n",
                scen$scenario_label, scen$n_repetitions,
                scen$n_repetitions - n_fail, n_fail))
  }
  out
}

# Run a set of scenarios (protocol section 26): parallelize across scenarios,
# keep repetitions sequential within a scenario so each worker only holds its
# own test population. If future/furrr are unavailable, run sequentially.
run_all_scenarios <- function(grid, master_seed, workers = 1L,
                              n_repetitions = NULL, n_test = NULL,
                              n_train_override = NULL,
                              scenario_ids = NULL,
                              raw_dir = sim_paths$raw, save = TRUE,
                              verbose = TRUE) {
  sub <- grid
  if (!is.null(scenario_ids)) {
    sub <- grid[grid$scenario_id %in% scenario_ids, , drop = FALSE]
  }
  if (nrow(sub) == 0L) stop("No scenarios selected.")

  scenario_list <- split(sub, seq_len(nrow(sub)))
  run_one <- function(sc) run_one_scenario(
    sc, master_seed,
    n_repetitions = n_repetitions, n_test = n_test,
    n_train_override = n_train_override,
    raw_dir = raw_dir, save = save, verbose = verbose
  )

  use_parallel <- workers > 1L &&
    requireNamespace("furrr", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE)

  if (use_parallel) {
    cat(sprintf("Parallel execution across %d scenarios with %d workers\n",
                length(scenario_list), workers))
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    out <- furrr::future_map(
      scenario_list, run_one,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    if (workers > 1L) {
      cat("future/furrr not available; falling back to sequential execution.\n")
    }
    out <- lapply(scenario_list, run_one)
  }
  invisible(out)
}
