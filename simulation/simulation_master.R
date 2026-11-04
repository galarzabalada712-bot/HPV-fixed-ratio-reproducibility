#!/usr/bin/env Rscript
# simulation_master.R ---------------------------------------------------------
# Master orchestrator for the simulation study.
#
# Usage (run from the project root):
#   Rscript simulation/simulation_master.R --stage 1              # function testing
#   Rscript simulation/simulation_master.R --stage 2              # pilot validation
#   Rscript simulation/simulation_master.R --stage 3              # reduced full grid (M=100)
#   Rscript simulation/simulation_master.R --stage 4              # final run (M=1000)
#
# Options:
#   --workers N        parallel workers across scenarios (default 1)
#   --M N              override number of repetitions per scenario
#   --n_test N         override test population size
#   --scenario ids     run only the given scenario ids (comma separated)
#   --n_train N        override training size (used by stage 1: n=500)
#   --force-calibrate  recompute gamma/alpha calibration
#   --qc               run QC battery (incl. prevalence check)
#   --no-summarize     skip aggregation after running
#   --no-tables        skip tables/figures
#
# Protocol section 30 (staged execution) presets:
#   Stage 1: 1 scenario, n=500, 5 reps, n_test=5000
#   Stage 2: 4 core scenarios (D=0/0.4 n=1000 rho=0.2, D=0/0.8 n=5000 rho=0.7),
#            M=100, n_test=20000
#   Stage 3: all 24 scenarios, M=100, n_test=100000
#   Stage 4: all 24 scenarios, M=1000, n_test=100000
#   Stage 5: independent verification — selected null/mild scenarios with a
#            DIFFERENT master seed and independently generated test
#            populations (protocol section 30, stage 5).

# --- Command-line arguments --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  idx <- match(name, args)
  if (is.na(idx) || idx >= length(args)) default else args[idx + 1]
}
has_arg <- function(name) name %in% args

stage          <- as.integer(get_arg("--stage", 1L))
workers        <- as.integer(get_arg("--workers", 1L))
M_override     <- as.integer(get_arg("--M", NA_integer_))
n_test_override <- as.integer(get_arg("--n_test", NA_integer_))
scenario_ids   <- get_arg("--scenario", NULL)
n_train_override <- as.integer(get_arg("--n_train", NA_integer_))
force_cal      <- has_arg("--force-calibrate")
run_qc         <- has_arg("--qc")
no_summarize   <- has_arg("--no-summarize")
no_tables      <- has_arg("--no-tables")

master_seed <- 20260731L
run_start <- Sys.time()

# --- Source modules -----------------------------------------------------------
sim_dir <- dirname(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]))
if (length(sim_dir) == 0 || sim_dir == ".") sim_dir <- "simulation"
source(file.path(sim_dir, "00_packages.R"))   # defines sim_paths
sim_paths$root <- sim_dir
sim_paths$raw <- file.path(sim_dir, "output", "raw")
sim_paths$summaries <- file.path(sim_dir, "output", "summaries")
sim_paths$tables <- file.path(sim_dir, "output", "tables")
sim_paths$figures <- file.path(sim_dir, "output", "figures")

modules <- c(
  "00_packages.R", "01_scenario_grid.R", "02_calibrate_parameters.R",
  "03_generate_data.R", "04_fit_models.R", "05_performance_metrics.R",
  "06_run_one_repetition.R", "07_run_scenarios.R", "08_summarize_results.R",
  "09_quality_control.R", "10_tables_figures.R"
)
for (m in modules) source(file.path(sim_dir, m))

# --- Load packages (protocol section 6) ---------------------------------------
load_simulation_packages()

# --- RNG setup (protocol section 7) -------------------------------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(master_seed)

# --- Scenario grid and calibration --------------------------------------------
grid <- build_scenario_grid()
grid <- calibrate_parameters(grid, master_seed, force = force_cal)
cat("Scenario grid:", nrow(grid), "scenarios\n")

# --- Stage presets (protocol section 30) --------------------------------------
stage1_ids <- grid$scenario_id[grid$scenario_label == "CORE_D0_N1000_RHO02"]
stage2_labels <- c("CORE_D0_N1000_RHO02", "CORE_D0_N5000_RHO07",
                   "CORE_D020_N1000_RHO02", "CORE_D040_N5000_RHO07")
stage2_ids <- grid$scenario_id[grid$scenario_label %in% stage2_labels]

preset <- switch(
  as.character(stage),
  "1" = list(scenario_ids = stage1_ids, M = 5L, n_test = 5000L, n_train = 500L),
  "2" = list(scenario_ids = stage2_ids, M = 100L, n_test = 20000L, n_train = NULL),
  "3" = list(scenario_ids = NULL,       M = 100L, n_test = 100000L, n_train = NULL),
  "4" = list(scenario_ids = NULL,       M = 1000L, n_test = 100000L, n_train = NULL),
  "5" = list(scenario_ids = NULL,       M = 2000L, n_test = 100000L, n_train = NULL),
  stop("Unknown --stage: ", stage)
)

# Stage 5 uses a DIFFERENT master seed and its own output directories, so the
# verification stream/test populations are independent of Stage 4, and its
# run log / session info / summaries never overwrite the Stage 4 artifacts.
stage5_labels <- c("CORE_D0_N1000_RHO02", "CORE_D0_N1000_RHO07",
                   "CORE_D005_N1000_RHO02", "CORE_D005_N5000_RHO07")
stage5_seed <- 20260801L
if (stage == 5L) {
  preset$scenario_ids <- grid$scenario_id[grid$scenario_label %in% stage5_labels]
  sim_paths$raw <- file.path(sim_dir, "output", "raw_verify")
  sim_paths$summaries <- file.path(sim_dir, "output", "verify")
  sim_paths$tables <- file.path(sim_dir, "output", "verify")
  sim_paths$figures <- file.path(sim_dir, "output", "verify")
  master_seed <- stage5_seed
}

if (!is.null(scenario_ids)) {
  preset$scenario_ids <- as.integer(strsplit(scenario_ids, ",")[[1]])
}
M <- if (!is.na(M_override)) M_override else preset$M
n_test <- if (!is.na(n_test_override)) n_test_override else preset$n_test
n_train <- if (!is.na(n_train_override)) n_train_override else preset$n_train

# --- Run the selected scenarios ------------------------------------------------
cat(sprintf("Stage %d: M=%d, n_test=%d, workers=%d%s\n",
            stage, M, n_test, workers,
            if (stage == 5L) sprintf(" (master seed %d, output: raw_verify/)", stage5_seed) else ""))
if (is.null(preset$scenario_ids)) {
  cat("Scenarios: all", nrow(grid), "\n")
} else {
  cat("Scenarios:", paste(grid$scenario_label[grid$scenario_id %in% preset$scenario_ids],
                          collapse = ", "), "\n")
}

run_all_scenarios(
  grid,
  master_seed = master_seed,
  workers = workers,
  n_repetitions = M,
  n_test = n_test,
  n_train_override = n_train,
  scenario_ids = preset$scenario_ids
)

# --- Summaries -----------------------------------------------------------------
if (!no_summarize) {
  summaries <- aggregate_results()
  cat("\n--- Main summary (inferential) ---\n")
  print(summaries$inferential[, c("scenario_id", "scenario_label", "delta",
                                  "rejection_rate", "bias_D", "coverage")])
}

# --- Quality control ------------------------------------------------------------
if (run_qc) {
  null_scenario <- if (any(grid$delta == 0)) {
    grid[grid$delta == 0, , drop = FALSE][1, , drop = FALSE]
  } else NULL
  qc <- run_quality_control(
    grid, master_seed,
    summaries = if (!no_summarize) summaries else NULL,
    null_scenario = null_scenario,
    include_prevalence = TRUE
  )
  if (!qc$pass) {
    warning("Quality control reported failures - inspect output/summaries/quality_control_log.txt")
  }
}

# --- Tables and figures ----------------------------------------------------------
if (!no_tables && !no_summarize) {
  make_main_table(summaries$inferential, summaries$predictive)
  make_supplementary_complete(summaries$inferential, summaries$predictive)
  make_main_figure(summaries$inferential, summaries$predictive)
  make_supplementary_figure(summaries$inferential, summaries$predictive)
  make_mechanistic_figure(summaries$inferential, summaries$predictive)
}

# --- Data dictionary (protocol §34: reproducibility outputs) ----------------------
write_column_dictionary()
cat("Saved column dictionary to", file.path(sim_paths$summaries, "column_dictionary.csv"), "\n")

# --- Reproducibility outputs (protocol section 34) --------------------------------
dir.create(sim_paths$summaries, recursive = TRUE, showWarnings = FALSE)
sink(file.path(sim_paths$summaries, "session_info.txt"))
cat("R version:\n"); print(R.version.string)
cat("\nPackages:\n"); print(sessionInfo())
sink()
cat("Saved sessionInfo to", file.path(sim_paths$summaries, "session_info.txt"), "\n")

# --- Centralized run log (protocol section 34) ------------------------------------
write_run_log(
  master_seed = master_seed, stage = stage, workers = workers,
  M = M, n_test = n_test,
  summaries = if (!no_summarize) summaries else NULL,
  raw_dir = sim_paths$raw, out_dir = sim_paths$summaries,
  start_time = run_start
)
cat("Saved run log to", file.path(sim_paths$summaries, "run_log.txt"), "\n")

cat("\nDone. Raw results: ", sim_paths$raw,
    "\nSummaries:        ", sim_paths$summaries,
    "\nTables:           ", sim_paths$tables,
    "\nFigures:          ", sim_paths$figures, "\n", sep = "")
