# 01_scenario_grid.R ----------------------------------------------------------
# Protocol sections 3 and 8: build the complete 24-scenario grid
# (16 core + 4 mechanistic sensitivity + 4 strong-baseline sensitivity).

make_scenario_label <- function(group, delta, n_train, rho) {
  dcode <- ifelse(delta == 0, "D0", sprintf("D%03d", round(delta * 100)))
  sprintf("%s_%s_N%d_RHO%02d", group, dcode, n_train, round(rho * 10))
}

# --- 16 core scenarios: theta = 0.5, delta in {0, .05, .2, .4},
#     n_train in {1000, 5000}, rho in {0.2, 0.7} -------------------------------
build_core_scenarios <- function(n_test = 100000L, n_repetitions = 1000L) {
  core <- expand.grid(
    delta   = c(0, 0.05, 0.20, 0.40),
    n_train = c(1000L, 5000L),
    rho     = c(0.2, 0.7),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  core$theta  <- 0.5
  core$beta_A <- core$theta + core$delta          # theta + delta
  core$beta_B <- -core$theta + core$delta         # -theta + delta
  core$D_true <- 2 * core$delta                   # beta_A + beta_B
  core$scenario_group    <- "CORE"
  core$baseline_strength <- "moderate"
  core$mechanism         <- "fixed-ratio design; D = 2*delta"
  core$scenario_label    <- make_scenario_label("CORE", core$delta, core$n_train, core$rho)
  core
}

# --- 4 mechanistic sensitivity scenarios (n = 5000) --------------------------
# One informative component:  (beta_A, beta_B) = (0.5, 0)   -> theta=0.25, delta=0.25, D=0.5
# Same-direction effects:     (beta_A, beta_B) = (0.5, 0.2) -> theta=0.15, delta=0.35, D=0.7
# These are mechanistic illustrations, NOT distinct forms of constraint violation.
build_mechanistic_scenarios <- function(n_test = 100000L, n_repetitions = 1000L) {
  configs <- list(
    list(tag = "ONECOMP", beta_A = 0.5, beta_B = 0.0, desc = "one informative component"),
    list(tag = "SAMEDIR", beta_A = 0.5, beta_B = 0.2, desc = "same-direction component effects")
  )
  rows <- vector("list", 4L)
  k <- 0L
  for (cfg in configs) {
    for (rho in c(0.2, 0.7)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        delta   = (cfg$beta_A + cfg$beta_B) / 2,
        n_train = 5000L,
        rho     = rho,
        theta   = (cfg$beta_A - cfg$beta_B) / 2,
        beta_A  = cfg$beta_A,
        beta_B  = cfg$beta_B,
        D_true  = cfg$beta_A + cfg$beta_B,
        scenario_group    = "MECH",
        baseline_strength = "moderate",
        mechanism         = cfg$desc,
        scenario_label    = sprintf("MECH_%s_RHO%02d", cfg$tag, round(rho * 10)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# --- 4 strong-baseline scenarios (n = 5000, covariate-only AUC ~ 0.85) -------
build_strong_scenarios <- function(n_test = 100000L, n_repetitions = 1000L) {
  strong <- expand.grid(
    delta = c(0.05, 0.20),
    rho   = c(0.2, 0.7),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  strong$theta  <- 0.5
  strong$beta_A <- strong$theta + strong$delta
  strong$beta_B <- -strong$theta + strong$delta
  strong$D_true <- 2 * strong$delta
  strong$scenario_group    <- "STRONG"
  strong$baseline_strength <- "strong"
  strong$mechanism         <- "strong baseline model (covariate-only AUC ~ 0.85)"
  strong$n_train           <- 5000L
  strong$scenario_label    <- make_scenario_label("STRONG", strong$delta, strong$n_train, strong$rho)
  strong
}

# --- Assemble the full grid --------------------------------------------------
build_scenario_grid <- function(n_test = 100000L, n_repetitions = 1000L) {
  grid <- rbind(
    build_core_scenarios(n_test, n_repetitions),
    build_mechanistic_scenarios(n_test, n_repetitions),
    build_strong_scenarios(n_test, n_repetitions)
  )
  grid$scenario_id      <- seq_len(nrow(grid))
  grid$prevalence_target <- 0.25
  grid$n_test           <- n_test
  grid$n_repetitions    <- n_repetitions
  # Standardized violation magnitude (protocol section 12):
  #   V_LP = |delta| * sqrt(Var(S | R, Z)) = |delta| * sqrt(2*(1+rho))
  # Equivalently V_LP = |D| * sqrt((1+rho)/2) since D = 2*delta
  #   (2*|delta|*sqrt((1+rho)/2) = |delta|*sqrt(2*(1+rho))).
  # This is the standardized effect of the omitted sum-direction signal
  # (expected log-odds change per 1 SD of S). Both expressions are computed
  # and checked for numeric agreement in QC 29.2.
  grid$V_LP  <- abs(grid$delta) * sqrt(2 * (1 + grid$rho))
  grid$V_LP2 <- abs(grid$D_true) * sqrt((1 + grid$rho) / 2)
  # Calibration columns (filled by 02_calibrate_parameters.R)
  grid$gamma_1 <- NA_real_
  grid$gamma_2 <- NA_real_
  grid$alpha   <- NA_real_
  grid[, c(
    "scenario_id", "scenario_label", "scenario_group", "mechanism",
    "baseline_strength", "theta", "delta", "beta_A", "beta_B", "D_true",
    "rho", "n_train", "n_test", "n_repetitions", "prevalence_target",
    "V_LP", "V_LP2", "gamma_1", "gamma_2", "alpha"
  )]
}
