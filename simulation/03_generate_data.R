# 03_generate_data.R ----------------------------------------------------------
# Protocol sections 9, 13 and 14: generate predictor variables, the fixed test
# population, and one training dataset for a given scenario + repetition seed.

# --- Generate (X_A, X_B) ~ MVN(0, Sigma), plus independent Z1, Z2 -------------
generate_predictors <- function(n, rho) {
  Sigma <- matrix(c(1, rho,
                    rho, 1), nrow = 2, byrow = TRUE)
  x <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  X_A <- x[, 1]
  X_B <- x[, 2]
  R <- X_A - X_B
  S <- X_A + X_B
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  data.frame(X_A = X_A, X_B = X_B, R = R, S = S, Z1 = Z1, Z2 = Z2)
}

# --- Attach true probabilities and binary outcomes ---------------------------
# logit{P(Y = 1)} = alpha + theta*R + delta*S + gamma_1*Z1 + gamma_2*Z2
generate_outcomes <- function(pred, theta, delta, gamma_1, gamma_2, alpha) {
  eta <- alpha + theta * pred$R + delta * pred$S +
         gamma_1 * pred$Z1 + gamma_2 * pred$Z2
  pred$p_true <- plogis(eta)
  pred$Y <- rbinom(nrow(pred), 1, pred$p_true)
  pred
}

# --- One training dataset (protocol section 14) -------------------------------
# Repetition-specific seed is set by the caller (run_one_repetition) via
# set.seed(seed_id) BEFORE this call, so the sequence is deterministic.
generate_train_data <- function(scenario, seed) {
  set.seed(seed)
  pred <- generate_predictors(scenario$n_train, scenario$rho)
  train <- generate_outcomes(pred, scenario$theta, scenario$delta,
                             scenario$gamma_1, scenario$gamma_2, scenario$alpha)
  # Immediate data checks (protocol section 14)
  stopifnot(all(is.finite(train$X_A)))
  stopifnot(all(is.finite(train$X_B)))
  stopifnot(all(train$p_true > 0 & train$p_true < 1))
  train
}

# --- Fixed test population (protocol section 13) ------------------------------
# One deterministic, scenario-specific test population of n_test subjects,
# generated from the same data-generating mechanism as training data.
generate_test_population <- function(scenario, seed) {
  set.seed(seed)
  pred <- generate_predictors(scenario$n_test, scenario$rho)
  test_pop <- generate_outcomes(pred, scenario$theta, scenario$delta,
                                scenario$gamma_1, scenario$gamma_2, scenario$alpha)
  stopifnot(all(is.finite(test_pop$X_A)))
  stopifnot(all(is.finite(test_pop$X_B)))
  stopifnot(all(test_pop$p_true > 0 & test_pop$p_true < 1))
  test_pop
}
