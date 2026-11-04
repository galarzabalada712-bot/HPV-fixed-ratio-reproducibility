# 04_fit_models.R -------------------------------------------------------------
# Protocol sections 15-17: fit the competing models, run the likelihood-ratio
# test, and estimate D = beta_A + beta_B with its uncertainty.

# Fixed-ratio model: imposes beta_A = -beta_B  (logit = a + lambda*R + gammas)
fit_ratio_model <- function(train) {
  glm(Y ~ R + Z1 + Z2, family = binomial(), data = train)
}

# Unconstrained-component model (nests the ratio model)
fit_component_model <- function(train) {
  glm(Y ~ X_A + X_B + Z1 + Z2, family = binomial(), data = train)
}

# Covariate-only model (documents actual baseline-model strength)
fit_covariate_model <- function(train) {
  glm(Y ~ Z1 + Z2, family = binomial(), data = train)
}

# --- Fit-success check (protocol section 24) ----------------------------------
# A fit is unsuccessful if: glm errors, fitted probabilities are exactly 0 or 1,
# coefficients / covariance are non-finite.
check_fit <- function(fit) {
  if (is.null(fit)) return("glm returned NULL (error)")
  if (!isTRUE(fit$converged)) return("glm did not converge")
  co <- coef(fit)
  if (any(!is.finite(co))) return("non-finite coefficients")
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(vc) || any(!is.finite(vc))) return("non-finite covariance matrix")
  fp <- fitted(fit)
  if (any(fp <= 1e-8 | fp >= 1 - 1e-8)) return("fitted probabilities at 0 or 1")
  NULL
}

# --- 1-df likelihood-ratio test (protocol section 16) -------------------------
compute_lrt <- function(fit_ratio, fit_component) {
  lrt <- anova(fit_ratio, fit_component, test = "LRT")
  statistic <- lrt$Deviance[2]
  df <- lrt$Df[2]
  p <- lrt[2, "Pr(>Chi)"]
  if (!identical(as.integer(df), 1L)) {
    stop("LRT degrees of freedom is not 1 (got ", df, ")")
  }
  list(lrt_statistic = statistic, lrt_df = df, lrt_p = p)
}

# --- Estimate D and its Wald CI (protocol section 17) -------------------------
estimate_D <- function(fit_component, D_true) {
  b_A_hat <- unname(coef(fit_component)["X_A"])
  b_B_hat <- unname(coef(fit_component)["X_B"])
  D_hat <- b_A_hat + b_B_hat
  V <- vcov(fit_component)
  var_D_hat <- V["X_A", "X_A"] + V["X_B", "X_B"] + 2 * V["X_A", "X_B"]
  se_D_hat <- sqrt(var_D_hat)
  D_lower <- D_hat - qnorm(0.975) * se_D_hat
  D_upper <- D_hat + qnorm(0.975) * se_D_hat
  cover_D <- as.integer(D_lower <= D_true && D_upper >= D_true)
  list(b_A_hat = b_A_hat, b_B_hat = b_B_hat, D_hat = D_hat,
       se_D_hat = se_D_hat, D_lower = D_lower, D_upper = D_upper,
       cover_D = cover_D)
}
