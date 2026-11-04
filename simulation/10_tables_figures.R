# 10_tables_figures.R ---------------------------------------------------------
# Protocol sections 31 and 32: main and supplementary tables, and the
# three-panel main figure.

# --- Scenario short codes -----------------------------------------------------
# C1..C16 (CORE), M1..M4 (MECH), S1..S4 (STRONG), in scenario_id order. Used as
# the first column of the main and complete supplementary tables; the original
# verbose labels are kept as a trailing Scenario_label column for auditing.
scenario_code <- function(scenario_group, scenario_id) {
  prefix <- c(CORE = "C", MECH = "M", STRONG = "S")
  out <- character(length(scenario_group))
  for (g in c("CORE", "MECH", "STRONG")) {
    idx <- which(scenario_group == g)
    idx <- idx[order(scenario_id[idx])]
    out[idx] <- sprintf("%s%d", prefix[[g]], seq_along(idx))
  }
  out
}

# --- Main table: streamlined headline results ----------------------------------
# Reviewer-guided column set, kept deliberately lean for the main text:
# scenario code (C1..S4), n / rho, violation (D / V_LP), inferential headline
# (rejection rate, bias of D, coverage), and the single headline predictive
# metric ΔAUC. ΔBrier, ΔLogLoss and all other metrics live in the complete
# supplementary table (make_supplementary_complete).
make_main_table <- function(inferential, predictive, out_dir = sim_paths$tables) {
  main <- merge(
    inferential[, c("scenario_id", "scenario_label", "scenario_group",
                    "n_train", "rho", "D_true", "V_LP",
                    "rejection_rate", "bias_D", "coverage")],
    predictive[, c("scenario_id", "delta_auc")],
    by = "scenario_id", sort = FALSE
  )
  main <- main[order(main$scenario_id), , drop = FALSE]
  tab <- data.frame(
    Scenario = scenario_code(main$scenario_group, main$scenario_id),
    n = main$n_train,
    rho = round(main$rho, 1),
    D = round(main$D_true, 1),
    V_LP = round(main$V_LP, 3),
    Rejection_rate = round(main$rejection_rate, 3),
    Bias_of_D = round(main$bias_D, 3),
    Coverage = round(main$coverage, 3),
    Delta_AUC = round(main$delta_auc, 4),
    Scenario_label = main$scenario_label
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(tab, file.path(out_dir, "main_table.csv"), row.names = FALSE)
  invisible(tab)
}

# --- Supplementary complete table (single table with ALL results) --------------
# One complete simulation results table for the supplement: every scenario
# (all 24) with the full inferential block (D estimation, SE behaviour,
# coverage, LRT rejection, failure rate) AND the full predictive block
# (baseline / ratio / component AUC, calibration, net benefit). Merged by
# scenario_id so each scenario contributes exactly one row.
make_supplementary_complete <- function(inferential, predictive,
                                        out_dir = sim_paths$tables) {
  complete <- merge(
    inferential[, c("scenario_id", "scenario_label", "scenario_group",
                    "D_true", "bias_D", "ese_D", "mean_model_se",
                    "se_ratio", "coverage", "rejection_rate",
                    "rejection_mcse", "failure_rate")],
    predictive[, c("scenario_id", "auc_covariate", "auc_ratio",
                   "auc_component", "delta_auc", "delta_brier",
                   "delta_logloss",
                   "cal_intercept_ratio", "cal_intercept_component",
                   "cal_slope_ratio", "cal_slope_component",
                   "delta_NB_010", "delta_NB_020", "delta_NB_030")],
    by = "scenario_id", sort = FALSE
  )
  complete <- complete[order(complete$scenario_id), , drop = FALSE]
  tab <- data.frame(
    Scenario = scenario_code(complete$scenario_group, complete$scenario_id),
    Group = complete$scenario_group,
    D_true = round(complete$D_true, 3),
    Mean_D_estimate = round(complete$D_true + complete$bias_D, 3),
    Bias = round(complete$bias_D, 3),
    Empirical_SE = round(complete$ese_D, 3),
    Mean_model_SE = round(complete$mean_model_se, 3),
    SE_ratio = round(complete$se_ratio, 3),
    Coverage = round(complete$coverage, 3),
    LRT_rejection = round(complete$rejection_rate, 3),
    Monte_Carlo_SE = round(complete$rejection_mcse, 3),
    Failure_rate = round(complete$failure_rate, 4),
    Covariate_only_AUC = round(complete$auc_covariate, 3),
    Ratio_AUC = round(complete$auc_ratio, 3),
    Component_AUC = round(complete$auc_component, 3),
    Delta_AUC = round(complete$delta_auc, 4),
    Delta_Brier = round(complete$delta_brier, 4),
    Delta_LogLoss = round(complete$delta_logloss, 4),
    Cal_intercept_ratio = round(complete$cal_intercept_ratio, 3),
    Cal_intercept_component = round(complete$cal_intercept_component, 3),
    Cal_slope_ratio = round(complete$cal_slope_ratio, 3),
    Cal_slope_component = round(complete$cal_slope_component, 3),
    Delta_NB_010 = round(complete$delta_NB_010, 4),
    Delta_NB_020 = round(complete$delta_NB_020, 4),
    Delta_NB_030 = round(complete$delta_NB_030, 4),
    Scenario_label = complete$scenario_label
  )
  write.csv(tab, file.path(out_dir, "supplementary_complete.csv"),
            row.names = FALSE)
  invisible(tab)
}

# --- Main figure: three panels ------------------------------------------------
# Reviewer guidance:
#   Panel A — constraint-test rejection rate vs V_LP (0.05 reference line)
#   Panel B — coverage of the Wald CI for D vs V_LP (0.95 reference line)
#   Panel C — mean test-population ΔAUC vs V_LP (0 reference line)
# ΔBrier, ΔLogLoss and ΔNB go to a supplementary figure so the main figure
# stays readable and the conclusions do not rest on a single metric.
#
# The main figure shows CORE and STRONG scenarios only: they share theta = 0.5,
# so V_LP (a standardized effect) is on a comparable scale and connecting them
# with lines is legitimate. MECH scenarios use different theta (0.25 / 0.15)
# and are displayed separately as points in make_mechanistic_figure() — they
# must NOT be joined to the CORE/STRONG trend lines, which would imply a
# continuous trend across different data-generating mechanisms.
make_main_figure <- function(inferential, predictive,
                             out_dir = sim_paths$figures) {
  df <- merge(
    inferential[inferential$scenario_group %in% c("CORE", "STRONG"),
                c("scenario_id", "scenario_label", "scenario_group",
                  "baseline_strength", "delta", "V_LP", "rho", "n_train",
                  "rejection_rate", "rejection_mcse", "coverage",
                  "coverage_mcse")],
    predictive[predictive$scenario_group %in% c("CORE", "STRONG"),
               c("scenario_id", "delta_auc", "delta_auc_mcse")],
    by = "scenario_id", sort = FALSE
  )
  df$n_label <- factor(df$n_train, levels = c(1000, 5000),
                       labels = c("n = 1000", "n = 5000"))
  df$rho_label <- factor(df$rho, levels = c(0.2, 0.7),
                         labels = c("rho = 0.2", "rho = 0.7"))
  df$strength_label <- factor(df$baseline_strength,
                              levels = c("moderate", "strong"))

  panel_theme <- theme_bw(base_size = 9) + theme(legend.position = "bottom")

  # Panel A: rejection rate vs standardized violation magnitude
  pA <- ggplot(df, aes(x = V_LP, y = rejection_rate,
                       color = n_label, linetype = strength_label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey40") +
    facet_wrap(~ rho_label) +
    labs(x = "Standardized violation magnitude (V_LP)",
         y = "Constraint-test rejection rate",
         color = NULL, linetype = NULL) +
    panel_theme

  # Panel B: coverage of the Wald CI for D
  pB <- ggplot(df, aes(x = V_LP, y = coverage,
                       color = n_label, linetype = strength_label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0.95, linetype = "dotted", color = "grey40") +
    facet_wrap(~ rho_label) +
    labs(x = "Standardized violation magnitude (V_LP)",
         y = "Coverage of 95% Wald CI for D",
         color = NULL, linetype = NULL) +
    panel_theme

  # Panel C: mean test-population delta AUC vs V_LP
  pC <- ggplot(df, aes(x = V_LP, y = delta_auc,
                       color = n_label, linetype = strength_label)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    facet_wrap(~ rho_label) +
    labs(x = "Standardized violation magnitude (V_LP)",
         y = expression(paste("Mean ", Delta, "AUC")),
         color = NULL, linetype = NULL) +
    panel_theme

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(out_dir, "main_figure.png"),
      width = 11, height = 4.2, units = "in", res = 300)
  # Compose the three panels with base grid (no gridExtra dependency)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 3)))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pC, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid::popViewport()
  dev.off()
  cat("Saved main figure to", file.path(out_dir, "main_figure.png"), "\n")
  invisible(list(panel_A = pA, panel_B = pB, panel_C = pC))
}

# --- Mechanistic sensitivity figure -------------------------------------------
# MECH scenarios (one informative component; same-direction effects) use
# theta = 0.25 / 0.15, different from the CORE/STRONG theta = 0.5. They are
# mechanistic illustrations rather than points on the same V_LP trend, so they
# are plotted as isolated points (no connecting lines), colored by mechanism.
make_mechanistic_figure <- function(inferential, predictive,
                                    out_dir = sim_paths$figures) {
  df <- merge(
    inferential[inferential$scenario_group == "MECH",
                c("scenario_id", "scenario_label", "mechanism",
                  "theta", "delta", "D_true", "V_LP", "rho",
                  "rejection_rate", "rejection_mcse")],
    predictive[predictive$scenario_group == "MECH",
               c("scenario_id", "delta_auc", "delta_auc_mcse")],
    by = "scenario_id", sort = FALSE
  )
  df$rho_label <- factor(df$rho, levels = c(0.2, 0.7),
                         labels = c("rho = 0.2", "rho = 0.7"))
  df$mech_label <- factor(df$mechanism,
                          levels = c("one informative component",
                                     "same-direction component effects"))

  panel_theme <- theme_bw(base_size = 9) + theme(legend.position = "bottom")

  pA <- ggplot(df, aes(x = V_LP, y = rejection_rate, color = mech_label)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey40") +
    facet_wrap(~ rho_label) +
    labs(x = "Standardized violation magnitude (V_LP)",
         y = "Constraint-test rejection rate",
         color = NULL) +
    panel_theme

  pB <- ggplot(df, aes(x = V_LP, y = delta_auc, color = mech_label)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    facet_wrap(~ rho_label) +
    labs(x = "Standardized violation magnitude (V_LP)",
         y = expression(paste("Mean ", Delta, "AUC")),
         color = NULL) +
    panel_theme

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(out_dir, "mechanistic_figure.png"),
      width = 8.5, height = 3.6, units = "in", res = 300)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(pA, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pB, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
  dev.off()
  cat("Saved mechanistic figure to", file.path(out_dir, "mechanistic_figure.png"), "\n")
  invisible(list(panel_A = pA, panel_B = pB))
}

# --- Supplementary figure: ΔBrier, ΔLogLoss, ΔNB vs V_LP ----------------------
# Also restricted to CORE + STRONG (theta = 0.5): connecting points with lines
# is only meaningful within a shared data-generating mechanism. MECH scenarios
# are shown as isolated points in mechanistic_figure.png.
make_supplementary_figure <- function(inferential, predictive,
                                      out_dir = sim_paths$figures) {
  df <- merge(
    inferential[inferential$scenario_group %in% c("CORE", "STRONG"),
                c("scenario_id", "scenario_label", "scenario_group",
                  "baseline_strength", "V_LP", "rho", "n_train")],
    predictive[predictive$scenario_group %in% c("CORE", "STRONG"),
               c("scenario_id", "delta_brier", "delta_logloss",
                 "delta_NB_010", "delta_NB_020", "delta_NB_030")],
    by = "scenario_id", sort = FALSE
  )
  df$n_label <- factor(df$n_train, levels = c(1000, 5000),
                       labels = c("n = 1000", "n = 5000"))
  df$rho_label <- factor(df$rho, levels = c(0.2, 0.7),
                         labels = c("rho = 0.2", "rho = 0.7"))
  df$strength_label <- factor(df$baseline_strength,
                              levels = c("moderate", "strong"))

  panel_theme <- theme_bw(base_size = 9) + theme(legend.position = "bottom")
  base <- function(yvar, ylab) {
    ggplot(df, aes(x = V_LP, y = .data[[yvar]],
                   color = n_label, linetype = strength_label)) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.8) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
      facet_wrap(~ rho_label) +
      labs(x = "Standardized violation magnitude (V_LP)", y = ylab,
           color = NULL, linetype = NULL) +
      panel_theme
  }
  pD <- base("delta_brier", expression(paste("Mean ", Delta, "Brier")))
  pE <- base("delta_logloss", expression(paste("Mean ", Delta, "log loss")))
  pF <- base("delta_NB_010", expression(paste("Mean ", Delta, "NB (p_t = 0.10)")))
  pG <- base("delta_NB_020", expression(paste("Mean ", Delta, "NB (p_t = 0.20)")))
  pH <- base("delta_NB_030", expression(paste("Mean ", Delta, "NB (p_t = 0.30)")))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(out_dir, "supplementary_figure.png"),
      width = 11, height = 8.4, units = "in", res = 300)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
  print(pD, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pE, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(pF, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(pG, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pH, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid::popViewport()
  dev.off()
  cat("Saved supplementary figure to", file.path(out_dir, "supplementary_figure.png"), "\n")
  invisible(list(p_D = pD, p_E = pE, p_F = pF, p_G = pG, p_H = pH))
}
