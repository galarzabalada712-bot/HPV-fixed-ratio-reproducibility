# 00_packages.R ---------------------------------------------------------------
# Protocol section 6: required R packages and shared utilities.
# Core model fitting uses base R glm(); extra packages are kept minimal.
# data.table is deliberately NOT required (the protocol lists it as part of a
# "minimal" set but also states: "Do not make the simulation depend on
# unnecessary packages"). future/furrr are optional: if they are installed the
# scenario loop parallelizes; otherwise it falls back to sequential execution.

sim_required_packages <- c("MASS", "pROC", "dplyr", "tidyr", "purrr", "ggplot2")
sim_optional_packages <- c("future", "furrr")

# Central output paths. Master script rewrites sim_paths$root before sourcing
# so that relative/absolute locations can be adjusted in one place.
if (!exists("sim_paths", inherits = FALSE)) {
  sim_paths <- new.env(parent = emptyenv())
  sim_paths$root      <- "simulation"
  sim_paths$raw       <- file.path("simulation", "output", "raw")
  sim_paths$summaries <- file.path("simulation", "output", "summaries")
  sim_paths$tables    <- file.path("simulation", "output", "tables")
  sim_paths$figures   <- file.path("simulation", "output", "figures")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Fast AUC via the Mann-Whitney U statistic (rank based, no ties assumed).
# Mathematically identical to the trapezoidal AUC for continuous scores.
# Used ONLY inside the calibration search; every reported AUC uses pROC.
fast_auc <- function(score, y) {
  n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

load_simulation_packages <- function() {
  missing <- setdiff(sim_required_packages, rownames(installed.packages()))
  if (length(missing) > 0) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "),
         ". Install them before running the simulation.", call. = FALSE)
  }
  for (p in sim_required_packages) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  available <- intersect(sim_optional_packages, rownames(installed.packages()))
  for (p in available) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(list(optional_available = available))
}
