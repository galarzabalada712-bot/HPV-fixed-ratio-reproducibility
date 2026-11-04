# Undercoverage diagnosis — CORE_D0_N1000_RHO07 (final, frozen)

**Date:** after Stage 4 (M=1000 × 24 scenarios) and a targeted diagnostic run
(2 scenarios × M=5000, L'Ecuyer-CMRG, master seed 20260731).
**Data:** `output/summaries/undercoverage_diagnosis.csv`.

## Question

Stage 4 showed, for the null scenario CORE_D0_N1000_RHO07 (n=1000, rho=0.7):
rejection 0.079, Wald coverage 0.921, SE ratio 0.935. Was this a real
structural problem (data generation / fitting / testing code) or a finite-
sample estimation effect?

## Diagnostic design

- Target: CORE_D0_N1000_RHO07 (n=1000, rho=0.7, null)
- Control: CORE_D0_N1000_RHO02 (n=1000, rho=0.2, null)
- M = 5000 repetitions per scenario, same fixed test population and the
  SAME RNG kind as the main run (RNGkind("L'Ecuyer-CMRG")).
- Per repetition: Z_D = (D_hat - D_true)/se_D_hat; Wald, LRT (anova, test="LRT")
  and score (Rao) tests for H0: D = 0; theta/delta reparameterization
  (theta=(bA-bB)/2, delta=(bA+bB)/2) with model SEs; vcov components
  Var(bA), Var(bB), Cov(bA,bB); information-matrix condition number kappa.
- Profile-likelihood CI coverage for a null scenario is exactly
  1 - LRT rejection rate (profile CI contains D_true=0 iff the LRT accepts),
  so it is reported as such.

## Results

| metric | control rho=0.2 | target rho=0.7 |
| --- | --- | --- |
| Wald rejection | 0.0504 (MCSE 0.0031) | 0.0576 (MCSE 0.0033) |
| LRT rejection | 0.0506 | 0.0584 |
| Score (Rao) rejection | 0.0506 | 0.0580 |
| profile coverage (=1-LRT) | 0.9494 | 0.9416 |
| SE ratio (mean model SE / empirical SD) | 0.994 | 0.969 |
| Z_D: mean / SD / skew / exkurt | 0.005 / 1.002 / -0.057 / 0.066 | -0.015 / 1.027 / 0.077 / 0.092 |
| P(Z<-1.96) / P(Z>1.96) | 0.027 / 0.024 | 0.030 / 0.028 |
| bias(D_hat) | 0.0006 | -0.0013 |
| theta: empirical SD / mean model SE | 0.0699 / 0.0690 | 0.1044 / 0.1063 |
| delta: empirical SD / mean model SE | 0.0529 / 0.0526 | 0.0448 / 0.0435 |
| Cor(theta_hat, delta_hat) | -0.004 | 0.002 |
| kappa (info matrix) mean / p95 | 2.7 / 3.1 | 6.5 / 7.4 |

## Conclusions (evidence-supported only)

1. **Not a Wald-specific problem.** Wald / LRT / score rejection rates are
   nearly identical (0.0576 / 0.0584 / 0.0580). Switching tests or using a
   profile-likelihood interval would not remove the mild deviation.
2. **Not estimation bias.** E(Z_D) ~ -0.015 and bias(D_hat) is negligible; no
   central shift that could explain the elevation.
3. **Not skewness or heavy tails.** SD(Z_D)=1.027, skewness 0.077, excess
   kurtosis 0.092; left/right tail probabilities 0.030 / 0.028 are symmetric.
   Consistent with a statistic whose scale is slightly larger than standard
   normal, not with a clearly non-normal distribution.
4. **Mainly a modest underestimation of the model-based SE.**
   mean(se_D_hat)/SD(D_hat) = 0.969, i.e. ~3.1% underestimation relative to
   the empirical SD. The two-sided rejection rate predicted from
   SD(Z_D)=1.027 is 2*Phi(-1.96/1.027) ≈ 0.057, matching the observed 0.058.
5. **The underestimation sits in the sum direction, not the nuisance
   direction.** theta_hat SE is accurate or slightly conservative; delta_hat
   SE is underestimated by ~3%; Cor(theta_hat, delta_hat) ~ 0; kappa = 6.5 is
   elevated but not pathological. The earlier working hypothesis ("instability
   of theta_hat at high rho drives D") is NOT supported by the decomposition.

## Reinterpretation of the Stage 4 numbers

- The diagnostic M=5000 estimate of rejection is ~0.058, about
  **0.8 percentage points** above nominal 0.05 (relative +16%).
- Stage 4's 0.079 (M=1000) exceeds the diagnostic estimate by ~0.021
  (≈ 2.5-2.8 MCSE under M=1000), consistent with an unfavorable Monte Carlo
  realization on top of a real but small elevation.
- The mechanism is described as "consistent with" a small finite-sample
  underestimation of the model-based variance in the sum direction — NOT
  asserted as a proven second-order information-matrix effect (no theoretical
  derivation was performed).

## Recommended manuscript wording

> In the null scenario with n=1,000 and rho=0.7, the Stage 4 main-table
> estimate of the constraint-test rejection rate was 0.079, with Wald
> coverage 0.921 and an SE ratio of 0.935. An extended 5,000-replicate
> diagnostic analysis estimated the underlying rejection rate at
> approximately 0.058 (MCSE 0.003; Wald, likelihood-ratio and score tests
> all ~0.058), and an independently seeded 2,000-replicate verification run
> gave 0.071 (MCSE 0.006), with Wald coverage 0.93-0.94 in both. The
> model-based standard error was modestly underestimated relative to the
> empirical standard deviation (approximately 3% in the diagnostic run and
> 6% in the verification run), and the standardized statistic remained
> approximately normally distributed with no material bias, skewness, heavy
> tails, or asymmetry between the two rejection tails. These findings were
> consistent with a mild finite-sample underestimation of variance in the
> sum direction rather than a failure of a particular testing procedure.

**Note on manuscript consistency:** the main table MUST keep the Stage 4
value (0.079) — it is the reported result of the primary run. The diagnostic
(0.058) and verification (0.071) numbers are supplementary evidence that the
elevation is partly Monte Carlo variation; they must never replace 0.079 in
the main table, and the text should attribute each number to its source run.

## Suggested short statement for the results section

> Type I error was generally close to the nominal 0.05 level. A mild
> elevation was observed in the null scenario with n=1,000 and rho=0.7
> (Stage 4 rejection rate 0.079). An extended 5,000-replicate diagnostic
> analysis and an independently seeded 2,000-replicate verification run
> estimated the underlying rejection rate at approximately 0.058 and 0.071,
> respectively, across the Wald, likelihood-ratio, and score tests. This
> deviation was accompanied by approximately 3-6% underestimation of the
> empirical standard error but no meaningful bias or departure from an
> approximately normal sampling distribution.

Detailed Z-distribution diagnostics, tail probabilities, condition numbers,
and the (theta, delta) decomposition are appropriate for supplementary
material.

## RNG note (audit)

The diagnostic module MUST run under RNGkind("L'Ecuyer-CMRG") to reproduce the
Stage 4 repetition stream (per-repetition set.seed calls are RNG-kind
dependent). A first diagnostic run under the default Mersenne-Twister produced
non-comparable numbers (rejection 0.048) and was discarded. The module forces
the RNG kind internally; repetition-level results were verified identical to
the Stage 4 raw file for rep 1 (D_hat = -0.04667061, lrt_p = 0.5900379).
