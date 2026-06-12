# Phase 8 — Liability- and risk-difference-scale conversion of the urate->gout
# cross-ancestry contrast. Addresses reviewer major point #2: "~3-fold" is the
# OR-ratio (exp(2.291)/exp(1.213) ~ 2.9x); the causal slope is the log-odds beta,
# whose ratio is ~1.9x; and log-odds is non-collapsible with a prevalence-dependent
# scaling. The reviewer asks: convert to a liability / risk-difference scale using
# ancestry-specific gout prevalence and show how much of the OR gap survives.
#
# Liability-threshold mapping. With population prevalence K, threshold t=Phi^{-1}(1-K)
# and density phi(t), a per-1-SD-exposure log-odds beta maps (to first order) to a
# per-1-SD liability-scale shift  delta = beta * K(1-K) / phi(t)  [probit slope].
# Risk-difference at baseline K:  RD = K' - K, K' = OR*K/(1-K) / (1 + OR*K/(1-K)),
# OR = exp(beta).
#
# Population gout prevalence is NOT in the GWAS (the harmonised prevalence.outcome
# field is a 0.1 placeholder, identical across arms), so K is SWEPT over literature-
# anchored ranges (EUR/Finnish gout ~1-3%; East-Asian gout generally reported higher,
# ~2-6%) and clearly labelled ASSUMED, not derived. The equal-K diagonal isolates the
# pure slope (delta-ratio == beta-ratio there).
#
# Writes: outputs/tables/liability_scale_urate_gout.csv

source(here::here("R", "utils.R"))
suppressPackageStartupMessages(library(data.table))

# MR-APSS primary urate->gout arms (outputs/tables/mr_apss_estimates.csv).
bE <- 1.21267213415297; bA <- 2.29088311690581
or_ratio  <- exp(bA) / exp(bE)
beta_ratio <- bA / bE
cat(sprintf("Scale-dependence of the magnitude (prevalence-independent):\n"))
cat(sprintf("  OR-ratio (the '~3-fold' headline) = %.2fx   [EUR OR %.2f, EAS OR %.2f]\n",
            or_ratio, exp(bE), exp(bA)))
cat(sprintf("  log-odds (causal-slope) ratio      = %.2fx   <- the defensible magnitude\n\n",
            beta_ratio))

g  <- function(K) K * (1 - K) / dnorm(qnorm(1 - K))   # log-odds -> liability factor
delta <- function(b, K) b * g(K)
rd <- function(b, K) {                                 # risk difference per 1-SD at K
  o0 <- K / (1 - K); o1 <- exp(b) * o0
  o1 / (1 + o1) - K
}

K_EUR <- c(0.01, 0.02, 0.03)
K_EAS <- c(0.02, 0.04, 0.06)
grid <- CJ(K_EUR = K_EUR, K_EAS = K_EAS)
grid[, `:=`(
  delta_EUR = delta(bE, K_EUR), delta_EAS = delta(bA, K_EAS),
  RD_EUR_pct = 100 * rd(bE, K_EUR), RD_EAS_pct = 100 * rd(bA, K_EAS))]
grid[, `:=`(delta_ratio = delta_EAS / delta_EUR, RD_ratio = RD_EAS_pct / RD_EUR_pct)]

# equal-prevalence reference (pure slope, removes any prevalence asymmetry)
eq <- data.table(K = K_EUR)[, .(K, delta_ratio_equalK = delta(bA, K) / delta(bE, K))]

fwrite(grid, here::here("outputs", "tables", "liability_scale_urate_gout.csv"))

cat("Liability-scale delta-ratio and risk-difference-ratio over ASSUMED gout prevalences:\n")
cat(sprintf("  %-9s %-9s | delta-ratio | RD_EUR%%  RD_EAS%%  RD-ratio\n","K_EUR","K_EAS"))
for (i in seq_len(nrow(grid))) with(grid[i], cat(sprintf(
  "  %-9.2f %-9.2f |   %.2fx    | %.2f    %.2f    %.2fx\n",
  K_EUR, K_EAS, delta_ratio, RD_EUR_pct, RD_EAS_pct, RD_ratio)))
cat(sprintf("\n  equal-prevalence diagonal (K_EUR=K_EAS): delta-ratio = beta-ratio = %.2fx\n", beta_ratio))
cat(sprintf("  delta-ratio range over the grid: %.2f - %.2fx\n",
            min(grid$delta_ratio), max(grid$delta_ratio)))

cat("\nReading:\n")
cat(" - The '~3-fold' is the OR-ratio; the causal log-odds ratio is ~1.9x; report that.\n")
cat(" - Baseline-prevalence NON-COLLAPSIBILITY works AGAINST the observed gap: when\n")
cat("   K_EAS > K_EUR the same liability shift yields a SMALLER log-odds, so a larger\n")
cat("   EAS log-odds implies an even larger liability-scale gap (delta-ratio >= beta-ratio).\n")
cat("   => the EAS excess is NOT a statistical-prevalence artefact.\n")
cat(" - A residual biological reading remains (a steeper near-MSU-saturation dose-response\n")
cat("   in higher-baseline-urate EAS) — effect-modification, not a measurement artefact —\n")
cat("   and should be stated as interpretation, separate from the non-collapsibility point.\n")
