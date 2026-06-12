# Phase 8 — Sample-overlap bias bound for the EAS urate->gout arm, WITHOUT routing
# through LDSC. Addresses reviewer major point #3: CKDGen trans-ethnic urate (exposure)
# includes BioBank Japan and the Major-2024 EAS gout meta (outcome) includes Japanese
# cohorts, so the EAS arm has uncertain sample overlap; the LDSC cross-trait intercept
# (|C_off_diag| ~ 0.025) is unreliable for an oligogenic trait dominated by SLC2A9/ABCG2,
# so "the LDSC overlap is small" cannot be load-bearing.
#
# Bound (does NOT use LDSC). Overlap biases two-sample MR toward the one-sample
# estimate, which is biased toward the confounded OBSERVATIONAL association with
# relative magnitude ~1/F per instrument (Burgess, Davies & Thompson 2016, Genet
# Epidemiol 40:597-608). With overlap proportion rho:
#       |bias|  <=  rho * |beta_obs - beta_MR| / Fbar
# We sweep rho up to the absurd worst case (rho=1, complete overlap) and a plausible
# range of the observational per-1-SD-urate log-odds for gout, and stress Fbar down to
# the median / a pessimistic floor. The bound is compared to the cross-ancestry gap
# |beta_EAS - beta_EUR|.
#
# Writes: outputs/tables/overlap_bias_bound.csv

source(here::here("R", "utils.R"))
suppressPackageStartupMessages(library(data.table))

bA <- 2.29088311690581           # EAS urate->gout MR-APSS log-odds per 1-SD
bE <- 1.21267213415297           # EUR
gap <- bA - bE                   # cross-ancestry difference to be explained (~0.82-1.08 incl. IVW)

# EAS urate instrument strength (computed in R/14 inputs: mean F 126.7, median 49.7).
Fbar_mean   <- 126.7
Fbar_median <- 49.7
Fbar_floor  <- 30.0              # pessimistic stress (near the per-SNP min)

# Observational per-1-SD-urate log-odds for gout. Strong but uncertain; sweep a wide,
# generous range (OR ~2x to ~6x per SD => log-odds ~0.69 to ~1.79).
beta_obs_grid <- c(0.69, 1.10, 1.39, 1.79)   # log(2), log(3), log(4), log(6)
rho_grid      <- c(0.10, 0.25, 0.50, 1.00)   # overlap proportion, up to complete-overlap worst case

bound <- function(rho, beta_obs, Fbar) rho * abs(beta_obs - bA) / Fbar

rows <- list()
for (Fb in c(Fbar_mean, Fbar_median, Fbar_floor))
  for (bo in beta_obs_grid)
    for (rho in rho_grid)
      rows[[length(rows) + 1]] <- data.table(
        Fbar = Fb, beta_obs = bo, rho = rho,
        max_bias = bound(rho, bo, Fb),
        pct_of_gap = 100 * bound(rho, bo, Fb) / gap)
res <- rbindlist(rows)
fwrite(res, here::here("outputs", "tables", "overlap_bias_bound.csv"))

worst <- res[which.max(max_bias)]
cat(sprintf("Cross-ancestry gap to explain: |beta_EAS - beta_EUR| = %.3f log-odds\n", gap))
cat(sprintf("EAS urate instrument strength: Fbar(mean)=%.0f, median=%.1f, floor stress=%.0f\n\n",
            Fbar_mean, Fbar_median, Fbar_floor))
cat("Max overlap-induced bias (LDSC-free; Burgess-Davies-Thompson 2016) over the sweep:\n")
cat("  Fbar    beta_obs  rho     |bias|   % of gap\n")
show <- res[Fbar == Fbar_median][order(-max_bias)]
for (i in seq_len(nrow(show))) with(show[i], cat(sprintf(
  "  %-6.1f  %-7.2f   %-5.2f   %.4f   %.1f%%\n", Fbar, beta_obs, rho, max_bias, pct_of_gap)))
cat(sprintf("\n  WORST CASE over ALL Fbar/beta_obs/rho (incl. rho=1, F-floor=%.0f):\n", Fbar_floor))
cat(sprintf("    |bias| = %.4f log-odds = %.1f%% of the cross-ancestry gap.\n",
            worst$max_bias, worst$pct_of_gap))
cat("\nReading: even under complete overlap (rho=1) and a pessimistic F floor, the\n")
cat("overlap-induced bias is a small fraction of the EAS-EUR gap, and it does NOT\n")
cat("route through LDSC. Overlap cannot generate the urate->gout cross-ancestry\n")
cat("exception. (Definitive follow-up: re-estimate the EAS arm on a Major-2024 gout\n")
cat("stratum with no BBJ contribution, if one is separable.)\n")
